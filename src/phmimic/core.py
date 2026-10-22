"""Configuration, provenance and file IO.

Every stage goes through this module so that:
  * there is exactly one place a parameter is read from (config.yaml);
  * the data files are hash-checked before any computation touches them;
  * reading a FASTA and hashing a file behave identically everywhere.
"""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import subprocess
from dataclasses import dataclass, field

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[2]


class ProvenanceError(RuntimeError):
    """Raised when a data file does not match the recorded manifest."""


def sha256_file(path: pathlib.Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def rel(path: pathlib.Path) -> str:
    """Path relative to the project root, or absolute if it lives outside it."""
    try:
        return str(pathlib.Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


@dataclass
class Config:
    raw: dict
    path: pathlib.Path
    sha256: str

    # -------------------------------------------------------------- accessors
    def get(self, dotted: str, default=None):
        node = self.raw
        for key in dotted.split("."):
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node

    def path_of(self, dotted: str) -> pathlib.Path:
        value = self.get(dotted)
        if value is None:
            raise KeyError(f"config key not set: {dotted}")
        p = pathlib.Path(value)
        return p if p.is_absolute() else ROOT / p

    @property
    def reports(self) -> pathlib.Path:
        return ROOT / "reports"

    @property
    def logs(self) -> pathlib.Path:
        return ROOT / "logs"


def load_config(path: pathlib.Path | str | None = None) -> Config:
    cfg_path = pathlib.Path(path) if path else ROOT / "config.yaml"
    text = cfg_path.read_text(encoding="utf-8")
    return Config(raw=yaml.safe_load(text), path=cfg_path, sha256=sha256_text(text))


# ------------------------------------------------------------------- manifest
def read_manifest(cfg: Config) -> dict:
    """file basename -> manifest row (dict of strings)."""
    manifest = cfg.path_of("data.manifest")
    if not manifest.exists():
        manifest = cfg.path_of("data.cache_dir") / "MANIFEST.tsv"
    rows = {}
    with open(manifest, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            rows[row["file"]] = row
    return rows


def verify_data(cfg: Config) -> dict:
    """Hash every declared input; raise ProvenanceError on any mismatch."""
    rows = read_manifest(cfg)
    results = {}
    problems = []
    for key in ("data.query", "data.proteome"):
        path = cfg.path_of(key)
        name = path.name
        if not path.exists():
            problems.append(f"{name}: file missing")
            continue
        digest = sha256_file(path)
        row = rows.get(name)
        if row is None:
            problems.append(f"{name}: not listed in manifest")
            continue
        ok = digest == row["sha256"]
        if not ok:
            problems.append(f"{name}: sha256 {digest} != manifest {row['sha256']}")
        results[key] = dict(
            file=rel(path),
            sha256=digest,
            manifest_sha256=row["sha256"],
            ok=ok,
            uniprot_release=row.get("uniprot_release", ""),
            release_date=row.get("release_date", ""),
            source_url=row.get("source_url", ""),
            n_seqs=int(row.get("n_seqs") or 0),
        )
    expected = cfg.get("data.uniprot_release_expected")
    releases = {r["uniprot_release"] for r in results.values()}
    if expected and releases and releases != {expected}:
        problems.append(f"UniProt release {sorted(releases)} != expected {expected}")
    if problems and cfg.get("data.require_sha256", True):
        raise ProvenanceError("; ".join(problems))
    return dict(ok=not problems, problems=problems, inputs=results)


def tool_versions(cfg: Config) -> dict:
    """Record the version of every external tool the pipeline can call."""
    out = {}
    venv = ROOT / "env" / "venv" / "bin" / "python"
    try:
        import sys
        out["python"] = sys.version.split()[0]
    except Exception:  # pragma: no cover
        out["python"] = "unknown"
    for name, argv in (
        ("mmseqs", [str(ROOT / "env" / "bin" / "mmseqs"), "version"]),
        ("mimic_scan", [str(ROOT / "env" / "bin" / "mimic_scan")]),
    ):
        try:
            p = subprocess.run(argv, capture_output=True, text=True, timeout=30)
            text = (p.stdout or p.stderr).strip().splitlines()
            out[name] = text[0] if text else f"rc={p.returncode}"
        except (FileNotFoundError, subprocess.SubprocessError):
            out[name] = "not installed"
    try:
        p = subprocess.run(["gcc", "--version"], capture_output=True, text=True, timeout=20)
        out["gcc"] = p.stdout.splitlines()[0]
    except (FileNotFoundError, subprocess.SubprocessError):
        out["gcc"] = "not installed"
    out["venv_python"] = str(venv.relative_to(ROOT)) if venv.exists() else "absent"
    return out


# ----------------------------------------------------------------------- fasta
def read_fasta(path: pathlib.Path) -> list[tuple[str, str, str]]:
    """Return [(accession, gene, sequence)] in file order."""
    records: list[tuple[str, str, str]] = []
    header = None
    buf: list[str] = []

    def flush():
        if header is None:
            return
        acc = header.split("|")[1] if "|" in header else header.split()[0]
        gene = next((t[3:] for t in header.split() if t.startswith("GN=")), "")
        records.append((acc, gene, "".join(buf).upper()))

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    flush()
    return records


def write_json(path: pathlib.Path, payload) -> pathlib.Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def read_json(path: pathlib.Path):
    return json.loads(path.read_text(encoding="utf-8"))
