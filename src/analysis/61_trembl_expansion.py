"""A2: does the negative conclusion survive a 9x larger library?

Everything in the dossier was measured against 20,431 reviewed human proteins.
The first question any reviewer asks is "and TrEMBL?". This script answers it
with a measurement rather than an argument, and it decides what to measure
*before* looking at the answer.

Pre-registered questions
------------------------
Q1  Scope. At the s61 parameters, how many EBNA1 alignments survive the
    enlargement, what is the best E-value, and does any alignment reach E < 1?
    This is the number that widens or narrows the negative claim.

Q2  Rate, not count. The enlarged library has ~10x more proteins, so raw counts
    must rise; a rise by itself is meaningless. The informative comparison is
    the real EBNA1 against composition-matched shuffles measured **in the same
    enlarged library** -- that is the null re-fit the ROADMAP requires. What
    carries over is the ratio, not the count.

Q3  Published positives. Are CRYAB / ANO2 / MEF2B / GlialCAM matched by any
    entry in the enlarged library, fragment or not?

Q4  Composition. Does the Gly-Ala block (90-330) still collect zero overlapping
    alignments, or do TrEMBL fragments manufacture GA-ish hits that reviewed
    Swiss-Prot does not contain?

Q5  Distinct genes. TrEMBL contains fragments of proteins that are already in
    Swiss-Prot, so an accession count overstates what was added. The honest
    statistic is how many *different genes* are hit in each library.

Policies fixed in advance
-------------------------
Dedupe      drop unreviewed entries whose sequence is byte-identical to a
            reviewed entry (pure redundancy), and count how many were dropped.
            Near-duplicates are NOT clustered: clustering is itself a modelling
            choice and it would change the search space in a way this script
            cannot defend.
Isoforms    UniProtKB canonical entries only. Alternatively spliced isoform
            sequences live in a separate database and are absent here, so
            "protein" throughout means "canonical entry".
Null        the same 40-aa windowed shuffle family used elsewhere (seeded per
            query), because that is the null the sequence channel was calibrated
            against; a different null here would make the two libraries
            incomparable.
Labels      every record is written to the working FASTA as ``>ACC|GENE`` with
            the accession taken from the ``|`` token that matches the UniProt
            accession pattern. Swiss-Prot cache files use ``>ACC|GENE`` while
            UniProt stream files use ``>tr|ACC|NAME_HUMAN``, and taking field 1
            blindly labels reviewed proteins with their *gene* name -- which
            silently destroys every gene-level comparison downstream.

Run from the repository root:

    PYTHONPATH=src env/venv/bin/python src/analysis/61_trembl_expansion.py
"""

from __future__ import annotations

import glob
import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
REPORTS = ROOT / "reports"
LOGS = ROOT / "logs"
WORK = ROOT / "work" / "trembl"

MM = "env/bin/mmseqs"
FORMAT = "query,target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend"
COLS = FORMAT.split(",")
GA_BLOCK = (90, 330)
WINDOW = 40
B = 200
SEED0 = 9000
SEARCH = ("-s 7.5 --max-seqs 20000 -e {e} --comp-bias-corr 1 "
          "--threads 14 -a")
# The reporting threshold has to be run twice. E-values scale with the size of
# the target database, so a 10x larger library pushes weak alignments above any
# fixed -e cutoff. Comparing a 20k library at -e 10000 with a 200k library at
# -e 10000 therefore confounds "what the library contains" with "what the
# threshold still admits"; the loose pass removes that confound.
CUTOFFS = (("e1e4", "10000"), ("e1e6", "1000000"))
PUBLISHED = {"CRYAB": "P02511", "ANO2": "Q9NQ90", "MEF2B": "Q02080",
             "GlialCAM": "Q14CZ8"}
ACC_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
                    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
                    r"|[A-Z0-9]{10})$")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str], log_path: Path) -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("$ " + " ".join(map(str, cmd)) + "\n")
        log.flush()
        proc = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed with {proc.returncode}; see {log_path}")


def _split_accession(header: str) -> tuple[str, str]:
    """Return (accession, gene) from a UniProt header, whatever its dialect."""
    body = header[1:]
    acc = ""
    for token in body.split("|"):
        if ACC_RE.match(token):
            acc = token
            break
    if not acc:
        acc = body.split()[0]
    match = re.search(r"GN=(\S+)", body)
    if match:
        gene = match.group(1)
    else:
        rest = [t for t in body.split("|") if not ACC_RE.match(t)]
        gene = rest[0].split()[0] if rest else acc
    return acc, gene


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    acc = gene = ""
    parts: list[str] = []

    def flush() -> None:
        if parts:
            records.append((acc, gene, "".join(parts)))

    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                parts = []
                acc, gene = _split_accession(line.rstrip("\n"))
            else:
                parts.append(line.strip())
    flush()
    return records


def fetch_trembl() -> tuple[Path, str]:
    out = CACHE / "human_trembl.fasta"
    if not out.exists() or out.stat().st_size == 0:
        run(["curl", "-fsSL", "--retry", "3", "-G",
             "https://rest.uniprot.org/uniprotkb/stream",
             "--data-urlencode", "query=organism_id:9606 AND reviewed:false",
             "--data-urlencode", "format=fasta",
             "-o", str(out)], LOGS / "61b_trembl_download.log")
    headers = subprocess.run(
        ["curl", "-s", "-G", "https://rest.uniprot.org/uniprotkb/search",
         "--data-urlencode", "query=organism_id:9606 AND reviewed:true",
         "--data-urlencode", "format=list", "--data-urlencode", "size=1",
         "-D", "-", "-o", "/dev/null"],
        cwd=ROOT, capture_output=True, text=True).stdout
    release = ""
    for line in headers.splitlines():
        if line.lower().startswith("x-uniprot-release:"):
            release = line.split(":", 1)[1].strip()
    return out, release


def dedupe(reviewed: list[tuple[str, str, str]],
           trembl: list[tuple[str, str, str]]) -> tuple[list, int]:
    seen = {seq for _, _, seq in reviewed}
    kept, dropped = [], 0
    for acc, gene, seq in trembl:
        if seq in seen:
            dropped += 1
            continue
        seen.add(seq)
        kept.append((acc, gene, seq))
    return kept, dropped


def write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    with path.open("w") as fh:
        for acc, gene, seq in records:
            # a query record has no gene, and ">name|" would then stop matching
            # the label the search engine reports back
            fh.write(f">{acc}|{gene}\n" if gene else f">{acc}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def windowed_shuffle(seq: str, width: int, rng: random.Random) -> str:
    out = list(seq)
    for start in range(0, len(seq), width):
        block = out[start:start + width]
        rng.shuffle(block)
        out[start:start + width] = block
    return "".join(out)


def _clean(prefix: Path) -> None:
    """Remove a mmseqs prefix: it may be a file, a directory, or a family."""
    for suffix in ("", ".dbtype", ".index", ".lookup", ".source", ".mapping",
                   "_h", "_h.dbtype", "_h.index", "_h_h"):
        for path in glob.glob(str(prefix) + suffix):
            p = Path(path)
            if p.is_dir():
                shutil.rmtree(p)
            else:
                p.unlink()


def build_target(fasta: Path, prefix: Path, log: Path) -> Path:
    _clean(prefix)
    run([MM, "createdb", str(fasta), str(prefix), "--shuffle", "0"], log)
    return prefix


def search(queries: Path, target: Path, out_m8: Path, log: Path,
           e_cutoff: str) -> None:
    qdb, res, tmp = WORK / "qdb", WORK / "res", WORK / "tmp"
    for p in (qdb, res, tmp):
        _clean(p)
    run([MM, "createdb", str(queries), str(qdb)], log)
    run([MM, "search", str(qdb), str(target), str(res), str(tmp)]
        + SEARCH.format(e=e_cutoff).split(), log)
    run([MM, "convertalis", str(qdb), str(target), str(res), str(out_m8),
         "--format-output", FORMAT], log)


def load_m8(path: Path, gene_of: dict[str, str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, names=COLS, dtype=str)
    for col in ("evalue", "bits", "alnlen", "pident", "qstart", "qend"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    acc = frame["target"].map(lambda s: _split_accession(">" + s)[0])
    frame["acc"] = acc
    frame["gene"] = acc.map(lambda a: gene_of.get(a, a))
    return frame


def cell(frame: pd.DataFrame, real_name: str) -> dict:
    real = frame[frame["query"] == real_name]
    null = frame[frame["query"] != real_name]
    best_null = null.groupby("query")["evalue"].min()
    best_real = float(real["evalue"].min()) if len(real) else None
    overlap = int(((real["qend"] >= GA_BLOCK[0]) &
                   (real["qstart"] <= GA_BLOCK[1])).sum()) if len(real) else 0
    n_ge = int((best_null <= best_real).sum()) if best_real is not None else 0
    return {
        "n_align": int(len(real)),
        "n_acc": int(real["acc"].nunique()),
        "n_gene": int(real["gene"].nunique()),
        "null_n_align_mean": round(float(len(null) / max(len(best_null), 1)), 2),
        "best_e": best_real,
        "n_below_1": int((real["evalue"] < 1).sum()) if len(real) else 0,
        "region_overlap_ga": overlap,
        "null_best_median": round(float(best_null.median()), 3) if len(best_null) else None,
        "null_best_iqr": [round(float(best_null.quantile(.25)), 3),
                          round(float(best_null.quantile(.75)), 3)] if len(best_null) else None,
        "n_null": int(len(best_null)),
        "n_null_ge_real": n_ge,
        "p": round((1 + n_ge) / (1 + len(best_null)), 4) if len(best_null) else None,
    }


def main() -> int:
    started = time.time()
    for d in (CACHE, WORK, REPORTS, LOGS):
        d.mkdir(parents=True, exist_ok=True)
    log = LOGS / "61b_trembl_expansion.log"
    log.write_text("")

    trembl_path, release = fetch_trembl()
    reviewed_path = CACHE / "human_swissprot.fasta"
    reviewed = read_fasta(reviewed_path)
    trembl = read_fasta(trembl_path)
    kept, dropped = dedupe(reviewed, trembl)
    gene_of = {acc: gene for acc, gene, _ in reviewed + kept}
    print(f"reviewed {len(reviewed)} | trembl {len(trembl)} "
          f"| dropped as exact duplicates {dropped} | kept {len(kept)}")

    combined_path = CACHE / "human_combined.fasta"
    write_fasta(combined_path, reviewed + kept)

    ebna1 = read_fasta(CACHE / "P03211.fasta")[0][2]
    rng = random.Random(SEED0)
    queries = [("real_EBNA1", ebna1)]
    for i in range(B):
        queries.append((f"w40_{i:04d}", windowed_shuffle(ebna1, WINDOW, rng)))
    queries_path = WORK / "queries.fasta"
    write_fasta(queries_path, [(name, "", seq) for name, seq in queries])

    cells, frames = {}, {}
    for tag, fasta in (("reviewed", reviewed_path), ("combined", combined_path)):
        tdb = build_target(fasta, WORK / f"tdb_{tag}", log)
        for cut_name, cut in CUTOFFS:
            key = f"{tag}_{cut_name}"
            m8 = WORK / f"{key}.m8"
            search(queries_path, tdb, m8, log, cut)
            frames[key] = load_m8(m8, gene_of)
            cells[key] = cell(frames[key], "real_EBNA1")
            print(f"[{key}] {cells[key]['n_align']} alignments, "
                  f"{cells[key]['n_gene']} genes, best E {cells[key]['best_e']}, "
                  f"E<1 {cells[key]['n_below_1']}, "
                  f"overlap 90-330 {cells[key]['region_overlap_ga']}, "
                  f"null {cells[key]['null_best_median']} p={cells[key]['p']}")

    def real_rows(key: str) -> pd.DataFrame:
        return frames[key][frames[key]["query"] == "real_EBNA1"]

    real_combined = real_rows("combined_e1e6")
    published = {}
    for name, acc in PUBLISHED.items():
        hits = real_combined[real_combined["acc"] == acc]
        published[name] = {
            "acc": acc,
            "in_library": acc in gene_of,
            "alignments_in_combined": int(len(hits)),
            "best_e": float(hits["evalue"].min()) if len(hits) else None,
        }

    genes = {}
    for cut_name, _ in CUTOFFS:
        g_rev = set(real_rows(f"reviewed_{cut_name}")["gene"])
        g_com = set(real_rows(f"combined_{cut_name}")["gene"])
        genes[cut_name] = {
            "hit_in_reviewed": len(g_rev),
            "hit_in_combined": len(g_com),
            "new_in_combined": len(g_com - g_rev),
            "lost_from_reviewed": len(g_rev - g_com),
            "examples_new": sorted(g_com - g_rev)[:20],
        }

    out = {
        "release": release,
        "db": {
            "reviewed_n": len(reviewed),
            "trembl_n": len(trembl),
            "trembl_exact_duplicates_dropped": dropped,
            "trembl_kept": len(kept),
            "combined_n": len(reviewed) + len(kept),
            "sha256": {
                "reviewed": sha256(reviewed_path),
                "trembl": sha256(trembl_path),
                "combined": sha256(combined_path),
            },
        },
        "parameters": {
            "mmseqs": SEARCH, "format": FORMAT, "window": WINDOW, "b": B,
            "seed": SEED0, "ga_block": GA_BLOCK,
            "e_cutoffs": dict(CUTOFFS),
        },
        "cells": cells,
        "published": published,
        "genes": genes,
        "note": ("exact-duplicate sequences between TrEMBL and Swiss-Prot are "
                 "dropped; near-duplicates are not clustered, so the enlarged "
                 "library is a union of canonical entries rather than a "
                 "non-redundant set. Isoform sequences are absent. Because "
                 "TrEMBL holds fragments of already-reviewed proteins, the "
                 "accession count overstates what was added; the gene count is "
                 "the honest one."),
        "seconds": round(time.time() - started, 1),
    }
    (REPORTS / "61b_trembl_expansion.json").write_text(json.dumps(out, indent=2))
    rows = [{"library": tag, **vals} for tag, vals in cells.items()]
    pd.DataFrame(rows).to_csv(REPORTS / "61b_trembl_expansion.tsv", sep="\t",
                              index=False)
    print(f"wrote reports/61b_trembl_expansion.json ({out['seconds']}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
