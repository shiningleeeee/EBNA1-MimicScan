"""Shared machinery for the unbiased channel: whole EBNA1 vs the human proteome.

Nothing in this module reads `seeds:` or `published:`, and that is the point.
The seed channel answers "what does a match to a validated peptide look like".
This channel answers the complementary question -- given all 641 residues of
EBNA1, where does the human proteome locally resemble it? -- and any influence
from a known epitope would destroy it. The published claims are audited
afterwards (see `s71_candidates`), never used as input.

Two properties the callers rely on:

  * every external call is retried and logged, because a silent network failure
    would drop a null-model sample and change a p-value between two runs that
    are supposed to be identical;
  * every table this module writes is sorted on a total order before it is
    written to disk, so mmseqs's thread-dependent output order cannot leak into
    a file hash.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from .core import ROOT

# ---------------------------------------------------------------- vocabularies
# Antibody reachability and CNS expression are read off UniProt text, so they
# are annotations about the protein, not about the segment. They are reported
# as columns and never combined into a score (see s68_layered).
SURFACE_TOKENS = (
    "cell membrane", "secreted", "cell surface", "extracellular",
    "gpi-anchor", "basement membrane",
)
CNS_TOKENS = (
    "brain", "cerebellum", "cerebral cortex", "hippocampus", "spinal cord",
    "neuron", "neuro", "synap", "myelin", "glia", "axon", "dendrit",
)

M8_COLUMNS = ["target", "evalue", "bits", "qstart", "qend", "tstart", "tend",
              "alnlen", "pident", "qaln", "taln"]


# ------------------------------------------------------------------------ fasta
def load_fasta(path) -> dict[str, str]:
    """accession -> sequence. Accession comes from the UniProt `sp|ACC|NAME` head."""
    out: dict[str, str] = {}
    acc, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    out[acc] = "".join(buf)
                head = line[1:]
                acc = head.split("|")[1] if "|" in head else head.split()[0]
                buf = []
            else:
                buf.append(line)
    if acc is not None:
        out[acc] = "".join(buf)
    return out


def first_record(path) -> tuple[str, str]:
    acc, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    break
                head = line[1:]
                acc = head.split("|")[1] if "|" in head else head.split()[0]
            else:
                seq.append(line)
    return acc, "".join(seq)


def write_fasta(path, records) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i: i + 60] + "\n")


# ------------------------------------------------------------------ composition
def entropy(seq: str) -> tuple[float, float]:
    """Shannon entropy in bits and the largest single-residue fraction."""
    if not seq:
        return 0.0, 1.0
    counts = Counter(seq)
    n = len(seq)
    h = -sum((v / n) * math.log2(v / n) for v in counts.values())
    return h, max(counts.values()) / n


def flag(text, tokens) -> bool:
    low = str(text).lower()
    return any(tok in low for tok in tokens)


def region_of(qstart: int, qend: int, ga_block=(90, 330)) -> str:
    """Where on EBNA1 an alignment sits -- descriptive, never a filter.

    The rule is a cascade with two different tests inside it, and both halves
    matter: the Gly/Ala block is decided by *overlap* (an alignment that touches
    it at all is called a Gly/Ala-block alignment, even when most of its length
    lies further C-terminal), while the SPR-rich and folded-core bands are
    decided by *start position*. That mix is what the region counts in the
    reports were computed with, so it is written down here rather than left to
    be inferred from the numbers.
    """
    lo, hi = ga_block
    if qend < lo:
        return "N-terminal"
    if qstart <= hi and qend >= lo:
        return "Gly/Ala block (%d-%d)" % (lo, hi)
    if 370 <= qstart <= 445:
        return "SPR-rich (370-445)"
    if 476 <= qstart <= 585:
        return "folded core (476-585)"
    return "C-terminal"


def overlaps(a1: int, a2: int, b1: int, b2: int) -> bool:
    return a1 <= b2 and b1 <= a2


# --------------------------------------------------------------- subprocess
def run_logged(cmd, log_path) -> None:
    """Run a command, appending the exact argv and its output to the stage log."""
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as fh:
        fh.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
        fh.flush()
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh,
                       stderr=subprocess.STDOUT)


def _mmseqs_binary() -> Path:
    return ROOT / "env" / "bin" / "mmseqs"


def pair_null_binary() -> Path:
    return ROOT / "env" / "bin" / "pair_null"


# ------------------------------------------------------------------ mmseqs
def _clean_prefix(prefix: Path) -> None:
    """Remove a mmseqs database, which is a *set* of files sharing a prefix.

    `createdb out` writes out, out.dbtype, out.index, ... so `out` itself is a
    regular file. Deleting it with rmtree raises NotADirectoryError, which is
    what a second run of the pipeline used to do.
    """
    if prefix.is_dir():
        shutil.rmtree(prefix)
    elif prefix.exists():
        prefix.unlink()
    for sibling in prefix.parent.glob(prefix.name + ".*"):
        if sibling.is_dir():
            shutil.rmtree(sibling)
        else:
            sibling.unlink()


def mmseqs_search(cfg, query_fa, out_m8, work_dir, log_path, name="q",
                  require_query_column=False) -> Path:
    """createdb + search + convertalis, with a reused target database.

    `-a` is mandatory: without it convertalis refuses with
    "Backtrace cigar is missing in the alignment result" and no qaln/taln exist.
    """
    mm = _mmseqs_binary()
    if not mm.exists():
        raise RuntimeError(f"mmseqs not built at {mm}; run `make` first")
    work = ROOT / work_dir
    work.mkdir(parents=True, exist_ok=True)
    tdb = work / "tdb"
    if not tdb.exists():
        run_logged([mm, "createdb", cfg.path_of("data.proteome"), tdb], log_path)
    qdb = work / f"{name}db"
    _clean_prefix(qdb)
    run_logged([mm, "createdb", query_fa, qdb], log_path)
    res = work / f"{name}res"
    _clean_prefix(res)
    tmp = work / f"{name}tmp"
    _clean_prefix(tmp)
    run_logged([
        mm, "search", qdb, tdb, res, tmp,
        "-s", cfg.get("mmseqs.sensitivity", 7.5),
        "--max-seqs", cfg.get("mmseqs.max_seqs", 20000),
        "-e", cfg.get("mmseqs.e_cutoff", 10000),
        "--comp-bias-corr", cfg.get("mmseqs.comp_bias_corr", 1),
        "-a", "--threads", cfg.get("mining.threads", 14),
    ], log_path)
    fmt = "query," + ",".join(M8_COLUMNS) if require_query_column else ",".join(M8_COLUMNS)
    out_m8 = ROOT / out_m8 if not Path(out_m8).is_absolute() else Path(out_m8)
    out_m8.parent.mkdir(parents=True, exist_ok=True)
    run_logged([mm, "convertalis", qdb, tdb, res, out_m8,
                "--format-output", fmt], log_path)
    return out_m8


def read_m8(path, require_query_column=False):
    import pandas as pd

    names = (["query"] + M8_COLUMNS) if require_query_column else M8_COLUMNS
    # Width is checked first: a search run without `-a` drops two columns, and
    # reading that file with the full column list would silently shift qaln/taln
    # into evalue/bits instead of failing.
    width = pd.read_csv(path, sep="\t", header=None, nrows=1).shape[1]
    if width != len(names):
        raise ValueError(
            f"{path}: expected {len(names)} columns, found {width}"
            + (" (missing -a in mmseqs search?)" if width < len(names) else ""))
    frame = pd.read_csv(path, sep="\t", names=names, dtype={"target": str,
                                                             "query": str})
    # A total order, so a multi-threaded search cannot change the file hash.
    order = [c for c in ("query", "evalue", "target", "qstart", "tstart") if c in frame.columns]
    return frame.sort_values(order, kind="mergesort").reset_index(drop=True)


def write_m8(frame, path) -> None:
    frame.to_csv(path, sep="\t", header=False, index=False)


# ------------------------------------------------------------------ pair null
def pair_null(cfg, query_fa, cands_tsv, out_tsv, b=None, seed=None, flank=None,
              log_path=None, quiet=True):
    """Score every candidate alignment against its own two-sided composition null."""
    import pandas as pd

    binary = pair_null_binary()
    if not binary.exists():
        raise RuntimeError(f"pair_null not built at {binary}; run `make` first")
    b = cfg.get("mining.pair_null.b", 10000) if b is None else b
    seed = cfg.get("mining.pair_null.seed", 42) if seed is None else seed
    flank = cfg.get("mining.pair_null.flank", 5) if flank is None else flank
    log_path = log_path or (cfg.logs / "pair_null.log")
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "--query", query_fa, "--db", cfg.path_of("data.proteome"),
           "--cands", cands_tsv, "--b", b, "--seed", seed, "--flank", flank,
           "--out", out_tsv]
    if quiet:
        with open(log_path, "a") as fh:
            fh.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
            fh.flush()
            subprocess.run([str(c) for c in cmd], check=True,
                           stdout=fh, stderr=subprocess.STDOUT)
    else:
        run_logged(cmd, log_path)
    frame = pd.read_csv(out_tsv, sep="\t")
    order = [c for c in ("target", "qstart", "tstart") if c in frame.columns]
    return frame.sort_values(order, kind="mergesort").reset_index(drop=True)


# ------------------------------------------------------ maturation annotation
def ensure_processing_tsv(cfg, log_path) -> Path:
    """Download the UniProt maturation features once and cache them.

    Signal peptide / propeptide / transit peptide are proteolytic removal
    annotations: a match inside one of those regions does not exist on the
    mature protein an antibody would meet.
    """
    path = cfg.path_of("mining.maturation.cache")
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    params = {
        "query": cfg.get("mining.maturation.query", "organism_id:9606 AND reviewed:true"),
        "format": "tsv",
        "fields": cfg.get("mining.maturation.fields", ""),
    }
    url = cfg.get("mining.maturation.url",
                  "https://rest.uniprot.org/uniprotkb/stream")
    query_url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(query_url, timeout=300) as fh:
                blob = fh.read()
            path.write_bytes(blob)
            with open(log_path, "a") as log:
                log.write(f"\n# downloaded {query_url} -> {len(blob)} bytes\n")
            return path
        except Exception as exc:  # noqa: BLE001 - retried, then surfaced
            last = exc
            time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not download maturation features: {last}")


def load_annotation(cfg) -> dict:
    """accession -> {cns, surface} read off the bulk UniProt annotation TSV."""
    import pandas as pd

    path = cfg.path_of("data.uniprot_annot")
    if not path.exists():
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    key = "Entry" if "Entry" in frame.columns else frame.columns[0]
    out = {}
    for row in frame.to_dict("records"):
        blob = " ".join(str(v).lower() for v in row.values())
        out[str(row[key])] = {
            "cns": flag(blob, CNS_TOKENS),
            "surface": flag(blob, SURFACE_TOKENS),
        }
    return out
