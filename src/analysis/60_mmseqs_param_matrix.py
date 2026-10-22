#!/usr/bin/env python3
"""A1 -- MMseqs2 low-complexity / composition-correction parameter matrix.

The unbiased channel reports that 0 of 363 alignments touch EBNA1 residues
90-330, a 240-residue stretch that is 99.6% Gly/Ala. Three mechanisms could
produce that and none of them was separated:

  (i)   MMseqs2's default low-complexity masking (tantan, --mask 1);
  (ii)  local composition correction (--comp-bias-corr 1);
  (iii) the intrinsic difficulty of scoring a near-homopolymer at all.

This script runs the same search in a grid where (i) and (ii) are switched
independently, and separately masks the *query* instead of the *database* --
those answer different questions: target-side masking asks whether human
low-complexity proteins are being excluded as sources of hits, query-side
masking asks whether EBNA1's own block can generate hits at all.

Every cell gets its own null. Raw hit counts are NOT comparable across cells
because masking changes the size of the effective search space, so the number
that matters per cell is where the real sequence sits inside its own null
distribution.

Null family: shuffle EBNA1 inside 40-aa windows (preserves local composition,
destroys order), B=200, same seed as 07_batch_null.py. For query-masked cells
the same span is masked in every shuffle, so the comparison stays like-for-like.

Usage: PYTHONPATH=src env/venv/bin/python src/analysis/60_mmseqs_param_matrix.py
"""
from __future__ import annotations

import json
import os
import random
import shutil
import subprocess

import numpy as np
import pandas as pd

MM = "env/bin/mmseqs"
FORMAT = "query,target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend"
GA_BLOCK = (90, 330)
WORK = "work/param_matrix"
OUT_TSV = "reports/60b_mmseqs_param_matrix.tsv"
OUT_JSON = "reports/60b_mmseqs_param_matrix.json"
B = 200
WINDOW = 40
SEED0 = 7000


def read_fasta_seq(path: str) -> str:
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def windowed_shuffle(seq: str, w: int, rng: random.Random) -> str:
    if w >= len(seq):
        out = list(seq)
        rng.shuffle(out)
        return "".join(out)
    out = list(seq)
    for s in range(0, len(seq), w):
        blk = list(seq[s:s + w])
        rng.shuffle(blk)
        out[s:s + w] = blk
    return "".join(out)


def mask_span(seq: str, span: tuple[int, int]) -> str:
    """Replace a 1-based closed interval with X, the way a user mask would."""
    lo, hi = span[0] - 1, span[1]
    return seq[:lo] + "X" * (hi - lo) + seq[hi:]


def run(cmd: list[str], log) -> None:
    log.write(" ".join(cmd) + "\n")
    log.flush()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[1]} failed:\n{proc.stdout}\n{proc.stderr}")


def build_target(mask: int) -> str:
    prefix = os.path.join(WORK, f"tdb_mask{mask}")
    if os.path.exists(prefix + ".dbtype"):
        return prefix
    run([MM, "createdb", "data/cache/human_swissprot.fasta", prefix,
         "--mask", str(mask), "--shuffle", "0"], log=open(os.path.join(WORK, "createdb.log"), "a"))
    return prefix


def search_cell(tag: str, queries: list[tuple[str, str]], tdb: str,
                mask: int, cbc: int, log) -> pd.DataFrame:
    """One batched search: every query in one FASTA, one target index."""
    cell = os.path.join(WORK, tag)
    shutil.rmtree(cell, ignore_errors=True)
    os.makedirs(cell, exist_ok=True)
    qfa = os.path.join(cell, "queries.fa")
    with open(qfa, "w") as fh:
        for name, seq in queries:
            fh.write(f">{name}\n{seq}\n")
    qdb = os.path.join(cell, "qdb")
    res = os.path.join(cell, "res")
    m8 = os.path.join(cell, "out.m8")
    run([MM, "createdb", qfa, qdb, "--shuffle", "0"], log=log)
    run([MM, "search", qdb, tdb, res, os.path.join(cell, "tmp"),
         "-s", "7.5", "--max-seqs", "20000", "-e", "10000",
         "--comp-bias-corr", str(cbc), "--mask", str(mask),
         "-a", "--threads", "14"], log=log)
    run([MM, "convertalis", qdb, tdb, res, m8, "--format-output", FORMAT,
         "--threads", "14"], log=log)
    frame = pd.read_csv(m8, sep="\t", header=None,
                        names=FORMAT.split(","), dtype={"query": str})
    shutil.rmtree(cell, ignore_errors=True)
    return frame


def cell_stats(frame: pd.DataFrame, real_name: str) -> dict:
    real = frame[frame["query"] == real_name]
    ev = real["evalue"].astype(float)
    ga = ((real["qstart"].astype(int) <= GA_BLOCK[1]) &
          (real["qend"].astype(int) >= GA_BLOCK[0])).sum()
    per_query = frame.groupby("query").size()
    null = per_query.drop(index=real_name, errors="ignore")
    best_per_query = frame.groupby("query")["evalue"].min()
    null_best = best_per_query.drop(index=real_name, errors="ignore")
    n_real = int(len(real))
    n_ge = int((null >= n_real).sum())
    real_best = float(ev.min()) if n_real else float("inf")
    return dict(
        n_align=n_real,
        best_evalue=round(real_best, 4) if np.isfinite(real_best) else None,
        best_target=(real.loc[ev.idxmin(), "target"] if n_real else ""),
        n_e_lt_1=int((ev < 1).sum()),
        n_e_lt_10=int((ev < 10).sum()),
        n_e_lt_100=int((ev < 100).sum()),
        n_overlapping_ga_block=int(ga),
        null_n_align_mean=round(float(null.mean()), 2),
        null_n_align_sd=round(float(null.std(ddof=1)), 2),
        null_n_align_max=int(null.max()) if len(null) else 0,
        null_best_evalue_median=round(float(null_best.median()), 2) if len(null_best) else None,
        p_n_align=round((1 + n_ge) / (1 + len(null)), 4),
        n_null=len(null),
    )


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    os.makedirs("reports", exist_ok=True)
    log = open("logs/60b_mmseqs_param_matrix.log", "w")
    seq = read_fasta_seq("data/cache/P03211.fasta")
    seq_masked = mask_span(seq, GA_BLOCK)

    rng = random.Random(SEED0)
    shuffles = [windowed_shuffle(seq, WINDOW, rng) for _ in range(B)]
    shuffles_masked = [mask_span(s, GA_BLOCK) for s in shuffles]

    cells = [
        # tag,                 query_masked, target mask, comp-bias-corr
        ("mask1_cbc1", False, 1, 1),
        ("mask1_cbc0", False, 1, 0),
        ("mask0_cbc1", False, 0, 1),
        ("mask0_cbc0", False, 0, 0),
        ("qmask_cbc1", True, 0, 1),
        ("qmask_cbc0", True, 0, 0),
    ]
    rows = []
    for tag, qmasked, tmask, cbc in cells:
        tdb = build_target(tmask)
        real_seq = seq_masked if qmasked else seq
        null_seqs = shuffles_masked if qmasked else shuffles
        queries = [("real", real_seq)] + [
            (f"n{i:04d}", s) for i, s in enumerate(null_seqs)
        ]
        print(f"[{tag}] query_masked={qmasked} target_mask={tmask} cbc={cbc}",
              flush=True)
        frame = search_cell(tag, queries, tdb, tmask, cbc, log)
        stats = cell_stats(frame, "real")
        stats["cell"] = tag
        stats["query_masked_90_330"] = qmasked
        stats["target_mask"] = tmask
        stats["comp_bias_corr"] = cbc
        rows.append(stats)
        print(f"    real n_align={stats['n_align']} bestE={stats['best_evalue']} "
              f"ga_overlap={stats['n_overlapping_ga_block']} "
              f"null={stats['null_n_align_mean']}+-{stats['null_n_align_sd']} "
              f"p={stats['p_n_align']}", flush=True)

    table = pd.DataFrame(rows)[[
        "cell", "query_masked_90_330", "target_mask", "comp_bias_corr",
        "n_align", "best_evalue", "best_target", "n_e_lt_1", "n_e_lt_10",
        "n_e_lt_100", "n_overlapping_ga_block", "null_n_align_mean",
        "null_n_align_sd", "null_n_align_max", "null_best_evalue_median",
        "p_n_align", "n_null",
    ]]
    table.to_csv(OUT_TSV, sep="\t", index=False)
    summary = dict(
        what="MMseqs2 low-complexity and composition-correction matrix",
        version=subprocess.run([MM, "version"], capture_output=True, text=True).stdout.strip(),
        fixed=dict(sensitivity=7.5, max_seqs=20000, e_cutoff=10000, gap_open=11,
                   database="data/cache/human_swissprot.fasta (reviewed, 20431)"),
        null=dict(family=f"shuffle inside {WINDOW}-aa windows", b=B, seed0=SEED0,
                  note=("each cell is calibrated against its own null; raw hit "
                        "counts are not comparable across cells because masking "
                        "changes the effective search space")),
        cells=rows,
        note=("query-side masking replaces EBNA1 90-330 with X in the real query "
              "and in every shuffle; target-side masking is baked into the "
              "database by createdb --mask and repeated in search --mask"),
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)
    log.close()
    print(f"\nwrote {OUT_TSV}")


if __name__ == "__main__":
    main()
