"""s50 -- the sequence-alignment channel, kept for continuity, not for evidence.

This runs the local C port of the reference implementation so the historical
numbers stay reproducible (seeded + exhaustive hit counts). It is NOT the basis
of any candidate call: a whole-proteome alignment count is dominated by
composition, which is exactly what stages s40/s60 measure.
"""
from __future__ import annotations

import subprocess

import pandas as pd

from ..core import ROOT, write_json

BINARY = ROOT / "env" / "bin" / "mimic_scan"


def _run(ctx, mode: str, out_name: str) -> pd.DataFrame:
    cmd = [
        str(BINARY),
        "--query", str(ctx.cfg.path_of("data.query")),
        "--db", str(ctx.cfg.path_of("data.proteome")),
        "--mode", mode,
        "--out", str(ctx.out(out_name)),
        "--kmer", str(ctx.cfg.get("alignment.kmer", 4)),
        "--xdrop", str(ctx.cfg.get("alignment.xdrop", 6)),
        "--slack", str(ctx.cfg.get("alignment.slack", 20)),
        "--gap-open", str(ctx.cfg.get("alignment.gap_open", 11)),
        "--min-len", str(ctx.cfg.get("alignment.min_len", 6)),
        "--min-score", str(ctx.cfg.get("alignment.min_score", 20)),
        "--n-iter", str(ctx.cfg.get("alignment.exhaustive_iter", 4)),
        "--threads", str(ctx.cfg.get("alignment.threads", 8)),
    ]
    ctx.log(f"    {' '.join(cmd[1:])}")
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        raise RuntimeError(f"mimic_scan {mode} failed rc={proc.returncode}: "
                           f"{proc.stderr.strip()[:400]}")
    for line in proc.stderr.strip().splitlines():
        ctx.log(f"      {line}")
    return pd.read_csv(ctx.out(out_name), sep="\t", dtype={"acc": str, "gene": str})


def run(ctx):
    if not ctx.cfg.get("alignment.enabled", True):
        ctx.log("    alignment channel disabled in config; skipped")
        return dict(outputs=[], facts=dict(skipped=True))
    if not BINARY.exists():
        ctx.log(f"    {BINARY} not built; skipped (run: make -C {ROOT.name})")
        return dict(outputs=[], facts=dict(skipped="binary missing"))

    block = ctx.cfg.get("composition.gar_block", [90, 330])
    facts = {}
    outputs = []
    for mode, name in (("seeded", "50_seeded.tsv"), ("exhaustive", "50_exhaustive.tsv")):
        frame = _run(ctx, mode, name)
        overlap = int(((frame["qend"] >= block[0]) & (frame["qstart"] <= block[1])).sum())
        facts[mode] = dict(n_hits=int(len(frame)), n_overlapping_block=overlap,
                           top_score=float(frame["gscore"].max()) if len(frame) else None)
        outputs.append(f"reports/{name}")
        ctx.log(f"    {mode}: {len(frame):,} hits, {overlap} overlapping EBNA1 "
                f"{block[0]}-{block[1]}")
    write_json(ctx.out("50_alignment_summary.json"), facts)
    outputs.append("reports/50_alignment_summary.json")
    return dict(outputs=outputs, facts=facts)
