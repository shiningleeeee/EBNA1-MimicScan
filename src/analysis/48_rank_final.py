#!/usr/bin/env python3
"""48_rank_final.py — build the hypothesis list AND calibrate the list itself.

Design contract
---------------
1. The ranking is a HYPOTHESIS-GENERATION ordering. It is not a detection
   claim. Sequence evidence here cannot establish cross-reactivity; what a
   ranking can do is put the entries most worth an experiment at the top and
   say exactly why each one is there.

2. Every axis is declared before looking at the result, and none of them is
   fitted to the three published mimics. The published claims are audited
   afterwards (§ benchmark), never used as input.

3. The composite is calibrated AS A WHOLE, not axis by axis: the identical
   procedure is applied to composition-preserving shuffles of EBNA1. If a
   shuffle reaches the same number of high-composite entries, the list is a
   property of the proteome, not of EBNA1.

Axes (1 point each, documented):
    order_signal   per-pair composition null margin_z >= 3.0
                   (the observed alignment beats 99.87% of its own
                    two-sided residue-order null)
    length_ok      aligned length >= 20
    identity_ok    percent identity >= 40
    target_normal  target window entropy >= 2.8 and max residue <= 30%
    query_normal   EBNA1 window entropy >= 2.8 and max residue <= 30%
    cns            brain / CNS expression evidence
    surface        plasma membrane / secreted / extracellular
    outside_gablock  EBNA1 window does not overlap 90-330

Usage:
  env/venv/bin/python src/48_rank_final.py --b 10000 --min-score 5
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
WORK = ROOT / "work" / "pair_null_batch"
PROTEOME = ROOT / "data" / "cache" / "human_swissprot.fasta"
QUERY = ROOT / "data" / "cache" / "P03211.fasta"
ANNOT = ROOT / "data" / "cache" / "uniprot_human_annot.tsv"
PAIR_NULL = ROOT / "env" / "bin" / "pair_null"
LOGS = ROOT / "logs"

CNS_TOKENS = (
    "brain", "cerebellum", "cerebral cortex", "hippocampus", "spinal cord",
    "neuron", "neuro", "synap", "myelin", "glia", "axon", "dendrit",
)
SURFACE_TOKENS = (
    "cell membrane", "secreted", "cell surface", "extracellular",
    "gpi-anchor", "basement membrane",
)
GA_BLOCK = (90, 330)
Z_MIN = 3.0
ENT_MIN = 2.8
MAXFRAC_MAX = 0.30


def load_fasta(path: Path) -> dict[str, str]:
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


def entropy(seq: str) -> tuple[float, float]:
    if not seq:
        return 0.0, 1.0
    c = Counter(seq)
    n = len(seq)
    h = -sum((v / n) * math.log2(v / n) for v in c.values())
    return h, max(c.values()) / n


def overlap(a1: int, a2: int, b1: int, b2: int) -> bool:
    return a1 <= b2 and b1 <= a2


def compose(row, qseq: str, tseq: str) -> dict:
    qw = qseq[row["qstart"] - 1 : row["qend"]]
    tw = tseq[row["tstart"] - 1 : row["tend"]]
    qh, qmf = entropy(qw)
    th, tmf = entropy(tw)
    return dict(q_win=qw, t_win=tw, q_ent=qh, q_maxfrac=qmf,
                t_ent=th, t_maxfrac=tmf)


def score_axes(row, ann: dict, extra: dict) -> dict:
    ident = float(row["pident"]) if row["pident"] != "-" else 0.0
    a = ann.get(str(row["target"]), {})
    axes = {
        "order_signal": float(row["margin_z"]) >= Z_MIN,
        "length_ok": int(row["alnlen"]) >= 20,
        "identity_ok": ident >= 40.0,
        "target_normal": extra["t_ent"] >= ENT_MIN and extra["t_maxfrac"] <= MAXFRAC_MAX,
        "query_normal": extra["q_ent"] >= ENT_MIN and extra["q_maxfrac"] <= MAXFRAC_MAX,
        "cns": bool(a.get("cns", False)),
        "surface": bool(a.get("surface", False)),
        "outside_gablock": not overlap(int(row["qstart"]), int(row["qend"]), *GA_BLOCK),
    }
    return axes


def run_pair_null(cands: Path, qfasta: Path, out: Path, b: int, seed: int) -> pd.DataFrame:
    subprocess.run(
        [str(PAIR_NULL), "--query", str(qfasta), "--db", str(PROTEOME),
         "--cands", str(cands), "--b", str(b), "--seed", str(seed),
         "--flank", "5", "--out", str(out)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return pd.read_csv(out, sep="\t")


def build(query_name: str, pn: pd.DataFrame, qseq: str,
          proteome: dict[str, str], annot: dict) -> pd.DataFrame:
    rows = []
    for _, r in pn.iterrows():
        tseq = proteome.get(r["target"])
        if not tseq:
            continue
        ex = compose(r, qseq, tseq)
        ax = score_axes(r, annot, ex)
        rows.append({
            "query": query_name,
            "acc": r["target"],
            "gene": r["gene"],
            "ebna1_start": int(r["qstart"]),
            "ebna1_end": int(r["qend"]),
            "t_start": int(r["tstart"]),
            "t_end": int(r["tend"]),
            "aln_len": int(r["alnlen"]),
            "identity": float(r["pident"]) if r["pident"] != "-" else 0.0,
            "obs_score": int(r["obs_score"]),
            "null_mean": float(r["null_mean"]),
            "null_sd": float(r["null_sd"]),
            "null_max": int(r["null_max"]),
            "margin_z": float(r["margin_z"]),
            "p_pair": float(r["p_pair"]),
            "q_win": ex["q_win"],
            "t_win": ex["t_win"],
            "q_ent": round(ex["q_ent"], 3),
            "q_maxfrac": round(ex["q_maxfrac"], 3),
            "t_ent": round(ex["t_ent"], 3),
            "t_maxfrac": round(ex["t_maxfrac"], 3),
            **{f"ax_{k}": v for k, v in ax.items()},
            "score": int(sum(ax.values())),
        })
    return pd.DataFrame(rows)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-score", type=int, default=5)
    args = ap.parse_args(argv)

    proteome = load_fasta(PROTEOME)
    qseq_real = list(load_fasta(QUERY).values())[0]
    print(f"[48] proteome {len(proteome)}, query {len(qseq_real)} aa")

    ann = {}
    if ANNOT.exists():
        df = pd.read_csv(ANNOT, sep="\t", dtype=str).fillna("")
        for _, r in df.iterrows():
            blob = " ".join(str(x).lower() for x in r.values)
            ann[r["Entry"]] = {
                "cns": any(t in blob for t in CNS_TOKENS),
                "surface": any(t in blob for t in SURFACE_TOKENS),
            }
    print(f"[48] annotation for {len(ann)} proteins")

    real = run_pair_null(REPORTS / "cands_for_null.tsv",
                         WORK / "q_real.fasta", WORK / "pn_real.tsv",
                         args.b, args.seed)
    real_df = build("real", real, qseq_real, proteome, ann)
    real_df = real_df.sort_values(["score", "margin_z", "aln_len"],
                                 ascending=[False, False, False])
    real_df.to_csv(REPORTS / "48_candidates_scored.tsv", sep="\t", index=False)
    print(f"[48] real: {len(real_df)} rows; score distribution "
          f"{dict(sorted(Counter(real_df['score']).items()))}")

    # ---- calibrate the composite against composition-preserving shuffles ----
    null_rows = []
    for cands in sorted(WORK.glob("cands_shuf*.tsv")):
        name = cands.stem.replace("cands_", "")
        qf = WORK / f"q_{name}.fasta"
        if not qf.exists():
            continue
        pn = run_pair_null(cands, qf, WORK / f"pn_{name}_z.tsv", args.b, args.seed)
        df = build(name, pn, list(load_fasta(qf).values())[0], proteome, ann)
        if df.empty:
            continue
        null_rows.append({
            "query": name, "n_rows": len(df),
            "max_score": int(df["score"].max()),
            **{f"n_ge_{k}": int((df["score"] >= k).sum()) for k in (4, 5, 6)},
        })
        print(f"[48] {name}: n={len(df)} max_score={null_rows[-1]['max_score']} "
              f"n>=5={null_rows[-1]['n_ge_5']}")
    null_df = pd.DataFrame(null_rows)
    if not null_df.empty:
        null_df.to_csv(REPORTS / "48_null_composite.csv", index=False)

    print("\n=== composite calibration (the list as a whole) ===")
    real_stats = {"max_score": int(real_df["score"].max()),
                  **{f"n_ge_{k}": int((real_df["score"] >= k).sum()) for k in (4, 5, 6)}}
    for k, v in real_stats.items():
        if null_df.empty:
            print(f"{k:12s}: real={v}  (no nulls)")
            continue
        vals = null_df[k].astype(float)
        ge = int((vals >= v).sum())
        p = (1 + ge) / (1 + len(vals))
        print(f"{k:12s}: null mean={vals.mean():7.3f} sd={vals.std():6.3f} "
              f"min={vals.min():5.0f} max={vals.max():5.0f} | real={v:5d} -> p={p:.3f}")

    top = real_df[real_df["score"] >= args.min_score]
    print(f"\n=== entries at score >= {args.min_score}: {len(top)} ===")
    cols = ["gene", "acc", "ebna1_start", "ebna1_end", "t_start", "t_end",
            "aln_len", "identity", "margin_z", "score"]
    if not top.empty:
        print(top[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
