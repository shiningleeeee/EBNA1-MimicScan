#!/usr/bin/env python3
"""Summarise motif_mine counts: real (run 0) vs window-shuffled nulls (runs 1..B).

Usage:  env/venv/bin/python src/11_motif_null_summary.py [reports/gated_counts.tsv]

Note: read columns with r['name'] -- r.L / r.count collide with pandas
Series attributes (.count is a method) and raise
    TypeError: int() argument must be ... not 'method'
"""
import sys
import numpy as np
import pandas as pd

path = sys.argv[1] if len(sys.argv) > 1 else "reports/gated_counts.tsv"
counts = pd.read_csv(path, sep="\t")

real = counts[counts["run"] == 0]
null = counts[counts["run"] > 0]
B = null["run"].nunique()
piv = null.pivot_table(index=["L", "mismatch"], columns="run",
                       values="count", aggfunc="first").fillna(0)

rows = []
print(f"{'L':>3} {'mm':>3} {'real':>8} {'null_mean':>10} {'null_max':>9} "
      f"{'#null>=real':>12} {'p':>7}")
for _, r in real.iterrows():
    L, mm, n = int(r["L"]), int(r["mismatch"]), int(r["count"])
    if (L, mm) in piv.index:
        v = piv.loc[(L, mm)].values.astype(float)
        ge = int((v >= n).sum())
        p = (ge + 1) / (len(v) + 1)
        print(f"{L:>3} {mm:>3} {n:>8} {v.mean():>10.1f} {v.max():>9.0f} "
              f"{ge:>12} {p:>7.3f}")
        rows.append((L, mm, n, v.mean(), v.max(), ge, p))
    else:
        # null side never produced this cell at all
        print(f"{L:>3} {mm:>3} {n:>8} {0.0:>10.1f} {0.0:>9.0f} {0:>12} "
              f"  n/a   <- null all zero")
        rows.append((L, mm, n, 0.0, 0.0, 0, np.nan))

out = pd.DataFrame(rows, columns=["L", "mm", "real", "null_mean", "null_max",
                                  "n_ge", "p"])
out.to_csv("reports/gated_summary.tsv", sep="\t", index=False)
print(f"\nB = {B} shuffles; wrote reports/gated_summary.tsv")
