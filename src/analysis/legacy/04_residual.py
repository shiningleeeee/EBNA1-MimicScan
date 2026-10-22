#!/usr/bin/env python3
"""What, if anything, survives the null once composition artefacts are removed?

Reads reports/null_stats_B200.tsv (E_raw / E_dens / E_band per hit) and applies
the filters the original package *defined* (tier A/B/C, GAr flag, low-complexity)
but never used for ranking.

Filters, all computed from the aligned strings themselves:
  gar_q      fraction of G/A/R in the aligned query segment (GAr tract drive)
  entropy_t  Shannon entropy (bits) of the aligned target segment
  maxfrac_t  largest single-residue fraction of the aligned target segment
  ga_t       G+A fraction of the aligned target segment
  spr_t      S+P+R fraction of the aligned target segment

Usage: 04_residual.py reports/null_stats_B200.tsv
"""
import sys

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"
BENCH = {"CRYAB", "HEPACAM", "MEF2B", "ANO2"}


def _clean(s):
    return (s if isinstance(s, str) else "").replace("-", "")


def comp(s, letters):
    s = _clean(s)
    if not s:
        return np.nan
    return sum(ch in letters for ch in s) / len(s)


def entropy(s):
    s = _clean(s)
    if not s:
        return np.nan
    v, c = np.unique(list(s), return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def maxfrac(s):
    s = _clean(s)
    if not s:
        return np.nan
    v, c = np.unique(list(s), return_counts=True)
    return float(c.max() / c.sum())


def main(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"[read] {path}: {len(df):,} hits, B={int(df['B'].iloc[0])}")
    B = int(df["B"].iloc[0])

    df["gar_q"] = df["gqseq"].map(lambda s: comp(s, set("GAR")))
    df["ga_t"] = df["gtseq"].map(lambda s: comp(s, set("GA")))
    df["spr_t"] = df["gtseq"].map(lambda s: comp(s, set("SPR")))
    df["entropy_t"] = df["gtseq"].map(entropy)
    df["maxfrac_t"] = df["gtseq"].map(maxfrac)
    df["ident"] = df["n_id"] / df["aln_len"].clip(lower=1)

    print("\n== survivor counts by statistic ==")
    for name in ["E_raw", "E_dens", "E_band"]:
        v = df[name].to_numpy()
        print(f"{name:7s} ==0 {np.sum(v==0):7,} ({(v==0).mean():6.2%})   "
              f"<B {np.sum(v<B):7,} ({(v<B).mean():6.2%})   median {np.median(v):12,.0f}")

    steps = [
        ("all hits", pd.Series(True, index=df.index)),
        ("aln_len >= 20", df["aln_len"] >= 20),
        ("query segment not GAr-rich (gar_q < 0.50)",
         df["gar_q"] < 0.50),
        ("target segment ordered (entropy >= 2.8 & maxfrac <= 0.30)",
         (df["entropy_t"] >= 2.8) & (df["maxfrac_t"] <= 0.30)),
        ("not a GA/SPR repeat target (ga_t < 0.40 & spr_t < 0.50)",
         (df["ga_t"] < 0.40) & (df["spr_t"] < 0.50)),
        ("beats its length-matched null (E_band < B)", df["E_band"] < B),
        ("identity >= 0.35", df["ident"] >= 0.35),
    ]
    mask = pd.Series(True, index=df.index)
    print("\n== cumulative filter funnel ==")
    for label, m in steps:
        mask &= m.fillna(False)
        n = int(mask.sum())
        print(f"  {label:58s} {n:7,}  ({n/len(df):6.3%})  "
              f"[step: {int(m.sum()):,}]")

    print(f"\n== survivors after the full funnel ==")
    surv = df[mask].sort_values(["E_band", "score"], ascending=[True, False])
    print(f"n = {len(surv)}")
    if len(surv) == 0:
        print("  nothing survives: every hit is either a composition artefact or "
              "not better than its length-matched null")
        return
    print(f"unique genes: {surv['gene'].nunique()}")
    print(surv[["gene", "acc", "mode", "score", "aln_len", "ident", "gar_q",
                "entropy_t", "E_band", "E_raw", "qstart", "qend", "tstart", "tend"]]
          .head(40).to_string(index=False))
    print(f"\nbenchmarks among survivors: "
          f"{sorted(set(surv['gene']) & BENCH) or 'none'}")

    print("\n== where do the published benchmarks sit on each axis? ==")
    for g in ["CRYAB", "HEPACAM", "MEF2B", "ANO2"]:
        sub = df[df["gene"] == g]
        if not len(sub):
            continue
        b = sub.sort_values("E_band").iloc[0]
        print(f"  {g:8s} best E_band {int(b['E_band']):>8,}  score {b['score']:5.0f} "
              f"len {b['aln_len']:4d} ident {b['ident']:.2f} gar_q {b['gar_q']:.2f} "
              f"ent {b['entropy_t']:.2f}")


if __name__ == "__main__":
    main(sys.argv[1])
