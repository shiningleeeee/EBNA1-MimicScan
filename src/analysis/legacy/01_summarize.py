#!/usr/bin/env python3
"""Summarise a mimic_scan run: composition of the hit set, top genes,
and presence of the three published EBNA1 cross-reactive benchmarks.

Usage: 01_summarize.py reports/seeded_raw.tsv reports/exhaustive_raw.tsv
"""
import sys
import pandas as pd

BENCH = {"CRYAB": "P02511", "HEPACAM": "Q14CZ8", "MEF2B": "Q02080", "ANO2": "Q9NQM4"}
GAR = set("GAR")


def gar_frac(s: str) -> float:
    s = s.replace("-", "")
    if not s:
        return 0.0
    return sum(ch in GAR for ch in s) / len(s)


def main(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t")
        print(f"[read] {p}: {len(df):,} rows")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    print(f"\n== totals ==\nrows: {len(df):,}")
    print(df["mode"].value_counts().to_string())
    print(f"unique acc: {df['acc'].nunique():,}   unique gene: {df['gene'].nunique():,}")

    df["gar_frac"] = df["gqseq"].fillna("").map(gar_frac)
    df["aln_ident"] = df["n_id"] / df["aln_len"]
    print("\n== aligned-query GAr fraction ==")
    bins = pd.cut(df["gar_frac"], [-0.01, 0.001, 0.25, 0.5, 0.75, 1.01])
    print(bins.value_counts().sort_index().to_string())
    print(f"\ngar_frac >= 0.75 : {(df['gar_frac'] >= 0.75).sum():,}"
          f"  ({(df['gar_frac'] >= 0.75).mean():.1%})")
    print(f"aln_len median   : {df['aln_len'].median():.0f}")
    print(f"identity median  : {df['aln_ident'].median():.3f}")

    print("\n== top 25 genes by best score ==")
    best = (df.assign(score=df["gscore"].fillna(df["uscore"]))
              .sort_values("score", ascending=False)
              .groupby("gene", as_index=False).first())
    best = best.sort_values("score", ascending=False)
    cols = ["gene", "acc", "score", "aln_len", "aln_ident", "gar_frac"]
    print(best[cols].head(25).to_string(index=False))

    print("\n== published benchmarks ==")
    for g, acc in BENCH.items():
        sub = df[(df["gene"] == g) | (df["acc"] == acc)]
        if len(sub) == 0:
            print(f"{g:8s} {acc}: ABSENT")
            continue
        b = sub.assign(score=sub["gscore"].fillna(sub["uscore"])).sort_values(
            "score", ascending=False).iloc[0]
        print(f"{g:8s} {acc}: {len(sub):3d} hits, best score {b['score']:4.0f}, "
              f"q {b['qstart']}-{b['qend']} / t {b['tstart']}-{b['tend']}, "
              f"aln_len {b['aln_len']}, ident {b['aln_ident']:.2f}, gar {b['gar_frac']:.2f}")

    top50 = set(best.head(50)["gene"])
    print(f"\ntop-50 genes GAr-laden: "
          f"{sum(best.set_index('gene').loc[g, 'gar_frac'] >= 0.75 for g in top50)}/50")


if __name__ == "__main__":
    main(sys.argv[1:])
