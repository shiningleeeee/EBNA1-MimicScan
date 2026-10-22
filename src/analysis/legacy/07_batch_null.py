#!/usr/bin/env python3
"""Batched local-composition null -- all shuffles in ONE MMseqs2 search.

Same question as 06, but the real sequence and every null variant go into a
single multi-FASTA query database, so they are searched in one run against the
same target index and compared like-for-like. This raises B to 200 per window
at negligible cost and gives p-values with 1/200 resolution.

Null family: shuffle inside windows of size w. w = 641 is the global shuffle
(preserves composition only); small w preserves progressively more local
composition. If the real sequence's advantage over a global shuffle is local
composition, it disappears as w shrinks.

Usage: 07_batch_null.py --B 200 --windows 10,20,40
"""
import argparse
import os
import subprocess

import numpy as np
import pandas as pd

MM = "env/bin/mmseqs"
FMT = "query,target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend"
SEARCH = ["-s", "7.5", "--max-seqs", "20000", "-e", "10000",
          "--comp-bias-corr", "1", "-a", "--threads", "14"]
WORK = "work/batch"


def read_query_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def windowed_shuffle(seq, w, rng):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--windows", default="10,20,40")
    ap.add_argument("--seed0", type=int, default=7000)
    args = ap.parse_args()

    seq = read_query_seq("data/cache/P03211.fasta")
    windows = [int(x) for x in args.windows.split(",")]
    os.makedirs(WORK, exist_ok=True)

    qfa = os.path.join(WORK, "queries.fa")
    meta = []
    with open(qfa, "w") as fh:
        def emit(name, s):
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
        emit("real", seq)
        meta.append(("real", "real", -1))
        for w in windows + [len(seq)]:
            rng = np.random.default_rng(args.seed0 + w)
            for b in range(args.B):
                nm = f"w{w}_{b:04d}"
                emit(nm, windowed_shuffle(seq, w, rng))
                meta.append((nm, f"w{w}", b))
    print(f"[batch] {len(meta)} queries ({len(windows)+1} families x B={args.B} + real)")

    for f in ("qdb", "res", "tmp"):
        p = os.path.join(WORK, f)
        subprocess.run(["rm", "-rf", p + "*"], check=False)
    subprocess.run([MM, "createdb", qfa, os.path.join(WORK, "qdb")],
                   check=True, capture_output=True)
    subprocess.run([MM, "search", os.path.join(WORK, "qdb"),
                    os.path.abspath("work/mm/tdb"), os.path.join(WORK, "res"),
                    os.path.join(WORK, "tmp")] + SEARCH,
                   check=True, capture_output=True)
    m8 = "reports/batch_null.m8"
    subprocess.run([MM, "convertalis", os.path.join(WORK, "qdb"),
                    os.path.abspath("work/mm/tdb"), os.path.join(WORK, "res"), m8,
                    "--format-output", FMT, "--threads", "14"],
                   check=True, capture_output=True)

    df = pd.read_csv(m8, sep="\t", header=None, names=FMT.split(","))
    meta_df = pd.DataFrame(meta, columns=["query", "family", "b"])
    print(f"[batch] {len(df):,} alignments over {df['query'].nunique():,} queries")

    g = df.groupby("query")
    per = pd.DataFrame({
        "bestE": g["evalue"].min(),
        "best_bits": g["bits"].max(),
        "n_E10": g["evalue"].apply(lambda s: int((s < 10).sum())),
        "n_E100": g["evalue"].apply(lambda s: int((s < 100).sum())),
        "n": g.size(),
    }).reset_index()
    per = meta_df.merge(per, on="query", how="left")
    per["bestE"] = per["bestE"].fillna(np.inf)
    per["n_E10"] = per["n_E10"].fillna(0)
    per["n_E100"] = per["n_E100"].fillna(0)
    per["n"] = per["n"].fillna(0)
    per.to_csv("reports/batch_null_per_query.csv", index=False)

    realE = float(per.loc[per["family"] == "real", "bestE"].iloc[0])
    print(f"\nreal EBNA1 best E = {realE:.4g}")

    print("\n== best E-value by null family (bias-corrected real-database search) ==")
    hdr = (f"{'family':>10} {'n':>4} {'median':>10} {'p10':>10} {'p90':>10} "
           f"{'min':>9} {'vs real':>9} {'P(null<=real)':>14} {'E<10':>6} {'E<100':>7}")
    print(hdr)
    print(f"{'real':>10} {1:>4} {realE:>10.3g} {'':>10} {'':>10} {'':>9} "
          f"{1.0:>9.2f} {'-':>14} {int(per.loc[per['family']=='real','n_E10'].iloc[0]):>6} "
          f"{int(per.loc[per['family']=='real','n_E100'].iloc[0]):>7}")
    for fam in [f"w{w}" for w in windows] + [f"w{len(seq)}"]:
        sub = per[per["family"] == fam]
        be = sub["bestE"].to_numpy(float)
        med = np.median(be)
        p = float((be <= realE).mean())
        print(f"{fam:>10} {len(be):>4} {med:>10.3g} {np.percentile(be,10):>10.3g} "
              f"{np.percentile(be,90):>10.3g} {be.min():>9.3g} {med/realE:>8.1f}x "
              f"{p:>14.3f} {sub['n_E10'].mean():>6.1f} {sub['n_E100'].mean():>7.1f}")

    print("\n== interpretation ==")
    print("A null family whose P(null<=real) is not small means: real EBNA1 is")
    print("indistinguishable from a sequence with the same *local* composition.")
    print("The family that first makes P small identifies which level of sequence")
    print("organisation carries the apparent signal.")


if __name__ == "__main__":
    main()
