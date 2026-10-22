#!/usr/bin/env python3
"""N1 null model: composition-preserving shuffles of EBNA1 + seeded rescan.

Faithful reproduction of the original package's N1 null, rebuilt so that the
per-shuffle hit counts are recorded (the shipped package never delivered the
raw hit tables, so its null could not be re-derived).

E_emp(h) = sum over shuffles b of #{null alignments in run b with score >= score(h)}

This is a COUNT, not a p-value: E_emp = 0 means "no null run produced a hit this
good", which is right-censored at 1/B of the observed mass -- exactly the
interpretation the original docs demand.

Usage: 02_null_seeded.py --real reports/seeded_raw.tsv[,more.tsv] \
                          --out reports/null_seeded_B200.tsv --B 200 --workers 8
"""
import argparse
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

BIN = "env/bin/mimic_scan"
QUERY = "data/cache/P03211.fasta"
DB = "data/cache/human_swissprot.fasta"
MAXSCORE = 800


def read_query_seq(path):
    seq, keep = [], False
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            keep = True
            continue
        if keep:
            seq.append(line)
    return "".join(seq)


def run_one(args):
    """Shuffle EBNA1 with seed b, rescan, return a count histogram of scores."""
    b, seed, tmpdir = args
    rng = np.random.default_rng(seed)
    letters = np.frombuffer(Q_SEQ.encode(), dtype=np.uint8)
    shuffled = rng.permutation(letters).tobytes().decode()
    qfa = os.path.join(tmpdir, f"shuf_{b}.fa")
    out = os.path.join(tmpdir, f"shuf_{b}.tsv")
    with open(qfa, "w") as fh:
        fh.write(f">shuffle_{b}\n")
        for i in range(0, len(shuffled), 60):
            fh.write(shuffled[i:i + 60] + "\n")
    r = subprocess.run(
        [BIN, "--query", qfa, "--db", DB, "--mode", "seeded", "--out", out],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"shuffle {b} failed: {r.stderr[-500:]}")
    df = pd.read_csv(out, sep="\t", usecols=["gscore", "uscore"], low_memory=False)
    sc = df["gscore"].fillna(df["uscore"]).to_numpy()
    hist = np.bincount(np.clip(sc, 0, MAXSCORE).astype(int), minlength=MAXSCORE + 1)
    os.remove(qfa)
    os.remove(out)
    return b, hist.astype(np.int64)


Q_SEQ = ""


def main():
    global Q_SEQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=42)
    args = ap.parse_args()

    Q_SEQ = read_query_seq(QUERY)
    print(f"[null] query {len(Q_SEQ)} aa, B={args.B}, workers={args.workers}")

    real_parts = []
    for p in args.real.split(","):
        df = pd.read_csv(p, sep="\t", low_memory=False)
        real_parts.append(df)
        print(f"[null] real hits from {p}: {len(df):,}")
    real = pd.concat(real_parts, ignore_index=True)
    real["score"] = real["gscore"].fillna(real["uscore"]).astype(float)

    tmpdir = tempfile.mkdtemp(prefix="nullB", dir="reports")
    jobs = [(b + 1, args.seed0 + b, tmpdir) for b in range(args.B)]
    hists = np.zeros((args.B, MAXSCORE + 1), dtype=np.int64)
    with Pool(args.workers) as pool:
        for k, (b, hist) in enumerate(pool.imap_unordered(run_one, jobs), 1):
            hists[b - 1] = hist
            if k % 25 == 0 or k == args.B:
                tot = hists[:k].sum(axis=1)
                print(f"  {k}/{args.B} shuffles done; null hits/shuffle "
                      f"mean={tot.mean():,.0f} min={tot.min():,} max={tot.max():,}",
                      flush=True)

    # suffix sum: cum[b, s] = #null hits in run b with score >= s
    cum = np.cumsum(hists[:, ::-1], axis=1)[:, ::-1]

    s_idx = np.clip(real["score"].to_numpy().astype(int), 0, MAXSCORE)
    E = cum[:, s_idx].sum(axis=0).astype(np.int64)
    real["E_emp"] = E
    real["B"] = args.B

    null_tot = hists.sum(axis=1)
    print(f"\n== null hit count per shuffle ==\n"
          f"mean {null_tot.mean():,.1f}  sd {null_tot.std(ddof=1):.1f}  "
          f"CV {null_tot.std(ddof=1)/null_tot.mean():.4f}")
    print(f"real seeded hits: {len(real):,}")

    print("\n== E_emp distribution over real hits ==")
    qs = [0, 1, 5, 10, 25, 50, 75, 90, 100]
    print(real["E_emp"].describe(percentiles=[q / 100 for q in qs]).to_string())
    print(f"E_emp == 0 : {(real['E_emp'] == 0).sum():,} "
          f"({(real['E_emp'] == 0).mean():.1%})")
    print(f"E_emp >= {args.B}: {(real['E_emp'] >= args.B).sum():,} "
          f"({(real['E_emp'] >= args.B).mean():.1%})  <- no better than null")

    for g in ["CRYAB", "HEPACAM", "MEF2B", "ANO2"]:
        sub = real[real["gene"] == g]
        if len(sub):
            b = sub.sort_values("score", ascending=False).iloc[0]
            print(f"benchmark {g:8s}: best score {b['score']:.0f}, E_emp {int(b['E_emp']):,}")

    print("\n== top 15 genes by E_emp (lowest = furthest from null) ==")
    best = (real.sort_values("E_emp")
                .groupby("gene", as_index=False).first()
                .sort_values("E_emp"))
    print(best[["gene", "acc", "score", "E_emp"]].head(15).to_string(index=False))

    real.to_csv(args.out, sep="\t", index=False)
    print(f"\n[out] {args.out}  ({len(real):,} rows with E_emp)")
    os.rmdir(tmpdir)


if __name__ == "__main__":
    main()
