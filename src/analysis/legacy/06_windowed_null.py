#!/usr/bin/env python3
"""Local-composition null: shuffle EBNA1 inside windows.

Motivation. A global shuffle preserves the amino-acid composition but destroys
the fact that EBNA1's Gly/Ala residues are clustered in a long tract. If the
real sequence's advantage over a global shuffle (best E 3.09 vs null median 518)
comes from "having a tract that can align to some human protein's tract", then
shuffling *within windows* -- which keeps every window's composition and only
rearranges locally -- reproduces that advantage, and the real sequence stops
being special.

Window size is the knob: w = L is the global shuffle of the previous stage,
w = 1 is the real sequence. Watching best-E rise (or not) as w shrinks tells us
whether the signal is tract composition or something order-dependent.

Usage: 06_windowed_null.py --B 25 --windows 20,40,80
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

MM = "env/bin/mmseqs"
FMT = "target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend"
SEARCH = ["-s", "7.5", "--max-seqs", "20000", "-e", "10000",
          "--comp-bias-corr", "1", "-a", "--threads", "2"]
REAL_BEST = 3.088

Q_SEQ = ""


def read_query_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def windowed_shuffle(seq, w, rng):
    out = list(seq)
    for s in range(0, len(seq), w):
        blk = list(seq[s:s + w])
        rng.shuffle(blk)
        out[s:s + w] = blk
    return "".join(out)


def run_search(qfa, workdir, tag, tdb):
    qdb = os.path.join(workdir, f"q_{tag}")
    out = os.path.join(workdir, f"r_{tag}")
    m8 = os.path.join(workdir, f"{tag}.m8")
    subprocess.run([MM, "createdb", qfa, qdb], check=True, capture_output=True)
    subprocess.run([MM, "search", qdb, tdb, out,
                    os.path.join(workdir, "t_" + tag)] + SEARCH,
                   check=True, capture_output=True)
    subprocess.run([MM, "convertalis", qdb, tdb, out, m8,
                    "--format-output", FMT, "--threads", "2"],
                   check=True, capture_output=True)
    if os.path.getsize(m8) == 0:
        return np.array([]), np.array([])
    df = pd.read_csv(m8, sep="\t", header=None, names=FMT.split(","))
    return df["evalue"].to_numpy(float), df["bits"].to_numpy(float)


def run_one(args):
    mode, w, b, seed, tdb = args
    workdir = tempfile.mkdtemp(prefix="w", dir="work/mm")
    rng = np.random.default_rng(seed)
    seq = windowed_shuffle(Q_SEQ, w, rng)
    qfa = os.path.join(workdir, "q.fa")
    with open(qfa, "w") as fh:
        fh.write(f">null_w{w}_{b}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    try:
        e, bits = run_search(qfa, workdir, "s", tdb)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if len(e) == 0:
        return mode, w, b, np.inf, 0, 0, np.nan
    return mode, w, b, float(e.min()), int((e < 10).sum()), int((e < 100).sum()), float(bits.max())


def main():
    global Q_SEQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=25)
    ap.add_argument("--windows", default="20,40,80")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--seed0", type=int, default=1000)
    args = ap.parse_args()

    Q_SEQ = read_query_seq("data/cache/P03211.fasta")
    tdb = os.path.abspath("work/mm/tdb")
    windows = [int(x) for x in args.windows.split(",")]

    jobs = []
    for w in windows:
        for b in range(args.B):
            jobs.append((f"w{w}", w, b + 1, args.seed0 + w * 1000 + b, tdb))
    print(f"[null] {len(jobs)} runs: windows={windows} x B={args.B}, "
          f"query {len(Q_SEQ)} aa", flush=True)

    rows = []
    with Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(run_one, jobs), 1):
            rows.append(r)
            if k % 20 == 0 or k == len(jobs):
                print(f"  {k}/{len(jobs)}", flush=True)

    df = pd.DataFrame(rows, columns=["mode", "w", "b", "bestE", "n_E10",
                                     "n_E100", "best_bits"])
    df.to_csv("reports/null_windowed_B.csv", index=False)

    print("\n== windowed-shuffle null: best E-value vs window size ==")
    print(f"{'window':>7} {'n':>4} {'median':>10} {'min':>10} {'p25':>10} "
          f"{'p75':>10} {'P(null<=real)':>14} {'E<10':>6} {'E<100':>7}")
    print(f"{'real':>7} {1:>4} {REAL_BEST:>10.3g} {'-':>10} {'-':>10} {'-':>10} "
          f"{'-':>14} {2:>6} {12:>7}")
    for w in windows:
        sub = df[df["w"] == w]
        be = sub["bestE"].to_numpy(float)
        frac = float((be <= REAL_BEST).mean())
        print(f"{w:>7} {len(be):>4} {np.median(be):>10.3g} {be.min():>10.3g} "
              f"{np.percentile(be,25):>10.3g} {np.percentile(be,75):>10.3g} "
              f"{frac:>14.3f} {sub['n_E10'].mean():>6.1f} {sub['n_E100'].mean():>7.1f}")

    print("\n== read-out ==")
    for w in windows:
        sub = df[df["w"] == w]
        med = np.median(sub["bestE"].to_numpy(float))
        print(f"w={w:>3d}: null median best-E {med:8.3g}  "
              f"= {med/REAL_BEST:6.1f}x worse than real")
    glob = pd.read_csv("reports/mmseqs_null_B.csv")
    gmed = np.median(glob["bestE"].to_numpy(float))
    print(f"w=641: null median best-E {gmed:8.3g}  = {gmed/REAL_BEST:6.1f}x worse than real")
    print("\nIf the ratio collapses toward 1 as w shrinks, EBNA1's advantage over a "
          "global shuffle is local composition (tract matching), not residue order.")


if __name__ == "__main__":
    main()
