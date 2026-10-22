#!/usr/bin/env python3
"""MMseqs2 real-database search + composition-matched shuffle null.

The point of this stage: MMseqs2 supplies a *calibrated* E-value (composition
bias correction on, real database size 20,431), so the question "is any
EBNA1 <-> human alignment better than chance?" becomes answerable without
inventing a null by hand.

The null is the same question asked of a random sequence with EBNA1's exact
amino-acid composition: shuffle, search the same database, record how good the
best hit is. If the real best E-value sits inside the null distribution, the
sequence channel carries no signal for this query.

Usage: 05_mmseqs_null.py --B 50 --workers 6
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
DB = "data/cache/human_swissprot.fasta"
QUERY = "data/cache/P03211.fasta"
FMT = "target,evalue,bits,qstart,qend,tstart,tend,alnlen,pident,qaln,taln"
SEARCH = ["-s", "7.5", "--max-seqs", "20000", "-e", "10000",
          "--comp-bias-corr", "1", "-a", "--threads", "2"]
BENCH = {"P02511": "CRYAB", "Q14CZ8": "HEPACAM", "Q02080": "MEF2B",
         "Q9NQM4": "ANO2"}

Q_SEQ = ""


def read_query_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def run_search(qfa, workdir, tag, tdb):
    """Run one query fasta through MMseqs2, return a DataFrame of hits."""
    qdb = os.path.join(workdir, f"q_{tag}")
    out = os.path.join(workdir, f"r_{tag}")
    m8 = os.path.join(workdir, f"{tag}.m8")
    subprocess.run([MM, "createdb", qfa, qdb], check=True,
                   capture_output=True)
    subprocess.run([MM, "search", qdb, tdb, out, os.path.join(workdir, "t_" + tag)]
                   + SEARCH, check=True, capture_output=True)
    subprocess.run([MM, "convertalis", qdb, tdb, out, m8,
                    "--format-output", FMT, "--threads", "2"],
                   check=True, capture_output=True)
    if os.path.getsize(m8) == 0:
        return pd.DataFrame(columns=FMT.split(","))
    df = pd.read_csv(m8, sep="\t", header=None, names=FMT.split(","))
    for f in (qdb, out, m8):
        for suf in ("", ".dbtype", ".index", ".lookup", ".source", "_h",
                    "_h.dbtype", "_h.index"):
            p = f + suf
            if os.path.exists(p):
                os.remove(p)
    return df


def run_shuffle(args):
    b, seed, tdb = args
    workdir = tempfile.mkdtemp(prefix=f"sh{b}_", dir="work/mm")
    rng = np.random.default_rng(seed)
    letters = np.frombuffer(Q_SEQ.encode(), dtype=np.uint8)
    seq = rng.permutation(letters).tobytes().decode()
    qfa = os.path.join(workdir, "q.fa")
    with open(qfa, "w") as fh:
        fh.write(f">shuffle_{b}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    try:
        df = run_search(qfa, workdir, "s", tdb)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if df.empty:
        return b, np.inf, 0, 0, np.nan
    e = df["evalue"].to_numpy(float)
    return b, float(e.min()), int((e < 1).sum()), int((e < 100).sum()), float(df["bits"].max())


def main():
    global Q_SEQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=50)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed0", type=int, default=42)
    args = ap.parse_args()

    Q_SEQ = read_query_seq(QUERY)
    print(f"[mmseqs] query {len(Q_SEQ)} aa; db 20,431 human reviewed")

    real = pd.read_csv("reports/mmseqs_ebna1.m8", sep="\t", header=None,
                       names=FMT.split(","))
    real["target"] = real["target"].astype(str)
    print(f"\n== real EBNA1 vs human (MMseqs2, comp-bias-corr on) ==")
    print(f"alignments reported : {len(real):,}")
    print(f"best E-value        : {real['evalue'].min():.4g}")
    print(f"best bits           : {real['bits'].max():.1f}")
    for thr in [0.001, 0.01, 0.1, 1, 10, 100]:
        n = int((real["evalue"] < thr).sum())
        print(f"  E < {thr:<8g}: {n:5,}")
    print("\ntop 12 by E-value:")
    cols = ["target", "evalue", "bits", "alnlen", "pident", "qstart", "qend",
            "tstart", "tend"]
    print(real.sort_values("evalue")[cols].head(12).to_string(index=False))
    print("\nbenchmarks:")
    for acc, g in BENCH.items():
        sub = real[real["target"] == acc]
        if sub.empty:
            print(f"  {g:8s} {acc}: not reported by MMseqs2 at all")
        else:
            b = sub.sort_values("evalue").iloc[0]
            print(f"  {g:8s} {acc}: E={b['evalue']:.4g} bits={b['bits']:.1f} "
                  f"len={b['alnlen']:.0f} pid={b['pident']:.1f} "
                  f"q{b['qstart']}-{b['qend']} t{b['tstart']}-{b['tend']}")

    # ---- null ----
    tdb = os.path.abspath("work/mm/tdb")
    jobs = [(b + 1, args.seed0 + b, tdb) for b in range(args.B)]
    print(f"\n== shuffle null: {args.B} composition-matched random EBNA1 ==")
    rows = []
    with Pool(args.workers) as pool:
        for k, (b, bestE, n1, n100, bits) in enumerate(
                pool.imap_unordered(run_shuffle, jobs), 1):
            rows.append((b, bestE, n1, n100, bits))
            if k % 10 == 0 or k == args.B:
                be = np.array([r[1] for r in rows])
                print(f"  {k}/{args.B}: null best-E  median {np.median(be):.3g}  "
                      f"min {be.min():.3g}  max {be.max():.3g}", flush=True)
    null = pd.DataFrame(rows, columns=["b", "bestE", "n_E1", "n_E100", "best_bits"])
    null.to_csv("reports/mmseqs_null_B.csv", index=False)

    be = null["bestE"].to_numpy(float)
    print(f"\n== null summary ({args.B} shuffles) ==")
    print(f"null best E-value : median {np.median(be):.4g}  "
          f"IQR [{np.percentile(be,25):.4g}, {np.percentile(be,75):.4g}]  "
          f"min {be.min():.4g}")
    print(f"null alignments with E<1   : mean {null['n_E1'].mean():.1f} "
          f"(sd {null['n_E1'].std(ddof=1):.1f})")
    print(f"null alignments with E<100 : mean {null['n_E100'].mean():.1f} "
          f"(sd {null['n_E100'].std(ddof=1):.1f})")

    real_best = real["evalue"].min()
    print(f"\n== verdict ==")
    print(f"real best E     = {real_best:.4g}")
    print(f"null best E med = {np.median(be):.4g}")
    frac = float((be <= real_best).mean())
    print(f"P(null best E <= real best E) = {frac:.3f}  over {len(be)} shuffles")
    print(f"real alignments E<1 = {int((real['evalue']<1).sum())}  vs null mean "
          f"{null['n_E1'].mean():.1f}")


if __name__ == "__main__":
    main()
