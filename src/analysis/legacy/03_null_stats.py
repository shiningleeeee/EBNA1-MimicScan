#!/usr/bin/env python3
"""Null model v2: three statistics side by side.

For each real hit we compute three "expected number of null hits at least this
good" counts, all over the same B composition-preserving EBNA1 shuffles:

  1. E_raw     null hits with score            >= observed score
               (this is what the original package used for ranking)
  2. E_dens    null hits with score/aln_len    >= observed score/aln_len
               (removes the trivial "longer alignment scores higher" effect)
  3. E_band    null hits *of the same alignment-length band* with score
               >= observed score
               (length-stratified: asks "for an alignment this long, is this
               score unusual?")

E_raw conflates alignment length with evidence: a 291-aa alignment scoring 542
beats a 20-aa alignment scoring 60, although per-residue the latter is far more
surprising. E_band is the comparison that actually holds length fixed.

E_* are COUNTS, not p-values. 0 is right-censored at the resolution 1/B.

Usage: 03_null_stats.py --real a.tsv[,b.tsv] --out reports/null_stats_B200.tsv \
                         --B 200 --workers 10
"""
import argparse
import os
import subprocess
import tempfile
from multiprocessing import Pool

import numpy as np
import pandas as pd

BIN = "env/bin/mimic_scan"
QUERY = "data/cache/P03211.fasta"
DB = "data/cache/human_swissprot.fasta"
MAXSCORE = 800
DENS_SCALE = 100
MAXDENS = 1500
BAND_EDGES = [6, 15, 30, 60, 120, 250, 500, 10 ** 9]
NBAND = len(BAND_EDGES) - 1

Q_SEQ = ""


def read_query_seq(path):
    seq = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">") or not line:
            continue
        seq.append(line)
    return "".join(seq)


def band_of(L):
    return np.searchsorted(BAND_EDGES, L, side="right") - 1


def run_one(args):
    b, seed, tmpdir = args
    rng = np.random.default_rng(seed)
    letters = np.frombuffer(Q_SEQ.encode(), dtype=np.uint8)
    shuffled = rng.permutation(letters).tobytes().decode()
    qfa = os.path.join(tmpdir, f"s{b}.fa")
    out = os.path.join(tmpdir, f"s{b}.tsv")
    with open(qfa, "w") as fh:
        fh.write(f">shuffle_{b}\n")
        for i in range(0, len(shuffled), 60):
            fh.write(shuffled[i:i + 60] + "\n")
    r = subprocess.run([BIN, "--query", qfa, "--db", DB, "--mode", "seeded",
                        "--out", out], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"shuffle {b}: {r.stderr[-400:]}")
    df = pd.read_csv(out, sep="\t", usecols=["gscore", "uscore", "aln_len"],
                     low_memory=False)
    sc = df["gscore"].fillna(df["uscore"]).to_numpy()
    L = df["aln_len"].to_numpy()
    ok = (L > 0)
    sc, L = sc[ok], L[ok]
    # raw score histogram
    h_raw = np.bincount(np.clip(sc, 0, MAXSCORE).astype(int), minlength=MAXSCORE + 1)
    # density histogram
    dens = np.clip((sc / L * DENS_SCALE).astype(int), 0, MAXDENS - 1)
    h_dens = np.bincount(dens, minlength=MAXDENS)
    # per-band score histograms
    h_band = np.zeros((NBAND, MAXSCORE + 1), dtype=np.int32)
    bi = band_of(L)
    for k in range(NBAND):
        m = bi == k
        if m.any():
            h_band[k] = np.bincount(np.clip(sc[m], 0, MAXSCORE).astype(int),
                                    minlength=MAXSCORE + 1)
    os.remove(qfa)
    os.remove(out)
    return b, h_raw.astype(np.int64), h_dens.astype(np.int64), h_band


def sufsum(h):
    return np.cumsum(h[..., ::-1], axis=-1)[..., ::-1]


def main():
    global Q_SEQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--seed0", type=int, default=42)
    args = ap.parse_args()

    Q_SEQ = read_query_seq(QUERY)
    print(f"[null] query {len(Q_SEQ)} aa  B={args.B}  workers={args.workers}")

    parts = []
    for p in args.real.split(","):
        d = pd.read_csv(p, sep="\t", low_memory=False)
        print(f"[null] {p}: {len(d):,} hits")
        parts.append(d)
    real = pd.concat(parts, ignore_index=True)
    real["score"] = real["gscore"].fillna(real["uscore"]).astype(float)
    real["dens"] = real["score"] / real["aln_len"].clip(lower=1)
    real["band"] = band_of(real["aln_len"].to_numpy())

    tmpdir = tempfile.mkdtemp(prefix="nul2", dir="reports")
    jobs = [(b + 1, args.seed0 + b, tmpdir) for b in range(args.B)]
    H_raw = np.zeros((args.B, MAXSCORE + 1), dtype=np.int64)
    H_dens = np.zeros((args.B, MAXDENS), dtype=np.int64)
    H_band = np.zeros((args.B, NBAND, MAXSCORE + 1), dtype=np.int64)

    with Pool(args.workers) as pool:
        for k, (b, hr, hd, hb) in enumerate(pool.imap_unordered(run_one, jobs), 1):
            H_raw[b - 1] = hr
            H_dens[b - 1] = hd
            H_band[b - 1] = hb
            if k % 50 == 0 or k == args.B:
                tot = H_raw[:k].sum(axis=1)
                print(f"  {k}/{args.B}: null hits/shuffle mean={tot.mean():,.0f} "
                      f"cv={tot.std(ddof=1)/tot.mean():.4f}", flush=True)

    C_raw = sufsum(H_raw)      # (B, S)
    C_dens = sufsum(H_dens)    # (B, D)
    C_band = sufsum(H_band)    # (B, NBAND, S)

    si = np.clip(real["score"].to_numpy().astype(int), 0, MAXSCORE)
    di = np.clip((real["dens"].to_numpy() * DENS_SCALE).astype(int), 0, MAXDENS - 1)
    bi = real["band"].to_numpy()

    real["E_raw"] = C_raw[:, si].sum(axis=0)
    real["E_dens"] = C_dens[:, di].sum(axis=0)
    real["E_band"] = C_band[np.arange(args.B)[:, None], bi[None, :], si[None, :]].sum(axis=0)
    real["B"] = args.B

    null_tot = H_raw.sum(axis=1)
    print(f"\n== N1 null (seeded, composition-preserving shuffle of EBNA1) ==")
    print(f"null hits / shuffle : mean {null_tot.mean():,.1f}  sd {null_tot.std(ddof=1):,.1f}"
          f"  CV {null_tot.std(ddof=1)/null_tot.mean():.4f}")
    print(f"real hits           : {len(real):,}")
    print(f"ratio real/null     : {len(real)/null_tot.mean():.2f}x")

    print("\n== how many hits survive each statistic ==")
    for name in ["E_raw", "E_dens", "E_band"]:
        v = real[name].to_numpy()
        print(f"{name:7s}  ==0: {(v == 0).sum():7,} ({(v==0).mean():6.2%})   "
              f"<B: {(v < args.B).sum():7,} ({(v<args.B).mean():6.2%})   "
              f"median {np.median(v):12,.0f}")

    print("\n== benchmarks ==")
    for g in ["CRYAB", "HEPACAM", "MEF2B", "ANO2"]:
        sub = real[real["gene"] == g]
        if not len(sub):
            continue
        for name in ["E_raw", "E_dens", "E_band"]:
            b = sub.sort_values(name).iloc[0]
            print(f"  {g:8s} {name:7s} {int(b[name]):>10,}  (score {b['score']:.0f}, "
                  f"len {b['aln_len']}, dens {b['dens']:.2f})")

    print("\n== top 20 genes, ranked by E_band (length-stratified) ==")
    best = (real.sort_values("E_band").groupby("gene", as_index=False).first()
                .sort_values(["E_band", "score"], ascending=[True, False]))
    print(best[["gene", "acc", "score", "aln_len", "dens", "E_band", "E_raw"]]
          .head(20).to_string(index=False))

    print("\n== ranking disagreement ==")
    mr = real["E_raw"].rank(method="average")
    mb = real["E_band"].rank(method="average")
    print(f"spearman(E_raw, E_band) over all {len(real):,} hits: "
          f"{np.corrcoef(mr, mb)[0,1]:.3f}")
    for name in ["E_raw", "E_dens", "E_band"]:
        top = set(real.sort_values(name).groupby("gene").first()
                      .sort_values(name).head(50).index)
        if name == "E_raw":
            keep = top
            continue
        print(f"top-50 gene overlap E_raw vs {name}: {len(keep & top)}/50")

    real.to_csv(args.out, sep="\t", index=False)
    print(f"\n[out] {args.out}")
    os.rmdir(tmpdir)


if __name__ == "__main__":
    main()
