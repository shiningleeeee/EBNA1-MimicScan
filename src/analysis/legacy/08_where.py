#!/usr/bin/env python3
"""Where on EBNA1 do the (few) real-database alignments actually sit?

The 363 MMseqs2 alignments are the only ones a calibrated search engine is
willing to report at all. This script asks two questions about them:

  1. Do they spread over EBNA1, or do they all pile onto one low-complexity
     region? (If the latter, the "mimicry" is one repeat tract matching many
     different human low-complexity regions.)
  2. Are the published cross-reactive benchmarks and the step-1 survivor
     (CBX8) even present?

Usage: 08_where.py
"""
import numpy as np
import pandas as pd

M8 = "reports/mmseqs_ebna1.m8"
DB = "data/cache/human_swissprot.fasta"
Q = "data/cache/P03211.fasta"
WATCH = {"P02511": "CRYAB*", "Q14CZ8": "HEPACAM*", "Q02080": "MEF2B*",
         "Q9NQM4": "ANO2*", "Q9HC52": "CBX8 (step-1 survivor)"}


def read_fa(path):
    recs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                recs[name] = "".join(buf)
            parts = line[1:].split("|")
            name = parts[1] if len(parts) > 2 else line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if name:
        recs[name] = "".join(buf)
    return recs


def gene_of(acc, db_genes):
    return db_genes.get(acc, acc)


def main():
    seq = read_fa(Q)["P03211"]
    L = len(seq)
    db_genes = {}
    nm = None
    for line in open(DB):
        if line.startswith(">"):
            f = line[1:].split("|")
            nm = f[1] if len(f) > 2 else line[1:].split()[0]
            g = [t[3:] for t in line.split() if t.startswith("GN=")]
            db_genes[nm] = g[0].split(";")[0] if g else nm
    print(f"[q] EBNA1 {L} aa; db genes mapped: {len(db_genes):,}")

    df = pd.read_csv(M8, sep="\t", header=None,
                     names=["target", "evalue", "bits", "qstart", "qend",
                            "tstart", "tend", "alnlen", "pident", "qaln", "taln"])
    df["gene"] = df["target"].map(lambda a: gene_of(a, db_genes))
    print(f"[m8] {len(df)} alignments, {df['target'].nunique()} distinct targets")
    print(f"     best E {df['evalue'].min():.4g};  median alnlen {df['alnlen'].median():.0f}; "
          f"median pident {df['pident'].median():.1f}%")

    print("\n== where on EBNA1 (50-aa bins, by alignment midpoint) ==")
    mid = ((df["qstart"] + df["qend"]) / 2).astype(int)
    bins = np.arange(0, L + 50, 50)
    hist, _ = np.histogram(mid, bins=bins)
    ga = np.array([sum(c in "GA" for c in seq[i:i + 50]) / 50 for i in bins[:-1]])
    gap = np.array([sum(c in "GAP" for c in seq[i:i + 50]) / 50 for i in bins[:-1]])
    print(f"{'bin':>9} {'n_hits':>7} {'bar':>22} {'G+A':>6} {'G+A+P':>6}")
    for i, cnt in enumerate(hist):
        bar = "#" * int(cnt / max(hist.max(), 1) * 20)
        print(f"{bins[i]:>4}-{bins[i]+49:<4} {cnt:>7} {bar:<22} {ga[i]:>6.2f} {gap[i]:>6.2f}")

    top = hist.max()
    if top:
        peak = bins[hist.argmax()]
        inside = np.histogram(mid, bins=[peak, peak + 50])[0][0]
        print(f"\npeak bin {peak}-{peak+49}: {inside}/{len(df)} = {inside/len(df):.1%} "
              f"of all alignments")
    cov = np.zeros(L, dtype=bool)
    for _, r in df.iterrows():
        cov[int(r.qstart) - 1:int(r.qend)] = True
    print(f"EBNA1 residues touched by ANY reported alignment: {cov.sum()}/{L} "
          f"({cov.sum()/L:.1%})")

    print("\n== composition of the region every hit piles onto ==")
    for lo, hi in [(0, 90), (90, 330), (330, 460), (460, 530), (530, L)]:
        s = seq[lo:hi]
        g = sum(c in "GA" for c in s) / len(s)
        p = sum(c in "GAP" for c in s) / len(s)
        r = sum(c == "R" for c in s) / len(s)
        print(f"  {lo:>3}-{hi:<3} len {len(s):>3}  G+A {g:.2f}  G+A+P {p:.2f}  R {r:.2f}")

    print("\n== the alignments themselves (top 10) ==")
    cols = ["gene", "target", "evalue", "bits", "alnlen", "pident",
            "qstart", "qend", "tstart", "tend"]
    print(df.sort_values("evalue")[cols].head(10).to_string(index=False))

    print("\n== do the benchmarks / step-1 survivor appear at all? ==")
    for acc, label in WATCH.items():
        sub = df[df["target"] == acc]
        if sub.empty:
            print(f"  {label:26s} {acc}: NOT reported by MMseqs2")
        else:
            b = sub.sort_values("evalue").iloc[0]
            print(f"  {label:26s} {acc}: E={b['evalue']:.4g} bits={b['bits']:.0f} "
                  f"len={b['alnlen']:.0f} pid={b['pident']:.0f}%")

    print("\n== targets: how many are distinct, and are they low-complexity? ==")
    def lowc(s):
        s = s.replace("-", "")
        if not s:
            return np.nan
        v, c = np.unique(list(s), return_counts=True)
        p = c / c.sum()
        return float(-(p * np.log2(p)).sum())
    df["ent_t"] = df["taln"].map(lowc)
    df["ga_t"] = df["taln"].map(
        lambda s: sum(c in "GAP" for c in s.replace("-", "")) / max(len(s.replace("-", "")), 1))
    print(f"  target segments with entropy < 2.8 : {int((df['ent_t']<2.8).sum())}/{len(df)}"
          f"  ({(df['ent_t']<2.8).mean():.1%})")
    print(f"  target segments with G+A+P >= 0.50 : {int((df['ga_t']>=0.5).sum())}/{len(df)}"
          f"  ({(df['ga_t']>=0.5).mean():.1%})")
    print(f"  median target-segment entropy      : {df['ent_t'].median():.2f} bits")


if __name__ == "__main__":
    main()
