#!/usr/bin/env python3
"""Mine EBNA1 mimicry candidates under a rule that recovers the known positives,
then ask how many candidates the same rule returns for composition-matched
shuffles of EBNA1.  The ratio is the only thing that says whether the list is
worth a bench's time.

Rule R was set so that all four experimentally validated cross-reactivities
(CRYAB, HEPACAM, MEF2B, ANO2) fall inside it -- that is the recall requirement
of a retrieval task, and it is stated as such rather than disguised as a
significance cut.

Usage:
  10_mine_candidates.py real      # apply R to reports/exhaustive_raw.tsv
  10_mine_candidates.py null N    # generate N window-shuffled EBNA1 queries,
                                  # scan them, apply the same R
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

SEQ = "data/cache/P03211.fasta"
DB = "data/cache/human_swissprot.fasta"
SCAN = "env/bin/mimic_scan"
BENCH = {"CRYAB", "HEPACAM", "MEF2B", "ANO2"}

# --- rule R (recall-oriented, calibrated below by a null rate) ---
R = dict(aln_min=20, aln_max=90, ident_min=0.33, nid_min=7,
         gar_q_max=0.70, ga_t_max=0.40, ent_t_min=2.4, mf_t_max=0.35)
GAR_LO, GAR_HI = 90, 330          # the Gly-Ala block: immune-evasive, excluded


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def clean(s):
    return (s if isinstance(s, str) else "").replace("-", "")


def comp(s, letters):
    s = clean(s)
    return np.nan if not s else sum(c in letters for c in s) / len(s)


def entropy(s):
    s = clean(s)
    if not s:
        return np.nan
    _, c = np.unique(list(s), return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def maxfrac(s):
    s = clean(s)
    if not s:
        return np.nan
    _, c = np.unique(list(s), return_counts=True)
    return c.max() / c.sum()


def annotate(df):
    df = df.copy()
    df["ident"] = df.n_id / df.aln_len
    df["gar_q"] = [comp(x, "GAR") for x in df.gqseq]
    df["ent_t"] = [entropy(x) for x in df.gtseq]
    df["mf_t"] = [maxfrac(x) for x in df.gtseq]
    df["ga_t"] = [comp(x, "GA") for x in df.gtseq]
    df["spr_t"] = [comp(x, "SPR") for x in df.gtseq]
    df["in_gar_block"] = (df.qstart <= GAR_HI) & (df.qend >= GAR_LO)
    return df


def apply_rule(df):
    m = ((df.aln_len >= R["aln_min"]) & (df.aln_len <= R["aln_max"])
         & (df.ident >= R["ident_min"]) & (df.n_id >= R["nid_min"])
         & (df.gar_q < R["gar_q_max"]) & (df.ga_t < R["ga_t_max"])
         & (df.ent_t >= R["ent_t_min"]) & (df.mf_t <= R["mf_t_max"])
         & (~df.in_gar_block))
    return df[m]


def dedup_genes(df):
    """One row per gene: its best-scoring candidate under R."""
    score_col = "score" if "score" in df.columns else "gscore"
    return (df.sort_values([score_col, "ident"], ascending=False)
              .drop_duplicates("gene", keep="first"))


def window_shuffle(seq, w, rng):
    out = []
    for i in range(0, len(seq), w):
        chunk = list(seq[i:i + w])
        rng.shuffle(chunk)
        out.extend(chunk)
    return "".join(out)


def run_null(n, w=20, seed=42):
    rng = np.random.default_rng(seed)
    seq = read_seq(SEQ)
    os.makedirs("work/mine", exist_ok=True)
    counts, genes = [], []
    for i in range(n):
        q = window_shuffle(seq, w, rng)
        qf = f"work/mine/q{i}.fasta"
        with open(qf, "w") as fh:
            fh.write(f">shuf{i}\n")
            for j in range(0, len(q), 60):
                fh.write(q[j:j + 60] + "\n")
        out = f"work/mine/s{i}.tsv"
        subprocess.run([SCAN, "--query", qf, "--db", DB, "--mode", "exhaustive",
                        "--n-iter", "4", "--threads", "8", "--out", out],
                       check=True, capture_output=True)
        d = annotate(pd.read_csv(out, sep="\t"))
        f = apply_rule(d)
        g = dedup_genes(f)
        counts.append(len(g))
        genes.append(set(g.gene))
        print(f"  shuffle {i}: {len(d):6d} hits -> rule R: {len(f):4d} rows, "
              f"{len(g):3d} genes")
    return counts, genes


def run_real():
    d = annotate(pd.read_csv("reports/exhaustive_raw.tsv", sep="\t"))
    f = apply_rule(d)
    g = dedup_genes(f)
    print(f"real exhaustive: {len(d)} hits -> rule R: {len(f)} rows, {len(g)} genes")
    hit = sorted(BENCH & set(g.gene))
    missing = sorted(BENCH - set(g.gene))
    print(f"  benchmarks recalled: {len(hit)}/4  {hit}")
    if missing:
        print(f"  benchmarks MISSED:   {missing}")
        for name in missing:
            s = d[d.gene == name]
            if len(s):
                r = s.sort_values("gscore", ascending=False).iloc[0]
                print(f"    {name}: len={r.aln_len:.0f} ident={r.ident:.3f} "
                      f"gar_q={r.gar_q:.2f} ent_t={r.ent_t:.2f} ga_t={r.ga_t:.2f} "
                      f"mf_t={r.mf_t:.2f} q={r.qstart:.0f}-{r.qend:.0f}")
    g.to_csv("reports/candidates_rule_R.csv", index=False)
    print("\n  top 40 candidate genes under R:")
    cols = ["gene", "acc", "gscore", "aln_len", "n_id", "ident", "qstart", "qend",
            "tstart", "tend", "gar_q", "ent_t"]
    print(g.head(40)[cols].to_string(index=False))
    return g


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "real"
    if mode == "real":
        run_real()
    else:
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 5
        print(f"null: {n} window-shuffled (w=20) EBNA1 queries, exhaustive n_iter=4")
        counts, genes = run_null(n)
        c = np.array(counts)
        print(f"\nnull genes under R: mean {c.mean():.1f}  sd {c.std(ddof=1):.2f}  "
              f"range [{c.min()},{c.max()}]")
        allg = set().union(*genes)
        print(f"null union of genes over {n} shuffles: {len(allg)}")
        np.save("work/mine/null_counts.npy", c)
