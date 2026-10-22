#!/usr/bin/env python3
"""12_candidate_cores.py -- seed candidates on the EXPERIMENTALLY VALIDATED
EBNA1 cross-reactive epitopes, not on genome-wide similarity.

Rationale
---------
Step 2 and Step 3 established that genome-wide similarity ranking is a
composition artefact: aggregate short-motif sharing, MMseqs2 E-values and
Biomni's E_emp all fail to separate EBNA1 from composition-matched shuffles.
What *is* discriminating is the opposite direction: start from the four
epitopes that antibodies/sera actually validated, and ask which human
proteins carry that exact motif.

Seeds (EBNA1 coordinates are 1-based, inclusive):
  C1  PPRRPP        394-399   MEF2B / GlialCAM (Lum 2026 PNAS)
  C2  RRPFF         402-406   CRYAB          (Thomas 2023 Sci Adv)
  C3  PGAIEQGPAD    431-440   ANO2           (Tengvall 2019 PNAS)
  C4  (20-mer)      386-405   HEPACAM        (Lanz 2022 Nature) - soft

Mismatch budget: 6-mer d<=1, 5-mer d<=1, 10-mer d<=3 (the published ANO2
alignment is 7/10, so d=3 is required or ANO2 itself is lost). C4 is a
composition-level assignment, so it is d<=6 and flagged soft.

Usage: env/venv/bin/python src/12_candidate_cores.py
Outputs: reports/candidates_cores.tsv, prints a ranked table.
"""
import json
import sys
from collections import defaultdict

import numpy as np
import requests

DB = "data/cache/human_swissprot.fasta"
QUERY = "data/cache/P03211.fasta"
UA = {"User-Agent": "ph-mimic/0.1 (candidate scoring)"}

CORES = [
    ("C1_PPRRPP",     "PPRRPP",             394, 399, 1, "hard"),
    ("C2_RRPFF",      "RRPFF",              402, 406, 1, "hard"),
    ("C3_PGAIEQGPAD", "PGAIEQGPAD",         431, 440, 3, "hard"),
    ("C4_HEPACAM20",  "SQSSSSGSPPRRPPPGRRPF", 386, 405, 6, "soft"),
]


def read_fasta(path):
    """Return list of (acc, gene, seq) -- UniProt-style sp|ACC|NAME headers."""
    out, hdr, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(buf)))
            hdr = line[1:]
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    recs = []
    for h, s in out:
        acc = h.split("|")[1] if "|" in h else h.split()[0]
        gene = ""
        for tok in h.split():
            if tok.startswith("GN="):
                gene = tok[3:]
                break
        recs.append((acc, gene, s.upper()))
    return recs


def scan_core(mat, offsets, core, d):
    """Vectorised Hamming scan of core against the concatenated proteome."""
    L = len(core)
    if len(mat) < L:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    win = np.lib.stride_tricks.sliding_window_view(mat, L)
    ref = np.frombuffer(core.encode(), dtype=np.uint8)
    mm = (win != ref).sum(axis=1)
    idx = np.flatnonzero(mm <= d)
    return idx, mm[idx].astype(np.int32)


def main():
    ebna = {a: s for a, g, s in read_fasta(QUERY)}["P03211"]
    recs = read_fasta(DB)
    accs = [r[0] for r in recs]
    genes = [r[1] for r in recs]
    seqs = [r[2] for r in recs]

    # concatenate with a sentinel that cannot match any residue
    blobs = []
    offsets = np.zeros(len(seqs) + 1, dtype=np.int64)
    for i, s in enumerate(seqs):
        offsets[i + 1] = offsets[i] + len(s) + 1
        blobs.append(s)
        blobs.append("*")
    mat = np.frombuffer("".join(blobs).encode(), dtype=np.uint8)
    print(f"[db] {len(recs)} proteins, {len(mat):,} residues (with sentinels)")

    # per-protein best match to each core
    best = defaultdict(dict)          # acc -> core_name -> (ident, tpos, tseq)
    carry = {}                        # core_name -> number of proteins matched
    for name, core, qs, qe, d, kind in CORES:
        idx, mm = scan_core(mat, offsets, core, d)
        hit_acc = np.searchsorted(offsets, idx, side="right") - 1
        keep = offsets[hit_acc + 1] - 1 - idx >= len(core)   # not spanning sentinel
        idx, mm, hit_acc = idx[keep], mm[keep], hit_acc[keep]
        prot_hit = set()
        for i, m, a in zip(idx, mm, hit_acc):
            acc = accs[a]
            pos = int(i - offsets[a]) + 1
            ident = len(core) - int(m)
            ts = "".join(chr(c) for c in mat[i:i + len(core)])
            prot_hit.add(acc)
            prev = best[acc].get(name)
            if prev is None or ident > prev[0]:
                best[acc][name] = (ident, pos, ts, len(core), kind)
        carry[name] = len(prot_hit)
        print(f"[{name}] {core} EBNA1 {qs}-{qe} d<={d} -> "
              f"{len(idx):,} matches in {len(prot_hit):,} proteins")

    # score: how many HARD cores are matched at >=80% identity
    rows = []
    for acc, hits in best.items():
        hard = [(n, v) for n, v in hits.items() if v[4] == "hard"]
        n80 = sum(1 for n, v in hard if v[0] / v[3] >= 0.80)
        ids = [v[0] / v[3] for n, v in hard]
        soft = sum(1 for n, v in hits.items() if v[4] == "soft")
        rows.append(dict(
            acc=acc,
            gene=genes[accs.index(acc)],
            prot_len=len(seqs[accs.index(acc)]),
            n_hard_cores=len(hard),
            n_hard_cores_80=n80,
            best_ident=round(max(ids), 3) if ids else 0.0,
            soft_hepacam=soft,
            cores=",".join(sorted(hits)),
        ))

    rows.sort(key=lambda r: (-r["n_hard_cores_80"], -r["best_ident"],
                             -r["n_hard_cores"], -r["soft_hepacam"], r["gene"]))
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv("reports/candidates_cores.tsv", sep="\t", index=False)
    print(f"\n[out] reports/candidates_cores.tsv  ({len(df)} proteins)")
    print(f"[rarity] proteins per core: {carry}")

    # detail for the top slice
    top = df.head(40)
    print(f"\n{'acc':<8} {'gene':<11} {'len':>5} {'hard':>4} {'>=80%':>5} "
          f"{'best_id':>7} {'soft':>4}  cores")
    for _, r in top.iterrows():
        print(f"{r['acc']:<8} {r['gene']:<11} {r['prot_len']:>5} "
              f"{r['n_hard_cores']:>4} {r['n_hard_cores_80']:>5} "
              f"{r['best_ident']:>7.3f} {r['soft_hepacam']:>4}  {r['cores']}")

    # per-candidate motif detail for the shortlist
    short = top[top["n_hard_cores_80"] >= 1]["acc"].tolist()[:30]
    print(f"\n[motif detail] {len(short)} proteins")
    for acc in short:
        for name, v in sorted(best[acc].items()):
            ident, pos, ts, ln, kind = v
            print(f"  {acc} {name:<15} {ts:<21} t@{pos:<6} {ident}/{ln} {kind}")

    json.dump({a: {n: v for n, v in d.items()} for a, d in best.items()
               if a in short}, open("reports/candidates_cores_detail.json", "w"), indent=1)


if __name__ == "__main__":
    main()
