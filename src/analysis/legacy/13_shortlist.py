#!/usr/bin/env python3
"""13_shortlist.py -- specificity-first shortlist.

Why this filter and not the previous ones
-----------------------------------------
Step 2 (significance) and Step 3 (aggregate short-motif null) both failed for the
same reason: whole-proteome reasoning is dominated by composition. The fix is to
invert the direction. The only thing in this problem that is experimentally
anchored is the EBNA1 cross-reactive region 386-440, where antibodies/sera
validated four epitopes. So:

  seed  = exact (d=0) match to a validated 5/6-mer core, or d<=3 to the 10-mer
  rank  = rarity of that motif in the human proteome, then orthogonal biology
          (CNS expression / surface accessibility / known autoantigen), which
          NO composition statistic can supply

A one-mismatch budget is deliberately NOT used for the short cores: d<=1 turns
PPRRPP from 6 proteins into 282 and RRPFF from 3 into 289 (see
reports/candidates_cores.tsv). That is the same artefact regime as Step 2, so
the mismatch budget is spent only where the core is long enough to afford it
(the 10-mer, where the published ANO2 alignment is itself 7/10).

Usage: env/venv/bin/python src/13_shortlist.py
Writes: reports/candidates_shortlist.tsv
"""
import json
import time

import numpy as np
import pandas as pd
import requests

DB = "data/cache/human_swissprot.fasta"
QUERY = "data/cache/P03211.fasta"
UA = {"User-Agent": "ph-mimic/0.1 (candidate shortlist)"}

CORES = [
    ("C1", "PPRRPP",     394, 399, 0),
    ("C2", "RRPFF",      402, 406, 0),
    ("C3", "PGAIEQGPAD", 431, 440, 3),
]

# biology that no sequence statistic can produce; filled from UniProt + literature
CNS_KW = ("brain", "neuron", "neuronal", "synap", "axon", "glial", "myelin",
          "cerebell", "cerebral", "spinal", "central nervous", "retina")
SURF_KW = ("Cell membrane", "Secreted", "Extracellular", "Cell surface",
           "Single-pass type I membrane protein", "GPI-anchor")


def read_fasta(path):
    out, hdr, buf = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(buf)))
            hdr, buf = line[1:], []
        else:
            buf.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    recs = []
    for h, s in out:
        acc = h.split("|")[1] if "|" in h else h.split()[0]
        gene = next((t[3:] for t in h.split() if t.startswith("GN=")), "")
        desc = h.split(" OS=")[0]
        recs.append((acc, gene, s.upper(), desc))
    return recs


def scan(mat, offsets, core, d):
    L = len(core)
    win = np.lib.stride_tricks.sliding_window_view(mat, L)
    ref = np.frombuffer(core.encode(), dtype=np.uint8)
    mm = (win != ref).sum(axis=1)
    idx = np.flatnonzero(mm <= d)
    return idx, mm[idx].astype(np.int32)


def fetch_uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    for _ in range(3):
        try:
            r = requests.get(url, headers=UA, timeout=30)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(1.5)
    return None


def main():
    ebna = {a: s for a, g, s, d in read_fasta(QUERY)}["P03211"]
    recs = read_fasta(DB)
    accs = [r[0] for r in recs]
    genes = {r[0]: r[1] for r in recs}
    descs = {r[0]: r[3] for r in recs}
    lens = {r[0]: len(r[2]) for r in recs}

    blobs, offsets = [], np.zeros(len(recs) + 1, dtype=np.int64)
    for i, r in enumerate(recs):
        offsets[i + 1] = offsets[i] + len(r[2]) + 1
        blobs.append(r[2])
        blobs.append("*")
    mat = np.frombuffer("".join(blobs).encode(), dtype=np.uint8)

    hits = {}
    for name, core, qs, qe, d in CORES:
        idx, mm = scan(mat, offsets, core, d)
        acc_i = np.searchsorted(offsets, idx, side="right") - 1
        ok = offsets[acc_i + 1] - 1 - idx >= len(core)
        idx, mm, acc_i = idx[ok], mm[ok], acc_i[ok]
        for i, m, a in zip(idx, mm, acc_i):
            acc = accs[a]
            hits.setdefault(acc, []).append(dict(
                core=name, ebna1=f"{qs}-{qe}", motif_t="".join(chr(c) for c in mat[i:i + len(core)]),
                tpos=int(i - offsets[a]) + 1, n_ident=len(core) - int(m), core_len=len(core)))
        print(f"[{name}] {core} d<={d}: {len(set(accs[a] for a in acc_i))} proteins")

    rows = []
    for acc, hs in hits.items():
        u = fetch_uniprot(acc)
        tissue, subcell, disease, function, feats = "", [], "", "", []
        kw, protname = [], descs.get(acc, "")
        if u:
            protname = u.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", protname)
            for c in u.get("comments", []):
                t = c.get("commentType")
                if t == "TISSUE SPECIFICITY":
                    tissue = " ".join(x.get("value", "") for x in c.get("texts", []))
                elif t == "SUBCELLULAR LOCATION":
                    subcell += [l.get("location", {}).get("value", "") for l in c.get("subcellularLocations", [])]
                elif t == "DISEASE":
                    disease = "; ".join(x.get("disease", {}).get("diseaseId", "") for x in c.get("disease", {}).get("involvement", []) if x)
                elif t == "FUNCTION":
                    function = " ".join(x.get("value", "") for x in c.get("texts", []))[:220]
            kw = [k.get("name", "") for k in u.get("keywords", [])]
            for f in u.get("features", []):
                feats.append((f.get("type", ""), f.get("location", {}).get("start", {}).get("value"),
                              f.get("location", {}).get("end", {}).get("value")))
        blob = " ".join([tissue, function, protname] + kw) + " " + " ".join(subcell)
        cns = any(k.lower() in blob.lower() for k in CNS_KW)
        surface = any(any(s.lower() in sl.lower() for sl in subcell) for s in SURF_KW)
        # is the motif inside an annotated structured domain, or in a linker?
        in_dom = any(ty in ("Domain", "Repeat", "Zinc finger", "Immunoglobulin domain",
                            "Coiled coil") and s and e and s <= h["tpos"] <= e
                     for h in hs for ty, s, e in feats)
        disord = any(ty == "Region" and s and e and s <= h["tpos"] <= e
                     for h in hs for ty, s, e in feats)
        best = max(hs, key=lambda h: (h["n_ident"] / h["core_len"], h["core_len"]))
        rows.append(dict(
            acc=acc, gene=genes.get(acc, ""), name=protname[:60], prot_len=lens.get(acc, 0),
            cores=",".join(f'{h["core"]}@{h["tpos"]}({h["n_ident"]}/{h["core_len"]})' for h in hs),
            n_exact=sum(1 for h in hs if h["n_ident"] == h["core_len"]),
            best_frac=round(best["n_ident"] / best["core_len"], 3),
            cns="Y" if cns else "", surface="Y" if surface else "",
            rel_pos=round(best["tpos"] / max(1, lens.get(acc, 1)), 2),
            subcell="; ".join(sorted(set(x for x in subcell if x)))[:70],
            tissue=tissue[:110], disease=disease[:60],
            keywords="; ".join(kw)[:90], function=function[:150],
        ))

    df = pd.DataFrame(rows)
    df["support"] = (df["cns"].eq("Y").astype(int) + df["surface"].eq("Y").astype(int)
                     + df["disease"].ne("").astype(int))
    df = df.sort_values(["n_exact", "best_frac", "support", "cns", "surface"],
                        ascending=[False, False, False, False, False])
    df.to_csv("reports/candidates_shortlist.tsv", sep="\t", index=False)
    print(f"\n[out] reports/candidates_shortlist.tsv ({len(df)} proteins)\n")
    for _, r in df.iterrows():
        print(f"{r['acc']:<8} {r['gene']:<10} {r['cores']:<26} exact={r['n_exact']} "
              f"id={r['best_frac']:.3f} cns={r['cns'] or '-':<1} surf={r['surface'] or '-':<1} "
              f"| {r['name'][:34]:<34} | {r['subcell'][:30]:<30} | {r['tissue'][:40]}")


if __name__ == "__main__":
    main()
