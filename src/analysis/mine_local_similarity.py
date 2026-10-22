"""Unbiased local-similarity mining: whole EBNA1 vs the human proteome.

No epitope prior. The query is the entire 641-aa EBNA1 and the only thing that
decides what appears is how well a stretch of it aligns to a human protein.

Two views are reported because they answer two different questions:

  * by E-value -- "what would a standard homology search report?", which
    favours long alignments whatever their identity
  * by identity -- "where is there a short, near-identical stretch?", which is
    closer to what an antibody actually sees

Composition is not filtered here; it is reported, so the reader can see which
alignments are driven by the Gly/Ala block and which are not.
"""
from __future__ import annotations

import pandas as pd

COLUMNS = ["target", "evalue", "bits", "qstart", "qend", "tstart", "tend",
           "alnlen", "pident", "qaln", "taln"]

SURFACE_TOKENS = ("cell membrane", "secreted", "cell surface", "extracellular",
                  "gpi-anchor", "basement membrane")
CNS_TOKENS = ("brain", "cerebellum", "cerebral cortex", "hippocampus",
              "spinal cord", "neuron", "neuro", "synap", "myelin", "glia")


def _region(qstart, qend):
    """Which part of EBNA1 the alignment sits on -- descriptive, not a filter."""
    if qend < 90:
        return "N-terminal"
    if qstart <= 330 and qend >= 90:
        return "Gly/Ala block (90-330)"
    if 370 <= qstart <= 445:
        return "SPR-rich (370-445)"
    if 476 <= qstart <= 585:
        return "folded core (476-585)"
    return "C-terminal"


def _flag(text, tokens):
    low = str(text).lower()
    return any(t in low for t in tokens)


def main():
    m = pd.read_csv("reports/mmseqs_ebna1.m8", sep="\t", names=COLUMNS)
    ann = pd.read_csv("data/cache/uniprot_human_annot.tsv", sep="\t",
                      dtype=str).set_index("Entry")

    rows = []
    for r in m.itertuples(index=False):
        acc = r.target
        gene = ""
        name = ""
        loc = ""
        tissue = ""
        kw = ""
        if acc in ann.index:
            gene = str(ann.loc[acc, "Gene Names"]).split(";")[0].split()[0] \
                if str(ann.loc[acc, "Gene Names"]) not in ("nan", "") else ""
            name = str(ann.loc[acc, "Protein names"])
            loc = str(ann.loc[acc, "Subcellular location [CC]"])
            tissue = str(ann.loc[acc, "Tissue specificity"])
            kw = str(ann.loc[acc, "Keywords"])
        rows.append(dict(
            acc=acc, gene=gene, protein=name[:70],
            ebna1_start=int(r.qstart), ebna1_end=int(r.qend),
            region=_region(int(r.qstart), int(r.qend)),
            t_start=int(r.tstart), t_end=int(r.tend),
            aln_len=int(r.alnlen), identity=round(float(r.pident), 1),
            evalue=float(r.evalue), bits=int(r.bits),
            cns=_flag(" ".join((tissue, kw)), CNS_TOKENS),
            surface=_flag(" ".join((loc, kw)), SURFACE_TOKENS),
            qaln=r.qaln, taln=r.taln,
        ))
    frame = pd.DataFrame(rows)
    frame = frame.sort_values(["identity", "aln_len"], ascending=[False, False])
    frame.insert(0, "rank", range(1, len(frame) + 1))
    frame.to_csv("reports/45_local_similarity.tsv", sep="\t", index=False)

    print(f"{len(frame)} alignments over {frame.acc.nunique()} human proteins\n")
    print("--- view A: highest local identity (what an antibody could see) ---")
    for r in frame.head(15).itertuples(index=False):
        print(f"{r.rank:>3} {r.identity:5.1f}% {r.aln_len:>3}aa E={r.evalue:9.3g} "
              f"EBNA1 {r.ebna1_start}-{r.ebna1_end} -> {r.gene:<10} {r.acc}  "
              f"[{r.region}] {'CNS ' if r.cns else ''}{'surface' if r.surface else ''}")
    print("\n--- view B: lowest E-value (what a homology search reports) ---")
    eb = frame.sort_values("evalue").head(10)
    for r in eb.itertuples(index=False):
        print(f"    {r.identity:5.1f}% {r.aln_len:>3}aa E={r.evalue:9.3g} "
              f"EBNA1 {r.ebna1_start}-{r.ebna1_end} -> {r.gene:<10} {r.acc}  [{r.region}]")
    print("\nregion distribution:")
    print(frame.region.value_counts().to_string())
    print(f"\nidentity>=40%: {int((frame.identity>=40).sum())} alignments, "
          f"{frame[frame.identity>=40].acc.nunique()} proteins")
    print(f"aln_len>=20 and identity>=40%: {int(((frame.aln_len>=20)&(frame.identity>=40)).sum())}")


if __name__ == "__main__":
    main()
