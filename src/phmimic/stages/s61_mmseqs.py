"""s61 -- the search that involves no prior at all.

The query is the whole 641-aa EBNA1. No seed, no masked region, no coordinate
from the literature: the only thing that decides what appears is how well a
stretch of EBNA1 aligns to a human protein under mmseqs2's own statistics.

Outputs
    reports/mmseqs_ebna1.m8        raw alignments (tab-separated, headerless)
    reports/45_local_similarity.tsv the same alignments with annotation
    reports/cands_for_null.tsv      the coordinates the pair-null stage consumes

Why the raw file has no header: it is the wire format of `mmseqs convertalis`,
kept verbatim so it can be diffed against a run made outside the pipeline. The
column order is fixed and is recorded in `mining.M8_COLUMNS`.
"""
from __future__ import annotations

import pandas as pd

from .. import mining
from ..core import write_json


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "61_mmseqs.log"
    query_fa = cfg.path_of("mining.query")

    m8_path = ctx.out("mmseqs_ebna1.m8")
    mining.mmseqs_search(cfg, query_fa, m8_path, "work/mmseqs", log, name="real")
    hits = mining.read_m8(m8_path)
    mining.write_m8(hits, m8_path)

    # The pair-null stage wants only the coordinates it must re-score; writing
    # them from the sorted table (not from the search) keeps the two in step.
    cands = hits[["target", "qstart", "qend", "tstart", "tend", "alnlen", "pident"]]
    cands.to_csv(ctx.out("cands_for_null.tsv"), sep="\t", index=False)

    # ---- annotated view: two orderings, because they answer two questions ----
    ann = {}
    ann_path = cfg.path_of("data.uniprot_annot")
    if ann_path.exists():
        ann = pd.read_csv(ann_path, sep="\t", dtype=str).fillna("").set_index("Entry")

    ga_block = cfg.get("composition.gar_block", [90, 330])
    rows = []
    for r in hits.itertuples(index=False):
        acc = str(r.target)
        gene = protein = loc = tissue = kw = ""
        if acc in ann.index:
            raw_gene = str(ann.loc[acc, "Gene Names"])
            gene = "" if raw_gene in ("nan", "") else raw_gene.split(";")[0].split()[0]
            protein = str(ann.loc[acc, "Protein names"])
            loc = str(ann.loc[acc, "Subcellular location [CC]"])
            tissue = str(ann.loc[acc, "Tissue specificity"])
            kw = str(ann.loc[acc, "Keywords"])
        rows.append(dict(
            acc=acc, gene=gene, protein=protein[:70],
            ebna1_start=int(r.qstart), ebna1_end=int(r.qend),
            region=mining.region_of(int(r.qstart), int(r.qend), tuple(ga_block)),
            t_start=int(r.tstart), t_end=int(r.tend),
            aln_len=int(r.alnlen), identity=round(float(r.pident), 1),
            evalue=float(r.evalue), bits=int(r.bits),
            cns=mining.flag(" ".join((tissue, kw)), mining.CNS_TOKENS),
            surface=mining.flag(" ".join((loc, kw)), mining.SURFACE_TOKENS),
            qaln=r.qaln, taln=r.taln,
        ))
    frame = pd.DataFrame(rows)
    order = [c for c in ("identity", "aln_len", "acc", "ebna1_start") if c in frame.columns]
    frame = frame.sort_values(order, ascending=[False, False, True, True],
                              kind="mergesort").reset_index(drop=True)
    frame.insert(0, "rank", range(1, len(frame) + 1))
    frame.to_csv(ctx.out("45_local_similarity.tsv"), sep="\t", index=False)

    best = frame.sort_values(["evalue", "acc"], kind="mergesort").iloc[0]
    regions = frame["region"].value_counts().to_dict()
    write_json(ctx.out("61_mmseqs_summary.json"), dict(
        n_alignments=int(len(frame)),
        n_targets=int(frame["acc"].nunique()),
        best_evalue=float(best["evalue"]),
        best_acc=str(best["acc"]),
        best_gene=str(best["gene"]),
        best_aln_len=int(best["aln_len"]),
        best_identity=float(best["identity"]),
        n_evalue_below_1=int((frame["evalue"] < 1.0).sum()),
        region_counts={str(k): int(v) for k, v in regions.items()},
        ga_block_overlap=int(sum(
            1 for r in frame.itertuples(index=False)
            if mining.overlaps(int(r.ebna1_start), int(r.ebna1_end), *ga_block)
        )),
        mmseqs=dict(
            sensitivity=cfg.get("mmseqs.sensitivity"),
            max_seqs=cfg.get("mmseqs.max_seqs"),
            e_cutoff=cfg.get("mmseqs.e_cutoff"),
            comp_bias_corr=cfg.get("mmseqs.comp_bias_corr"),
        ),
        note=("mmseqs2 default low-complexity masking and --comp-bias-corr are "
              "both active. Their contributions are separated in "
              "src/analysis/60_mmseqs_param_matrix.py: with --comp-bias-corr 1 "
              "the Gly-Ala block collects 0 overlapping alignments whether or "
              "not tantan masking is on (363 vs 431 alignments overall), while "
              "with --comp-bias-corr 0 it collects 1254 (mask 1) / 1696 "
              "(mask 0) overlapping hits and each such cell fails its own "
              "windowed null (p=0.995 / 0.881). The zero overlap is therefore "
              "attributable to composition correction, not to masking."),
    ))

    ctx.log(f"    {len(frame)} alignments over {frame['acc'].nunique()} proteins; "
            f"best E={float(best['evalue']):.3g} ({best['gene']}); "
            f"E<1: {int((frame['evalue'] < 1.0).sum())}")
    return dict(
        outputs=["reports/mmseqs_ebna1.m8", "reports/45_local_similarity.tsv",
                 "reports/cands_for_null.tsv", "reports/61_mmseqs_summary.json"],
        facts=dict(n_alignments=int(len(frame)),
                   n_targets=int(frame["acc"].nunique()),
                   best_evalue=float(best["evalue"]),
                   n_evalue_below_1=int((frame["evalue"] < 1.0).sum())),
    )
