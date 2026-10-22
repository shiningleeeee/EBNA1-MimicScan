"""s72 -- score every screening axis against the published positives.

The whole project reduces to one methodological claim:

    A sequence-similarity axis may be used to rank candidates only if it has
    survived an examination by the experimentally established cases.

Every axis built in earlier rounds was checked that way and every statistical
one failed. That check lived in prose. This stage turns it into a table the
pipeline produces, so the verdicts are computed rather than asserted.

Declared beforehand (config `validation:`):
    min_positives_passing   how many published pairs an axis must keep
    max_pool_retained       how much of the pool it may keep and still count
                            as discriminating

Three kinds of axis, and the difference is the point:
    hard_constraint   mechanistic necessity, not a statistic. A match inside a
                      propeptide that is cleaved away cannot be an antibody
                      target whatever its score. Fails only if it removes a
                      published positive.
    ranking_axis      statistical selection. Must keep the positives AND be
                      selective, or it may not order candidates.
    annotation        displayed beside a candidate, never used to filter.

The published set deliberately retains the two pairs a linear screen cannot see
(GlialCAM: epitope defined by crystal structure; MBP: conformational), so no
axis is credited for recovering only the easy cases.
"""
from __future__ import annotations

import pandas as pd

from .. import mining, structure
from ..core import rel, write_json
from .s64_rank import _build
from .s66_maturation import _classify, _parse


def _structure_call(acc, start, end, cache, params):
    """-> (plddt, rsasa, call, folded_exposed). A missing model is not a failure."""
    empty = ("", "", "not_modelled", False)
    try:
        atom = structure.parse_pdb(structure.af_pdb(acc, cache))
    except Exception:                                   # noqa: BLE001
        return empty
    stats = structure.residue_stats(atom, structure.sasa_per_atom(atom))
    span = structure.segment(atom, stats, int(start), int(end))
    if not span.get("n"):
        return empty
    plddt, rsasa = span["plddt_mean"], span["rsasa_mean"]
    folded = (plddt >= params["folded_plddt"] and rsasa >= params["exposed_rsasa"])
    return round(plddt, 1), round(rsasa, 3), structure.classify(span), folded


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "72_axis_validation.log"
    params = dict(
        folded_plddt=float(cfg.get("structure.folded_plddt", 70)),
        exposed_rsasa=float(cfg.get("structure.exposed_rsasa", 0.25)),
    )
    cache = cfg.path_of("structure.cache_dir")
    query_acc = str(cfg.get("structure.query_accession", "P03211"))

    # ---- published positives, run through the SAME axis builder as the pool
    pairing = pd.read_csv(ctx.out("55_published_pair_null.tsv"), sep="\t")
    proteome = mining.load_fasta(cfg.path_of("mining.query_db")
                                 if cfg.get("mining.query_db") else
                                 cfg.path_of("data.proteome"))
    query_seq = mining.load_fasta(cfg.path_of("data.query"))[query_acc]
    annotation = mining.load_annotation(cfg)
    pos_scored = _build("real", pairing, query_seq, proteome, annotation, cfg)
    # `gene` comes from the pairing row via _build -- do NOT reassign it
    # positionally: _build sorts by score, so a positional list would be shifted.

    features = _parse(mining.ensure_processing_tsv(cfg, log))
    published = cfg.get("published") or []
    meta_by_acc = {str(e["acc"]): e for e in published}

    extra = {}
    for row in pos_scored.to_dict("records"):
        entry = meta_by_acc.get(str(row["acc"]), {})
        acc = str(row["acc"])
        t_start, t_end = entry.get("region", [row["t_start"], row["t_end"]])
        q_start, q_end = entry.get("ebna1_region", [row["ebna1_start"],
                                                    row["ebna1_end"]])
        meta = features.get(acc)
        call = ("no_uniprot_entry" if meta is None
                else _classify(meta["features"], int(t_start), int(t_end))[0])
        q_plddt, q_rsasa, q_call, q_folded = _structure_call(
            query_acc, q_start, q_end, cache, params)
        t_plddt, t_rsasa, t_call, t_folded = _structure_call(
            acc, t_start, t_end, cache, params)
        extra[acc] = dict(
            retained=call.startswith("retained"), processing_call=call,
            query_folded_exposed=bool(q_folded),
            target_folded_exposed=bool(t_folded),
            both_folded_exposed=bool(t_folded and q_folded),
            comp_accessible=bool(annotation.get(acc, {}).get("surface")),
            reachable=bool(entry.get("reachable", False)),
            query_struct=q_call, target_struct=t_call,
            query_plddt=q_plddt, target_plddt=t_plddt,
        )
    pos = pos_scored.copy()
    for key in ("retained", "processing_call", "query_folded_exposed",
                "target_folded_exposed", "both_folded_exposed",
                "comp_accessible", "reachable", "query_struct",
                "target_struct", "query_plddt", "target_plddt"):
        pos[key] = [extra.get(str(a), {}).get(key) for a in pos["acc"]]
    keep = ["gene", "acc", "ebna1_start", "ebna1_end", "t_start", "t_end",
            "aln_len", "identity", "margin_z", "score", "processing_call",
            "retained", "query_plddt", "query_struct", "query_folded_exposed",
            "target_plddt", "target_struct", "target_folded_exposed",
            "both_folded_exposed", "comp_accessible", "reachable"]
    pos[keep].to_csv(ctx.out("72_axis_positives.tsv"), sep="\t", index=False)

    # ---- the pool, axis by axis
    scored = pd.read_csv(ctx.out("48_candidates_scored.tsv"), sep="\t")
    layered = pd.read_csv(ctx.out("53_layered_table.tsv"), sep="\t", dtype=str)
    maturated = pd.read_csv(ctx.out("51_processing_state.tsv"), sep="\t", dtype=str)

    z_final = float(cfg.get("mining.final.min_margin_z", 3.0))
    z_iedb = float(cfg.get("mining.iedb.min_margin_z", 6.0))
    min_score = float(cfg.get("mining.composite.min_score", 6))
    min_pos = int(cfg.get("validation.min_positives_passing", 3))
    max_keep = float(cfg.get("validation.max_pool_retained", 0.25))
    min_len = float(cfg.get("mining.final.min_aln_len", 20))
    min_ident = float(cfg.get("mining.final.min_identity", 40))

    def frac(series):
        return round(float(series.mean()), 4) if len(series) else None

    axes = [
        dict(axis="maturation", kind="hard_constraint",
             rule="matched span must survive proteolytic maturation",
             hits=pos["retained"].astype(bool),
             pool=frac(maturated["retained"] == "yes"), base="maturated (363)"),
        dict(axis="min_alignment_length", kind="ranking_axis",
             rule=f"aligned length >= {min_len:g} aa",
             hits=pos["aln_len"] >= min_len,
             pool=frac(scored["aln_len"] >= min_len), base="scored (363)"),
        dict(axis="min_identity", kind="ranking_axis",
             rule=f"percent identity >= {min_ident:g}",
             hits=pos["identity"] >= min_ident,
             pool=frac(scored["identity"] >= min_ident), base="scored (363)"),
        dict(axis="margin_z_final", kind="ranking_axis",
             rule=f"composition-null margin z >= {z_final:g}",
             hits=pos["margin_z"] >= z_final,
             pool=frac(scored["margin_z"] >= z_final), base="scored (363)"),
        dict(axis="margin_z_iedb", kind="ranking_axis",
             rule=f"composition-null margin z >= {z_iedb:g}",
             hits=pos["margin_z"] >= z_iedb,
             pool=frac(scored["margin_z"] >= z_iedb), base="scored (363)"),
        dict(axis="target_normal", kind="ranking_axis",
             rule="target window entropy >= ent_min and max residue <= maxfrac",
             hits=pos["ax_target_normal"].astype(bool),
             pool=frac(scored["ax_target_normal"].astype(bool)),
             base="scored (363)"),
        dict(axis="query_folded_exposed", kind="ranking_axis",
             rule="EBNA1-side span folded and exposed in the AlphaFold model",
             hits=pos["query_folded_exposed"].astype(bool),
             pool=frac(layered["q_call"] == "surface_of_folded_domain"),
             base="layered (55)"),
        dict(axis="both_folded_exposed", kind="ranking_axis",
             rule="both spans folded and exposed",
             hits=pos["both_folded_exposed"].astype(bool),
             pool=frac(layered["both_folded_exposed"] == "yes"),
             base="layered (55)"),
        dict(axis="comp_accessible", kind="annotation",
             rule="target annotated secreted / membrane / extracellular",
             hits=pos["comp_accessible"].astype(bool),
             pool=frac(layered["comp_accessible"] == "yes"),
             base="layered (55)"),
        dict(axis="composite_score", kind="ranking_axis",
             rule=f"declared eight-axis composite score >= {min_score:g}",
             hits=pos["score"] >= min_score,
             pool=frac(scored["score"] >= min_score), base="scored (363)"),
    ]

    rows = []
    for spec in axes:
        hits = spec["hits"]
        n_pass = int(hits.sum())
        lost = ", ".join(pos.loc[~hits, "gene"].tolist()) or "-"
        kept = spec["pool"]
        kind = spec["kind"]
        if kind == "hard_constraint":
            verdict = ("pass" if n_pass == len(pos)
                       else "fail: removes a published positive")
        elif kind == "annotation":
            verdict = "annotation_only"
        elif n_pass < min_pos:
            verdict = "fail: loses positives"
        elif kept is not None and kept > max_keep:
            verdict = "fail: not selective"
        else:
            verdict = "pass"
        consequence = {
            "pass": "may order candidates" if kind == "ranking_axis"
                    else "applied as a hard filter",
            "annotation_only": "may be displayed, never used to filter or order",
        }.get(verdict, "may not be used to filter or order candidates")
        rows.append(dict(
            axis=spec["axis"], kind=kind, rule=spec["rule"],
            positives_passing=n_pass, positives_total=int(len(pos)),
            positives_lost=lost, pool_base=spec["base"],
            pool_retained_fraction=kept, verdict=verdict,
            consequence=consequence,
        ))

    table = pd.DataFrame(rows)
    table.to_csv(ctx.out("72_axis_validation.tsv"), sep="\t", index=False)

    n_rank = int((table["kind"] == "ranking_axis").sum())
    n_rank_pass = int(((table["kind"] == "ranking_axis")
                       & (table["verdict"] == "pass")).sum())
    payload = dict(
        declared_before_measuring=dict(
            min_positives_passing=min_pos,
            max_pool_retained=max_keep,
            positives=[e["gene"] for e in published],
            positives_a_linear_screen_cannot_see=[
                e["gene"] for e in published if not e.get("reachable", False)],
        ),
        n_axes=int(len(table)),
        n_ranking_axes=n_rank,
        n_ranking_axes_passed=n_rank_pass,
        n_hard_constraints=int((table["kind"] == "hard_constraint").sum()),
        verdicts={k: int(v) for k, v in table["verdict"].value_counts().items()},
        note=("A ranking axis keeping fewer than "
              f"{min_pos} of {len(pos)} published pairs, or keeping more than "
              f"{max_keep:.0%} of the pool, cannot discriminate and is "
              "forbidden from ordering candidates. Positive values are computed "
              "by the same axis builder used for the pool; no threshold was "
              "tuned after seeing them."),
    )
    write_json(cfg.reports / "72_axis_validation_summary.json", payload)
    print(f"72_axis_validation: {len(table)} axes, "
          f"{n_rank_pass}/{n_rank} ranking axes pass", flush=True)
    for row in table.to_dict("records"):
        kept = ("-" if row["pool_retained_fraction"] is None
                else f"{row['pool_retained_fraction']:.3f}")
        print(f"  {row['axis']:<22} {row['kind']:<16} "
              f"{row['positives_passing']}/{row['positives_total']} "
              f"pool {kept:>6} -> {row['verdict']}", flush=True)
    return dict(
        outputs=["reports/72_axis_positives.tsv",
                 "reports/72_axis_validation.tsv",
                 "reports/72_axis_validation_summary.json"],
        facts=dict(n_axes=int(len(table)), n_ranking_axes=n_rank,
                   n_ranking_axes_passed=n_rank_pass,
                   n_positives=int(len(pos)),
                   verdicts=payload["verdicts"]),
    )
