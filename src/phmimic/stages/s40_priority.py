"""s40 -- score every seed and calibrate the real ones against their nulls.

This is the stage that decides whether anything here is a result. A real seed's
best reachable candidate is only interesting if the peptides that share its
composition -- or arbitrary peptides of the same length -- cannot reach the same
score. If they can, the ranking below is a description of the human proteome's
composition, not evidence of mimicry.
"""
from __future__ import annotations

import pandas as pd

from ..annotate import ensure_annotation
from ..core import write_json
from ..score import empirical_p, score_seed


def _tuples(sub: pd.DataFrame):
    return list(zip(sub["acc"], sub["tpos"].astype(int),
                    sub["n_ident"].astype(int), sub["matched"]))


def run(ctx):
    seeds = pd.read_csv(ctx.out("10_seeds.tsv"), sep="\t", dtype=str).fillna("")
    matches = pd.read_csv(ctx.out("20_motif_matches.tsv"), sep="\t",
                          dtype={"acc": str, "motif": str, "matched": str,
                                 "seed_id": str, "family": str, "parent": str})
    annotation = ensure_annotation(ctx.cfg)

    by_seed = {sid: sub for sid, sub in matches.groupby("seed_id")}
    scored: dict[str, dict] = {}
    for seed in seeds.itertuples(index=False):
        sub = by_seed.get(seed.seed_id)
        result = score_seed(_tuples(sub) if sub is not None else [], annotation, ctx.cfg)
        result.update(family=seed.family, parent=seed.parent, motif=seed.motif,
                      length=len(seed.motif), max_mismatch=int(seed.max_mismatch))
        scored[seed.seed_id] = result

    seed_rows = []
    for sid, r in scored.items():
        seed_rows.append(dict(
            seed_id=sid, family=r["family"], parent=r["parent"], motif=r["motif"],
            length=r["length"], max_mismatch=r["max_mismatch"],
            n_matches=r["n_matches"], n_proteins=r["n_proteins"],
            best_score=r["best_score"], best_acc=r["best_acc"], best_gene=r["best_gene"],
            best_ident=r["best_ident"],
        ))
    seed_frame = pd.DataFrame(seed_rows)
    seed_frame.to_csv(ctx.out("40_seed_scores.tsv"), sep="\t", index=False)

    nb_ids = [s for s in scored if scored[s]["family"] == "N_B"]
    nb_len = {scored[s]["length"] for s in nb_ids}
    calibration = []
    for sid, r in scored.items():
        if r["family"] != "real":
            continue
        na_ids = [s for s in scored
                  if scored[s]["family"] == "N_A" and scored[s]["parent"] == sid]
        na_scores = [scored[s]["best_score"] for s in na_ids]
        p_na = empirical_p(r["best_score"], na_scores)

        p_nb = None
        if nb_ids and r["length"] in nb_len:
            nb_scores = [scored[s]["best_score"] for s in nb_ids]
            p_nb = empirical_p(r["best_score"], nb_scores)

        calibration.append(dict(
            seed_id=sid, motif=r["motif"], length=r["length"],
            best_score=r["best_score"], best_gene=r["best_gene"], best_acc=r["best_acc"],
            n_proteins=r["n_proteins"],
            null_A=p_na,
            null_B=p_nb,
            null_B_note="" if p_nb else "N_B seeds differ in length; not comparable",
        ))
        ctx.log(f"    [{sid}] {r['motif']} best={r['best_score']} ({r['best_gene'] or '-'}) | "
                f"N_A n={p_na['n_null']} max={p_na['null_max']} "
                f"p_uncond={p_na['p_unconditional']} p_cond={p_na['p_conditional']}"
                + (f" | N_B p_cond={p_nb['p_conditional']}" if p_nb else ""))

    write_json(ctx.out("40_null_calibration.json"),
               dict(calibration=calibration,
                    note="null_A permutes each real seed's own residues; "
                         "null_B is unrelated random 6-mers of a different length "
                         "for C2/C3 and is therefore reported but not used there"))

    cand_frames = []
    for sid, r in scored.items():
        if r["family"] != "real" or "candidates" not in r:
            continue
        frame = r["candidates"].copy()
        frame.insert(0, "seed_id", sid)
        frame.insert(1, "seed_motif", r["motif"])
        cand_frames.append(frame)
    if cand_frames:
        cands = pd.concat(cand_frames, ignore_index=True)
        cands = cands.sort_values(["score", "cns", "surface", "ident", "acc"],
                                  ascending=[False, False, False, False, True])
    else:
        cands = pd.DataFrame()
    if len(cands):
        published = {str(p["acc"]): str(p.get("gene", ""))
                     for p in (ctx.cfg.get("published") or [])}
        cands["published"] = [published.get(str(a), "") for a in cands["acc"]]
    cands.to_csv(ctx.out("40_candidates_ranked.tsv"), sep="\t", index=False)

    top = cands.head(5)
    for row in top.itertuples(index=False):
        ctx.log(f"    {row.seed_id} {row.acc} {row.gene:<9} score={row.score} "
                f"ident={row.ident} cns={row.cns} surf={row.surface}")
    return dict(
        outputs=["reports/40_seed_scores.tsv", "reports/40_candidates_ranked.tsv",
                 "reports/40_null_calibration.json"],
        facts=dict(n_seeds_scored=len(scored), n_candidates=int(len(cands)),
                   calibration=calibration),
    )
