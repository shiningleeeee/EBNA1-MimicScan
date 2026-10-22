"""s65 -- structural accessibility of the matched segment.

Same question as every other stage, asked in a different space: given that the
peptide matches, can an antibody physically reach it? A motif buried in a folded
core cannot be a serological mimic even when it is sequence-identical.

AlphaFold models are used as a proxy and labelled as such. pLDDT is a model
confidence, not a measurement of disorder; rSASA is computed from a predicted
monomer, so it ignores partners, membranes and the Gly/Ala block's own
flexibility. Nothing here feeds the priority score -- it either removes a
candidate from further wet-lab consideration or leaves it where it was.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import structure as st
from ..core import write_json


def _measure(acc, spans, cache, params):
    """pLDDT / rSASA over one model, for every span asked about."""
    try:
        meta = st.af_metadata(acc, cache)
        atom = st.parse_pdb(st.af_pdb(acc, cache))
    except Exception as exc:  # noqa: BLE001 - network or 404, both are "no model"
        return None, f"{type(exc).__name__}"
    sasa = st.sasa_per_atom(atom, probe=params["probe_radius"],
                            n_points=params["sphere_points"])
    stats = st.residue_stats(atom, sasa)
    out = {}
    for key, (start, end) in spans.items():
        seg = st.segment(atom, stats, start, end)
        # An AlphaFold model can be shorter than the UniProt entry (fragments,
        # removed signal peptides), so a span may fall outside it entirely.
        if seg.get("n", 0) == 0:
            out[key] = dict(n=0, plddt_mean=float("nan"), plddt_min=float("nan"),
                            rsasa_mean=float("nan"), frac_rsasa_ge_025=float("nan"),
                            call="not_modelled")
            continue
        seg["call"] = st.classify(seg)
        out[key] = seg
    return dict(meta=meta, spans=out, stats=stats), None


def run(ctx):
    params = dict(
        probe_radius=float(ctx.cfg.get("structure.probe_radius", 1.4)),
        sphere_points=int(ctx.cfg.get("structure.sphere_points", 92)),
        folded_plddt=float(ctx.cfg.get("structure.folded_plddt", 70)),
        exposed_rsasa=float(ctx.cfg.get("structure.exposed_rsasa", 0.25)),
        null_per_family=int(ctx.cfg.get("structure.null_per_family", 40)),
    )
    cache = ctx.cfg.path_of("structure.cache_dir")
    cache.mkdir(parents=True, exist_ok=True)

    ranked_path = ctx.out("40_candidates_ranked.tsv")
    if not ranked_path.exists():
        ctx.log("    no ranking; skipped")
        return dict(outputs=[], facts=dict(skipped="no ranking"))
    ranked = pd.read_csv(ranked_path, sep="\t", dtype={"acc": str})

    # ---- the query protein itself: where are the validated epitopes?
    query_acc = ctx.cfg.get("data.query_accession", "P03211")
    epitopes = {f"{s['id']}_{s['ebna1'][0]}_{s['ebna1'][1]}": tuple(s["ebna1"])
                for s in ctx.cfg.get("seeds.real")}
    epitopes["block"] = tuple(ctx.cfg.get("composition.gar_block"))
    query_measure, query_err = _measure(query_acc, epitopes, cache, params)
    query_facts = dict(acc=query_acc, error=query_err)
    if query_measure:
        query_facts.update(
            global_plddt=query_measure["meta"].get("globalMetricValue"),
            fraction_plddt_very_low=query_measure["meta"].get("fractionPlddtVeryLow"),
            segments={k: dict(plddt_mean=v["plddt_mean"], rsasa_mean=v["rsasa_mean"],
                              call=v["call"])
                      for k, v in query_measure["spans"].items()},
        )
        ctx.log(f"    query {query_acc}: global pLDDT "
                f"{query_facts['global_plddt']:.1f}, "
                f"{query_facts['fraction_plddt_very_low']:.0%} of residues very low")

    # ---- candidates
    rows = []
    excluded = []
    for row in ranked.itertuples(index=False):
        span = (int(row.pos), int(row.pos) + len(str(row.matched)) - 1)
        measured, err = _measure(row.acc, {"m": span}, cache, params)
        if measured is None:
            excluded.append(dict(acc=row.acc, gene=row.gene, family="real",
                                 reason=err))
            rows.append(dict(acc=row.acc, gene=row.gene, seed_id=row.seed_id,
                             motif=row.matched, pos=span[0], span=span[1] - span[0] + 1,
                             plddt=np.nan, plddt_min=np.nan, rsasa=np.nan,
                             frac_exposed=np.nan, call=f"no_model:{err}",
                             global_plddt=np.nan, folded_exposed=False))
            ctx.log(f"    {row.gene:<9} no structural model ({err})")
            continue
        seg = measured["spans"]["m"]
        fe = bool(seg["plddt_mean"] >= params["folded_plddt"]
                  and seg["rsasa_mean"] >= params["exposed_rsasa"])
        rows.append(dict(acc=row.acc, gene=row.gene, seed_id=row.seed_id,
                         motif=row.matched, pos=span[0], span=span[1] - span[0] + 1,
                         plddt=seg["plddt_mean"], plddt_min=seg["plddt_min"],
                         rsasa=seg["rsasa_mean"], frac_exposed=seg["frac_rsasa_ge_025"],
                         call=seg["call"],
                         global_plddt=measured["meta"].get("globalMetricValue"),
                         folded_exposed=fe))
        ctx.log(f"    {row.gene:<9} @{span[0]:<5d} pLDDT {seg['plddt_mean']:5.1f} "
                f"rSASA {seg['rsasa_mean']:.2f} {seg['call']}")
    real = pd.DataFrame(rows)
    real.to_csv(ctx.out("65_structure_real.tsv"), sep="\t", index=False)

    # ---- the same measurement on null matches: does a random match land on a
    # folded surface just as often? Same code path, same parameters.
    null_rows = []
    matches_path = ctx.out("20_motif_matches.tsv")
    if matches_path.exists():
        matches = pd.read_csv(matches_path, sep="\t", dtype={"acc": str})
        matches = matches[matches["family"].isin(["N_A", "N_B"])]
        best = (matches.sort_values(["n_ident", "seed_id"], ascending=[False, True])
                .drop_duplicates("seed_id"))
        rng = np.random.default_rng(int(ctx.cfg.get("seeds.random_state", 42)))
        picked = []
        for family in ("N_A", "N_B"):
            pool = best[best["family"] == family]
            if len(pool) == 0:
                continue
            take = min(params["null_per_family"], len(pool))
            picked.append(pool.iloc[rng.choice(len(pool), size=take, replace=False)])
        if picked:
            sample = pd.concat(picked)
            n_unmodelled = 0
            for row in sample.itertuples(index=False):
                span = (int(row.tpos),
                        int(row.tpos) + len(str(row.matched)) - 1)
                measured, err = _measure(row.acc, {"m": span}, cache, params)
                if measured is None:
                    excluded.append(dict(acc=row.acc, seed_id=row.seed_id,
                                         family=row.family, reason=err))
                    n_unmodelled += 1
                    continue
                seg = measured["spans"]["m"]
                if seg["call"] == "not_modelled":
                    excluded.append(dict(acc=row.acc, seed_id=row.seed_id,
                                         family=row.family,
                                         reason="span outside model"))
                    n_unmodelled += 1
                    continue
                null_rows.append(dict(acc=row.acc, seed_id=row.seed_id,
                                      family=row.family, pos=span[0],
                                      plddt=seg["plddt_mean"], rsasa=seg["rsasa_mean"],
                                      call=seg["call"],
                                      folded_exposed=bool(
                                          seg["plddt_mean"] >= params["folded_plddt"]
                                          and seg["rsasa_mean"] >= params["exposed_rsasa"])))
            if n_unmodelled:
                ctx.log(f"    null matches without a usable model: {n_unmodelled}")
            ctx.log(f"    null matches measured: {len(null_rows)}")
    null = pd.DataFrame(null_rows)
    if len(null):
        null.to_csv(ctx.out("65_structure_null.tsv"), sep="\t", index=False)

    # ---- calibration of the structural claim
    n_real = int(real["folded_exposed"].sum()) if len(real) else 0
    calib = {}
    if len(null):
        ge = int(null["folded_exposed"].sum())
        calib = dict(
            n_null=len(null),
            n_null_folded_exposed=ge,
            p_unconditional=round((1 + ge) / (1 + len(null)), 4),
            null_rate=round(ge / len(null), 4),
            real_rate=round(n_real / len(real), 4) if len(real) else None,
            n_real_folded_exposed=n_real,
        )
        ctx.log(f"    folded+exposed: real {n_real}/{len(real)} vs null {ge}/{len(null)} "
                f"-> p={calib['p_unconditional']}")

    facts = dict(
        parameters=params,
        query=query_facts,
        n_candidates=int(len(real)),
        n_folded_exposed=n_real,
        calls=real["call"].value_counts().to_dict() if len(real) else {},
        null_calls=null["call"].value_counts().to_dict() if len(null) else {},
        calibration=calib,
        excluded=excluded,
        caveat="pLDDT is model confidence and rSASA is computed on a predicted "
               "monomer; both are proxies for accessibility, not measurements. "
               "This stage is deterministic once the model cache is complete; a "
               "match whose model could not be retrieved is listed under "
               "'excluded' rather than counted as unexposed.",
    )
    write_json(ctx.out("65_structure_summary.json"), facts)
    return dict(outputs=["reports/65_structure_real.tsv",
                         "reports/65_structure_summary.json"]
                + (["reports/65_structure_null.tsv"] if len(null) else []),
                facts=facts)
