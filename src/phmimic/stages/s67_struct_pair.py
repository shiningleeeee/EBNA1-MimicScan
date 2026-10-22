"""s67 -- is a folded, reachable surface available on BOTH sides?

An antibody cross-reaction needs a paratope that can sit on both partners. If
one side's matched stretch is a disordered loop, the comparison is between a
conformation the antibody never sees and one it might. This stage measures
AlphaFold models at the matched spans and marks the pairs where both sides are
a folded, solvent-exposed surface.

Read it as the weakest of the annotation axes, not as evidence. pLDDT is
model confidence, not truth; rSASA is computed on the single predicted
conformer; neither knows about the actual binding-competent conformation, PTMs,
binding partners, or the conformational ensemble. And the axis has no positive
control: all four published EBNA1 epitopes are in disordered stretches, so
"folded and exposed on the EBNA1 side" keeps 0 of them. It can sort a shortlist
for a structural follow-up. It cannot filter, and it cannot rank.

Modelling costs one AlphaFold fetch plus ~5 s of SASA per protein, so only the
filtered set is measured; the query-side spans are deduplicated first.
"""
from __future__ import annotations

import time

import pandas as pd

from .. import mining, structure
from ..core import rel, write_json


def _call(summary, folded, exposed):
    """`structure.classify` with the thresholds taken from the config.

    The shared helper hardcodes 70 / 0.25 because that is what its own tests pin;
    here the numbers come from `mining.structure_pair.*` so a run is described by
    its config rather than by a constant two modules away. Same branches, same
    wording -- the only difference is where the cut-offs come from.
    """
    if summary.get("n", 0) == 0:
        return "not_modelled"
    is_exposed = summary["rsasa_mean"] >= exposed
    is_folded = summary["plddt_mean"] >= folded
    if not is_folded:
        return "disordered_or_loop_exposed" if is_exposed else "disordered_low_confidence"
    if is_exposed:
        return "surface_of_folded_domain"
    return "partly_buried" if summary["rsasa_mean"] >= 0.10 else "buried"


def _measure(acc, spans, cache, params, log):
    """One model fetch, many spans: returns {span: summary-or-None}."""
    spans = sorted(set(spans))
    try:
        pdb = structure.af_pdb(acc, cache)
    except Exception as exc:                      # 404, truncated, network
        log(f"    [{acc}] no model ({type(exc).__name__}: {exc})")
        return {s: None for s in spans}
    atom = structure.parse_pdb(pdb)
    sasa = structure.sasa_per_atom(atom, probe=params["probe_radius"],
                                   n_points=params["sphere_points"])
    stats = structure.residue_stats(atom, sasa)
    out = {}
    for span in spans:
        seg = structure.segment(atom, stats, *span)
        if seg.get("n", 0) == 0:
            out[span] = None                      # outside the modelled range
        else:
            out[span] = dict(
                seg,
                call=_call(seg, params["folded_plddt"], params["exposed_rsasa"]),
            )
    return out


def run(ctx):
    cfg = ctx.cfg
    log_path = cfg.logs / "67_struct_pair.log"
    log = ctx.log

    def _log(msg):
        log(msg)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(msg.split("    ", 1)[-1] + "\n")

    src = pd.read_csv(ctx.out("51_processing_state.tsv"), sep="\t", dtype=str)
    keep = (
        (src["retained"] == "yes")
        & (src["ax_length_ok"] == "True")
        & (src["ax_identity_ok"] == "True")
        & (src["ax_outside_gablock"] == "True")
        & (src["margin_z"].astype(float) >= float(cfg.get("mining.final.min_margin_z", 3.0)))
    )
    sel = src[keep].copy().reset_index(drop=True)
    _log(f"    selecting on maturation + length + identity + outside_gablock + "
         f"margin_z: {len(src)} -> {len(sel)}")

    cache = cfg.path_of("mining.structure_pair.cache_dir")
    query_acc = str(cfg.get("structure.query_accession", "P03211"))
    params = dict(
        probe_radius=float(cfg.get("mining.structure_pair.probe_radius", 1.4)),
        sphere_points=int(cfg.get("mining.structure_pair.sphere_points", 92)),
        folded_plddt=float(cfg.get("mining.structure_pair.folded_plddt", 70.0)),
        exposed_rsasa=float(cfg.get("mining.structure_pair.exposed_rsasa", 0.25)),
    )

    q_spans = [(int(a), int(b)) for a, b in zip(sel["ebna1_start"], sel["ebna1_end"])]
    t_spans = {}
    for acc, s, e in zip(sel["acc"], sel["t_start"], sel["t_end"]):
        t_spans.setdefault(str(acc), []).append((int(s), int(e)))

    t0 = time.time()
    q_meas = _measure(query_acc, q_spans, cache, params, _log)
    t_meas = {}
    for i, (acc, spans) in enumerate(sorted(t_spans.items()), 1):
        t_meas[acc] = _measure(acc, spans, cache, params, _log)
        if i % 10 == 0:
            _log(f"    modelled {i}/{len(t_spans)} targets ({time.time() - t0:.0f}s)")

    rows, calls = [], {"yes": 0, "no": 0, "not_modelled": 0}
    for row in sel.to_dict("records"):
        qspan = (int(row["ebna1_start"]), int(row["ebna1_end"]))
        tspan = (int(row["t_start"]), int(row["t_end"]))
        q = q_meas.get(qspan)
        t = t_meas.get(str(row["acc"]), {}).get(tspan)
        if q is None or t is None:
            flag = "not_modelled"
            note = "no usable model for query or target span"
        else:
            both = (q["call"] == "surface_of_folded_domain"
                    and t["call"] == "surface_of_folded_domain")
            flag = "yes" if both else "no"
            note = ""
        calls[flag] += 1
        rows.append(dict(
            row,
            q_plddt=round(q["plddt_mean"], 2) if q else "",
            q_rsasa=round(q["rsasa_mean"], 3) if q else "",
            q_call=q["call"] if q else "not_modelled",
            t_plddt=round(t["plddt_mean"], 2) if t else "",
            t_rsasa=round(t["rsasa_mean"], 3) if t else "",
            t_call=t["call"] if t else "not_modelled",
            both_folded_exposed=flag,
            struct_note=note,
        ))
    out = pd.DataFrame(rows)
    out.to_csv(ctx.out("52_biology_first.tsv"), sep="\t", index=False)

    payload = dict(
        n_in=int(len(src)),
        n_selected=int(len(sel)),
        selection=["maturation retained", "aln_len >= 20", "identity >= 40",
                   "outside Gly/Ala block", "margin_z >= min_margin_z"],
        min_margin_z=float(cfg.get("mining.final.min_margin_z", 3.0)),
        structure_params=params,
        n_both_folded_exposed=int(calls["yes"]),
        calls=calls,
        n_target_models_attempted=len(t_spans),
        caveat=("pLDDT is model confidence, not truth; rSASA is computed on the "
                "single predicted conformer; conformational ensembles, PTMs and "
                "binding partners are not represented. This axis has no positive "
                "control -- all four published EBNA1 epitopes are disordered -- so "
                "it annotates and sorts, it does not filter."),
    )
    write_json(ctx.out("67_struct_pair_summary.json"), payload)
    _log(f"    both folded+exposed: {calls['yes']}/{len(sel)} "
         f"(not modelled {calls['not_modelled']})")
    return dict(
        outputs=["reports/52_biology_first.tsv", "reports/67_struct_pair_summary.json"],
        facts=dict(n_selected=int(len(sel)),
                   n_both_folded_exposed=int(calls["yes"]),
                   n_not_modelled=int(calls["not_modelled"])),
    )
