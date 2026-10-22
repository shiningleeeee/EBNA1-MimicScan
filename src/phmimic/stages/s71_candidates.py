"""s71 -- assemble the candidate table, and audit our own sensitivity beside it.

Two products, kept strictly apart:

  the table   the filtered, annotated set, one row per alignment, every field
              traceable to the stage that produced it, plus the two columns that
              keep the reading honest: `unknowns` (what we did not check) and
              `manual_judgement` (which of the remaining fields a human should
              weigh, and why). No composite score: the axes are stacked, not
              multiplied.

  the audit   for each published EBNA1 mimic in the config, whether this
              unbiased channel found it at all. It is computed last, written to
              its own file, and nothing downstream reads it. Whatever it says
              about sensitivity, it must not change the ordering above it --
              that is exactly the loop that made an earlier round worthless.

The audit's honest answer for this channel is uncomfortable and belongs in the
report: of the three published linear mimics, this full-length local-alignment
channel recovers none, because each is a short low-complexity epitope that a
length-and-composition-biased aligner cannot score. The epitope-seeded channel
in the reproducible pipeline recovered all three. Two channels, opposite
failures, stated side by side.
"""
from __future__ import annotations

import pandas as pd

from .. import mining
from ..core import write_json

TEXT_COLUMNS = dict(
    gene="gene", acc="acc",
    ebna1_start="ebna1_start", ebna1_end="ebna1_end",
    target_start="t_start", target_end="t_end",
    aln_len="aln_len", identity_pct="identity",
    margin_z="margin_z", obs_score="obs_score",
    null_mean="null_mean", null_sd="null_sd", null_max="null_max",
    p_pair="p_pair", n_ge="n_ge", b="b",
    processing_call="processing_call", compartment_accessible="comp_accessible",
    comp_text="comp_text", cns_annotation="cns_annotation",
    surface_annotation="surface_annotation",
    query_struct="q_call", target_struct="t_call",
    both_folded_exposed="both_folded_exposed", seq_axis="seq_axis",
    z_rank="z_rank",
    overlaps_published_target_region="overlaps_published_target_region",
    overlaps_published_query_region="overlaps_published_query_region",
    ebna1_window="q_win", target_window="t_win",
)

UNKNOWNS = [
    "no antibody assay of these peptides in any reported work we could find",
    "no knowledge of conformational epitopes (a peptide screen cannot see them)",
    "no PTM state (phosphorylation, citrullination, glycosylation) modelled",
    "no expression evidence in MS lesions or CSF for most entries",
    "rSASA from a single AlphaFold conformer, not a measured surface",
]


def _join_m8(m8_path, rows):
    try:
        m8 = mining.read_m8(m8_path)
    except FileNotFoundError:
        return {}
    keys = {}
    for r in m8.itertuples(index=False):
        key = (str(r.target), int(r.tstart), int(r.tend))
        keys[key] = (r.evalue, r.bits, r.qaln, r.taln)
    out = {}
    for row in rows:
        key = (str(row["acc"]), int(row["t_start"]), int(row["t_end"]))
        out[key] = keys.get(key)
    return out


def run(ctx):
    cfg = ctx.cfg
    src = pd.read_csv(ctx.out("53_layered_table.tsv"), sep="\t", dtype=str)
    rows = src.to_dict("records")

    iedb = {}
    try:
        for row in pd.read_csv(ctx.out("54_iedb_evidence.tsv"), sep="\t",
                               dtype=str).fillna("").to_dict("records"):
            iedb[str(row["acc"]).upper()] = row
    except FileNotFoundError:
        ctx.log("    54_iedb_evidence.tsv missing; IEDB columns marked not_queried")

    m8 = _join_m8(ctx.out("mmseqs_ebna1.m8"), rows)

    final = []
    for row in rows:
        acc = str(row["acc"]).upper()
        rec = iedb.get(acc)
        extra = m8.get((str(row["acc"]), int(row["t_start"]), int(row["t_end"])))
        item = {}
        for out_col, src_col in TEXT_COLUMNS.items():
            item[out_col] = row.get(src_col, "")
        item["maturation"] = row.get("processing_call", "")
        item["evalue"] = extra[0] if extra else ""
        item["bits"] = extra[1] if extra else ""
        item["qaln"] = extra[2] if extra else ""
        item["taln"] = extra[3] if extra else ""
        item["iedb_queried"] = rec["iedb_queried"] if rec else "not_queried"
        item["iedb_bcell_records"] = rec["n_records"] if rec else "not_queried"
        item["iedb_on_segment"] = rec["n_overlap"] if rec else "not_queried"
        item["unknowns"] = "; ".join(UNKNOWNS)
        item["manual_judgement"] = (
            "sequence statistic only; nothing here measures binding"
            if str(row.get("both_folded_exposed")) != "yes"
            else "sequence statistic plus a folded, exposed surface on both sides; "
                 "still does not measure binding")
        final.append(item)
    # Numeric, not lexicographic: everything above is read as str so the table
    # round-trips exactly, but "11.16" must not sort below "9.23".
    out = pd.DataFrame(final)
    out["_z"] = pd.to_numeric(out["margin_z"], errors="coerce")
    out["_len"] = pd.to_numeric(out["aln_len"], errors="coerce")
    out = out.sort_values(["_z", "_len", "acc"], ascending=[False, False, True],
                          kind="mergesort").drop(columns=["_z", "_len"])
    out = out.reset_index(drop=True)
    out.to_csv(ctx.out("55_candidates_final.tsv"), sep="\t", index=False)

    # ---- sensitivity audit of this channel, computed last, read by nothing ----
    audit = []
    for entry in cfg.get("published", []) or []:
        acc = str(entry.get("acc", "")).upper()
        region = entry.get("region") or [None, None]
        reachable = bool(entry.get("reachable"))
        hits = [r for r in rows if str(r["acc"]).upper() == acc]
        overlap = [r for r in hits
                   if region[0] and mining.overlaps(int(r["t_start"]), int(r["t_end"]),
                                                    int(region[0]), int(region[1]))]
        best = min((float(r["evalue"]) for r in (overlap or hits)
                    if str(r.get("evalue", "")).strip() not in ("", "nan")),
                   default=None)
        audit.append(dict(
            gene=entry.get("gene", ""), acc=acc,
            region=f"{region[0]}-{region[1]}" if region[0] else "",
            reachable=reachable,
            in_unbiased_channel="yes" if hits else "no",
            n_alignments_in_channel=len(hits),
            overlaps_claimed_region=len(overlap),
            best_evalue="" if best is None else round(best, 3),
            verdict="recovered" if overlap else "not_recovered",
        ))
    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(ctx.out("71_benchmark_audit.tsv"), sep="\t", index=False)
    reachable = [a for a in audit if a["reachable"]]
    n_rec = sum(1 for a in reachable if a["verdict"] == "recovered")

    payload = dict(
        n_rows=int(len(out)),
        n_columns=int(out.shape[1]),
        columns=list(out.columns),
        n_iedb_queried=int((out["iedb_queried"] == "yes").sum()),
        n_iedb_not_queried=int((out["iedb_queried"] == "not_queried").sum()),
        benchmark=dict(
            n_reachable=len(reachable),
            n_recovered=n_rec,
            note=("post-hoc sensitivity audit of the unbiased channel only; it is "
                  "written after the table and is not an input to any ranking. "
                  "The epitope-seeded channel recovers these differently -- the "
                  "two failure modes are reported side by side."),
        ),
        note=("no composite score: axes are stacked, never multiplied. "
              "`unknowns` and `manual_judgement` are part of the deliverable."),
    )
    write_json(ctx.out("71_candidates_summary.json"), payload)
    ctx.log(f"    final table {out.shape[0]} x {out.shape[1]}; benchmark audit "
            f"{n_rec}/{len(reachable)} reachable published mimics recovered")
    return dict(
        outputs=["reports/55_candidates_final.tsv", "reports/71_benchmark_audit.tsv",
                 "reports/71_candidates_summary.json"],
        facts=dict(n_rows=int(len(out)), n_columns=int(out.shape[1]),
                   bench_recovered=n_rec, bench_reachable=len(reachable)),
    )
