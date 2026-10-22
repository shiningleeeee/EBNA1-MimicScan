"""s68 -- annotate the selected set on independent axes, without fusing them.

Everything here is a separate annotation, deliberately not combined into one
number: a compartment call from UniProt, a crude band on the sequence statistic,
and whether the matched target span coincides with a region that someone has
already published as an EBNA1 mimic. The last one is a post-hoc audit of our own
sensitivity -- it is read-only, it never feeds the ordering, and it is reported
with the caveat that a hit there means we re-found known biology, not that we
found something new.

Note the direction of the compartment axis: it excludes CRYAB, which is
cytosolic and is one of the three published mimics. A filter that drops a known
positive is a filter that cannot be trusted to filter, so this axis is recorded
as an annotation and never used to remove a row.
"""
from __future__ import annotations

from collections import Counter

import pandas as pd

from .. import mining
from ..core import write_json


def _published_ranges(cfg):
    """acc -> [(target span, ebna1 span, reachable, gene)] from the config."""
    out = {}
    for entry in cfg.get("published", []) or []:
        try:
            start, end = int(entry["region"][0]), int(entry["region"][1])
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        try:
            q_start, q_end = int(entry["ebna1_region"][0]), int(entry["ebna1_region"][1])
        except (KeyError, IndexError, TypeError, ValueError):
            q_start, q_end = 0, 0
        out.setdefault(str(entry.get("acc", "")).upper(), []).append(
            (start, end, q_start, q_end, bool(entry.get("reachable")),
             str(entry.get("gene", ""))))
    return out


def run(ctx):
    cfg = ctx.cfg
    published = _published_ranges(cfg)

    # Same annotation file stage 30 used, so `comp_accessible` here and the
    # surface axis there cannot disagree about the same protein.
    annot = mining.load_annotation(cfg)

    src = pd.read_csv(ctx.out("52_biology_first.tsv"), sep="\t", dtype=str)
    z_sorted = src.assign(_z=src["margin_z"].astype(float)).sort_values(
        "_z", ascending=False, kind="mergesort")
    z_rank = {str(a): i + 1 for i, a in enumerate(z_sorted["acc"])}

    rows, axes = [], Counter()
    for row in src.to_dict("records"):
        acc = str(row["acc"])
        t_start, t_end = int(row["t_start"]), int(row["t_end"])
        q_start, q_end = int(row["ebna1_start"]), int(row["ebna1_end"])
        on_target, on_query = [], []
        for s, e, qs, qe, reachable, gene in published.get(acc.upper(), []):
            if reachable and mining.overlaps(t_start, t_end, s, e):
                on_target.append(gene)
        for entry_acc, entries in published.items():
            for s, e, qs, qe, reachable, gene in entries:
                # No `reachable` test here on purpose: the query-side question is
                # only "did this span land on a validated EBNA1 stretch". Whether
                # our channel can reach that pair through the *target* protein is
                # a different question, and it is what the target column records.
                if qs and mining.overlaps(q_start, q_end, qs, qe):
                    on_query.append(gene)
        z = float(row["margin_z"])
        seq_axis = "z>=8" if z >= 8 else ("z>=6" if z >= 6 else "z>=3")
        rows.append(dict(
            row,
            comp_accessible="yes" if annot.get(acc, {}).get("surface") else "no",
            comp_text=row.get("", ""),
            seq_axis=seq_axis,
            z_rank=z_rank.get(acc, ""),
            overlaps_published_target_region=";".join(sorted(set(on_target))),
            overlaps_published_query_region=";".join(sorted(set(on_query))),
        ))
        axes[row["both_folded_exposed"]] += 1
    out = pd.DataFrame(rows)
    out.to_csv(ctx.out("53_layered_table.tsv"), sep="\t", index=False)

    payload = dict(
        n_rows=int(len(out)),
        both_folded_exposed={k: int(v) for k, v in axes.most_common()},
        comp_accessible=dict(Counter(out["comp_accessible"])),
        seq_axis=dict(Counter(out["seq_axis"])),
        n_overlapping_published_target=int(
            (out["overlaps_published_target_region"] != "").sum()),
        n_overlapping_published_query=int(
            (out["overlaps_published_query_region"] != "").sum()),
        query_region_hits=dict(Counter(
            g for v in out["overlaps_published_query_region"] for g in v.split(";") if g)),
        published_entries_in_config=sum(len(v) for v in published.values()),
        note=("axes are reported side by side and never multiplied into a score. "
              "comp_accessible excludes CRYAB, a published mimic, so it is an "
              "annotation and never a filter. Both overlap columns are post-hoc "
              "sensitivity audits of our own run: the target column says we "
              "re-found a known pair, the query column says we landed on a "
              "validated EBNA1 stretch. Neither influences the ordering."),
    )
    write_json(ctx.out("68_layered_summary.json"), payload)
    ctx.log(f"    {len(out)} rows; accessible {payload['comp_accessible']}; "
            f"published overlap target={payload['n_overlapping_published_target']} "
            f"query={payload['n_overlapping_published_query']}")
    return dict(
        outputs=["reports/53_layered_table.tsv", "reports/68_layered_summary.json"],
        facts=dict(n_rows=int(len(out)),
                   n_accessible=int(payload["comp_accessible"].get("yes", 0)),
                   n_both_folded_exposed=int(axes.get("yes", 0)),
                   n_overlapping_published_query=payload["n_overlapping_published_query"]),
    )
