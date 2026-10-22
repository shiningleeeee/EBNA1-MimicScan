"""s62 -- per-alignment two-sided composition null.

For each reported alignment, EBNA1's side and the target's side are each
reshuffled keeping their exact residue composition, and the same SW kernel is
re-run B times. `margin_z` is how many null standard deviations the observed
score sits above the null mean.

What this does and does not establish:

  * it asks "would random orderings of these two compositions have matched this
    well?" -- so a high margin_z is evidence that residue ORDER matters here;
  * it does NOT ask "is this alignment unusual among everything the search
    could have reported". The pairs tested here were already selected as the
    best a whole-proteome search found, so the selection-aware null in s63 is
    the calibration that actually applies to them.

`n_ge` saturates for short alignments (a 10-aa window's null often never
reaches the observation), so `margin_z` is the quantity to carry forward, not
the count.
"""
from __future__ import annotations

import pandas as pd

from .. import mining
from ..core import write_json


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "62_pair_null.log"
    res = mining.pair_null(
        cfg,
        cfg.path_of("mining.query"),
        ctx.out("cands_for_null.tsv"),
        ctx.out("46_pair_null.tsv"),
        log_path=log,
    )
    z = res["margin_z"].astype(float)
    summary = dict(
        n_alignments=int(len(res)),
        b=int(cfg.get("mining.pair_null.b")),
        seed=int(cfg.get("mining.pair_null.seed")),
        flank=int(cfg.get("mining.pair_null.flank")),
        margin_z_max=round(float(z.max()), 3),
        margin_z_median=round(float(z.median()), 3),
        margin_z_min=round(float(z.min()), 3),
        n_ge_zero=int((res["n_ge"] == 0).sum()),
        n_ge_le_2=int((res["n_ge"] <= 2).sum()),
        length_bias=dict(
            note=("n_ge==0 rate falls monotonically with alignment length "
                  "(short windows saturate), which is the opposite direction "
                  "from the E-value length bias"),
            rate_by_length_bin={
                str(k): round(float(v), 4)
                for k, v in _n_ge_zero_rate_by_length(res).items()
            },
        ),
    )
    write_json(ctx.out("62_pair_null_summary.json"), summary)
    ctx.log(f"    {len(res)} alignments; margin_z max {summary['margin_z_max']} "
            f"median {summary['margin_z_median']}; n_ge==0: {summary['n_ge_zero']}")
    return dict(
        outputs=["reports/46_pair_null.tsv", "reports/62_pair_null_summary.json"],
        facts=dict(n_alignments=summary["n_alignments"],
                   margin_z_max=summary["margin_z_max"],
                   margin_z_median=summary["margin_z_median"],
                   n_ge_zero=summary["n_ge_zero"]),
    )


def _n_ge_zero_rate_by_length(frame: pd.DataFrame) -> dict:
    bins = [(10, 15), (15, 20), (20, 30), (30, 50), (50, 10 ** 6)]
    out = {}
    for lo, hi in bins:
        sub = frame[(frame["alnlen"] >= lo) & (frame["alnlen"] < hi)]
        if len(sub):
            label = f"{lo}-{hi}" if hi < 10 ** 6 else f"{lo}+"
            out[label] = float((sub["n_ge"] == 0).mean())
    return out
