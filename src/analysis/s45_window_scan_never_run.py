"""s45 -- unbiased whole-sequence window scan.

Every other candidate in this pipeline starts from a peptide that somebody else
already validated. That is a real limitation: searching outward from a known
epitope can only ever rediscover the region the epitope came from, and it hides
whatever else the protein might be doing. This stage removes the prior entirely.

The query protein is tiled end to end into overlapping windows. Each window is
scanned against the human proteome, and then the *same window with its own
residues permuted* -- identical composition, identical length -- is scanned the
same number of times. A window counts as interesting only if it finds more human
proteins than its own permutations do, which is the only honest way to ask
"is there order-specific information here?" for a protein that is 38% Gly/Ala.

Nothing in this stage knows where the published epitopes are. Their position is
read afterwards, as a blind test of the method, and is reported separately.
"""
from __future__ import annotations

import random

import numpy as np
import pandas as pd

from ..core import read_fasta, write_json
from ..scan import Proteome


def _count(proteome, motif, max_mismatch):
    """Distinct human proteins containing the motif within the mismatch budget."""
    hits = proteome.find(motif, max_mismatch)
    return len({h[0] for h in hits}), hits


def run(ctx):
    length = int(ctx.cfg.get("window_scan.length", 10))
    step = int(ctx.cfg.get("window_scan.step", 5))
    max_mm = int(ctx.cfg.get("window_scan.max_mismatch", 2))
    n_shuffles = int(ctx.cfg.get("window_scan.n_shuffles", 30))
    seed = int(ctx.cfg.get("window_scan.seed", 4242))
    blind = tuple(ctx.cfg.get("window_scan.blind_region", [386, 440]))

    sequence = read_fasta(ctx.cfg.path_of("data.query"))[0][2]
    proteome = Proteome(read_fasta(ctx.cfg.path_of("data.proteome")))
    ctx.log(f"    {len(sequence)} aa tiled into {length}-aa windows, step {step}; "
            f"{n_shuffles} composition shuffles each")

    rng = random.Random(seed)
    rows = []
    for start in range(0, len(sequence) - length + 1, step):
        window = sequence[start:start + length]
        n_real, hits = _count(proteome, window, max_mm)
        null = []
        for _ in range(n_shuffles):
            perm = list(window)
            rng.shuffle(perm)
            null.append(_count(proteome, "".join(perm), max_mm)[0])
        null_arr = np.asarray(null)
        ge = int((null_arr >= n_real).sum())
        rows.append(dict(
            start=start + 1, end=start + length, window=window,
            n_real=n_real,
            null_mean=round(float(null_arr.mean()), 2),
            null_max=int(null_arr.max()),
            null_median=float(np.median(null_arr)),
            fold_vs_null=round(n_real / null_arr.mean(), 3) if null_arr.mean() else None,
            p_unconditional=round((1 + ge) / (1 + n_shuffles), 4),
            p_conditional=round((1 + ge) / (1 + n_shuffles), 4),
            n_null_ge=ge,
            beats_all_nulls=bool(n_real > null_arr.max()),
            _hits=hits,
        ))
        if (start // step) % 20 == 0:
            ctx.log(f"      {start + 1}/{len(sequence)}")

    frame = pd.DataFrame([{k: v for k, v in r.items() if k != "_hits"} for r in rows])
    frame = frame.sort_values(["p_unconditional", "n_real", "start"],
                             ascending=[True, False, True]).reset_index(drop=True)
    frame.insert(0, "rank", np.arange(1, len(frame) + 1))
    frame.to_csv(ctx.out("45_window_scan.tsv"), sep="\t", index=False)

    # ---- global test: how many windows beat every one of their own shuffles?
    n_windows = len(frame)
    fl = 1.0 / (1 + n_shuffles)
    n_top = int((frame["p_unconditional"] <= fl + 1e-12).sum())
    expected_top = n_windows * fl
    # Poisson tail for the count under the null, reported rather than asserted
    from math import exp, factorial
    tail = 1.0 - sum(exp(-expected_top) * expected_top ** k / factorial(k)
                     for k in range(n_top))

    # ---- what the top windows actually matched, so the output is a list
    best = []
    for row in rows:
        if row["n_real"] == 0:
            continue
        if not row["beats_all_nulls"]:
            continue
        targets = sorted({h[0] for h in row["_hits"]})
        best.append(dict(start=row["start"], end=row["end"], window=row["window"],
                         n_real=row["n_real"], null_mean=row["null_mean"],
                         null_max=row["null_max"], targets=targets[:50],
                         n_targets=len(targets)))

    # ---- blind test: were we ever told where to look?
    in_blind = frame[(frame["start"] <= blind[1]) & (frame["end"] >= blind[0])]
    blind_facts = dict(
        region=list(blind),
        n_windows_in_region=int(len(in_blind)),
        best_rank=int(in_blind["rank"].min()) if len(in_blind) else None,
        best_percentile=round(100.0 * in_blind["rank"].min() / n_windows, 1)
        if len(in_blind) else None,
        max_n_real=int(in_blind["n_real"].max()) if len(in_blind) else None,
        any_beats_null=bool(in_blind["beats_all_nulls"].any()) if len(in_blind) else None,
    )

    facts = dict(
        parameters=dict(length=length, step=step, max_mismatch=max_mm,
                        n_shuffles=n_shuffles, seed=seed),
        n_windows=int(n_windows),
        n_windows_beating_all_shuffles=n_top,
        expected_by_chance=round(expected_top, 2),
        poisson_tail_p=round(max(tail, 0.0), 4),
        p_floor=round(fl, 4),
        top_windows=[dict(start=b["start"], window=b["window"], n_real=b["n_real"],
                          null_mean=b["null_mean"], null_max=b["null_max"],
                          n_targets=b["n_targets"], targets=b["targets"][:10])
                     for b in best[:10]],
        blind_test=blind_facts,
        note="windows are scored against permutations of themselves only; no "
             "published epitope position enters the scan. The blind test is an "
             "evaluation of the method, not an input to it.",
    )
    write_json(ctx.out("45_window_summary.json"), facts)

    ctx.log(f"    {n_top}/{n_windows} windows beat all {n_shuffles} of their own "
            f"shuffles (expected by chance {expected_top:.1f}, Poisson tail "
            f"p={tail:.3f})")
    if blind_facts["best_rank"]:
        ctx.log(f"    blind test: published region {blind[0]}-{blind[1]} ranks "
                f"{blind_facts['best_rank']}/{n_windows}")
    for b in best[:5]:
        ctx.log(f"      {b['start']}-{b['end']} {b['window']}: {b['n_real']} proteins "
                f"(null {b['null_mean']} max {b['null_max']})")

    return dict(outputs=["reports/45_window_scan.tsv", "reports/45_window_summary.json"],
                facts=facts)
