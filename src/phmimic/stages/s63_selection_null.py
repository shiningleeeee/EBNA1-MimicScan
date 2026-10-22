"""s63 -- selection-aware null: run the WHOLE procedure on shuffled EBNA1.

The per-pair null in s62 is biased when applied to what a whole-proteome search
already chose to report: those pairs were selected for being a maximum, so of
course they beat their own shuffles. The bias is removed by doing to a shuffled
EBNA1 exactly what was done to the real one -- search the proteome, take every
reported alignment, score each against its own composition null -- and asking
where the real query lands in that distribution.

Two levels of randomization, which must never be conflated:

  outer  n_shuffles  composition-preserving shuffles of the whole 641-aa EBNA1,
                     each searched independently. Answers "is the NUMBER of
                     alignments EBNA1 reports unusual?"
  inner  b           shuffles per alignment inside pair_null. Answers "is THIS
                     one alignment better than its own two-sided order null?"
                     Resolution floor 1/(b+1).

The resolution floor is a real limit, not a footnote: with n_shuffles=200 the
smallest attainable selection p is 1/201 = 0.005.
"""
from __future__ import annotations

import random

import pandas as pd

from .. import mining
from ..core import ROOT, write_json

WORK = "work/selection_null"


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "63_selection_null.log"
    n_shuffles = int(cfg.get("mining.selection_null.n_shuffles", 200))
    seed = int(cfg.get("mining.selection_null.seed", 42))
    b = int(cfg.get("mining.selection_null.b", 2000))

    acc, seq = mining.first_record(cfg.path_of("mining.query"))
    rng = random.Random(seed)
    records = [("real", seq)]
    for i in range(1, n_shuffles + 1):
        chars = list(seq)
        rng.shuffle(chars)
        records.append((f"shuf{i:03d}", "".join(chars)))

    wdir = ROOT / WORK
    wdir.mkdir(parents=True, exist_ok=True)
    queries = wdir / "queries.fasta"
    mining.write_fasta(queries, records)

    m8 = mining.mmseqs_search(cfg, queries, wdir / "all.m8", WORK, log,
                              name="all", require_query_column=True)
    hits = mining.read_m8(m8, require_query_column=True)
    mining.write_m8(hits, m8)
    ctx.log(f"    {len(hits)} alignments over {hits['query'].nunique()} queries "
            f"({n_shuffles} shuffles, inner b={b})")

    by_query = {q: g for q, g in hits.groupby("query")}
    rows = []
    for name, seq_i in records:
        sub = by_query.get(name)
        if sub is None or sub.empty:
            rows.append(dict(query=name, n_align=0, n_p_lt_1e4=0, n_p_lt_1e3=0,
                             n_p_lt_1e2=0, best_p=float("nan"), med_p=float("nan")))
            continue
        cands = wdir / f"cands_{name}.tsv"
        sub[["target", "qstart", "qend", "tstart", "tend", "alnlen", "pident"]].to_csv(
            cands, sep="\t", index=False)
        qfa = wdir / f"q_{name}.fasta"
        mining.write_fasta(qfa, [(name, seq_i)])
        pn = mining.pair_null(cfg, qfa, cands, wdir / f"pn_{name}.tsv",
                              b=b, seed=seed, log_path=log)
        rows.append(dict(query=name, n_align=int(len(pn)),
                         n_p_lt_1e4=int((pn["p_pair"] < 1e-4).sum()),
                         n_p_lt_1e3=int((pn["p_pair"] < 1e-3).sum()),
                         n_p_lt_1e2=int((pn["p_pair"] < 1e-2).sum()),
                         best_p=float(pn["p_pair"].min()),
                         med_p=float(pn["p_pair"].median())))

    summary_frame = pd.DataFrame(rows).sort_values("query", kind="mergesort")
    summary_frame.to_csv(ctx.out("47_pair_null_batch.csv"), index=False)

    real = summary_frame[summary_frame["query"] == "real"].iloc[0]
    nulls = summary_frame[summary_frame["query"] != "real"]
    stats = {}
    for col in ("n_align", "n_p_lt_1e4", "n_p_lt_1e3", "n_p_lt_1e2", "best_p"):
        vals = nulls[col].astype(float)
        observed = float(real[col])
        ge = int((vals <= observed).sum()) if col == "best_p" else int((vals >= observed).sum())
        stats[col] = dict(
            real=observed,
            null_mean=round(float(vals.mean()), 4),
            null_sd=round(float(vals.std()), 4),
            null_min=round(float(vals.min()), 4),
            null_max=round(float(vals.max()), 4),
            n_null_ge=ge,
            p=(1 + ge) / (1 + len(vals)),
        )
    payload = dict(
        n_shuffles=n_shuffles,
        inner_b=b,
        seed=seed,
        resolution_floor=round(1.0 / (1 + n_shuffles), 6),
        inner_resolution_floor=round(1.0 / (1 + b), 8),
        real_n_align=int(real["n_align"]),
        null_n_align_mean=stats["n_align"]["null_mean"],
        null_n_align_sd=stats["n_align"]["null_sd"],
        p_selection_n_align=stats["n_align"]["p"],
        stats=stats,
        note=("best_p saturates at the inner resolution floor for both the real "
              "query and the shuffles, so it carries no information in the tail."),
    )
    write_json(ctx.out("63_selection_summary.json"), payload)
    ctx.log(f"    real {int(real['n_align'])} alignments vs null "
            f"{stats['n_align']['null_mean']} +- {stats['n_align']['null_sd']} "
            f"-> p={stats['n_align']['p']:.4f}")
    return dict(
        outputs=["reports/47_pair_null_batch.csv", "reports/63_selection_summary.json"],
        facts=dict(real_n_align=int(real["n_align"]),
                   null_mean=stats["n_align"]["null_mean"],
                   null_sd=stats["n_align"]["null_sd"],
                   p_selection=stats["n_align"]["p"],
                   n_shuffles=n_shuffles, inner_b=b),
    )
