"""s64 -- score every alignment on eight declared axes, then calibrate the LIST.

The ordering is hypothesis generation, not detection. It exists to put the
entries most worth an experiment at the top and to say exactly why each one is
there; it cannot establish cross-reactivity, and the calibration below is what
keeps anyone from reading it as if it could.

The eight axes are declared before the result is seen and none is fitted to the
three published mimics (which are audited in s71, never used as input):

    order_signal    per-pair composition null margin_z >= z_min
    length_ok       aligned length >= 20
    identity_ok     percent identity >= 40
    target_normal   target window entropy >= ent_min and max residue <= maxfrac_max
    query_normal    same, on the EBNA1 window
    cns             brain / CNS expression evidence (annotation)
    surface         antibody-reachable compartment (annotation)
    outside_gablock EBNA1 window does not overlap the Gly/Ala block

Two calibrations are reported, because they answer different questions:

  quality  does the real query's BEST entry beat the shuffles' best entry?
  quantity does the real query put a larger FRACTION of its alignments at or
           above a given composite than the shuffles do?

Quantity must be compared as a rate: the real EBNA1 reports roughly 18x more
alignments than a shuffle, so comparing raw counts inflates the effect ~19-fold.
"""
from __future__ import annotations

from collections import Counter

import pandas as pd

from .. import mining
from ..core import ROOT, write_json

WORK = ROOT / "work" / "selection_null"


def _build(query_name, pn, qseq, proteome, annotation, cfg):
    ga_block = tuple(cfg.get("composition.gar_block", [90, 330]))
    z_min = float(cfg.get("mining.composite.z_min", 3.0))
    ent_min = float(cfg.get("mining.composite.ent_min", 2.8))
    mf_max = float(cfg.get("mining.composite.maxfrac_max", 0.30))
    rows = []
    for r in pn.itertuples(index=False):
        target = str(r.target)
        tseq = proteome.get(target)
        if not tseq:
            continue
        qwin = qseq[int(r.qstart) - 1: int(r.qend)]
        twin = tseq[int(r.tstart) - 1: int(r.tend)]
        q_ent, q_mf = mining.entropy(qwin)
        t_ent, t_mf = mining.entropy(twin)
        ann = annotation.get(target, {})
        ident = float(r.pident) if r.pident != "-" else 0.0
        axes = dict(
            order_signal=float(r.margin_z) >= z_min,
            length_ok=int(r.alnlen) >= 20,
            identity_ok=ident >= 40.0,
            target_normal=t_ent >= ent_min and t_mf <= mf_max,
            query_normal=q_ent >= ent_min and q_mf <= mf_max,
            cns=bool(ann.get("cns", False)),
            surface=bool(ann.get("surface", False)),
            outside_gablock=not mining.overlaps(int(r.qstart), int(r.qend), *ga_block),
        )
        rows.append(dict(
            query=query_name, acc=target, gene=str(r.gene),
            ebna1_start=int(r.qstart), ebna1_end=int(r.qend),
            t_start=int(r.tstart), t_end=int(r.tend),
            aln_len=int(r.alnlen), identity=ident,
            obs_score=int(r.obs_score),
            null_mean=round(float(r.null_mean), 4),
            null_sd=round(float(r.null_sd), 4),
            null_max=int(r.null_max),
            margin_z=round(float(r.margin_z), 4),
            p_pair=float(r.p_pair),
            n_ge=int(r.n_ge), b=int(r.b),
            q_win=qwin, t_win=twin,
            q_ent=round(q_ent, 3), q_maxfrac=round(q_mf, 3),
            t_ent=round(t_ent, 3), t_maxfrac=round(t_mf, 3),
            **{f"ax_{k}": v for k, v in axes.items()},
            score=int(sum(axes.values())),
        ))
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(["score", "margin_z", "aln_len", "acc"],
                             ascending=[False, False, False, True],
                             kind="mergesort").reset_index(drop=True)


def _calibrate(real_df, null_df, min_score):
    out = {}
    keys = ["max_score", "n_ge_4", "n_ge_5", "n_ge_6"]
    real_stats = {"max_score": int(real_df["score"].max()),
                  **{f"n_ge_{k}": int((real_df['score'] >= k).sum()) for k in (4, 5, 6)}}
    real_stats["frac_ge_min_score"] = float((real_df["score"] >= min_score).mean())
    for key, observed in real_stats.items():
        if null_df.empty or key not in null_df.columns:
            out[key] = dict(real=observed, null_mean=None, p=None)
            continue
        vals = null_df[key].astype(float)
        ge = int((vals >= observed).sum())
        out[key] = dict(
            real=observed,
            null_mean=round(float(vals.mean()), 4),
            null_sd=round(float(vals.std()), 4),
            null_min=round(float(vals.min()), 4),
            null_max=round(float(vals.max()), 4),
            n_null_ge=ge,
            p=round((1 + ge) / (1 + len(vals)), 4),
        )
    return out


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "64_rank_composite.log"
    composite_b = int(cfg.get("mining.composite.b", 10000))
    min_score = int(cfg.get("mining.composite.min_score", 6))
    interval_b = int(cfg.get("mining.pair_null.b_interval", 1000000))

    proteome = mining.load_fasta(cfg.path_of("data.proteome"))
    query_fa = cfg.path_of("mining.query")
    qseq = list(mining.load_fasta(query_fa).values())[0]
    annotation = mining.load_annotation(cfg)
    ctx.log(f"    proteome {len(proteome)}; query {len(qseq)} aa; "
            f"annotation {len(annotation)}")

    real_pn = mining.pair_null(cfg, query_fa, ctx.out("cands_for_null.tsv"),
                               ROOT / "work" / "selection_null" / "pn_real_composite.tsv",
                               b=composite_b, log_path=log)
    real_df = _build("real", real_pn, qseq, proteome, annotation, cfg)
    real_df.to_csv(ctx.out("48_candidates_scored.tsv"), sep="\t", index=False)
    ctx.log(f"    real: {len(real_df)} rows; score distribution "
            f"{dict(sorted(Counter(real_df['score']).items()))}")

    null_rows = []
    for cands in sorted(WORK.glob("cands_shuf*.tsv")):
        name = cands.stem.replace("cands_", "")
        qfa = WORK / f"q_{name}.fasta"
        if not qfa.exists():
            continue
        pn = mining.pair_null(cfg, qfa, cands,
                              WORK / f"pn_{name}_composite.tsv",
                              b=composite_b, log_path=log)
        frame = _build(name, pn, list(mining.load_fasta(qfa).values())[0],
                       proteome, annotation, cfg)
        if frame.empty:
            continue
        null_rows.append(dict(
            query=name, n_rows=int(len(frame)),
            max_score=int(frame["score"].max()),
            **{f"n_ge_{k}": int((frame["score"] >= k).sum()) for k in (4, 5, 6)},
            frac_ge_min_score=round(float((frame["score"] >= min_score).mean()), 4),
        ))
    null_df = pd.DataFrame(null_rows)
    if not null_df.empty:
        null_df = null_df.sort_values("query", kind="mergesort")
        null_df.to_csv(ctx.out("48_null_composite.csv"), index=False)

    calibration = _calibrate(real_df, null_df, min_score)

    # ---- per-candidate interval for the top band (B far above the table's) ----
    top = real_df[real_df["score"] == int(real_df["score"].max())]
    interval_rows = []
    if not top.empty:
        top_cands = top[["acc", "ebna1_start", "ebna1_end", "t_start", "t_end",
                         "aln_len", "identity"]].rename(
            columns={"acc": "target", "ebna1_start": "qstart", "ebna1_end": "qend",
                     "t_start": "tstart", "t_end": "tend", "aln_len": "alnlen",
                     "identity": "pident"})
        top_path = ROOT / "work" / "selection_null" / "cands_top.tsv"
        top_cands.to_csv(top_path, sep="\t", index=False)
        pn_top = mining.pair_null(cfg, query_fa, top_path,
                                  ctx.out("46b_top_interval.tsv"), b=interval_b,
                                  log_path=log)
        for r in pn_top.itertuples(index=False):
            n = int(r.n_ge)
            interval_rows.append(dict(
                target=str(r.target), gene=str(r.gene),
                qstart=int(r.qstart), qend=int(r.qend),
                obs_score=int(r.obs_score), null_max=int(r.null_max),
                margin_z=round(float(r.margin_z), 4),
                n_ge=n, b=int(r.b),
                p_pair=float(r.p_pair),
                # Clopper-Pearson-free binomial interval: exact order statistics
                ci_low=round(n / (1 + int(r.b)), 8) if n else 0.0,
                ci_high=round((n + 1) / (1 + int(r.b)), 8),
            ))

    payload = dict(
        composite_b=composite_b,
        min_score=min_score,
        interval_b=interval_b,
        axes=[c for c in real_df.columns if c.startswith("ax_")],
        score_distribution={str(k): int(v)
                            for k, v in sorted(Counter(real_df["score"]).items())},
        n_at_or_above_min_score=int((real_df["score"] >= min_score).sum()),
        calibration=calibration,
        interval_top_n=int(len(top)),
        note=("max_score calibration is the honest one for 'is the best entry "
              "special'; the n_ge_* rows must be read as FRACTIONS, not counts, "
              "because the real query reports ~18x more alignments than a shuffle "
              "and therefore has more chances to reach any threshold."),
    )
    write_json(ctx.out("64_composite_summary.json"), payload)

    ctx.log(f"    composite calibration: max_score real="
            f"{calibration['max_score']['real']} vs null "
            f"{calibration['max_score']['null_mean']} +- {calibration['max_score']['null_sd']} "
            f"-> p={calibration['max_score']['p']}")
    return dict(
        outputs=["reports/48_candidates_scored.tsv", "reports/48_null_composite.csv",
                 "reports/46b_top_interval.tsv", "reports/64_composite_summary.json"],
        facts=dict(n_rows=int(len(real_df)),
                   max_score=int(real_df["score"].max()),
                   n_at_or_above_min_score=payload["n_at_or_above_min_score"],
                   p_max_score=calibration["max_score"]["p"],
                   p_frac=calibration["frac_ge_min_score"]["p"]),
    )
