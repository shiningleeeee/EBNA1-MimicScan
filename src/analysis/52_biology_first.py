#!/usr/bin/env python3
"""Biology-first ranking: sequence signal is a filter, not the ranking axis.

Two questions are asked of each candidate, and reported separately:

  sequence axis  -- how surprising is the match under a composition-preserving
                    null (margin_z from the per-pair null model)
  biology axis   -- could one antibody plausibly see both segments?
                    * retained after maturation (no signal/propeptide removal)
                    * EBNA1 segment folded and exposed on the model
                    * target segment folded and exposed on the model

A match that fails the biology axis can still be a real sequence match; it just
cannot support an antibody cross-reactivity hypothesis. The two orderings are
deliberately not combined into one score.

Inputs
    reports/51_processing_state.tsv
    data/cache/alphafold/*.pdb (downloaded on demand)

Output
    reports/52_biology_first.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phmimic import structure  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
QUERY_ACC = "P03211"
CACHE = ROOT / "data/cache/alphafold"
FOLDED_PLDDT = 70.0
EXPOSED_RSASA = 0.25


def measure(acc, spans):
    """Classify several spans of one accession from a single parsed model."""
    try:
        path = structure.af_pdb(acc, CACHE)
    except (RuntimeError, OSError) as exc:
        # AlphaFoldDB has no model for every accession (404), and the API is
        # occasionally flaky. Absence of a model is a missing observation, not
        # a negative structural result, so it is recorded as such.
        print(f"    no model for {acc}: {exc}", flush=True)
        return {start: None for start, end in spans}
    if path is None:
        return {start: None for start, end in spans}
    atom = structure.parse_pdb(path)
    sasa = structure.sasa_per_atom(atom)
    stats = structure.residue_stats(atom, sasa)
    out = {}
    for start, end in spans:
        summary = structure.segment(atom, stats, start, end)
        summary["call"] = structure.classify(summary)
        out[(start, end)] = summary
    return out


def main():
    src = ROOT / "reports/51_processing_state.tsv"
    dst = ROOT / "reports/52_biology_first.tsv"

    rows = []
    with open(src, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["retained"] != "yes":
                continue
            if row["ax_length_ok"] != "True" or row["ax_identity_ok"] != "True":
                continue
            if row["ax_outside_gablock"] != "True":
                continue
            if float(row["margin_z"]) < 3.0:
                continue
            rows.append(row)

    print(f"sequence-filtered candidates: {len(rows)}")

    q_spans = sorted({(int(r["ebna1_start"]), int(r["ebna1_end"])) for r in rows})
    print(f"query spans to model: {len(q_spans)}")
    q_stats = measure(QUERY_ACC, q_spans)

    by_acc = {}
    for r in rows:
        by_acc.setdefault(r["acc"], set()).add((int(r["t_start"]), int(r["t_end"])))

    t_stats = {}
    for i, (acc, spans) in enumerate(sorted(by_acc.items()), 1):
        print(f"  [{i}/{len(by_acc)}] {acc} ...", flush=True)
        t_stats[acc] = measure(acc, sorted(spans))

    out_rows = []
    for r in rows:
        qs = q_stats.get((int(r["ebna1_start"]), int(r["ebna1_end"])))
        ts = t_stats.get(r["acc"], {}).get((int(r["t_start"]), int(r["t_end"])))
        q_call = qs["call"] if qs else "not_modelled"
        t_call = ts["call"] if ts else "not_modelled"
        q_ok = bool(qs) and qs.get("n", 0) > 0
        t_ok = bool(ts) and ts.get("n", 0) > 0
        both = q_call == "surface_of_folded_domain" and t_call == "surface_of_folded_domain"
        out_rows.append(
            {
                **r,
                "q_plddt": f"{qs['plddt_mean']:.1f}" if q_ok else "",
                "q_rsasa": f"{qs['rsasa_mean']:.2f}" if q_ok else "",
                "q_call": q_call,
                "t_plddt": f"{ts['plddt_mean']:.1f}" if t_ok else "",
                "t_rsasa": f"{ts['rsasa_mean']:.2f}" if t_ok else "",
                "t_call": t_call,
                "both_folded_exposed": "yes" if both else "no",
            }
        )

    fieldnames = list(out_rows[0].keys())
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    out_rows.sort(key=lambda r: -float(r["margin_z"]))
    print(f"\n=== sequence priority (margin_z) ===")
    for r in out_rows[:12]:
        print(
            f"  {float(r['margin_z']):5.2f} {r['gene']:10s} {r['acc']} "
            f"{r['aln_len']:>3s}aa {r['identity']:>5s}%  "
            f"q {r['q_call'][:14]:14s} t {r['t_call'][:14]:14s} "
            f"{'BOTH' if r['both_folded_exposed'] == 'yes' else ''}"
        )

    both = [r for r in out_rows if r["both_folded_exposed"] == "yes"]
    print(f"\n=== biology priority: both segments folded and exposed ({len(both)}) ===")
    for r in sorted(both, key=lambda r: -float(r["margin_z"])):
        print(
            f"  {float(r['margin_z']):5.2f} {r['gene']:10s} {r['acc']} "
            f"EBNA1 {r['ebna1_start']}-{r['ebna1_end']} vs {r['t_start']}-{r['t_end']} "
            f"{r['aln_len']:>3s}aa {r['identity']:>5s}% "
            f"qPLDDT {r['q_plddt']:>5s} qrSASA {r['q_rsasa']} "
            f"tPLDDT {r['t_plddt']:>5s} trSASA {r['t_rsasa']}"
        )


if __name__ == "__main__":
    sys.exit(main())
