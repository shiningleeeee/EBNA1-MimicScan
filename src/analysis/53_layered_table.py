#!/usr/bin/env python3
"""Layered candidate table: every axis reported separately, nothing summed.

Axes
    seq      margin_z from the per-pair composition null
    mature   the matched target segment survives proteolytic maturation
    comp     the target compartment admits an antibody (secreted / plasma membrane)
    struct   both segments folded and exposed on the AlphaFold model

The union is NOT a score. In particular the `struct` axis is known to be wrong
for this problem: all four literature EBNA1 epitopes (386-405, 394-399, 402-406,
431-440) sit in disordered regions, so a filter requiring a folded EBNA1 segment
would have discarded every validated positive. `struct` is therefore reported
for information only, and the report says so.

Output
    reports/53_layered_table.tsv
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

ANTIBODY_COMPARTMENTS = (
    "secreted",
    "cell membrane",
    "cell surface",
    "extracellular",
    "gpi-anchor",
    "basement membrane",
)

# EBNA1 segments with published cross-reactive evidence, per config.yaml
# `published:` (query side only; the human-side coordinates are irrelevant here)
PUBLISHED_QUERY = {
    (386, 405): "HEPACAM/GlialCAM (structural, 7K7R)",
    (394, 399): "MEF2B (PPRRPP motif)",
    (402, 406): "CRYAB (RRPFF core)",
    (431, 440): "ANO2 (PGAIEQGPAD)",
}


def load_locations(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["Entry"]] = (
                (row.get("Subcellular location [CC]") or "")
                + " "
                + (row.get("Keywords") or "")
            ).lower()
    return out


def main():
    locations = load_locations(ROOT / "data/cache/uniprot_processing.tsv")
    with open(ROOT / "reports/52_biology_first.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    for row in rows:
        text = locations.get(row["acc"], "")
        row["comp_accessible"] = (
            "yes" if any(tok in text for tok in ANTIBODY_COMPARTMENTS) else "no"
        )
        row["seq_axis"] = "high" if float(row["margin_z"]) >= 6.0 else "moderate"
        qs, qe = int(row["ebna1_start"]), int(row["ebna1_end"])
        hit = [v for (s, e), v in PUBLISHED_QUERY.items() if s <= qe and qs <= e]
        row["overlaps_published_epitope"] = hit[0] if hit else ""

    dst = ROOT / "reports/53_layered_table.tsv"
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    def count(pred):
        return sum(1 for r in rows if pred(r))

    print(f"rows (survive seq+mature+length+identity filters): {len(rows)}")
    print(f"  mature retained                : {count(lambda r: r['retained']=='yes')}")
    print(f"  compartment accessible         : {count(lambda r: r['comp_accessible']=='yes')}")
    print(f"  both folded+exposed            : {count(lambda r: r['both_folded_exposed']=='yes')}")
    print(
        "  accessible AND both-folded     : "
        f"{count(lambda r: r['comp_accessible']=='yes' and r['both_folded_exposed']=='yes')}"
    )

    print("\n=== compartment-accessible candidates (antibody could in principle reach) ===")
    acc = [r for r in rows if r["comp_accessible"] == "yes"]
    for r in sorted(acc, key=lambda r: -float(r["margin_z"])):
        print(
            f"  {float(r['margin_z']):5.2f} {r['gene']:10s} {r['acc']} "
            f"EBNA1 {r['ebna1_start']}-{r['ebna1_end']} vs {r['t_start']}-{r['t_end']} "
            f"{r['aln_len']:>3s}aa {r['identity']:>5s}% "
            f"q {r['q_call'][:16]:16s} t {r['t_call'][:16]:16s}"
        )

    print("\n=== same list, keyed by which EBNA1 immunodominant region is used ===")
    from collections import Counter

    buckets = Counter()
    for r in acc:
        qs = int(r["ebna1_start"])
        if qs < 90:
            buckets["N-terminal 1-89"] += 1
        elif qs <= 330:
            buckets["Gly/Ala block 90-330"] += 1
        elif qs < 445:
            buckets["SPR-rich 331-444"] += 1
        elif qs < 476:
            buckets["445-475"] += 1
        elif qs <= 585:
            buckets["folded core 476-585"] += 1
        else:
            buckets["C-terminal 586-641"] += 1
    for k, v in buckets.most_common():
        print(f"  {k:26s} {v}")


if __name__ == "__main__":
    sys.exit(main())
