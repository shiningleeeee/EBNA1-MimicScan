#!/usr/bin/env python3
"""Classify every candidate target span against UniProt maturation features.

Rationale: a sequence match is only a candidate for antibody cross-reactivity if
the matched segment is actually present in the protein species an antibody can
meet. UniProt annotates signal peptides, propeptides and transit peptides that
are proteolytically removed during maturation. A match lying inside a removed
region cannot be recognised on the mature protein (MMP19 43-66 is the case that
forced this filter).

Inputs
    data/cache/uniprot_processing.tsv   UniProt TSV with maturation features
    reports/candidates_v2_scored.tsv    all sequence-level candidates

Output
    reports/51_processing_state.tsv
"""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FEATURE_RE = re.compile(r"\b(SIGNAL|PROPEP|CHAIN|TRANSIT|INIT_MET|PEPTIDE) (\d+)\.\.(\d+)")
REMOVING = ("SIGNAL", "PROPEP", "TRANSIT", "INIT_MET")


def parse_features(tsv_path):
    """-> {acc: {'length': int, 'features': [(kind, start, end)], 'tissue': str}}"""
    out = {}
    with open(tsv_path, newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            acc = row["Entry"]
            feats = []
            for col in (
                "Signal peptide",
                "Propeptide",
                "Chain",
                "Transit peptide",
                "Initiator methionine",
                "Peptide",
            ):
                text = row.get(col) or ""
                for kind, start, end in FEATURE_RE.findall(text):
                    feats.append((kind, int(start), int(end)))
            out[acc] = {
                "length": int(row["Length"]) if row.get("Length") else None,
                "features": feats,
                "tissue": (row.get("Tissue specificity") or "").strip(),
                "location": (row.get("Subcellular location [CC]") or "").strip(),
                "keywords": (row.get("Keywords") or "").strip(),
            }
    return out


def removed_residues(features, start, end):
    """Residues in [start, end] covered by any maturation-removed feature."""
    covered = set()
    kinds = set()
    for kind, fstart, fend in features:
        if kind not in REMOVING:
            continue
        lo = max(start, fstart)
        hi = min(end, fend)
        if lo <= hi:
            covered.update(range(lo, hi + 1))
            kinds.add(kind)
    return len(covered), kinds


def classify(features, start, end, length):
    n_removed, kinds = removed_residues(features, start, end)
    span = end - start + 1
    if n_removed == 0:
        # Is the span inside an annotated mature chain? If no chain annotation
        # exists we cannot claim maturation, but absence of a removed region is
        # still the permissive case.
        in_chain = any(
            k == "CHAIN" and s <= start and end <= e for k, s, e in features
        )
        return ("retained_in_chain" if in_chain else "retained_unannotated"), 0
    if n_removed == span:
        kind = sorted(kinds)[0]
        return f"fully_in_{kind.lower()}", n_removed
    return "spans_removed_boundary", n_removed


def main():
    features = parse_features(ROOT / "data/cache/uniprot_processing.tsv")
    src = ROOT / "reports/candidates_v2_scored.tsv"
    dst = ROOT / "reports/51_processing_state.tsv"

    rows_out = []
    calls = Counter()
    with open(src, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["acc"]
            meta = features.get(acc)
            ts, te = int(row["t_start"]), int(row["t_end"])
            if meta is None:
                call, n_removed = "no_uniprot_entry", 0
                tissue = ""
                length = ""
            else:
                call, n_removed = classify(meta["features"], ts, te, meta["length"])
                tissue = meta["tissue"]
                length = meta["length"]
            calls[call] += 1
            rows_out.append(
                {
                    **row,
                    "target_length": length,
                    "processing_call": call,
                    "removed_residues_in_span": n_removed,
                    "span_len": te - ts + 1,
                    "retained": "yes" if call.startswith("retained") else "no",
                    "tissue_specificity": tissue,
                }
            )

    fieldnames = list(rows_out[0].keys())
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"classified {len(rows_out)} candidate spans -> {dst.name}")
    for call, n in calls.most_common():
        print(f"  {call:32s} {n:4d}")

    retained = [r for r in rows_out if r["retained"] == "yes"]
    print(f"\nretained (no maturation-removed residue in span): {len(retained)}")

    good = [
        r
        for r in retained
        if r["ax_length_ok"] == "True"
        and r["ax_identity_ok"] == "True"
        and float(r["margin_z"]) >= 3.0
    ]
    print(f"retained + len>=20 + ident>=40% + margin_z>=3: {len(good)}")
    for r in sorted(good, key=lambda r: -float(r["margin_z"])):
        flag = "CNS" if r["ax_cns"] == "True" else "   "
        surf = "SURF" if r["ax_surface"] == "True" else "    "
        print(
            f"  {float(r['margin_z']):5.2f} {r['gene']:12s} {r['acc']} "
            f"{r['ebna1_start']}-{r['ebna1_end']} vs {r['t_start']}-{r['t_end']} "
            f"{r['aln_len']:>3s}aa {r['identity']:>5s}% {flag} {surf} {r['processing_call']}"
        )


if __name__ == "__main__":
    sys.exit(main())
