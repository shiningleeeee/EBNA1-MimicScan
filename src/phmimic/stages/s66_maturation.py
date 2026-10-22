"""s66 -- does the matched segment exist on the protein an antibody can meet?

A sequence match is only a candidate for antibody cross-reactivity if the
matched residues are present in the species of the protein that circulates or
sits on a membrane. UniProt annotates the stretches that are proteolytically
removed during maturation -- signal peptide, propeptide, transit peptide, the
initiator methionine. A match entirely inside one of those cannot be recognised
on the mature protein.

This is the only filter in the whole project that passed a positive-control
check: it keeps every published linear mimic (CRYAB 11-16, ANO2 140-149,
MEF2B 272-277, and GlialCAM 370-389 on the far side of its signal peptide) and
still removes 12 of 363 alignments. MMP19 43-66, which an earlier round had put
first, is one of the 12: it sits inside PROPEP 19-97 and is cut away.

States, and what each means:
    retained_in_chain       span inside an annotated mature chain
    retained_unannotated    no removal feature covers it, but no chain either
                            -- NOT confirmed, reported as unknown
    fully_in_<kind>         entirely inside a removed feature
    spans_removed_boundary  partly removed; the fragment is not intact
"""
from __future__ import annotations

import re
from collections import Counter

import pandas as pd

from .. import mining
from ..core import rel, write_json

FEATURE_RE = re.compile(r"\b(SIGNAL|PROPEP|CHAIN|TRANSIT|INIT_MET|PEPTIDE) (\d+)\.\.(\d+)")
REMOVING = ("SIGNAL", "PROPEP", "TRANSIT", "INIT_MET")

FEATURE_COLUMNS = ("Signal peptide", "Propeptide", "Chain", "Transit peptide",
                   "Initiator methionine", "Peptide")


def _parse(path):
    out = {}
    for row in pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict("records"):
        feats = []
        for col in FEATURE_COLUMNS:
            for kind, start, end in FEATURE_RE.findall(row.get(col, "")):
                feats.append((kind, int(start), int(end)))
        out[str(row["Entry"])] = dict(
            length=int(row["Length"]) if row.get("Length") else None,
            features=feats,
            tissue=(row.get("Tissue specificity") or "").strip(),
        )
    return out


def _classify(features, start, end):
    covered, kinds = set(), set()
    for kind, fstart, fend in features:
        if kind not in REMOVING:
            continue
        lo, hi = max(start, fstart), min(end, fend)
        if lo <= hi:
            covered.update(range(lo, hi + 1))
            kinds.add(kind)
    span = end - start + 1
    if not covered:
        in_chain = any(k == "CHAIN" and s <= start and end <= e
                       for k, s, e in features)
        return ("retained_in_chain" if in_chain else "retained_unannotated"), 0
    if len(covered) == span:
        return f"fully_in_{sorted(kinds)[0].lower()}", len(covered)
    return "spans_removed_boundary", len(covered)


def run(ctx):
    cfg = ctx.cfg
    log = cfg.logs / "66_maturation.log"
    path = mining.ensure_processing_tsv(cfg, log)
    from ..core import sha256_file

    features = _parse(path)
    src = pd.read_csv(ctx.out("48_candidates_scored.tsv"), sep="\t", dtype=str)
    rows, calls = [], Counter()
    for row in src.to_dict("records"):
        acc = str(row["acc"])
        meta = features.get(acc)
        start, end = int(row["t_start"]), int(row["t_end"])
        if meta is None:
            call, removed, length, tissue = "no_uniprot_entry", 0, "", ""
        else:
            call, removed = _classify(meta["features"], start, end)
            length, tissue = meta["length"], meta["tissue"]
        calls[call] += 1
        rows.append(dict(
            row,
            target_length=length,
            processing_call=call,
            removed_residues_in_span=removed,
            span_len=end - start + 1,
            retained="yes" if call.startswith("retained") else "no",
            tissue_specificity=tissue,
        ))
    out = pd.DataFrame(rows)
    out.to_csv(ctx.out("51_processing_state.tsv"), sep="\t", index=False)

    n_unannotated = calls.get("retained_unannotated", 0)
    payload = dict(
        source=rel(path),
        source_sha256=sha256_file(path),
        n_classified=int(len(out)),
        calls={k: int(v) for k, v in calls.most_common()},
        n_retained=int((out["retained"] == "yes").sum()),
        n_confirmed=int((out["processing_call"] == "retained_in_chain").sum()),
        n_unannotated=n_unannotated,
        n_removed=int((out["retained"] == "no").sum()),
        caveat=("retained_unannotated means no removal feature covers the span, "
                "but no mature chain is annotated either; it is reported as "
                "state-unknown and is NOT counted as a confirmed pass."),
    )
    write_json(ctx.out("66_maturation_summary.json"), payload)
    ctx.log(f"    {payload['n_classified']} candidates: "
            f"confirmed {payload['n_confirmed']}, "
            f"unannotated {n_unannotated}, removed {payload['n_removed']}")
    return dict(
        outputs=["reports/51_processing_state.tsv", "reports/66_maturation_summary.json"],
        facts=dict(n_classified=payload["n_classified"],
                   n_confirmed=payload["n_confirmed"],
                   n_unannotated=n_unannotated,
                   n_removed=payload["n_removed"]),
    )
