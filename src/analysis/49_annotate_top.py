#!/usr/bin/env python3
"""49_annotate_top.py — add structural state and literature evidence to the
shortlisted alignments, for BOTH sides of every match.

For each shortlisted entry we ask two independent questions that the sequence
axes cannot answer:

  1. Structural state — is the matched segment, on the EBNA1 side and on the
     human side, inside a folded domain (pLDDT >= 70) or in a disordered
     region?  This matters because an antibody can only reach a folded,
     exposed surface, while a disordered segment is the natural substrate
     for linear-epitope cross-reactivity.  The two published linear mimics
     sit in disordered regions; if our unbiased top hits sit in folded cores
     instead, that is a mechanistic difference, not a better candidate.

  2. Literature — has the human protein ever been reported in an
     autoantibody / autoantigen / multiple-sclerosis / EBV context?  Counts
     are recorded as context, never as a score.

Neither number enters the ranking. They are reported so a reader can decide
which entry is worth an experiment.

Outputs: reports/49_top_annotated.tsv (and prints a markdown table).
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
REPORTS = ROOT / "reports"
CACHE = ROOT / "data" / "cache"

from phmimic import structure  # noqa: E402


def pubmed_count(gene: str) -> int:
    terms = f'({gene}[tiab]) AND (autoantibod*[tiab] OR autoantigen*[tiab] OR "multiple sclerosis"[tiab] OR "Epstein-Barr"[tiab])'
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
           "&retmode=json&retmax=0&term=" + urllib.parse.quote(terms))
    try:
        with urllib.request.urlopen(url, timeout=30) as fh:
            return int(json.load(fh)["esearchresult"]["count"])
    except Exception as exc:  # noqa: BLE001
        print(f"[49] pubmed {gene}: {exc}")
        return -1


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--min-score", type=int, default=7)
    args = ap.parse_args(argv)

    df = pd.read_csv(REPORTS / "48_candidates_scored.tsv", sep="\t")
    df = df.sort_values(["score", "margin_z", "aln_len"],
                        ascending=[False, False, False])
    sel = df[df["score"] >= args.min_score].head(args.top).copy()
    print(f"[49] annotating {len(sel)} entries (score >= {args.min_score})")

    af_cache = CACHE / "alphafold"
    af_cache.mkdir(parents=True, exist_ok=True)

    models: dict[str, dict | None] = {}

    def model_for(acc: str):
        if acc in models:
            return models[acc]
        try:
            pdb = structure.af_pdb(acc, af_cache)
            atom = structure.parse_pdb(pdb)
            sasa = structure.sasa_per_atom(atom)
            stats = structure.residue_stats(atom, sasa)
            models[acc] = {"atom": atom, "stats": stats}
        except Exception as exc:  # noqa: BLE001
            print(f"[49] {acc}: no model ({exc})")
            models[acc] = None
        return models[acc]

    rows = []
    for _, r in sel.iterrows():
        m = model_for(str(r["acc"]))
        if m is None:
            t_sum = {"n": 0}
        else:
            t_sum = structure.segment(m["atom"], m["stats"],
                                      int(r["t_start"]), int(r["t_end"]))
        qm = model_for("P03211")
        q_sum = structure.segment(qm["atom"], qm["stats"],
                                  int(r["ebna1_start"]), int(r["ebna1_end"]))
        rows.append({
            **r.to_dict(),
            "ebna1_plddt": round(q_sum.get("plddt_mean", float("nan")), 1),
            "ebna1_rsasa": round(q_sum.get("rsasa_mean", float("nan")), 3),
            "ebna1_state": structure.classify(q_sum),
            "target_plddt": round(t_sum.get("plddt_mean", float("nan")), 1),
            "target_rsasa": round(t_sum.get("rsasa_mean", float("nan")), 3),
            "target_state": structure.classify(t_sum),
            "pubmed_auto": pubmed_count(str(r["gene"])),
        })
        print(f"[49] {r['gene']:10s} q {r['ebna1_start']}-{r['ebna1_end']} "
              f"pLDDT {rows[-1]['ebna1_plddt']:5.1f} {rows[-1]['ebna1_state']:26s} | "
              f"t {rows[-1]['target_state']:26s} | pubmed {rows[-1]['pubmed_auto']}")

    out = pd.DataFrame(rows)
    out.to_csv(REPORTS / "49_top_annotated.tsv", sep="\t", index=False)

    cols = ["gene", "acc", "ebna1_start", "ebna1_end", "aln_len", "identity",
            "margin_z", "score", "ebna1_state", "target_state", "pubmed_auto"]
    print("\n" + out[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
