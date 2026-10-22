"""Quick probe: does the structural axis carry any information at all?

Run before integrating anything into the pipeline: fetch AlphaFold models for
EBNA1 and the top candidates, and look at pLDDT and rSASA over the motif.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phmimic import structure as st  # noqa: E402
from phmimic.core import load_config, read_fasta  # noqa: E402

cfg = load_config()
cache = Path("data/cache/alphafold")
cache.mkdir(parents=True, exist_ok=True)

ranked = pd.read_csv("reports/40_candidates_ranked.tsv", sep="\t", dtype={"acc": str})

print("=== EBNA1 (P03211) ===")
ebna = st.af_metadata("P03211", cache)
print(f"  global pLDDT {ebna['globalMetricValue']}  very_low {ebna['fractionPlddtVeryLow']}")
atom = st.parse_pdb(st.af_pdb("P03211", cache))
sasa = st.sasa_per_atom(atom)
stats = st.residue_stats(atom, sasa)
for lo, hi, label in [(90, 330, "Gly/Ala block"), (386, 405, "HEPACAM/PPRRPP region"),
                      (394, 399, "C1 PPRRPP"), (402, 406, "C2 RRPFF"),
                      (431, 440, "C3 PGAIEQGPAD")]:
    s = st.segment(atom, stats, lo, hi)
    print(f"  {label:24s} {lo}-{hi}  pLDDT {s['plddt_mean']:.1f}  "
          f"rSASA {s['rsasa_mean']:.2f}  {st.classify(s)}")

print()
print("=== candidates ===")
rows = []
for row in ranked.itertuples(index=False):
    acc, gene, pos, motif = row.acc, row.gene, int(row.pos), row.matched
    t0 = time.time()
    try:
        meta = st.af_metadata(acc, cache)
        atom = st.parse_pdb(st.af_pdb(acc, cache))
        sasa = st.sasa_per_atom(atom)
        stats = st.residue_stats(atom, sasa)
        seg = st.segment(atom, stats, pos, pos + len(motif) - 1)
        call = st.classify(seg)
    except Exception as exc:  # noqa: BLE001
        print(f"  {acc} {gene}: FAILED {type(exc).__name__}: {exc}")
        continue
    rows.append(dict(acc=acc, gene=gene, motif=motif, pos=pos,
                     plddt=seg["plddt_mean"], rsasa=seg["rsasa_mean"],
                     frac_exposed=seg["frac_rsasa_ge_025"], call=call,
                     global_plddt=meta.get("globalMetricValue")))
    print(f"  {gene:9s} @{pos:<5d} pLDDT {seg['plddt_mean']:5.1f}  "
          f"rSASA {seg['rsasa_mean']:.2f}  exposed {seg['frac_rsasa_ge_025']:.2f}  "
          f"{call:28s} ({time.time() - t0:.1f}s)")

pd.DataFrame(rows).to_csv("reports/probe_structure.tsv", sep="\t", index=False)
print(f"\nwrote reports/probe_structure.tsv ({len(rows)} rows)")
