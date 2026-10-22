"""Probe 2: is "the motif sits on a folded, exposed surface" special?

For every null seed (a composition-permuted real seed, or an arbitrary peptide),
take its best match in the proteome and compute the same pLDDT / rSASA numbers
the real candidates got. If the null lands on folded surfaces just as often,
the structural axis is another description of composition.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from phmimic import structure as st  # noqa: E402

CACHE = Path("data/cache/alphafold")
CACHE.mkdir(parents=True, exist_ok=True)

matches = pd.read_csv("reports/20_motif_matches.tsv", sep="\t", dtype={"acc": str})
matches = matches[matches["family"].isin(["N_A", "N_B"])]

# best match per null seed: most identities, then longest
matches = matches.sort_values(["n_ident", "seed_id"], ascending=[False, True])
best = matches.drop_duplicates("seed_id")
print(f"null seeds with a hit: {len(best)}")

rng = np.random.default_rng(7)
sample = []
for family, n in (("N_A", 40), ("N_B", 40)):
    pool = best[best["family"] == family]
    take = min(n, len(pool))
    idx = rng.choice(len(pool), size=take, replace=False)
    sample.append(pool.iloc[idx])
sample = pd.concat(sample)

rows = []
for i, row in enumerate(sample.itertuples(index=False), start=1):
    acc, pos, ident = row.acc, int(row.tpos), int(row.n_ident)
    try:
        atom = st.parse_pdb(st.af_pdb(acc, CACHE))
        sasa = st.sasa_per_atom(atom)
        stats = st.residue_stats(atom, sasa)
        seg = st.segment(atom, stats, pos, pos + ident - 1)
        rows.append(dict(seed_id=row.seed_id, family=row.family, acc=acc,
                         pos=pos, n_ident=ident,
                         plddt=seg["plddt_mean"], rsasa=seg["rsasa_mean"],
                         call=st.classify(seg)))
    except Exception as exc:  # noqa: BLE001
        rows.append(dict(seed_id=row.seed_id, family=row.family, acc=acc,
                         pos=pos, n_ident=ident,
                         plddt=np.nan, rsasa=np.nan, call=f"error:{type(exc).__name__}"))
    if i % 10 == 0:
        print(f"  {i}/{len(sample)}")

frame = pd.DataFrame(rows)
frame.to_csv("reports/probe_structure_null.tsv", sep="\t", index=False)

ok = frame.dropna(subset=["plddt", "rsasa"])
folded = ok[(ok["plddt"] >= 70) & (ok["rsasa"] >= 0.25)]
print()
print(f"null matches modelled: {len(ok)}/{len(frame)}")
for family, part in ok.groupby("family"):
    f = part[(part["plddt"] >= 70) & (part["rsasa"] >= 0.25)]
    print(f"  {family}: n={len(part):3d}  folded+exposed {len(f):3d} ({len(f)/len(part):.1%})  "
          f"median pLDDT {part['plddt'].median():.1f}  median rSASA {part['rsasa'].median():.2f}")
print(f"  all: {len(folded)}/{len(ok)} = {len(folded)/max(1,len(ok)):.1%}")
print()
print(ok["call"].value_counts().to_string())
