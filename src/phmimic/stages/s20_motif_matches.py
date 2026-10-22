"""s20 -- match every seed (real and null) against the human proteome.

Real seeds and null seeds go through exactly the same code path and land in the
same table. If they did not, any difference between them would be a property of
the analysis rather than of the peptide.
"""
from __future__ import annotations

import pandas as pd

from ..core import read_fasta
from ..scan import Proteome


def run(ctx):
    proteome_records = read_fasta(ctx.cfg.path_of("data.proteome"))
    proteome = Proteome(proteome_records)
    seeds = pd.read_csv(ctx.out("10_seeds.tsv"), sep="\t", dtype=str).fillna("")
    ctx.log(f"    proteome: {len(proteome)} proteins, {len(proteome.mat):,} bytes packed")

    rows = []
    for seed in seeds.itertuples(index=False):
        hits = proteome.find(seed.motif, int(seed.max_mismatch))
        for acc, pos, n_ident, text in hits:
            rows.append(
                dict(
                    seed_id=seed.seed_id,
                    family=seed.family,
                    parent=seed.parent,
                    motif=seed.motif,
                    max_mismatch=int(seed.max_mismatch),
                    acc=acc,
                    tpos=pos,
                    n_ident=n_ident,
                    matched=text,
                )
            )
        if seed.family == "real":
            ctx.log(f"    [{seed.seed_id}] {seed.motif}: {len({h[0] for h in hits})} proteins")

    frame = pd.DataFrame(rows, columns=["seed_id", "family", "parent", "motif",
                                        "max_mismatch", "acc", "tpos", "n_ident", "matched"])
    frame.to_csv(ctx.out("20_motif_matches.tsv"), sep="\t", index=False)

    with_hits = frame["seed_id"].nunique() if len(frame) else 0
    real = frame[frame["family"] == "real"]
    ctx.log(f"    {len(frame)} matches over {with_hits}/{len(seeds)} seeds "
            f"({real['seed_id'].nunique()} real seeds hit)")
    return dict(
        outputs=["reports/20_motif_matches.tsv"],
        facts=dict(n_matches=int(len(frame)), n_seeds_with_hits=int(with_hits),
                   n_seeds=int(len(seeds))),
    )
