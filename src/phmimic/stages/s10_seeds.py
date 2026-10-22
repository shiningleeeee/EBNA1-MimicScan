"""s10 -- build the real seeds and both null seed families.

The seed table is written to disk and hashed, so the null that a candidate is
judged against is itself an artefact with a fixed identity.
"""
from __future__ import annotations

import pandas as pd

from ..core import read_fasta
from ..seeds import build_seed_table, seed_table_sha256


def run(ctx):
    proteome = read_fasta(ctx.cfg.path_of("data.proteome"))
    rows = build_seed_table(ctx.cfg, [r[2] for r in proteome])
    frame = pd.DataFrame(rows)
    path = ctx.out("10_seeds.tsv")
    frame.to_csv(path, sep="\t", index=False)

    by_family = frame["family"].value_counts().to_dict()
    ctx.log(f"    seeds: {len(frame)} total {by_family}")
    for family in ("N_A", "N_B"):
        sub = frame[frame["family"] == family]
        if len(sub):
            ctx.log(f"    {family}: {len(sub)} seeds, e.g. {', '.join(sub['motif'].head(6))}")
    return dict(
        outputs=["reports/10_seeds.tsv"],
        facts=dict(n_seeds=len(frame), by_family=by_family,
                   seed_table_sha256=seed_table_sha256(rows)),
    )
