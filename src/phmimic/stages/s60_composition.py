"""s60 -- the composition profile that explains everything else.

EBNA1 90-330 is a pure Gly/Ala block. Any method that compares composition --
shuffled nulls, short-motif counting, alignment scores -- will put this block's
matches at the top, because the block is 240 residues drawn from a two-letter
alphabet. This stage writes that profile out so the claim is inspectable rather
than asserted.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..core import read_fasta, write_json


def run(ctx):
    query = read_fasta(ctx.cfg.path_of("data.query"))[0]
    seq = query[2]
    window = int(ctx.cfg.get("composition.window", 20))
    block = ctx.cfg.get("composition.gar_block", [90, 330])

    rows = []
    for start in range(0, max(1, len(seq) - window + 1)):
        chunk = seq[start : start + window]
        if len(chunk) < window:
            break
        ga = sum(1 for aa in chunk if aa in "GA") / window
        counts = pd.Series(list(chunk)).value_counts()
        rows.append(dict(
            start=start + 1,
            end=start + window,
            ga_fraction=round(ga, 3),
            n_distinct=int(len(counts)),
            top_residue=counts.index[0],
            top_fraction=round(float(counts.iloc[0]) / window, 3),
        ))
    profile = pd.DataFrame(rows)
    profile.to_csv(ctx.out("60_composition_profile.tsv"), sep="\t", index=False)

    inner = profile[(profile["start"] >= block[0]) & (profile["end"] <= block[1])]
    block_seq = seq[block[0] - 1 : block[1]]
    facts = dict(
        query=query[0],
        query_length=len(seq),
        block=block,
        block_length=len(block_seq),
        block_ga_fraction=round(
            sum(1 for aa in block_seq if aa in "GA") / max(1, len(block_seq)), 6),
        block_n_distinct=len(set(block_seq)),
        block_non_ga=[dict(pos=block[0] + i - 1, aa=aa)
                      for i, aa in enumerate(block_seq, start=1) if aa not in "GA"],
        window=window,
        windows_fully_inside_block=int(len(inner)),
        block_windows_ga_min=round(float(inner["ga_fraction"].min()), 6) if len(inner) else None,
        block_windows_ga_max=round(float(inner["ga_fraction"].max()), 6) if len(inner) else None,
        profile_max_ga_outside_block=round(
            float(profile[~profile.index.isin(inner.index)]["ga_fraction"].max()), 6),
    )
    write_json(ctx.out("60_composition_summary.json"), facts)
    ctx.log(f"    EBNA1 {len(seq)} aa; block {block[0]}-{block[1]} is "
            f"{facts['block_ga_fraction']:.0%} Gly/Ala over {facts['block_length']} aa "
            f"({facts['block_n_distinct']} distinct residues)")
    return dict(outputs=["reports/60_composition_profile.tsv",
                         "reports/60_composition_summary.json"], facts=facts)
