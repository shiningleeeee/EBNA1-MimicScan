"""Seed definition and null-seed generation.

The seed is the only experimentally anchored object in this problem, so the
null must be built around *it*: what would a peptide with the same residues in a
different order -- or a random peptide of the same length -- have found?

Two null families, both reported:

  N_A  every distinct arrangement of the seed's own residue multiset.
       Exhaustive when the multiset has at most `enumerate_max` arrangements,
       otherwise a fixed-seed sample. Composition is held exactly constant, so
       any difference between the seed and N_A is attributable to ORDER alone.

  N_B  random windows drawn from the human proteome. Composition is free, so
       N_B answers "what does an arbitrary hexamer find".
"""
from __future__ import annotations

import itertools
import random

import numpy as np


def distinct_arrangements(motif: str, limit: int, rng: random.Random) -> list[str]:
    """All distinct strings with the same multiset of residues as `motif`.

    Enumerated exhaustively when the multiset is small enough (this is the case
    for every 5/6-mer seed here, so there is no sampling variance to argue
    about); otherwise `limit` draws without replacement.
    """
    total = _n_arrangements(motif)
    if total <= limit:
        return sorted({"".join(perm) for perm in set(itertools.permutations(motif))})
    seen: set[str] = set()
    residues = list(motif)
    guard = 0
    while len(seen) < limit and guard < limit * 200:
        rng.shuffle(residues)
        seen.add("".join(residues))
        guard += 1
    return sorted(seen)


def _counts(motif: str) -> list[tuple[str, int]]:
    order: list[str] = []
    seen: dict[str, int] = {}
    for ch in motif:
        if ch not in seen:
            seen[ch] = 0
            order.append(ch)
        seen[ch] += 1
    return [(ch, seen[ch]) for ch in order]


def _n_choose_k(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    out = 1
    for i in range(k):
        out = out * (n - i) // (i + 1)
    return out


def _n_arrangements(motif: str) -> int:
    import math

    total = math.factorial(len(motif))
    for _, count in _counts(motif):
        total //= math.factorial(count)
    return total


def proteome_random_windows(
    sequences: list[str], length: int, n: int, rng: random.Random
) -> list[str]:
    """`n` distinct length-`length` windows sampled uniformly from the proteome."""
    usable = [s for s in sequences if len(s) >= length]
    if not usable:
        return []
    picked: set[str] = set()
    guard = 0
    while len(picked) < n and guard < n * 100:
        seq = usable[rng.randrange(len(usable))]
        start = rng.randrange(len(seq) - length + 1)
        picked.add(seq[start : start + length])
        guard += 1
    return sorted(picked)


def build_seed_table(cfg, sequences: list[str]) -> list[dict]:
    """Real seeds plus both null families, as flat rows ready for scanning."""
    rng = random.Random(cfg.get("seeds.random_state", 42))
    rows: list[dict] = []

    for spec in cfg.get("seeds.real", []):
        rows.append(
            dict(
                seed_id=spec["id"],
                family="real",
                parent="",
                motif=spec["motif"],
                length=len(spec["motif"]),
                max_mismatch=int(spec.get("max_mismatch", 0)),
                ebna1_start=spec["ebna1"][0],
                ebna1_end=spec["ebna1"][1],
                provenance=spec.get("provenance", ""),
            )
        )

    if cfg.get("seeds.null.shuffle_real", True):
        limit = int(cfg.get("seeds.null.enumerate_max", 200))
        for spec in cfg.get("seeds.real", []):
            arrangements = distinct_arrangements(spec["motif"], limit, rng)
            for i, motif in enumerate(arrangements):
                if motif == spec["motif"]:
                    continue  # the seed itself is not part of its own null
                rows.append(
                    dict(
                        seed_id=f"N_A.{spec['id']}.{i:03d}",
                        family="N_A",
                        parent=spec["id"],
                        motif=motif,
                        length=len(motif),
                        max_mismatch=int(spec.get("max_mismatch", 0)),
                        ebna1_start="",
                        ebna1_end="",
                        provenance=f"permutation of {spec['motif']}",
                    )
                )

    n_random = int(cfg.get("seeds.null.proteome_random_n", 200))
    if n_random:
        length = int(cfg.get("seeds.null.proteome_random_len", 6))
        for i, motif in enumerate(
            proteome_random_windows(sequences, length, n_random, rng)
        ):
            rows.append(
                dict(
                    seed_id=f"N_B.{i:03d}",
                    family="N_B",
                    parent="",
                    motif=motif,
                    length=length,
                    max_mismatch=0,
                    ebna1_start="",
                    ebna1_end="",
                    provenance="random proteome window",
                )
            )
    return rows


def seed_table_sha256(rows: list[dict]) -> str:
    import hashlib

    blob = "\n".join(f"{r['seed_id']}\t{r['motif']}\t{r['max_mismatch']}" for r in rows)
    return hashlib.sha256(blob.encode()).hexdigest()
