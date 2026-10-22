"""The null must be the null it claims to be.

If the permutation family were sampled instead of enumerated, or if a seed were
allowed into its own null, every p-value in the report would change. These are
cheap invariants that catch exactly that.
"""
from __future__ import annotations

import random

from phmimic.seeds import (
    _n_arrangements,
    build_seed_table,
    distinct_arrangements,
    proteome_random_windows,
)


def test_arrangement_counts_are_exact():
    assert _n_arrangements("PPRRPP") == 15          # 6! / (4! 2!)
    assert _n_arrangements("RRPFF") == 30           # 5! / (2! 2!)
    assert _n_arrangements("PGAIEQGPAD") == 453600  # 10! / (2! 2! 2!)


def test_small_spaces_are_enumerated_not_sampled():
    rng = random.Random(42)
    arrangements = distinct_arrangements("PPRRPP", 200, rng)
    assert len(arrangements) == 15
    assert "PPRRPP" in arrangements
    assert all(sorted(a) == sorted("PPRRPP") for a in arrangements)


def test_large_spaces_are_sampled_to_the_limit():
    rng = random.Random(42)
    arrangements = distinct_arrangements("PGAIEQGPAD", 200, rng)
    assert len(arrangements) == 200
    assert len(set(arrangements)) == 200
    assert all(sorted(a) == sorted("PGAIEQGPAD") for a in arrangements)


def test_seed_table_is_deterministic_and_excludes_self(cfg, proteome_records):
    sequences = [r[2] for r in proteome_records]
    first = build_seed_table(cfg, sequences)
    second = build_seed_table(cfg, sequences)
    assert [r["seed_id"] for r in first] == [r["seed_id"] for r in second]
    assert [r["motif"] for r in first] == [r["motif"] for r in second]

    real = {r["seed_id"]: r["motif"] for r in first if r["family"] == "real"}
    for row in first:
        if row["family"] == "N_A":
            assert row["motif"] != real[row["parent"]], "a seed leaked into its own null"


def test_random_windows_respect_length(proteome_records):
    sequences = [r[2] for r in proteome_records]
    windows = proteome_random_windows(sequences, 6, 50, random.Random(42))
    assert len(windows) == 50
    assert all(len(w) == 6 for w in windows)
