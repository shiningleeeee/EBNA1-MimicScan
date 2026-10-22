"""The p-value arithmetic, pinned to hand-computed values.

A p-value that is off by an off-by-one is the difference between "worth a
Western blot" and "indistinguishable from noise", so the formula is tested
directly rather than trusted.
"""
from __future__ import annotations

from phmimic.score import empirical_p


def test_counts_and_both_p_values():
    result = empirical_p(4, [None, 3, 4, 2])
    assert result["n_null"] == 4
    assert result["n_null_present"] == 3      # the None seed matched nothing
    assert result["n_null_ge"] == 1           # only the 4
    assert result["p_unconditional"] == 0.4   # (1+1)/(1+4)
    assert result["p_conditional"] == 0.5     # (1+1)/(1+3)
    assert result["null_max"] == 4


def test_observation_above_every_null():
    result = empirical_p(5, [1, 2, 3])
    assert result["n_null_ge"] == 0
    assert result["p_unconditional"] == 0.25  # (1+0)/(1+3)


def test_observation_matching_everything_is_not_significant():
    """A null that always reaches the observed score leaves p at 1, not 0."""
    result = empirical_p(0, [1, 2, 3])
    assert result["n_null_ge"] == 3
    assert result["p_unconditional"] == 1.0


def test_missing_observation_is_not_significant():
    result = empirical_p(None, [2, 3, None])
    assert result["n_null_ge"] == 3
    assert result["p_unconditional"] == 1.0


def test_empty_null_returns_none_not_one():
    result = empirical_p(3, [])
    assert result["n_null"] == 0
    assert result["p_unconditional"] is None
    assert result["p_conditional"] is None
