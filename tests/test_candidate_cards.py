"""Invariants for the deliverable itself (s73).

These are not statistical tests. They exist because the candidate list is the piece
of this project a reader is most likely to quote, and every claim on it has to be
recoverable from a table rather than from a sentence someone wrote once.
"""

import json

import pandas as pd
import pytest

CARDS = "reports/73_candidate_cards.tsv"
SUMMARY = "reports/73_candidate_cards_summary.json"
MARKDOWN = "reports/CANDIDATES.md"
FINAL = "reports/55_candidates_final.tsv"

pytestmark = pytest.mark.skipif(
    not __import__("pathlib").Path(CARDS).exists(),
    reason="run the pipeline first: reports/73_candidate_cards.tsv missing",
)


@pytest.fixture(scope="module")
def cards():
    return pd.read_csv(CARDS, sep="\t")


@pytest.fixture(scope="module")
def summary():
    return json.loads(open(SUMMARY).read())


def test_cards_are_a_subset_of_the_pool(cards):
    """Every card must be a row of the 55, not a separately assembled list."""
    pool = pd.read_csv(FINAL, sep="\t")
    assert len(pool) == 55
    assert set(cards["acc"]) <= set(pool["acc"])
    assert cards["acc"].is_unique


def test_null_tail_flag_matches_its_definition(cards):
    """null_cleared must be exactly null_max < observed score, for every row."""
    for row in cards.to_dict("records"):
        assert bool(row["null_cleared"]) == (
            float(row["null_max"]) < float(row["obs_score"])
        ), row["gene"]


def test_cleared_rows_sort_before_uncleared(cards):
    """The declared ordering, checked where it applies: inside the shortlist a
    pair that cannot be defended never outranks one that can, and the whole file
    keeps shortlisted rows first."""
    assert list(cards["shortlist"]) == sorted(cards["shortlist"], reverse=True)
    shortlist = cards[cards["shortlist"]]
    flags = list(shortlist["null_cleared"])
    assert flags == sorted(flags, reverse=True)
    for _, group in shortlist.groupby("null_cleared"):
        assert list(group["margin_z"]) == sorted(group["margin_z"], reverse=True)


def test_low_complexity_flag_matches_its_definition(cards, cfg):
    ent_min = float(cfg.get("mining.cards.low_complexity_entropy", 3.0))
    maxfrac = float(cfg.get("mining.cards.low_complexity_maxfrac", 0.30))
    for row in cards.to_dict("records"):
        expected = float(row["query_entropy"]) < ent_min or float(
            row["query_maxfrac"]
        ) > maxfrac
        assert bool(row["query_low_complexity"]) == expected, (
            row["gene"], row["query_entropy"], row["query_maxfrac"])


def test_no_card_carries_a_composite_score(cards):
    """The list must not contain a fused score: no axis here passed validation."""
    forbidden = {"score", "priority", "rank_score", "composite", "p_value", "pvalue"}
    assert not (forbidden & set(cards.columns))


def test_every_card_states_its_unknowns(cards):
    assert "unknowns" in cards.columns
    assert cards["unknowns"].astype(str).str.len().min() > 10


def test_summary_counts_match_the_table(cards, summary):
    facts = summary.get("facts", summary)
    assert facts["n_candidates"] == len(cards)
    assert facts["n_shortlist"] == int(cards["shortlist"].sum())
    assert facts["n_normal_composition"] == int((~cards["query_low_complexity"]).sum())
    assert facts["n_query_low_complexity"] == int(cards["query_low_complexity"].sum())
    assert facts["n_null_tail_cleared"] == int(cards["null_cleared"].sum())


def test_markdown_quotes_the_same_numbers(cards, summary):
    """The prose render and the tables must not drift apart."""
    text = open(MARKDOWN).read()
    facts = summary.get("facts", summary)
    assert f"{facts['n_candidates']} candidates" in text or "55" in text
    assert f"**{facts['n_normal_composition']} of {len(cards)}**" in text
    for row in cards.head(3).to_dict("records"):
        assert row["gene"] in text
        assert f"{row['margin_z']:.2f}" in text


def test_low_complexity_subset_is_stated_not_silent(cards):
    """The document must say that the composition filter also removes the
    published linear mimics; otherwise a reader will read the subset as an
    improvement."""
    text = open(MARKDOWN).read()
    if (cards["query_low_complexity"]).any():
        assert "PPRRPP" in text and "RRPFF" in text and "PGAIEQGPAD" in text
        assert "removes them too" in text
