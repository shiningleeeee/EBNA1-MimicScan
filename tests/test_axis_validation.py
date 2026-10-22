"""The axis verdicts are a deliverable, so they are pinned like one.

s72 turns the project's central methodological claim into a table: a statistical
axis may order candidates only if it keeps the published positives and stays
selective. If that table could say anything, the claim would be decoration.

What is checked here is the verdict logic itself, recomputed from the table's
own columns, plus the two properties that make the table worth having:
that the declared thresholds are the ones actually applied, and that a
hard mechanistic constraint is never confused with a statistical axis.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest


def _table(reports):
    path = reports / "72_axis_validation.tsv"
    if not path.exists():
        pytest.skip("72_axis_validation.tsv not produced in this workspace")
    return pd.read_csv(path, sep="\t", dtype=str)


def _summary(reports):
    path = reports / "72_axis_validation_summary.json"
    if not path.exists():
        pytest.skip("72_axis_validation_summary.json not produced")
    return json.loads(path.read_text(encoding="utf-8"))


def test_every_axis_has_a_verdict(reports):
    table = _table(reports)
    assert len(table) >= 5
    assert table["verdict"].notna().all()
    assert table["consequence"].notna().all()
    assert set(table["kind"]) <= {"hard_constraint", "ranking_axis", "annotation"}


def test_ranking_verdicts_match_the_declared_rule(reports):
    table = _table(reports)
    summary = _summary(reports)
    declared = summary["declared_before_measuring"]
    min_pos = declared["min_positives_passing"]
    max_keep = declared["max_pool_retained"]

    for row in table.to_dict("records"):
        if row["kind"] != "ranking_axis":
            continue
        n_pass = int(row["positives_passing"])
        kept = row["pool_retained_fraction"]
        kept = None if pd.isna(kept) else float(kept)
        if row["verdict"] == "pass":
            assert n_pass >= min_pos, f"{row['axis']} passed while losing positives"
            assert kept is not None and kept <= max_keep, (
                f"{row['axis']} passed while keeping {kept} of the pool")
        else:
            assert (n_pass < min_pos) or (kept is not None and kept > max_keep), (
                f"{row['axis']} failed without a declared reason")


def test_a_hard_constraint_is_not_a_ranking_axis(reports):
    table = _table(reports)
    hard = table[table["kind"] == "hard_constraint"]
    assert len(hard) >= 1
    for row in hard.to_dict("records"):
        # a mechanistic necessity is either applied or it removes a positive
        assert row["verdict"] in {"pass", "fail: removes a published positive"}
        assert "may not be used" not in row["consequence"]


def test_annotation_axes_never_filter(reports):
    table = _table(reports)
    for row in table[table["kind"] == "annotation"].to_dict("records"):
        assert row["verdict"] == "annotation_only"
        assert "never" in row["consequence"]


def test_summary_counts_agree_with_the_table(reports):
    table = _table(reports)
    summary = _summary(reports)
    ranking = table[table["kind"] == "ranking_axis"]
    assert summary["n_axes"] == len(table)
    assert summary["n_ranking_axes"] == len(ranking)
    assert summary["n_ranking_axes_passed"] == int(
        (ranking["verdict"] == "pass").sum())


def test_positives_table_keeps_both_hard_cases(reports):
    """The positive set must retain the pairs a linear screen cannot see."""
    path = reports / "72_axis_positives.tsv"
    if not path.exists():
        pytest.skip("72_axis_positives.tsv not produced in this workspace")
    pos = pd.read_csv(path, sep="\t", dtype=str)
    assert len(pos) >= 4
    unreachable = pos[pos["reachable"].str.lower() == "false"]
    assert len(unreachable) >= 1, (
        "no reachable=false positive: an axis could be credited for the easy cases only")
    # identity is a number, and gene labels follow their own accession
    assert pd.to_numeric(pos["identity"]).notna().all()
    assert pos["gene"].notna().all()
    assert pos["acc"].nunique() == len(pos)
