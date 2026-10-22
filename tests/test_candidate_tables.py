"""Invariants of the delivered tables.

These tests exist because the earlier rounds failed in a way that no unit test
caught: the tables were fine, the numbers were fine, and the reading was wrong.
So what is pinned here is the shape of the deliverable -- that the ordering is
numeric, that no fused score exists, that the columns which exist to stop
over-reading are populated, and that every headline number in a summary can be
recomputed from the table it summarises.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest


def _read(reports, name):
    path = reports / name
    if not path.exists():
        pytest.skip(f"{name} not produced in this workspace")
    return pd.read_csv(path, sep="\t", dtype=str)


def _json(reports, name):
    path = reports / name
    if not path.exists():
        pytest.skip(f"{name} not produced in this workspace")
    return json.loads(path.read_text(encoding="utf-8"))


def test_candidate_table_has_no_composite_score(reports):
    table = _read(reports, "55_candidates_final.tsv")
    forbidden = {"score", "total", "total_score", "priority", "rank_score",
                 "composite", "composite_score"}
    assert not (forbidden & {c.lower() for c in table.columns}), (
        "axes must stay stacked; a fused score would hide which axis carried it")


def test_candidate_table_is_sorted_numerically(reports):
    table = _read(reports, "55_candidates_final.tsv")
    z = pd.to_numeric(table["margin_z"], errors="coerce")
    assert z.is_monotonic_decreasing, "margin_z must descend numerically, not lexically"
    assert z.notna().all()


def test_candidate_table_states_what_is_unknown(reports):
    table = _read(reports, "55_candidates_final.tsv")
    assert table["unknowns"].fillna("").str.len().gt(0).all()
    assert table["manual_judgement"].fillna("").str.len().gt(0).all()


def test_unqueried_iedb_is_not_recorded_as_zero(reports):
    table = _read(reports, "55_candidates_final.tsv")
    allowed = {"yes", "no", "not_queried", "failed"}
    assert set(table["iedb_queried"].unique()) <= allowed
    blank = table[table["iedb_queried"] == "not_queried"]["iedb_bcell_records"]
    assert (blank == "not_queried").all(), (
        "an unqueried protein must not look like a queried protein with 0 records")


def test_maturation_summary_matches_its_table(reports):
    table = _read(reports, "51_processing_state.tsv")
    summary = _json(reports, "66_maturation_summary.json")
    assert int(summary["n_classified"]) == len(table)
    calls = table["processing_call"].value_counts().to_dict()
    assert int(summary["calls"]["retained_in_chain"]) == calls.get("retained_in_chain", 0)
    # An unannotated span is not a confirmed one; the two counts stay separate.
    n_confirmed = int(summary["n_confirmed"])
    n_unannotated = int(summary["n_unannotated"])
    assert calls.get("retained_unannotated", 0) == n_unannotated
    assert n_confirmed + n_unannotated <= len(table)


def test_layered_summary_matches_its_table(reports):
    table = _read(reports, "53_layered_table.tsv")
    summary = _json(reports, "68_layered_summary.json")
    assert int(summary["n_rows"]) == len(table)
    assert int(summary["comp_accessible"].get("yes", 0)) == int(
        (table["comp_accessible"] == "yes").sum())
    assert int(summary["both_folded_exposed"].get("yes", 0)) == int(
        (table["both_folded_exposed"] == "yes").sum())
    hits = table["overlaps_published_query_region"].fillna("")
    assert int(summary["n_overlapping_published_query"]) == int((hits != "").sum())


def test_overlap_columns_are_named_for_their_side(reports):
    table = _read(reports, "53_layered_table.tsv")
    assert {"overlaps_published_target_region",
            "overlaps_published_query_region"} <= set(table.columns)
    assert "overlaps_published_epitope" not in table.columns, (
        "an unqualified name hid which protein the overlap was measured on")


def test_benchmark_audit_is_a_separate_table(reports):
    audit = _read(reports, "71_benchmark_audit.tsv")
    assert len(audit) == 5
    assert set(audit["verdict"]) <= {"recovered", "not_recovered"}
    table = _read(reports, "55_candidates_final.tsv")
    assert "verdict" not in table.columns, (
        "the sensitivity audit must not live in the candidate table, where it "
        "would look like an input to the ranking")


def test_strategy_of_negative_results_is_recorded(reports):
    summary = _json(reports, "61_mmseqs_summary.json")
    assert "n_evalue_below_1" in summary
    assert "ga_block_overlap" in summary
    note = summary["note"].lower()
    assert "not been audited" in note or "unaudited" in note or "未审计" in note, (
        "the low-complexity attribution must stay flagged as unaudited")


def test_selection_null_reports_its_resolution_floor(reports):
    summary = _json(reports, "63_selection_summary.json")
    n_shuffles = int(summary["n_shuffles"])
    assert float(summary["resolution_floor"]) == pytest.approx(
        1.0 / (1 + n_shuffles), rel=1e-3)
    # A p value sitting on the floor is not evidence of anything sharper; the
    # summary has to say so itself rather than leave it to a reader's memory.
    if float(summary["p_selection_n_align"]) <= float(summary["resolution_floor"]) * 1.001:
        assert summary.get("note")


def test_composite_calibration_reports_both_max_and_rate(reports):
    summary = _json(reports, "64_composite_summary.json")
    cal = summary["calibration"]
    assert {"max_score", "frac_ge_min_score"} <= set(cal)
    assert cal["max_score"]["real"] <= cal["max_score"]["null_max"] or \
        float(cal["max_score"]["p"]) < 1.0
    note = (summary.get("note") or "").lower()
    assert "fraction" in note or "比例" in note, (
        "the count columns must be flagged as fractions, since the real query "
        "reports roughly 18x more alignments than a shuffle")


def test_published_entries_carry_both_sides(reports, cfg):
    for entry in cfg.get("published", []):
        assert "region" in entry
        assert "ebna1_region" in entry, (
            "the target-side span alone cannot say whether a candidate landed on "
            "a validated EBNA1 stretch")
        assert "reachable" in entry
