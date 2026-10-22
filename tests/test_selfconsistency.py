"""The report must not contain a number that the tables do not contain.

This is the test that makes the report auditable: every headline figure is
re-derived from a delivered table and compared with what was written into
REPORT.md and summary.json.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def summary(reports):
    path = reports / "summary.json"
    if not path.exists():
        pytest.skip("pipeline has not been run")
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def report_text(reports):
    path = reports / "REPORT.md"
    if not path.exists():
        pytest.skip("pipeline has not been run")
    return path.read_text(encoding="utf-8")


def test_candidate_counts_match_the_table(reports, summary):
    ranked = pd.read_csv(reports / "40_candidates_ranked.tsv", sep="\t")
    assert summary["n_matches"] == len(ranked)
    assert summary["n_candidate_proteins"] == ranked["acc"].nunique()


def test_seed_counts_match_the_seed_table(reports, summary):
    seeds = pd.read_csv(reports / "10_seeds.tsv", sep="\t", dtype=str).fillna("")
    counts = seeds["family"].value_counts().to_dict()
    assert summary["n_real_seeds"] == counts.get("real", 0)
    assert summary["n_null_A"] == counts.get("N_A", 0)
    assert summary["n_null_B"] == counts.get("N_B", 0)


def test_every_headline_number_appears_in_the_report(report_text, summary):
    for key in ("n_matches", "n_candidate_proteins", "n_real_seeds"):
        assert str(summary[key]) in report_text, f"{key} missing from REPORT.md"


def test_calibration_best_score_is_recomputable(reports, summary):
    seed_scores = pd.read_csv(reports / "40_seed_scores.tsv", sep="\t")
    real = seed_scores[seed_scores["family"] == "real"]["best_score"].dropna()
    assert summary["best_real_score"] == int(real.max())


def test_composition_block_is_recomputable(cfg, summary):
    """Recompute the Gly/Ala fraction of the block from the raw query FASTA."""
    from phmimic.core import read_fasta

    sequence = read_fasta(cfg.path_of("data.query"))[0][2]
    start, end = summary["composition"]["block"]
    block = sequence[start - 1 : end]
    fraction = sum(1 for aa in block if aa in "GA") / len(block)
    assert summary["composition"]["block_length"] == len(block)
    assert abs(summary["composition"]["block_ga_fraction"] - fraction) < 1e-6


def test_config_hash_in_summary_matches_config_file(cfg, summary):
    from phmimic.core import sha256_text

    assert summary["config_sha256"] == sha256_text(
        cfg.path.read_text(encoding="utf-8")
    )


@pytest.fixture(scope="module")
def structure(reports):
    path = reports / "65_structure_summary.json"
    if not path.exists():
        pytest.skip("structure stage has not been run")
    return json.loads(path.read_text(encoding="utf-8"))


def test_structure_counts_match_the_tables(reports, structure):
    """Every structure number must be re-derivable from the two delivered tables."""
    real = pd.read_csv(reports / "65_structure_real.tsv", sep="\t")
    assert structure["n_candidates"] == len(real)
    assert structure["n_folded_exposed"] == int(real["folded_exposed"].sum())
    assert structure["calls"] == real["call"].value_counts().to_dict()

    null_path = reports / "65_structure_null.tsv"
    if null_path.exists():
        null = pd.read_csv(null_path, sep="\t")
        cal = structure["calibration"]
        assert cal["n_null"] == len(null)
        assert cal["n_null_folded_exposed"] == int(null["folded_exposed"].sum())
        expected = (1 + cal["n_null_folded_exposed"]) / (1 + cal["n_null"])
        assert cal["p_unconditional"] == pytest.approx(expected, abs=5e-5)


def test_structure_call_agrees_with_its_own_thresholds(reports, structure):
    """A segment labelled folded+exposed must pass the published thresholds."""
    params = structure["parameters"]
    real = pd.read_csv(reports / "65_structure_real.tsv", sep="\t")
    for row in real.itertuples(index=False):
        if row.call == "not_modelled":
            continue
        passes = (row.plddt >= params["folded_plddt"]
                  and row.rsasa >= params["exposed_rsasa"])
        assert bool(row.folded_exposed) == bool(passes), row.gene
        expected = "surface_of_folded_domain" if passes else None
        if expected:
            assert row.call == expected, row.gene
