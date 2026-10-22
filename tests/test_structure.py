"""Unit tests for the structure channel.

No network: the PDB parser and the SASA calculation are fed a hand-written
three-residue model, so the arithmetic is checked without depending on
AlphaFold being reachable.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phmimic import structure as st

# Two alanines and a glycine; the alanines sit 1.4 A apart so they bury each
# other, the glycine is 30 A away and must come out fully exposed.
PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.460   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       3.200   1.400   0.000  1.00 90.00           O
ATOM      5  CB  ALA A   1       0.700  -0.700   1.000  1.00 90.00           C
ATOM      6  N   ALA A   2       1.400   1.400   0.000  1.00 40.00           N
ATOM      7  CA  ALA A   2       1.460   2.860   0.000  1.00 40.00           C
ATOM      8  C   ALA A   2       2.900   3.200   0.000  1.00 40.00           C
ATOM      9  N   GLY A   3      30.000   0.000   0.000  1.00 95.00           N
ATOM     10  CA  GLY A   3      31.460   0.000   0.000  1.00 95.00           C
ATOM     11  C   GLY A   3      32.000   1.400   0.000  1.00 95.00           C
ENDMDL
ATOM     12  N   ALA A   9       0.000   0.000   0.000  1.00 10.00           N
"""


@pytest.fixture()
def atom_file(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text(PDB)
    return path


def test_parse_pdb_stops_at_first_model(atom_file):
    atom = st.parse_pdb(atom_file)
    assert len(atom["xyz"]) == 11, "the atom after ENDMDL must be ignored"
    assert atom["resnum"].tolist() == [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3]
    assert atom["resname"][8] == "G", "residue names are kept as one-letter codes"


def test_parse_pdb_reads_plddt_from_bfactor(atom_file):
    atom = st.parse_pdb(atom_file)
    assert atom["plddt"][0] == pytest.approx(90.0)
    assert atom["plddt"][8] == pytest.approx(95.0)


def test_isolated_atom_has_full_sphere_sasa():
    """A lone atom exposes 4*pi*(r+probe)^2 of accessible surface."""
    atom = dict(xyz=np.zeros((1, 3)), element=["C"], resnum=np.array([1]),
                resname=np.array(["ALA"]), plddt=np.array([90.0]))
    area = st.sasa_per_atom(atom, probe=1.4, n_points=920)[0]
    expected = 4 * math.pi * (st.VDW["C"] + 1.4) ** 2
    assert area == pytest.approx(expected, rel=0.02)


def test_burial_reduces_sasa(atom_file):
    atom = st.parse_pdb(atom_file)
    area = st.sasa_per_atom(atom, probe=1.4, n_points=920)
    gly_total = area[8:].sum()
    ala_total = area[:8].sum()
    # three glycine atoms at 1.44 A spacing are mostly buried; the two alanines
    # are packed too, so both clusters expose far less than an equal number of
    # free atoms would.
    free = 4 * math.pi * (st.VDW["N"] + 1.4) ** 2
    assert gly_total < 2.0 * free
    assert ala_total < 5.0 * free
    assert area[0] > 0


def test_more_sphere_points_converge(atom_file):
    atom = st.parse_pdb(atom_file)
    coarse = st.sasa_per_atom(atom, probe=1.4, n_points=92)
    fine = st.sasa_per_atom(atom, probe=1.4, n_points=920)
    assert coarse.shape == fine.shape
    # the coarse lattice must not be off by more than a few percent per atom
    assert np.allclose(coarse.sum(), fine.sum(), rtol=0.05)


def test_residue_stats_are_relative_to_max_asa(atom_file):
    atom = st.parse_pdb(atom_file)
    sasa = st.sasa_per_atom(atom, probe=1.4, n_points=920)
    stats = st.residue_stats(atom, sasa)
    assert sorted(stats) == [1, 2, 3]
    res1 = stats[1]
    assert res1["resname"] == "A"
    assert res1["rsasa"] == pytest.approx(res1["sasa"] / st.MAX_ASA["A"])
    assert stats[3]["rsasa"] > stats[1]["rsasa"], "the isolated residue is the exposed one"


def test_segment_averages_only_residues_present(atom_file):
    atom = st.parse_pdb(atom_file)
    sasa = st.sasa_per_atom(atom, probe=1.4, n_points=920)
    stats = st.residue_stats(atom, sasa)
    seg = st.segment(atom, stats, 1, 3)
    assert seg["n"] == 3
    assert seg["plddt_mean"] == pytest.approx((90 + 40 + 95) / 3)
    assert seg["plddt_min"] == pytest.approx(40)
    assert st.segment(atom, stats, 5, 8)["n"] == 0


def test_classify_prefers_disorder_over_burial_when_unconfident():
    folded_surface = dict(n=10, plddt_mean=85.0, rsasa_mean=0.40)
    folded_buried = dict(n=10, plddt_mean=85.0, rsasa_mean=0.05)
    unconfident = dict(n=10, plddt_mean=40.0, rsasa_mean=0.70)
    missing = dict(n=0)
    assert st.classify(folded_surface) == "surface_of_folded_domain"
    assert st.classify(folded_buried) == "buried"
    assert st.classify(unconfident) == "disordered_or_loop_exposed"
    assert st.classify(missing) == "not_modelled"


def test_classify_flags_low_confidence_but_buried_separately():
    """A segment with no confident fold and no exposure is neither surface nor
    buried -- saying 'buried' would over-claim."""
    assert st.classify(dict(n=10, plddt_mean=45.0, rsasa_mean=0.12)) == \
        "disordered_low_confidence"
    assert st.classify(dict(n=10, plddt_mean=85.0, rsasa_mean=0.15)) == "partly_buried"
