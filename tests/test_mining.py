"""Tests for the helpers the unbiased chain is built on.

These are cheap and they cover the two places where a silent bug would not
raise: composition arithmetic (an entropy that is quietly wrong changes which
window counts as "normal") and m8 ordering (a file whose row order follows
thread scheduling cannot be hash-compared between runs).
"""
from __future__ import annotations

import pandas as pd
import pytest

from phmimic import mining


def test_entropy_of_a_homopolymer_is_zero():
    h, maxfrac = mining.entropy("AAAAAAAAAA")
    assert h == pytest.approx(0.0, abs=1e-9)
    assert maxfrac == pytest.approx(1.0)


def test_entropy_of_an_even_alphabet():
    h, maxfrac = mining.entropy("ACDEFGHIKLMNPQRSTVWY")
    assert h == pytest.approx(4.3219, abs=1e-3)   # log2(20)
    assert maxfrac == pytest.approx(0.05)


def test_entropy_sees_two_residues_as_low_complexity():
    # EBNA1 394-399 in one letter: two residue types, so this must land well
    # under the 2.8 threshold the pipeline uses for "normal".
    h, _ = mining.entropy("PPRRPP")
    assert h < 1.5


def test_region_of_covers_the_declared_bands():
    assert mining.region_of(1, 20) == "N-terminal"
    assert mining.region_of(200, 260) == "Gly/Ala block (90-330)"
    assert mining.region_of(380, 400) == "SPR-rich (370-445)"
    assert mining.region_of(480, 500) == "folded core (476-585)"
    assert mining.region_of(600, 620) == "C-terminal"


def test_region_of_is_a_cascade_not_a_majority_vote():
    # Gly/Ala is decided by overlap, SPR-rich and folded core by start position.
    # Both halves are pinned: the reported region counts depend on this mix, not
    # on anyone's intuition about where an alignment "mostly" is.
    assert mining.region_of(300, 500) == "Gly/Ala block (90-330)"   # touches block
    assert mining.region_of(390, 500) == "SPR-rich (370-445)"       # starts in band
    assert mining.region_of(331, 400) == "C-terminal"               # starts in neither


def test_overlaps_is_inclusive_and_symmetric():
    assert mining.overlaps(1, 10, 10, 20)
    assert mining.overlaps(10, 20, 1, 10)
    assert not mining.overlaps(1, 9, 10, 20)


def test_flag_matches_tokens_case_insensitively():
    assert mining.flag("Secreted, Cell Membrane", mining.SURFACE_TOKENS)
    assert not mining.flag("Nucleus", mining.SURFACE_TOKENS)


def test_m8_roundtrip_preserves_values(tmp_path):
    frame = pd.DataFrame([
        dict(target="B", evalue=1e-5, bits=40.0, qstart=1, qend=10,
             tstart=5, tend=14, alnlen=10, pident=90.0, qaln="A" * 10,
             taln="A" * 10),
        dict(target="A", evalue=3.0, bits=20.0, qstart=2, qend=11,
             tstart=7, tend=16, alnlen=10, pident=60.0, qaln="C" * 10,
             taln="C" * 10),
    ])
    path = tmp_path / "x.m8"
    mining.write_m8(frame, path)
    back = mining.read_m8(path, require_query_column=False)
    assert list(back.columns) == mining.M8_COLUMNS
    assert len(back) == 2
    # read_m8 orders by evalue, so the best hit comes first regardless of the
    # order the rows were written in.
    assert list(back["target"]) == ["B", "A"]
    assert back.iloc[0]["evalue"] == pytest.approx(1e-5)
    assert back.iloc[1]["evalue"] == pytest.approx(3.0)


def test_read_m8_is_order_insensitive(tmp_path):
    """Row order must depend on the data, never on which thread finished first.

    `write_m8` is a plain writer and keeps whatever order it is given; the
    guarantee lives in `read_m8`, which sorts, so two files holding the same
    alignments in different orders produce identical frames.
    """
    rows = [
        dict(target=t, evalue=e, bits=1.0, qstart=1, qend=5, tstart=1, tend=5,
             alnlen=5, pident=50.0, qaln="AAAAA", taln="AAAAA")
        for t, e in [("C", 5.0), ("A", 0.5), ("B", 2.0)]
    ]
    forward = tmp_path / "fwd.m8"
    reverse = tmp_path / "rev.m8"
    mining.write_m8(pd.DataFrame(rows), forward)
    mining.write_m8(pd.DataFrame(list(reversed(rows))), reverse)
    assert mining.read_m8(forward).equals(mining.read_m8(reverse))
    assert list(mining.read_m8(forward)["target"]) == ["A", "B", "C"]


def test_read_m8_requires_query_column_when_asked(tmp_path):
    path = tmp_path / "q.m8"
    pd.DataFrame([dict(target="A", evalue=1.0, bits=1.0, qstart=1, qend=5,
                       tstart=1, tend=5, alnlen=5, pident=50.0, qaln="AAAAA",
                       taln="AAAAA")]).to_csv(path, sep="\t", index=False,
                                            header=False)
    with pytest.raises(ValueError):
        mining.read_m8(path, require_query_column=True)
