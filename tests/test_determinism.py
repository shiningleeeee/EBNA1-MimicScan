"""The scan must be deterministic and must not match across protein boundaries."""
from __future__ import annotations

import hashlib

from phmimic.scan import Proteome


def test_packing_is_deterministic(proteome_records):
    a = Proteome(proteome_records)
    b = Proteome(proteome_records)
    ha = hashlib.sha256(a.mat.tobytes()).hexdigest()
    hb = hashlib.sha256(b.mat.tobytes()).hexdigest()
    assert ha == hb


def test_no_cross_boundary_matches(proteome):
    """A motif built from the tail of one protein and the head of the next must
    never be reported, even with a mismatch budget."""
    first, second = proteome.accs[0], proteome.accs[1]
    tail = proteome.sequence(first)[-3:]
    head = proteome.sequence(second)[:3]
    if "*" in tail + head:
        return
    hits = proteome.find(tail + head, max_mismatch=0)
    assert hits == []


def test_find_is_repeatable(proteome):
    first = proteome.find("PPRRPP", 0)
    second = proteome.find("PPRRPP", 0)
    assert first == second
    assert len(first) > 0


def test_mismatch_budget_is_monotone(proteome):
    exact = proteome.find("RRPFF", 0)
    loose = proteome.find("RRPFF", 1)
    assert len(loose) >= len(exact)
    assert {h[0] for h in exact} <= {h[0] for h in loose}


def test_reported_positions_are_inside_proteins(proteome):
    for acc, pos, n_ident, text in proteome.find("PPRRPP", 0):
        assert 1 <= pos <= proteome.lengths[acc]
        assert pos + len(text) - 1 <= proteome.lengths[acc]
        assert text == "PPRRPP" and n_ident == 6
