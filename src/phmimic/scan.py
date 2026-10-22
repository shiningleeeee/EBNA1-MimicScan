"""Exact / near-exact motif matching over a reference proteome.

One flat byte array holds every protein, separated by `*`. Any window that
reaches across a boundary necessarily contains a `*`, and every window is
additionally required to lie inside a single protein, so no cross-protein match
can survive.

Mismatch counting is accumulated one column at a time rather than materialised
as a (positions x length) boolean matrix -- at 10.8M positions that matrix would
be >100 MB per seed and the pipeline scans several hundred seeds.
"""
from __future__ import annotations

import numpy as np


class Proteome:
    """Byte-packed concatenation of all sequences, with per-protein offsets."""

    def __init__(self, records: list[tuple[str, str, str]]):
        self.accs = [r[0] for r in records]
        self.genes = {r[0]: r[1] for r in records}
        self.lengths = {r[0]: len(r[2]) for r in records}
        self.offsets = np.zeros(len(records) + 1, dtype=np.int64)
        blobs: list[str] = []
        for i, rec in enumerate(records):
            self.offsets[i + 1] = self.offsets[i] + len(rec[2]) + 1
            blobs.append(rec[2])
            blobs.append("*")
        self.mat = np.frombuffer("".join(blobs).encode("ascii"), dtype=np.uint8)

    def __len__(self) -> int:
        return len(self.accs)

    def sequence(self, acc: str) -> str:
        i = self.accs.index(acc)
        start = int(self.offsets[i])
        return self.mat[start : start + self.lengths[acc]].tobytes().decode()

    def find(self, motif: str, max_mismatch: int = 0) -> list[tuple[str, int, int, str]]:
        """Return [(accession, 1-based target position, n_identical, matched_text)]."""
        L = len(motif)
        n = len(self.mat) - L + 1
        if n <= 0:
            return []
        ref = np.frombuffer(motif.encode("ascii"), dtype=np.uint8)
        window = np.lib.stride_tricks.sliding_window_view(self.mat, L)
        mismatch = np.zeros(n, dtype=np.int16)
        for k in range(L):
            mismatch += window[:, k] != ref[k]
        idx = np.flatnonzero(mismatch <= max_mismatch)
        if idx.size == 0:
            return []
        mm = mismatch[idx]
        acc_i = np.searchsorted(self.offsets, idx, side="right") - 1
        # window must end before the '*' separator that follows its protein
        inside = self.offsets[acc_i + 1] - 1 - idx >= L
        idx, acc_i, mm = idx[inside], acc_i[inside], mm[inside]
        out = []
        for i, a, mismatches in zip(idx.tolist(), acc_i.tolist(), mm.tolist()):
            text = self.mat[i : i + L].tobytes().decode("ascii")
            out.append((self.accs[a], int(i - self.offsets[a]) + 1, L - int(mismatches), text))
        return out


def cns_tokens() -> tuple[str, ...]:
    return (
        "brain", "neuron", "neuronal", "synap", "axon", "glial", "myelin",
        "oligodendro", "astrocyt", "microglia", "schwann", "cerebell",
        "cerebral", "spinal", "central nervous", "retina", "nerve",
    )


SURFACE_TOKENS = (
    "cell membrane",
    "secreted",
    "cell surface",
    "extracellular",
    "gpi-anchor",
)
