# Changelog

## 1.0.0 — 2026-10-22

First public version. Data pinned to UniProt release 2026_03 (2026-09-02);
human proteome restricted to reviewed entries (20,431).

**Pipeline.** 21 stages, `bash run.sh`, ~375 s with warm caches. Every input is
hash-verified against `data/MANIFEST.tsv` before use, and every stage records the
SHA256 of each table it writes into `reports/run_manifest.json`.

**Two channels, deliberately separate.**

- Anchored (`s10`–`s40`): seeds taken from the antibody-validated EBNA1 386-440
  region, plus two generated null seed families.
- Unbiased (`s61`–`s73`): all 641 EBNA1 residues against the whole proteome, no
  region excluded, every alignment given a both-sides composition null.

**Calibration.** Four levels of null — whole-proteome composition shuffle,
windowed shuffle at 10/20/40 aa, per-alignment two-sided shuffle, and exact
enumeration where the space is small enough to enumerate. Every empirical p value
is reported next to the resolution floor `1/(1+b)` that produced it.

**Axis validation (`s72`).** Each candidate-ranking axis is judged against the
four published mimics before it is allowed to rank anything. 0 of 8 statistical
axes pass: every one either loses a published positive or retains more than a
quarter of the pool. The candidate list therefore carries a declared,
non-significant order and a per-row list of what is unknown.

**Deliverable.** `docs/results/CANDIDATES.md`: 55 composition-controlled local
matches, of which 3 survive the stricter test that no shuffle of the same
residues reaches the observed score. 46 of the 55 have an EBNA1 window that is
not itself low-complexity.

**Known limits, stated in the documents rather than inferred by the reader.**
Conformational epitopes are invisible to this method; binding is not measured
anywhere; the negative result is not novel (shared pentapeptides cover essentially
every human protein, as shown by Kanduc and colleagues in 2010) — what this
repository adds is the composition control, applied consistently.

**Engineering.** 75 tests, fresh-clone test run green (30 passed / 45 skipped
without data), byte-identical output across two full runs, and a CI job that
builds the C scanners, runs the tests and checks that every stage module is
registered.
