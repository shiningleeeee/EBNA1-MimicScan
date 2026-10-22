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

## Post-1.0.0 analyses — 2026-10-22

Three items from `ROADMAP.md` section A, each an *analysis* rather than a pipeline
change: they add frozen tables under `docs/results/` and edit prose, but no stage
and no axis threshold changes, so the version stays 1.0.0 and the existing frozen
tables remain valid.

**A1 — composition correction, measured.** `src/analysis/60_mmseqs_param_matrix.py`.
Six cells over {database masking on/off} × {`--comp-bias-corr` 1/0}, each with its
own 40-aa windowed null. Turning composition correction off lets the Gly-Ala block
collect 1,254 and 1,696 alignments at E = 0.0, and every such cell fails its own
null (p = 0.995, 0.881). The zero overlap reported in 1.0.0 is that correction
working — not tantan masking, which changes totals without changing overlap.
`DOSSIER.md` §16 and `METHOD.md` now attribute the effect instead of calling the
cause undetermined.

**A2 — the negative claim widened.** `src/analysis/61_trembl_expansion.py`.
205,766 canonical human entries (20,431 reviewed + 185,335 unreviewed after
dropping 4,940 exact duplicates; no near-duplicate clustering, since clustering is
itself a choice that changes the search space). No alignment reaches E < 1; the
best is 23.36. A tenfold library multiplies alignments by 7.3 and
genes-with-any-similarity by 1.6, and the real EBNA1's standing against its own
shuffles gets *worse* (p 0.1493 → 0.2736). Reported at two `-e` cutoffs because a
fixed E cutoff compares thresholds, not contents, once library size changes.

**A3 — a fifth axis, this time from structural data.** `src/analysis/62_paratope_contact_identity.py`.
At the 10 EBNA1 residues that MS39p2w174 actually contacts in 7K7R (262 heavy-atom
contacts), the crystallised binder GlialCAM `Q14CZ8` shares 5 of 10. Candidates
COBL shares 7 of 10, MAPK8IP2 5 of 6, PPP2R3B 2 of 3 — all far above the random
same-register mean of 0.064. The verified binder is not the highest scorer, so
contact identity is real but not discriminating: the fifth axis to fail a positive
control, and the first to fail one built on measured structure rather than
sequence. The wet half (7K7R Fab × MAPK8IP2 peptide) is blocked: MS39p2w174 is not
a catalogue antibody.

