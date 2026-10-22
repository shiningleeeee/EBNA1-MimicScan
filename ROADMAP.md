# Roadmap

What is left, what it would cost, and — for each item — **which sentence in the
documents it would change**. An item that changes no sentence is not worth doing.

Version 1.0.0 (`config.yaml: version`) marks a *contract*: these data (UniProt
2026_03), this method, these tables. Rerunning the pipeline tomorrow reproduces
1.0.0 byte for byte; changing any input or any axis threshold makes it 1.1.0 and
the frozen tables in `docs/results/` must be regenerated with `make archive`.
Numbers in this file refer to the state at 1.0.0.

## Three things not to do

These are not preferences; each one converts a defensible result into an
indefensible one.

1. **Do not add randomisations to make a p value smaller.** Every empirical p is
   printed next to its resolution floor `1/(1+B)`. Going from B=200 to B=20000
   moves the floor, not the answer — the two-level calibration is already
   saturated at its floor (`p=0.005` = `1/201`).
2. **Do not retune a threshold to recover a known positive.** The positives are
   in the file on purpose so that this is possible to do and visible if done.
   `s72` exists to make "the axis passes" a checkable claim rather than an
   assertion.
3. **Do not treat an uncalibrated co-folding prediction as binding evidence.**
   AF3/Boltz will produce a confident-looking interface for any two sequences.
   Until there is a calibration set of *measured* peptide–Fab affinities, a
   predicted interface is a hypothesis about a hypothesis.

---

## A. Doable now, no new science required

**Status, 2026-10-22.** A1 and A2 are done, with outputs frozen in
`docs/results/`. A3's computational half is done and came back negative; A3's wet
half and all of A4 are blocked on reagents rather than on method.

| item | state | what it changed |
|---|---|---|
| A1 | **done** | `DOSSIER.md` §16 and `METHOD.md` no longer say "cause not determined": composition correction is the mechanism, measured over six cells (`docs/results/60b_mmseqs_param_matrix.tsv`) |
| A2 | **done** | the negative claim now covers 205,766 canonical human entries instead of 20,431 reviewed ones, and got *weaker* rather than stronger (`docs/results/61b_trembl_expansion.tsv`) |
| A3 | half done, **negative** | a fifth statistic — identity at the residues the Fab actually contacts — puts the crystallised binder *inside* the candidate range (`docs/results/62b_paratope_contact_identity.tsv`) |
| A4 | **not started** | wet lab; see the peptide list below |

Three numbers worth carrying, because each one is the kind that gets quoted back:

- **A1.** With `--comp-bias-corr 0` the Gly-Ala block collects 1,254 / 1,696
  alignments at E = 0.0, and *every such cell fails its own windowed null*
  (p = 0.995 / 0.881). The zero overlap in the published run is composition
  correction working, and the uncorrected alternative is exactly the artefact
  this project exists to remove.
- **A2.** A tenfold library multiplies alignments by 7.3 and genes-with-any-
  similarity by 1.6, produces no alignment below E = 1, and makes the real
  EBNA1's standing against its own shuffles worse: p 0.1493 → 0.2736.
- **A3.** At the 10 EBNA1 positions MS39p2w174 actually touches, GlialCAM — a
  verified binder with a crystal structure — shares 5 of 10. Candidates COBL,
  MAPK8IP2 and PPP2R3B share 7/10, 5/6 and 2/3. Contact identity therefore does
  not separate a real binder from an unvalidated match, which is the same verdict
  the sequence statistics earned, now from structural data.

### A1. MMseqs2 low-complexity / composition-correction matrix

**Status: done.** `src/analysis/60_mmseqs_param_matrix.py`,
`docs/results/60b_mmseqs_param_matrix.tsv`.

**Why.** `docs/DOSSIER.md` §16 says the Gly-Ala block is hit by 0 of 363
alignments and that the contributions of (i) MMseqs2's default low-complexity
handling, (ii) `--comp-bias-corr 1`, and (iii) the intrinsic difficulty of
scoring a 240-aa near-homopolymer are **not separated**. That "not separated" is
an admitted hole in an otherwise quantitative document, and it is the most
obvious thing a reader will attack.

**Do.** Four cells, holding everything else fixed: {low-complexity masking on,
off} × {`--comp-bias-corr 1`, `0`}. Separately, mask the *database* side
(`mmseqs makepaddedseqdb`-style target masking) versus the *query* side
(pre-mask EBNA1 90-330 before searching) — these answer different questions and
must not be conflated. **Each cell needs its own null model**; raw hit counts
across cells are not comparable and comparing them is the standard error this
project exists to avoid.

**Cost.** One afternoon. ~20 searches × 6.5 s plus nulls.

**Changes.** §16 "cause not yet determined" → a determined sentence with per-cell
numbers; `61_mmseqs_summary.json`'s `note` field loses its "have not been
audited" clause.

**Acceptance.** A table with one row per cell: hits, best E, E distribution,
fraction overlapping 90-330, and the matched null calibration for that cell.

### A2. TrEMBL expansion

**Why.** The negative result currently covers 20,431 reviewed human proteins. The
first question any reviewer asks is "TrEMBL is 30× larger". Right now the honest
answer is "we did not look".

**Do.** Not to find candidates — to bound the negative claim. Requirements before
starting: record the TrEMBL release, a redundancy-reduction policy (MMseqs2
`cluster` at 90% or similar), an isoform policy, and disk/time budget. **After
expanding, the null models must be re-fitted**, because the search space changed
and every E-value and every null distribution is a function of database size.

**Cost.** 2–3 hours of compute plus the recalibration; ~1 GB disk.

**Done, and differently from the plan above.** `src/analysis/61_trembl_expansion.py`,
`docs/results/61b_trembl_expansion.tsv`. Near-duplicates were **not** clustered:
clustering is itself a modelling choice that changes the search space, so the
library is a union of canonical entries with exact duplicates removed (4,940
dropped, 205,766 kept). Both a tight and a loose `-e` cutoff are reported, because
E-values scale with library size and a single cutoff compares thresholds rather
than contents. The scope clause in `METHOD.md` §8b now reads reviewed + TrEMBL
rather than "reviewed only".

**Changes.** `METHOD.md` §0 scope clause ("reviewed only") → "reviewed +
redundancy-reduced TrEMBL"; §8's "0 of 363" numbers get re-derived.

**Acceptance.** The unbiased channel's per-protein rate is reported against both
databases, and the conclusion either holds or is weakened explicitly.

### A3. The one structural test with a real basis: 7K7R Fv × MAPK8IP2 peptide

**Why.** `MAPK8IP2` (`Q13387`, EBNA1 394-429 ↔ 309-341, 54.0%) is the only
candidate that covers a stretch where two published epitopes sit (`PPRRPP`
394-399, `RRPFF` 402-406), and 7K7R gives the actual paratope that binds there.
This is the one place where "would the same antibody bind both?" is answerable
with a structure rather than a guess.

**Precondition — do this first, it may kill the item.** Verify the antibody is
obtainable: IEDB returns **0** B-cell records for `Q13387`, and the MS39p2w174
reagent's availability was never confirmed. If there is no accessible reagent and
no reliable model, this is a modelling exercise, and see "three things not to do"
item 3.

**Cost.** Half a day modelling; weeks if a reagent must be sourced.

**Result, computational half only.** `src/analysis/62_paratope_contact_identity.py`,
`docs/results/62b_paratope_contact_identity.tsv`. The precondition is not met —
MS39p2w174 is not a catalogue antibody, so the only route is recombinant Fab from
the chains deposited in 7K7R / 7UQC. What could be done without a reagent was to
ask whether candidates keep the residues this Fab actually touches: of the 262
heavy-atom contacts, 10 EBNA1 residues carry them all (R396 56, P395 44, P398 38,
R397 36, then single digits).

| pair | contacts covered | identical | simple | weighted |
|---|---|---|---|---|
| GlialCAM `Q14CZ8` (**verified binder**) | 10 | 5 | 0.500 | 0.668 |
| COBL `O75128` (candidate) | 10 | 7 | 0.700 | 0.741 |
| MAPK8IP2 `Q13387` (candidate) | 6 | 5 | 0.833 | 0.818 |
| PPP2R3B `Q9Y5P8` (candidate) | 3 | 2 | 0.667 | 0.513 |

Random human peptides at the same register average 0.064. So every one of these
is far above random, and the verified binder is **not** the highest — COBL is.
The statistic is therefore real but not discriminating, which makes it the fifth
axis to fail a positive control and the first to fail one built on experimental
structural data rather than on sequence. A candidate scoring above random here is
not evidence of binding; it is evidence that this paratope's contacts can be
satisfied by composition. That verdict is the item's outcome, and it does not
justify buying an antibody.

**Changes.** `CANDIDATES.md` MAPK8IP2 card's `unknowns` field; possibly a
`docs/results/74_*` table with the interface energy comparison.

**Acceptance.** Either a calibrated comparison against the 7K7R/7UQC pair as
positive controls, or a documented decision not to claim anything.

### A4. Wet lab: three peptides

**Why.** Everything in this repository is a statement about sequences and
predicted structures. Nothing has been measured. Three synthetic peptides turn
the three `null tail: cleared` candidates into testable hypotheses.

**Do.** ELISA against MS serum/CSF for `PRKAR1B` 20-43, `RHBDF2` 339-360, and
`MAPK8IP2` 309-341, each with a competition step against the corresponding EBNA1
peptide (564-587, 541-562, 394-429). The competition step is what separates "this
peptide is antigenic" from "the same antibodies recognise both".

**Cost.** ~3 peptides + one ELISA plate; days, not weeks, if serum is at hand.

**Changes.** It would be the first entry in this project that is not a
computation. A positive result raises the corresponding card; a negative result
is equally publishable and is what the calibration in `s72` predicts is more
likely.

**Acceptance.** Pre-register the read-out and the decision rule before running.

---

## B. Needs the field to move first

### B1. A bigger calibration set of validated mimics

The method's weakest joint is that four positives define the standard, all four
are short linear stretches chosen because they were easy to demonstrate, and one
of them (GlialCAM) depends on phosphorylation that the sequence channel cannot
see. `s72`'s two thresholds (`min_positives_passing: 3`,
`max_pool_retained: 0.25`) are calibrated on n=4. With n≈20 the thresholds could
be set by data instead of by judgement — and the axes that fail today might fail
for a defensible reason rather than an anecdotal one.

### B2. PTM-aware matching

GlialCAM 370-389 is bound as a phospho-peptide (SEP376/377 in 7UQC), and the
whole screen treats `S` and `pS` as the same residue. Any future version needs
modification-aware alignment and a null model that shuffles modified positions
consistently. Until then, this class of mimic is invisible *by construction*, and
that limitation should be read as applying to every result here.

### B3. Affinity as the quantity, not similarity

Every published analysis of this system — including the two MD studies — points
at the same gap: EBNA1 386-405 and GlialCAM 370-389 are read by the same paratope
with different affinities, and the difference is carried by conformational
adjustment, not by the alignment. A screen built on predicted binding free energy,
calibrated against measured Kd for a panel of peptide–Fab pairs, would be a
different kind of tool. Nothing in this repository is a step toward it except the
demonstration that the sequence-level tool cannot get there.

---

## C. Known engineering gaps

- **CI does not verify that the frozen tables match the code.** CI builds, tests
  and checks stage registration; it does not run the pipeline, because that needs
  ~193 MB of cached structures. A partial fix that needs no cache: run
  `s00`–`s64` in CI (UniProt FASTA + MMseqs2 = ~14 MB) and diff only the tables
  those stages produce. That would cover `40_*`, `46_*`, `48_*`, `51_*`, `61_*`,
  `62_*`, `63_*`, `64_*` — everything except the AlphaFold, PDB and IEDB layers.
- **`src/phmimic/stages/s45_window_scan.py` does not exist any more** by design
  (it was never registered and never run; it lives in
  `src/analysis/s45_window_scan_never_run.py`). `make check-stages` enforces that
  this stays true.
- **`src/analysis/` is not covered by tests.** Two of those scripts overturned
  conclusions in `docs/history/`, so they are load-bearing for the record even
  though they are not part of the pipeline.
- **`src/c/motif_mine.c` has one dead variable** (`TABS defined but not used`) and
  no test beyond the pipeline's own use of it. It produced the step-3 result, so
  it is load-bearing; it is not clean.

## D. Rerun recipes

**Refresh the data (new UniProt release).** Update the two `source_url`s and the
expected release in `data/MANIFEST.tsv`, run `bash scripts/01_fetch_data.sh`, then
bump `version` in `config.yaml` and `bash run.sh`. `s00` will fail loudly if the
release does not match — that failure is the feature.

**Change the pathogen.** `config.yaml: data.query` points at one FASTA;
`seeds.real`, `published` and `composition.gar_block` are all EBNA1-specific and
must be replaced. Keep at least one positive that a linear method cannot find,
or `s72` tests nothing. Full checklist in `METHOD.md` §9.
