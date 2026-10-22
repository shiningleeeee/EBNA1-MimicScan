# Handoff — state of this project, and what to do next

Written 2026-10-22, at version 1.0.0. This file is the entry point for anyone
(including a future session of the agent that built it) picking the project up
cold. It is deliberately short and points at the documents that are not.

Read in this order: **this file** → `docs/results/CANDIDATES.md` (the deliverable)
→ `METHOD.md` (how any claim here is allowed to be made) → `docs/DOSSIER.md`
(everything, including the errors) → `ROADMAP.md` (what is left).

## What this repository is

A reproducible screen for **pathogen-to-human molecular mimicry**, exercised on
Epstein–Barr nuclear antigen 1 (EBNA1, `P03211`) against the human proteome. It
exists to answer one question honestly: *given that short peptide sharing between
unrelated proteins is nearly universal, is there anything left to see in EBNA1's
similarity to human proteins once composition is controlled for?*

The answer this repository delivers is **no, not by sequence, and not yet by
structure** — and the more useful half of that answer is *why*: every statistic
that looked selective turned out to also throw away the four published,
antibody-validated EBNA1 mimics. That failure is the finding, and it is measured,
not asserted.

## Short feature plate

- 21 stages, one entry point (`bash run.sh`), ~375 s with warm caches.
- 75 tests. Every stage module must be registered or `make check-stages` fails.
- Byte-identical outputs across two full runs; the only deliberately
  non-deterministic field is the IEDB fetch timestamp, and `run_manifest.json`
  records where that leaks into a hash.
- Data pinned to UniProt release **2026_03** (2026-09-02) and hash-verified
  against `data/MANIFEST.tsv` before any stage runs.
- Four levels of null model, each empirical p printed next to the resolution
  floor `1/(1+b)` that produced it.

## Where the science stands

Three claims survive scrutiny; everything else in the history folder is either
superseded or withdrawn.

1. **The published EBNA1 mimics are not special by sequence.** The one pair with
   a crystal structure (EBNA1 386-405 ↔ GlialCAM 370-389, PDB 7K7R / 7UQC, same
   Fab, Fv superposition 0.63 Å) sits at `margin_z` **0.865**, p = **0.20**
   against a composition null; over 10⁶ shuffles the null reaches 69 where the
   observation is 35. CRYAB and MEF2B score exactly what their best shuffle
   scores. The three verified mimics rank at the 31st, 38th and 15th percentile
   of the 363-alignment pool.
2. **The candidate list is a hypothesis set, not a hit list.** 55 local matches
   survive composition and maturation filters; 3 of them have a null tail that is
   actually cleared (PRKAR1B `P31321`, RHBDF2 `Q6PJF5`, MAPK8IP2 `Q13387`). The
   top of the list is not special: in 200 whole-proteome shuffles, 27 reach the
   same maximum score (p = 0.139). What *is* in excess is quantity, not quality
   (46.8% vs 20.2% of alignments passing the score cut, p = 0.025).
3. **The axes that looked sensible are the ones that fail.** 0 of 8 statistical
   ranking axes pass a positive control (`docs/results/72_axis_validation.tsv`);
   so do a structural axis and a paratope-contact axis (`62b`). No axis is
   allowed to rank anything until it has been run against the published positives
   — that rule is `s72`, and it is the one piece of machinery in here worth
   reusing elsewhere.

## What is done, and what is open

| item | state |
|---|---|
| A1 composition correction, measured (6 cells × own nulls) | **done**, `docs/results/60b_*` |
| A2 negative claim widened to 205,766 canonical human entries | **done**, `docs/results/61b_*` |
| A3 paratope-contact identity, computational half | **done, negative**, `docs/results/62b_*` |
| A3 wet half (7K7R Fab × MAPK8IP2 peptide) | **blocked**: MS39p2w174 is not a catalogue antibody; needs recombinant Fab from the deposited chains |
| A4 wet lab (PRKAR1B 20-43, RHBDF2 339-360, MAPK8IP2 309-341, each with EBNA1 peptide competition) | **not started**, needs a lab |
| B: larger positive-control set, PTM-aware epitopes, affinity as the discriminating quantity | **needs field progress**, not compute |
| C: CI does not verify the frozen tables; `src/analysis/` has no tests; `motif_mine.c` has a dead variable | **open engineering**, none of it affects a claim |

If you are a future session asked to "keep improving this": the honest options are
A3-wet, A4, or item B. **Not** more randomisation, not a threshold retune, and not
an uncalibrated co-folding prediction — see the three prohibitions at the top of
`ROADMAP.md`.

## Resume from nothing

```sh
git clone <this repo> && cd EBNA1-MimicScan
python3 -m venv env/venv && env/venv/bin/pip install -r requirements.txt
make all                      # three C scanners, no mmseqs2 needed
bash scripts/01_fetch_data.sh # UniProt 2026_03, hash-checked, ~15 MB
# mmseqs2: static avx2 build, commit fcf516ab…, place at env/bin/mmseqs
bash run.sh                   # 21 stages; first run downloads ~350 MB of
                              # AlphaFold/PDB/IEDB caches, ~40 min cold
make test                     # 75 passed with data, 30 passed / 45 skipped without
make archive                  # refresh docs/results/ from reports/
```

The previous working copy lived on a disposable GPU pod at
`/workspace/PH-Mimic`. Nothing in it is unique: `data/cache/` re-downloads and
every table that matters is committed under `docs/results/`. Pods get released;
this repository does not.

## Traps, so they are not re-entered

- **Do not read a smaller p as progress.** The two-level calibration is saturated
  at its floor (`p = 0.005` *is* `1/201`). Raising B moves the floor, not the
  answer.
- **The four published positives are in `config.yaml` on purpose**, so that
  retuning a threshold to recover them is possible and visible. Don't.
- **Two claims in the history folder are withdrawn**, and are kept rather than
  deleted: that the antibody interface concentrates on one arginine (it was a
  dict-key collision counting half a paratope — 100 contacts instead of 262), and
  that all 363 unbiased alignments fall in EBNA1 350-641 (seven start before 90;
  the correct statement is that none overlaps the Gly-Ala block).
- **A negative result here is not novel and does not claim to be.** Shared
  pentapeptides cover essentially every human protein (Kanduc and colleagues,
  2010). What this repository adds is the composition control, applied
  consistently, and the positive-control discipline above.
- **E-values are not comparable across library sizes.** A fixed `-e` cutoff on a
  10× library compares thresholds, not contents — which is why `61b` reports at
  two cutoffs.
