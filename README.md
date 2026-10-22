# EBNA1-MimicScan

A reproducibility-first screen for human proteins that the immune system could
mistake for EBNA1, and the negative results that constrain it.

**Read this first:** `docs/results/CANDIDATES.md` (the candidate list),
`METHOD.md` (the rule every ranking axis must pass), `docs/results/REPORT.md`
(the numbers, rendered from the tables).

Screen for human proteins that the immune system could mistake for EBNA1.

Two channels run here and they are deliberately separate. The anchored one
starts from the EBNA1 386-440 region, where antibodies and sera have validated
cross-reactive epitopes. The unbiased one takes all 641 EBNA1 residues against
20,431 reviewed human proteins and excludes no region. The deliverable is
`reports/CANDIDATES.md`: 55 composition-controlled local matches, ranked with the
statistical weakness printed on every card. Every claim is checked against what a
shuffle of the same residues — or of the same length — would have found.

## Install and run

```bash
git clone https://github.com/shiningleeeee/EBNA1-MimicScan.git
cd EBNA1-MimicScan

python -m venv env/venv
env/venv/bin/pip install -r requirements.txt   # pinned versions
make all                                       # gcc: three C scanners

bash scripts/01_fetch_data.sh   # the two sequence inputs, hash-verified
bash run.sh                     # 21 stages, ~375 s with warm caches
make test                       # pytest; skips the data-dependent tests if
                                # the inputs have not been fetched
```

The search backend (MMseqs2, static build) is **not** vendored. Stages from `s61`
on need it in `env/bin/mmseqs`; `make check-mmseqs` verifies it is there before a
run fails halfway through. The build this project was developed against is
`mmseqs` commit `fcf516aba34be4569b3eb5bda2b9db9e26a8f278` (static, avx2,
`https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz`, binary SHA256
`79e4c3722ae4ba75ea5d46a67cc6d56415b59cec6f910445c49a961ee6ab0753`).

```bash
bash run.sh --list              # stage list
bash run.sh --from s61_mmseqs   # resume from a stage (full stage id, not a prefix)
bash run.sh --only s73_candidate_cards
make archive                    # refresh docs/results/ from reports/
```

Everything the pipeline writes lands in `reports/`, which is not committed: the
published tables are frozen under `docs/results/` by `make archive`.

Outputs land in `reports/`. `reports/run_manifest.json` records the config
hash, the verified hash of every input, and the sha256 of every table the report
was rendered from — so a report can always be traced back to exactly the bytes
that produced it.

## Stages

| stage | what it does |
|---|---|
| `s00_provenance` | hash-checks every input against `data/MANIFEST.tsv`, records tool versions; refuses to continue on mismatch |
| `s10_seeds` | real seeds from published epitopes + two null seed families |
| `s20_motif_matches` | exact / near-exact scan of every seed over the proteome |
| `s30_annotate` | one bulk UniProt download supplies the CNS and surface axes |
| `s40_priority` | scores every seed, calibrates the real ones against their nulls |
| `s50_alignment` | alignment channel, kept for continuity — not evidence |
| `s60_composition` | EBNA1 composition profile and the Gly/Ala block |
| `s61_mmseqs` | **unbiased channel**: all 641 residues vs the proteome, no seed |
| `s62_pair_null` | per-alignment two-sided composition null, gaps and matrix unchanged |
| `s63_selection_null` | is the *number* of alignments itself above a composition shuffle? |
| `s64_rank` | eight declared axes, plus calibration of the whole list |
| `s65_structure` | AlphaFold proxy for accessibility of each matched segment |
| `s66_maturation` | does the matched span survive signal/propeptide removal? |
| `s67_struct_pair` | folded, exposed surface on both sides of the pair |
| `s68_layered` | independent annotations, stacked and never fused |
| `s69_iedb` | antibody-level evidence per accession, with references for scale |
| `s70_confirm` | HPA + PubMed confirmation layer for finalists (never scored) |
| `s71_candidates` | assembles the candidate table + the sensitivity audit |
| `s72_axis_validation` | **the rule, as a table**: every axis against the published positives |
| `s73_candidate_cards` | **the deliverable**: `CANDIDATES.md`, one card per shortlisted pair |
| `s80_report` | renders `REPORT.md` and `summary.json` from the tables |

Two channels, deliberately separate. `s10`–`s40` are anchored on the validated
386-440 peptide; `s61`–`s73` search the whole protein with no anchor and no
region excluded. Neither can validate the other — they fail in opposite
directions, and `REPORT.md` §6c states both.

`s72` is the methodological gate and the reason this repository is worth reading
even though it found nothing. An axis may **order** candidates only if it both
keeps the published positives and stays selective; the two thresholds are
declared in `config.yaml` before anything is measured
(`validation.min_positives_passing: 3`, `validation.max_pool_retained: 0.25`).
Of ten axes, **eight statistical ones fail** and one annotation axis is barred
from filtering at all:

| axis | kind | keeps positives | keeps pool | verdict |
|---|---|---|---|---|
| `maturation` | hard constraint | 4/4 | 0.967 | **pass** |
| `min_alignment_length` (≥20) | ranking | 1/4 | 0.645 | fail: loses positives |
| `min_identity` (≥40%) | ranking | 3/4 | 0.493 | fail: not selective |
| `margin_z_final` (≥3.0) | ranking | 3/4 | 0.920 | fail: not selective |
| `margin_z_iedb` (≥6.0) | ranking | 0/4 | 0.416 | fail: loses positives |
| `target_normal` | ranking | 0/4 | 0.807 | fail: loses positives |
| `query_folded_exposed` | ranking | 0/4 | 0.236 | fail: loses positives |
| `both_folded_exposed` | ranking | 0/4 | 0.182 | fail: loses positives |
| `composite_score` (≥6) | ranking | 0/4 | 0.468 | fail: loses positives |
| `comp_accessible` | annotation | 3/4 | 0.218 | annotation only |

The full method — scope and refusals, the four-rung null ladder with its
resolution arithmetic, quantity-versus-quality separation, the provenance
contract, and what would change the answer — is in [`METHOD.md`](METHOD.md).

## Result in one line

The prioritisation is not significant: compositions of the same residues reach
the same top score (p = 0.867 for `PPRRPP`, 0.267 for `RRPFF`, 0.343 for
`PGAIEQGPAD`), so the candidate ordering describes what the proteome contains,
not what this virus does. `REPORT.md` §3 is the statistic, §4 the ranking, §4b
what it can and cannot recover.

The structure channel returns the same verdict from the other direction
(`REPORT.md` §6b): EBNA1 is itself unfolded where it matters — mean pLDDT 41.4,
82% of residues in the very-low band, and all three validated epitopes in
disordered stretches — while 1 of 14 candidates and 18 of 75 null matches land
on a folded, exposed surface (p = 0.25). Reachability does not separate the list
either.

The unbiased channel (`REPORT.md` §6c) agrees for a third time: 363 alignments
across 363 human proteins, **none with E < 1**; the number of alignments beats a
composition shuffle (363 vs 19.4 ± 9.3, p at the 1/201 floor) but the best entry
does not (score 8 vs a null maximum of 8, p = 0.139). What does survive is
maturation: 12 of 363 matched spans sit in a signal peptide, propeptide or
transit peptide and are gone before the protein reaches the surface, and the
published mimics all survive that filter (3/3) while they fail the alignment one
(0/3).

The interface channel (`reports/REPORT_step9_paratope.md`) supplies the most
useful single number in the project. The one antibody with two solved peptide
structures — MS39p2w174 against EBNA1 386-405 (7K7R) and phospho-GlialCAM 370-389
(7UQC), same paratope at 0.63 Å — binds them at a clean register (+16) with 6/20
identity. That confirmed pair, measured by our composition null, gives
**margin z 0.865, p = 0.20**: indistinguishable from a shuffle of its own
residues. Sequence similarity and antibody cross-reactivity are not the same
quantity, and this is the only place where that is shown with experimental
structures rather than argued.

## Reproducibility

Verified on the reference container (Ubuntu 22.04, gcc 11.4, CPython 3.10.12):

- full pipeline, **21 stages, ~375 s** wall clock with the AlphaFold and IEDB
  caches warm (`bash run.sh`); the first run also downloads ~75 MB of models and
  takes roughly four minutes longer
- **bit-for-bit deterministic**: two independent full runs produce identical
  sha256 for every output file and for the seed table — 90 files compared
  pairwise, byte for byte. Two differences are expected and both are accounted
  for. `run_manifest.json` differs in its per-stage durations, never in the hash
  of an input or of any table whose contents are deterministic. And
  `69_iedb_summary.json` records `queried_utc`, the wall-clock time the IEDB API
  was hit — provenance, not a result. That single field propagates: the manifest
  also stores the hash **of that file**, so the recorded hash for
  `69_iedb_summary.json` differs between runs for the same reason the file does,
  and for no other.
- **74 tests** (`make test`), covering input hashes, release pinning, seed-table
  determinism, null p-value arithmetic, the SASA arithmetic, the m8 reader's
  ordering and width checks, the candidate table's invariants (no fused score,
  numeric ordering, `unknowns` populated, unqueried IEDB not recorded as zero),
  the axis verdicts in `s72` (each verdict must be recomputable from the declared
  thresholds), and a self-consistency check that re-derives every headline number
  in `REPORT.md` from the delivered tables
- every empirical p value is reported next to the resolution floor that produced
  it (`1/(1+b)`), so a number sitting on the floor cannot be read as sharper
  evidence

## Environment notes

- `pypi.org` may be unreachable; the venv was built with
  `uv pip install --python env/venv/bin/python` against a mirror.
- `apt-get install` is very slow here; the two C helpers are compiled with the
  system `gcc` (`make all`) instead of pulling a toolchain.
- The venv has no `pip`; install additions through `uv`.


## Why the nulls are the point

Any whole-proteome statistic on this protein is dominated by the fact that
EBNA1 residues 90-330 are a 241-residue block that is 99.6% Gly/Ala (a single
Arg at 329). Three independent null designs in this project (whole-proteome
composition shuffle, windowed shuffle against MMseqs2, and short-motif counting
under a window shuffle) all returned p > 0.11 for the sequence channel. What
survives is a filter, not a detection: the ranking is only as meaningful as the
ability of a composition-matched peptide to reach the same score, which is what
`s40` measures and `REPORT.md` §3 reports.

## Layout

```
config.yaml            every parameter; its hash is recorded in the run manifest
METHOD.md              the rule: which axes may rank, and how a null is chosen
data/MANIFEST.tsv      accession/url -> sha256 + UniProt release (verified at s00)
data/cache/            fetched inputs and HTTP caches (rebuildable, not committed)
src/phmimic/           the pipeline package (core, mining, scan, seeds, score,
                       structure, annotate, runner + stages/s00..s73)
src/c/                 mimic_scan.c, motif_mine.c, pair_null.c
src/analysis/          the exploratory scripts that produced docs/history/;
                       superseded by the stages, kept so the history is re-runnable
tests/                 provenance, determinism, SASA arithmetic, m8 reader,
                       candidate-table invariants, axis verdicts, self-consistency
docs/results/          frozen copies of the published tables (committed)
docs/history/          earlier reports, including the ones later revised
docs/figures/          figures
reports/               all generated; safe to delete and rebuild (not committed)
```

## Scope, and what this repository does not claim

- It does not claim EBNA1 does or does not cross-react with any of these
  proteins. Nothing here measures binding.
- It does not produce a "significant" candidate list. The ranking axes were
  tested against the four published mimics first (`docs/results/72_axis_validation.tsv`)
  and **0 of 8** survived; the list is a hypothesis set with a declared,
  non-significant order.
- It cannot see conformational epitopes (the GlialCAM/MBP mechanism), PTM state,
  or MHC presentation.
- The negative result here is not novel: Kanduc and colleagues showed in 2010
  that shared pentapeptides cover essentially every human protein, so mimicry
  screens need a composition control — that is what this project implements, not
  something it discovered.

## Licence and data

The code, documentation and derived tables are MIT (`LICENSE`); the upstream
data each table is derived from keeps its own licence and is not redistributed
here — see `DATA-LICENCES.md`.

Code, documentation and derived tables: MIT (see `LICENSE`). No upstream
sequence, structure or annotation data is redistributed: `data/MANIFEST.tsv`
records where each input came from and its SHA256, and `make data` refetches it
from the providers under their own licences (UniProt CC BY 4.0, AlphaFold DB
CC BY 4.0, RCSB PDB CC0, Human Protein Atlas CC BY 4.0).
