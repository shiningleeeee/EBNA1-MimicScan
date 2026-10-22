# METHOD — a calibrated screen for molecular mimicry, specified end to end

This document fixes the method. It is written so that a second person can run
it on a different pathogen protein and get a result whose meaning is already
defined, and so that a reviewer can attack it at named places instead of
guessing at intent.

The method exists because a sequence-similarity screen for molecular mimicry is
easy to write and almost impossible to interpret. In this project the screen
found 363 human proteins matching EBNA1 at E < 10 000, and **every subsequent
test showed the ranking carried no information beyond local composition**. The
value of the method is not that it finds mimics. It is that it decides, in
advance and in code, what would count as evidence — and then reports what
happened when the answer was no.

---

## 0. Scope

**Answers.** Given a pathogen protein and a host proteome, which host proteins
share a local stretch with the pathogen protein, how surprising is each match
under a null model that preserves local composition, and which of the declared
screening axes have demonstrated that they can discriminate anything at all.

**Refuses to answer.** Whether a listed match is a real cross-reaction. That
requires affinity, processing or antibody-binding data. The method produces a
ranked hypothesis set with its own refutations attached; it does not produce
"candidates" in the sense of validated binders.

**Explicitly out of scope.** Conformational epitopes, post-translational
modifications other than those carried as non-standard residues in a deposited
structure, MHC presentation, and any inference that statistical
non-significance means the biology is absent.

---

## 1. The one rule that makes the rest meaningful

> **A statistical axis may order candidates only after it has been shown to
> keep the experimentally established cases and to remain selective.**

Everything else in this document is machinery for applying that rule honestly.
It is enforced mechanically in stage `s72_axis_validation`, whose verdict table
is generated, not written by hand.

Each axis is declared one of three kinds:

| kind | meaning | may it rank? |
|---|---|---|
| `hard_constraint` | a mechanistic necessity — a match inside a propeptide that is cleaved away cannot be an antibody target, whatever its score | yes, as a filter |
| `ranking_axis` | a statistical selection; must pass the positive-control test | only if `pass` |
| `annotation` | informative context, never selective | never |

Two numbers are declared in `config.yaml` **before** the stage measures
anything, so a verdict cannot be tuned into existence:

```yaml
validation:
  min_positives_passing: 3     # of the published pairs
  max_pool_retained: 0.25      # an axis keeping more than a quarter is not selective
```

### 1.1 The result on this screen (generated: `reports/72_axis_validation.tsv`)

| axis | kind | positives kept | pool kept | verdict |
|---|---|---|---|---|
| maturation (span survives proteolytic processing) | hard_constraint | 4/4 | 0.967 | **pass** |
| aligned length ≥ 20 aa | ranking_axis | 1/4 | 0.645 | fail: loses positives |
| identity ≥ 40 % | ranking_axis | 3/4 | 0.493 | fail: not selective |
| margin z ≥ 3.0 | ranking_axis | 3/4 | 0.920 | fail: not selective |
| margin z ≥ 6.0 | ranking_axis | 0/4 | 0.416 | fail: loses positives |
| target window composition normal | ranking_axis | 0/4 | 0.807 | fail: loses positives |
| query span folded and exposed | ranking_axis | 0/4 | 0.236 | fail: loses positives |
| both spans folded and exposed | ranking_axis | 0/4 | 0.182 | fail: loses positives |
| compartment accessible | annotation | 3/4 | 0.218 | annotation only |
| eight-axis composite ≥ 6 | ranking_axis | 0/4 | 0.468 | fail: loses positives |

**Eight statistical axes, none survives.** The four published pairs score
5, 4, 4, 4 out of 8 on the composite that was originally used to rank the
candidates; the composite threshold that produced the headline list (6) keeps
none of them.

The published set deliberately keeps the two pairs a linear screen cannot see —
GlialCAM, whose epitope is defined by a crystal structure rather than by
identity, and MBP, which is conformational — so that no axis can be credited
for recovering only the easy cases. Per-pair evidence:
`reports/72_axis_positives.tsv`.

### 1.2 What follows from a failed axis

A failed axis may still be **displayed**. It may not filter, order, or enter a
composite score. This is why the deliverable is a table of aligned candidates
with columns of unweighted evidence and an `unknowns` column, not a single
"priority score": any weighting of axes at least one of which is unvalidated
launders noise into apparent precision.

---

## 2. The null ladder

Four nested nulls, each preserving more of the query than the last. The number
of matches a screen reports falls as the null gets stricter, and that fall is
the measurement — not an inconvenience.

| level | what it preserves | what it destroys | question it answers |
|---|---|---|---|
| 1. global shuffle | amino-acid composition of the whole query | all positional information | is the match count more than a composition-matched protein would give? |
| 2. window shuffle (10 / 20 / 40 aa) | composition *within each window* | order inside the window | does local arrangement matter? |
| 3. per-pair, both sides | composition of the query span *and* the target span | everything else | is this one alignment better than two composition-matched stretches? |
| 4. exact enumeration | the full permutation multiset | nothing | the same question, without Monte-Carlo error, when the space is small enough |

Measured on this screen: level 1 gives 1.23× (seeded channel) and 144.8×
(global, MMseqs2); level 2 gives p = 0.176 / 0.170 / 0.141 for 10 / 20 / 40 aa
windows; level 3 with B = 10⁶ bounds each published pair at
margin z = 0.87–5.01; level 4 was used for the six-residue seeds where the
permutation space is 15 and 30 arrangements.

**Reading rule.** A large ratio at level 1 with a null result at level 2 means
the signal was local composition. That is the finding here, and it is why the
global-shuffle ratio is reported next to the windowed p and never alone.

### 2.1 Resolution accounting

Every empirical p is reported with its floor `1/(1+B)`. With B = 20 the
smallest attainable p is 0.0476 — three "p = 0.048" results in an earlier round
were three copies of the floor, not three findings. A p value quoted below
`1/(1+B)` is a defect. `B ≥ 1/α − 1` for the smallest α claimed.

### 2.2 What a null does not do

Shuffling preserves composition and destroys *every* order-dependent feature,
including any real short linear motif. A window null that kills a three-residue
motif has not shown the motif to be unreal — it has shown that the test cannot
see it. Every null in the table above is therefore reported with the biological
question it can and cannot address, and no null result is written as
"the residues carry no signal".

---

## 3. Separation of quantity from quality

A screen that reports more alignments than a shuffled query is not thereby
showing better alignments. The two are measured separately:

- **selection**: matches reported, real vs shuffled. Here 363 vs 19.4 ± 9.3,
  p = 0.005 (floor, 1/201).
- **quality**: best match, real vs shuffled. Here composite 8 vs 6.61 ± 0.83
  over 200 shuffles, of which **27 also reach 8** → p = 0.139.

Comparing counts across a screen where the real query reports 18× more
alignments than a null inflates any apparent effect by roughly that factor.
The honest comparison for a threshold is the **fraction** of reported
alignments above it: 46.8 % vs 20.2 % ± 12.4 %, p = 0.025. An earlier version
of this document reported the count comparison as a 44× enrichment; the
fraction comparison gives 2.3×.

---

## 4. Composition diagnostics come first, not last

Before any alignment is scored, the query is described:

- per-window entropy and maximum single-residue fraction, sliding window 20 aa;
- the longest low-complexity block, its exact composition and length;
- a census of the exact motifs of interest in the host proteome.

For EBNA1 that block is 241 residues at 99.59 % Gly/Ala, three distinct residue
types, the only non-Gly/Ala residue being R329. **Zero** of the 363 reported
alignments overlap it, and the reason for that is *not* established: MMseqs2's
default low-complexity masking and `--comp-bias-corr 1` both act on it, and
their separate contributions have not been audited. The report says exactly
that, and no more.

The census that decides how a motif may be read:

| motif | proteins (of 20 431) | occurrences |
|---|---|---|
| `PPRRPP` | 6 | 6 |
| `RRPFF` | 3 | 3 |
| `PGAIEQGPAD` | 0 | 0 |
| Hecker 17-mer `GSPPRRPPPGRRPFFHP` | 0 | 0 |
| `P.R.P` (3-position abstraction) | 2 841 (13.9 %) | 3 719 |

A motif's rarity says nothing about its specificity: `P.R.P` is common because
it is short and class-poor, and the exact 17-mer is absent because it is long.
Neither number is evidence about binding.

---

## 5. Provenance and determinism contract

**Provenance.** Every input is pinned by release and SHA-256 in
`data/MANIFEST.tsv`; the pipeline refuses to run on a hash or release mismatch
(`ProvenanceError`). In this run: UniProt 2026_03 (2026-09-02), EBNA1 `P03211`
641 aa `56d0ecac…`, human reviewed 20 431 sequences `bf1bc7e1…`.

**Determinism.** Two full runs must agree byte for byte on every output except
two fields, both of them deliberately not reproducible:

- `run_manifest.json` → `duration_s` (wall clock, recorded because it is part
  of the run, not because it should be stable);
- `69_iedb_summary.json` → `queried_utc` (when the external database was read;
  a re-run legitimately reads at a different time).

Everything else — including every hash in the manifest, which covers the
config, every input, and every output — must match. Multi-threaded search
output is sorted by `(query, evalue, target, qstart, tstart)` before writing so
that thread scheduling cannot leak into a file hash.

**Declared parameters.** Every threshold that can affect a ranking lives in
`config.yaml` with a comment giving the reason. `seed_score_min: 14` and
`gap_ext: 1` are carried and marked **NO-OP** and **DEAD** respectively,
because they are inert in the reference implementation and pretending otherwise
would misrepresent what was run.

---

## 6. The stage graph

One command: `bash run.sh`. Twenty stages, each writing to `reports/` and
recording its duration and output hashes into `reports/run_manifest.json`.

| stage | what it produces |
|---|---|
| `s00_provenance` | hash/release verification, tool versions |
| `s10_seeds` | declared seeds and their null arrangements |
| `s20_motif_matches` | exact / mismatch-bounded motif matches |
| `s30_annotate` | one bulk UniProt annotation fetch (CNS, surface) |
| `s40_priority` | motif-level scores with permutation calibration |
| `s50_alignment` | reference SW implementation, seeded + exhaustive |
| `s60_composition` | low-complexity profile of the query |
| `s61_mmseqs` | unbiased full-length search, all matches kept |
| `s62_pair_null` | per-pair composition null, B = 10⁴ |
| `s63_selection_null` | quantity vs quality, outer 200 × inner 2 000 |
| `s64_rank` | eight axes, composite, whole-list calibration |
| `s65_structure` | AlphaFold accessibility of the query and of null matches |
| `s66_maturation` | does the matched span survive proteolytic processing |
| `s67_struct_pair` | both-sides accessibility for the filtered pool |
| `s68_layered` | evidence columns stacked, never fused |
| `s69_iedb` | external B-cell records, reported as coverage only |
| `s70_confirm` | expression and literature context, not scored |
| `s71_candidates` | the candidate table + a read-only benchmark audit |
| `s72_axis_validation` | the verdict table of §1.1 |
| `s73_candidate_cards` | `CANDIDATES.md` -- the candidate list itself |
| `s80_report` | `REPORT.md` and `summary.json`, interpolated from tables |

`tests/` holds 61 assertions; the ones that matter are the self-consistency
tests, which recompute every headline number from the delivered tables rather
than trusting the report that quotes it. A stale `config.sha256` in
`summary.json` is a test failure, by design.

---

## 7. Decision rules

1. A host protein enters the candidate table if it shares an aligned stretch
   with the query, that stretch survives maturation, and the alignment is
   longer than the declared minimum. These are mechanics, not evidence.
2. No candidate is ranked by any axis whose verdict in §1.1 is not `pass`.
   With the current table, **no statistical ranking is permitted at all**, so
   the table is ordered by margin z with the admission that margin z has not
   passed.
3. Benchmark recall is a sensitivity audit and appears in its own read-only
   table (`reports/71_benchmark_audit.tsv`). It never reorders candidates, and
   a low recall is never reported as a pipeline defect: on this screen the
   unbiased channel recovers 0 of 3 linear published mimics while the
   epitope-seeded channel recovers 3 of 3 — the difference measures the seeds,
   not the quality.
4. Every candidate carries `unknowns`, naming what was not measured. An empty
   `unknowns` is a defect.
5. A negative result is reported as "no candidate survived the declared
   calibration", never as "there is no mimicry".

---

## 8. Worked result on this screen

- Unbiased full-length search: 363 alignments, 363 distinct human proteins,
  best E = 3.088, **none below E = 1**.
- Global shuffle: 144.8× more alignments than a composition-matched query.
- Local-composition nulls (10 / 20 / 40 aa): p = 0.176 / 0.170 / 0.141. The
  advantage disappears when local composition is preserved.
- Per-pair null (B = 10⁶) on the only crystal-confirmed pair
  (EBNA1 386–405 ↔ GlialCAM 370–389): margin z = 0.865, p = 0.20.
  **A pair proven to cross-react is not distinguishable from a
  composition-matched permutation.**
- Whole-list calibration: best composite 8 vs 6.61 ± 0.83, 27/200 shuffles also
  reach 8, p = 0.139.
- Axis validation: 0 of 8 statistical axes pass; 1 mechanistic constraint passes.
- External correlates of a mechanistic story we proposed were tested and
  **refuted**: across 50 PDB peptide–antibody complexes the dominant residue
  carries a median 24.4 % of contacts, and the two calibration crystals sit at
  21.4 % and 19.9 % — below median, not above. The "the antibody only reads one
  arginine" reading came from a `{**light, **heavy}` dict merge that silently
  discarded the heavy chain; corrected, the claim does not hold.

---

## 9. Reuse

To run this on another pathogen protein:

1. Replace `data.query` in `config.yaml` and record the accession, release and
   SHA-256 in `data/MANIFEST.tsv`.
2. Replace `published:` with the experimentally established pairs for that
   protein, including **at least one case a linear screen is expected to miss**.
   Without that, §1 tests nothing.
3. Re-derive `composition.gar_block` from stage `s60` — do not copy EBNA1's.
4. Run `bash run.sh`. Read `reports/72_axis_validation.tsv` **before** reading
   any candidate list. If no ranking axis passes, the candidate list is a
   hypothesis set and must be reported as one.
5. Report the null ladder, not the global ratio. Report the selection and
   quality statistics separately. Report the resolution floor of every p.

The method's transferable output is the verdict table, not the candidates.

---

## 10. What would change the answer

- A valid ranking axis: one that keeps ≥ 3 of 4 published pairs and keeps
  ≤ 25 % of the pool. None exists yet; finding one is the open problem.
- Independent evidence: antibody-binding, affinity, or processing data for a
  specific candidate. Statistics cannot substitute, and the sequence channel is
  exhausted for this query.
- A different evidence type: interface thermodynamics, or a paratope with a
  measured affinity for both peptides. Note that the cheapest version of this —
  co-folding predicted structures — is explicitly *not* used, because without
  calibration data it produces a plausible complex for anything.
