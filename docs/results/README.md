# Published tables

These files are frozen copies of the pipeline's output, taken from `reports/`
after the last full run. They are committed here so that every number quoted in
the documentation can be checked without downloading 200 MB of inputs and
re-running a six-minute pipeline.

To regenerate them from scratch:

```bash
make data        # fetch and hash-verify the two sequence inputs
bash run.sh      # 21 stages, ~375 s with warm caches
make archive     # copy the tables below back into this directory
```

| file | what it is |
|---|---|
| `CANDIDATES.md` | **the deliverable**: 55 composition-controlled local matches, one card per shortlisted pair |
| `REPORT.md` | the pipeline's own report, rendered from the tables (never hand-edited) |
| `summary.json` | every headline number, as data |
| `55_candidates_final.tsv` | the full candidate table: 55 rows x 40 columns, one column per piece of evidence |
| `73_candidate_cards.tsv` | the same pairs with the card fields (screen, entropy, null tail) |
| `72_axis_validation.tsv` | each ranking axis judged against the published positives; 0 of 8 statistical axes pass |
| `72_axis_positives.tsv` | the four published mimics as measured by this pipeline |
| `40_candidates_ranked.tsv` | the anchored channel's candidates |
| `45_local_similarity.tsv`, `46_pair_null.tsv` | the unbiased channel: every alignment MMseqs2 reports, with its composition null |
| `65_structure_real.tsv` | AlphaFold accessibility of each matched span |
| `panel_paratope.tsv` | 50 peptide-Fab complexes, used to test whether an interface concentrates on one residue |
| `60b_mmseqs_param_matrix.tsv`, `.json` | six-cell MMseqs2 audit separating composition correction from low-complexity masking (from `src/analysis/60_mmseqs_param_matrix.py`) |
| `61b_trembl_expansion.tsv`, `.json` | library-scope audit: reviewed Swiss-Prot vs + human TrEMBL, each with a re-fitted windowed null (from `src/analysis/61_trembl_expansion.py`) |
| `*_summary.json` | per-stage facts (counts, parameters, p values, caveats) |

Figures live in `docs/figures/`.
