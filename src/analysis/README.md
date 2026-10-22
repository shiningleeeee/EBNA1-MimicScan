# Exploratory analysis scripts

These are the scripts that produced the material in `docs/history/`. They are
kept because that history is part of the record — including the parts that were
later revised — not because they are the recommended entry point.

**The pipeline supersedes most of them.** Everything needed to reproduce the
candidate table, the null calibrations, the maturation filter, the structure
channel and the IEDB layer now runs as numbered stages via `bash run.sh`:

| script | superseded by |
|---|---|
| `45_*` (unbiased MMseqs2 search) | `s61_mmseqs` |
| `47_pair_null_batch.py` | `s63_selection_null` |
| `48_rank_final.py`, `49_annotate_top.py` | `s64_rank` |
| `51_processing_state.py` | `s66_maturation` |
| `52_biology_first.py` | `s67_struct_pair` |
| `53_layered_table.py` | `s68_layered` |
| `54_iedb_evidence.py` | `s69_iedb` |
| `60_mmseqs_param_matrix.py` | **not superseded** — this is the parameter audit that `s61_mmseqs` cites in its `note`; its output is frozen in `docs/results/60b_mmseqs_param_matrix.tsv`. It answers a question the pipeline deliberately leaves open (which mechanism suppresses the Gly-Ala block) rather than reproducing a stage |
| `mine_local_similarity.py` | `s61_mmseqs` |
| `panel_paratope.py`, `probe_paratope.py`, `probe_read_motif.py` | not in the pipeline: the 50-complex paratope panel is a one-off check, archived in `docs/results/panel_paratope.tsv` |
| `56_almulla_null.py`, `57_almulla_specificity.py` | not in the pipeline: their input is a published paper, not this project's candidate set |
| `61_trembl_expansion.py` | **not superseded** — the library-scope audit: reruns the unbiased search against reviewed Swiss-Prot plus human TrEMBL (205,766 canonical entries) and re-fits the windowed null in the enlarged library. Output frozen in `docs/results/61b_trembl_expansion.tsv` |
| `legacy/01..14` | the first five rounds; kept for provenance only |
| `s45_window_scan_never_run.py` | **written, never registered, never run.** A windowed-shuffle scan of all 641 EBNA1 residues that was designed and then abandoned; the pipeline's unbiased channel is `s61`-`s63`. It lives here rather than under `stages/` because `make check-stages` fails on any stage module the runner does not know about, and registering code that was never run would misrepresent it. |

One of these closed a caveat rather than a claim: `60_mmseqs_param_matrix.py`
turned `s61`'s "the reason is not established" into a measured attribution
(composition correction, not masking), after which the stage's `note` and
`METHOD.md` were rewritten to say so.

Two more were used to **overturn** a claim in `docs/history/`, which is the
reason they are worth keeping: `08_where.py` produced the "all alignments lie in
350-641" statement that a later recomputation contradicted, and
`panel_paratope.py` refuted the "the interface concentrates on one arginine"
reading that `probe_paratope.py` had produced from a mis-keyed dictionary.

Run them from the repository root with `PYTHONPATH=src env/venv/bin/python ...`.
