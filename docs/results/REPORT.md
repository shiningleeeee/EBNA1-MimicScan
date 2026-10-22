# EBNA1 molecular mimicry -- prioritised candidate screen

`ph-mimic` v1.0.0 -- config sha256 `8b1de6ed868bf37d`

## 0. What this is, and what it is not

A screen for human proteins that could be mistaken for EBNA1 by the immune system. The only experimentally anchored object is a short peptide from the EBNA1 386-440 region, where antibodies and sera have validated cross-reactive epitopes; every candidate here is a match to that peptide, ranked by biology that no sequence statistic can supply.

It is **not** a list of confirmed cross-reactive proteins. Nothing in this report is a measurement; it is a prioritisation. The score is a filter, and the null below states what a filter of this kind finds when the peptide means nothing.

## 1. Provenance

| input | file | UniProt release | sha256 |
|---|---|---|---|
| `data.proteome` | `data/cache/human_swissprot.fasta` | 2026_03 | `bf1bc7e188e55199...` |
| `data.query` | `data/cache/P03211.fasta` | 2026_03 | `56d0ecacb15063e7...` |

| tool | version |
|---|---|
| gcc | `gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0` |
| mimic_scan | `usage: mimic_scan --query Q.fa --db DB.fa --mode seeded|exhaustive --out O.tsv` |
| mmseqs | `fcf516aba34be4569b3eb5bda2b9db9e26a8f278` |
| python | `3.10.12` |
| venv_python | `env/venv/bin/python` |

## 2. Seeds and the nulls they are judged against

- **real seeds (3)**: `PPRRPP` (EBNA1 394-399, d<=0), `RRPFF` (EBNA1 402-406, d<=0), `PGAIEQGPAD` (EBNA1 431-440, d<=3)
- **N_A (243)**: every distinct arrangement of each seed's own residue multiset, enumerated exhaustively. Composition is held exactly constant, so anything separating a real seed from N_A is due to order alone.
- **N_B (200)**: random 6-aa windows drawn from the human proteome. Composition is free; this answers "what would an arbitrary peptide have found".

Real and null seeds are scanned by the same code and scored by the same rule.

## 3. Primary result: is the seed better than its own null?

Score = exactness (2 for a full-length match, 1 for >= 0.70 identity) + CNS expression (1) + antibody-reachable localisation (1). Maximum 4.

| seed | motif | best score | best candidate | proteins matched | N_A n | N_A max | p (unconditional) | p (conditional) | N_B p |
|---|---|---|---|---|---|---|---|---|---|
| `C1` | `PPRRPP` | 4 | CNTN2 | 6 | 14 | 4 | 0.867 | 0.867 | 0.149 |
| `C2` | `RRPFF` | 4 | MAGEE1 | 3 | 29 | 4 | 0.267 | 0.276 | n/a |
| `C3` | `PGAIEQGPAD` | 3 | ANO2 | 5 | 200 | 3 | 0.343 | 0.350 | n/a |

- *p (unconditional)* treats the whole null seed space as the reference, including null peptides that occur nowhere in the human proteome.
- *p (conditional)* restricts the reference to null peptides that match at least one human protein, which is the fairer comparison when the question is "given that this peptide exists in us, does its order matter".
- *N_B p* is blank where the real seed length differs from the N_B length; lengths are not comparable and no p is reported rather than a misleading one.

## 4. Candidates

14 protein-seed pairs pass the sequence filter; 14 distinct proteins. Top 14 by score:

| # | protein | gene | seed | matched | ident | CNS | surface | score | published |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `Q02246` | CNTN2 | `PPRRPP` | `PPRRPP` @911 | 1.000 | Y | Y | 4 | - |
| 2 | `Q9HCI5` | MAGEE1 | `RRPFF` | `RRPFF` @710 | 1.000 | Y | Y | 4 | - |
| 3 | `Q9NQ90` | ANO2 | `PGAIEQGPAD` | `PGDIELGPLD` @140 | 0.700 | Y | Y | 3 | yes |
| 4 | `Q02080` | MEF2B | `PPRRPP` | `PPRRPP` @272 | 1.000 | Y | - | 3 | yes |
| 5 | `P02511` | CRYAB | `RRPFF` | `RRPFF` @11 | 1.000 | - | Y | 3 | yes |
| 6 | `Q96BY6` | DOCK10 | `PGAIEQGPAD` | `PSAFFQGPAD` @1550 | 0.700 | Y | - | 2 | - |
| 7 | `O43493` | TGOLN2 | `PGAIEQGPAD` | `SGAEEQGPID` @224 | 0.700 | - | Y | 2 | - |
| 8 | `A1X283` | SH3PXD2B | `PPRRPP` | `PPRRPP` @758 | 1.000 | - | - | 2 | - |
| 9 | `P48634` | PRRC2A | `PPRRPP` | `PPRRPP` @1706 | 1.000 | - | - | 2 | - |
| 10 | `Q5T089` | MORN1 | `PPRRPP` | `PPRRPP` @18 | 1.000 | - | - | 2 | - |
| 11 | `Q8N201` | INTS1 | `PPRRPP` | `PPRRPP` @841 | 1.000 | - | - | 2 | - |
| 12 | `Q9H4Z2` | ZNF335 | `RRPFF` | `RRPFF` @720 | 1.000 | - | - | 2 | - |
| 13 | `Q07864` | POLE | `PGAIEQGPAD` | `PGAIRDGPAT` @1304 | 0.700 | - | - | 1 | - |
| 14 | `Q96L42` | KCNH8 | `PGAIEQGPAD` | `TGAHEQNPAD` @973 | 0.700 | - | - | 1 | - |

`CNS` and `surface` come from a single bulk UniProt download pinned in the run manifest, not from per-protein lookups, so the table and the null distribution share one annotation snapshot.

## 4b. Sensitivity against the published claims

This screen recovers 3 of 5 published EBNA1 mimicry claims recorded in `config.yaml`:

- recovered: `P02511` CRYAB -- alpha B-crystallin, RRPFF core, linear (Thomas 2023 Sci Adv 10.1126/sciadv.adg3032)
- recovered: `Q9NQ90` ANO2 -- anoctamin-2, PGAIEQGPAD/PGDIELGPLD 7 of 10, linear (Tengvall 2019; Lum 2026 PNAS Table 1)
- recovered: `Q02080` MEF2B -- MEF2B, shared PPRRPP motif, cross-reactivity of a clonal antibody (Lum 2026 PNAS 10.1073/pnas.2519445123 (PPRRPP section, not Table 1))
- **not recoverable**: `Q14CZ8` HEPACAM -- only 6 of 20 residues identical to EBNA1 386-405 (SSPPRA vs PPRRPP); the epitope is defined structurally, not by linear identity
- **not recoverable**: `P02686` MBP -- no linear homology claimed by the source; a peptide screen cannot see a conformational epitope

Recall is reported as a sensitivity check on the screen, never as a measure of its reliability: a screen can recover a known mimic by composition alone, and failing to recover one can mean the epitope is conformational rather than that the code is wrong. The null in section 3 is the statistic; this list is an audit.

## 5. Confirmation layer (not scored)

| protein | gene | HPA brain nTPM | HPA brain cluster | HPA single-cell cluster | PubMed (MS/EBV/auto) |
|---|---|---|---|---|---|
| `Q02246` | CNTN2 | 157.3 | Cluster 11: White matter - Signal transduction | Cluster 80: Oligodendrocytes - Myelination | 4 |
| `Q9HCI5` | MAGEE1 | 25.0 | Cluster 1: Neurons - Mixed function | Cluster 98: Neurons - Neurotransmission & signaling | 0 |
| `Q02080` | MEF2B | - | Cluster 11: White matter - Signal transduction | Cluster 67: Plasma cells - Plasma cell differentiation | 5 |
| `P02511` | CRYAB | - | Cluster 28: Oligodendrocytes - Mixed function | Cluster 61: Schwann cells & Melanocytes - Mixed function | 27 |
| `Q9NQ90` | ANO2 | - | Cluster 47: Endothelial cells - Vasculature | Cluster 12: Endothelial cells - Angiogenesis & vascular home | 12 |
| `A1X283` | SH3PXD2B | - | Cluster 6: Non-specific - Immune response | Cluster 103: Non-specific - Basic cellular processes | 0 |
| `P48634` | PRRC2A | - | Cluster 51: Neurons - Mixed function | Cluster 102: Non-specific - Mixed function | 2 |
| `Q5T089` | MORN1 | 24.4 | Cluster 2: Neurons - Nucleosome | Cluster 107: Ciliated cells - Cilia assembly & function | 0 |
| `Q8N201` | INTS1 | - | Cluster 32: Cerebral cortex - Mixed function | Cluster 65: Non-specific - Cell growth & division | 1 |
| `Q9H4Z2` | ZNF335 | - | Cluster 15: Non-specific - Transcription | Cluster 65: Non-specific - Cell growth & division | 0 |
| `Q96BY6` | DOCK10 | 34.2 | Cluster 11: White matter - Signal transduction | Cluster 80: Oligodendrocytes - Myelination | 2 |
| `O43493` | TGOLN2 | - | Cluster 28: Oligodendrocytes - Mixed function | Cluster 72: Macrophages - Phagocytosis & lysosomal degradati | 0 |
| `Q07864` | POLE | - | Cluster 15: Non-specific - Transcription | Cluster 37: Neutrophil progenitors - Granulopoiesis | 121 |
| `Q96L42` | KCNH8 | 10.2 | Cluster 11: White matter - Signal transduction | Cluster 59: Pituitary endocrine cells - Neuroendocrine signa | 0 |

The two HPA cluster columns answer different questions and can disagree: the brain cluster is a bulk-tissue cluster, the single-cell cluster names the cell type that actually expresses the protein. Where they disagree, the single-cell cluster is the one that should drive a CNS claim.

Literature counts are **not** part of the score. A protein with many papers has been looked at more, which is not evidence that it is a mimic.

## 6. Why composition explains the sequence channel

EBNA1 (641 aa) residues 90-330 are 99.6% Gly/Ala (241 residues, 3 distinct residue types; the only residues that are not Gly or Ala are R329). Outside that block the highest 20-aa Gly/Ala fraction is 0.95.

A block of this length over a two-letter alphabet will match collagen Gly-X-Y repeats, keratin and other low-complexity human proteins under almost any composition-sensitive statistic. That is why this pipeline never ranks on an alignment score.

- alignment channel, exhaustive: 80,255 hits; 16617 overlap the Gly/Ala block.
- alignment channel, seeded: 59,108 hits; 194 overlap the Gly/Ala block.

The alignment channel is retained for continuity with earlier work and is reported, but it does not contribute to any candidate call.

## 6b. Structure channel: can an antibody reach the matched segment?

The query protein itself is not folded. AlphaFold's mean pLDDT for EBNA1 is 41.4 and 82% of its residues are in the very-low-confidence band. The three experimentally validated linear epitopes sit in that band:

| span | mean pLDDT | mean rSASA | call |
|---|---|---|---|
| C1 394 399 | 30.4 | 0.90 | disordered_or_loop_exposed |
| C2 402 406 | 27.0 | 0.77 | disordered_or_loop_exposed |
| C3 431 440 | 26.0 | 0.81 | disordered_or_loop_exposed |
| block | 33.3 | 0.94 | disordered_or_loop_exposed |

Every validated mimicry epitope lies in a disordered, highly exposed stretch, and the only confidently folded parts of EBNA1 (roughly 476-585) are the DNA-binding / dimerisation core, which carries no known epitope. A structure-based mimicry search over this protein has nothing folded to compare in the epitope region.

The same measurement, applied to each candidate's matched segment and to null matches from both families with the same code path:

| call | candidates |
|---|---|
| disordered_or_loop_exposed | 11 |
| partly_buried | 2 |
| surface_of_folded_domain | 1 |

Folded and exposed (mean pLDDT >= 70.0, mean rSASA >= 0.25): **1 of 14 candidates** against **19 of 76** null matches, giving p = 0.2597. Accessibility is therefore not a selective property here: an arbitrary matching hexamer in a human protein lands on a folded surface about a quarter of the time.

This channel changes no score and adds no candidate. It removes none either -- what it removes is the idea that reachability distinguishes the list. pLDDT is model confidence and rSASA is computed on a predicted monomer; both are proxies for accessibility, not measurements. This stage is deterministic once the model cache is complete; a match whose model could not be retrieved is listed under 'excluded' rather than counted as unexposed.

## 6c. The unbiased channel: full EBNA1 against the whole proteome

The channel above is anchored on the validated peptide. This one asks the same question without an anchor: all 641 residues searched against all 20431 reviewed human proteins, no seed, no region excluded.

It reports 363 alignments over 363 proteins. The best E-value is 3.088 (PDE11A), and **0 alignments have E < 1**. Not one human protein is a credible homologue of EBNA1 on this evidence; that is a statement about local similarity, not about mimicry.

| EBNA1 region | alignments |
|---|---|
| folded core (476-585) | 276 |
| C-terminal | 68 |
| SPR-rich (370-445) | 12 |
| N-terminal | 7 |

Alignments overlapping the Gly/Ala block: **0**. The reason is **not** established: mmseqs2 default low-complexity handling and --comp-bias-corr are both active; their separate contributions to the absence of 90-330 hits have not been audited (parameter matrix outstanding).

Is the *number* of alignments itself surprising? Shuffling the query 200 times while preserving its composition gives 19.4 +- 9.3018 alignments; the real query reports 363. That is p = 0.0050, which is the floor of the test (1/201) and not a measurement below it. Read as: the real sequence does line up with the proteome more often than its own anagram does, by a factor of about 18.7.

For each individual alignment, a two-sided composition null (b = 10000 rearrangements of both partners, gaps and matrix identical to the reference implementation) gives a margin over the null ranging 0.602 to 11.16 (median 5.669). That statistic is length-dependent (r = -0.40), so it orders short alignments first by construction.

Calibrating the whole list against 200 composition shuffles: the best entry scores 8 against a null maximum of 6.605 +- 0.8261 (null max 8.0), p = 0.139. **The best candidate is not special.** The *proportion* of alignments clearing the mid threshold is 46.8% against 20.2% +- 12.4%, p = 0.025: a real but modest enrichment, roughly 2.3x. Counting entries rather than proportions would have inflated this ~18x, because the real query simply reports more alignments.

Before any of that, one filter is applied that has nothing to do with significance: does the matched span survive maturation? 350 of 363 spans are confirmed present on the mature chain, 1 cannot be decided from UniProt, and 12 lie in a signal peptide, propeptide or transit peptide and are removed before the protein reaches the surface. This is the only axis here that a known positive also passes (3/3 of the published linear mimics survive it).

The candidate table has **55 rows** and **40 columns**, one column per stage. No composite score: axes are stacked, never multiplied. `unknowns` and `manual_judgement` are part of the deliverable.

Antibody-level evidence, queried per accession: 23 proteins asked, 22 returned no IEDB record at all, 0 had a record overlapping the matched span. For scale, the same query returns CRYAB 57, EBNA1 507, MEF2B 0, HEPACAM 156, ANO2 70. Coverage tracks who has been studied; it is recorded, never used as a filter, and an empty result is not evidence of absence.

Finally, our own sensitivity, computed after the table and read by nothing: this channel recovers **0 of 3** of the reachable published mimics, because each is a short low-complexity epitope that an aligner scoring length and identity cannot rank. The seeded channel above recovers all three. Two channels, opposite failure modes, stated together rather than picking the flattering one.

## 6d. The rule, as a table: no statistical axis passes its own controls

Every axis in this pipeline is scored against the published positives before it is allowed to order anything. Two thresholds were declared in `config.yaml` before the measurement: keep at least **3** of the published pairs, and keep at most **0.25** of the candidate pool. An axis that misses either one cannot discriminate, and using it as a ranking criterion would be a choice dressed up as a result.

| axis | kind | keeps positives | keeps pool | verdict |
|---|---|---|---|---|
| `maturation` | hard_constraint | 4 | 0.967 | pass |
| `min_alignment_length` | ranking_axis | 1 | 0.645 | fail: loses positives |
| `min_identity` | ranking_axis | 3 | 0.493 | fail: not selective |
| `margin_z_final` | ranking_axis | 3 | 0.920 | fail: not selective |
| `margin_z_iedb` | ranking_axis | 0 | 0.416 | fail: loses positives |
| `target_normal` | ranking_axis | 0 | 0.807 | fail: loses positives |
| `query_folded_exposed` | ranking_axis | 0 | 0.236 | fail: loses positives |
| `both_folded_exposed` | ranking_axis | 0 | 0.182 | fail: loses positives |
| `comp_accessible` | annotation | 3 | 0.218 | annotation_only |
| `composite_score` | ranking_axis | 0 | 0.468 | fail: loses positives |

**0 of 8 ranking axes pass.** The only axis that does is a hard mechanistic constraint (does the matched span survive signal-peptide or propeptide removal), which is not a statistical choice at all — it is a fact about where the protein is cut. The practical consequence is stated in `METHOD.md` §7: with no passing ranking axis, no column in the candidate table is allowed to be presented as a statistically justified ordering. The table is ordered by margin z with that failure printed next to it.

The positives this is measured on are deliberately not all easy. Two of the five are pairs a linear screen cannot see at all — GlialCAM, whose epitope is defined by a crystal structure and shares only 6 of 20 residues with EBNA1, and MBP, which is conformational. An axis set that was only tested on the three short linear epitopes would have been tested on exactly the cases it was built for.

## 7. What would change the conclusion

- A real seed reaching a score that no member of N_A reaches (p_conditional below the multiple-testing threshold) would make the ordering informative.
- A candidate whose protein is CNS-expressed, antibody-reachable and carries the motif in a disordered or surface-exposed segment would be worth wet-lab follow-up regardless of the p-value.
- The opposite reading is also live: if N_A peptides reach the same scores, this screen is a description of the proteome's composition and the honest output is the negative result, not a candidate list.

## 8. Reproduce

```bash
bash run.sh                 # full pipeline
env/venv/bin/python -m phmimic.runner --list
env/venv/bin/python -m pytest tests -q
```

`reports/run_manifest.json` records the config hash, every verified input hash, and the sha256 of every table this report was rendered from.
