# Unbiased local-similarity mining: EBNA1 → human proteome

*Generated from `reports/45_local_similarity.tsv`; reproducible with
`env/venv/bin/python src/mine_local_similarity.py`.*

## What was done

The query is the entire 641-aa EBNA1 (UniProt `P03211`, release 2026_03). It was
searched against all 20,431 human reviewed Swiss-Prot entries with MMseqs2
(`-s 7.5 --max-seqs 20000 -e 10000 --comp-bias-corr 1`). No part of EBNA1 was
masked, no epitope was used as a seed, and nothing about the three published
mimics influenced which proteins appear. What decides inclusion is only how well
a stretch of EBNA1 aligns to a human protein.

## Result

363 human proteins are reported, each by exactly one alignment. **Every one has
E > 1** — the best is E = 3.09. By the ordinary standard for a database search
of this size, no human protein is significantly similar to EBNA1. This table is
therefore a *hypothesis ordering*, not a set of detections.

Two views are reported because they answer different questions:

* **by E-value** — what a homology search would report. Dominated by long,
  low-identity alignments: PDE11A (E = 3.09, 29.4% over 163 aa, EBNA1 464–619),
  MAP1S (E = 7.08, 36.9% over 46 aa), PLEKHH3 (E = 37.2, 27.9% over 86 aa).
* **by local identity** — where a short, near-identical stretch sits, which is
  closer to what an antibody sees. But the very top of this view (88–91% over
  8–11 aa, all inside EBNA1 541–560) is proline/glycine-rich and is a
  composition artifact, not mimicry.

## The working set

Requiring both a stretch long enough to be non-trivial and an identity high
enough to matter — **alignment length ≥ 20 aa and identity ≥ 40%** — leaves
**61 alignments**. Twenty-six of them sit on a protein that is either
CNS-expressed or surface-exposed, i.e. reachable by an antibody. Full table:
`reports/45_candidate_table.tsv`.

Highest local identity among the reachable set:

| identity | length | E | EBNA1 | human protein | acc | EBNA1 region |
|---:|---:|---:|---|---|---|---|
| 59.0% | 22 aa | 9.3e3 | 383–404 | COBL, cordon-bleu | O75128 | SPR-rich |
| 55.0% | 20 aa | 339 | 533–552 | GPC2, glypican-2 | Q8N158 | folded core |
| 54.0% | 37 aa | 64.6 | 394–429 | MAPK8IP2, JIP-2 | Q13387 | SPR-rich |
| 50.0% | 26 aa | 257 | 530–555 | FBLIM1, migfilin | Q8WUP2 | folded core |
| 48.2% | 29 aa | 4.1e3 | 52–78 | RPH3A, rabphilin-3A | Q9Y2J0 | N-terminal |
| 47.6% | 21 aa | 112 | 540–560 | TSEN54 | Q7Z6J9 | folded core |
| 47.6% | 21 aa | 5.4e3 | 538–558 | LARGE2 | Q8N3Y3 | folded core |
| 46.4% | 28 aa | 148 | 534–558 | KNDC1 | Q76NI1 | folded core |
| 45.4% | 22 aa | 777 | 535–556 | TMEM163 | Q8TC26 | folded core |
| 45.4% | 22 aa | 4.1e3 | 541–562 | RHBDF2, iRhom2 | Q6PJF5 | folded core |
| 45.0% | 20 aa | 2.4e3 | 576–595 | PRKN, parkin | O60260 | folded core |
| 43.7% | 32 aa | 1.4e3 | 527–558 | SYNE1, nesprin-1 | Q8NF91 | folded core |
| 42.8% | 42 aa | 777 | 506–546 | DPH2 | Q9BQC3 | folded core |

## Where in EBNA1 the alignments sit

| EBNA1 region | alignments |
|---|---:|
| folded core 476–585 | 276 |
| C-terminal | 68 |
| SPR-rich 370–445 | 12 |
| N-terminal | 7 |
| Gly/Ala block 90–330 | **0** |

The Gly/Ala block contributes nothing here — not because it was masked, but
because MMseqs2's composition-bias correction removes it. The block is still the
explanation for the *other* search methods (see §6b of `REPORT.md`), but a
composition-corrected search never reaches for it.

## What this list is not

1. **Not significant.** All E > 1. A 37-aa alignment at 54% identity (MAPK8IP2)
   is what a permissive search finds by chance in a 20k-protein database.
2. **Not composition-free at the top.** View A's leaders are 8–11 aa runs inside
   a proline/glycine stretch; they are listed for completeness, not as findings.
3. **Not predictions of antibody cross-reactivity.** A shared peptide is
   necessary, never sufficient — the published GlialCAM mimic required a crystal
   structure (PDB 7K7R) to be defined.
4. **Not ranked by biology.** The order is identity, then length. CNS expression
   and surface exposure are annotated, not scored.

## What would make one of these interesting

A shared stretch in a region of EBNA1 that is actually accessible in the intact
protein, on a host protein that is reachable by extracellular antibodies, with
conserved flanking sequence that could reproduce the peptide's conformation.
MAPK8IP2 is the single alignment that satisfies the first two conditions on its
face: 37 aa at 54% in the SPR-rich region, neuronal, and predicted
surface-exposed (pLDDT 85.1, rSASA 0.26 — the only candidate segment in the
folded-and-exposed class).
