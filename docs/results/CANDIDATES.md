# EBNA1 mimicry candidates

Fifty-five human proteins carry a composition-controlled local match to EBNA1, out of 20,431 reviewed human proteins screened end to end. This document is the list. The statistics that limit how far any ordering can be pushed are in `METHOD.md`; they are repeated on every card below because they apply to every card below.

**What a row means.** A stretch of EBNA1 and a stretch of a human protein align better than chance composition would predict, after both stretches are shuffled with their own residues held fixed. That is all it means. It does not mean an antibody binds both, and nothing in this project measures binding.

**Ordering.** Shortlist = length >= 20 aa, identity >= 40%, and the human protein is CNS-expressed, antibody-accessible, or both. Inside the shortlist the first key is whether the composition null's tail was cleared at all (`null_max < observed score` over 10,000 both-sides shuffles), and only then `margin_z`. `margin_z` on its own is a convenience, not evidence: `margin_z >= 6` retains 41% of the whole pool while recovering none of the four published mimics (`reports/72_axis_validation.tsv`).

## The shortlist

| # | gene | accession | EBNA1 | human | len | identity | longest run | margin z | null tail | low-complexity | accessible | CNS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **PRKAR1B** | `P31321` | 564-587 | 20-43 | 24 | 41.6% | 4 | 10.87 | cleared | - | yes | - |
| 2 | **RHBDF2** | `Q6PJF5` | 541-562 | 339-360 | 22 | 45.4% | 3 | 7.28 | cleared | - | yes | - |
| 3 | **MAPK8IP2** | `Q13387` | 394-429 | 309-341 | 37 | 54.0% | 3 | 6.95 | cleared | - | yes | - |
| 4 | **PCDH8** | `O95206` | 530-553 | 868-891 | 24 | 41.6% | 2 | 6.48 | NOT cleared | - | yes | - |
| 5 | **CACNA1H** | `O95180` | 521-548 | 2321-2348 | 28 | 42.8% | 3 | 5.88 | NOT cleared | - | yes | - |
| 6 | **COBL** | `O75128` | 383-404 | 1000-1021 | 22 | 59.0% | 3 | 5.79 | NOT cleared | YES | yes | - |
| 7 | **SLC4A5** | `Q9BY07` | 544-576 | 914-949 | 36 | 41.6% | 3 | 5.47 | NOT cleared | - | yes | - |
| 8 | **MUC17** | `Q685J3` | 587-606 | 4440-4461 | 22 | 40.9% | 3 | 5.16 | NOT cleared | - | yes | - |
| 9 | **SCRIB** | `Q14160` | 535-556 | 1532-1553 | 22 | 40.9% | 2 | 4.86 | NOT cleared | YES | yes | - |
| 10 | **SNAP91** | `O60641` | 535-557 | 862-888 | 27 | 44.4% | 3 | 4.45 | NOT cleared | YES | yes | - |
| 11 | **TMEM163** | `Q8TC26` | 535-556 | 21-42 | 22 | 45.4% | 4 | 4.21 | NOT cleared | YES | yes | - |
| 12 | **RPH3A** | `Q9Y2J0` | 52-78 | 292-320 | 29 | 48.2% | 4 | 3.81 | NOT cleared | YES | yes | - |

`null tail` says whether the strongest of 10,000 same-composition shuffles scored **below** the observation. When it does not, the pair carries no information beyond its own residue composition — no p-value at any sample size will change that, because the null is already at the observation. **9 of the 12 cards below are in that state**, and they are ordered after the ones that are not: a pair whose null tail is not cleared cannot be defended, so it should not outrank one that can.

### 1. PRKAR1B — `P31321`

*cAMP-dependent protein kinase type I-beta regulatory subunit* (genes: PRKAR1B)

```
EBNA1   564 VFLQTHIFAEVLKDAIKDLVMTKP
            +++| |    |||| |  | ++||   | = identical, + = same residue class
PRKAR1B  20 LYVQLHGIQQVLKDCIVHLCISKP
```

- **24 aa, 41.6% identical, longest identical run 4.** Observed BLOSUM62 score 51 against a both-sides composition null of 15.56 ± 3.26; the strongest of 10,000 shuffles scored 36, i.e. 71% of the observation (margin z = **10.87**, p = 9.999e-05).
- **EBNA1 side:** partly buried, mean pLDDT 80.4; window entropy 3.574 bits, most frequent residue 12.5%.
- **Human side:** surface of folded domain, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000269|PubMed:23115245, ECO:0000269|PubMed:27028580}.
- **Tissue:** TISSUE SPECIFICITY: Four types of regulatory chains are found: I-alpha, I-beta, II-alpha, and II-beta. Their expression varies among tissues and is in some cases constitutive and in others inducible.
- **Keywords:** 3D-structure;cAMP;cAMP-binding;Cell membrane;Disease variant;Disulfide bond;Intellectual disability;Membrane;Methylation;Nitration;Nucleotide-binding;Phosphoprotein;Proteomics identification;Reference proteome;Repeat
- **HPA:** brain nTPM 319.6; `Cluster 1: Neurons - Mixed function`; single-cell `Cluster 47: Platelets - Hemostasis`.
- **Literature:** 0 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; 0 IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 15.56.

### 2. RHBDF2 — `Q6PJF5`

*Inactive rhomboid protein 2* (genes: RHBDF2, IRHOM2, RHBDL5, RHBDL6)

```
EBNA1   541 FGMAPGPGPQPGPLRESIVCYF
            +| || |||  |    | |  |   | = identical, + = same residue class
RHBDF2  339 YGRAPVPGPRRGKRIASKVKHF
```

- **22 aa, 45.4% identical, longest identical run 3.** Observed BLOSUM62 score 44 against a both-sides composition null of 17.41 ± 3.65; the strongest of 10,000 shuffles scored 42, i.e. 95% of the observation (margin z = **7.28**, p = 9.999e-05).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 51.7; window entropy 3.477 bits, most frequent residue 22.7%.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Endoplasmic reticulum membrane {ECO:0000250|UniProtKB:Q80WQ6}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q80WQ6}. Cell membrane {ECO:0000269|PubMed:22265016, ECO:0000269|PubMed:29897333}.
- **Tissue:** TISSUE SPECIFICITY: Found in the epidermis and esophageal epithelium. {ECO:0000269|PubMed:22265016}.
- **Keywords:** 3D-structure;Alternative splicing;Cell membrane;Disease variant;Endoplasmic reticulum;Growth factor binding;Membrane;Palmoplantar keratoderma;Phosphoprotein;Protein transport;Proteomics identification;Reference proteome;
- **HPA:** `Cluster 20: Macrophages & Microglia - Immune response`; single-cell `Cluster 74: Microglia - Neuroimmunity & neuroinflammation`.
- **Literature:** 0 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; 0 IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 17.41.

### 3. MAPK8IP2 — `Q13387`

*C-Jun-amino-terminal kinase-interacting protein 2* (genes: MAPK8IP2, IB2, JIP2, PRKM8IPL)

```
EBNA1   394 PPRRPPPGRRPFFHPVGEAD-YFEYHQEGGPDGEPDV
            ||| ||  ||| | |||  |   ||  | | + |||+   | = identical, + = same residue class
MAPK8IP2 309 PPREPP--RRPAFLPVGPDDTNSEY--ESGSESEPDL
```

- **37 aa, 54.0% identical, longest identical run 3.** Observed BLOSUM62 score 65 against a both-sides composition null of 27.83 ± 5.35; the strongest of 10,000 shuffles scored 58, i.e. 89% of the observation (margin z = **6.95**, p = 9.999e-05).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 27.7; window entropy 3.18 bits, most frequent residue 25.0%. **This segment overlaps a validated epitope region.**
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:10490659}. Note=Accumulates in cell surface projections. {ECO:0000269|PubMed:10490659}.
- **Tissue:** TISSUE SPECIFICITY: Expressed mainly in the brain and pancreas, including insulin-secreting cells. In the nervous system, more abundantly expressed in the cerebellum, pituitary gland, occipital lobe and the amygdala. Als
- **Keywords:** 3D-structure;Alternative splicing;Cytoplasm;Phosphoprotein;Proteomics identification;Reference proteome;SH3 domain
- **HPA:** brain nTPM 103.2; `Cluster 24: Neurons - Mixed function`; single-cell `Cluster 2: Cerebellar granular cells - Neuronal signaling`.
- **Literature:** 0 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; 0 IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 27.83.

### 4. PCDH8 — `O95206`

*Protocadherin-8* (genes: PCDH8)

```
EBNA1   530 QCRLTPLSRLPFGMAPGPGPQPGP
            | || ++   |+|+ || |  | |   | = identical, + = same residue class
PCDH8   868 QQRLRGAHAEPYGASPGFGKEPAP
```

- **24 aa, 41.6% identical, longest identical run 2.** Observed BLOSUM62 score 48 against a both-sides composition null of 20.6 ± 4.23; the strongest of 10,000 shuffles scored 50, **at or above the observation** — the null tail is not cleared at all (margin z = **6.48**, p = 0.00029997).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 57.1; window entropy 3.068 bits, most frequent residue 29.2%.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000250}; Single-pass type I membrane protein {ECO:0000250}. Cell projection, dendrite {ECO:0000250}. Presynaptic cell membrane. Postsynaptic cell membrane. Note=Also expressed in neuronal cell bodies. Localized to excitatory, but not with inhibitory, synapse
- **Keywords:** 3D-structure;Alternative splicing;Calcium;Cell adhesion;Cell membrane;Cell projection;Glycoprotein;Membrane;Phosphoprotein;Postsynaptic cell membrane;Proteomics identification;Reference proteome;Repeat;Signal;Synapse;Tra
- **HPA:** brain nTPM 12.3; `Cluster 51: Neurons - Mixed function`; single-cell `Cluster 35: Neurons - Synaptic function & neurogenesis`.
- **Literature:** 0 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; 0 IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 20.6.

### 5. CACNA1H — `O95180`

*Voltage-dependent T-type calcium channel subunit alpha-1H* (genes: CACNA1H)

```
EBNA1   521 RRGTALAIPQCRLTPLSRLPFGMAPGPG
            |||  | +||| |          |||+|   | = identical, + = same residue class
CACNA1H2321 RRGLYLTVPQCPLEKPGSPSATPAPGGG
```

- **28 aa, 42.8% identical, longest identical run 3.** Observed BLOSUM62 score 48 against a both-sides composition null of 21.96 ± 4.43; the strongest of 10,000 shuffles scored 49, **at or above the observation** — the null tail is not cleared at all (margin z = **5.88**, p = 0.00029997).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 68.6; window entropy 3.294 bits, most frequent residue 17.9%.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000269|PubMed:25907736, ECO:0000269|PubMed:27149520, ECO:0000269|PubMed:27729216, ECO:0000269|PubMed:9670923, ECO:0000269|PubMed:9930755}; Multi-pass membrane protein {ECO:0000305}. Note=Interaction with STAC increases expression at the cell membrane. {ECO:0
- **Tissue:** TISSUE SPECIFICITY: Expressed in the adrenal glomerulosa (at protein level) (PubMed:25907736, PubMed:27729216). In nonneuronal tissues, the highest expression levels are found in the kidney, liver, and heart. In the brai
- **Keywords:** 3D-structure;Alternative splicing;Calcium;Calcium channel;Calcium transport;Cell membrane;Disease variant;Epilepsy;Glycoprotein;Ion channel;Ion transport;Membrane;Metal-binding;Proteomics identification;Reference proteom
- **HPA:** `Cluster 33: Neurons - Mixed function`; single-cell `Cluster 110: Mural cells - ECM organization`.
- **Literature:** 1 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 21.96.

### 6. COBL — `O75128`

*Protein cordon-bleu* (genes: COBL, KIAA0633)

```
EBNA1   383 SPSSQSSSSGSPPRRPPPGRRP
            | |||  || | ||| | |  |   | = identical, + = same residue class
COBL   1000 STSSQEASSASEPRRAPDGTDP
```

- **22 aa, 59.0% identical, longest identical run 3.** Observed BLOSUM62 score 60 against a both-sides composition null of 29.74 ± 5.23; the strongest of 10,000 shuffles scored 60, **at or above the observation** — the null tail is not cleared at all (margin z = **5.79**, p = 0.00019998).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 28.9; window entropy 2.021 bits, most frequent residue 36.4% — **low complexity: this match is partly a composition effect**. **This segment overlaps a validated epitope region.**
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000250}; Peripheral membrane protein {ECO:0000250}; Cytoplasmic side {ECO:0000250}. Cytoplasm, cytoskeleton {ECO:0000250}. Cell projection, ruffle {ECO:0000250}. Cytoplasm {ECO:0000250}. Note=Recruited to the cell membrane via interaction with PACSIN1. Coloc
- **Keywords:** 3D-structure;Actin-binding;Alternative splicing;Cell membrane;Cell projection;Cytoplasm;Cytoskeleton;Membrane;Phosphoprotein;Proteomics identification;Reference proteome;Repeat
- **HPA:** `Cluster 11: White matter - Signal transduction`; single-cell `Cluster 51: Secretory epithelial cells - Mixed function`.
- **Literature:** 4 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 29.74.

### 7. SLC4A5 — `Q9BY07`

*Electrogenic sodium bicarbonate cotransporter 4* (genes: SLC4A5, NBC4)

```
EBNA1   544 APGPGPQPGPLRESIVCYFMVFLQTHI---FAEVLK
            |||  ||  ++||+ |   +||+ | |    | +||   | = identical, + = same residue class
SLC4A5  914 APGEQPQFLGVREQRVTGIIVFILTGISVFLAPILK
```

- **36 aa, 41.6% identical, longest identical run 3.** Observed BLOSUM62 score 42 against a both-sides composition null of 20.08 ± 4.0; the strongest of 10,000 shuffles scored 44, **at or above the observation** — the null tail is not cleared at all (margin z = **5.47**, p = 0.00069993).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 65.5; window entropy 3.874 bits, most frequent residue 15.2%.
- **Human side:** partly buried, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Apical cell membrane {ECO:0000269|PubMed:17715183, ECO:0000269|PubMed:26447209}; Multi-pass membrane protein {ECO:0000255}. Basolateral cell membrane {ECO:0000250|UniProtKB:Q6RI88}; Multi-pass membrane protein {ECO:0000255}. Note=Expressed in the apical plasma membrane domain o
- **Tissue:** TISSUE SPECIFICITY: Highest expression observed in liver, spleen and testis; moderate expression in the choroid plexus, hippocampus, cerebrum and cerebellum of brain, and in kidney cortex and kidney medulla. Also observe
- **Keywords:** Alternative splicing;Anion exchange;Cell membrane;Ion transport;Membrane;Proteomics identification;Reference proteome;Sodium;Sodium transport;Transmembrane;Transmembrane helix;Transport
- **HPA:** `Cluster 36: Choroid plexus - Mixed function`; single-cell `Cluster 101: Choroid plexus epithelial cells - Ion transport`.
- **Literature:** 0 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 20.08.

### 8. MUC17 — `Q685J3`

*Mucin-17* (genes: MUC17, MUC3)

```
EBNA1   587 PAPTC--NIRVTVCSFDDGVDL
            |||    ||   +|+ || + |   | = identical, + = same residue class
MUC17  4440 PAPGTFQNIGFDICQDDDSIHL
```

- **22 aa, 40.9% identical, longest identical run 3.** Observed BLOSUM62 score 34 against a both-sides composition null of 16.11 ± 3.47; the strongest of 10,000 shuffles scored 36, **at or above the observation** — the null tail is not cleared at all (margin z = **5.16**, p = 0.00059994).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 66.2; window entropy 3.546 bits, most frequent residue 15.0%.
- **Human side:** not modelled, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: [Isoform 1]: Cell membrane; Single-pass membrane protein.; SUBCELLULAR LOCATION: [Isoform 2]: Secreted {ECO:0000305}. Cell membrane.
- **Tissue:** TISSUE SPECIFICITY: Expressed almost exclusively in the intestine. Expression is especially high in both the duodenum and transverse colon. Expressed in mature absorptive cells of the small intestinal villi. No expressio
- **Keywords:** Alternative splicing;Cell membrane;Disulfide bond;Glycoprotein;Membrane;Proteomics identification;Reference proteome;Repeat;Secreted;Signal;Transmembrane;Transmembrane helix
- **HPA:** single-cell `Cluster 1: Enterocytes - Absorption and Digestion`; secretome: Intracellular and membrane.
- **Literature:** 1 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 16.11.

### 9. SCRIB — `Q14160`

*Protein scribble homolog* (genes: SCRIB, CRIB1, KIAA0147, LAP4, SCRB1, VARTUL)

```
EBNA1   535 PLSRLPFGMAPGPGPQPGPLRE
            || ||     | | |+| |+ +   | = identical, + = same residue class
SCRIB  1532 PLERLAEAPSPAPTPSPTPVED
```

- **22 aa, 40.9% identical, longest identical run 2.** Observed BLOSUM62 score 49 against a both-sides composition null of 25.17 ± 4.9; the strongest of 10,000 shuffles scored 52, **at or above the observation** — the null tail is not cleared at all (margin z = **4.86**, p = 0.00069993).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 51.3; window entropy 2.896 bits, most frequent residue 31.8% — **low complexity: this match is partly a composition effect**.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000269|PubMed:19386766, ECO:0000269|PubMed:21849460, ECO:0000269|PubMed:28169360}; Peripheral membrane protein. Cell junction {ECO:0000269|PubMed:27380321}. Cell junction, adherens junction. Cell projection, lamellipodium. Cytoplasm {ECO:0000269|PubMed:19386
- **Tissue:** TISSUE SPECIFICITY: Expressed in kidney, skeletal muscles, liver, lung, breast, intestine, placenta and skin mainly in epithelial cells (at protein level). {ECO:0000269|PubMed:15806148}.
- **Keywords:** 3D-structure;Alternative splicing;Cell junction;Cell membrane;Cell projection;Coiled coil;Cytoplasm;Developmental protein;Differentiation;Disease variant;Leucine-rich repeat;Lipoprotein;Membrane;Palmitate;Phosphoprotein;
- **HPA:** `Cluster 13: Oligodendrocytes - Mixed function`; single-cell `Cluster 28: Astrocytes - Homeostasis & neuron support`.
- **Literature:** 2 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 25.17.

### 10. SNAP91 — `O60641`

*Clathrin coat assembly protein AP180* (genes: SNAP91, KIAA0656)

```
EBNA1   535 PLSRLPFGMAPGPG----PQPGPLRES
            |+ | |||+|  ||    |+| |+  |   | = identical, + = same residue class
SNAP91  862 PMMRPPFGAAAVPGTQLSPSPTPASQS
```

- **27 aa, 44.4% identical, longest identical run 3.** Observed BLOSUM62 score 49 against a both-sides composition null of 26.56 ± 5.04; the strongest of 10,000 shuffles scored 59, **at or above the observation** — the null tail is not cleared at all (margin z = **4.45**, p = 0.00069993).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 51.7; window entropy 2.941 bits, most frequent residue 30.4% — **low complexity: this match is partly a composition effect**.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cell membrane {ECO:0000250}. Membrane, coated pit {ECO:0000250}; Peripheral membrane protein {ECO:0000250}; Cytoplasmic side {ECO:0000250}. Note=Component of the coat surrounding the cytoplasmic face of coated vesicles in the plasma membrane. {ECO:0000250}.
- **Keywords:** Alternative splicing;Cell membrane;Coated pit;Glycoprotein;Membrane;Methylation;Phosphoprotein;Protein transport;Proteomics identification;Reference proteome;Transport
- **HPA:** brain nTPM 105.4; `Cluster 52: Neurons & Synapses - Synaptic function`; single-cell `Cluster 59: Pituitary endocrine cells - Neuroendocrine signaling`.
- **Literature:** 1 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 26.56.

### 11. TMEM163 — `Q8TC26`

*Transmembrane protein 163* (genes: TMEM163)

```
EBNA1   535 PLSRLPFGMAPGPGPQPGPLRE
            |  ++|  +|||| |   |+||   | = identical, + = same residue class
TMEM163  21 PRGHAPPAAAPGPAPLSSPVRE
```

- **22 aa, 45.4% identical, longest identical run 4.** Observed BLOSUM62 score 52 against a both-sides composition null of 29.37 ± 5.37; the strongest of 10,000 shuffles scored 57, **at or above the observation** — the null tail is not cleared at all (margin z = **4.21**, p = 0.00129987).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 51.3; window entropy 2.896 bits, most frequent residue 31.8% — **low complexity: this match is partly a composition effect**.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cytoplasmic vesicle, secretory vesicle, synaptic vesicle membrane {ECO:0000250|UniProtKB:Q8C996}; Multi-pass membrane protein {ECO:0000255}. Early endosome membrane {ECO:0000250|UniProtKB:A9CMA6}; Multi-pass membrane protein {ECO:0000255}. Late endosome membrane {ECO:0000269|Pu
- **Tissue:** TISSUE SPECIFICITY: Widely expressed. High expression is detected in brain, lung and testis. {ECO:0000269|PubMed:25130899}.
- **Keywords:** Alternative splicing;Cell membrane;Cytoplasmic vesicle;Disease variant;Endosome;Leukodystrophy;Lysosome;Membrane;Phosphoprotein;Proteomics identification;Reference proteome;Synapse;Transmembrane;Transmembrane helix;Trans
- **HPA:** brain nTPM 44.4; `Cluster 38: Sub-cortical - Mixed function`; single-cell `Cluster 14: Oligodendrocyte progenitor cells - Development & myelination`.
- **Literature:** 2 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 29.37.

### 12. RPH3A — `Q9Y2J0`

*Rabphilin-3A* (genes: RPH3A, KIAA0985)

```
EBNA1    52 GRPGAPGGS--GSGPRHRDGVRRPQKRPS
            | || ||||  | ||  | +  +|   ||   | = identical, + = same residue class
RPH3A   292 GQPGTPGGSRPGPGPAGRFPDQKPEVAPS
```

- **29 aa, 48.2% identical, longest identical run 4.** Observed BLOSUM62 score 53 against a both-sides composition null of 31.23 ± 5.71; the strongest of 10,000 shuffles scored 67, **at or above the observation** — the null tail is not cleared at all (margin z = **3.81**, p = 0.00279972).
- **EBNA1 side:** disordered or loop exposed, mean pLDDT 34.1; window entropy 2.847 bits, most frequent residue 25.9% — **low complexity: this match is partly a composition effect**.
- **Human side:** disordered or loop exposed, maturation `retained_in_chain`, antibody-accessible compartment.
- **Location:** SUBCELLULAR LOCATION: Cytoplasmic vesicle, secretory vesicle, synaptic vesicle membrane {ECO:0000250|UniProtKB:P47709}. Cell projection, dendritic spine {ECO:0000250|UniProtKB:P47709}. Postsynaptic cell membrane {ECO:0000250|UniProtKB:P47709}. Membrane {ECO:0000250|UniProtKB:P47709}; Peripheral memb
- **Keywords:** Alternative splicing;Calcium;Cell membrane;Cell projection;Cytoplasmic vesicle;Lipid-binding;Membrane;Metal-binding;Methylation;Phosphoprotein;Postsynaptic cell membrane;Protein transport;Proteomics identification;Refere
- **HPA:** brain nTPM 115.4; `Cluster 38: Sub-cortical - Mixed function`; single-cell `Cluster 48: Neutrophils - Phagocytosis & degranulation`.
- **Literature:** 1 PubMed records for this gene together with autoantibody / autoantigen / multiple sclerosis / EBV terms; not_queried IEDB B-cell records overlapping the matched segment.
- **Unknowns:** no antibody has been tested against these two peptides; no PTM state modelled; no expression data in MS lesions or CSF; disorder calls come from one AlphaFold conformer.
- **The experiment:** order both peptides, run a competition ELISA for the human peptide against whatever EBNA1 reactivity is available, or immunise and test cross-reactivity. A negative result is informative because the pair is a pure composition case: the null already says a shuffle of these same residues scores 31.23.

## The full list of 55

| gene | accession | EBNA1 | human | len | identity | longest run | margin z | mature | accessible | CNS | low-complexity |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PRKAR1B | `P31321` | 564-587 | 20-43 | 24 | 41.6% | 4 | 10.87 | yes | yes | - | - |
| RHBDF2 | `Q6PJF5` | 541-562 | 339-360 | 22 | 45.4% | 3 | 7.28 | yes | yes | - | - |
| MAPK8IP2 | `Q13387` | 394-429 | 309-341 | 37 | 54.0% | 3 | 6.95 | yes | yes | - | - |
| PCDH8 | `O95206` | 530-553 | 868-891 | 24 | 41.6% | 2 | 6.48 | yes | yes | - | - |
| CACNA1H | `O95180` | 521-548 | 2321-2348 | 28 | 42.8% | 3 | 5.88 | yes | yes | - | - |
| COBL | `O75128` | 383-404 | 1000-1021 | 22 | 59.0% | 3 | 5.79 | yes | yes | - | YES |
| SLC4A5 | `Q9BY07` | 544-576 | 914-949 | 36 | 41.6% | 3 | 5.47 | yes | yes | - | - |
| MUC17 | `Q685J3` | 587-606 | 4440-4461 | 22 | 40.9% | 3 | 5.16 | yes | yes | - | - |
| SCRIB | `Q14160` | 535-556 | 1532-1553 | 22 | 40.9% | 2 | 4.86 | yes | yes | - | YES |
| SNAP91 | `O60641` | 535-557 | 862-888 | 27 | 44.4% | 3 | 4.45 | yes | yes | - | YES |
| TMEM163 | `Q8TC26` | 535-556 | 21-42 | 22 | 45.4% | 4 | 4.21 | yes | yes | - | YES |
| RPH3A | `Q9Y2J0` | 52-78 | 292-320 | 29 | 48.2% | 4 | 3.81 | yes | yes | - | YES |
| ZNF347 | `Q96SE7` | 566-586 | 230-250 | 21 | 42.8% | 3 | 9.23 | yes | - | - | - |
| RNF228 | `A0A7I2V3R4` | 63-83 | 207-227 | 21 | 52.3% | 3 | 9.10 | yes | - | - | - |
| R3HDM4 | `Q96D70` | 478-503 | 169-193 | 26 | 46.1% | 3 | 8.70 | yes | - | - | - |
| CARD17P | `Q5XLA6` | 578-601 | 56-80 | 25 | 44.0% | 2 | 8.59 | yes | - | - | - |
| RFLNB | `Q8N5W9` | 588-614 | 78-104 | 27 | 40.7% | 2 | 8.40 | yes | - | - | - |
| DALRD3 | `Q5D0E6` | 401-421 | 424-444 | 21 | 42.8% | 2 | 8.34 | yes | - | - | - |
| ACACA | `Q13085` | 488-511 | 1938-1958 | 24 | 41.6% | 2 | 8.24 | yes | - | - | - |
| BLOC1S1 | `P78537` | 516-539 | 10-33 | 24 | 45.8% | 4 | 8.15 | yes | - | - | - |
| HECTD3 | `Q5T447` | 584-606 | 342-363 | 23 | 43.4% | 3 | 7.73 | yes | - | - | - |
| TMEM179 | `Q6ZVK1` | 477-497 | 72-92 | 21 | 42.8% | 3 | 7.47 | yes | - | - | - |
| TSEN54 | `Q7Z6J9` | 540-560 | 257-277 | 21 | 47.6% | 3 | 7.28 | yes | - | - | - |
| ACACB | `O00763` | 488-514 | 2049-2072 | 27 | 40.7% | 3 | 6.99 | yes | - | - | - |
| CCHCR1 | `Q8TD31` | 482-502 | 143-164 | 22 | 54.5% | 6 | 6.34 | yes | - | - | - |
| FBLIM1 | `Q8WUP2` | 530-555 | 135-160 | 26 | 50.0% | 6 | 6.30 | yes | - | - | - |
| DPH2 | `Q9BQC3` | 506-546 | 215-244 | 42 | 42.8% | 4 | 6.02 | yes | - | - | - |
| SYNE1 | `Q8NF91` | 527-558 | 5267-5297 | 32 | 43.7% | 2 | 7.52 | yes | - | - | - |
| IRF3 | `Q14653` | 580-606 | 360-393 | 35 | 40.0% | 3 | 6.36 | yes | - | - | - |
| WDFY3 | `Q8IZQ1` | 557-583 | 311-336 | 27 | 40.7% | 4 | 6.20 | yes | - | - | - |
| PRKN | `O60260` | 576-595 | 408-425 | 20 | 45.0% | 3 | 6.12 | yes | - | - | - |
| PTDSS2 | `Q9BVG9` | 416-450 | 5-40 | 37 | 43.2% | 3 | 5.71 | yes | - | - | - |
| OPA1 | `O60313` | 592-615 | 151-174 | 24 | 41.6% | 2 | 5.63 | yes | - | - | - |
| RPRD1B | `Q9NQG5` | 476-501 | 75-99 | 26 | 42.3% | 3 | 5.57 | yes | - | - | - |
| CBX2 | `Q14781` | 537-556 | 334-350 | 20 | 55.0% | 5 | 5.55 | yes | - | - | YES |
| NFKB1 | `P19838` | 479-502 | 696-720 | 25 | 40.0% | 2 | 5.44 | yes | - | - | - |
| KNDC1 | `Q76NI1` | 534-558 | 723-750 | 28 | 46.4% | 6 | 5.39 | yes | - | - | - |
| CMBL | `Q96DG6` | 505-524 | 150-171 | 22 | 54.5% | 3 | 5.20 | yes | - | - | - |
| ERICH2 | `A1L162` | 474-493 | 57-82 | 26 | 46.1% | 4 | 5.17 | yes | - | - | - |
| GALT | `P07902` | 544-565 | 106-131 | 26 | 42.3% | 3 | 5.16 | yes | - | - | - |
| ABTB3 | `A6QL63` | 573-594 | 766-786 | 22 | 40.9% | 2 | 5.10 | yes | - | - | - |
| TPRX1 | `Q8N7U7` | 529-554 | 80-103 | 26 | 42.3% | 5 | 5.06 | yes | - | - | YES |
| NCAPH2 | `Q6IBW4` | 531-554 | 230-251 | 24 | 45.8% | 3 | 5.02 | yes | - | - | - |
| PTGS1 | `P23219` | 574-603 | 557-589 | 33 | 42.4% | 3 | 4.93 | yes | - | - | - |
| ACTL9 | `Q8TC94` | 482-506 | 356-384 | 29 | 41.3% | 4 | 4.83 | yes | - | - | - |
| AHCYL2 | `Q96HN2` | 425-448 | 48-71 | 24 | 41.6% | 2 | 4.82 | yes | - | - | YES |
| ATAD3B | `Q5T9A4` | 602-624 | 14-36 | 23 | 43.4% | 2 | 4.71 | yes | - | - | - |
| ULK1 | `O75385` | 521-552 | 648-682 | 35 | 42.8% | 3 | 4.69 | yes | - | - | - |
| ATAD3A | `Q9NVI7` | 602-624 | 14-36 | 23 | 43.4% | 2 | 4.62 | yes | - | - | - |
| TACC2 | `O95359` | 487-525 | 2699-2732 | 39 | 41.0% | 3 | 4.56 | yes | - | - | - |
| ZBTB5 | `O15062` | 531-553 | 438-458 | 23 | 43.4% | 3 | 4.02 | yes | - | - | YES |
| PPP2R3B | `Q9Y5P8` | 397-417 | 542-565 | 24 | 45.8% | 4 | 3.89 | yes | - | - | - |
| KDM6B | `O15054` | 534-557 | 627-652 | 26 | 46.1% | 4 | 3.74 | yes | - | - | - |
| LARGE2 | `Q8N3Y3` | 538-558 | 36-56 | 21 | 47.6% | 4 | 3.74 | yes | - | - | - |
| TTLL10 | `Q6ZVT0` | 528-556 | 42-71 | 30 | 43.3% | 2 | 3.57 | yes | - | - | - |

## The normal-composition subset

Removing every pair whose EBNA1 window is itself low complexity (entropy < 3.0 bits or one residue above 30%) leaves **46 of 55**. This is the subset where the match is not explained by repeated residues on our side of the pair.

| gene | accession | EBNA1 | human | len | identity | margin z | accessible | CNS |
|---|---|---|---|---|---|---|---|---|
| PRKAR1B | `P31321` | 564-587 | 20-43 | 24 | 41.6% | 10.87 | yes | - |
| ZNF347 | `Q96SE7` | 566-586 | 230-250 | 21 | 42.8% | 9.23 | - | - |
| RNF228 | `A0A7I2V3R4` | 63-83 | 207-227 | 21 | 52.3% | 9.10 | - | - |
| R3HDM4 | `Q96D70` | 478-503 | 169-193 | 26 | 46.1% | 8.70 | - | - |
| CARD17P | `Q5XLA6` | 578-601 | 56-80 | 25 | 44.0% | 8.59 | - | - |
| RFLNB | `Q8N5W9` | 588-614 | 78-104 | 27 | 40.7% | 8.40 | - | - |
| DALRD3 | `Q5D0E6` | 401-421 | 424-444 | 21 | 42.8% | 8.34 | - | - |
| ACACA | `Q13085` | 488-511 | 1938-1958 | 24 | 41.6% | 8.24 | - | - |
| BLOC1S1 | `P78537` | 516-539 | 10-33 | 24 | 45.8% | 8.15 | - | - |
| HECTD3 | `Q5T447` | 584-606 | 342-363 | 23 | 43.4% | 7.73 | - | - |
| SYNE1 | `Q8NF91` | 527-558 | 5267-5297 | 32 | 43.7% | 7.52 | - | - |
| TMEM179 | `Q6ZVK1` | 477-497 | 72-92 | 21 | 42.8% | 7.47 | - | - |
| TSEN54 | `Q7Z6J9` | 540-560 | 257-277 | 21 | 47.6% | 7.28 | - | - |
| RHBDF2 | `Q6PJF5` | 541-562 | 339-360 | 22 | 45.4% | 7.28 | yes | - |
| ACACB | `O00763` | 488-514 | 2049-2072 | 27 | 40.7% | 6.99 | - | - |
| MAPK8IP2 | `Q13387` | 394-429 | 309-341 | 37 | 54.0% | 6.95 | yes | - |
| PCDH8 | `O95206` | 530-553 | 868-891 | 24 | 41.6% | 6.48 | yes | - |
| IRF3 | `Q14653` | 580-606 | 360-393 | 35 | 40.0% | 6.36 | - | - |
| CCHCR1 | `Q8TD31` | 482-502 | 143-164 | 22 | 54.5% | 6.34 | - | - |
| FBLIM1 | `Q8WUP2` | 530-555 | 135-160 | 26 | 50.0% | 6.30 | - | - |
| WDFY3 | `Q8IZQ1` | 557-583 | 311-336 | 27 | 40.7% | 6.20 | - | - |
| PRKN | `O60260` | 576-595 | 408-425 | 20 | 45.0% | 6.12 | - | - |
| DPH2 | `Q9BQC3` | 506-546 | 215-244 | 42 | 42.8% | 6.02 | - | - |
| CACNA1H | `O95180` | 521-548 | 2321-2348 | 28 | 42.8% | 5.88 | yes | - |
| PTDSS2 | `Q9BVG9` | 416-450 | 5-40 | 37 | 43.2% | 5.71 | - | - |

**A tension to state rather than hide.** All three published linear mimics are low-complexity pairs on our side as well (`PPRRPP`, `RRPFF`, `PGAIEQGPAD`), so this filter removes them too. That is not an argument against the filter: it is the reason the published cases could not be recovered by any composition-controlled statistic. But it means the two lists answer different questions — the full 55 asks *what does EBNA1 resemble*, the normal-composition subset asks *where does the resemblance survive removing repeated residues*. Pick the one that matches the experiment you intend to run.

## What would make this list smaller and better

Nothing in the pipeline currently measures binding, so no ordering here can be improved by more sequence analysis. The three things that would change the list are: a peptide-binding assay for any of these pairs; MS/CSF reactivity data against the human peptides; and an antibody against one of the human proteins, which would let the EBNA1 peptide be tested as a competitor. The only pair in this project with solved structures on both sides of a proven cross-reaction is EBNA1 386-405 / GlialCAM 370-389 (7K7R, 7UQC), and that pair scores `margin_z` 0.865 under the same null.
