# Step 9 — the paratope channel: one antibody, two peptides, and what it actually reads

**Headline (revised after §4b).** Two crystal structures of the same antibody
bound to two different peptides — the only experimentally proven molecular
mimicry case in this project — show that the peptides are **30% identical at
their true binding register (+16)** and that the antibody's contact mass is
**spread over about ten peptide residues**, the largest single contributor being
an arginine carrying ~20% of the contacts. A panel of 50 independent
peptide–antibody complexes puts that concentration **below the median**.

An earlier version of this report claimed the interface sat "almost entirely on
one arginine" (51% of contacts). That number was produced by merging the light
and heavy chain dictionaries with `{**light, **heavy}`, which silently discards
the heavy chain wherever residue numbers collide — 212 of 212 residues here.
Corrected: 56 of 262 contacts, not 51 of 100. The claim is retracted in §4b;
what survives is the register, the identity at that register, and the census.

---

## 0. Why this and not something else

The user's two options were **A** (change the evidence type: affinity/antibody
layer) and **C** (stop and write it up). I took A, but the cheap and falsifiable
end of it rather than a GPU co-folding run.

Rationale: the question that all previous work could not touch is *"can the same
paratope bind both peptides?"* Every axis we built measured sequence similarity,
and every one of them was falsified by its own positive controls (§8 of the
dossier). Installing AlphaFold3/Boltz would let us *predict* a complex, but with
no calibration data it would generate plausible structures for anything. The
crystal structures, by contrast, are measurements.

And a decisive piece of luck: the antibody **MS39p2w174** has been crystallised
twice, against **both** of its known peptide ligands.

| file | entry | content | resolution | sha256 |
|---|---|---|---|---|
| `data/cache/pdb/7k7r.pdb` | 7K7R | Fab MS39p2w174 + EBNA1 386-405 | 2.50 Å | `da75e257f2c993ecf9f7307086818f44fa87d183ba0cc63b570bb6c905f8edab` |
| `data/cache/pdb/7uqc.pdb` | 7UQC | Fab MS39p2w174 + phospho-GlialCAM 370-389 | 2.65 Å | `5e786188d000b9ad91a960f063c53d0686cc297b2fd9e370244137659f3d1bc3` |

Fetched from `https://files.rcsb.org/download/{ID}.pdb` on 2026-10-22T10:22:24Z.
Scripts: `src/probe_paratope.py`, `src/probe_read_motif.py`.

---

## 1. Result: it is the same paratope

Superposing the Fv (chains A+B, 429 CA atoms) of the two crystals:

```
Fv CA superposition: 429 residues, RMSD 0.63 A
```

0.63 Å over 429 residues means the two structures present essentially the same
combining site. Any comparison of the two peptides is therefore a comparison
within one paratope, not between two antibodies.

## 2. Result: the binding register is +16

Scanning every offset, after superposing the two Fabs:

| offset | mean CA distance | n |
|---|---|---|
| +14 | 7.88 Å | 14 |
| +15 | 5.34 Å | 14 |
| **+16** | **2.48 Å** | **14** |
| +17 | 3.21 Å | 14 |
| +18 | 5.16 Å | 14 |

A clean, unambiguous minimum: **GlialCAM residue *r* sits where EBNA1 residue
*r+16* sits.** The register is a property of the structures, not of an alignment
program — which matters, because a sequence aligner has no way to know it.

Note also what the crystals actually model: 7K7R resolves **EBNA1 389-402 only**.
The published "386-405" epitope has six residues (386-388, 403-405) unmodelled,
so the epitope's ordered core is 14 residues, not 20.

## 3. Result: 30% identity is enough

At the true register the two bound peptides are:

```
GliaCAM  R373/S389  T374/S390  H375/S391  S376/G392  S377/S393=  P378/P394=
         P379/P395=  R380/R396=  A381/R397  P382/P398=  S383/P399  S384/P400
         P385/G401   G386/R402
```

**5 of 14 modelled positions identical (36%); 6 of 20 over the full peptides
(30%)** — the latter reproduces the figure in the source paper (Lanz 2022).

The published pairing (EBNA1 386-405 with GlialCAM 370-389) **is** register
+16, so our composition-matched pair null was measuring the right pairing all
along: `margin_z = 0.865, p = 0.200`, obs 35 versus a null mean of 30.32 (dossier
§9, `reports/55_published_pair_null.tsv`). The proven cross-reactive pair is
**indistinguishable from a composition shuffle** on sequence grounds. Both
statements are true at once, and that is the whole problem in one line.

## 4. Result: the interface is distributed, not concentrated

Number of Fab heavy atoms within 4.5 Å of each peptide residue, counted against
**all** antibody chains (light and heavy, every copy), at the shared register:

| position | EBNA1 | GlialCAM |
|---|---|---|
| +3 | S389: 0 | R373: 40 |
| +4 | S390: 16 | T374: 2 |
| +5 | S391: 21 | H375: 8 |
| +6 | G392: 10 | S376: 33 |
| +7 | S393: 17 | S377: 22 |
| +8 | P394: 20 | P378: 19 |
| +9 | P395: 44 | P379: 35 |
| **+10** | **R396: 56** | **R380: 59** |
| +11 | R397: 36 | A381: 28 |
| +12 | P398: 38 | P382: 38 |
| +13 | P399: 4 | S383: 2 |
| +0…+2, +14…+20 | 0 (or unmodelled) | 0 |

Total contacts: **EBNA1 262, GlialCAM 290**.

The largest single contributor is an arginine — R396 at **21.4%** of the EBNA1
interface, R380 at **20.3%** of the GlialCAM one — but it is one of ten
contacting residues, and the next three (P395 16.8%, P398 14.5%, R397 13.7% in
EBNA1) are comparable. The source paper's "two central arginines" survives as
"both arginines are contacted"; it does not survive as "one arginine carries
the interface", which is what an earlier version of this section said.

## 4b. The concentration claim, tested against 50 other complexes, and retracted

If "one residue carries the interface" were a property of paratopes, it should
be visible in other peptide–antibody crystals. Panel: every RCSB entry holding
a 5–25 aa polymer entity alongside a SAbDab-annotated antibody — 52 candidates,
**50 usable** complexes, none of them EBNA1 (`src/panel_paratope.py`,
`reports/panel_paratope.tsv`, thresholds declared in the script and applied
uniformly).

| statistic | panel value |
|---|---|
| peptide residues in contact (median) | 10 of 12 (83%) |
| total contacts (median) | 296 |
| **dominant-residue share (median)** | **0.244** (range 0.159–0.617) |
| complexes with dominant share ≥ 0.30 | 6 / 50 |
| 7K7R (EBNA1) | 0.214 → **14th percentile** |
| 7UQC (GlialCAM) | 0.199 → **10th percentile** |

**The hypothesis is refuted and its reverse is what the data show.** Interfaces
are spread over most of the peptide, and our two calibration crystals are
*less* concentrated than typical, not more.

Three consequences, stated in order of how much they cost:

1. **Retracted:** "the paratope reads one arginine"; "the second arginine can be
   replaced by alanine"; "the contact mass is concentrated on one or two
   residues". None of these survive correct counting.
2. **Downgraded:** the "3-position read `P.R.P`" of §5. With ten contacting
   residues, a three-position abstraction is an operational summary I chose, not
   a description of what the antibody reads. The census numbers are unchanged
   and remain facts; their interpretation as "the read" does not.
3. **Unchanged:** the Fv RMSD (0.63 Å), the register (+16), the identity at that
   register (6/20), and — most importantly — the composition-null result for the
   published pair (margin z 0.865, p = 0.20).

The honest summary of this round is therefore not "here is why the sequence axes
failed" but "I proposed an explanation for why the sequence axes failed, tested
it in advance against an independent panel, and it did not hold". The negative
result stands on measurement; the mechanism I offered to explain it does not.

## 5. Result: what the antibody reads is not rare

If the paratope's determinant is a proline–arginine–proline pattern at spacing
+0/+2/+4, then its usefulness as evidence is decided by how often that pattern
occurs. Census over all 20,431 reviewed human proteins (11,418,237 residues):

| pattern | occurrences | proteins with ≥1 | share of proteome |
|---|---|---|---|
| `P.R.P` (the shared read) | 3,719 | 2,841 | **13.9%** |
| `PP.R.P` (read + one flank) | 449 | 405 | 2.0% |
| `PPRRPP` (our published seed) | 6 | 6 | 0.03% |
| `RRPFF` (our other seed) | 3 | 3 | 0.01% |
| `GSPPRRPPPGRRPFFHP` (Hecker 17-mer) | 0 | 0 | 0% |
| `AAAAA` (size control) | 2,443 | 691 | 3.4% |
| `WWWWW` (rarity control) | 0 | 0 | 0% |

**What these numbers do and do not support.** After §4b, the only defensible
reading of this table is narrower than the one an earlier draft gave it. The
census itself is exact and unchanged: the patterns that make our two seeds look
specific in sequence (`PPRRPP` 6 proteins, `RRPFF` 3 proteins, Hecker 17-mer 0)
are far rarer than the proline/arginine signals the structures actually show in
contact. What is **not** supported is calling `P.R.P` "the pattern the antibody
reads": with ten contacting residues per peptide (§4b median 10 of 12), no
three-position abstraction describes the interface. `P.R.P` was my operational
summary, chosen because R396 is the single largest contributor; it is not a
mechanism.

The durable point survives in weaker form: **the residues that create apparent
specificity in a linear alignment are not the residues that carry the protein–
protein contact**, and this is what all four sequence axes tripped over when
tested against their own positives (dossier §8). Why exactly they trip is not
answered by this step — the mechanism I proposed for that, concentration on one
side chain, was refuted in §4b.

## 6. What this does and does not change for our 55 candidates

**Applies to very few.** This paratope is defined for one epitope, EBNA1
386-405. Of the 55 candidates only **MAPK8IP2** (`Q13387`, EBNA1 394-429) spans
that region; the top-ranked candidates match elsewhere in EBNA1 (PRKAR1B:
564-587; ZNF347: 566-586) and cannot be threaded onto this paratope at all. For
them a paratope test would need an antibody that does not exist.

**Does not change the interpretation of the ranking as much as it first
appeared.** The mechanism it pointed to — a real mimic rescued by one
well-placed side chain that no alignment score can express — is exactly the
claim §4b refutes. What remains is the weaker and still useful statement: a
linear score sums over residues the antibody weights very unevenly, so the
mapping from score to binding is not merely noisy, it is misaligned.

**Does not touch the candidates themselves.** No score, no rank and no filter in
the pipeline is modified by this step.

## 7. Limits, stated plainly

- **Contact count is not energy.** A large arginine makes more 4.5 Å contacts
  than an alanine by geometry alone. These numbers say where the interface is,
  not what it costs to break.
- **Two calibration complexes, but a 50-complex panel for the claim.** The
  register and the identity (n=2) are descriptive; the concentration hypothesis
  was tested against 50 independent complexes and failed there. Note the panel
  is itself selected — complexes that could be crystallised, in a database, with
  a SAbDAb annotation — so it is a reference distribution, not a random sample
  of all paratopes.
- **A single static conformer.** No ensemble, no induced fit — and the MD study
  on this exact antibody (Zhang & Yu 2026) reports induced fit and
  conformational selection as the binding mechanism, i.e. exactly what a static
  crystal cannot see.
- **No side-chain modelling.** I did not rebuild or repack side chains, so this
  is not a threading score; it is a contact census of the experimental
  structures.
- **The `P.R.P` census ignores context.** A `P.R.P` in a folded protein core is
  not an accessible epitope. The 13.9% figure is an upper bound on how often
  such a pattern occurs at all, not a set of candidates — and after §4b it is
  not claimed to be the read pattern either.

## 8. What would actually move this forward

1. **Affinity, not similarity.** Take the three validated peptides plus the top
   candidates and measure/compute ΔΔG against one antibody. Requires AF3 or
   Boltz-2 plus a scoring function, or SPR/ELISA.
2. **Test the one testable candidate.** `MAPK8IP2` 309-341 is the only candidate
   that overlaps this epitope; competing peptides for MS39p2w174 would be a
   single, cheap, decisive experiment — if the antibody is obtainable (IEDB
   returns 0 records for `Q13387`, so availability is unverified).
3. **Ask what the paratope requires in 3D, not in sequence.** The honest version
   of this channel is "does the candidate's side chain fit the R396 pocket",
   which needs side-chain modelling or co-folding, not a contact census.

---

*Every number above is reproduced by `env/venv/bin/python src/probe_paratope.py`
and `env/venv/bin/python src/probe_read_motif.py` on the two hashed PDB files.
The census in §5 uses the same `data/cache/human_swissprot.fasta` as the rest of
the pipeline. Nothing in this step alters a pipeline output.*
