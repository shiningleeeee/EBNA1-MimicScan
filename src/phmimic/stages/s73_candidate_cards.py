"""Candidate cards: the deliverable, as opposed to the statistics about it.

Everything upstream of this file answers "is the ordering significant". This one
answers the user's actual question: *which human proteins, and what do I do with
them*. It therefore deliberately does not sort by a p-value. It sorts by how much
is known about each candidate and how cheaply it could be tested, and it prints
the statistical weakness on every card rather than once in a footnote.

Declared ordering (fixed here, before looking at the output):

  1. shortlist = rows with aln_len >= 20 AND identity >= 40 AND the human protein
     is either CNS-expressed or antibody-accessible. Both thresholds are the ones
     already used for the table; neither was chosen here.
  2. within the shortlist, order by whether the composition null's tail was
     cleared at all (null_max < observed score over 10,000 both-sides shuffles),
     then by margin_z descending. margin_z alone is NOT selective (s72: keeps 92%
     of the pool at z >= 3, 41% at z >= 6), so it is the second key, not the
     first. Every card says so. A pair whose null tail is not cleared cannot be
     defended at any sample size, so it must not outrank one that can.
  3. report the longest run of identical residues for each pair, because that is
     what a synthetic peptide would be ordered from -- and note that it does not
     track validity: CRYAB's validated core is a 5-residue run, ANO2's is not.
  4. flag every pair whose EBNA1 window is itself low complexity (entropy < 3.0
     bits, or one residue above 30%), because for those the match is partly a
     property of our own sequence rather than of the human protein. The same
     filter also removes all three published linear mimics -- that is stated in
     the document rather than silently applied.

Outputs: reports/73_candidate_cards.tsv and reports/CANDIDATES.md.
"""
from __future__ import annotations

import textwrap

import pandas as pd

from ..core import rel, write_json
from ..mining import entropy
from . import s70_confirm as confirm

CNS_OR_SURFACE = ("cns_annotation", "surface_annotation")
SHORTLIST_MIN_LEN = 20
SHORTLIST_MIN_IDENTITY = 40.0
N_CARDS = 12
LOW_COMPLEXITY_ENTROPY = 3.0
LOW_COMPLEXITY_MAXFRAC = 0.30
CLASSES = {
    frozenset("AVLIM"): "hydrophobic",
    frozenset("FYW"): "aromatic",
    frozenset("STNQ"): "polar",
    frozenset("DE"): "acidic",
    frozenset("KRH"): "basic",
    frozenset("GP"): "turn",
    frozenset("C"): "cys",
}
CLASS_OF = {aa: name for group, name in CLASSES.items() for aa in group}


def _class(aa: str) -> str:
    return CLASS_OF.get(aa, "other")


def render(q: str, t: str) -> tuple[str, str, int]:
    """Return (marker string, similarity string, longest identical run).

    The marker is what gets printed under the alignment: ``|`` identical,
    ``+`` conservative (same residue class, different residue), space elsewhere.
    The second return value is kept for the table's ``sim`` column and is empty
    on purpose -- similarity is already carried by the marker, and printing a
    second offset line only makes the alignment harder to read.
    """
    marks, run, best = [], 0, 0
    for a, b in zip(q, t):
        if a == b:
            marks.append("|")
            run += 1
            best = max(best, run)
        else:
            marks.append("+" if a != "-" and b != "-" and _class(a) == _class(b) else " ")
            run = 0
    return "".join(marks), "", best


def _query_plddt(span, models) -> float | None:
    """Mean pLDDT over an EBNA1 span, from the cached AlphaFold model."""
    if not models:
        return None
    start, end = span
    values = [models.get(r) for r in range(start, end + 1)]
    values = [v for v in values if v is not None]
    return round(sum(values) / len(values), 1) if values else None


def run(ctx):
    table = ctx.out("55_candidates_final.tsv")
    if not table.exists():
        raise FileNotFoundError(f"{table} missing -- run s71_candidates first")
    rows = pd.read_csv(table, sep="\t", dtype=str)
    rows["margin_z_num"] = pd.to_numeric(rows["margin_z"])
    rows["identity_num"] = pd.to_numeric(rows["identity_pct"])
    rows["len_num"] = pd.to_numeric(rows["aln_len"])

    # EBNA1 disorder, per residue, from the cached model (one download at most).
    models = None
    try:
        from pathlib import Path

        from ..structure import af_pdb, parse_pdb

        cache = ctx.cfg.path_of("data.cache_dir") / "alphafold"
        pdb = af_pdb(ctx.cfg.get("structure.query_accession", "P03211"), cache)
        atom = parse_pdb(Path(pdb))
        models = dict(zip(atom["resnum"], atom["plddt"]))
    except Exception as exc:  # a missing model must not stop a text report
        ctx.log(f"  (EBNA1 model unavailable: {exc})")

    cache_dir = ctx.cfg.path_of("data.cache_dir")
    annot = pd.read_csv(ctx.cfg.path_of("data.uniprot_annot"), sep="\t",
                        dtype=str, keep_default_na=False).set_index("Entry")

    out_rows = []
    for row in rows.to_dict("records"):
        mark, sim, run_len = render(row["qaln"], row["taln"])
        plddt = _query_plddt((int(row["ebna1_start"]), int(row["ebna1_end"])), models)
        q_ungapped = row["qaln"].replace("-", "")
        q_h, q_maxfrac = entropy(q_ungapped)
        # Round first, then flag: the number printed on the card and the boolean
        # that drives the subset have to be the same number, or a row can display
        # "3.00 bits" and still be filed as low complexity.
        q_h = round(q_h, 3)
        q_maxfrac = round(q_maxfrac, 3)
        low_complex = q_h < LOW_COMPLEXITY_ENTROPY or q_maxfrac > LOW_COMPLEXITY_MAXFRAC
        accessible = str(row.get("compartment_accessible", "")).lower() == "yes"
        cns = str(row.get("cns_annotation", "")).lower() == "true"
        surface = str(row.get("surface_annotation", "")).lower() == "true"
        shortlist = (row["len_num"] >= SHORTLIST_MIN_LEN
                     and row["identity_num"] >= SHORTLIST_MIN_IDENTITY
                     and (cns or surface or accessible))
        out_rows.append(dict(
            gene=row["gene"], acc=row["acc"],
            ebna1_start=int(row["ebna1_start"]), ebna1_end=int(row["ebna1_end"]),
            target_start=int(row["target_start"]), target_end=int(row["target_end"]),
            aln_len=int(row["len_num"]), identity=row["identity_num"],
            margin_z=row["margin_z_num"], longest_run=run_len,
            evalue=row["evalue"], obs_score=row["obs_score"],
            null_mean=row["null_mean"], null_sd=row["null_sd"],
            null_max=row["null_max"], p_pair=row["p_pair"],
            maturation=row["maturation"],
            query_plddt=plddt,
            query_struct=row["query_struct"], target_struct=row["target_struct"],
            accessible=accessible, cns=cns, surface=surface,
            overlaps_published=bool(str(row.get("overlaps_published_query_region", ""))
                                    .strip() not in ("", "nan", "False", "false")),
            iedb_on_segment=row.get("iedb_on_segment", ""),
            null_cleared=float(row["null_max"]) < float(row["obs_score"]),
            query_entropy=q_h, query_maxfrac=q_maxfrac,
            query_low_complexity=low_complex,
            shortlist=shortlist,
            unknowns=str(row.get("unknowns", "")),
            manual_judgement=str(row.get("manual_judgement", "")),
            qaln=row["qaln"], taln=row["taln"], mark=mark, sim=sim,
        ))

    frames = pd.DataFrame(out_rows).sort_values(
        ["shortlist", "null_cleared", "margin_z"], ascending=[False, False, False])
    cards = frames.head(N_CARDS)
    n_not_cleared = int((~cards["null_cleared"]).sum())

    lines = []
    add = lines.append
    add("# EBNA1 mimicry candidates")
    add("")
    add(f"Fifty-five human proteins carry a composition-controlled local match to "
        f"EBNA1, out of 20,431 reviewed human proteins screened end to end. "
        f"This document is the list. The statistics that limit how far any ordering "
        f"can be pushed are in `METHOD.md`; they are repeated on every card below "
        f"because they apply to every card below.")
    add("")
    add("**What a row means.** A stretch of EBNA1 and a stretch of a human protein "
        "align better than chance composition would predict, after both stretches "
        "are shuffled with their own residues held fixed. That is all it means. It "
        "does not mean an antibody binds both, and nothing in this project measures "
        "binding.")
    add("")
    add("**Ordering.** Shortlist = length >= 20 aa, identity >= 40%, and the human "
        "protein is CNS-expressed, antibody-accessible, or both. Inside the "
        "shortlist the first key is whether the composition null's tail was cleared "
        "at all (`null_max < observed score` over 10,000 both-sides shuffles), and "
        "only then `margin_z`. `margin_z` on its own is a convenience, not evidence: "
        "`margin_z >= 6` retains 41% of the whole pool while recovering none of the "
        "four published mimics (`reports/72_axis_validation.tsv`).")
    add("")
    add("## The shortlist")
    add("")
    add("| # | gene | accession | EBNA1 | human | len | identity | longest run | margin z | null tail | low-complexity | accessible | CNS |")
    add("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, row in enumerate(cards.to_dict("records"), 1):
        add(f"| {i} | **{row['gene']}** | `{row['acc']}` | "
            f"{row['ebna1_start']}-{row['ebna1_end']} | "
            f"{row['target_start']}-{row['target_end']} | {row['aln_len']} | "
            f"{row['identity']:.1f}% | {row['longest_run']} | {row['margin_z']:.2f} | "
            f"{'cleared' if row['null_cleared'] else 'NOT cleared'} | "
            f"{'YES' if row['query_low_complexity'] else '-'} | "
            f"{'yes' if row['accessible'] else '-'} | {'yes' if row['cns'] else '-'} |")
    add("")
    add("`null tail` says whether the strongest of 10,000 same-composition shuffles "
        "scored **below** the observation. When it does not, the pair carries no "
        "information beyond its own residue composition — no p-value at any sample "
        f"size will change that, because the null is already at the observation. "
        f"**{n_not_cleared} of the {len(cards)} cards below are in that state**, and "
        f"they are ordered after the ones that are not: a pair whose null tail is "
        f"not cleared cannot be defended, so it should not outrank one that can.")
    add("")

    for i, row in enumerate(cards.to_dict("records"), 1):
        acc = row["acc"]
        entry = annot.loc[acc] if acc in annot.index else None
        name = (entry["Protein names"].split(" (")[0] if entry is not None else "")
        genes = (entry["Gene Names"].replace(" ", ", ") if entry is not None else "")
        location = (entry["Subcellular location [CC]"][:300] if entry is not None else "")
        tissue = (entry["Tissue specificity"][:220] if entry is not None else "")
        keywords = (entry["Keywords"][:220] if entry is not None else "")

        ensg = confirm._ensg_of(acc, cache_dir)
        hpa = confirm._hpa(ensg, cache_dir) or {} if ensg else {}
        ntpms = hpa.get("RNA tissue specific nTPM") or {}
        brain = ntpms.get("brain", "")
        cluster = hpa.get("Brain expression cluster", "")
        single = hpa.get("Single cell expression cluster", "")
        secretome = hpa.get("Secretome location", "")
        pubmed = confirm._pubmed(row["gene"], cache_dir)

        add(f"### {i}. {row['gene']} — `{acc}`")
        add("")
        if name:
            add(f"*{name}*" + (f" (genes: {genes})" if genes else ""))
            add("")
        add("```")
        add(f"EBNA1  {row['ebna1_start']:>4} {row['qaln']}")
        add(f"            {row['mark']}   | = identical, + = same residue class")
        add(f"{row['gene']:<7}{row['target_start']:>4} {row['taln']}")
        add("```")
        add("")
        ratio = float(row["null_max"]) / float(row["obs_score"])
        if ratio >= 1.0:
            tail = (f"the strongest of 10,000 shuffles scored {row['null_max']}, "
                    f"**at or above the observation** — the null tail is not "
                    f"cleared at all")
        else:
            tail = (f"the strongest of 10,000 shuffles scored {row['null_max']}, "
                    f"i.e. {100 * ratio:.0f}% of the observation")
        add(f"- **{row['aln_len']} aa, {row['identity']:.1f}% identical, "
            f"longest identical run {row['longest_run']}.** Observed BLOSUM62 score "
            f"{row['obs_score']} against a both-sides composition null of "
            f"{row['null_mean']} ± {row['null_sd']}; {tail} "
            f"(margin z = **{row['margin_z']:.2f}**, p = {row['p_pair']}).")
        add(f"- **EBNA1 side:** {row['query_struct'].replace('_', ' ')}, mean pLDDT "
            f"{row['query_plddt'] if row['query_plddt'] is not None else 'n/a'}; "
            f"window entropy {row['query_entropy']} bits, most frequent residue "
            f"{100 * row['query_maxfrac']:.1f}%"
            + (" — **low complexity: this match is partly a composition effect**"
               if row["query_low_complexity"] else "")
            + (". **This segment overlaps a validated epitope region.**"
               if row["overlaps_published"] else "."))
        add(f"- **Human side:** {row['target_struct'].replace('_', ' ')}, maturation "
            f"`{row['maturation']}`"
            + (", antibody-accessible compartment" if row["accessible"] else "")
            + (", CNS expression annotation" if row["cns"] else "") + ".")
        if location:
            add(f"- **Location:** {location}")
        if tissue:
            add(f"- **Tissue:** {tissue}")
        if keywords:
            add(f"- **Keywords:** {keywords}")
        hpa_bits = []
        if brain not in ("", None):
            hpa_bits.append(f"brain nTPM {brain}")
        if cluster:
            hpa_bits.append(f"`{cluster}`")
        if single:
            hpa_bits.append(f"single-cell `{single}`")
        if secretome:
            hpa_bits.append(f"secretome: {secretome}")
        if hpa_bits:
            add(f"- **HPA:** " + "; ".join(hpa_bits) + ".")
        add(f"- **Literature:** {pubmed if pubmed is not None else 'n/a'} PubMed "
            f"records for this gene together with autoantibody / autoantigen / "
            f"multiple sclerosis / EBV terms; {row['iedb_on_segment']} IEDB B-cell "
            f"records overlapping the matched segment.")
        add(f"- **Unknowns:** no antibody has been tested against these two peptides; "
            f"no PTM state modelled; no expression data in MS lesions or CSF; "
            f"disorder calls come from one AlphaFold conformer.")
        add(f"- **The experiment:** order both peptides, run a competition ELISA for "
            f"the human peptide against whatever EBNA1 reactivity is available, or "
            f"immunise and test cross-reactivity. A negative result is informative "
            f"because the pair is a pure composition case: the null already says a "
            f"shuffle of these same residues scores {row['null_mean']}.")
        add("")

    add("## The full list of 55")
    add("")
    add("| gene | accession | EBNA1 | human | len | identity | longest run | margin z | mature | accessible | CNS | low-complexity |")
    add("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for row in frames.to_dict("records"):
        add(f"| {row['gene']} | `{row['acc']}` | {row['ebna1_start']}-{row['ebna1_end']} | "
            f"{row['target_start']}-{row['target_end']} | {row['aln_len']} | "
            f"{row['identity']:.1f}% | {row['longest_run']} | {row['margin_z']:.2f} | "
            f"{'yes' if row['maturation'].startswith('retained') else 'NO'} | "
            f"{'yes' if row['accessible'] else '-'} | {'yes' if row['cns'] else '-'} | "
            f"{'YES' if row['query_low_complexity'] else '-'} |")
    add("")
    add("## The normal-composition subset")
    add("")
    normal = frames[~frames["query_low_complexity"]].sort_values("margin_z",
                                                                 ascending=False)
    add(f"Removing every pair whose EBNA1 window is itself low complexity "
        f"(entropy < {LOW_COMPLEXITY_ENTROPY:.1f} bits or one residue above "
        f"{LOW_COMPLEXITY_MAXFRAC:.0%}) leaves **{len(normal)} of {len(frames)}**. "
        f"This is the subset where the match is not explained by repeated residues "
        f"on our side of the pair.")
    add("")
    add("| gene | accession | EBNA1 | human | len | identity | margin z | accessible | CNS |")
    add("|---|---|---|---|---|---|---|---|---|")
    for row in normal.head(25).to_dict("records"):
        add(f"| {row['gene']} | `{row['acc']}` | {row['ebna1_start']}-{row['ebna1_end']} | "
            f"{row['target_start']}-{row['target_end']} | {row['aln_len']} | "
            f"{row['identity']:.1f}% | {row['margin_z']:.2f} | "
            f"{'yes' if row['accessible'] else '-'} | {'yes' if row['cns'] else '-'} |")
    add("")
    add("**A tension to state rather than hide.** All three published linear mimics "
        "are low-complexity pairs on our side as well (`PPRRPP`, `RRPFF`, "
        "`PGAIEQGPAD`), so this filter removes them too. That is not an argument "
        "against the filter: it is the reason the published cases could not be "
        "recovered by any composition-controlled statistic. But it means the two "
        "lists answer different questions — the full 55 asks *what does EBNA1 "
        "resemble*, the normal-composition subset asks *where does the resemblance "
        "survive removing repeated residues*. Pick the one that matches the "
        "experiment you intend to run.")
    add("")
    add("## What would make this list smaller and better")
    add("")
    add("Nothing in the pipeline currently measures binding, so no ordering here can "
        "be improved by more sequence analysis. The three things that would change "
        "the list are: a peptide-binding assay for any of these pairs; MS/CSF "
        "reactivity data against the human peptides; and an antibody against one of "
        "the human proteins, which would let the EBNA1 peptide be tested as a "
        "competitor. The only pair in this project with solved structures on both "
        "sides of a proven cross-reaction is EBNA1 386-405 / GlialCAM 370-389 "
        "(7K7R, 7UQC), and that pair scores `margin_z` 0.865 under the same null.")
    add("")

    path = ctx.out("CANDIDATES.md")
    path.write_text("\n".join(lines), encoding="utf-8")
    tsv = ctx.out("73_candidate_cards.tsv")
    frames.drop(columns=["qaln", "taln", "mark", "sim"]).to_csv(tsv, sep="\t", index=False)

    shortlist = frames[frames["shortlist"]]
    normal = frames[~frames["query_low_complexity"]]
    facts = dict(
        n_candidates=int(len(frames)),
        n_shortlist=int(len(shortlist)),
        n_normal_composition=int(len(normal)),
        n_query_low_complexity=int(frames["query_low_complexity"].sum()),
        n_cards=min(N_CARDS, len(frames)),
        top_gene=str(cards.iloc[0]["gene"]),
        top_margin_z=float(cards.iloc[0]["margin_z"]),
        top_longest_run=int(cards.iloc[0]["longest_run"]),
        n_null_tail_cleared=int(frames["null_cleared"].sum()),
        ordering="shortlist first, then null_cleared, then margin_z descending",
        ordering_caveat="margin_z is not selective (s72 keeps 0.92 of the pool at "
                        "z>=3 and 0.42 at z>=6, recovering 0 of 4 published mimics); "
                        "null_cleared is a hard ceiling, not a p-value",
        low_complexity_rule=f"EBNA1 window entropy < {LOW_COMPLEXITY_ENTROPY} bits "
                            f"or one residue above {LOW_COMPLEXITY_MAXFRAC:.0%}",
        shortlist_rule=f"aln_len >= {SHORTLIST_MIN_LEN} and identity >= "
                       f"{SHORTLIST_MIN_IDENTITY} and (CNS or surface or accessible)",
    )
    summary = ctx.out("73_candidate_cards_summary.json")
    write_json(summary, facts)
    ctx.log(f"  {len(frames)} candidates; {len(shortlist)} shortlisted; "
            f"{len(cards)} cards; top {facts['top_gene']} z={facts['top_margin_z']:.2f}")
    return dict(outputs=[rel(path), rel(tsv), rel(summary)], facts=facts)
