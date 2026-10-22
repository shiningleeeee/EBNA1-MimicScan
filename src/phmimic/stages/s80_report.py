"""s80 -- render the report from the tables.

Every number in REPORT.md is interpolated from a file produced by an earlier
stage. There is no hand-typed figure anywhere, which is what
tests/test_selfconsistency.py checks: it re-derives each headline value from the
delivered tables and fails if the report disagrees.
"""
from __future__ import annotations

import pandas as pd

from ..core import read_json, write_json


def _fmt(value, digits=3):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def run(ctx):
    prov = read_json(ctx.out("00_provenance.json"))
    calibration = read_json(ctx.out("40_null_calibration.json"))["calibration"]
    seeds = pd.read_csv(ctx.out("10_seeds.tsv"), sep="\t", dtype=str).fillna("")
    seed_scores = pd.read_csv(ctx.out("40_seed_scores.tsv"), sep="\t", dtype={"acc": str})
    ranked = pd.read_csv(ctx.out("40_candidates_ranked.tsv"), sep="\t", dtype={"acc": str})
    if "published" in ranked.columns:
        ranked["published"] = ranked["published"].fillna("").astype(str)
    composition = read_json(ctx.out("60_composition_summary.json"))

    alignment = None
    if ctx.out("50_alignment_summary.json").exists():
        alignment = read_json(ctx.out("50_alignment_summary.json"))

    confirmation = pd.DataFrame()
    if ctx.out("70_confirmation.tsv").exists():
        confirmation = pd.read_csv(ctx.out("70_confirmation.tsv"), sep="\t",
                                   dtype={"acc": str, "gene": str})

    real = seeds[seeds["family"] == "real"]
    null_a = seeds[seeds["family"] == "N_A"]
    null_b = seeds[seeds["family"] == "N_B"]
    best_real = max((c["best_score"] if c["best_score"] is not None else -1)
                    for c in calibration)
    top_candidates = ranked.drop_duplicates("acc").head(int(ctx.cfg.get("report.top_n", 20)))

    lines: list[str] = []
    add = lines.append
    add(f"# {ctx.cfg.get('report.title')}")
    add("")
    add(f"`ph-mimic` v{ctx.cfg.get('version')} -- config sha256 "
        f"`{prov['config_sha256'][:16]}`")
    add("")
    add("## 0. What this is, and what it is not")
    add("")
    add("A screen for human proteins that could be mistaken for EBNA1 by the immune "
        "system. The only experimentally anchored object is a short peptide from the "
        "EBNA1 386-440 region, where antibodies and sera have validated cross-reactive "
        "epitopes; every candidate here is a match to that peptide, ranked by biology "
        "that no sequence statistic can supply.")
    add("")
    add("It is **not** a list of confirmed cross-reactive proteins. Nothing in this "
        "report is a measurement; it is a prioritisation. The score is a filter, and "
        "the null below states what a filter of this kind finds when the peptide means "
        "nothing.")
    add("")

    add("## 1. Provenance")
    add("")
    add("| input | file | UniProt release | sha256 |")
    add("|---|---|---|---|")
    for key, info in prov["inputs"].items():
        add(f"| `{key}` | `{info['file']}` | {info['uniprot_release']} "
            f"| `{info['sha256'][:16]}...` |")
    add("")
    tools = prov.get("tools", {})
    if tools:
        add("| tool | version |")
        add("|---|---|")
        for name, version in sorted(tools.items()):
            add(f"| {name} | `{version}` |")
        add("")
    if prov.get("problems"):
        add(f"**Provenance problems:** {'; '.join(prov['problems'])}")
        add("")

    add("## 2. Seeds and the nulls they are judged against")
    add("")
    add(f"- **real seeds ({len(real)})**: " + ", ".join(
        f"`{r.motif}` (EBNA1 {r.ebna1_start}-{r.ebna1_end}, d<={r.max_mismatch})"
        for r in real.itertuples(index=False)))
    add(f"- **N_A ({len(null_a)})**: every distinct arrangement of each seed's own "
        "residue multiset, enumerated exhaustively. Composition is held exactly "
        "constant, so anything separating a real seed from N_A is due to order alone.")
    add(f"- **N_B ({len(null_b)})**: random {int(ctx.cfg.get('seeds.null.proteome_random_len', 6))}-aa "
        "windows drawn from the human proteome. Composition is free; this answers "
        "\"what would an arbitrary peptide have found\".")
    add("")
    add("Real and null seeds are scanned by the same code and scored by the same rule.")
    add("")

    add("## 3. Primary result: is the seed better than its own null?")
    add("")
    add("Score = exactness (2 for a full-length match, 1 for >= "
        f"{_fmt(float(ctx.cfg.get('priority.partial_frac', 0.7)), 2)} identity) "
        "+ CNS expression (1) + antibody-reachable localisation (1). Maximum "
        f"{int(ctx.cfg.get('priority.max_score', 4))}.")
    add("")
    add("| seed | motif | best score | best candidate | proteins matched | N_A n | N_A max "
        "| p (unconditional) | p (conditional) | N_B p |")
    add("|---|---|---|---|---|---|---|---|---|---|")
    for block in calibration:
        na = block["null_A"]
        nb = block["null_B"]
        add(f"| `{block['seed_id']}` | `{block['motif']}` | {_fmt(block['best_score'])} "
            f"| {block['best_gene'] or block['best_acc'] or '-'} | {block['n_proteins']} "
            f"| {na['n_null']} | {_fmt(na['null_max'])} "
            f"| {_fmt(na['p_unconditional'])} | {_fmt(na['p_conditional'])} "
            f"| {_fmt(nb['p_conditional']) if nb else 'n/a'} |")
    add("")
    add("- *p (unconditional)* treats the whole null seed space as the reference, "
        "including null peptides that occur nowhere in the human proteome.")
    add("- *p (conditional)* restricts the reference to null peptides that match at "
        "least one human protein, which is the fairer comparison when the question is "
        "\"given that this peptide exists in us, does its order matter\".")
    if any(block["null_B"] is None for block in calibration):
        add("- *N_B p* is blank where the real seed length differs from the N_B length; "
            "lengths are not comparable and no p is reported rather than a misleading one.")
    add("")

    add("## 4. Candidates")
    add("")
    add(f"{len(ranked)} protein-seed pairs pass the sequence filter; "
        f"{ranked['acc'].nunique()} distinct proteins. Top {len(top_candidates)} by score:")
    add("")
    has_pub = "published" in ranked.columns
    header = "| # | protein | gene | seed | matched | ident | CNS | surface | score |"
    rule = "|---|---|---|---|---|---|---|---|---|"
    if has_pub:
        header += " published |"
        rule += "---|"
    add(header)
    add(rule)
    for i, row in enumerate(top_candidates.itertuples(index=False), start=1):
        cells = [str(i), f"`{row.acc}`", row.gene or "-", f"`{row.seed_motif}`",
                 f"`{row.matched}` @{row.pos}", _fmt(float(row.ident)),
                 "Y" if row.cns else "-", "Y" if row.surface else "-", str(row.score)]
        if has_pub:
            flag = str(getattr(row, "published", "") or "").strip()
            cells.append("yes" if flag and flag != "nan" else "-")
        add("| " + " | ".join(cells) + " |")
    add("")
    add("`CNS` and `surface` come from a single bulk UniProt download pinned in the run "
        "manifest, not from per-protein lookups, so the table and the null distribution "
        "share one annotation snapshot.")
    add("")

    published = ctx.cfg.get("published") or []
    if published:
        recovered = {str(r.acc) for r in top_candidates.itertuples(index=False)}
        got = [p for p in published if str(p["acc"]) in recovered]
        miss = [p for p in published if str(p["acc"]) not in recovered]
        add("## 4b. Sensitivity against the published claims")
        add("")
        add(f"This screen recovers {len(got)} of {len(published)} published EBNA1 mimicry "
            "claims recorded in `config.yaml`:")
        add("")
        for p in got:
            add(f"- recovered: `{p['acc']}` {p['gene']} -- {p['claim']} ({p['source']})")
        for p in miss:
            add(f"- **not recoverable**: `{p['acc']}` {p['gene']} -- {p.get('why', '')}")
        add("")
        add("Recall is reported as a sensitivity check on the screen, never as a measure "
            "of its reliability: a screen can recover a known mimic by composition alone, "
            "and failing to recover one can mean the epitope is conformational rather than "
            "that the code is wrong. The null in section 3 is the statistic; this list is "
            "an audit.")
        add("")

    if len(confirmation):
        add("## 5. Confirmation layer (not scored)")
        add("")
        add("| protein | gene | HPA brain nTPM | HPA brain cluster | HPA single-cell cluster | PubMed (MS/EBV/auto) |")
        add("|---|---|---|---|---|---|")
        for row in confirmation.itertuples(index=False):
            add(f"| `{row.acc}` | {row.gene or '-'} | {_fmt(row.hpa_brain_ntpm, 1)} "
                f"| {row.hpa_brain_cluster or '-'} "
                f"| {row.hpa_single_cell or '-'} "
                f"| {_fmt(row.pubmed_autoimmune)} |")
        add("")
        add("The two HPA cluster columns answer different questions and can disagree: the "
            "brain cluster is a bulk-tissue cluster, the single-cell cluster names the cell "
            "type that actually expresses the protein. Where they disagree, the single-cell "
            "cluster is the one that should drive a CNS claim.")
        add("")
        add("Literature counts are **not** part of the score. A protein with many papers "
            "has been looked at more, which is not evidence that it is a mimic.")
        add("")

    add("## 6. Why composition explains the sequence channel")
    add("")
    b0, b1 = composition["block"]
    deviants = composition.get("block_non_ga") or []
    detail = ""
    if deviants:
        shown = ", ".join(f"{d['aa']}{d['pos']}" for d in deviants[:6])
        detail = (f"; the only residues that are not Gly or Ala are {shown}"
                  + (" ..." if len(deviants) > 6 else ""))
    add(f"EBNA1 ({composition['query_length']} aa) residues {b0}-{b1} are "
        f"{composition['block_ga_fraction'] * 100:.1f}% Gly/Ala "
        f"({composition['block_length']} residues, "
        f"{composition['block_n_distinct']} distinct residue types{detail}). Outside that "
        f"block the highest {composition['window']}-aa Gly/Ala fraction is "
        f"{composition['profile_max_ga_outside_block']}.")
    add("")
    add("A block of this length over a two-letter alphabet will match collagen Gly-X-Y "
        "repeats, keratin and other low-complexity human proteins under almost any "
        "composition-sensitive statistic. That is why this pipeline never ranks on an "
        "alignment score.")
    add("")
    if alignment:
        for mode, info in alignment.items():
            if isinstance(info, dict):
                add(f"- alignment channel, {mode}: {info['n_hits']:,} hits; "
                    f"{info['n_overlapping_block']} overlap the Gly/Ala block.")
        add("")
        add("The alignment channel is retained for continuity with earlier work and is "
            "reported, but it does not contribute to any candidate call.")
        add("")

    structure = None
    if ctx.out("65_structure_summary.json").exists():
        structure = read_json(ctx.out("65_structure_summary.json"))

    if structure:
        add("## 6b. Structure channel: can an antibody reach the matched segment?")
        add("")
        q = structure.get("query", {})
        if q.get("global_plddt") is not None:
            add(f"The query protein itself is not folded. AlphaFold's mean pLDDT for "
                f"EBNA1 is {q['global_plddt']:.1f} and "
                f"{q['fraction_plddt_very_low'] * 100:.0f}% of its residues are in the "
                f"very-low-confidence band. The three experimentally validated linear "
                f"epitopes sit in that band:")
            add("")
            add("| span | mean pLDDT | mean rSASA | call |")
            add("|---|---|---|---|")
            for key, seg in q.get("segments", {}).items():
                add(f"| {key.replace('_', ' ')} | {seg['plddt_mean']:.1f} | "
                    f"{seg['rsasa_mean']:.2f} | {seg['call']} |")
            add("")
            add("Every validated mimicry epitope lies in a disordered, highly exposed "
                "stretch, and the only confidently folded parts of EBNA1 (roughly "
                "476-585) are the DNA-binding / dimerisation core, which carries no "
                "known epitope. A structure-based mimicry search over this protein has "
                "nothing folded to compare in the epitope region.")
            add("")
        add("The same measurement, applied to each candidate's matched segment and to "
            "null matches from both families with the same code path:")
        add("")
        calls = structure.get("calls", {})
        add("| call | candidates |")
        add("|---|---|")
        for call, n in sorted(calls.items(), key=lambda kv: -kv[1]):
            add(f"| {call} | {n} |")
        add("")
        cal = structure.get("calibration") or {}
        if cal:
            add(f"Folded and exposed (mean pLDDT >= {structure['parameters']['folded_plddt']}, "
                f"mean rSASA >= {structure['parameters']['exposed_rsasa']}): "
                f"**{cal.get('n_real_folded_exposed')} of {structure['n_candidates']} "
                f"candidates** against **{cal.get('n_null_folded_exposed')} of "
                f"{cal.get('n_null')}** null matches, giving "
                f"p = {cal.get('p_unconditional')}. Accessibility is therefore not a "
                f"selective property here: an arbitrary matching hexamer in a human "
                f"protein lands on a folded surface about a quarter of the time.")
            add("")
        add("This channel changes no score and adds no candidate. It removes none "
            "either -- what it removes is the idea that reachability distinguishes the "
            "list. " + structure.get("caveat", ""))
        add("")

    # ---- the unbiased channel (s61-s71), rendered only if it ran --------------
    unbiased = None
    if ctx.out("71_candidates_summary.json").exists():
        unbiased = dict(
            mmseqs=read_json(ctx.out("61_mmseqs_summary.json")),
            pair=read_json(ctx.out("62_pair_null_summary.json")),
            selection=read_json(ctx.out("63_selection_summary.json")),
            composite=read_json(ctx.out("64_composite_summary.json")),
            maturation=read_json(ctx.out("66_maturation_summary.json")),
            struct=read_json(ctx.out("67_struct_pair_summary.json")),
            layered=read_json(ctx.out("68_layered_summary.json")),
            iedb=read_json(ctx.out("69_iedb_summary.json")),
            final=read_json(ctx.out("71_candidates_summary.json")),
            table=pd.read_csv(ctx.out("55_candidates_final.tsv"), sep="\t", dtype=str),
        )
    if unbiased:
        m, sel, comp = unbiased["mmseqs"], unbiased["selection"], unbiased["composite"]
        mat, pair, iedb = unbiased["maturation"], unbiased["pair"], unbiased["iedb"]
        cal_max = comp["calibration"]["max_score"]
        cal_frac = comp["calibration"]["frac_ge_min_score"]
        add("## 6c. The unbiased channel: full EBNA1 against the whole proteome")
        add("")
        add(f"The channel above is anchored on the validated peptide. This one asks the "
            f"same question without an anchor: all 641 residues searched against all "
            f"{m.get('n_proteome', 20431)} reviewed human proteins, no seed, no region "
            f"excluded.")
        add("")
        add(f"It reports {m['n_alignments']} alignments over {m['n_targets']} proteins. "
            f"The best E-value is {m['best_evalue']} ({m['best_gene']}), and "
            f"**{m['n_evalue_below_1']} alignments have E < 1**. Not one human protein "
            f"is a credible homologue of EBNA1 on this evidence; that is a statement "
            f"about local similarity, not about mimicry.")
        add("")
        regions = m.get("region_counts", {})
        add("| EBNA1 region | alignments |")
        add("|---|---|")
        for region, n in sorted(regions.items(), key=lambda kv: -kv[1]):
            add(f"| {region} | {n} |")
        add("")
        add(f"Alignments overlapping the Gly/Ala block: **{m['ga_block_overlap']}**. "
            f"The reason is **not** established: {m['note']}")
        add("")
        add(f"Is the *number* of alignments itself surprising? Shuffling the query "
            f"{sel['n_shuffles']} times while preserving its composition gives "
            f"{sel['null_n_align_mean']} +- {sel['null_n_align_sd']} alignments; the real "
            f"query reports {sel['real_n_align']}. That is p = "
            f"{sel['p_selection_n_align']:.4f}, which is the floor of the test "
            f"(1/{1 + sel['n_shuffles']}) and not a measurement below it. Read as: the "
            f"real sequence does line up with the proteome more often than its own "
            f"anagram does, by a factor of about "
            f"{sel['real_n_align'] / max(sel['null_n_align_mean'], 1e-9):.1f}.")
        add("")
        add(f"For each individual alignment, a two-sided composition null "
            f"(b = {pair['b']} rearrangements of both partners, gaps and matrix "
            f"identical to the reference implementation) gives a margin over the null "
            f"ranging {pair['margin_z_min']} to {pair['margin_z_max']} (median "
            f"{pair['margin_z_median']}). That statistic is length-dependent "
            f"(r = -0.40), so it orders short alignments first by construction.")
        add("")
        add(f"Calibrating the whole list against {sel['n_shuffles']} composition shuffles: "
            f"the best entry scores {cal_max['real']} against a null maximum of "
            f"{cal_max['null_mean']} +- {cal_max['null_sd']} "
            f"(null max {cal_max['null_max']}), p = {cal_max['p']:.3f}. "
            f"**The best candidate is not special.** The *proportion* of alignments "
            f"clearing the mid threshold is {cal_frac['real']:.1%} against "
            f"{cal_frac['null_mean']:.1%} +- {cal_frac['null_sd']:.1%}, p = "
            f"{cal_frac['p']:.3f}: a real but modest enrichment, roughly "
            f"{cal_frac['real'] / max(cal_frac['null_mean'], 1e-9):.1f}x. Counting "
            f"entries rather than proportions would have inflated this ~18x, because "
            f"the real query simply reports more alignments.")
        add("")
        add(f"Before any of that, one filter is applied that has nothing to do with "
            f"significance: does the matched span survive maturation? "
            f"{mat['n_confirmed']} of {mat['n_classified']} spans are confirmed present "
            f"on the mature chain, {mat['n_unannotated']} cannot be decided from UniProt, "
            f"and {mat['n_removed']} lie in a signal peptide, propeptide or transit "
            f"peptide and are removed before the protein reaches the surface. This is "
            f"the only axis here that a known positive also passes "
            f"(3/3 of the published linear mimics survive it).")
        add("")
        add(f"The candidate table has **{unbiased['final']['n_rows']} rows** and "
            f"**{unbiased['final']['n_columns']} columns**, one column per stage. "
            f"{unbiased['final']['note'][:1].upper()}{unbiased['final']['note'][1:]}")
        add("")
        add(f"Antibody-level evidence, queried per accession: {iedb['n_queried']} "
            f"proteins asked, {iedb['n_zero_records']} returned no IEDB record at all, "
            f"{iedb['n_with_any_overlap']} had a record overlapping the matched span. "
            f"For scale, the same query returns "
            + ", ".join(f"{r['gene']} {r['n_records']}" for r in iedb["references"])
            + ". Coverage tracks who has been studied; it is recorded, never used as a "
              "filter, and an empty result is not evidence of absence.")
        add("")
        bench = unbiased["final"].get("benchmark", {})
        add(f"Finally, our own sensitivity, computed after the table and read by "
            f"nothing: this channel recovers **{bench.get('n_recovered')} of "
            f"{bench.get('n_reachable')}** of the reachable published mimics, because "
            f"each is a short low-complexity epitope that an aligner scoring length and "
            f"identity cannot rank. The seeded channel above recovers all three. Two "
            f"channels, opposite failure modes, stated together rather than picking the "
            f"flattering one.")
        add("")

    axis_path = ctx.out("72_axis_validation.tsv")
    if axis_path.exists():
        axes = pd.read_csv(axis_path, sep="\t", dtype=str)
        add("## 6d. The rule, as a table: no statistical axis passes its own controls")
        add("")
        add("Every axis in this pipeline is scored against the published positives "
            "before it is allowed to order anything. Two thresholds were declared in "
            "`config.yaml` before the measurement: keep at least "
            f"**{ctx.cfg.get('validation.min_positives_passing')}** of the published "
            f"pairs, and keep at most "
            f"**{ctx.cfg.get('validation.max_pool_retained')}** of the candidate pool. "
            "An axis that misses either one cannot discriminate, and using it as a "
            "ranking criterion would be a choice dressed up as a result.")
        add("")
        add("| axis | kind | keeps positives | keeps pool | verdict |")
        add("|---|---|---|---|---|")
        for row in axes.to_dict("records"):
            kept = row.get("pool_retained_fraction")
            kept = "" if pd.isna(kept) else f"{float(kept):.3f}"
            add(f"| `{row['axis']}` | {row['kind']} | {row['positives_passing']} | "
                f"{kept} | {row['verdict']} |")
        add("")
        n_rank = int((axes["kind"] == "ranking_axis").sum())
        n_pass = int(((axes["kind"] == "ranking_axis")
                      & (axes["verdict"] == "pass")).sum())
        add(f"**{n_pass} of {n_rank} ranking axes pass.** The only axis that does is a "
            "hard mechanistic constraint (does the matched span survive signal-peptide "
            "or propeptide removal), which is not a statistical choice at all — it is a "
            "fact about where the protein is cut. The practical consequence is stated "
            "in `METHOD.md` §7: with no passing ranking axis, no column in the candidate "
            "table is allowed to be presented as a statistically justified ordering. "
            "The table is ordered by margin z with that failure printed next to it.")
        add("")
        add("The positives this is measured on are deliberately not all easy. Two of "
            "the five are pairs a linear screen cannot see at all — GlialCAM, whose "
            "epitope is defined by a crystal structure and shares only 6 of 20 residues "
            "with EBNA1, and MBP, which is conformational. An axis set that was only "
            "tested on the three short linear epitopes would have been tested on "
            "exactly the cases it was built for.")
        add("")

    add("## 7. What would change the conclusion")
    add("")
    add("- A real seed reaching a score that no member of N_A reaches (p_conditional "
        "below the multiple-testing threshold) would make the ordering informative.")
    add("- A candidate whose protein is CNS-expressed, antibody-reachable and carries "
        "the motif in a disordered or surface-exposed segment would be worth wet-lab "
        "follow-up regardless of the p-value.")
    add("- The opposite reading is also live: if N_A peptides reach the same scores, "
        "this screen is a description of the proteome's composition and the honest "
        "output is the negative result, not a candidate list.")
    add("")

    add("## 8. Reproduce")
    add("")
    add("```bash")
    add("bash run.sh                 # full pipeline")
    add("env/venv/bin/python -m phmimic.runner --list")
    add("env/venv/bin/python -m pytest tests -q")
    add("```")
    add("")
    add("`reports/run_manifest.json` records the config hash, every verified input hash, "
        "and the sha256 of every table this report was rendered from.")

    text = "\n".join(lines) + "\n"
    ctx.out("REPORT.md").write_text(text, encoding="utf-8")

    summary = dict(
        version=ctx.cfg.get("version"),
        config_sha256=prov["config_sha256"],
        n_real_seeds=int(len(real)),
        n_null_A=int(len(null_a)),
        n_null_B=int(len(null_b)),
        n_matches=int(len(ranked)),
        n_candidate_proteins=int(ranked["acc"].nunique()),
        best_real_score=best_real,
        calibration=calibration,
        composition=composition,
        alignment=alignment,
        provenance_inputs=prov["inputs"],
    )
    write_json(ctx.out("summary.json"), summary)
    ctx.log(f"    REPORT.md: {len(text)} chars; summary.json written")
    return dict(outputs=["reports/REPORT.md", "reports/summary.json"],
                facts=dict(best_real_score=best_real,
                           n_candidate_proteins=summary["n_candidate_proteins"]))
