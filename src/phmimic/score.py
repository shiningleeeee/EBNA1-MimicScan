"""Scoring a seed by the best candidate it can reach.

The score deliberately mixes one sequence term and two biology terms, and the
biology terms are the ones that actually move the ranking. A perfect 6/6 match
to a protein nobody can put an antibody on, that is not expressed where the
disease happens, is not a candidate -- it is a string match.

Literature support is NOT scored: whether a protein has been studied is
confounded with how much attention the field paid it, so it is reported next to
the ranking and can only ever downgrade a call, never create one.
"""
from __future__ import annotations

import pandas as pd


def score_seed(matches: list[tuple[str, int, int, str]], annotation, cfg) -> dict:
    """Reduce one seed's raw matches to its best achievable candidate."""
    exact_pts = int(cfg.get("priority.exact_points", 2))
    partial_pts = int(cfg.get("priority.partial_points", 1))
    partial_frac = float(cfg.get("priority.partial_frac", 0.70))
    cns_pts = int(cfg.get("priority.cns_points", 1))
    surf_pts = int(cfg.get("priority.surface_points", 1))

    if not matches:
        return dict(n_matches=0, n_proteins=0, best_score=None, best_acc="",
                    best_gene="", best_motif="", best_pos="", best_ident="")

    best_by_acc: dict[str, tuple[int, int, str]] = {}
    for acc, pos, n_ident, text in matches:
        prev = best_by_acc.get(acc)
        if prev is None or n_ident > prev[0]:
            best_by_acc[acc] = (n_ident, pos, text)

    rows = []
    length = len(matches[0][3])
    for acc, (n_ident, pos, text) in best_by_acc.items():
        frac = n_ident / length
        if n_ident == length:
            seq_pts = exact_pts
        elif frac >= partial_frac:
            seq_pts = partial_pts
        else:
            continue
        cns = bool(annotation.at[acc, "cns"]) if acc in annotation.index else False
        surf = bool(annotation.at[acc, "surface"]) if acc in annotation.index else False
        rows.append(
            dict(
                acc=acc,
                gene=annotation.at[acc, "genes"].split()[0] if acc in annotation.index
                and annotation.at[acc, "genes"] else "",
                name=(annotation.at[acc, "name"][:60] if acc in annotation.index else ""),
                pos=pos,
                n_ident=n_ident,
                motif_len=length,
                ident=round(frac, 3),
                seq_pts=seq_pts,
                cns=int(cns),
                surface=int(surf),
                score=seq_pts + cns_pts * int(cns) + surf_pts * int(surf),
                matched=text,
            )
        )
    if not rows:
        return dict(n_matches=len(matches), n_proteins=len(best_by_acc), best_score=None,
                    best_acc="", best_gene="", best_motif="", best_pos="", best_ident="")
    frame = pd.DataFrame(rows).sort_values(
        ["score", "ident", "acc"], ascending=[False, False, True]
    )
    top = frame.iloc[0]
    return dict(
        n_matches=len(matches),
        n_proteins=len(best_by_acc),
        best_score=int(top["score"]),
        best_acc=top["acc"],
        best_gene=top["gene"],
        best_motif=top["matched"],
        best_pos=int(top["pos"]),
        best_ident=float(top["ident"]),
        candidates=frame,
    )


def empirical_p(observed: int | None, null_scores: list[int | None]) -> dict:
    """One-sided empirical p-value against a null seed family.

    `null_scores` may contain None for null seeds that match nothing at all; a
    seed that matches nothing cannot out-score anything, so it counts against
    the observation. Both the unconditional and the conditional (restricted to
    null seeds that matched something) values are returned: the first treats the
    null seed space as the reference, the second treats "seeds that occur in the
    human proteome at all" as the reference. They answer different questions and
    the difference between them is itself informative.
    """
    values = [(-1 if s is None else int(s)) for s in null_scores]
    n_all = len(values)
    present = [v for v in values if v >= 0]
    obs = -1 if observed is None else int(observed)
    ge_all = sum(1 for v in values if v >= obs)
    ge_present = sum(1 for v in present if v >= obs)
    return dict(
        n_null=n_all,
        n_null_present=len(present),
        n_null_ge=ge_all,
        p_unconditional=round((1 + ge_all) / (1 + n_all), 4) if n_all else None,
        p_conditional=round((1 + ge_present) / (1 + len(present)), 4) if present else None,
        null_max=max(values) if values else None,
        null_median_present=float(pd.Series(present).median()) if present else None,
    )
