"""s30 -- one bulk UniProt download supplies both orthogonal axes.

Downloaded once and cached; its sha256 goes into the run manifest so the
annotation snapshot behind every score is pinned.
"""
from __future__ import annotations

from ..annotate import annotation_provenance, ensure_annotation, annotation_sha256


def run(ctx):
    frame = ensure_annotation(ctx.cfg)
    n = len(frame)
    n_cns = int(frame["cns"].sum())
    n_surf = int(frame["surface"].sum())
    both = int((frame["cns"] & frame["surface"]).sum())
    ctx.log(f"    annotated {n:,} proteins: CNS {n_cns:,} ({n_cns / n:.1%}), "
            f"surface {n_surf:,} ({n_surf / n:.1%}), both {both:,}")
    return dict(
        outputs=[],
        facts=dict(n_proteins=n, n_cns=n_cns, n_surface=n_surf, n_both=both,
                   annotation=annotation_provenance(ctx.cfg)),
    )
