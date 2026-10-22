"""s00 -- prove the inputs before computing anything.

Nothing downstream is allowed to run on a file whose hash does not match the
manifest, and the UniProt release is checked as well as the bytes: a re-issued
proteome would change every count in the report even with identical filenames.
"""
from __future__ import annotations

from ..core import ProvenanceError, rel, tool_versions, verify_data, write_json


def run(ctx):
    try:
        verification = verify_data(ctx.cfg)
    except ProvenanceError as exc:
        raise RuntimeError(f"data provenance check failed: {exc}") from exc

    payload = dict(
        config_file=rel(ctx.cfg.path),
        config_sha256=ctx.cfg.sha256,
        inputs=verification["inputs"],
        problems=verification["problems"],
        tools=tool_versions(ctx.cfg),
    )
    write_json(ctx.out("00_provenance.json"), payload)
    ctx.log(f"    inputs verified: {len(payload['inputs'])}; tools recorded: {len(payload['tools'])}")
    return dict(
        outputs=["reports/00_provenance.json"],
        facts=dict(n_inputs=len(payload["inputs"]),
                   n_problems=len(payload["problems"]),
                   tools=payload["tools"]),
    )
