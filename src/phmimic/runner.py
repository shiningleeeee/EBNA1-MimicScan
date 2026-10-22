"""Stage runner: the only entry point.

Every stage declares what it produces, and the runner records the sha256 of each
output plus the config hash and the verified input hashes in
`reports/run_manifest.json`. A run therefore carries its own proof of what it
was: re-running the whole pipeline from the same config must reproduce both the
table contents and that manifest.

Usage:
    python -m phmimic.runner                 # run everything
    python -m phmimic.runner --from s40      # run from a stage onwards
    python -m phmimic.runner --only s40_priority
    python -m phmimic.runner --list
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
import traceback
from types import SimpleNamespace

from .core import ROOT, load_config, read_json, sha256_file, write_json

STAGES: list[tuple[str, str]] = [
    ("s00_provenance", "verify data hashes, record tool versions"),
    ("s10_seeds", "real seeds + two null seed families"),
    ("s20_motif_matches", "exact / near-exact motif scan of every seed"),
    ("s30_annotate", "bulk UniProt annotation (CNS + surface axes)"),
    ("s40_priority", "score every seed, calibrate against its null family"),
    ("s50_alignment", "sequence-alignment channel (seeded + exhaustive)"),
    ("s60_composition", "EBNA1 composition profile and Gly/Ala block"),
    ("s61_mmseqs", "unbiased local similarity: full EBNA1 vs the proteome"),
    ("s62_pair_null", "per-alignment composition null for every hit"),
    ("s63_selection_null", "is the number of hits itself above a shuffle?"),
    ("s64_rank", "eight declared axes + calibration of the whole list"),
    ("s65_structure", "structural accessibility of matched segments (AlphaFold)"),
    ("s66_maturation", "is the matched span present on the mature protein?"),
    ("s67_struct_pair", "folded, exposed surface on both sides of the pair"),
    ("s68_layered", "independent annotations, stacked and never fused"),
    ("s69_iedb", "antibody-level evidence, if any, with its references"),
    ("s70_confirm", "HPA + PubMed confirmation layer for finalists"),
    ("s71_candidates", "assemble the candidate table + sensitivity audit"),
    ("s72_axis_validation", "score every axis against the published positives"),
    ("s73_candidate_cards", "the candidate list itself, one card per pair"),
    ("s80_report", "generate REPORT.md and summary.json from the tables"),
]
STAGE_IDS = [sid for sid, _ in STAGES]


class Context:
    def __init__(self, cfg, quiet: bool = False):
        self.cfg = cfg
        self.quiet = quiet
        self.reports = cfg.reports
        self.logs = cfg.logs
        self.reports.mkdir(parents=True, exist_ok=True)
        self.logs.mkdir(parents=True, exist_ok=True)
        self.state: dict = {}

    def log(self, msg: str) -> None:
        if not self.quiet:
            print(msg, flush=True)

    def out(self, name: str):
        return self.reports / name


def run_stage(stage_id: str, ctx: Context) -> dict:
    module = importlib.import_module(f"phmimic.stages.{stage_id}")
    started = time.time()
    ctx.log(f"\n=== {stage_id} ===")
    result = module.run(ctx) or {}
    duration = round(time.time() - started, 2)
    outputs = {}
    for rel in result.get("outputs", []):
        path = ROOT / rel
        outputs[rel] = sha256_file(path) if path.exists() else "missing"
    ctx.log(f"    {duration}s  ->  {len(outputs)} output(s)")
    return dict(duration_s=duration, outputs=outputs, facts=result.get("facts", {}))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="ph-mimic pipeline")
    parser.add_argument("--config", default=None)
    parser.add_argument("--from", dest="start", default=None, choices=STAGE_IDS)
    parser.add_argument("--to", dest="stop", default=None, choices=STAGE_IDS)
    parser.add_argument("--only", default=None, choices=STAGE_IDS)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.list:
        for sid, description in STAGES:
            print(f"{sid:<18} {description}")
        return 0

    cfg = load_config(args.config)
    ctx = Context(cfg, quiet=args.quiet)

    if args.only:
        selected = [args.only]
    else:
        lo = STAGE_IDS.index(args.start) if args.start else 0
        hi = STAGE_IDS.index(args.stop) if args.stop else len(STAGE_IDS) - 1
        selected = STAGE_IDS[lo : hi + 1]

    manifest_path = ctx.out("run_manifest.json")
    manifest = read_json(manifest_path) if manifest_path.exists() else {}
    manifest.setdefault("project", cfg.get("project"))
    manifest.setdefault("version", cfg.get("version"))
    manifest["config"] = dict(
        file=str(cfg.path.relative_to(ROOT)), sha256=cfg.sha256
    )
    manifest.setdefault("stages", {})

    started_all = time.time()
    for stage_id in selected:
        try:
            manifest["stages"][stage_id] = run_stage(stage_id, ctx)
        except Exception as exc:  # noqa: BLE001 - a stage failure must stop the run
            manifest["stages"][stage_id] = dict(
                failed=type(exc).__name__, message=str(exc),
                traceback=traceback.format_exc().splitlines()[-6:],
            )
            write_json(manifest_path, manifest)
            print(f"\nSTAGE FAILED: {stage_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
            return 1
        write_json(manifest_path, manifest)

    manifest["total_duration_s"] = round(time.time() - started_all, 2)
    write_json(manifest_path, manifest)
    ctx.log(f"\nrun_manifest.json written ({manifest['total_duration_s']}s total)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
