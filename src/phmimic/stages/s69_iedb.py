"""s69 -- is there any antibody-level evidence attached to these proteins?

IEDB holds curated B-cell assay records. Querying it is the only step here that
leaves our own modelling and touches other people's experiments, so the result
has to be read carefully in both directions:

  a hit  means someone assayed antibody binding to this protein -- the assay
         may be negative, may be in another disease, may be at another site;
  a miss means nobody recorded one -- which says nothing about biology and
         everything about which proteins people have bothered to study.

The reference proteins make that second point measurable: EBNA1 itself has
hundreds of records, the published mimics have dozens, and the candidates have
almost none. Any use of "IEDB coverage" as a filter would therefore be a filter
on research attention. It is recorded as an annotation only.

Records are cached on disk so a re-run costs nothing and the report can name
exactly what was asked and when.
"""
from __future__ import annotations

import time
from collections import Counter

import pandas as pd
import requests

from .. import mining
from ..core import rel, write_json


def _fetch(acc, api, limit, cache_dir, attempts=3, pause=2.0):
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"bcell_{acc.upper()}.json"
    if path.exists():
        try:
            import json
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            pass
    url = f"{api}/bcell_search"
    params = {"curated_source_antigen->>accession": f"like.{acc.upper()}*",
              "limit": str(limit)}
    last = None
    for attempt in range(attempts):
        try:
            resp = requests.get(url, params=params, timeout=120)
            resp.raise_for_status()
            data = resp.json()
            import json
            path.write_text(json.dumps(data), encoding="utf-8")
            return data
        except Exception as exc:                  # noqa: BLE001 - retried
            last = exc
            time.sleep(pause * (attempt + 1))
    raise RuntimeError(f"IEDB query failed for {acc}: {last}")


def _summarise(records, span):
    start, end = span
    n_coords = 0
    overlap, peptides, pmids, positive = 0, set(), set(), 0
    for rec in records:
        s, e = rec.get("starting_position"), rec.get("ending_position")
        if s is None or e is None:
            continue
        n_coords += 1
        try:
            s, e = int(s), int(e)
        except (TypeError, ValueError):
            continue
        if not mining.overlaps(s, e, start, end):
            continue
        overlap += 1
        if rec.get("linear_sequence"):
            peptides.add(str(rec["linear_sequence"]))
        for ref in rec.get("reference", []) or []:
            if ref.get("pubmed_id"):
                pmids.add(str(ref["pubmed_id"]))
        if str(rec.get("assay_group", "")).lower().find("positive") >= 0 or \
                str(rec.get(" qualitative_measurement", "")).strip():
            positive += 1
    return dict(
        n_records=len(records), n_with_coords=n_coords, n_overlap=overlap,
        n_overlap_positive=positive,
        overlap_peptides=";".join(sorted(peptides)[:8]),
        overlap_pmids=";".join(sorted(pmids)[:8]),
    )


def run(ctx):
    cfg = ctx.cfg
    api = str(cfg.get("mining.iedb.api", "https://query-api.iedb.org")).rstrip("/")
    limit = int(cfg.get("mining.iedb.limit", 5000))
    min_z = float(cfg.get("mining.iedb.min_margin_z", 6.0))
    cache_dir = cfg.path_of("mining.iedb.cache_dir")

    src = pd.read_csv(ctx.out("53_layered_table.tsv"), sep="\t", dtype=str)
    subset = src[src["margin_z"].astype(float) >= min_z]
    ctx.log(f"    querying {len(subset)} candidates with margin_z >= {min_z}")

    rows, failed = [], []
    for row in subset.to_dict("records"):
        acc = str(row["acc"]).upper()
        span = (int(row["t_start"]), int(row["t_end"]))
        try:
            records = _fetch(acc, api, limit, cache_dir)
            stat = _summarise(records, span)
            queried = "yes"
        except Exception as exc:                  # noqa: BLE001 - recorded
            ctx.log(f"    [{acc}] {type(exc).__name__}: {exc}")
            stat = dict(n_records="", n_with_coords="", n_overlap="",
                        n_overlap_positive="", overlap_peptides="",
                        overlap_pmids="")
            queried = "failed"
            failed.append(acc)
        rows.append(dict(row, iedb_queried=queried, **stat))
    out = pd.DataFrame(rows)

    refs = []
    for entry in cfg.get("mining.iedb.references", []) or []:
        acc = str(entry.get("acc", "")).upper() if isinstance(entry, dict) else str(entry).upper()
        gene = entry.get("gene", "") if isinstance(entry, dict) else ""
        try:
            records = _fetch(acc, api, limit, cache_dir)
            refs.append(dict(acc=acc, gene=gene, n_records=len(records)))
        except Exception as exc:                  # noqa: BLE001
            refs.append(dict(acc=acc, gene=gene, n_records="failed",
                             error=type(exc).__name__))
    out.to_csv(ctx.out("54_iedb_evidence.tsv"), sep="\t", index=False)

    zero = int((out["n_records"].astype(str) == "0").sum())
    payload = dict(
        api=api, limit=limit, min_margin_z=min_z, queried_utc=time.strftime(
            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        n_queried=int(len(subset)),
        n_zero_records=zero,
        n_failed=len(failed),
        failed=failed,
        n_with_any_overlap=int((pd.to_numeric(
            out["n_overlap"], errors="coerce").fillna(0) > 0).sum()),
        references=sorted(refs, key=lambda r: r["acc"]),
        queried_counts=dict(Counter(out["iedb_queried"])),
        note=("a hit is somebody else's assay, not our evidence; a miss reports "
              "research attention, not biology. Never used as a filter."),
    )
    write_json(ctx.out("69_iedb_summary.json"), payload)
    ctx.log(f"    {payload['n_queried']} queried, {zero} with zero records, "
            f"{payload['n_with_any_overlap']} with any overlapping record")
    return dict(
        outputs=["reports/54_iedb_evidence.tsv", "reports/69_iedb_summary.json"],
        facts=dict(n_queried=payload["n_queried"], n_zero_records=zero,
                   n_with_any_overlap=payload["n_with_any_overlap"]),
    )
