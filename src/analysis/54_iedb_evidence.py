#!/usr/bin/env python3
"""External-evidence layer: what has anyone actually measured on these proteins?

The question this answers is narrow and deliberately so: for each candidate
protein, does the Immune Epitope Database hold B-cell (antibody) assay records,
and if so, do any of them cover the segment our alignment points at?

Three outcomes mean different things:

    no records at all          nobody has assayed this protein with antibodies
                               -- absence of evidence, not evidence of absence
    records, but not on the    the protein is assayed, our segment is not
    matched segment
    records on the segment     a measured antibody response exists there; the
                               assay outcome (positive/negative) is quoted

Inputs
    reports/53_layered_table.tsv   candidate list with matched segments
Output
    reports/54_iedb_evidence.tsv
    data/cache/iedb/bcell_<ACC>.json
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data/cache/iedb"
API = "https://query-api.iedb.org/bcell_search"

# Reference proteins: what "good" looks like, so a zero is interpretable
REFERENCES = {
    "P03211": "EBNA1 (the query itself)",
    "P02511": "CRYAB (validated cross-reactive target)",
    "Q9NQ90": "ANO2 (validated)",
    "Q14CZ8": "GlialCAM/HEPACAM (validated)",
    "Q02080": "MEF2B (PPRRPP cross-reactivity)",
}


def fetch(acc, attempts=3):
    """All B-cell records for one accession, cached on disk."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"bcell_{acc}.json"
    if path.exists():
        return json.loads(path.read_text())
    params = {
        "curated_source_antigen->>accession": f"like.{acc}*",
        "limit": "5000",
    }
    url = f"{API}?{urllib.parse.urlencode(params)}"
    last = None
    for i in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                data = json.load(fh)
            path.write_text(json.dumps(data))
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2.0 * (i + 1))
    print(f"  ! {acc}: {last}", flush=True)
    return []


def summarise(records, start, end):
    """Overlap between IEDB epitope coordinates and the matched segment."""
    coords = [
        r
        for r in records
        if isinstance(r.get("curated_source_antigen"), dict)
        and r["curated_source_antigen"].get("starting_position") is not None
        and r["curated_source_antigen"].get("ending_position") is not None
    ]
    hits = []
    for r in coords:
        s = r["curated_source_antigen"]["starting_position"]
        e = r["curated_source_antigen"]["ending_position"]
        if s <= end and start <= e:
            hits.append(r)
    positive = [
        r for r in hits if (r.get("qualitative_measure") or "").lower() == "positive"
    ]
    pmids = sorted({str(r.get("pubmed_id")) for r in hits if r.get("pubmed_id")})
    titles = sorted({t for r in hits for t in (r.get("reference_titles") or [])})
    return {
        "n_records": len(records),
        "n_with_coords": len(coords),
        "n_overlap": len(hits),
        "n_overlap_positive": len(positive),
        "overlap_peptides": "|".join(
            sorted({r.get("linear_sequence") or "?" for r in hits})[:6]
        ),
        "overlap_pmids": ",".join(pmids[:8]),
        "overlap_title": (titles[0][:110] if titles else ""),
    }


def main():
    with open(ROOT / "reports/53_layered_table.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cands = [r for r in rows if float(r["margin_z"]) >= 6.0]
    cands.sort(key=lambda r: -float(r["margin_z"]))

    targets = [(r["acc"], r["gene"], int(r["t_start"]), int(r["t_end"])) for r in cands]
    extra = [(acc, gene, 0, 0) for acc, gene in REFERENCES.items()]

    out = []
    for acc, gene, start, end in targets + extra:
        print(f"  querying {gene} ({acc}) ...", flush=True)
        records = fetch(acc)
        info = summarise(records, start, end)
        info.update(acc=acc, gene=gene, seg=f"{start}-{end}" if start else "-")
        out.append(info)
        time.sleep(0.5)

    dst = ROOT / "reports/54_iedb_evidence.tsv"
    with open(dst, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out)

    print(f"\nwrote {dst.name}\n")
    print(f"{'gene':11s}{'acc':12s}{'seg':>12s}{'recs':>7s}{'coords':>8s}{'on_seg':>8s}{'pos':>5s}")
    for r in out:
        seg = r["seg"]
        print(
            f"{r['gene']:11s}{r['acc']:12s}{seg:>12s}"
            f"{r['n_records']:>7d}{r['n_with_coords']:>8d}{r['n_overlap']:>8d}{r['n_overlap_positive']:>5d}"
        )

    print("\n=== proteins with a measured antibody response on the matched segment ===")
    any_hit = False
    for r in out:
        if r["n_overlap"]:
            any_hit = True
            print(f"  {r['gene']} {r['acc']} seg={r['seg']}")
            print(f"     peptides: {r['overlap_peptides']}")
            print(f"     positive: {r['n_overlap_positive']}/{r['n_overlap']}  PMIDs: {r['overlap_pmids']}")
            print(f"     ref: {r['overlap_title']}")
    if not any_hit:
        print("  none")

    zero = [r["gene"] for r in out if r["n_records"] == 0]
    print(f"\nproteins with no IEDB B-cell records at all ({len(zero)}): {', '.join(zero)}")


if __name__ == "__main__":
    sys.exit(main())
