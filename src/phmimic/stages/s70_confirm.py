"""s70 -- confirmation layer: higher-resolution expression and prior literature.

This stage does NOT feed the score. It exists to answer, for the handful of
proteins that already scored well, two questions the bulk annotation cannot:
where in the brain is it expressed (HPA single-nucleus clusters rather than a
keyword), and has anyone already looked at it in an EBV/MS/autoimmune context.

Literature counts are reported and never scored: a protein with 27 papers is
not more likely to be a mimic than one with none -- it is more likely to have
been looked at.
"""
from __future__ import annotations

import time
import urllib.parse

import pandas as pd
import requests

from ..core import write_json

UA = {"User-Agent": "ph-mimic/1.0 (confirmation)"}


def _get_json(url: str, timeout: int = 45):
    for _ in range(3):
        try:
            resp = requests.get(url, headers=UA, timeout=timeout)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(1.5)
    return None


def _fetch_ensg(acc: str, cache):
    path = cache / "uniprot" / f"{acc}.json"
    if path.exists():
        payload = path.read_text(encoding="utf-8")
        return payload
    resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                        headers=UA, timeout=45)
    if resp.status_code != 200:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(resp.text, encoding="utf-8")
    return resp.text


def _ensg_of(acc: str, cache) -> str:
    import json

    text = _fetch_ensg(acc, cache)
    if not text:
        return ""
    data = json.loads(text)
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "Ensembl":
            continue
        for prop in xref.get("properties", []):
            value = str(prop.get("value", ""))
            if prop.get("key") == "GeneId" and value.startswith("ENSG"):
                # UniProt carries a versioned id (ENSG00000184144.13); the atlas
                # only resolves the unversioned form.
                return value.split(".")[0]
    return ""


def _hpa(ensg: str, cache):
    import json

    if not ensg:
        return None
    path = cache / "hpa" / f"{ensg}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    data = _get_json(f"https://www.proteinatlas.org/{ensg}.json")
    if data:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data), encoding="utf-8")
    return data


def _pubmed(gene: str, cache) -> int | None:
    import json

    if not gene:
        return None
    path = cache / "pubmed" / f"{gene}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8")).get("count")
    query = (f'{gene}[tiab] AND (autoantibod*[tiab] OR autoantigen*[tiab] '
             f'OR "multiple sclerosis"[tiab] OR "Epstein-Barr"[tiab] OR EBNA1[tiab])')
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
           f"&term={urllib.parse.quote(query)}&retmode=json&retmax=3")
    data = _get_json(url)
    count = None
    if data:
        count = int(data.get("esearchresult", {}).get("count", 0))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(dict(count=count, query=query)), encoding="utf-8")
        time.sleep(0.35)  # E-utilities allows 3 requests/second
    return count


def run(ctx):
    ranked_path = ctx.out("40_candidates_ranked.tsv")
    if not ranked_path.exists():
        ctx.log("    40_candidates_ranked.tsv missing; skipped")
        return dict(outputs=[], facts=dict(skipped="no ranking"))
    ranked = pd.read_csv(ranked_path, sep="\t", dtype={"acc": str})
    top_n = int(ctx.cfg.get("report.top_n", 20))
    finalists = (ranked.sort_values(["score", "ident"], ascending=[False, False])
                 .drop_duplicates("acc").head(top_n))

    cache = ctx.cfg.path_of("data.cache_dir")
    rows = []
    for row in finalists.itertuples(index=False):
        ensg = _ensg_of(row.acc, cache)
        hpa = _hpa(ensg, cache) or {}
        brain_tpm = ((hpa.get("RNA tissue specific nTPM") or {}) or {}).get("brain")
        try:
            brain_tpm = float(brain_tpm)
        except (TypeError, ValueError):
            brain_tpm = None
        rows.append(dict(
            seed_id=row.seed_id,
            acc=row.acc,
            gene=row.gene,
            score=row.score,
            ident=row.ident,
            matched=row.matched,
            tpos=row.pos,
            ensg=ensg,
            hpa_brain_ntpm=brain_tpm,
            hpa_brain_cluster=str(hpa.get("Brain expression cluster") or "")[:60],
            hpa_single_cell=str(hpa.get("Single cell expression cluster") or "")[:60],
            hpa_secretome=str(hpa.get("Secretome location") or "")[:40],
            hpa_subcellular=str(hpa.get("Subcellular main location") or "")[:40],
            pubmed_autoimmune=_pubmed(row.gene, cache),
        ))
        ctx.log(f"    {row.acc} {row.gene:<9} brain={brain_tpm} "
                f"cluster={str(hpa.get('Brain expression cluster'))[:34]}")

    frame = pd.DataFrame(rows)
    frame.to_csv(ctx.out("70_confirmation.tsv"), sep="\t", index=False)
    write_json(ctx.out("70_confirmation.json"), rows)
    return dict(
        outputs=["reports/70_confirmation.tsv", "reports/70_confirmation.json"],
        facts=dict(n_finalists=int(len(frame)),
                   with_brain_expression=int(frame["hpa_brain_ntpm"].notna().sum())),
    )
