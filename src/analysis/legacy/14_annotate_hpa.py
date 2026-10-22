#!/usr/bin/env python3
"""14_annotate_hpa.py -- orthogonal annotation for the Step-4 shortlist.

Two axes that no sequence statistic can supply:
  (a) is the gene expressed in the CNS?   -> Human Protein Atlas RNA (brain nTPM,
      brain single-nucleus cluster, brain expression cluster)
  (b) is the protein reachable by an antibody from outside the cell?
      -> HPA secretome location + subcellular main location
  (c) has the protein already been reported as an autoantigen / in MS / in EBV
      work? -> NCBI E-utilities count on a fixed query

Usage: env/venv/bin/python src/14_annotate_hpa.py
Writes: reports/candidates_annotated.tsv
"""
import json
import sys
import time
import urllib.parse

import pandas as pd
import requests

UA = {"User-Agent": "ph-mimic/0.1 (annotation)"}
SL = pd.read_csv("reports/candidates_shortlist.tsv", sep="\t")


def get(url, tries=3, timeout=40):
    for _ in range(tries):
        try:
            r = requests.get(url, headers=UA, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(1.5)
    return None


def ensembl_of(acc):
    """Return the Ensembl GENE id (ENSG...). UniProt's Ensembl xref id itself is
    a transcript (ENST...), which 404s on proteinatlas.org -- the gene id lives
    in the cross-reference properties under key 'GeneId'."""
    u = get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    if not u:
        return None
    for x in u.get("uniProtKBCrossReferences", []):
        if x.get("database") != "Ensembl":
            continue
        for p in x.get("properties", []):
            if p.get("key") == "GeneId" and str(p.get("value", "")).startswith("ENSG"):
                # UniProt carries a versioned id (ENSG00000184144.13);
                # proteinatlas.org only resolves the unversioned form.
                return p["value"].split(".")[0]
    return None


def hpa(ensg):
    return get(f"https://www.proteinatlas.org/{ensg}.json")


def pubmed_hits(gene):
    q = (f'{gene}[tiab] AND (autoantibod*[tiab] OR autoantigen*[tiab] '
         f'OR "multiple sclerosis"[tiab] OR "Epstein-Barr"[tiab] OR EBNA1[tiab])')
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
           f"&term={urllib.parse.quote(q)}&retmode=json&retmax=3")
    j = get(url)
    if not j:
        return None, []
    res = j.get("esearchresult", {})
    return int(res.get("count", 0)), res.get("idlist", [])


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


rows = []
for _, r in SL.iterrows():
    acc, gene = r["acc"], (r["gene"] if isinstance(r["gene"], str) else "")
    ensg = ensembl_of(acc)
    h = hpa(ensg) if ensg else None
    brain_tpm = brain_cluster = sc_cluster = None
    secret = mainloc = None
    tissue_spec = None
    if h:
        bt = h.get("RNA tissue specific nTPM") or {}
        brain_tpm = num(bt.get("brain")) if isinstance(bt, dict) else None
        brain_cluster = h.get("Brain expression cluster")
        sc_cluster = h.get("Single cell expression cluster")
        secret = h.get("Secretome location")
        mainloc = h.get("Subcellular main location")
        tissue_spec = h.get("RNA tissue specificity")
    cnt, pms = pubmed_hits(gene) if gene else (None, [])
    cns = bool((brain_tpm and brain_tpm >= 20) or
               (brain_cluster and "brain" in str(brain_cluster).lower()) or
               (isinstance(sc_cluster, str) and any(k in sc_cluster for k in
                    ("Astrocyte", "Oligodendrocyte", "Neuron", "Microglia",
                     "Excitatory", "Inhibitory", "Glial"))))
    rows.append(dict(
        acc=acc, gene=gene, name=r["name"][:40], cores=r["cores"],
        n_exact=r["n_exact"], best_frac=r["best_frac"],
        brain_nTPM=brain_tpm, brain_cluster=(brain_cluster or "")[:38],
        brain_cell=(sc_cluster or "")[:34],
        secretome=(secret or "")[:26], subcell=(mainloc or "")[:30],
        tissue_spec=tissue_spec or "",
        cns_data="Y" if cns else "",
        pm_auto=cnt, pm_top=",".join(pms or []),
    ))
    print(f"{acc:<8} {gene:<9} brain_nTPM={str(brain_tpm):>7} | {str(sc_cluster)[:30]:<30} "
          f"| secret={str(secret)[:22]:<22} | PM_auto={cnt}")
    time.sleep(0.3)

out = pd.DataFrame(rows)
out.to_csv("reports/candidates_annotated.tsv", sep="\t", index=False)
print(f"\n[out] reports/candidates_annotated.tsv ({len(out)} rows)")
print(out[["gene", "n_exact", "best_frac", "brain_nTPM", "cns_data", "pm_auto"]]
      .to_string(index=False))
