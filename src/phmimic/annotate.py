"""Orthogonal annotation: is a matched protein CNS-expressed and reachable?

Both axes are read from ONE bulk UniProt download rather than one HTTP call per
protein. That matters for reproducibility as much as for speed: a per-protein
loop would make the score of every candidate depend on which requests happened
to succeed that evening, and the null distribution would then be built from a
different annotation snapshot than the real list.
"""
from __future__ import annotations

import time
import urllib.parse

import pandas as pd
import requests

from .core import ROOT, sha256_file

UNIPROT_FIELDS = ",".join(
    [
        "accession",
        "gene_names",
        "protein_name",
        "length",
        "cc_subcellular_location",
        "cc_tissue_specificity",
        "keyword",
    ]
)
UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=" + urllib.parse.quote("organism_id:9606 AND reviewed:true")
    + "&format=tsv&fields=" + UNIPROT_FIELDS
)
USER_AGENT = {"User-Agent": "ph-mimic/1.0 (annotation)"}

CNS_TOKENS = (
    "brain", "neuron", "neuronal", "synap", "axon", "glial", "myelin",
    "oligodendro", "astrocyt", "microglia", "schwann", "cerebell",
    "cerebral", "spinal", "central nervous", "retina", "nerve",
)
SURFACE_TOKENS = ("cell membrane", "secreted", "cell surface", "extracellular", "gpi-anchor")


def annotation_path(cfg):
    return cfg.path_of("data.uniprot_annot")


def ensure_annotation(cfg, force: bool = False) -> pd.DataFrame:
    """Download once, then reuse. Returns a frame indexed by accession."""
    path = annotation_path(cfg)
    if force or not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        last = None
        for attempt in range(4):
            try:
                resp = requests.get(UNIPROT_URL, headers=USER_AGENT, timeout=300)
                if resp.status_code == 200 and resp.text.startswith("Entry"):
                    path.write_bytes(resp.content)
                    break
                last = f"HTTP {resp.status_code}, {len(resp.content)} bytes"
            except requests.RequestException as exc:  # pragma: no cover
                last = str(exc)
            time.sleep(3 * (attempt + 1))
        else:
            raise RuntimeError(f"could not download UniProt annotation: {last}")
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    frame = frame.rename(
        columns={
            "Entry": "acc",
            "Gene Names": "genes",
            "Protein names": "name",
            "Length": "length",
            "Subcellular location [CC]": "subcell",
            "Tissue specificity": "tissue",
            "Keywords": "keywords",
        }
    )
    frame["cns"] = frame.apply(lambda r: _has(r, CNS_TOKENS), axis=1)
    frame["surface"] = frame.apply(lambda r: _has(r, SURFACE_TOKENS), axis=1)
    return frame.set_index("acc", drop=False)


def _has(row, tokens) -> bool:
    blob = " ".join(
        str(row.get(col, "")) for col in ("subcell", "tissue", "keywords", "name")
    ).lower()
    return any(tok in blob for tok in tokens)


def annotation_sha256(cfg) -> str:
    path = annotation_path(cfg)
    return sha256_file(path) if path.exists() else "absent"


def annotation_provenance(cfg) -> dict:
    return dict(
        file=str(annotation_path(cfg).relative_to(ROOT)),
        url=UNIPROT_URL,
        fields=UNIPROT_FIELDS,
        sha256=annotation_sha256(cfg),
    )
