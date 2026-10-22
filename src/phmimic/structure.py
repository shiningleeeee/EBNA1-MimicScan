"""Structural accessibility of a linear motif, from AlphaFold models.

The question this answers is narrow and mechanical: a residue that is buried in
a folded core cannot be bound by an antibody, no matter how identical its
sequence is to a viral epitope. So for each candidate we ask what the predicted
model says about the matched segment -- how confident the model is there
(pLDDT, read from the B-factor column) and how much of that segment is
solvent-exposed (rSASA, computed here from the coordinates).

Both numbers are proxies, not measurements. pLDDT is a confidence, not an
order/disorder call, and a predicted monomer is not the cellular context. They
are reported as such and are never folded into the priority score.
"""
from __future__ import annotations

import json
import math
import pathlib
import time
import urllib.request

import numpy as np

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# van der Waals radii (A). Hydrogens are absent from AlphaFold models.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90,
       "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "ZN": 1.39, "MG": 1.73,
       "CA": 1.95, "NA": 2.27, "K": 2.75}

# Theoretical maximum solvent-accessible surface per residue, Tien et al. 2013,
# Ala-X-Ala tripeptide, in A^2. Used only to turn absolute SASA into a fraction.
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
           "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 228,
           "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


def _fetch_json(url: str, attempts: int = 3, pause: float = 2.0) -> dict:
    """GET a JSON document, retrying transient failures.

    Without the retry a dropped connection silently removes a row from the null
    sample, which changes the p-value between two otherwise identical runs.
    """
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001 - network errors vary by platform
            last = exc
            if attempt + 1 < attempts:
                time.sleep(pause * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def af_metadata(acc: str, cache_dir: pathlib.Path) -> dict:
    """AlphaFold DB prediction metadata for one accession (cached on disk)."""
    cache_dir = pathlib.Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{acc}.json"
    if path.exists():
        return json.loads(path.read_text())
    payload = _fetch_json(AF_API.format(acc=acc))
    entry = payload[0] if isinstance(payload, list) else payload
    path.write_text(json.dumps(entry))
    return entry


def af_pdb(acc: str, cache_dir: pathlib.Path) -> pathlib.Path:
    """Download and cache the AlphaFold PDB for one accession."""
    cache_dir = pathlib.Path(cache_dir)
    pdb = cache_dir / f"{acc}.pdb"
    if pdb.exists() and pdb.stat().st_size > 0:
        return pdb
    url = af_metadata(acc, cache_dir)["pdbUrl"]
    last: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                body = resp.read()
            if not body:
                raise RuntimeError("empty body")
            pdb.write_bytes(body)
            return pdb
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt < 2:
                time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def parse_pdb(path: pathlib.Path) -> dict:
    """Atoms of the first model: coords, element, residue number, residue name, pLDDT."""
    coords, elements, resnums, resnames, bfac = [], [], [], [], []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip().upper()
        name = line[12:16].strip()
        if not element:
            element = name[0] if name[0].isalpha() else name[1]
            element = element.upper()
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elements.append(element)
        resnums.append(int(line[22:26]))
        resnames.append(_THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
        bfac.append(float(line[60:66]))
    return dict(
        xyz=np.asarray(coords, dtype=np.float64),
        element=np.asarray(elements),
        resnum=np.asarray(resnums, dtype=np.int32),
        resname=np.asarray(resnames),
        plddt=np.asarray(bfac, dtype=np.float64),
    )


_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def _sphere_points(n: int) -> np.ndarray:
    """Roughly uniform points on the unit sphere (golden spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = math.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


def sasa_per_atom(atom: dict, probe: float = 1.4, n_points: int = 92) -> np.ndarray:
    """Shrake-Rupley solvent-accessible surface per atom, in A^2."""
    from scipy.spatial import cKDTree

    xyz = atom["xyz"]
    radii = np.array([VDW.get(e, 1.8) + probe for e in atom["element"]], dtype=np.float64)
    if len(xyz) == 0:
        return np.zeros(0)
    tree = cKDTree(xyz)
    points = _sphere_points(n_points)
    areas = np.zeros(len(xyz))
    sphere_area = 4.0 * math.pi
    for i in range(len(xyz)):
        neighbours = tree.query_ball_point(xyz[i], radii[i] + radii.max())
        if not neighbours:
            areas[i] = sphere_area * radii[i] ** 2
            continue
        nb = np.asarray([j for j in neighbours if j != i], dtype=np.int64)
        if nb.size == 0:
            areas[i] = sphere_area * radii[i] ** 2
            continue
        test = xyz[i] + points * radii[i]
        d = np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2)
        buried = (d < radii[nb][None, :]).any(axis=1)
        areas[i] = sphere_area * radii[i] ** 2 * (1.0 - buried.mean())
    return areas


def residue_stats(atom: dict, sasa: np.ndarray) -> dict:
    """Per-residue pLDDT and relative SASA."""
    if len(atom["resnum"]) == 0:
        return {}
    order = np.argsort(atom["resnum"], kind="stable")
    resnum = atom["resnum"][order]
    sasa_o = sasa[order]
    plddt_o = atom["plddt"][order]
    resname_o = atom["resname"][order]
    uniq, starts = np.unique(resnum, return_index=True)
    counts = np.diff(np.append(starts, len(resnum)))
    sasa_sum = np.add.reduceat(sasa_o, starts)
    plddt_sum = np.add.reduceat(plddt_o, starts)
    out = {}
    for i, res in enumerate(uniq.tolist()):
        name = str(resname_o[starts[i]])
        cap = MAX_ASA.get(name, 200)
        out[int(res)] = dict(resname=name,
                             sasa=float(sasa_sum[i]),
                             rsasa=float(sasa_sum[i]) / cap,
                             plddt=float(plddt_sum[i] / counts[i]))
    return out


def segment(atom: dict, stats: dict, start: int, end: int) -> dict:
    """Summarise a 1-based inclusive segment of the model."""
    picked = [stats[p] for p in range(start, end + 1) if p in stats]
    if not picked:
        return dict(n=0)
    rsasa = np.array([p["rsasa"] for p in picked])
    plddt = np.array([p["plddt"] for p in picked])
    return dict(
        n=len(picked),
        rsasa_mean=float(rsasa.mean()),
        rsasa_max=float(rsasa.max()),
        frac_rsasa_ge_025=float((rsasa >= 0.25).mean()),
        plddt_mean=float(plddt.mean()),
        plddt_min=float(plddt.min()),
        frac_plddt_ge_70=float((plddt >= 70).mean()),
    )


def classify(summary: dict) -> str:
    """Coarse accessibility call, stated in words rather than a score.

    A low-pLDDT segment has no confidenct fold, so its relative SASA says little
    about whether an antibody can reach it -- those are reported as disordered,
    with the exposed ones flagged separately.
    """
    if summary.get("n", 0) == 0:
        return "not_modelled"
    exposed = summary["rsasa_mean"] >= 0.25
    folded = summary["plddt_mean"] >= 70
    if not folded:
        return "disordered_or_loop_exposed" if exposed else "disordered_low_confidence"
    if exposed:
        return "surface_of_folded_domain"
    if summary["rsasa_mean"] >= 0.10:
        return "partly_buried"
    return "buried"
