"""Is "one residue carries the interface" a property of one crystal, or of paratopes?

probe_paratope.py read two crystals of a single Fab (MS39p2w174 + EBNA1 386-405
in 7K7R, + phospho-GlialCAM 370-389 in 7UQC). In both, ~half of the interfacial
contact mass sits on ONE arginine, flanked by prolines at -2 and +2, and that
three-position read (P.R.P) occurs in 13.9% of human proteins.

n = 2 cannot tell whether that is a general property of peptide-antibody
recognition or a peculiarity of this antibody. This script assembles an
independent panel and re-measures.

Panel: every experimental entry that contains (a) a short polymer entity
(5-25 aa) and (b) at least one entity annotated by SAbDab as an antibody.
Selection uses no knowledge of the answer -- it is a structural-class filter,
not a candidate filter.

For each entry: identify the antibody chains (>= 100 modelled residues), take
the peptide chain with the most antibody contacts, and measure
  - how concentrated the contact mass is (dominant residue share),
  - what three-position read at (-2, 0, +2) around the dominant residue is,
  - how often that read occurs in the human proteome.

Threshold choices are declared once, here, and applied to every entry.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))

from probe_paratope import parse_pdb, one  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "pdb"
REPORTS = ROOT / "reports"
LOGS = ROOT / "logs"

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# --- declared thresholds -----------------------------------------------------
CUTOFF = 4.5          # A, heavy-atom contact distance (same as probe_paratope)
PEPTIDE_MIN = 5       # modelled residues
PEPTIDE_MAX = 25      # modelled residues
FAB_MIN = 100         # modelled residues for a chain to count as antibody
MIN_CONTACTS = 20     # total peptide-antibody contacts for an entry to be usable
DOM_NEIGHBOUR = 2     # the read is taken at (-2, 0, +2) around the dominant

QUERY = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_sample_sequence_length",
                "operator": "range",
                "value": {"from": PEPTIDE_MIN, "to": PEPTIDE_MAX,
                          "include_lower": True, "include_upper": True}}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_feature_summary.type",
                "operator": "exact_match",
                "value": "SABDAB_ANTIBODY_NAME"}},
        ],
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {"start": 0, "rows": 10000},
        "results_content_type": ["experimental"],
    },
}

# Census patterns must not be able to match the '.' separator itself.
ATOMS = set("ACDEFGHIKLMNPQRSTVWY")


def search_entries():
    r = requests.post(SEARCH_URL, json=QUERY, timeout=120)
    r.raise_for_status()
    payload = r.json()
    return [hit["identifier"] for hit in payload.get("result_set", [])]


def fetch_pdb(pdb_id, force=False):
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / f"{pdb_id.lower()}.pdb"
    if out.exists() and not force:
        return out
    r = requests.get(DOWNLOAD_URL.format(pdb_id=pdb_id.upper()), timeout=120)
    if r.status_code != 200:
        return None
    out.write_text(r.text)
    return out


def load_proteome():
    sys.path.insert(0, str(ROOT / "src"))
    from phmimic.mining import load_fasta
    return load_fasta(ROOT / "data" / "cache" / "human_swissprot.fasta")


def census(pattern, records):
    """Occurrences of a 'A.B.C' pattern, '.' = any residue."""
    regex = re.compile("(?=(" + pattern + "))")
    occ = 0
    prot = 0
    for seq in records.values():
        n = len(regex.findall(seq))
        if n:
            occ += n
            prot += 1
    return occ, prot


def analyse_entry(pdb_id, path):
    chains = parse_pdb(path)
    if not chains:
        return None
    sizes = {c: len(res) for c, res in chains.items()}
    fab_chains = [c for c, n in sizes.items() if n >= FAB_MIN]
    if len(fab_chains) < 2:
        return None                     # no antibody pair
    peptides = [c for c, n in sizes.items()
                if PEPTIDE_MIN <= n <= PEPTIDE_MAX]
    if not peptides:
        return None
    fab_atoms = np.array([xyz for c in fab_chains
                          for r in chains[c].values()
                          for _n, _e, xyz in r["atoms"]])
    best = None
    for pc in peptides:
        peptide = chains[pc]
        keys = sorted(peptide)
        seq = "".join(one(peptide, r) for r in keys)
        prof = {}
        for r in keys:
            A = np.array([xyz for _n, _e, xyz in peptide[r]["atoms"]])
            d = np.linalg.norm(A[:, None, :] - fab_atoms[None, :, :], axis=2)
            prof[r] = int((d < CUTOFF).sum())
        total = sum(prof.values())
        if best is None or total > best["total"]:
            best = {"chain": pc, "keys": keys, "seq": seq,
                    "prof": prof, "total": total}
    if best is None or best["total"] < MIN_CONTACTS:
        return None
    seq, keys, prof = best["seq"], best["keys"], best["prof"]
    contact_pos = [i for i, r in enumerate(keys) if prof[r] > 0]
    if not contact_pos:
        return None
    dom_i = max(range(len(keys)), key=lambda i: prof[keys[i]])
    dom_contacts = prof[keys[dom_i]]
    read = None
    lo, hi = dom_i - DOM_NEIGHBOUR, dom_i + DOM_NEIGHBOUR
    if 0 <= lo and hi < len(seq):
        letters = [seq[lo], seq[dom_i], seq[hi]]
        if all(ch in ATOMS for ch in letters):
            read = f"{letters[0]}.{letters[1]}.{letters[2]}"
    penta = None
    if dom_i - 4 >= 0 and dom_i + 4 < len(seq):
        cand = seq[dom_i - 4:dom_i + 5]
        if all(ch in ATOMS for ch in cand):
            penta = cand
    return {
        "pdb": pdb_id.upper(),
        "peptide_chain": best["chain"],
        "peptide_len": len(keys),
        "peptide_seq": seq,
        "n_fab_chains": len(fab_chains),
        "total_contacts": best["total"],
        "n_contact_positions": len(contact_pos),
        "dominant_index": dom_i + 1,        # 1-based within the peptide
        "dominant_aa": seq[dom_i],
        "dominant_contacts": dom_contacts,
        "dominant_share": round(dom_contacts / best["total"], 4),
        "read_pattern": read,
        "pentapeptide": penta,
    }


def main():
    LOGS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    entries = search_entries()
    print(f"[panel] {len(entries)} candidate entries", flush=True)

    rows = []
    skipped = {"no_pdb": 0, "no_fab": 0, "no_peptide_contacts": 0}
    for i, pdb_id in enumerate(entries, 1):
        path = fetch_pdb(pdb_id)
        if path is None:
            skipped["no_pdb"] += 1
            continue
        row = analyse_entry(pdb_id, path)
        if row is None:
            skipped["no_peptide_contacts"] += 1
            continue
        rows.append(row)
        if i % 10 == 0:
            print(f"[panel] {i}/{len(entries)} -> {len(rows)} usable", flush=True)

    print(f"[panel] usable complexes: {len(rows)}", flush=True)
    records = load_proteome()
    print(f"[panel] proteome: {len(records)} proteins", flush=True)

    cache: dict[str, tuple[int, int]] = {}
    for row in rows:
        for field in ("read_pattern", "pentapeptide"):
            pat = row[field]
            if not pat:
                row[f"{field}_occ"] = None
                row[f"{field}_proteins"] = None
                continue
            if pat not in cache:
                cache[pat] = census(pat, records)
            occ, prot = cache[pat]
            row[f"{field}_occ"] = occ
            row[f"{field}_proteins"] = prot

    cols = ["pdb", "peptide_chain", "peptide_len", "peptide_seq",
            "n_fab_chains", "total_contacts", "n_contact_positions",
            "dominant_index", "dominant_aa", "dominant_contacts",
            "dominant_share", "read_pattern", "read_pattern_occ",
            "read_pattern_proteins", "pentapeptide", "pentapeptide_occ",
            "pentapeptide_proteins"]
    out = REPORTS / "panel_paratope.tsv"
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join("" if row.get(c) is None else str(row[c])
                               for c in cols) + "\n")

    def vals(key):
        return sorted(r[key] for r in rows if r.get(key) is not None)

    shares = vals("dominant_share")
    n_prot = len(records)
    read_rate = [r["read_pattern_proteins"] / n_prot for r in rows
                 if r.get("read_pattern_proteins") is not None]
    penta_rate = [r["pentapeptide_proteins"] / n_prot for r in rows
                  if r.get("pentapeptide_proteins") is not None]
    summary = {
        "n_entries_considered": len(entries),
        "n_usable_complexes": len(rows),
        "skipped": skipped,
        "parameters": {"cutoff_a": CUTOFF, "peptide_min": PEPTIDE_MIN,
                       "peptide_max": PEPTIDE_MAX, "fab_min": FAB_MIN,
                       "min_contacts": MIN_CONTACTS,
                       "dominant_neighbour": DOM_NEIGHBOUR},
        "dominant_share": {
            "min": shares[0] if shares else None,
            "median": shares[len(shares) // 2] if shares else None,
            "max": shares[-1] if shares else None,
            "n_ge_0_30": sum(1 for s in shares if s >= 0.30),
            "n_ge_0_40": sum(1 for s in shares if s >= 0.40),
            "n_total": len(shares),
        },
        "read_pattern_proteome_fraction": {
            "min": min(read_rate) if read_rate else None,
            "median": sorted(read_rate)[len(read_rate) // 2] if read_rate else None,
            "max": max(read_rate) if read_rate else None,
            "n_proteome": n_prot,
        },
        "pentapeptide_proteome_fraction": {
            "min": min(penta_rate) if penta_rate else None,
            "median": sorted(penta_rate)[len(penta_rate) // 2] if penta_rate else None,
            "max": max(penta_rate) if penta_rate else None,
        },
        "note": ("Panel is PDB peptide-antibody complexes: entries that can be "
                 "crystallised are a biased sample and this panel contains no "
                 "EBNA1 complex. It tests whether 'one residue carries the "
                 "interface' is typical of paratopes, not whether EBNA1 "
                 "mimicry is special."),
    }
    (REPORTS / "panel_paratope_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
