"""A3, the half that does not need an antibody in a tube.

ROADMAP A3 was "7K7R Fv x MAPK8IP2 peptide, precondition: verify reagent
availability". The precondition fails: MS39p2w174 is not a catalogue antibody, and
the only route to it is recombinant expression from the Fab sequence deposited in
7K7R / 7UQC. So the wet half cannot start here. The computational half can, and it
is worth doing because it asks a question the candidate table never asked:

    does a candidate keep the residues this antibody actually touches?

Every ranking axis in this project compared whole windows. But the paratope does
not read whole windows: in 7K7R it contacts 14 of the 20 EBNA1 residues, with a
long tail (R396 56 heavy-atom contacts, P395 44, P398 38, R397 36, then single
digits, then nothing at S389). Two peptides can be 40% identical overall and
disagree at every contacting position, or 30% identical overall and agree at all
of them. This script measures which, for the three shortlisted pairs.

Statistic
---------
``contact_identity`` = fraction of antibody-contacting EBNA1 positions that the
partner peptide shares. ``weighted`` = the same, weighted by heavy-atom contact
count, so the arginine at 396 counts 56 times and P399 counts 4. Both are reported
because the unweighted version flatters peptides that match only the cheap
positions.

Calibration
-----------
A number like "6 of 14" means nothing on its own, so the same statistic is
computed for random human peptides placed at the same register and length --
2,000 draws per pair. What matters is where a candidate sits in that distribution,
not the raw count.

Register
--------
Candidates are compared through the alignment already stored in
``reports/55_candidates_final.tsv`` (columns ``qaln`` / ``taln``), not through a
new alignment, so the numbers describe the pairs the project actually reported.

Two traps this script is written around, both paid for earlier in the project:
``phmimic.structure.parse_pdb`` does not return chain identifiers, so the chain
column is parsed here; and residues must be keyed by (chain, resseq) anywhere a
dictionary is used, because both Fab chains number from 1 and merging them by
resseq alone silently deletes one chain (100 contacts instead of 262 in 7K7R).

Run from the repository root:

    PYTHONPATH=src env/venv/bin/python src/analysis/62_paratope_contact_identity.py
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
REPORTS = ROOT / "reports"
CUTOFF = 4.5
SHORTLIST = ["P31321", "Q6PJF5", "Q13387"]      # PRKAR1B, RHBDF2, MAPK8IP2
REFERENCE = "Q14CZ8"                            # GlialCAM, the crystallised partner
N_RANDOM = 2000
SEED = 11

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def read_pdb(path: Path) -> dict:
    """First model, with the chain identifier that parse_pdb does not return."""
    xyz, chains, resnums, residues = [], [], [], []
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper() or (
            name[0] if name[0].isalpha() else name[1]).upper()
        if element == "H":
            continue                      # heavy atoms only
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        chains.append(line[21])
        resnums.append(int(line[22:26]))
        residues.append(_THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
    return {"xyz": np.asarray(xyz, dtype=float), "chain": np.asarray(chains),
            "resnum": np.asarray(resnums, dtype=int),
            "resname": np.asarray(residues)}


def paratope(pdb: Path) -> dict[int, int]:
    """Heavy-atom contacts per EBNA1 residue against both Fab chains."""
    atom = read_pdb(pdb)
    by_chain: dict[str, np.ndarray] = {}
    for ch in sorted(set(atom["chain"].tolist())):
        by_chain[ch] = np.where(atom["chain"] == ch)[0]
    peptide_chain = None
    for ch, idx in by_chain.items():
        residues = set(atom["resnum"][idx].tolist())
        if 5 <= len(residues) <= 30:
            peptide_chain = ch
            break
    if peptide_chain is None:
        raise RuntimeError(f"{pdb}: no peptide chain found")
    fab = np.concatenate([idx for ch, idx in by_chain.items()
                          if ch != peptide_chain])
    fab_xyz = atom["xyz"][fab]
    per_residue: dict[int, int] = {}
    for i in by_chain[peptide_chain]:
        dist2 = ((fab_xyz - atom["xyz"][i]) ** 2).sum(axis=1)
        resnum = int(atom["resnum"][i])
        per_residue[resnum] = per_residue.get(resnum, 0) + int(
            (dist2 <= CUTOFF ** 2).sum())
    return per_residue


def score(query_aln: str, target_aln: str, query_start: int,
          profile: dict[int, int]) -> dict:
    qi = query_start - 1
    pairs: list[tuple[int, str, str]] = []
    for q, t in zip(query_aln, target_aln):
        if q != "-":
            qi += 1
        if q != "-" and t != "-":
            pairs.append((qi, q.upper(), t.upper()))
    touching = [(pos, q, t) for pos, q, t in pairs if profile.get(pos, 0) > 0]
    if not touching:
        return {"n_contact_positions": 0, "n_identical": 0,
                "contact_identity": None, "weighted": None}
    identical = sum(1 for _, q, t in touching if q == t)
    total_w = sum(profile[pos] for pos, _, _ in touching)
    ident_w = sum(profile[pos] for pos, q, t in touching if q == t)
    return {
        "n_contact_positions": len(touching),
        "n_identical": identical,
        "contact_identity": round(identical / len(touching), 4),
        "weighted": round(ident_w / total_w, 4),
        "differing_positions": [(pos, q, t) for pos, q, t in touching if q != t],
    }


def load_sequences(path: Path) -> list[str]:
    seqs, parts = [], []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    seqs.append("".join(parts))
                parts = []
            else:
                parts.append(line.strip())
    if parts:
        seqs.append("".join(parts))
    return seqs


def null_distribution(ebna1: str, profile: dict[int, int], span: tuple[int, int],
                      proteome: list[str], n: int, rng: random.Random) -> dict:
    """Random human peptides placed at the same register and length."""
    positions = [p for p in sorted(profile) if profile[p] > 0
                 and span[0] <= p <= span[1]]
    if not positions:
        return {}
    length = span[1] - span[0] + 1
    expected = {p: ebna1[p - 1] for p in positions}
    weights = [profile[p] for p in positions]
    simple, weighted = [], []
    for _ in range(n):
        seq = proteome[rng.randrange(len(proteome))]
        if len(seq) < length:
            continue
        start = rng.randrange(0, len(seq) - length + 1)
        peptide = seq[start:start + length].upper()
        got = [peptide[p - span[0]] for p in positions]
        ident = [g == expected[p] for g, p in zip(got, positions)]
        simple.append(sum(ident) / len(ident))
        weighted.append(sum(w for w, ok in zip(weights, ident) if ok) / sum(weights))
    return {
        "n": len(simple),
        "simple_mean": round(float(np.mean(simple)), 4),
        "simple_p95": round(float(np.percentile(simple, 95)), 4),
        "weighted_mean": round(float(np.mean(weighted)), 4),
        "weighted_p95": round(float(np.percentile(weighted, 95)), 4),
    }


def main() -> int:
    rng = random.Random(SEED)
    pdb = CACHE / "pdb" / "7k7r.pdb"
    if not pdb.exists():
        print(f"missing {pdb}: fetch 7K7R first", file=sys.stderr)
        return 1
    profile = paratope(pdb)
    touching = {k: v for k, v in profile.items() if v > 0}
    print(f"paratope touches {len(touching)} EBNA1 residues / "
          f"{sum(touching.values())} heavy-atom contacts")

    ebna1 = load_sequences(CACHE / "P03211.fasta")[0]
    proteome = load_sequences(CACHE / "human_swissprot.fasta")
    table = pd.read_csv(REPORTS / "55_candidates_final.tsv", sep="\t",
                        dtype=str).fillna("")

    rows, nulls = [], {}
    for acc in SHORTLIST:
        hit = table[table["acc"] == acc]
        if hit.empty:
            print(f"  {acc}: not in the candidate table")
            continue
        row = hit.iloc[0]
        qaln, taln = row["qaln"], row["taln"]
        if not qaln or not taln:
            print(f"  {acc}: no stored alignment")
            continue
        span = (int(row["ebna1_start"]), int(row["ebna1_end"]))
        res = score(qaln, taln, span[0], profile)
        differing = res.pop("differing_positions", [])
        null = null_distribution(ebna1, profile, span, proteome, N_RANDOM, rng)
        nulls[acc] = null
        rows.append({"acc": acc, "gene": row["gene"],
                     "ebna1_span": f"{span[0]}-{span[1]}",
                     "target_span": f"{row['target_start']}-{row['target_end']}",
                     "aln_len": row["aln_len"], "identity_pct": row["identity_pct"],
                     "null_simple_mean": null.get("simple_mean"),
                     "null_simple_p95": null.get("simple_p95"),
                     "null_weighted_mean": null.get("weighted_mean"),
                     **res, "differing": differing})
        print(f"  {row['gene']:9s} {res['n_identical']}/"
              f"{res['n_contact_positions']} contacting positions identical "
              f"(simple {res['contact_identity']} vs random "
              f"{null.get('simple_mean')}; weighted {res['weighted']} vs random "
              f"{null.get('weighted_mean')})")

    # Every candidate that covers a contacting residue, so a single striking
    # number can be read against the other alignments of the same region rather
    # than against random peptides, which is not the population it came from.
    covering = []
    for _, row in table.iterrows():
        qaln, taln = row.get("qaln", ""), row.get("taln", "")
        if not qaln or not taln:
            continue
        span = (int(row["ebna1_start"]), int(row["ebna1_end"]))
        if span[1] < min(touching) or span[0] > max(touching):
            continue
        res = score(qaln, taln, span[0], profile)
        if res["n_contact_positions"] == 0:
            continue
        covering.append({"acc": row["acc"], "gene": row["gene"],
                         "ebna1_span": f"{span[0]}-{span[1]}",
                         "n_contact_positions": res["n_contact_positions"],
                         "n_identical": res["n_identical"],
                         "contact_identity": res["contact_identity"],
                         "weighted": res["weighted"]})
    covering.sort(key=lambda r: (-r["n_identical"], -r["n_contact_positions"]))
    print(f"\n{len(covering)} candidates cover at least one contacting residue")

    # The crystallised partner, compared position by position. The published
    # pairing EBNA1 386-405 <-> GlialCAM 370-389 is a straight 20-mer register,
    # so no alignment is needed to make this honest.
    glialcam = load_sequences(CACHE / "Q14CZ8.fasta")
    reference = None
    if glialcam:
        partner = glialcam[0][370 - 1:389]
        window = ebna1[386 - 1:405]
        res = score(window, partner, 386, profile)
        reference = {"acc": REFERENCE, "gene": "GlialCAM",
                     "ebna1_span": "386-405", "target_span": "370-389",
                     "pairing": "straight 20-mer register, no gaps",
                     **{k: res[k] for k in ("n_contact_positions", "n_identical",
                                            "contact_identity", "weighted")}}
        print(f"  GlialCAM (crystallised partner) "
              f"{res['n_identical']}/{res['n_contact_positions']} identical")

    out = {
        "paratope": {"pdb": "7K7R", "cutoff_angstrom": CUTOFF,
                     "n_touching_residues": len(touching),
                     "total_contacts": int(sum(touching.values())),
                     "profile": {str(k): v for k, v in sorted(touching.items())}},
        "pairs": rows,
        "null": nulls,
        "covering_candidates": covering,
        "reference_binder": reference,
        "note": ("contact_identity compares only EBNA1 positions the antibody "
                 "touches in 7K7R; weighted weights each position by its "
                 "heavy-atom contact count. Whole-window identity and "
                 "contact-level identity are different questions, and only the "
                 "second one is about this antibody. A candidate above the null "
                 "p95 is worth an assay; a candidate above the null p95 while a "
                 "known binder sits far below it is evidence that the statistic "
                 "does not track binding, not evidence that the candidate binds."),
        "reagent_status": ("MS39p2w174 is not a catalogue antibody; the only "
                           "route is recombinant expression of the Fab chains "
                           "deposited in 7K7R / 7UQC."),
    }
    (REPORTS / "62b_paratope_contact_identity.json").write_text(
        json.dumps(out, indent=2, default=str))
    flat = pd.DataFrame(rows).drop(columns=["differing"], errors="ignore")
    flat.to_csv(REPORTS / "62b_paratope_contact_identity.tsv", sep="\t",
                index=False)
    print("wrote reports/62b_paratope_contact_identity.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
