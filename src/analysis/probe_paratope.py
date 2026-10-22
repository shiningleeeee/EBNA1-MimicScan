"""Crystal paratope: register, contacts and what the antibody actually reads.

7K7R = Fab MS39p2w174 + EBNA1 peptide 386-405    (X-ray 2.50 A, modelled 389-402)
7UQC = Fab MS39p2w174 + phospho-GlialCAM 370-389 (X-ray 2.65 A)

One antibody, two peptides, both true binders (Lanz 2022). That makes the
GlialCAM complex a blind positive for anything derived from the EBNA1 complex.

Three questions, answered from the coordinates alone:
  1. Do the two Fabs present the same paratope? (CA RMSD over Fv)
  2. At what register are the two 20-mers bound -- and what identity does that
     register show? (the pair is the experimentally proven case, so this number
     is the ground truth for "how much sequence identity is required")
  3. Which peptide positions does the antibody actually contact, and is the
     contact pattern identical between the two peptides?
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "PTR": "Y", "SEP": "S", "TPO": "T",
}
SKIP = {"HOH", "DOD", "WAT"}


def parse_pdb(path):
    """-> {chain: {resseq: {"resname", "atoms": [(name, element, xyz)]}}}"""
    out: dict[str, dict[int, dict]] = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[16] not in (" ", "A"):
            continue
        resname = line[17:20].strip()
        if resname in SKIP:
            continue
        name = line[12:16].strip()
        chain = line[21]
        resseq = int(line[22:26])
        xyz = np.array((float(line[30:38]), float(line[38:46]),
                        float(line[46:54])))
        element = (line[76:78].strip() or name[0]).upper()
        out.setdefault(chain, {}).setdefault(
            resseq, {"resname": resname, "atoms": []})["atoms"].append(
                (name, element, xyz))
    return out


def one(chain, r):
    return AA3.get(chain[r]["resname"], "X")


def ca(chain, r):
    for name, _e, xyz in chain[r]["atoms"]:
        if name == "CA":
            return xyz
    return None


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _s, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def contact_profile(peptide, fab, cutoff=4.5):
    F = np.array([xyz for r in fab.values() for _n, _e, xyz in r["atoms"]])
    prof = {}
    for r in sorted(peptide):
        A = np.array([xyz for _n, _e, xyz in peptide[r]["atoms"]])
        d = np.linalg.norm(A[:, None, :] - F[None, :, :], axis=2)
        prof[r] = int((d < cutoff).sum())
    return prof


def main():
    d = Path(sys.argv[1] if len(sys.argv) > 1
             else "/workspace/PH-Mimic/data/cache/pdb")
    e = parse_pdb(d / "7k7r.pdb")     # EBNA1
    g = parse_pdb(d / "7uqc.pdb")     # GlialCAM

    pep_e, pep_g = e["C"], g["C"]
    key_e, key_g = sorted(pep_e), sorted(pep_g)
    seq_e = "".join(one(pep_e, r) for r in key_e)
    seq_g = "".join(one(pep_g, r) for r in key_g)
    print(f"7K7R EBNA1     {len(seq_e)} modelled residues "
          f"{key_e[0]}-{key_e[-1]}:  {seq_e}")
    print(f"7UQC GlialCAM  {len(seq_g)} modelled residues "
          f"{key_g[0]}-{key_g[-1]}:  {seq_g}")

    # ---- 1. same paratope?
    common = []
    for ch in ("A", "B"):
        for r in sorted(set(e[ch]) & set(g[ch])):
            pa, pb = ca(e[ch], r), ca(g[ch], r)
            if pa is not None and pb is not None:
                common.append((pa, pb))
    P = np.array([c[0] for c in common])
    Q = np.array([c[1] for c in common])
    R, t = kabsch(P, Q)                       # 7K7R frame -> 7UQC frame
    moved = (R @ P.T).T + t
    rmsd = float(np.sqrt(((moved - Q) ** 2).sum(1).mean()))
    print(f"\nFv CA superposition: {len(common)} residues, RMSD {rmsd:.2f} A "
          f"-> same paratope")

    # ---- 2. register: scan offsets of EBNA1 numbering vs GlialCAM numbering
    print("\nregister scan (mean CA distance after superposing the two Fabs):")
    best = None
    for off in range(8, 22):                  # EBNA1 residue = GlialCAM + off
        ds = []
        for rg in key_g:
            re_ = rg + off
            if re_ not in pep_e:
                continue
            a, b_ = ca(pep_e, re_), ca(pep_g, rg)
            if a is None or b_ is None:
                continue
            ds.append(float(np.linalg.norm((R @ a + t) - b_)))
        if len(ds) >= 10:
            m = float(np.mean(ds))
            if best is None or m < best[1]:
                best = (off, m, len(ds))
            print(f"  offset +{off:2d}: n={len(ds):2d} mean CA d = {m:5.2f} A")
    off, mdist, n = best
    print(f"\nbest register: GlialCAM residue r  <->  EBNA1 residue r+{off} "
          f"(mean CA distance {mdist:.2f} A over {n} positions)")

    pairs = [(rg, rg + off) for rg in key_g if rg + off in pep_e]
    ident = [(rg, re_) for rg, re_ in pairs if one(pep_g, rg) == one(pep_e, re_)]
    print(f"identity at that register: {len(ident)}/{len(pairs)} = "
          f"{len(ident)/len(pairs):.0%}")
    print("  " + "  ".join(f"{one(pep_g, rg)}{rg}/{one(pep_e, re_)}{re_}"
                           f"{'=' if one(pep_g, rg) == one(pep_e, re_) else ''}"
                           for rg, re_ in pairs))

    # ---- 3. what does the antibody read?
    pe = contact_profile(pep_e, {**e["A"], **e["B"]})
    pg = contact_profile(pep_g, {**g["A"], **g["B"]})
    print("\nparatope contacts (Fab heavy atoms within 4.5 A per peptide residue)")
    print("  position        EBNA1          GlialCAM")
    for rg, re_ in pairs:
        print(f"  {re_:>4d}/{rg:<4d}   {one(pep_e, re_)}{re_:>5d}:{pe[re_]:<4d}"
              f"   {one(pep_g, rg)}{rg:>5d}:{pg[rg]}")
    tot_e = sum(pe.values())
    tot_g = sum(pg.values())
    top_e = sorted(pe.items(), key=lambda kv: -kv[1])[:3]
    top_g = sorted(pg.items(), key=lambda kv: -kv[1])[:3]
    print(f"\n  total contacts: EBNA1 {tot_e}, GlialCAM {tot_g}")
    print(f"  top EBNA1 residues: "
          f"{[(one(pep_e, r) + str(r), n) for r, n in top_e]}")
    print(f"  top GlialCAM residues: "
          f"{[(one(pep_g, r) + str(r), n) for r, n in top_g]}")
    unmodelled = [r for r in range(386, 406) if r not in pep_e]
    print(f"\n  EBNA1 386-405 residues NOT modelled in 7K7R: {unmodelled}")


if __name__ == "__main__":
    main()
