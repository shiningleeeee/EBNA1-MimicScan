"""How rare is the motif the paratope actually reads?

probe_paratope.py established, from two crystal structures of the same Fab:
  - the two peptides bind in a fixed register (GlialCAM r <-> EBNA1 r+16),
  - at that register they are 5/14 identical (6/20 over the published 20-mers),
  - the paratope's contact mass sits almost entirely on ONE arginine
    (EBNA1 R396: 51 of 100 contacts; GlialCAM R380: 52 of 83),
    flanked by prolines at +-2 and +-4.

So the shared read across both true binders is a three-position pattern:
P at register 0, R at register +2, P at register +4  ->  "P.R.P".

This script asks the only question that decides whether that pattern is usable
as evidence: how many places in the human proteome present it?

A pattern that occurs in most proteins cannot distinguish a mimic from anything
else, no matter how well it is anchored in a crystal structure.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from phmimic.mining import load_fasta  # noqa: E402

PATTERN = re.compile(r"(?=(P.R.P))")   # overlapping
CONTROLS = {
    "full read P.R.P": re.compile(r"(?=(P.R.P))"),
    "wider read PP.R.P": re.compile(r"(?=(PP.R.P))"),
    "published seed PPRRPP": re.compile(r"(?=PPRRPP)"),
    "published seed RRPFF": re.compile(r"(?=RRPFF)"),
    "Hecker 17-mer GSPPRRPPPGRRPFFHP": re.compile(r"(?=GSPPRRPPPGRRPFFHP)"),
    "control AAAAA": re.compile(r"(?=AAAAA)"),
    "control WWWWW": re.compile(r"(?=WWWWW)"),
}


def main():
    db = Path(sys.argv[1] if len(sys.argv) > 1
              else "data/cache/human_swissprot.fasta")
    records = load_fasta(db)          # {accession: sequence}
    n_prot = len(records)
    total_aa = sum(len(s) for s in records.values())
    print(f"proteome: {n_prot} proteins, {total_aa} residues")

    for label, rx in CONTROLS.items():
        hits = 0
        prot_with = 0
        lengths = []
        for seq in records.values():
            n = len(rx.findall(seq))
            if n:
                hits += n
                prot_with += 1
                lengths.append(n)
        per_prot = prot_with / n_prot
        print(f"\n{label}")
        print(f"  occurrences      : {hits}")
        print(f"  proteins with >=1: {prot_with} / {n_prot} ({per_prot:.1%})")
        if lengths:
            lengths.sort()
            print(f"  median per protein: {lengths[len(lengths)//2]}")

    # where do the two known binders sit?
    for acc, lo, hi, name in (("P03211", 394, 398, "EBNA1 394-398"),
                              ("Q14CZ8", 378, 382, "GlialCAM 378-382")):
        seq = records.get(acc)
        if seq is None:
            print(f"\n{name}: {acc} not in the human proteome (expected for EBNA1)")
            continue
        sub = seq[lo - 1:hi]
        print(f"\n{name}: {sub}  (P.R.P satisfied: "
              f"{bool(re.fullmatch('P.R.P', sub))})")


if __name__ == "__main__":
    main()
