#!/usr/bin/env python
"""
57_almulla_specificity.py — second falsification test.

Almulla et al. (PMID 42125999) frame their 91 shared pentapeptides as evidence
for a *virus -> CNS autoimmunity* link. That is a claim about the CNS antigen set
being special. This script tests it directly.

Test: are the 10 CNS antigens more pentapeptide-similar to the 33 EBV/HHV-6
proteins than length- and composition-matched human proteins?

Control set construction (both confounds controlled):
  * length within +/-15% of the CNS antigen
  * amino-acid composition within the closest 20% of that length pool
  * sampled uniformly from that restricted pool, 200 replicates

Statistic: total number of distinct shared pentapeptides between the 10 drawn
human proteins and the 33 viral proteins (observed: 99 for the real CNS set).
"""
from __future__ import annotations

import csv
import json
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "almulla"
HUMAN = ROOT / "data" / "cache" / "human_swissprot.fasta"
REPORTS = ROOT / "reports"
K = 5
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    acc = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc:
                out[acc] = "".join(buf).upper()
            acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        out[acc] = "".join(buf).upper()
    return out


def pentamers(seq: str) -> set[str]:
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def compvec(seq: str) -> Counter:
    c = Counter(seq)
    n = max(1, len(seq))
    return {a: c.get(a, 0) / n for a in AA}


def l1(a: dict, b: dict) -> float:
    return sum(abs(a[k] - b[k]) for k in AA)


def main(argv: list[str]) -> int:
    reps = 200
    seed = 42
    for i, a in enumerate(argv):
        if a == "--reps":
            reps = int(argv[i + 1])
        elif a == "--seed":
            seed = int(argv[i + 1])

    viral = list(csv.DictReader(open(CACHE / "table1_viral.tsv"), delimiter="\t"))
    cns = list(csv.DictReader(open(CACHE / "table1_cns.tsv"), delimiter="\t"))
    vs = {}
    for r in viral:
        lines = (CACHE / f"{r['acc']}.fasta").read_text().splitlines()
        vs[r["acc"]] = "".join(l for l in lines[1:] if not l.startswith(">")).upper()
    vset = {a: pentamers(s) for a, s in vs.items()}

    cns_seq = {}
    for r in cns:
        lines = (CACHE / f"{r['acc']}.fasta").read_text().splitlines()
        cns_seq[r["acc"]] = "".join(l for l in lines[1:] if not l.startswith(">")).upper()
    cset = {a: pentamers(s) for a, s in cns_seq.items()}
    observed = sum(len(vset[v] & cset[c]) for v in vset for c in cset)

    human = read_fasta(HUMAN)
    hset = {a: pentamers(s) for a, s in human.items()}
    hcomp = {a: compvec(s) for a, s in human.items()}
    hlen = {a: len(s) for a, s in human.items()}
    print(f"human proteome: {len(human)} proteins; viral {len(vset)}; CNS {len(cset)}")
    print(f"observed total shared pentapeptides (real CNS set vs viral): {observed}")

    rng = random.Random(seed)
    pools: dict[str, list[str]] = {}
    for acc, seq in cns_seq.items():
        L = len(seq)
        lo, hi = int(0.85 * L), int(1.15 * L)
        cand = [a for a in human if lo <= hlen[a] <= hi and a not in cns_seq]
        cv = compvec(seq)
        cand.sort(key=lambda a: l1(cv, hcomp[a]))
        keep = cand[: max(1, int(0.20 * len(cand)))]
        pools[acc] = keep
        print(f"  pool[{acc}] L={L:4d} len-matched={len(cand):5d} comp+len matched={len(keep):5d}")

    null = []
    for _ in range(reps):
        tot = 0
        for acc in cns_seq:
            pick = rng.choice(pools[acc])
            tot += sum(len(vset[v] & hset[pick]) for v in vset)
        null.append(tot)
    null.sort()
    ge = sum(1 for x in null if x >= observed)
    mean = sum(null) / len(null)
    sd = (sum((x - mean) ** 2 for x in null) / (len(null) - 1)) ** 0.5

    print(f"\nnull ({reps} draws of 10 length+composition-matched human proteins):")
    print(f"  mean {mean:.1f}  sd {sd:.1f}  min {null[0]}  max {null[-1]}  p95 {null[int(0.95*(len(null)-1))]}")
    print(f"  observed {observed}  |  n_null_ge {ge}  p_one_sided {(1+ge)/(1+reps):.4f}")

    out = {
        "question": "are the 10 CNS antigens more pentapeptide-similar to EBV/HHV-6 proteins than length+composition-matched human proteins?",
        "paper": "Almulla AF et al., Immun Inflamm Dis 2026;14(5):e70445, PMID 42125999",
        "k": K, "reps": reps, "seed": seed,
        "observed_total_shared": observed,
        "null_mean": mean, "null_sd": sd, "null_min": null[0], "null_max": null[-1],
        "null_p95": null[int(0.95 * (len(null) - 1))],
        "n_null_ge": ge, "p_one_sided": (1 + ge) / (1 + reps),
        "cns_lengths": {a: len(s) for a, s in cns_seq.items()},
    }
    (REPORTS / "57_almulla_specificity.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("wrote reports/57_almulla_specificity.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
