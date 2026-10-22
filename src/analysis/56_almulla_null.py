#!/usr/bin/env python
"""
56_almulla_null.py — external falsification test.

Question: Almulla et al. (Immun Inflamm Dis 2026;14(5):e70445, PMID 42125999,
PMC13169152) report 91 shared pentapeptides between 33 EBV/HHV-6 proteins and
10 human CNS antigens, and read them as evidence of molecular mimicry.

Test: are 91 shared pentapeptides surprising, once each protein's amino-acid
composition is held fixed?

Null models (each preserves composition, destroys only residue ORDER):
  A  global shuffle      — uniform random permutation of each protein
  B  window shuffle w=10 — permute within non-overlapping 10-aa windows
  C  window shuffle w=20 — permute within non-overlapping 20-aa windows

Statistics compared against the null:
  total_pairs_shared   number of the 330 (33x10) pairs sharing >=1 pentapeptide
  total_shared_types   sum over pairs of |pentamers(a) & pentamers(b)|
  max_pair_shared      largest per-pair shared count
  reported_pairs_hit   how many of the pairs actually listed in Table 2 are
                       reproduced by the shuffled sequences

Inputs (all cached, provenance in data/cache/almulla/sequences_manifest.json):
  data/cache/almulla/table1_viral.tsv
  data/cache/almulla/table1_cns.tsv
  data/cache/almulla/table2_pentapeptides.tsv
  data/cache/almulla/<ACC>.fasta
"""
from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "almulla"
REPORTS = ROOT / "reports"
K = 5


def read_fasta_seq(acc: str) -> str:
    lines = (CACHE / f"{acc}.fasta").read_text().splitlines()
    return "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()


def pentamers(seq: str) -> set[str]:
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def shuffle_global(seq: str, rng: random.Random) -> str:
    a = list(seq)
    rng.shuffle(a)
    return "".join(a)


def shuffle_windows(seq: str, w: int, rng: random.Random) -> str:
    out = []
    for i in range(0, len(seq), w):
        chunk = list(seq[i:i + w])
        rng.shuffle(chunk)
        out.extend(chunk)
    return "".join(out)


def statistics(sets_v: dict[str, set[str]], sets_c: dict[str, set[str]],
               pairs: list[tuple[str, str]]) -> tuple[int, int, int, int]:
    n_pairs_shared = 0
    total_types = 0
    max_pair = 0
    for va, ca in pairs:
        n = len(sets_v[va] & sets_c[ca])
        if n:
            n_pairs_shared += 1
            total_types += n
            max_pair = max(max_pair, n)
    return n_pairs_shared, total_types, max_pair


def main(argv: list[str]) -> int:
    b = 1000
    seed = 42
    for i, a in enumerate(argv):
        if a == "--b":
            b = int(argv[i + 1])
        elif a == "--seed":
            seed = int(argv[i + 1])

    viral = list(csv.DictReader(open(CACHE / "table1_viral.tsv"), delimiter="\t"))
    cns = list(csv.DictReader(open(CACHE / "table1_cns.tsv"), delimiter="\t"))
    reported = list(csv.DictReader(open(CACHE / "table2_verified.tsv"), delimiter="\t"))

    vs = {r["acc"]: read_fasta_seq(r["acc"]) for r in viral}
    cs = {r["acc"]: read_fasta_seq(r["acc"]) for r in cns}
    pairs = [(v["acc"], c["acc"]) for v in viral for c in cns]
    assert len(pairs) == 330, len(pairs)

    # ---- observed ----
    obs_v = {a: pentamers(s) for a, s in vs.items()}
    obs_c = {a: pentamers(s) for a, s in cs.items()}
    o_pairs, o_types, o_max = statistics(obs_v, obs_c, pairs)

    reported_pairs = {(r["viral_acc"], r["cns_acc"]) for r in reported}
    missing = [p for p in reported_pairs if p[0] not in obs_v or p[1] not in obs_c]
    assert not missing, missing
    o_reportedpairs_hit = sum(
        1 for va, ca in reported_pairs if obs_v[va] & obs_c[ca]
    )

    # per-pair reported vs observed detail (for the record)
    detail = []
    for r in reported:
        va, ca = r["viral_acc"], r["cns_acc"]
        detail.append({
            "virus": r["virus"], "viral_protein": r["viral_protein"], "viral_acc": va,
            "cns_antigen": r["cns_antigen"], "cns_acc": ca,
            "pentapeptide": r["pentapeptide"],
            "in_viral": r["in_viral"], "in_cns": r["in_cns"],
        })

    # ---- null ----
    results = {}
    for label, fn in (
        ("global", lambda s, rng: shuffle_global(s, rng)),
        ("window10", lambda s, rng: shuffle_windows(s, 10, rng)),
        ("window20", lambda s, rng: shuffle_windows(s, 20, rng)),
    ):
        rng = random.Random(seed)
        null_pairs, null_types, null_max, null_rep = [], [], [], []
        for _ in range(b):
            nv = {a: pentamers(fn(s, rng)) for a, s in vs.items()}
            nc = {a: pentamers(fn(s, rng)) for a, s in cs.items()}
            p, t, m = statistics(nv, nc, pairs)
            null_pairs.append(p)
            null_types.append(t)
            null_max.append(m)
            null_rep.append(sum(1 for va, ca in reported_pairs if nv[va] & nc[ca]))

        def summ(obs, arr):
            arr_s = sorted(arr)
            ge = sum(1 for x in arr_s if x >= obs)
            return {
                "observed": obs,
                "null_mean": sum(arr_s) / len(arr_s),
                "null_sd": (sum((x - sum(arr_s) / len(arr_s)) ** 2 for x in arr_s) / (len(arr_s) - 1)) ** 0.5,
                "null_min": arr_s[0],
                "null_max": arr_s[-1],
                "null_p95": arr_s[int(0.95 * (len(arr_s) - 1))],
                "n_null_ge": ge,
                "p_one_sided": (1 + ge) / (1 + len(arr_s)),
            }

        results[label] = {
            "b": b,
            "seed": seed,
            "pairs_with_shared": summ(o_pairs, null_pairs),
            "total_shared_types": summ(o_types, null_types),
            "max_pair_shared": summ(o_max, null_max),
            "reported_pairs_reproduced": summ(o_reportedpairs_hit, null_rep),
        }
        r = results[label]
        print(f"[{label:8s}] pairs_with_shared obs={o_pairs:5d} null={r['pairs_with_shared']['null_mean']:8.1f} p={r['pairs_with_shared']['p_one_sided']:.4f} | "
              f"total_shared obs={o_types:6d} null={r['total_shared_types']['null_mean']:9.1f} p={r['total_shared_types']['p_one_sided']:.4f} | "
              f"max_pair obs={o_max:3d} null_max={r['max_pair_shared']['null_max']:3d} p={r['max_pair_shared']['p_one_sided']:.4f}",
              flush=True)

    out = {
        "question": "are 91 shared pentapeptides between EBV/HHV-6 and human CNS proteins surprising under composition-matched nulls?",
        "paper": "Almulla AF et al., Immun Inflamm Dis 2026;14(5):e70445, PMID 42125999, PMC13169152",
        "k": K,
        "n_viral": len(viral), "n_cns": len(cns), "n_pairs": len(pairs),
        "observed": {
            "pairs_with_shared": o_pairs,
            "total_shared_types": o_types,
            "max_pair_shared": o_max,
            "reported_pairs_listed": len(reported_pairs),
            "reported_pairs_reproduced": o_reportedpairs_hit,
            "reported_pentapeptides": len(reported),
            "reported_verified_in_both": sum(1 for r in detail if r["in_viral"] == "1" and r["in_cns"] == "1"),
        },
        "nulls": results,
        "detail": detail,
    }
    (REPORTS / "56_almulla_null.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("\nobserved: pairs_with_shared=%d total_shared=%d max_pair=%d (of %d pairs)"
          % (o_pairs, o_types, o_max, len(pairs)))
    print("wrote reports/56_almulla_null.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
