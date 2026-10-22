#!/usr/bin/env python3
"""47_pair_null_batch.py — selection-aware per-pair composition null.

Why this exists
---------------
A per-pair composition null (pair_null) answers "given these two local
compositions, how often would random orderings match this well?".  Applied to
the alignments that a whole-proteome search already picked as *best*, it is
biased: those pairs were selected for being the maximum of a large search, so
of course they beat their own shuffles.

The bias is removed by running the SAME full procedure -- search the proteome,
take every reported alignment, score each one against its own composition
null -- on composition-preserving shuffles of EBNA1.  Whatever the shuffles
report is the selection-aware expectation; the real query is only interesting
where it stands outside that.

Design
------
  1. real EBNA1 + N composition-preserving shuffles in ONE multi-FASTA query
     DB, searched in ONE mmseqs invocation (removes batch effects entirely).
  2. every reported alignment of every query is scored by pair_null against
     its own two-sided composition null (B shuffles through the exact same
     SW kernel).
  3. per query summary: n_align, n_p_lt_1e-4, n_p_lt_1e-3, best/median p.
  4. the real query is compared against the empirical distribution of the
     shuffles -> p_selection.

Usage:
  env/venv/bin/python src/47_pair_null_batch.py --n-shuffles 20 --b 2000
"""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "pair_null_batch"
REPORTS = ROOT / "reports"
LOGS = ROOT / "logs"
MMSEQS = ROOT / "env" / "bin" / "mmseqs"
PAIR_NULL = ROOT / "env" / "bin" / "pair_null"
QUERY = ROOT / "data" / "cache" / "P03211.fasta"
PROTEOME = ROOT / "data" / "cache" / "human_swissprot.fasta"


def read_first_fasta(path: Path) -> tuple[str, str]:
    acc, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    break
                acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
            else:
                seq.append(line)
    return acc, "".join(seq)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def run(cmd: list[str], log: Path) -> None:
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
        fh.flush()
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh, stderr=subprocess.STDOUT)


def mmseqs_search(query_fa: Path, out_m8: Path) -> None:
    work = WORK / "mm"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    qdb = work / "qdb"
    tdb = work / "tdb"
    res = work / "res"
    run([MMSEQS, "createdb", query_fa, qdb], LOGS / "47_batch_null.log")
    run([MMSEQS, "createdb", PROTEOME, tdb], LOGS / "47_batch_null.log")
    run(
        [MMSEQS, "search", qdb, tdb, res, work / "tmp",
         "-s", "7.5", "--max-seqs", "20000", "-e", "10000",
         "--comp-bias-corr", "1", "-a", "--threads", "14"],
        LOGS / "47_batch_null.log",
    )
    run(
        [MMSEQS, "convertalis", qdb, tdb, res, out_m8,
         "--format-output",
         "query,target,evalue,bits,qstart,qend,tstart,tend,alnlen,pident,qaln,taln"],
        LOGS / "47_batch_null.log",
    )


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-shuffles", type=int, default=20)
    ap.add_argument("--b", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    WORK.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    acc, seq = read_first_fasta(QUERY)
    print(f"[47] query {acc} {len(seq)} aa; {args.n_shuffles} shuffles; B={args.b}")

    rng = random.Random(args.seed)
    records: list[tuple[str, str]] = [("real", seq)]
    for i in range(1, args.n_shuffles + 1):
        chars = list(seq)
        rng.shuffle(chars)
        records.append((f"shuf{i:02d}", "".join(chars)))

    qfa = WORK / "queries.fasta"
    write_fasta(qfa, records)

    m8 = WORK / "all.m8"
    if not m8.exists():
        mmseqs_search(qfa, m8)
    cols = ["query", "target", "evalue", "bits", "qstart", "qend",
            "tstart", "tend", "alnlen", "pident", "qaln", "taln"]
    allhits = pd.read_csv(m8, sep="\t", names=cols, dtype={"query": str, "target": str})
    print(f"[47] {len(allhits)} alignments over {allhits['query'].nunique()} queries")
    print(allhits.groupby("query").size().sort_values(ascending=False).to_string())

    rows = []
    by_query: dict[str, pd.DataFrame] = {q: g for q, g in allhits.groupby("query")}
    for name, seq_i in records:
        sub = by_query.get(name)
        if sub is None or sub.empty:
            rows.append(dict(query=name, n_align=0, n_p_lt_1e4=0, n_p_lt_1e3=0,
                             n_p_lt_1e2=0, best_p=float("nan"), med_p=float("nan")))
            print(f"[47] {name}: no alignments")
            continue
        cands = WORK / f"cands_{name}.tsv"
        sub[["target", "qstart", "qend", "tstart", "tend", "alnlen", "pident"]].to_csv(
            cands, sep="\t", index=False
        )
        qonly = WORK / f"q_{name}.fasta"
        write_fasta(qonly, [(name, seq_i)])
        outp = WORK / f"pn_{name}.tsv"
        run([PAIR_NULL, "--query", qonly, "--db", PROTEOME, "--cands", cands,
             "--b", args.b, "--seed", args.seed, "--flank", "5", "--out", outp],
            LOGS / "47_batch_null.log")
        pn = pd.read_csv(outp, sep="\t")
        rows.append(dict(query=name, n_align=len(pn),
                         n_p_lt_1e4=int((pn["p_pair"] < 1e-4).sum()),
                         n_p_lt_1e3=int((pn["p_pair"] < 1e-3).sum()),
                         n_p_lt_1e2=int((pn["p_pair"] < 1e-2).sum()),
                         best_p=float(pn["p_pair"].min()),
                         med_p=float(pn["p_pair"].median())))
        print(f"[47] {name}: n={len(pn)} p<1e-4={rows[-1]['n_p_lt_1e4']} "
              f"p<1e-3={rows[-1]['n_p_lt_1e3']} best={rows[-1]['best_p']:.2e}")

    summ = pd.DataFrame(rows)
    summ.to_csv(REPORTS / "47_pair_null_batch.csv", index=False)

    real = summ[summ["query"] == "real"].iloc[0]
    nulls = summ[summ["query"] != "real"]
    print("\n=== selection-aware comparison ===")
    print(f"real  : n_align={real.n_align}  p<1e-4={real.n_p_lt_1e4}  "
          f"p<1e-3={real.n_p_lt_1e3}  best_p={real.best_p:.2e}")
    for col in ("n_align", "n_p_lt_1e4", "n_p_lt_1e3", "n_p_lt_1e2", "best_p"):
        vals = nulls[col].astype(float)
        obs = float(real[col])
        if col == "best_p":
            ge = int((vals <= obs).sum())
        else:
            ge = int((vals >= obs).sum())
        p = (1 + ge) / (1 + len(vals))
        print(f"{col:12s}: null mean={vals.mean():10.4g} sd={vals.std():10.4g} "
              f"min={vals.min():10.4g} max={vals.max():10.4g} | real={obs:10.4g} -> p={(p):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
