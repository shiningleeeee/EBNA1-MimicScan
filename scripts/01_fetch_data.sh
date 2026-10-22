#!/usr/bin/env bash
# Fetch the two sequence inputs, then verify them against data/MANIFEST.tsv.
# Everything else (annotations, AlphaFold models, PDB entries, IEDB records) is
# downloaded lazily by the stage that needs it and cached under data/cache/.
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1

BASE=$(awk -F'\t' 'NR>1 && $1=="P03211.fasta"{print $2}' data/MANIFEST.tsv)
PROT=$(awk -F'\t' 'NR>1 && $1=="human_swissprot.fasta"{print $2}' data/MANIFEST.tsv)
[ -n "$BASE" ] && [ -n "$PROT" ] || { echo "data/MANIFEST.tsv: cannot read source URLs" >&2; exit 1; }

mkdir -p data/cache
for pair in "P03211.fasta|$BASE" "human_swissprot.fasta|$PROT"; do
    name=${pair%%|*}; url=${pair#*|}
    if [ -s "data/cache/$name" ]; then
        echo "have data/cache/$name"
    else
        echo "fetching $name"
        curl -fsSL --retry 3 --retry-delay 2 "$url" -o "data/cache/$name"
    fi
done

echo "== verifying against data/MANIFEST.tsv =="
( cd data/cache && awk -F'\t' 'NR>1{print $6"  "$1}' ../MANIFEST.tsv | sha256sum -c - )
echo "ok: inputs match the manifest"
