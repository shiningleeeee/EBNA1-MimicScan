#!/usr/bin/env bash
# Verify environment + data provenance. Exit non-zero on any mismatch.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
rc=0
echo "== python =="
env/venv/bin/python - <<'PY' || rc=1
import sys, numpy, pandas, scipy, Bio, yaml, requests, tqdm, pytest
print("python", sys.version.split()[0])
for m in (numpy, pandas, scipy, Bio, yaml, requests, tqdm, pytest):
    print(f"  {m.__name__:10s} {getattr(m,'__version__','?')}")
PY
echo "== mmseqs =="
env/bin/mmseqs version || rc=1
echo "== data checksums =="
( cd data/cache && sha256sum -c <(awk -F'\t' 'NR>1{print $6"  "$1}' MANIFEST.tsv) ) || rc=1
echo "== done (rc=$rc) =="
exit $rc
