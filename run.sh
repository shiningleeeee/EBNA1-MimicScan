#!/usr/bin/env bash
# Single entry point. Everything the report contains is produced by this script.
#
#   bash run.sh                  full pipeline
#   bash run.sh --from s40       resume from a stage
#   bash run.sh --only s40_priority
#   bash run.sh --list
set -euo pipefail
cd "$(dirname "$0")"

export PYTHONPATH=src
PY=env/venv/bin/python

if [ ! -x "$PY" ]; then
    echo "error: $PY not found -- run scripts/00_check_env.sh first" >&2
    exit 2
fi

mkdir -p reports logs
make -s all

exec "$PY" -m phmimic.runner "$@"
