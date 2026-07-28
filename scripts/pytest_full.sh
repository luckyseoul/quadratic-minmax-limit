#!/usr/bin/env bash
# Full residual suite with mandatory multi-worker fan-out.
# NEVER run the full suite on 1 core. See evidence/E1_FAILURE_GRAPH.md F17.
set -euo pipefail
cd "$(dirname "$0")/.."
W="${PYTEST_WORKERS:-$(( $(nproc) - 2 ))}"
if [[ "${W}" -lt 2 ]]; then W=2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "pytest_full: W=${W}  nproc=$(nproc)" >&2
exec python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py \
  -n "${W}" -q --tb=line "$@"
