#!/usr/bin/env bash
# Full residual suite with mandatory multi-worker fan-out.
# NEVER run the full suite on 1 core. See evidence/E1_FAILURE_GRAPH.md F17.
set -euo pipefail
cd "$(dirname "$0")/.."

# Prefer the largest reliable CPU count (cgroup nproc can lie low; Python
# os.cpu_count and /proc/cpuinfo are fallbacks). Always leave 2 cores free.
_n_nproc="$(nproc 2>/dev/null || echo 0)"
_n_py="$(python3 -c 'import os; print(os.cpu_count() or 0)' 2>/dev/null || echo 0)"
_n_procfile="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 0)"
_n_max="${_n_nproc}"
for _c in "${_n_py}" "${_n_procfile}"; do
  if [[ "${_c}" -gt "${_n_max}" ]]; then _n_max="${_c}"; fi
done
if [[ "${_n_max}" -lt 4 ]]; then
  echo "pytest_full: FATAL — only ${_n_max} CPUs visible; refusing single-core suite (F17)" >&2
  exit 2
fi
W="${PYTEST_WORKERS:-$(( _n_max - 2 ))}"
if [[ "${W}" -lt 2 ]]; then W=2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "pytest_full: W=${W}  cpus(nproc=${_n_nproc},py=${_n_py},getconf=${_n_procfile},max=${_n_max})" >&2
exec python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py \
  -n "${W}" -q --tb=line "$@"
