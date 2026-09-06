#!/usr/bin/env bash
set -euo pipefail
cd /tmp/original-mo-final-docs-gate.uCXyc4/repo
sha256sum --check ../source.sha256
gate_status=0
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -n 14 tests/test_original_mo_status.py tests/test_main_chain_docs.py::test_main_chain_L_open_and_docs_ok > ../pytest_verified.log 2>&1 || gate_status=$?
cat ../pytest_verified.log
exit "$gate_status"
