#!/usr/bin/env bash
# Resume p=13 k=6 gauged GPU enum on whatever CUDA GPU is here (V100, H100, …).
#
# Progress is GAUGE_OUT/orb{tvidx}.npy — existing files are skipped.
# Fast orbit skip is LOAD_TASKS pickle (k6_tasks.pkl).
# Single status file: RESUME_JSON (k6_resume.json).
#
# Copy onto a rental box:
#   k6_code/          (these scripts)
#   k6_tasks.pkl
#   k6_resume.json
#   k6_gpu_out/       (the orb npy files)
#
# Then:
#   PALEY_P=13 GAUGE_OUT=.../k6_gpu_out LOAD_TASKS=.../k6_tasks.pkl \
#     GPU_WORKERS=2 ./run_k6_p13_resume.sh
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
export PALEY_P="${PALEY_P:-13}"
export GAUGE_OUT="${GAUGE_OUT:-/mnt/storage/e1work/maxplus_p13/k6_gpu_out}"
export LOAD_TASKS="${LOAD_TASKS:-/mnt/storage/e1work/maxplus_p13/k6_tasks.pkl}"
export RESUME_JSON="${RESUME_JSON:-/mnt/storage/e1work/maxplus_p13/k6_resume.json}"
export GPU_WORKERS="${GPU_WORKERS:-2}"
export GPU_MEM_FRAC="${GPU_MEM_FRAC:-0.70}"
export GEN_CAP="${GEN_CAP:-40000000}"
export ENUM_SHARD_MOD="${ENUM_SHARD_MOD:-1}"
export ENUM_SHARD_REM="${ENUM_SHARD_REM:-0}"
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export RESUME_K=6
nw="${GPU_WORKERS}"
# H100 80GiB: raise GEN_CAP / workers if you want. Kernel is arch-agnostic.
if [[ ! -f "$LOAD_TASKS" ]]; then
  echo "missing LOAD_TASKS pickle: $LOAD_TASKS" >&2
  exit 1
fi
mkdir -p "$GAUGE_OUT"
echo "resume $(date -Is) GAUGE_OUT=$GAUGE_OUT orbs=$(ls "$GAUGE_OUT"/orb*.npy 2>/dev/null | wc -l) workers=$nw"
cd "$DIR"
exec python3 -u run_kgauged.py 6 "$nw"
