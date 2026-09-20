#!/usr/bin/env bash
set -euo pipefail
SSH=(ssh -F /home/nick/.ssh/mesh.config -o ConnectTimeout=5 -o BatchMode=yes nuka)
DIR=/home/nick/quadratic-minmax-limit/evidence/k50_nuka_gpu_20260919
LOG=$DIR/k50_gpu.log
PIDF=$DIR/k50_gpu.pid
diag=/tmp/watch_k50_gpu_$$.log
while :; do
  if ! "${SSH[@]}" "PID=\$(cat $PIDF 2>/dev/null || echo 0); echo PID:\$PID; kill -0 \$PID 2>/dev/null; echo alive:\$?; tail -5 $LOG" >"$diag" 2>&1; then
    echo "FAILED: ssh nuka failed"
    exit 1
  fi
  if grep -q '"phi":' "$diag"; then
    echo "DONE: $(grep '"phi":' "$diag" | tail -1)"
    exit 0
  fi
  if grep -q 'alive:1' "$diag"; then
    echo "FAILED: pid dead without phi. $(tail -3 $diag | tr '\n' ' ')"
    exit 1
  fi
  sleep 30
done
