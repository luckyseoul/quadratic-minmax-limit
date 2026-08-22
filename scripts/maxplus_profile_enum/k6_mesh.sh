#!/usr/bin/env bash
# k6-mesh — p=13 k=6 orbit mesh. See k6-mesh.1 (k6-mesh man).
# Soft stop always finishes the current orbit(s); it never SIGKILLs an in-flight
# outer. Start of one worker never stops another.
set -euo pipefail
ROOT="${K6_ROOT:-/mnt/storage/e1work/maxplus_p13}"
# Default to this script's directory so dash/start still work after mesh
# scripts left main (they live on mesh/k6-p13-enum).
CODE="${K6_CODE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SSH=(ssh -F /home/nick/.ssh/mesh.config)
PIDS="$ROOT/mesh_pids"
GAUGE="$ROOT/k6_gpu_out"
PICKLE="$ROOT/k6_tasks.pkl"
STOP="$ROOT/mesh_stop"
NUKA_MNT=/home/nick/mnt/maxplus_p13
NUKA_PY=/home/nick/.venvs/rocm72/bin/python
ORIN_MNT=/home/nick/mnt/maxplus_p13
JF_MNT=/home/nick/mnt/maxplus_p13
WORKERS=(v100 nuka orin a380 cpu44 dash)
export PYTHONPATH="$CODE${PYTHONPATH:+:$PYTHONPATH}"

usage() {
  echo "usage: k6-mesh start|stop|status|man|dash [NAME]  (NAME: all ${WORKERS[*]})" >&2
  echo "       stop is always soft (finish current orbit). stop --hard NAME for kill." >&2
  exit 2
}

is_worker() {
  local w
  for w in "${WORKERS[@]}"; do [[ "$1" == "$w" ]] && return 0; done
  [[ "$1" == "all" ]] && return 0
  return 1
}

pidfile() { echo "$PIDS/$1.pid"; }

k6py() {
  python3 -c "import k6_control as c; $1"
}

alive() {
  local pf name=$1 pid
  pf=$(pidfile "$name")
  pid=""
  [[ -f "$pf" ]] && pid=$(cat "$pf" 2>/dev/null || true)
  case "$name" in
    v100)
      pgrep -f 'run_kgauged.py 6 2' >/dev/null
      ;;
    cpu44)
      [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null
      ;;
    dash)
      [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null
      ;;
    nuka)
      "${SSH[@]}" nuka "pgrep -f '[r]un_kgauged.py' >/dev/null" 2>/dev/null
      ;;
    orin)
      "${SSH[@]}" orin "pgrep -f '[r]un_kgauged.py' >/dev/null" 2>/dev/null
      ;;
    a380)
      "${SSH[@]}" jellyfin "pgrep -f '[r]un_kgauged.py' >/dev/null" 2>/dev/null
      ;;
    *) return 1 ;;
  esac
}

rsync_code() {
  local host=$1
  rsync -e "ssh -F /home/nick/.ssh/mesh.config" -a \
    --include='*.py' --include='*.sh' --include='*.cpp' --include='*.so' --include='*.1' --exclude='*' \
    "$CODE"/ "$host:/tmp/e1work/"
}

ensure_sshfs() {
  local host=$1 mnt=$2
  # jellyfin (and any new box) may lack mesh.config; sshfs needs Host soulkiller.
  if ! "${SSH[@]}" "$host" "test -f /home/nick/.ssh/mesh.config" 2>/dev/null; then
    scp -F /home/nick/.ssh/mesh.config /home/nick/.ssh/mesh.config \
      "$host:/home/nick/.ssh/mesh.config" >/dev/null
    "${SSH[@]}" "$host" "chmod 600 /home/nick/.ssh/mesh.config"
  fi
  "${SSH[@]}" "$host" "bash -s" <<EOF
set -e
mkdir -p $mnt /tmp/e1work
if ! mountpoint -q $mnt; then
  sshfs -F /home/nick/.ssh/mesh.config -o reconnect,ServerAliveInterval=15,idmap=user \
    soulkiller:/mnt/storage/e1work/maxplus_p13 $mnt
fi
mountpoint -q $mnt
touch $mnt/.wtest && rm -f $mnt/.wtest
EOF
}

prepare_start() {
  mkdir -p "$STOP" "$PIDS"
  k6py "c.prepare_start('$1')"
}

# --- start ---
start_v100() {
  alive v100 && { echo "v100 already running (pgrep run_kgauged.py 6 2)"; return 0; }
  prepare_start v100
  mkdir -p "$PIDS" "$GAUGE"
  [[ -f "$PICKLE" ]] || { echo "missing $PICKLE" >&2; return 1; }
  (
    cd "$CODE"
    export PALEY_P=13 GAUGE_OUT="$GAUGE" LOAD_TASKS="$PICKLE"
    export RESUME_JSON="$ROOT/k6_resume.json" K6_ROOT="$ROOT" K6_STOP_DIR="$STOP"
    export GPU_WORKERS=2 GPU_MEM_FRAC=0.70 GEN_CAP=40000000
    export ENUM_SHARD_MOD=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
    export K6_HOST=v100 K6_BACKEND=cuda PYTHONPATH="$CODE${PYTHONPATH:+:$PYTHONPATH}"
    echo "===== mesh V100 $(date -Is) =====" >> "$ROOT/enum_p13_k6.log"
    nohup python3 -u run_kgauged.py 6 2 >> "$ROOT/enum_p13_k6.log" 2>&1 &
    echo $! > "$(pidfile v100)"
  )
  echo "v100 started pid $(cat "$(pidfile v100)")"
}

remote_pid_from() {
  printf '%s\n' "$1" | sed -n 's/^K6PID=//p' | tail -1
}

start_nuka() {
  alive nuka && { echo "nuka already running"; return 0; }
  prepare_start nuka
  rsync_code nuka
  ensure_sshfs nuka "$NUKA_MNT"
  # soulkiller pickle is numpy 2; nuka venv numpy cannot load it.
  out=$("${SSH[@]}" nuka "bash -s" <<EOF
set -e
export PALEY_P=13 GAUGE_OUT=$NUKA_MNT/k6_gpu_out
unset LOAD_TASKS || true
export RESUME_JSON=$NUKA_MNT/k6_resume.json K6_ROOT=$NUKA_MNT K6_STOP_DIR=$NUKA_MNT/mesh_stop
export GPU_WORKERS=1 GPU_MEM_FRAC=0.70 GEN_CAP=40000000
export ENUM_SHARD_MOD=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export K6_HOST=nuka K6_BACKEND=cuda PYTHONPATH=/tmp/e1work
cd /tmp/e1work
echo "===== mesh nuka \$(date -Is) =====" >> $NUKA_MNT/enum_p13_k6_nuka.log
nohup $NUKA_PY -u run_kgauged.py 6 1 >> $NUKA_MNT/enum_p13_k6_nuka.log 2>&1 &
echo K6PID=\$!
EOF
)
  nuka_pid=$(remote_pid_from "$out")
  echo "$nuka_pid" > "$(pidfile nuka)"
  echo "nuka started pid $nuka_pid"
}

start_orin() {
  alive orin && { echo "orin already running"; return 0; }
  prepare_start orin
  rsync_code orin
  ensure_sshfs orin "$ORIN_MNT"
  local load=""
  [[ -f "$PICKLE" ]] && load="export LOAD_TASKS=$ORIN_MNT/k6_tasks.pkl"
  out=$("${SSH[@]}" orin "bash -s" <<EOF
set -e
export PALEY_P=13 GAUGE_OUT=$ORIN_MNT/k6_gpu_out
$load
export RESUME_JSON=$ORIN_MNT/k6_resume.json K6_ROOT=$ORIN_MNT K6_STOP_DIR=$ORIN_MNT/mesh_stop
export GPU_WORKERS=1 GPU_MEM_FRAC=0.45 GEN_CAP=8000000
export ENUM_SHARD_MOD=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export K6_HOST=orin K6_BACKEND=cuda PYTHONPATH=/tmp/e1work
cd /tmp/e1work
echo "===== mesh orin \$(date -Is) =====" >> $ORIN_MNT/enum_p13_k6_orin.log
nohup python3 -u run_kgauged.py 6 1 >> $ORIN_MNT/enum_p13_k6_orin.log 2>&1 &
echo K6PID=\$!
EOF
)
  opid=$(remote_pid_from "$out")
  echo "$opid" > "$(pidfile orin)"
  echo "orin started pid $opid"
}

start_a380() {
  alive a380 && { echo "a380 already running"; return 0; }
  prepare_start a380
  rsync_code jellyfin
  ensure_sshfs jellyfin "$JF_MNT"
  "${SSH[@]}" jellyfin "bash -s" <<'EOF'
set -e
cd /tmp/e1work
# shellcheck disable=SC1091
source /opt/intel/oneapi/setvars.sh >/dev/null
# lib prefix so the .so does not shadow gpu_gen_sycl.py
if [[ ! -f libgpu_gen_sycl.so ]] || [[ gpu_gen_sycl.cpp -nt libgpu_gen_sycl.so ]]; then
  icpx -fsycl -O3 -shared -fPIC -o libgpu_gen_sycl.so gpu_gen_sycl.cpp -lsycl
fi
rm -f gpu_gen_sycl.so
EOF
  local load=""
  [[ -f "$PICKLE" ]] && load="export LOAD_TASKS=$JF_MNT/k6_tasks.pkl"
  out=$("${SSH[@]}" jellyfin "bash -s" <<EOF
set -e
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
# intelpython's libsycl is older than compiler/2025.2 — force compiler .so + system python.
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2025.2/lib:\${LD_LIBRARY_PATH:-}
export LIBSYCL_SO=/opt/intel/oneapi/compiler/2025.2/lib/libsycl.so.8
export PALEY_P=13 GAUGE_OUT=$JF_MNT/k6_gpu_out
$load
export RESUME_JSON=$JF_MNT/k6_resume.json K6_ROOT=$JF_MNT K6_STOP_DIR=$JF_MNT/mesh_stop
export GPU_WORKERS=1 GEN_CAP=20000000
export ENUM_SHARD_MOD=1 PYTHONUNBUFFERED=1
export NUMBA_NUM_THREADS=14 OMP_NUM_THREADS=14
export NUMBA_THREADING_LAYER=workqueue
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export K6_HOST=a380 K6_BACKEND=sycl GPU_GEN_SYCL_SO=/tmp/e1work/libgpu_gen_sycl.so
export ONEAPI_DEVICE_SELECTOR=level_zero:gpu
export ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE
export PYTHONPATH=/tmp/e1work
cd /tmp/e1work
echo "===== mesh a380 \$(date -Is) =====" >> $JF_MNT/enum_p13_k6_a380.log
nohup /usr/bin/python3 -u run_kgauged.py 6 1 >> $JF_MNT/enum_p13_k6_a380.log 2>&1 &
echo K6PID=\$!
EOF
)
  jpid=$(remote_pid_from "$out")
  echo "$jpid" > "$(pidfile a380)"
  echo "a380 started pid $jpid"
}

start_cpu44() {
  alive cpu44 && { echo "cpu44 already running pid $(cat "$(pidfile cpu44)")"; return 0; }
  prepare_start cpu44
  mkdir -p "$PIDS" "$GAUGE"
  [[ -f "$PICKLE" ]] || { echo "missing $PICKLE" >&2; return 1; }
  (
    cd "$CODE"
    export PALEY_P=13 GAUGE_OUT="$GAUGE" LOAD_TASKS="$PICKLE"
    export RESUME_JSON="$ROOT/k6_resume.json" K6_ROOT="$ROOT" K6_STOP_DIR="$STOP"
    # 44 independent orbit processes. Inner OpenMP/numba stays 1: NumpyTester
    # is serial numpy, so one process × 44 OMP threads pegs ~1 core.
    export GPU_WORKERS=44 GEN_CAP=8000000
    export ENUM_SHARD_MOD=1 PYTHONUNBUFFERED=1
    export K6_HOST=cpu44 K6_BACKEND=cpu PYTHONPATH="$CODE${PYTHONPATH:+:$PYTHONPATH}"
    export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
    export NUMBA_THREADING_LAYER=workqueue
    export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    export CUDA_VISIBLE_DEVICES=
    echo "===== mesh cpu44 $(date -Is) =====" >> "$ROOT/enum_p13_k6_cpu44.log"
    nohup python3 -u run_kgauged.py 6 44 >> "$ROOT/enum_p13_k6_cpu44.log" 2>&1 &
    echo $! > "$(pidfile cpu44)"
  )
  echo "cpu44 started pid $(cat "$(pidfile cpu44)") (44 orbit processes, OMP=1, no CUDA)"
}

start_dash() {
  if alive dash; then
    echo "dash already running pid $(cat "$(pidfile dash)")"
    return 0
  fi
  mkdir -p "$PIDS"
  (
    cd "$CODE"
    export K6_ROOT="$ROOT" K6_STOP_DIR="$STOP" K6_CODE="$CODE" PYTHONPATH="$CODE${PYTHONPATH:+:$PYTHONPATH}"
    nohup python3 -u k6_dashboard.py --serve 8765 --loop 8 >> "$ROOT/k6_dashboard.log" 2>&1 &
    echo $! > "$(pidfile dash)"
  )
  echo "dashboard http://127.0.0.1:8765/k6_dashboard.html  pid $(cat "$(pidfile dash)")"
  echo "ssh -L 8765:127.0.0.1:8765 nick@192.168.1.113"
}

mark_if_needed() {
  mkdir -p "$PIDS" "$GAUGE"
  n0=$(ls "$GAUGE"/orb*.npy 2>/dev/null | wc -l)
  python3 - <<PY
import json, time, os
p="$ROOT/k6_mesh_mark.json"
if not os.path.isfile(p):
    json.dump({"t0": time.time(), "n_done": int("$n0")}, open(p,"w"))
    print("mark n_done=$n0")
else:
    print("mark exists")
PY
}

start_one() {
  mark_if_needed
  case "$1" in
    v100) start_v100 ;;
    nuka) start_nuka ;;
    orin) start_orin ;;
    a380) start_a380 ;;
    cpu44) start_cpu44 ;;
    dash) start_dash ;;
    all)
      start_dash
      start_v100
      start_nuka
      start_orin
      start_a380
      start_cpu44
      ;;
    *) echo "unknown worker $1" >&2; return 2 ;;
  esac
}

# --- stop (always soft) ---
soft_stop_legacy_starve() {
  k6py "import json; print(json.dumps(c.starve_unclaimed('$GAUGE')))"
}

drop_dummy_locks() {
  k6py "print('dropped starve locks', c.drop_starve_locks('$GAUGE'))"
}

wait_pid_gone() {
  local name=$1 timeout=${2:-3600} i=0
  while alive "$name"; do
    if (( i % 15 == 0 )); then
      echo "waiting $name to finish current orbit(s) (${i}s)"
    fi
    sleep 1
    i=$((i + 1))
    if (( i >= timeout )); then
      echo "timeout waiting for $name (still finishing an orbit; not SIGKILL)" >&2
      return 1
    fi
  done
  rm -f "$(pidfile "$name")"
  echo "$name stopped (soft)"
}

stop_one() {
  local name=$1 hard=${2:-0}
  mkdir -p "$STOP"
  cd "$CODE"
  if [[ "$hard" == 1 ]]; then
    echo "HARD stop $name (SIGTERM then SIGKILL) — not the default" >&2
    hard_kill "$name"
    return
  fi
  if [[ "$name" == "all" ]]; then
    k6py "c.request_stop('ALL')"
    # Old workers ignore the flag; starve new orbs so current locks still finish.
    # New workers also honour the flag and will not mkdir.
    soft_stop_legacy_starve
    for w in v100 nuka orin a380 cpu44; do
      wait_pid_gone "$w" 3600 || true
    done
    drop_dummy_locks
    echo "k6-mesh all soft-stopped"
    return
  fi
  k6py "c.request_stop('$name')"
  if [[ "$name" == "dash" ]]; then
    if [[ -f "$(pidfile dash)" ]]; then
      kill -TERM "$(cat "$(pidfile dash)")" 2>/dev/null || true
      sleep 1
      kill -KILL "$(cat "$(pidfile dash)")" 2>/dev/null || true
      rm -f "$(pidfile dash)"
    fi
    echo "dash stopped"
    return
  fi
  # Per-worker stop is flag-only so the rest of the mesh keeps claiming orbs.
  # Dummy starve is mesh-wide and would halt everyone — only used for `stop all`.
  wait_pid_gone "$name" 3600 || true
}

kill_tree_by_pid() {
  local pid=$1 sig=${2:-TERM} kids k
  [[ -n "$pid" ]] || return 0
  kill -0 "$pid" 2>/dev/null || return 0
  kids=$(pgrep -P "$pid" 2>/dev/null || true)
  for k in $kids; do
    kill_tree_by_pid "$k" "$sig"
  done
  kill -"$sig" "$pid" 2>/dev/null || true
}

hard_kill() {
  local name=$1 pf pid
  pf=$(pidfile "$name")
  pid=$(cat "$pf" 2>/dev/null || true)
  case "$name" in
    v100)
      pgrep -f 'run_kgauged.py 6 2' | xargs -r kill -TERM 2>/dev/null || true
      sleep 2
      pgrep -f 'run_kgauged.py 6 2' | xargs -r kill -KILL 2>/dev/null || true
      ;;
    cpu44)
      # Kill the Pool parent and descendants by PID tree (never pkill -f:
      # that matches this wrapper). Also reap leftover 6 1 / 6 44 workers.
      kill_tree_by_pid "$pid" TERM
      sleep 2
      kill_tree_by_pid "$pid" KILL
      for wp in $(pgrep -f '[r]un_kgauged.py 6 (1|44)' || true); do
        # v100 is '6 2'; only reap cpu44 argv shapes.
        kill -KILL "$wp" 2>/dev/null || true
      done
      ;;
    dash)
      [[ -n "$pid" ]] && kill -TERM "$pid" 2>/dev/null || true
      kill -KILL "$pid" 2>/dev/null || true
      ;;
    nuka)
      "${SSH[@]}" nuka "pkill -TERM -f 'run_kgauged.py' 2>/dev/null; sleep 2; pkill -KILL -f 'run_kgauged.py' 2>/dev/null" || true
      ;;
    orin)
      "${SSH[@]}" orin "pkill -TERM -f 'run_kgauged.py' 2>/dev/null; sleep 2; pkill -KILL -f 'run_kgauged.py' 2>/dev/null" || true
      ;;
    a380)
      "${SSH[@]}" jellyfin "pkill -TERM -f 'run_kgauged.py' 2>/dev/null; sleep 2; pkill -KILL -f 'run_kgauged.py' 2>/dev/null" || true
      ;;
  esac
  rm -f "$pf"
}

status() {
  python3 "$CODE/k6_dashboard.py" || true
  echo "--- pids ---"
  for w in "${WORKERS[@]}"; do
    if alive "$w"; then
      echo "$w RUNNING pid $(cat "$(pidfile "$w")" 2>/dev/null || echo '?')"
    else
      echo "$w stopped"
    fi
  done
}

cmd="${1:-}"
shift || true
hard=0
if [[ "${1:-}" == "--hard" ]]; then
  hard=1
  shift
fi
name="${1:-all}"
case "$cmd" in
  start)
    is_worker "$name" || usage
    start_one "$name"
    ;;
  stop)
    is_worker "$name" || usage
    stop_one "$name" "$hard"
    ;;
  status) status ;;
  dash) python3 "$CODE/k6_dashboard.py" --serve 8765 ;;
  man)
    if command -v man >/dev/null && [[ -f "$CODE/k6-mesh.1" ]]; then
      man -l "$CODE/k6-mesh.1"
    else
      cat "$CODE/k6-mesh.1"
    fi
    ;;
  *) usage ;;
esac
