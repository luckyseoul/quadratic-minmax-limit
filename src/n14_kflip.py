#!/usr/bin/env python3
"""Exact min Φ at Hamming distance k from Paley n=14 (q=13)."""
from __future__ import annotations
import json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from minmax_quadratic import paley_conference_matrix, is_conference_matrix
from e1_paley_gap_campaign import exact_phi_gray

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

C = paley_conference_matrix(13)
assert is_conference_matrix(C)
EDGES = [(i, j) for i in range(14) for j in range(i + 1, 14)]
C_LIST = C.tolist()


def _shard(args):
    for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_k] = "1"
    C_list, combos = args
    A0 = np.array(C_list, dtype=np.float64)
    best = 1e18
    for combo in combos:
        A = A0.copy()
        for t in combo:
            i, j = EDGES[t]
            A[i, j] *= -1
            A[j, i] *= -1
        v = exact_phi_gray(A)
        if v < best:
            best = v
    return best


def main():
    workers = int(os.environ.get("WORKERS", 86))
    k_max = int(os.environ.get("K_MAX", "4"))
    Phi = exact_phi_gray(C)
    out = {"n": 14, "Phi": Phi, "k_table": []}
    print(f"Paley n=14 Φ={Phi}", flush=True)
    for k in range(0, k_max + 1):
        if k == 0:
            row = {"k": 0, "min_phi": Phi, "n_combos": 1}
        else:
            combos = list(combinations(range(len(EDGES)), k))
            shard = max(1, (len(combos) + workers * 4 - 1) // (workers * 4))
            jobs = []
            for s in range(0, len(combos), shard):
                jobs.append((C_LIST, combos[s : s + shard]))
            t0 = time.time()
            best = Phi
            with ProcessPoolExecutor(max_workers=workers) as ex:
                for b in ex.map(_shard, jobs):
                    best = min(best, b)
            row = {
                "k": k,
                "min_phi": best,
                "n_combos": len(combos),
                "time_s": time.time() - t0,
                "undercuts": best < Phi - 1e-9,
            }
        out["k_table"].append(row)
        print(row, flush=True)
    path = Path(os.environ.get("OUT", str(ROOT / "evidence"))) / "n14_kflip.json"
    path.write_text(json.dumps(out, indent=2, default=lambda o: bool(o) if type(o).__name__=="bool_" else float(o) if hasattr(o, "item") else str(o)))
    print("wrote", path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
