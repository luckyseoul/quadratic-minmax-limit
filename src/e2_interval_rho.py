#!/usr/bin/env python3
"""
E(2) campaign: cube/sphere ratio ρ of Paley conference via interval boolean vectors.

For Paley C of order n=q+1 (q prime ≡1 mod 4), the interval signing
  x_∞ = +1,  x_t = +1 iff t ∈ {0,1,...,⌊q/2⌋}  (field elements 0..q-1 as vertices 1..q)
gives a constructive lower bound ρ_int ≤ ρ(C).

If ρ_int → 1 then ρ(C) → 1 (since ρ ≤ 1), closing E(2) with ρ_*=1.

This script certifies ρ_int at many q by exact matrix-vector (not SA), and records
the asymptotic (1-ρ_int)√n for the conjecture ρ = 1 - Θ(n^{-1/2}).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    is_conference_matrix,
    paley_conference_matrix,
    form_Q,
)

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_1_mod_4(limit: int) -> List[int]:
    out = []
    for n in range(5, limit + 1):
        if n % 4 == 1 and _is_prime(n):
            out.append(n)
    return out


def interval_boolean_vector(q: int) -> np.ndarray:
    """
    Interval signing on F_q ∪ {∞}: vertices 0=∞, 1..q ↔ field 0..q-1.
    x_∞=+1; x_t = +1 for t=0..⌊q/2⌋ else -1.
    """
    n = q + 1
    x = np.ones(n, dtype=np.float64)
    half = q // 2  # floor(q/2); include 0..half → (half+1) pluses among field
    for t in range(q):
        x[t + 1] = 1.0 if t <= half else -1.0
    return x


def rho_from_vector(C: np.ndarray, x: np.ndarray) -> float:
    n = C.shape[0]
    # ρ = |x^T C x| / (n √(n-1))
    return abs(float(x @ C @ x)) / (n * math.sqrt(n - 1))


def _job(q: int) -> Dict:
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    t0 = time.time()
    C = paley_conference_matrix(q)
    assert is_conference_matrix(C)
    n = C.shape[0]
    x = interval_boolean_vector(q)
    # also try shifted intervals: x_t = +1 iff (t-s) mod q ≤ q/2
    best = 0.0
    best_s = 0
    # sample shifts (full for small q, stride for large)
    step = 1 if q <= 500 else max(1, q // 200)
    for s in range(0, q, step):
        xs = np.ones(n, dtype=np.float64)
        half = q // 2
        for t in range(q):
            u = (t - s) % q
            xs[t + 1] = 1.0 if u <= half else -1.0
        r = rho_from_vector(C, xs)
        if r > best:
            best = r
            best_s = s
    return {
        "q": q,
        "n": n,
        "rho_int": best,
        "best_shift": best_s,
        "one_minus_rho_times_sqrt_n": (1.0 - best) * math.sqrt(n),
        "Phi_LB": 0.5 * n * math.sqrt(n - 1) * best,
        "Phi_UB_spherical": 0.5 * n * math.sqrt(n - 1),
        "time_s": time.time() - t0,
        "shift_step": step,
    }


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    q_limit = int(os.environ.get("Q_LIMIT", "2500"))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    qs = primes_1_mod_4(q_limit)
    # also include a few larger for asymptotics
    extra = [3001, 4001, 5003, 6007, 8009, 10007, 12007, 15013]
    for q in extra:
        if q % 4 == 1 and _is_prime(q) and q not in qs:
            qs.append(q)
    qs = sorted(set(qs))

    print(f"E2 interval ρ: {len(qs)} primes, workers={workers}", flush=True)
    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_job, q): q for q in qs}
        done = 0
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            done += 1
            if done % 20 == 0 or row["q"] > 2000:
                print(
                    f"  [{done}/{len(qs)}] q={row['q']} ρ_int={row['rho_int']:.6f} "
                    f"(1-ρ)√n={row['one_minus_rho_times_sqrt_n']:.4f}",
                    flush=True,
                )

    rows.sort(key=lambda r: r["q"])
    best = max(rows, key=lambda r: r["rho_int"])
    gap_scaled = [r["one_minus_rho_times_sqrt_n"] for r in rows]
    report = {
        "campaign": "E2_interval_rho",
        "workers": workers,
        "q_limit": q_limit,
        "n_primes": len(rows),
        "best_rho_int": best["rho_int"],
        "best_q": best["q"],
        "best_n": best["n"],
        "min_rho_int": min(r["rho_int"] for r in rows),
        "rho_int_at_largest_n": rows[-1]["rho_int"],
        "largest_n": rows[-1]["n"],
        "one_minus_rho_sqrt_n_median": float(np.median(gap_scaled)),
        "one_minus_rho_sqrt_n_mean": float(np.mean(gap_scaled)),
        "total_time_s": time.time() - t0,
        "rows": rows,
        "e2_status": {
            "proved_rho_to_1": False,
            "constructive_limsup_rho_ge": best["rho_int"],
            "note": (
                "ρ_int is a lower bound on ρ(C). Growth of ρ_int toward 1 supports "
                "E(2) with ρ_*=1 but is not a full proof for all Paley (need analytic "
                "Weil-type bound uniform in q). Limit existence still requires E(1)."
            ),
            "limit_still_open": True,
        },
    }
    out = out_dir / "e2_interval_rho.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}", flush=True)
    print(
        f"best ρ_int={best['rho_int']:.6f} at n={best['n']}; "
        f"largest n={rows[-1]['n']} ρ_int={rows[-1]['rho_int']:.6f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
