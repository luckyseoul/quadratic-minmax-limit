#!/usr/bin/env python3
"""
Probe: random perfect-matching flips of Paley C_26 (ρ=1 family, Φ=65).

At n=10, 144/945 perfect-matching flips achieve the exact optimum m_10=13
(absolute gap 2 from Paley Φ=15). Here we sample matchings at n=26 and
score exact Φ (numba half-cube) to see whether a constant/small gap reappears.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import (
    is_conference_matrix,
    halfspace_boolean_vector,
    paley_conference_prime_power,
    rho_from_boolean_evec,
    spherical_half_bound,
)

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False


if HAS_NUMBA:

    @njit(cache=True)
    def exact_phi_numba(A: np.ndarray) -> float:
        """
        Exact max |Q| over x in {±1}^n with x[0]=+1, via Gray code.

        Cost O(n 2^{n-1}): when coordinate i flips, Q changes by
        -2 x_i (Ax)_i and Ax -= 2 x_i A[:, i].
        """
        n = A.shape[0]
        x = np.ones(n, dtype=np.float64)
        Ax = A @ x
        q2 = float(x @ Ax)  # = x^T A x = 2Q
        best = abs(q2) * 0.5
        gray = 0
        npat = 1 << (n - 1)
        for i in range(1, npat):
            ng = i ^ (i >> 1)
            bit = 0
            diff = gray ^ ng
            while (diff >> bit) & 1 == 0:
                bit += 1
            gray = ng
            # flip coordinate bit+1 (x[0] fixed)
            j = bit + 1
            xj = x[j]
            # delta (x^T A x) = -4 xj (Ax)_j ? 
            # new x' = x - 2 xj e_j; x'^T A x' = x^T A x - 4 xj (Ax)_j
            # because diag A = 0.
            delta = -4.0 * xj * Ax[j]
            q2 = q2 + delta
            # update Ax: Ax' = Ax - 2 xj A[:, j]
            for r in range(n):
                Ax[r] = Ax[r] - 2.0 * xj * A[r, j]
            x[j] = -xj
            aq = abs(q2) * 0.5
            if aq > best:
                best = aq
        return best

else:

    def exact_phi_numba(A: np.ndarray) -> float:
        n = A.shape[0]
        x = np.ones(n, dtype=np.float64)
        Ax = A @ x
        q2 = float(x @ Ax)
        best = abs(q2) * 0.5
        gray = 0
        npat = 1 << (n - 1)
        for i in range(1, npat):
            ng = i ^ (i >> 1)
            diff = gray ^ ng
            bit = diff.bit_length() - 1
            gray = ng
            j = bit + 1
            xj = x[j]
            q2 += -4.0 * xj * Ax[j]
            Ax -= 2.0 * xj * A[:, j]
            x[j] = -xj
            aq = abs(q2) * 0.5
            if aq > best:
                best = aq
        return best


def random_perfect_matching(n: int, rng: np.random.Generator):
    verts = list(range(n))
    rng.shuffle(verts)
    return [(verts[2 * i], verts[2 * i + 1]) for i in range(n // 2)]


def flip_matching(C: np.ndarray, matching):
    A = C.copy()
    for i, j in matching:
        A[i, j] *= -1.0
        A[j, i] *= -1.0
    return A


def _worker(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    C_list, seeds = args
    C = np.array(C_list, dtype=np.float64)
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        m = random_perfect_matching(C.shape[0], rng)
        A = flip_matching(C, m)
        t0 = time.time()
        ph = float(exact_phi_numba(A))
        rows.append({"seed": int(seed), "phi": ph, "time_s": time.time() - t0})
    return rows


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    n_samples = int(os.environ.get("N_SAMPLES", "86"))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    p = 5
    n = p * p + 1
    print(f"Building Paley C_{n} (p={p}) ...", flush=True)
    C = paley_conference_prime_power(p)
    assert is_conference_matrix(C)
    x = halfspace_boolean_vector(p)
    assert abs(rho_from_boolean_evec(C, x) - 1.0) < 1e-8
    phi_paley = 0.5 * n * p  # ρ=1
    assert abs(phi_paley - spherical_half_bound(n)) < 1e-9

    # Warmup numba on tiny
    if HAS_NUMBA:
        print("Warmup numba exact_phi on n=10 ...", flush=True)
        from minmax_quadratic import paley_conference_prime_power as pp

        C10 = pp(3)
        t0 = time.time()
        ph10 = exact_phi_numba(C10)
        print(f"  Φ(C_10)={ph10} (expect 15) t={time.time()-t0:.3f}s", flush=True)

    print(
        f"Sampling {n_samples} random perfect-matching flips; workers={workers} ...",
        flush=True,
    )
    seeds = list(range(50_000, 50_000 + n_samples))
    # one seed per worker job for load balance when n_samples≈workers
    jobs = []
    C_list = C.tolist()
    # batch seeds
    batch = max(1, (n_samples + workers - 1) // workers)
    for i in range(0, n_samples, batch):
        jobs.append((C_list, seeds[i : i + batch]))

    t_all = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        done = 0
        for fut in as_completed(futs):
            part = fut.result()
            rows.extend(part)
            done += len(part)
            best = min(r["phi"] for r in rows)
            print(
                f"  progress {done}/{n_samples}  running_best_Φ={best}  "
                f"last_t={part[-1]['time_s']:.2f}s",
                flush=True,
            )

    phis = [r["phi"] for r in rows]
    report = {
        "n": n,
        "p": p,
        "Phi_Paley": phi_paley,
        "has_numba": HAS_NUMBA,
        "n_samples": n_samples,
        "workers": workers,
        "phi_min": min(phis),
        "phi_max": max(phis),
        "phi_mean": float(np.mean(phis)),
        "phi_median": float(np.median(phis)),
        "n_undercut_paley": sum(1 for p_ in phis if p_ < phi_paley - 1e-9),
        "n_at_or_below_paley": sum(1 for p_ in phis if p_ <= phi_paley + 1e-9),
        "gap_min": phi_paley - min(phis),
        "total_time_s": time.time() - t_all,
        "samples": rows,
        "note": (
            "Random perfect-matching flips of Paley-26. At n=10 this construction "
            "hits m_10. Here exact Φ via half-cube (numba if available)."
        ),
    }
    out = out_dir / "n26_matching_probe.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"wrote {out}", flush=True)
    print(
        f"RESULT minΦ={report['phi_min']} Paley={phi_paley} "
        f"gap={report['gap_min']} undercuts={report['n_undercut_paley']}/{n_samples}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
