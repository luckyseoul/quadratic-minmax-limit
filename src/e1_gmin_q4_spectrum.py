#!/usr/bin/env python3
"""
Exact Q / Q4 spectrum on zero-diag ∩ V_+ with multiplicities and formula fit.
Also: relate image Gdisj Rayleigh to (p,d,N) closed form.

Multi-worker. Writes evidence/e1_gmin_q4_spectrum.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    raise ValueError(p)


def _A_from_vec(v: np.ndarray, d: int) -> np.ndarray:
    A = np.zeros((d, d))
    k = 0
    for i in range(d):
        for j in range(i, d):
            A[i, j] = A[j, i] = v[k]
            k += 1
    return A


def _null_zero_diag(Vp: np.ndarray) -> np.ndarray:
    n, d = Vp.shape
    dimS = d * (d + 1) // 2
    M = np.zeros((n, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        M[:, k] = np.diag(Vp @ _A_from_vec(e, d) @ Vp.T)
    _u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    return vh[rank:].T


def _worker_spectrum(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    null = _null_zero_diag(Vp)

    # ONB of zero-diag ∩ V_+ matrices
    Bs = []
    for t in range(null.shape[1]):
        A = _A_from_vec(null[:, t], d)
        B = Vp @ A @ Vp.T
        for Bb in Bs:
            B = B - np.sum(B * Bb) * Bb
        nrm = np.linalg.norm(B)
        if nrm > 1e-12:
            Bs.append(B / nrm)
    m = len(Bs)

    # Full Q Gram: M_st = sum_y (y Bs y)(y Bt y)
    # Vectorized: build Qvals N x m
    Qvals = np.stack([np.einsum("ai,ij,aj->a", Y, B, Y) for B in Bs], axis=1)
    Mq = Qvals.T @ Qvals  # m x m
    evals_Q = np.linalg.eigvalsh(0.5 * (Mq + Mq.T))

    # Q4 = Q - 6N on diagonal (for unit B): M4 = Mq - 6N I
    M4 = Mq - 6.0 * N * np.eye(m)
    evals_Q4 = np.linalg.eigvalsh(0.5 * (M4 + M4.T))

    def pack(evals):
        # round to rational-ish
        rounded = [round(float(e), 8) for e in evals]
        ctr = Counter(rounded)
        items = []
        for val, mult in sorted(ctr.items(), reverse=True):
            f = Fraction(val / N).limit_denominator(500)
            items.append(
                {
                    "eigenvalue": val,
                    "eig_over_N": val / N,
                    "frac_over_N": f"{f.numerator}/{f.denominator}",
                    "frac_err": abs(val / N - float(f)),
                    "mult": mult,
                }
            )
        return items

    q_pack = pack(evals_Q)
    q4_pack = pack(evals_Q4)

    # Image Gdisj rays = Q4/(2N) for unit-B eigenspaces
    rays = []
    for item in q4_pack:
        ray = item["eigenvalue"] / (2.0 * N)
        f = Fraction(ray).limit_denominator(500)
        rays.append(
            {
                "ray": ray,
                "frac": f"{f.numerator}/{f.denominator}",
                "frac_err": abs(ray - float(f)),
                "mult": item["mult"],
                "Q4": item["eigenvalue"],
            }
        )

    # Formula candidates for λ_max(Q)/N and λ_max(Q4)/N and ray_max
    lamQ_N = float(evals_Q[-1]) / N
    lamQ4_N = float(evals_Q4[-1]) / N
    ray_max = lamQ4_N / 2.0

    formulas = {
        "lamQ_N": lamQ_N,
        "lamQ4_N": lamQ4_N,
        "ray_max": ray_max,
        "16": 16.0,
        "16*(d-2)/d": 16 * (d - 2) / d,
        "2*(d+1)/(d-1)": 2 * (d + 1) / (d - 1) if d > 1 else None,
        # ray candidates
        "ray_5": 5.0,
        "ray_49_over_d": 49 / d,
        "ray_(p+2)": p + 2,
        "ray_2p/(p-1)": 2 * p / (p - 1) if p > 1 else None,
        "ray_(p^2-1)/(2p)": (p * p - 1) / (2 * p),
        "ray_4d/(d+2)": 4 * d / (d + 2),
        "ray_2*(d-1)/d": 2 * (d - 1) / d,
        "ray_(d+1)/4": (d + 1) / 4,
        "ray_d/(d-2)": d / (d - 2) if d > 2 else None,
        "ray_(d-1)/4": (d - 1) / 4,
        "ray_4*(d-1)/(d+1)": 4 * (d - 1) / (d + 1),
        # from p=5: 49/13; 49=(2d-3)^2 / something? 2*13-3=23 no. 7^2. 7=(p+2)
        "ray_(p+2)^2/d": (p + 2) ** 2 / d,
        "ray_(p+2)^2/(2d)": (p + 2) ** 2 / (2 * d),
        "ray_p^2/d": p * p / d,
        "ray_(2p-1)^2/(2d)": (2 * p - 1) ** 2 / (2 * d),
        "ray_(p^2+1)/(2p)": (p * p + 1) / (2 * p),  # n/(2p)
        "ray_n/(2p)": n / (2 * p),
        "ray_(p+1)/2 * something": (p + 1) / 2,
    }

    # Exact match check for ray
    ray_matches = {}
    for name, val in formulas.items():
        if name.startswith("ray_") and val is not None and name != "ray_max":
            ray_matches[name] = abs(ray_max - val)

    best_ray = min(ray_matches.items(), key=lambda kv: kv[1]) if ray_matches else None

    # dim checks
    # binom(d-1,2) = d(d-1)/2 - (d-1) wait binom(d-1,2)=(d-1)(d-2)/2
    expected_dim = (d - 1) * (d - 2) // 2

    return {
        "p": int(p),
        "n": int(n),
        "d": int(d),
        "N": int(N),
        "dim_space": int(m),
        "expected_binom_d1_2": int(expected_dim),
        "Q_spectrum": q_pack,
        "Q4_spectrum": q4_pack,
        "Gdisj_image_rays": rays,
        "lamQ_over_N": lamQ_N,
        "lamQ4_over_N": lamQ4_N,
        "ray_max": ray_max,
        "thr_ray_5": 5.0,
        "ray_below_5": bool(ray_max <= 5.0 + 1e-9),
        "best_ray_formula": {"name": best_ray[0], "abs_err": best_ray[1]} if best_ray else None,
        "all_ray_formula_errs": {k: v for k, v in sorted(ray_matches.items(), key=lambda x: x[1])[:15]},
        "ratio_to_16N": lamQ_N / 16.0,
        "bound_16N_ok": bool(lamQ_N <= 16.0 + 1e-9),
    }


def _worker_p7_ratio_exact() -> dict:
    """High-precision ratio at p=7 via many power starts; guess fraction."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    p = 7
    Y = _load_Y_safe7()
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(7)
    best_Q = -1.0
    for trial in range(30):
        A = rng.standard_normal((d, d))
        A = 0.5 * (A + A.T)
        A -= np.trace(A) / d * np.eye(d)
        A /= np.linalg.norm(A) + 1e-30
        for _ in range(100):
            q = np.einsum("bi,ij,bj->b", U, A, U)
            Ph = (U * q[:, None]).T @ U
            Ph = 0.5 * (Ph + Ph.T)
            Ph -= np.trace(Ph) / d * np.eye(d)
            A = Ph / (np.linalg.norm(Ph) + 1e-30)
        B = Vp @ A @ Vp.T
        bf2 = float(np.sum(B * B))
        ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
        Q = float(np.sum(ytBy**2)) / bf2
        if Q > best_Q:
            best_Q = Q
    qn = best_Q / N
    q4n = qn - 6.0
    ray = q4n / 2.0
    # fraction guesses
    guesses = []
    for den in range(1, 400):
        for num in range(1, 20 * den):
            if abs(ray - num / den) < 1e-9:
                guesses.append(f"{num}/{den}")
            if abs(qn - num / den) < 1e-9:
                guesses.append(f"qn={num}/{den}")
    f_ray = Fraction(ray).limit_denominator(1000)
    f_qn = Fraction(qn).limit_denominator(1000)
    f_r16 = Fraction(best_Q / (16 * N)).limit_denominator(1000)
    return {
        "p": 7,
        "N": N,
        "d": d,
        "Q_max": best_Q,
        "qn": qn,
        "q4n": q4n,
        "ray_max": ray,
        "frac_ray": f"{f_ray.numerator}/{f_ray.denominator}",
        "frac_ray_err": abs(ray - float(f_ray)),
        "frac_qn": f"{f_qn.numerator}/{f_qn.denominator}",
        "frac_qn_err": abs(qn - float(f_qn)),
        "frac_ratio16": f"{f_r16.numerator}/{f_r16.denominator}",
        "frac_ratio16_err": abs(best_Q / (16 * N) - float(f_r16)),
        "bound_ok": best_Q <= 16 * N + 1e-3,
    }


def _load_Y_safe7():
    return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"q4_spectrum: workers={W}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {
            ex.submit(_worker_spectrum, 3): 3,
            ex.submit(_worker_spectrum, 5): 5,
            ex.submit(_worker_p7_ratio_exact): 7,
        }
        # fan out p=7 power more
        for s in range(W - 3):
            futs[ex.submit(_worker_p7_one, 10000 + s)] = ("p7s", s)

        p7_starts = []
        for fut in as_completed(futs):
            tag = futs[fut]
            r = fut.result()
            if isinstance(tag, int) and tag in (3, 5):
                rows.append(r)
                print(
                    f"  p={r['p']}: dim={r['dim_space']} exp={r['expected_binom_d1_2']} "
                    f"ray_max={r['ray_max']:.8f} best={r['best_ray_formula']} "
                    f"Q4spec={[ (x['frac_over_N'], x['mult']) for x in r['Q4_spectrum'] ]}",
                    flush=True,
                )
            elif tag == 7:
                rows.append(r)
                print(
                    f"  p=7 exact: ray={r['ray_max']:.10f}={r['frac_ray']} "
                    f"qn={r['qn']:.10f}={r['frac_qn']} r16={r['frac_ratio16']} "
                    f"err_ray={r['frac_ray_err']:.2e}",
                    flush=True,
                )
            else:
                p7_starts.append(r)

    if p7_starts:
        best = max(p7_starts, key=lambda r: r["Q"])
        print(
            f"  p=7 fanout best Q/N={best['Q']/best['N']:.10f} ray={best['ray']:.10f} "
            f"frac={best['frac_ray']}",
            flush=True,
        )
        rows.append({"p": 7, "fanout_best": best, "n_starts": len(p7_starts)})

    out = {
        "results": rows,
        "status": (
            "Q4/Gdisj image spectrum with multiplicities p=3,5; p=7 high-prec ratio. "
            "Need closed form ray_max(p)≤5 for proof."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_q4_spectrum.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


def _worker_p7_one(seed: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    C = paley_conference_prime_power(7).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(90):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    bf2 = float(np.sum(B * B))
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2)) / bf2
    qn = Q / N
    q4n = qn - 6.0
    ray = q4n / 2.0
    f = Fraction(ray).limit_denominator(2000)
    return {
        "seed": seed,
        "Q": Q,
        "N": N,
        "qn": qn,
        "q4n": q4n,
        "ray": ray,
        "frac_ray": f"{f.numerator}/{f.denominator}",
        "frac_err": abs(ray - float(f)),
    }


if __name__ == "__main__":
    main()
