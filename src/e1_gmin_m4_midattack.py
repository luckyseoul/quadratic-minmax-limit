#!/usr/bin/env python3
"""
P0 attack: prove max_{|κ|=1}|m4| ≤ M_mid(p)=(p-2)/(2p(p+1)) for primes p≥5.

Strategy:
  1) Algebra already: M_mid ≤ L_abs < T_abs (Prop 15.70).
  2) GPU: E[dot^4], Tr(G²) exact for Max+; fit closed forms in p.
  3) Multi-worker: moduli nullity + Tr pin at p=5,7; compare selected g_min to -M_mid.
  4) Boolean kurtosis: E[dot^4] ≤ Wick - δ candidate ⇒ g_min > T; check if enough for mid/L.
  5) Atomic evidence write.

F17: GPU for dense dots; no serial sparse T. File script only.
Writes evidence/e1_gmin_m4_midattack.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def T_abs(p: int) -> float:
    return (p - 2) / (p * (2.0 * p - 1))


def algebra() -> dict:
    rows = {}
    for p in range(5, 60, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        mid, La, Ta = M_mid(p), L_abs(p), T_abs(p)
        rows[str(p)] = {
            "M_mid": mid,
            "L_abs": La,
            "T_abs": Ta,
            "mid_le_L": mid <= La + 1e-15,
            "L_lt_T": La < Ta,
            "mid_over_L": mid / La,
        }
    return {
        "proved": all(r["mid_le_L"] and r["L_lt_T"] for r in rows.values()),
        "by_p": rows,
        "note": "max|m4|≤M_mid ⇒ g_min≥-M_mid≥L>T ⇒ bi-tight empty",
    }


def _gpu_dot_moments(p: int) -> dict:
    """E[dot²], E[dot^4], Tr(G²) for Max+ via GPU (pairwise dots in batches)."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import require_gpu, gpu_info
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    cp = require_gpu()
    Y = np.ascontiguousarray(load_maxplus_array(p, mmap=True), dtype=np.float64)
    N, n = Y.shape
    # 2-design formula (Max+-free)
    p_conf = int(round(np.sqrt(n - 1)))
    E_dot2_formula = n + n * (n - 1) / (p_conf * p_conf)

    # GPU: Y is N×n of ±1. Gram G = Y @ Y.T is N×N of dots (including self).
    # For large N=11452, G is 11452² × 8 ≈ 1.05 GB — fits V100 16GB.
    t0 = time.perf_counter()
    Yg = cp.asarray(Y)
    # dots matrix (N,N)
    G = Yg @ Yg.T  # float64
    # E[dot²] over all ordered pairs including diagonal
    # mean of G_ij²
    E_dot2 = float(cp.mean(G * G).get())
    E_dot4 = float(cp.mean(G * G * G * G).get())
    # Also exclude diagonal for "distinct pairs" variant
    # mask diagonal
    eye = cp.eye(N, dtype=cp.float64)
    off = 1.0 - eye
    n_off = N * (N - 1)
    E_dot2_off = float(cp.sum((G * G) * off).get() / n_off)
    E_dot4_off = float(cp.sum((G * G * G * G) * off).get() / n_off)
    wall = time.perf_counter() - t0

    # Tr(G_feat²) via K_ab = ((dot_ab)² - n)/2 on unordered pairs a<b of Max+ vectors
    # Tr(G_edge²) identity: (1/4)(E[dot^4]-2n E[dot²]+n²) uses expectation over Max+ pairs
    # Use all ordered pairs including self for the shipped identity check
    Tr_from_mom = 0.25 * (E_dot4 - 2 * n * E_dot2 + n * n)
    Tr_from_mom_off = 0.25 * (E_dot4_off - 2 * n * E_dot2_off + n * n)

    # Wick comparison: Σ = I + C/p, Gaussian 4th moment
    C = paley_conference_prime_power(p).astype(np.float64)
    Sigma = np.eye(n) + C / p_conf
    # For two independent N(0,Σ) vectors u,v: dot=u·v
    # Actually the identity uses Max+ vectors as rows — dots of design points
    # Wick for design is approximate via Σ

    # Closed form candidates for E_dot4
    candidates = {
        "3*(E2)^2": 3 * E_dot2**2,  # Gaussian-like for mean0 — wrong scale
        "n^2 + 6*n*(n-1)/p2 + ...": None,
        "2design_E2": E_dot2_formula,
    }

    # Fit simple polynomials in n,p
    # Record fractions if close
    from fractions import Fraction

    def near_frac(x, lim=10**6):
        try:
            return str(Fraction(x).limit_denominator(lim))
        except Exception:
            return None

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "wall_s": float(wall),
        "gpu": gpu_info(),
        "E_dot2_allpairs": E_dot2,
        "E_dot4_allpairs": E_dot4,
        "E_dot2_offdiag": E_dot2_off,
        "E_dot4_offdiag": E_dot4_off,
        "E_dot2_2design_formula": E_dot2_formula,
        "E_dot2_matches_2design": abs(E_dot2 - E_dot2_formula) < 1e-6
        or abs(E_dot2_off - E_dot2_formula) < 1e-6,
        "TrG2_from_allpairs": Tr_from_mom,
        "TrG2_from_offdiag": Tr_from_mom_off,
        "E_dot4_frac": near_frac(E_dot4),
        "E_dot4_off_frac": near_frac(E_dot4_off),
        "TrG2_frac": near_frac(Tr_from_mom_off),
        "M_mid": M_mid(p),
        "L_abs": L_abs(p),
        # boolean vs wick kurtosis ratio
        "kurtosis_ratio_E4_over_3E2sq": E_dot4_off / (3 * E_dot2_off**2 + 1e-30),
    }


def _worker_moduli_gmin_vs_mid(p: int) -> dict:
    """Run shipped cbound/moduli path if available; else load evidence."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    mid = M_mid(p)
    La = L_abs(p)
    # Prefer live recompute via e1_gmin_cbound if main/export exists
    try:
        # Import module and call main pieces
        import e1_gmin_cbound as cb

        # Try common entry points
        if hasattr(cb, "run_p"):
            r = cb.run_p(p)
            gmin = r.get("selected_g_min_true_Tr") or r.get("g_min_at_c_true")
            return {
                "p": p,
                "via": "e1_gmin_cbound.run_p",
                "g_min": gmin,
                "ge_neg_mid": gmin is not None and gmin >= -mid - 1e-12,
                "ge_L": gmin is not None and gmin >= -La - 1e-12,
                "raw": {k: r[k] for k in r if k in (
                    "selected_g_min_true_Tr", "g_min_at_c_true", "true_root_ge_L",
                    "Tr_G2_true", "Tr_G2_wick", "nullity",
                )},
            }
    except Exception as e:
        err = str(e)

    path = ROOT / "evidence" / "e1_gmin_cbound.json"
    if path.is_file():
        data = json.loads(path.read_text())
        r = data.get("results", {}).get(str(p), {})
        gmin = r.get("selected_g_min_true_Tr")
        return {
            "p": p,
            "via": "evidence",
            "g_min": gmin,
            "ge_neg_mid": gmin is not None and gmin >= -mid - 1e-12,
            "ge_L": r.get("true_root_ge_L"),
            "Tr_true": r.get("Tr_G2_true"),
            "Tr_wick": r.get("Tr_G2_wick"),
            "wick_selected": r.get("selected_g_min_wick_Tr"),
        }
    return {"p": p, "error": "no moduli"}


def _worker_boolean_dot4_ub(p: int) -> dict:
    """
    Boolean kurtosis bound attempt (Max+-free structure only):
    For z_i = y_i y'_i on coordinates, z_i=±1, sum z = dot.
    E[dot^4] expansion via power sums of z under design constraints.
    """
    _blas1()
    # Pure algebra from E[yy^T]=I+C/p only — Wick is the Gaussian upper envelope
    # Boolean: for mean-zero pairwise, Maclaurin / hypercontractivity
    # |E[prod of 4 coords]| ≤ ... 
    # Record Wick E_dot4 as UB candidate
    n = p * p + 1
    E2 = n + n * (n - 1) / (p * p)
    # For Gaussian vectors with cov Σ=I+C/p, E[(u·v)^4] for u,v iid N(0,Σ):
    # This is more subtle. Shipped wick uses 3||Σ||_F^4 + 6 Tr(Σ^4) for a different form.
    # Here: naive Gaussian on the scalar dot with variance E2 would give 3 E2²
    wick_scalar = 3 * E2 * E2
    return {
        "p": p,
        "E_dot2_2design": E2,
        "wick_scalar_3E2sq": wick_scalar,
        "note": "scalar Gaussian 3E2² is usually an upper envelope for boolean; verify vs GPU",
    }


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    print(f"midattack: W={W} gpu={snap['gpu'].get('available')}", flush=True)

    out: dict = {
        "workers": W,
        "compute": snap,
        "parts": {"algebra": algebra()},
    }
    print(f"  algebra proved={out['parts']['algebra']['proved']}", flush=True)

    # Boolean UB formulas (tiny, parent)
    out["parts"]["boolean_ub"] = {
        str(p): _worker_boolean_dot4_ub(p) for p in (5, 7, 11, 13)
    }

    # GPU dot moments for p=5 and p=7 sequentially (one CUDA context; each is dense N×N)
    if prefer_gpu_for("dense_matmul"):
        for p in (5, 7):
            print(f"  GPU dot moments p={p}...", flush=True)
            out["parts"][f"gpu_dots_p{p}"] = _gpu_dot_moments(p)
            r = out["parts"][f"gpu_dots_p{p}"]
            print(
                f"    E4_off={r['E_dot4_offdiag']:.6f} E2={r['E_dot2_offdiag']:.6f} "
                f"kurt={r['kurtosis_ratio_E4_over_3E2sq']:.4f} wall={r['wall_s']:.2f}s",
                flush=True,
            )
    else:
        out["parts"]["gpu_dots"] = {"skipped": "GPU not preferred/available"}

    # Multi-worker moduli compare (p=5,7 and evidence)
    print("  moduli g_min vs mid (process pool)...", flush=True)
    with ProcessPoolExecutor(max_workers=min(W, 4)) as ex:
        futs = {ex.submit(_worker_moduli_gmin_vs_mid, p): p for p in (5, 7)}
        for fut in as_completed(futs):
            r = fut.result()
            out["parts"][f"moduli_p{r['p']}"] = r
            print(f"    moduli p={r['p']}: g_min={r.get('g_min')} ge_mid={r.get('ge_neg_mid')}", flush=True)

    # Compare GPU E4 to scalar wick
    synth = {}
    for p in (5, 7):
        g = out["parts"].get(f"gpu_dots_p{p}", {})
        b = out["parts"]["boolean_ub"].get(str(p), {})
        if g and b:
            synth[str(p)] = {
                "E4_off": g.get("E_dot4_offdiag"),
                "wick_3E2sq": b.get("wick_scalar_3E2sq"),
                "E4_le_wick_scalar": (
                    g.get("E_dot4_offdiag") is not None
                    and g["E_dot4_offdiag"] <= b["wick_scalar_3E2sq"] + 1e-6
                ),
                "kurtosis_ratio": g.get("kurtosis_ratio_E4_over_3E2sq"),
                "g_min_census_ge_neg_mid": out["parts"].get(f"moduli_p{p}", {}).get("ge_neg_mid"),
            }
    out["synthesis"] = synth
    out["status"] = (
        "Algebra M_mid≤L proved. GPU E[dot^4] measured p=5,7. "
        "OPEN: Max+-free closed form / UB on E[dot^4] tight enough to force "
        "moduli pin g_min≥-M_mid for all primes p≥5 (Wick-scalar alone too loose for mid)."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_midattack.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
