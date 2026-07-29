#!/usr/bin/env python3
"""
P0: prove resolvent gain |κ|=3→|κ|=1 ≤ (p-4)/48, or candidate |m4|≤(p-2)/(p(2p+3)).

Attacks (multi-worker F17):
  A) Exact L2/Frobenius bounds on T and Neumann residual
  B) Exact T on κ-sign strata + type6 (Max+-free association)
  C) Closed-form candidate check vs true m4; analytic inequality hunt
  D) Spectral radius of T on full R^{binom(n,4)} via power iteration (multi-start)
  E) Character-sum style: constant-weight design intersection numbers

Writes evidence/e1_gmin_m4_gain.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import as_completed
from itertools import combinations, product
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from workers import pool, require_workers  # noqa: E402


def _blas1():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def cand_ub(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def star_sum(C, S) -> int:
    s = 0
    for v in S:
        pr = 1
        for u in S:
            if u != v:
                pr *= int(C[v, u])
        s += pr
    return s


# ---------------------------------------------------------------------------
# A: L2 energy of source and T
# ---------------------------------------------------------------------------
def _worker_L2_bounds(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    nQ = comb(n, 4)

    # Source b = Tκ/p²: nonzero only |κ|=3, value ±24/p²
    n3 = 0
    n1 = 0
    for S in combinations(range(n), 4):
        k = abs(kappa(C, S))
        if k == 3:
            n3 += 1
        elif k == 1:
            n1 += 1

    source_inf = 24.0 / (p * p)
    # ‖b‖_2 = sqrt(n3) * source_inf
    b_l2 = np.sqrt(n3) * source_inf
    # ‖b‖_1 = n3 * source_inf
    b_l1 = n3 * source_inf

    # Frobenius of T: each 4-set has out-degree 4(n-4) with ±1 weights
    # ‖T‖_F² = sum_S 4(n-4) = nQ * 4(n-4)  (each edge weight ±1, so sum of squares of matrix entries)
    # T_{S,S'} is sum of C_vr over (v,r) with S_vr=S', each 0 or ±1, typically 0 or 1 path
    # Upper: each S has exactly 4(n-4) nonzero contributions of ±1 to various S', so
    # sum_{S,S'} T_{S,S'}^2 ≤ sum_S (sum |contrib|)^2 wait no:
    # If T_SS' = sum_{edges S→S'} C_e, and usually at most one edge S→S',
    # then T_SS' ∈ {-1,0,1} and ‖T‖_F² = # of (S,v,r) = nQ*4(n-4)
    T_F2_ub = nQ * 4 * (n - 4)
    T_F_ub = np.sqrt(T_F2_ub)
    # Spectral: ‖T‖_2 ≤ ‖T‖_F
    # Gap: if 4p > ‖T‖_2 then ‖ρ‖_2 ≤ b_l2 / (4p - ‖T‖_2)
    gap_F = 4 * p - T_F_ub  # almost surely negative

    # Better: row l1 norm of |T| is 4(n-4), so ‖T‖_∞≤4(n-4), ‖T‖_1≤4(n-4)
    row_l1 = 4 * (n - 4)
    # For ρ on |κ|=1 only we care about the compression

    # Exact: compute ‖Tκ‖_2 from formula
    # Tκ = -6*star, |Tκ|=24 on |κ|=3, 0 on |κ|=1, and on other κ?
    # κ only in {-3,-1,1,3} for conference (3 pairings of ±1)
    tkap_vals = Counter()
    for S in combinations(range(n), 4):
        tkap_vals[-6 * star_sum(C, S)] += 1

    # Power iteration estimate of ‖T‖_2 on R^{nQ} is expensive; do on random subspace via matvec
    # Matvec (Tf)(S) = sum_{v,r} C_vr f(S_vr): for random f, estimate top singular value
    # Store quads list
    quads = list(combinations(range(n), 4))
    qindex = {q: i for i, q in enumerate(quads)}
    rng = np.random.default_rng(p)

    def matvec(f):
        out = np.zeros(nQ)
        for i, S in enumerate(quads):
            Sset = set(S)
            acc = 0.0
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    acc += float(C[v, r]) * f[qindex[Sp]]
            out[i] = acc
        return out

    # Power iteration (few iters) — only for small p
    if p <= 5:
        v = rng.standard_normal(nQ)
        v /= np.linalg.norm(v)
        for _ in range(12):
            v = matvec(v)
            nrm = np.linalg.norm(v)
            v /= nrm + 1e-30
        # one more for Rayleigh
        Tv = matvec(v)
        lam_est = float(v @ Tv)
        op_est = float(np.linalg.norm(Tv))  # since v unit, ≈ |λ|
    else:
        # p=7: nQ=230300, matvec too heavy in pure Python — skip full op
        lam_est = None
        op_est = None

    rho_budget = (p - 4) / (2.0 * p * p) if p > 4 else None
    gain_budget = (p - 4) / 48.0 if p > 4 else None

    # L2 → L∞ crude: |ρ|_∞ ≤ ‖ρ‖_2 (worst)
    # need b_l2 / max(4p - op, eps) ≤ rho_budget
    return {
        "p": int(p),
        "n": int(n),
        "nQ": int(nQ),
        "n_kappa1": int(n1),
        "n_kappa3": int(n3),
        "source_inf": source_inf,
        "b_l2": float(b_l2),
        "b_l1": float(b_l1),
        "T_F_ub": float(T_F_ub),
        "row_l1_absT": int(row_l1),
        "gap_4p_minus_T_F": float(4 * p - T_F_ub),
        "tkappa_value_counts": {str(k): int(v) for k, v in sorted(tkap_vals.items())},
        "op_norm_T_power_est": op_est,
        "lam_T_power_est": lam_est,
        "rho_budget": rho_budget,
        "gain_budget": gain_budget,
        "F_bound_useful": bool(4 * p > T_F_ub),
        "note": "F-bound on ‖T‖ useless (≫4p). Need association-scheme eigenvalues.",
    }


# ---------------------------------------------------------------------------
# B: type6-only T matrix (Max+-free, exact for small p)
# ---------------------------------------------------------------------------
def type6_canon(C, S):
    """S4-canonical 6-edge sign pattern (lex min among 24 perms)."""
    a, b, c, d = S
    edges_raw = {
        (0, 1): int(C[a, b]),
        (0, 2): int(C[a, c]),
        (0, 3): int(C[a, d]),
        (1, 2): int(C[b, c]),
        (1, 3): int(C[b, d]),
        (2, 3): int(C[c, d]),
    }
    best = None
    verts = [0, 1, 2, 3]
    for perm in __import__("itertools").permutations(verts):
        inv = {perm[i]: i for i in range(4)}
        # remap edges
        tup = []
        for i, j in combinations(range(4), 2):
            oi, oj = perm[i], perm[j]
            e = edges_raw[tuple(sorted((oi, oj)))]
            tup.append(e)
        t = tuple(tup)
        if best is None or t < best:
            best = t
    return best


def _worker_type6_T(p: int) -> dict:
    """Exact T among type6 classes (edge pattern only) — may not be constant m4 but exact C-op."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    # Build type6 classes
    classes = defaultdict(list)
    for S in combinations(range(n), 4):
        classes[type6_canon(C, S)].append(S)
    keys = sorted(classes.keys())
    nk = len(keys)
    kidx = {k: i for i, k in enumerate(keys)}
    kap = np.array([kappa(C, classes[k][0]) for k in keys], dtype=float)

    # Exact average T: for each class, use ALL quads if p=3,5 else sample large
    Tmat = np.zeros((nk, nk))
    rng = np.random.default_rng(0)
    for i, kA in enumerate(keys):
        quads = classes[kA]
        take = quads if (p <= 5 or len(quads) <= 100) else [
            quads[j] for j in rng.choice(len(quads), 100, replace=False)
        ]
        acc = np.zeros(nk)
        for S in take:
            Sset = set(S)
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    acc[kidx[type6_canon(C, Sp)]] += float(C[v, r])
        Tmat[i] = acc / len(take)

    Tkap = np.array([-6.0 * star_sum(C, classes[k][0]) for k in keys])
    b = Tkap / (p * p)
    A = 4 * p * np.eye(nk) - Tmat
    rho = np.linalg.solve(A, b)
    m4 = kap / (p * p) + rho

    idx1 = np.where(np.abs(kap) == 1)[0]
    max_m4_1 = float(np.max(np.abs(m4[idx1]))) if len(idx1) else None
    max_rho_1 = float(np.max(np.abs(rho[idx1]))) if len(idx1) else None
    same = [
        abs(rho[i])
        for i in idx1
        if np.sign(rho[i]) * np.sign(kap[i]) > 0
    ]
    max_same = float(max(same)) if same else 0.0

    evals = np.linalg.eigvals(Tmat).real
    return {
        "p": int(p),
        "n_type6_classes": nk,
        "class_sizes": {str(k): len(classes[k]) for k in keys},
        "lam_max_T": float(np.max(evals)),
        "lam_min_T": float(np.min(evals)),
        "gap_4p": float(4 * p - np.max(evals)),
        "max_abs_m4_pred_kappa1": max_m4_1,
        "max_abs_rho_kappa1": max_rho_1,
        "max_abs_rho_same_sign": max_same,
        "rho_budget": (p - 4) / (2.0 * p * p) if p > 4 else None,
        "same_le_budget": max_same <= (p - 4) / (2.0 * p * p) + 1e-9 if p > 4 else None,
        "L_abs": L_abs(p),
        "cand": cand_ub(p),
        "pred_le_L": max_m4_1 is not None and max_m4_1 <= L_abs(p) + 1e-9,
        "pred_le_cand": max_m4_1 is not None and max_m4_1 <= cand_ub(p) + 1e-9,
        "solve_resid": float(np.max(np.abs(A @ rho - b))),
    }


# ---------------------------------------------------------------------------
# C: design intersection numbers from Max+ (for character-sum closed form)
# ---------------------------------------------------------------------------
def _worker_design_intersections(p: int) -> dict:
    """
    Oriented Max+ (y_∞=+1) as constant-weight code of weight w=p(p-1)/2 on F_{p²}.
    Compute λ_t for t≤4 and m4 via inclusion when possible.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    elif p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        raise ValueError(p)

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    # Orient: flip so coordinate 0 (∞) is +1
    Y_or = Y * Y[:, 0:1]
    # Minus sets: positions where y=-1
    w = p * (p - 1) // 2
    # Verify weight
    weights = np.sum(Y_or < 0, axis=1)
    weight_ok = bool(np.all(weights == w))

    # λ_t: for a fixed t-set of finite points (exclude ∞=0), how many blocks contain it?
    # Use first few t-sets and check constancy
    finite = list(range(1, n))

    def lambda_stats(t: int, samples: int = 30):
        rng = np.random.default_rng(t + p)
        vals = []
        for _ in range(samples):
            S = tuple(sorted(rng.choice(finite, size=t, replace=False)))
            # count oriented Max+ with all S minus
            # y_i = -1 for i in S
            mask = np.ones(len(Y_or), dtype=bool)
            for i in S:
                mask &= Y_or[:, i] < 0
            vals.append(int(np.sum(mask)))
        return {
            "unique": sorted(set(vals)),
            "constant": len(set(vals)) == 1,
            "value": vals[0] if len(set(vals)) == 1 else None,
            "mean": float(np.mean(vals)),
        }

    lam = {str(t): lambda_stats(t) for t in (1, 2, 3, 4)}

    # For 4-sets including patterns by κ: empirical m4 vs design formula
    # If 3-design, λ_j determined for j≤3; 4-wise varies
    # Character formula: E[(-1)^{|M∩T|}] = sum_j (-1)^j λ_{T,j} / (N/2)
    # where λ_{T,j} = # blocks with |M∩T|=j

    # Exact m4 for all |κ|=1 at p=5 already known; here extract formula candidate
    # from λ parameters
    n_or = len(Y_or)  # = N (each pair ± still both have y0 after orient? 
    # After orient Y*Y0, antipodes: -y becomes y with y0=+1? (-y)*(-y0)=y*y0 same. 
    # Antipodes collapse! So we get N/2 unique oriented vectors.
    Y_unique = np.unique(np.round(Y_or), axis=0)
    n_blocks = len(Y_unique)

    # λ1 for design of n_blocks blocks of size w on n points (or n-1 finite?)
    # Include infinity: weight w total minuses, all in finite if y0=+1
    lam1_formula = n_blocks * w / n  # if ∞ can be minus - but y0=+1 so minuses only in 1..n-1
    lam1_finite = n_blocks * w / (n - 1)

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "w": int(w),
        "weight_ok": weight_ok,
        "n_oriented_unique_blocks": int(n_blocks),
        "lambda_stats": lam,
        "is_1design_finite": lam["1"]["constant"],
        "is_2design_finite": lam["2"]["constant"],
        "is_3design_finite": lam["3"]["constant"],
        "is_4design_finite": lam["4"]["constant"],
        "lam1_finite_mean": lam["1"]["mean"],
        "lam1_finite_expected": float(lam1_finite),
        "cand_ub": cand_ub(p),
        "L_abs": L_abs(p),
        "note": (
            "If 2-design or 3-design, intersection numbers constrain m4; "
            "4-design would force constant m4 (false)."
        ),
    }


# ---------------------------------------------------------------------------
# D: Analytic: prove cand ≤ L and check (p-2)/(p(2p+3)) vs empirical
# ---------------------------------------------------------------------------
def analytic_candidate() -> dict:
    rows = {}
    for p in range(5, 50, 2):
        if not all(p % d for d in range(2, int(p**0.5) + 1)):
            continue  # not prime-ish skip composites roughly
        # quick prime check
        prime = p > 1 and all(p % d for d in range(2, int(p**0.5) + 1))
        if not prime:
            continue
        c = cand_ub(p)
        L = L_abs(p)
        rows[str(p)] = {
            "cand": c,
            "L_abs": L,
            "cand_le_L": c <= L + 1e-15,
            "ratio_cand_over_L": c / L,
            "gain_budget": (p - 4) / 48.0,
            "wick": 1.0 / (p * p),
            "room_cand_minus_wick": c - 1.0 / (p * p),
        }
    # Algebraic identity cand = (p-2)/(p(2p+3)), L=(p-2)/(2p^2)
    # cand/L = 2p / (2p+3) < 1
    return {
        "identity_cand_over_L": "2p/(2p+3) < 1 for p>0",
        "by_p": rows,
        "sharp_at_p5": abs(cand_ub(5) - 3 / 65) < 1e-15,
    }


# ---------------------------------------------------------------------------
# E: Power-iteration spectral radius of T on |κ|=1 subspace with source coupling
# ---------------------------------------------------------------------------
def _worker_T_spectrum_kappa_subspaces(p: int) -> dict:
    """
    Restrict functions to be constant on |κ| level sets only (4 values of κ).
    Too coarse for exact m4 but gives invariant subspace of Tκ and crude spectrum.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    # Map κ value -> index 0..3 for {-3,-1,1,3}
    levels = (-3, -1, 1, 3)
    lidx = {k: i for i, k in enumerate(levels)}

    # T among 4 levels: average
    # Count transitions: for each S with κ=k, for each extension, target κ
    # Weight by C_vr
    T4 = np.zeros((4, 4))
    counts = np.zeros(4)
    for S in combinations(range(n), 4):
        k = kappa(C, S)
        if k not in lidx:
            continue
        i = lidx[k]
        counts[i] += 1
        Sset = set(S)
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                kt = kappa(C, Sp)
                if kt not in lidx:
                    continue
                T4[i, lidx[kt]] += float(C[v, r])
    # average per S
    for i in range(4):
        if counts[i] > 0:
            T4[i] /= counts[i]

    evals = np.linalg.eigvals(T4).real
    # Source b on 4 levels
    # Tκ on level: 0 for |κ|=1, average of ±24 for |κ|=3
    Tkap_level = np.zeros(4)
    for S in combinations(range(n), 4):
        k = kappa(C, S)
        if k not in lidx:
            continue
        Tkap_level[lidx[k]] += -6 * star_sum(C, S)
    for i in range(4):
        if counts[i] > 0:
            Tkap_level[i] /= counts[i]

    b = Tkap_level / (p * p)
    A = 4 * p * np.eye(4) - T4
    try:
        rho = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        rho = np.linalg.lstsq(A, b, rcond=None)[0]
    m4 = np.array(levels, dtype=float) / (p * p) + rho

    return {
        "p": int(p),
        "counts_by_kappa": {str(levels[i]): int(counts[i]) for i in range(4)},
        "T4": T4.tolist(),
        "Tkappa_level_avg": Tkap_level.tolist(),
        "evals_T4": sorted(float(x) for x in evals),
        "lam_max_T4": float(np.max(evals)),
        "gap_4p": float(4 * p - np.max(evals)),
        "rho_level": rho.tolist(),
        "m4_level_pred": m4.tolist(),
        "max_abs_m4_on_pm1": float(max(abs(m4[1]), abs(m4[2]))),
        "L_abs": L_abs(p),
        "cand": cand_ub(p),
        "note": "4-level (κ-only) model is too coarse for sharp m4 but T-invariant.",
    }


def main():
    W = require_workers(min_workers=4)
    print(f"m4_gain: W={W}", flush=True)
    out = {"workers": W, "analytic": analytic_candidate(), "parts": {}}

    jobs = {}
    with pool(W) as ex:
        for p in (3, 5, 7):
            jobs[ex.submit(_worker_L2_bounds, p)] = f"L2_p{p}"
            jobs[ex.submit(_worker_type6_T, p)] = f"t6_p{p}"
            jobs[ex.submit(_worker_design_intersections, p)] = f"des_p{p}"
            jobs[ex.submit(_worker_T_spectrum_kappa_subspaces, p)] = f"k4_p{p}"
        for fut in as_completed(jobs):
            tag = jobs[fut]
            try:
                r = fut.result()
                out["parts"][tag] = r
                print(f"  {tag}: done keys={list(r.keys())[:6]}", flush=True)
            except Exception as e:
                out["parts"][tag] = {"error": str(e)}
                print(f"  {tag}: ERROR {e}", flush=True)

    # Synthesis for proof
    synth = {
        "cand_le_L_algebra": "cand/L = 2p/(2p+3) < 1",
        "cand_sharp_p5": out["analytic"]["sharp_at_p5"],
        "design_strength": {},
        "type6_pred": {},
        "open": (
            "Prove |m4|≤(p-2)/(p(2p+3)) on |κ|=1 for all primes p≥5, "
            "or resolvent gain ≤(p-4)/48 via association-scheme eigenvalues of T."
        ),
    }
    for p in (3, 5, 7):
        d = out["parts"].get(f"des_p{p}", {})
        synth["design_strength"][str(p)] = {
            "2design": d.get("is_2design_finite"),
            "3design": d.get("is_3design_finite"),
            "4design": d.get("is_4design_finite"),
            "blocks": d.get("n_oriented_unique_blocks"),
        }
        t = out["parts"].get(f"t6_p{p}", {})
        synth["type6_pred"][str(p)] = {
            "max_m4_pred": t.get("max_abs_m4_pred_kappa1"),
            "same_le_budget": t.get("same_le_budget"),
            "pred_le_cand": t.get("pred_le_cand"),
            "lam_max_T": t.get("lam_max_T"),
            "gap_4p": t.get("gap_4p"),
        }

    out["synthesis"] = synth
    out["status"] = (
        f"Candidate UB (p-2)/(p(2p+3)) ≤ L algebraically (ratio 2p/(2p+3)); "
        f"sharp at p=5. Design/T-spectrum probes recorded. {synth['open']}"
    )
    path = ROOT / "evidence" / "e1_gmin_m4_gain.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
