#!/usr/bin/env python3
"""
Prove hypothesis H: ray ≤ (p+2)^2/d  (⇔ λ_cycle ≤ 3+(p+2)^2/d ≤ 8).

Approach:
  1) Closed form for λ2(K) on Max+ via association-scheme / harmonic analysis
  2) Direct bound λ_cycle ≤ 3+(p+2)^2/d using V_+ row structure of Vp
  3) Verify exact formula candidates against p=3,5,7

K_{yy'} = ((y·y')^2 - n)/(2N),  λ_max(K)=max(n/2, λ_cycle),
λ_cycle = 3 + ray_max on cycle space.

Multi-worker. Writes evidence/e1_gmin_H_proof.json
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
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _worker_scheme_params(p: int) -> dict:
    """Association scheme parameters of Max+ (dot-product classes) and K eigs."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    D = Y @ Y.T
    # Unique dots
    dots = sorted({round(float(x), 8) for x in np.unique(D)})
    # Degrees from first row
    row0 = D[0]
    deg = {}
    for v in dots:
        deg[v] = int(np.sum(np.abs(row0 - v) < 1e-6))

    # K matrix eigs via FF^T route (stable)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if N * len(edges) > 5e8:
        # p=7: use power for top eigs of K only
        return _worker_p7_K_eigs(p)

    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    # K = F F.T / N
    # For large N use chunked or eig of smaller if N small
    if N <= 300:
        K = (F @ F.T) / N
        evals = np.linalg.eigvalsh(K)
    else:
        # randomized SVD top k
        return _worker_p7_K_eigs(p)

    evals = np.sort(evals)[::-1]
    evals_pos = evals[evals > 1e-8]
    ctr = Counter(np.round(evals_pos, 10))
    eigs = []
    for val, mult in sorted(ctr.items(), reverse=True):
        f = Fraction(val).limit_denominator(500)
        eigs.append(
            {
                "value": float(val),
                "frac": f"{f.numerator}/{f.denominator}",
                "frac_err": abs(val - float(f)),
                "mult": int(mult),
            }
        )

    # Identify n/2 and λ_cycle
    nhalf = n / 2.0
    lam_cycle = None
    for e in eigs:
        if abs(e["value"] - nhalf) < 1e-6:
            continue
        lam_cycle = e["value"]
        break

    H = (p + 2) ** 2 / d
    H_plus_3 = 3.0 + H
    ray = (lam_cycle - 3.0) if lam_cycle is not None else None

    # Intersection numbers for class of antipodes and min class
    # p_ij^k rough: from first row only if distance-homogeneous
    # Check distance-homogeneous: all rows have same deg multiset
    degs_ok = True
    for r in range(min(N, 20)):
        ctr_r = Counter(np.round(D[r], 8))
        if dict(ctr_r) != {round(k, 8): v for k, v in deg.items()}:
            # compare as multiset of values
            pass
    # Actually check all rows have same degree sequence
    ref = tuple(sorted(np.round(D[0], 8).tolist()))
    for r in range(1, min(N, 50)):
        if tuple(sorted(np.round(D[r], 8).tolist())) != ref:
            degs_ok = False
            break

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "d": int(d),
        "dots": dots,
        "degrees": {str(k): int(v) for k, v in deg.items()},
        "distance_homogeneous_sample": degs_ok,
        "K_eigs": eigs,
        "n_half": nhalf,
        "lam_cycle": lam_cycle,
        "ray": ray,
        "H": H,
        "H_plus_3": H_plus_3,
        "H_holds": bool(ray is not None and ray <= H + 1e-8),
        "cycle_le_8": bool(lam_cycle is not None and lam_cycle <= 8 + 1e-8),
        "gap_ok": bool(lam_cycle is not None and lam_cycle <= nhalf + 1e-8),
    }


def _worker_p7_K_eigs(p: int = 7) -> dict:
    """Power iteration for top eigenvalues of K at p=7 without full matrix."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(7)
    C = paley_conference_prime_power(7).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    # K x = (1/N) F (F.T x) with F_ye = y_i y_j
    # F.T x = sum_y x_y f_e(y) = Be-like
    # (F.T x)_e for e=ij: sum_y x_y y_i y_j = T_ij
    # K x = (1/N) * (row action): (Kx)_y = (1/N) sum_e f_e(y) (F.T x)_e
    #      = (1/N) sum_{i<j} y_i y_j T_ij = (1/(2N)) (y^T T y)
    # with T zero-diag symmetric, T_ij = sum_{y'} x_{y'} y'_i y'_j

    def apply_K(x):
        T = (Y * x[:, None]).T @ Y  # n x n
        np.fill_diagonal(T, 0.0)  # numerical
        # (Kx)_y = (1/(2N)) y^T T y   but wait FF^T/N:
        # (FF^T x)_y = sum_e f_e(y) (F.T x)_e = sum_{i<j} y_i y_j T_ij
        # = (1/2) (y^T T y - sum_i T_ii y_i^2) = (1/2) y^T T y if diag T=0
        ytTy = np.einsum("ai,ij,aj->a", Y, T, Y)
        return ytTy / (2.0 * N)

    # Also verify: K 1 should be (n/2) 1
    ones = np.ones(N)
    K1 = apply_K(ones)
    nhalf_emp = float(np.mean(K1))

    # Power iteration in 1^⊥ for λ_cycle (largest in orthogonal complement)
    rng = np.random.default_rng(7)
    # Deflate n/2: work in 1^⊥
    def top_in_perp(n_iter=80):
        x = rng.standard_normal(N)
        x = x - x.mean()
        x /= np.linalg.norm(x) + 1e-30
        for _ in range(n_iter):
            y = apply_K(x)
            y = y - y.mean()
            nrm = np.linalg.norm(y)
            x = y / (nrm + 1e-30)
        # Rayleigh
        y = apply_K(x)
        return float(np.dot(x, y)), x

    # Several restarts
    best_lam = -1e300
    best_x = None
    for s in range(15):
        rng = np.random.default_rng(7 + s)
        x = rng.standard_normal(N)
        x = x - x.mean()
        x /= np.linalg.norm(x) + 1e-30
        for _ in range(100):
            y = apply_K(x)
            y = y - y.mean()
            x = y / (np.linalg.norm(y) + 1e-30)
        lam = float(np.dot(x, apply_K(x)))
        if lam > best_lam:
            best_lam = lam
            best_x = x.copy()

    # Second eigenvalue in perp via deflation of best_x
    lam2 = None
    for s in range(10):
        rng = np.random.default_rng(100 + s)
        x = rng.standard_normal(N)
        x = x - x.mean()
        x = x - np.dot(x, best_x) * best_x
        x /= np.linalg.norm(x) + 1e-30
        for _ in range(80):
            y = apply_K(x)
            y = y - y.mean()
            y = y - np.dot(y, best_x) * best_x
            x = y / (np.linalg.norm(y) + 1e-30)
        lam = float(np.dot(x, apply_K(x)))
        if lam2 is None or lam > lam2:
            lam2 = lam

    H = (7 + 2) ** 2 / d
    ray = best_lam - 3.0
    f = Fraction(best_lam).limit_denominator(2000)
    return {
        "p": 7,
        "N": int(N),
        "n": int(n),
        "d": int(d),
        "n_half": n / 2.0,
        "n_half_emp": nhalf_emp,
        "lam_cycle": best_lam,
        "lam_cycle_frac": f"{f.numerator}/{f.denominator}",
        "lam_cycle_frac_err": abs(best_lam - float(f)),
        "lam2_cycle_space": lam2,
        "ray": ray,
        "H": H,
        "H_plus_3": 3.0 + H,
        "H_holds": bool(ray <= H + 1e-6),
        "cycle_le_8": bool(best_lam <= 8 + 1e-6),
        "gap_ok": bool(best_lam <= n / 2.0 + 1e-6),
        "method": "power_K_matvec",
    }


def _worker_bound_attempt(p: int) -> dict:
    """
    Analytical bound attempts for λ_cycle.

    Key identity: for unit zero-diag B on V_+,
      E[(yBy)^2] = 6 + 2 ray  (with ‖B‖F=1)
      yBy = s^T A s, s=Vp^T y, ‖s‖^2=n, A unit Tr0, r_i^T A r_i=0.

    Attempt A: expand E[(s^T A s)^2] using only 2-design of {s} plus
    the zero-diag constraint to kill some terms.

    Attempt B: Gershgorin / row bounds on K.
    Attempt C: |y·y'| ≤ p+1 type bounds (likely fail).
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    rows = Vp  # n x d

    # 2-design moment of s = Vp^T y: E[s s^T] = (1/N) Y? 
    # s_y = Vp^T y, E[s s^T] = Vp^T E[yy^T] Vp = Vp^T (2 P_+) Vp = 2 I_d
    # So E[s s^T] = 2 I, ‖s‖^2 = n = 2d, consistent with E[‖s‖^2]=2d.

    # For unit A Tr0, sphere/2-design formula for E[(sAs)^2]:
    # If s were Gaussian N(0,2I): E[(sAs)^2] = 8 ‖A‖F² * something?
    # s = √2 g with g ~ isotropic; or E[(sAs)^2] = 2*4 Tr(A²) + ... Wick
    # Wick: E[(sAs)^2] = 2 Tr((Σ A)^2) + (Tr Σ A)^2 with Σ=2I
    # = 2 Tr(4 A²) + 0 = 8 ‖A‖F² for Tr A=0.
    # So E_Gauss = 8, Q_Gauss = 8N, which is Wick.

    # Boolean residual δ = E_Max+[(sAs)^2] - 8
    # H: E ≤ 6 + 2(p+2)^2/d = 6 + 4(p+2)^2/n
    # At unit A.

    # Bound attempt using ‖A‖_op under zero-diag constraints.
    # rows r_i: r_i^T A r_i = 0.
    # The quadratic forms r_i^T A r_i = 0 are n constraints.

    # CS: |s^T A s| = |sum_{αβ} A_αβ s_α s_β|
    # With A orthogonal to all r_i r_i^T in the sense of Frobenius after projection.

    # Compute max |r_i · r_j| = max |P_+_ij| = max |C_ij|/(2p) = 1/(2p) for i≠j
    # since |C_ij|=1.

    # Try: represent A in the basis of the Jacobi / harmonic polynomials.

    # Practical bound from coherence:
    # s = sum_i y_i r_i, s^T A s = sum_{i,j} y_i y_j (r_i^T A r_j)
    # = sum_{i≠j} y_i y_j H_ij where H_ij = r_i^T A r_j (H_ii=0 by constraint)
    # H is exactly B! So circular.

    # Gershgorin on K: for each y,
    # sum_{y'} |K_yy'| = sum_{y'} |((y·y')^2 - n)/(2N)|
    dots = Y @ Y[0]
    gersh = float(np.sum(np.abs((dots**2 - n) / (2 * N))))
    # This bounds all |λ| ≤ gersh — check if ≤ H+3 or 8
    K_diag = (n**2 - n) / (2 * N)

    # Better: the Seidel / conference bound on the second eigenvalue
    # of the graph with weights (y·y')^2.

    # Compute exact for maximizer E[(sAs)^2]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(p)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(70):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    B /= np.linalg.norm(B)
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    E_sAs2 = float(np.mean(ytBy**2))  # = Q/N for unit B

    H_qn = 6 + 2 * (p + 2) ** 2 / d
    return {
        "p": int(p),
        "N": int(N),
        "d": int(d),
        "K_diag": K_diag,
        "gershgorin_row_abs_sum": gersh,
        "gersh_le_8": bool(gersh <= 8 + 1e-9),
        "gersh_le_H3": bool(gersh <= 3 + (p + 2) ** 2 / d + 1e-9),
        "E_sAs2_max": E_sAs2,
        "H_qn": H_qn,
        "wick_8": 8.0,
        "residual_vs_wick": E_sAs2 - 8.0,
        "H_holds_qn": bool(E_sAs2 <= H_qn + 1e-6),
    }


def _worker_closed_form_fit(args: tuple) -> dict:
    """Fit λ_cycle as rational function of p,d,n,N."""
    _blas1()
    p, lam_cycle, ray, N, n, d = args
    H = (p + 2) ** 2 / d
    candidates = {
        "3+(p+2)^2/d": 3 + H,
        "3+(p+2)^2/(d+1)": 3 + (p + 2) ** 2 / (d + 1),
        "n/2": n / 2,
        "8": 8.0,
        "4+4/p": 4 + 4 / p,
        "2+2*p/d": 2 + 2 * p / d,
        "2*n/N * something": None,
        # from p=5: 88/13; p=3: 8; p=7: ?
        "4*d/(d-1)": 4 * d / (d - 1) if d > 1 else None,
        "2*(d+1)/(d-1)": 2 * (d + 1) / (d - 1) if d > 1 else None,
        "(p^2+7)/4": (p * p + 7) / 4,  # p=5: 32/4=8 too big
        "2*p": 2 * p,
        "p+1": p + 1,
        "4*(p+1)/p": 4 * (p + 1) / p,
        # λ_cycle = 2N * λ2(P⊙P); maybe λ2(P⊙P)=1/(2d) etc
        "N/d": N / d,
        "2N/d^2 * (p+2)^2": 2 * N * (p + 2) ** 2 / (d**2),  # unlikely
        # exact match attempts for known values
        "88/13_only_p5": 88 / 13 if p == 5 else None,
    }
    errs = {}
    for name, val in candidates.items():
        if val is None:
            continue
        errs[name] = abs(lam_cycle - val)
    best = min(errs.items(), key=lambda kv: kv[1])
    # Fraction of measured
    f = Fraction(lam_cycle).limit_denominator(500)
    fr = Fraction(ray).limit_denominator(500)
    return {
        "p": p,
        "lam_cycle": lam_cycle,
        "ray": ray,
        "H_plus_3": 3 + H,
        "frac_lam": f"{f.numerator}/{f.denominator}",
        "frac_lam_err": abs(lam_cycle - float(f)),
        "frac_ray": f"{fr.numerator}/{fr.denominator}",
        "frac_ray_err": abs(ray - float(fr)),
        "best_candidate": {"name": best[0], "err": best[1]},
        "all_errs": errs,
        "equals_H_plus_3": bool(abs(lam_cycle - (3 + H)) < 1e-6),
    }


def _worker_direct_H_inequality(p: int) -> dict:
    """
    Direct proof ingredients for ray ≤ (p+2)^2/d.

    Idea: for unit Be in image (‖Be‖^2=2 for unit B),
      Be^T G_disj Be = E[(yBy)^2] - 6.

    Write yBy = 2 sum_{i<j} B_ij y_i y_j.
    Use orthogonal decomposition of the edge space into
    stars + image(V_+) + orthogonal complement.

    From Prop 15.56: stars have G-Rayleigh 1 (or map to eigenvalue n/2 for all-ones combination).

    New attempt: bound Be^T G Be using only ‖Be‖ and the fact that
    Be comes from a rank-≤d matrix (B = Vp A Vp.T).

    Be_ij = 2 (Vp A Vp.T)_ij = 2 sum_{αβ} (Vp)_{iα} A_αβ (Vp)_{jβ}.

    This means the matrix B_off = (Be as matrix) has rank ≤ d and factors through Vp.

    Then E[(⟨Be, f⟩)^2] = E[(sum_{i<j} Be_ij y_i y_j)^2].

    Equivalent to Frobenius: E[‖ P_{sym off} (y y^T) projected ‖ related].
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]

    # For each Max+ y, the matrix yy^T - 2 P_+ has Frobenius norm^2 = n(n-2) as before
    # T(x) for x = yBy maximizer direction

    # Claim to test: λ_cycle = max_{‖x‖=1,x⊥1} ‖T(x)‖F^2 / (2N)
    # and ‖T(x)‖F^2 ≤ 2N * (3 + (p+2)^2/d) ‖x‖^2

    # Construct the operator on the d-dimensional space of "linear" features
    # φ_a(y) = (Vp^T y)_a = s_a, and quadratic features s_a s_b.

    # Quadratic features Q_ab(y) = s_a(y) s_b(y), a≤b
    # The space of even functions spanned by {s_a s_b} has dimension d(d+1)/2
    # yBy = s^T A s = sum_{ab} A_ab Q_ab, so image of B→yBy is inside this span.
    # Dimension of zero-diag constraint reduces it to ~ binom(d-1,2).

    S = Y @ Vp  # N x d, rows are s^T
    # Build feature matrix for a≤b: φ_{ab} = S_a * S_b
    feats = []
    pairs = []
    for a in range(d):
        for b in range(a, d):
            feats.append(S[:, a] * S[:, b])
            pairs.append((a, b))
    Phi = np.stack(feats, axis=1)  # N x dimS
    dimS = Phi.shape[1]
    # Center (orthogonal to 1)
    Phi = Phi - Phi.mean(axis=0, keepdims=True)

    # Gram Gφ = Phi.T Phi / N : eigenvalues give E[(s^T A s)^2] for A in Sym
    # For A corresponding to coefficient vector c (with off-diag doubled),
    # s^T A s = c · φ, E[(sAs)^2] = c^T Gφ c.
    Gφ = (Phi.T @ Phi) / N
    evals = np.linalg.eigvalsh(Gφ)
    evals = evals[evals > 1e-10]

    # Now impose zero ambient diag: r_i^T A r_i = 0
    # r_i = Vp row i, (r_i^T A r_i) = sum_{ab} A_ab (Vp_{ia} Vp_{ib}) = 0
    # In c-coordinates: for each i, sum_{a≤b} c_{ab} * factor * Vp_ia Vp_ib = 0
    # A_aa = c_{aa}, A_ab = c_{ab}/2 for a<b if we use s^T A s = sum_a A_aa s_a^2 + 2 sum_{a<b} A_ab s_a s_b
    # Standard: s^T A s = sum_{a≤b} M_{ab} φ_ab with M_aa=A_aa, M_ab=2 A_ab for a<b.

    # Build constraint matrix on c (dimS):
    # For each vertex i: sum_a A_aa Vp_ia^2 + 2 sum_{a<b} A_ab Vp_ia Vp_ib = 0
    cons = []
    for i in range(n):
        row = np.zeros(dimS)
        k = 0
        for a in range(d):
            for b in range(a, d):
                if a == b:
                    row[k] = Vp[i, a] ** 2
                else:
                    row[k] = 2.0 * Vp[i, a] * Vp[i, b]
                k += 1
        cons.append(row)
    Cons = np.stack(cons, axis=0)  # n x dimS
    # Also Tr A = 0: sum_a A_aa = 0
    tr = np.zeros(dimS)
    k = 0
    for a in range(d):
        for b in range(a, d):
            if a == b:
                tr[k] = 1.0
            k += 1
    Cons = np.vstack([Cons, tr])

    # Nullspace of Cons
    _u, s, vh = np.linalg.svd(Cons, full_matrices=True)
    rank = int(np.sum(s > 1e-8))
    null = vh[rank:].T  # dimS x nullity
    nullity = null.shape[1]

    # Restrict Gφ to nullspace: max Rayleigh of null.T @ Gφ @ null
    # But need unit A in Frobenius, not unit c.
    # Map c → A → ‖A‖F
    # c corresponds to M: M_aa=c_aa, M_ab=c_ab for a<b meaning 2A_ab=c_ab, A_ab=c_ab/2
    def c_to_A(c):
        A = np.zeros((d, d))
        k = 0
        for a in range(d):
            for b in range(a, d):
                if a == b:
                    A[a, a] = c[k]
                else:
                    A[a, b] = A[b, a] = c[k] / 2.0
                k += 1
        return A

    # Build metric L such that ‖A‖F^2 = c^T L c
    L = np.zeros((dimS, dimS))
    for k1 in range(dimS):
        e1 = np.zeros(dimS)
        e1[k1] = 1.0
        A1 = c_to_A(e1)
        for k2 in range(dimS):
            e2 = np.zeros(dimS)
            e2[k2] = 1.0
            A2 = c_to_A(e2)
            L[k1, k2] = np.sum(A1 * A2)

    # On nullspace: max c=null @ z,  (c^T Gφ c) / (c^T L c)
    # Generalized eigenproblem: null.T Gφ null z = λ null.T L null z
    Gn = null.T @ Gφ @ null
    Ln = null.T @ L @ null
    # Symmetrize
    Gn = 0.5 * (Gn + Gn.T)
    Ln = 0.5 * (Ln + Ln.T)
    # Solve gen eig
    # Ln should be SPD on null
    try:
        # whitening
        w, V = np.linalg.eigh(Ln)
        w = np.maximum(w, 1e-15)
        W12 = V * np.sqrt(w)
        Winv = V * (1.0 / np.sqrt(w))
        M = Winv.T @ Gn @ Winv
        mevals = np.linalg.eigvalsh(0.5 * (M + M.T))
        lam_max_E = float(mevals[-1])  # = max E[(sAs)^2] = max Q/N for unit B
    except Exception as ex:
        return {"p": p, "error": str(ex)}

    H_qn = 6 + 2 * (p + 2) ** 2 / d
    ray = (lam_max_E - 6.0) / 2.0
    H = (p + 2) ** 2 / d
    return {
        "p": int(p),
        "dimS": int(dimS),
        "nullity": int(nullity),
        "expected_binom": (d - 1) * (d - 2) // 2,
        "lam_max_E_sAs2": lam_max_E,
        "ray": ray,
        "H": H,
        "H_qn": H_qn,
        "H_holds": bool(ray <= H + 1e-6),
        "qn_le_16": bool(lam_max_E <= 16 + 1e-6),
        "Gφ_max_unconstrained": float(evals[-1]) if len(evals) else None,
        "note": "exact max E[(sAs)^2] on zero-diag∩Tr0 via gen.eig",
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"H_proof: workers={W}", flush=True)

    scheme, bounds, direct = [], [], []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5, 7):
            futs[ex.submit(_worker_scheme_params, p)] = ("sch", p)
            futs[ex.submit(_worker_bound_attempt, p)] = ("bnd", p)
            futs[ex.submit(_worker_direct_H_inequality, p)] = ("dir", p)

        for fut in as_completed(futs):
            kind, p = futs[fut]
            r = fut.result()
            if kind == "sch":
                scheme.append(r)
                print(
                    f"  scheme p={p}: λ_cycle={r.get('lam_cycle')} ray={r.get('ray')} "
                    f"H={r.get('H')} holds={r.get('H_holds')} "
                    f"eigs={[e.get('frac') for e in r.get('K_eigs', [])[:5]]}",
                    flush=True,
                )
            elif kind == "bnd":
                bounds.append(r)
                print(
                    f"  bound p={p}: Emax={r['E_sAs2_max']:.6f} H_qn={r['H_qn']:.6f} "
                    f"gersh={r['gershgorin_row_abs_sum']:.4f} gersh_le8={r['gersh_le_8']}",
                    flush=True,
                )
            else:
                direct.append(r)
                print(
                    f"  direct p={p}: ray={r.get('ray')} H={r.get('H')} "
                    f"holds={r.get('H_holds')} nullity={r.get('nullity')} "
                    f"Emax={r.get('lam_max_E_sAs2')}",
                    flush=True,
                )

    # Closed form fit from scheme results
    fits = []
    for r in scheme:
        if r.get("lam_cycle") is not None:
            fits.append(
                _worker_closed_form_fit(
                    (
                        r["p"],
                        r["lam_cycle"],
                        r["ray"],
                        r["N"],
                        r["n"],
                        r["d"],
                    )
                )
            )
            print(
                f"  fit p={r['p']}: lam={fits[-1]['frac_lam']} "
                f"eq_H+3={fits[-1]['equals_H_plus_3']} best={fits[-1]['best_candidate']}",
                flush=True,
            )

    scheme.sort(key=lambda r: r["p"])
    bounds.sort(key=lambda r: r["p"])
    direct.sort(key=lambda r: r["p"])

    H_ok = all(r.get("H_holds") for r in scheme if r.get("H_holds") is not None)
    H_ok = H_ok and all(r.get("H_holds") for r in direct if r.get("H_holds") is not None)

    # Algebra already proved in Prop 15.63
    out = {
        "hypothesis_H": "ray ≤ (p+2)^2/d for unit zero-diag B on V_+",
        "algebra_H_le_5": "proved in Prop 15.63.1",
        "scheme": scheme,
        "bounds": bounds,
        "direct_exact": direct,
        "fits": fits,
        "H_certified": H_ok,
        "status": (
            ("H CERTIFIED p=3,5,7 via exact gen.eig + scheme; " if H_ok else "H FAIL; ")
            + "uniform proof OPEN unless closed form locks. L OPEN."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_H_proof.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
