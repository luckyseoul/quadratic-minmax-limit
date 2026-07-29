#!/usr/bin/env python3
"""
Prop 15.63 attack: prove Q_4 ≤ 10 N ‖B‖_F² for zero-diag B on V_+.

Approaches probed (multi-worker):
  A) Closed-form candidates for λ_max(Q) / N
  B) Exact Q via moments: Tr(B^4), Tr((CB)^2), row structure
  C) Image of edge map: Rayleigh of Gu_disj vs combinatorial UB
  D) Cumulant residual R_ijkl vs Wick; bound on image
  E) Irrep spectrum pattern (p=3,5) for formula guess

Writes evidence/e1_gmin_q4_proof.json
"""
from __future__ import annotations

import json
import os
import sys
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


def _maximizer_B(Y, Vp, d, seed=0):
    rng = np.random.default_rng(seed)
    U = (Y @ Vp) / np.sqrt(Y.shape[1])
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(80):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    return B / (np.linalg.norm(B) + 1e-30), A


def _worker_invariants(p: int) -> dict:
    """Compute polynomial invariants of maximizer B and test closed forms."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    B, A = _maximizer_B(Y, Vp, d, seed=p + 7)
    bf2 = float(np.sum(B * B))
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))
    Q4 = Q - 6.0 * N * bf2

    # Invariants
    trB2 = bf2
    trB4 = float(np.trace(B @ B @ B @ B))
    trCB2 = float(np.trace(C @ B @ B))  # = p ‖B‖F² on V+
    trCB4 = float(np.trace(C @ B @ B @ B @ B))
    # row structure
    row_n2 = np.sum(B * B, axis=1)
    # By moments
    By = Y @ B  # N x n? wait Y is N x n, B is n x n → Y@B is N x n
    # actually By for each y: B @ y, so (B @ Y.T).T
    BY = (B @ Y.T).T  # N x n
    By_n2 = np.sum(BY * BY, axis=1)
    yBy = ytBy
    # alignment cos^2 = (yBy)^2 / (n * ‖By‖^2)
    cos2 = (yBy**2) / (n * By_n2 + 1e-30)

    # Candidate formulas for Q/N (unit B)
    qn = Q / (N * bf2)
    candidates = {
        "16": 16.0,
        "16*(d-2)/d": 16.0 * (d - 2) / d,
        "16*d/(d+2)": 16.0 * d / (d + 2),  # 2×sphere full
        "8": 8.0,  # Wick
        "8+8*(d-2)/(d+2)": 8.0 + 8.0 * (d - 2) / (d + 2),
        "2*(p+1)": 2.0 * (p + 1),
        "4*p": 4.0 * p,
        "2*n/d": 2.0 * n / d,  # = 4
        "12+4*(d-2)/d": 12.0 + 4.0 * (d - 2) / d,
        "6+10*(d-2)/d": 6.0 + 10.0 * (d - 2) / d,  # if Q4=10N*(d-2)/d
        "6+98/(d)": 6.0 + 98.0 / d,  # p=5: 6+98/13=6+7.538=13.538
        "6+4*(p-1)": 6.0 + 4.0 * (p - 1),
        "6+2*(p+1)": 6.0 + 2.0 * (p + 1),
        "2*(d+3)/(d-1)*something": None,
    }
    # Fit Q4/N = ?
    q4n = Q4 / (N * bf2)
    # Exact fractions
    def near_frac(x, max_den=200):
        f = Fraction(x).limit_denominator(max_den)
        return f"{f.numerator}/{f.denominator}", float(f), abs(x - float(f))

    qn_f, qn_fv, qn_err = near_frac(qn)
    q4n_f, q4n_fv, q4n_err = near_frac(q4n)
    ratio16 = Q / (16.0 * N * bf2)
    r16_f, r16_fv, r16_err = near_frac(ratio16)

    # More candidates involving p only
    extra = {
        "qn_actual": qn,
        "q4n_actual": q4n,
        "ratio16_actual": ratio16,
        "qn_frac": qn_f,
        "q4n_frac": q4n_f,
        "ratio16_frac": r16_f,
        "qn_frac_err": qn_err,
        "q4n_frac_err": q4n_err,
        "ratio16_frac_err": r16_err,
    }
    # Test: Q4/N = 2(p^2-1)/p ? etc
    trial_q4n = {
        "10": 10.0,
        "10*(d-2)/d": 10.0 * (d - 2) / d,
        "10*(d-4)/(d-2)": 10.0 * (d - 4) / (d - 2) if d > 2 else None,
        "8*(d-2)/(d+2)": 8.0 * (d - 2) / (d + 2),
        "2*(p+1)": 2.0 * (p + 1),
        "4*(p-2)": 4.0 * (p - 2),
        "2*p": 2.0 * p,
        "p+3": p + 3.0,
        "98/d": 98.0 / d,
        "2*(n-4)/(d)": 2.0 * (n - 4) / d,
        "4*(d-1)/d": 4.0 * (d - 1) / d,
        "8*(d-1)/(d+1)": 8.0 * (d - 1) / (d + 1),
        # from 16N*(d-2)/d - 6N = N(16(d-2)/d - 6) = N(16d-32-6d)/d = N(10d-32)/d
        "(10d-32)/d": (10.0 * d - 32.0) / d,
        # sphere-based residual budget
        "16d/(d+2)-6": 16.0 * d / (d + 2) - 6.0,
    }

    matches = {}
    for name, val in trial_q4n.items():
        if val is None:
            continue
        err = abs(q4n - val)
        matches[name] = {"value": val, "abs_err": err, "rel_err": err / max(abs(q4n), 1e-12)}

    best = min(matches.items(), key=lambda kv: kv[1]["abs_err"])

    return {
        "p": int(p),
        "n": int(n),
        "d": int(d),
        "N": int(N),
        "Q": Q,
        "Q4": Q4,
        "qn": qn,
        "q4n": q4n,
        "ratio16": ratio16,
        "trB4": trB4,
        "trCB2": trCB2,
        "trCB2_over_p": trCB2 / p,
        "trCB4": trCB4,
        "row_n2_mean": float(np.mean(row_n2)),
        "row_n2_max": float(np.max(row_n2)),
        "By2_mean": float(np.mean(By_n2)),
        "cos2_mean": float(np.mean(cos2)),
        "cos2_max": float(np.max(cos2)),
        "yBy2_mean": float(np.mean(yBy**2)),
        "fracs": extra,
        "q4n_trials": matches,
        "best_q4n_trial": {"name": best[0], **best[1]},
        "diagmax": float(np.max(np.abs(np.diag(B)))),
    }


def _worker_algebra_expand(p: int) -> dict:
    """
    Expand Q for zero-diag B on V_+ using only Σ and C.
    Test whether Q = a‖B‖F² + b Tr(B^4) + c something closed.
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
    null = _null_zero_diag(Vp)
    nullity = null.shape[1]
    rng = np.random.default_rng(p * 17)

    # Sample many B, collect features for linear regression: Q ~ c0*N*bf2 + c1*N*trB4 + ...
    rows = []
    for s in range(60):
        coef = rng.standard_normal(nullity)
        v = null @ coef
        A = _A_from_vec(v, d)
        B = Vp @ A @ Vp.T
        nrm = np.linalg.norm(B)
        if nrm < 1e-14:
            continue
        B = B / nrm
        bf2 = 1.0
        ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
        Q = float(np.sum(ytBy**2))
        trB4 = float(np.trace(np.linalg.matrix_power(B, 4)))
        trB3 = float(np.trace(B @ B @ B))  # 0 for? 
        # Hadamard / Schur with C
        Bc = B * C  # entrywise
        trBc2 = float(np.sum(Bc * Bc))
        # sum_ij B_ij^4
        b4 = float(np.sum(B**4))
        # sum_i (sum_j B_ij^2)^2 = sum_i ‖row_i‖^4
        row2 = np.sum(B * B, axis=1)
        sum_row2_sq = float(np.sum(row2**2))
        # sum_{i≠j≠k} B_ij^2 B_ik^2 etc already in row
        rows.append(
            {
                "Q": Q,
                "Nbf2": N * bf2,
                "NtrB4": N * trB4,
                "Nb4": N * b4,
                "Nsum_row2sq": N * sum_row2_sq,
                "NtrBc2": N * trBc2,
                "trB4": trB4,
                "b4": b4,
                "sum_row2sq": sum_row2_sq,
            }
        )

    # Also maximizer
    B, _ = _maximizer_B(Y, Vp, d, seed=1)
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Qm = float(np.sum(ytBy**2))
    trB4m = float(np.trace(np.linalg.matrix_power(B, 4)))
    b4m = float(np.sum(B**4))
    row2 = np.sum(B * B, axis=1)
    sum_row2_sq_m = float(np.sum(row2**2))

    # Linear least squares: Q = a*N + b*N*trB4 + c*N*b4 + d*N*sum_row2sq
    # (unit B)
    X = []
    y = []
    for r in rows:
        X.append([r["Nbf2"], r["NtrB4"], r["Nb4"], r["Nsum_row2sq"], r["NtrBc2"]])
        y.append(r["Q"])
    X = np.asarray(X)
    y = np.asarray(y)
    coef, residuals, rank, svals = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    max_abs_resid = float(np.max(np.abs(pred - y)))
    rms = float(np.sqrt(np.mean((pred - y) ** 2)))

    # Test exact identity candidates from graph counting
    # For boolean Rademacher chaos of order 2, there is a formula in terms of
    # ‖B‖F², sum_row2sq, trB4, and 4-cycle signed weights.
    # Hypothesis: for B on V_+ zero-diag of conference,
    # Q = α N ‖B‖F² + β N sum_i ‖row_i‖^4 + γ N tr(B^4) + δ N * (signed 4-cycles)
    #
    # If residual of lstsq is ~0 with these features, identity is in this span.

    # Check if Q = 2N*bf2 + 4N*bf2 + something*trB4
    # i.e. already know 6N + Q4, Q4 related to trB4?
    q4s = y - 6 * N
    trB4s = X[:, 1] / N
    # correlation
    if len(q4s) > 2:
        corr = float(np.corrcoef(q4s, trB4s)[0, 1])
    else:
        corr = None

    return {
        "p": int(p),
        "n_samples": len(rows),
        "lstsq_coef_Nbf2_NtrB4_Nb4_Nrow2sq_NtrBc2": [float(c) for c in coef],
        "lstsq_max_abs_resid": max_abs_resid,
        "lstsq_rms": rms,
        "lstsq_rank": int(rank),
        "Q4_trB4_corr": corr,
        "maximizer": {
            "Q": Qm,
            "trB4": trB4m,
            "b4": b4m,
            "sum_row2sq": sum_row2_sq_m,
            "Q4": Qm - 6 * N,
            "Q4_over_N": (Qm - 6 * N) / N,
            "trB4_over": trB4m,
        },
        "identity_in_span": bool(max_abs_resid < 1e-4 * max(float(np.max(y)), 1.0)),
    }


def _worker_disj_m4_image(p: int) -> dict:
    """
    Bound attempt: Q4 = N * Be^T M4_disj Be.
    On image, compare to ‖Be‖² * (max Rayleigh).
    Also test: |m4| ≤ (p-2)/(p(p+1)) or known L(p) style bounds.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    if p == 7:
        return {"p": 7, "skipped": True, "reason": "E large"}

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    G = (F.T @ F) / N  # E[f f']
    share = np.zeros((E, E), dtype=bool)
    for a, (i, j) in enumerate(edges):
        for b in range(a + 1, E):
            k, l = edges[b]
            if len({i, j, k, l}) == 3:
                share[a, b] = share[b, a] = True
    disj = ~share
    np.fill_diagonal(disj, False)
    Gdisj = G * disj

    # m4 stats already known; compute crude CS bound
    # |Be^T Gdisj Be| ≤ ‖Gdisj‖_F ‖Be Be^T‖_F = ‖Gdisj‖_F ‖Be‖²
    fro = float(np.linalg.norm(Gdisj))
    # Better: spectral on full, then image
    evals = np.linalg.eigvalsh(Gdisj)
    lam_max_full = float(evals[-1])
    lam_min_full = float(evals[0])

    # Build image basis of Be from zero-diag∩V_+
    null = _null_zero_diag(Vp)
    Bes = []
    for t in range(null.shape[1]):
        A = _A_from_vec(null[:, t], d)
        B = Vp @ A @ Vp.T
        Be = np.array([2.0 * B[i, j] for i, j in edges])
        for Bb in Bes:
            Be = Be - np.dot(Be, Bb) * Bb
        nrm = np.linalg.norm(Be)
        if nrm > 1e-12:
            Bes.append(Be / nrm)
    m = len(Bes)
    BesM = np.stack(Bes, axis=1)  # E x m
    # Rayleigh matrix: BesM.T @ Gdisj @ BesM  (for unit Be in image, Q4/N = b^T R b)
    # But ‖Be‖² = 2 for unit B, so if we use ONB of Be-space with ‖Be‖=1,
    # Q4 = N * (Be^T Gdisj Be) and for unit B, Be has norm √2, so
    # max Q4 = N * 2 * λ_max(R) where R is Rayleigh on unit sphere in Be-image.
    R = BesM.T @ Gdisj @ BesM
    revals = np.linalg.eigvalsh(0.5 * (R + R.T))
    # max Q4 for unit B: Be = √2 * u, u unit in image ⇒ Q4 = N * 2 * u^T Gdisj u = 2N λ_max(R)
    q4_max = 2.0 * N * float(revals[-1])
    thr = 10.0 * N

    # Sign pattern: is Gdisj positive on image?
    return {
        "p": int(p),
        "skipped": False,
        "E": int(E),
        "image_dim": int(m),
        "Gdisj_fro": fro,
        "Gdisj_lam_max": lam_max_full,
        "Gdisj_lam_min": lam_min_full,
        "image_ray_max": float(revals[-1]),
        "image_ray_min": float(revals[0]),
        "image_ray_unique_r6": sorted({round(float(x), 6) for x in revals}),
        "Q4_max_from_image": q4_max,
        "thr_10N": thr,
        "bound_ok": bool(q4_max <= thr + 1e-4 * thr),
        "ratio": q4_max / thr,
        # crude Frobenius bound for unit B: Q4 ≤ N * fro * ‖Be‖² = N*fro*2
        "fro_bound_Q4": 2.0 * N * fro,
        "fro_bound_ok": bool(2.0 * N * fro <= thr),
    }


def _worker_proof_attempt_c4(p: int) -> dict:
    """
    Attempt: for zero-diag B on V_+,
      Q4 = 2 N sum_{i<j,k<l disjoint} 4 B_ij B_kl m4
    Use that B = P_+ B P_+ and conference:
      m4(i,j,k,l) = (1/p^2) * poly(C) + κ residual.

    Direct identity test:
      Q = sum_y (2 sum_{i<j} B_ij y_i y_j)^2
        = 4 sum_{i<j} B_ij^2 * N  + 8 sum_wedge ... + 8 sum_disj B_ij B_kl * sum yiyjykyl
    Already have typeA=2N (in Be language) = 4 * sum_{i<j} B_ij^2 * N with Be=2B.

    New: relate Q4 to ‖B⊙B‖ and cycles in C.
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

    B, A = _maximizer_B(Y, Vp, d, seed=3)
    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    Q = float(np.sum(ytBy**2))
    Q4 = Q - 6 * N

    # Cycle form: for Seidel, tr(B^4) counts weighted closed walks
    trB4 = float(np.trace(np.linalg.matrix_power(B, 4)))
    # sum of B_ij B_jk B_kl B_li over distinct 4-cycles
    # tr(B^4) = sum_{ijkl} B_ij B_jk B_kl B_li
    # includes non-simple walks

    # Boolean identity attempt from literature / chaos:
    # For x in {±1}^n, (x^T B x)^2 with B zero-diag symmetric:
    # (x^T B x)^2 = 2‖B‖F² * n? no only sum over cube.
    # On Max+ only.

    # Use conditional: since Cy=py, the coordinates are dependent.
    # y_∞ = 1 (affine), etc.

    # Compute Q via double count of pairs of edges
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    Be = np.array([2.0 * B[i, j] for i, j in edges])
    F = np.stack([Y[:, i] * Y[:, j] for i, j in edges], axis=1)
    Gu = F.T @ F
    Q_from_Be = float(Be @ Gu @ Be)
    assert abs(Q_from_Be - Q) < 1e-4

    # Wick on full Gaussian with Σ=2P_+: E[(yBy)^2]=2‖Σ^{1/2} B Σ^{1/2}‖F² + (Tr ΣB)^2
    # = 2‖2B‖F²? ΣB = 2B on V_+, Tr(ΣB)=2 Tr B=0
    # E = 2 * 4 ‖B‖F² = 8 ‖B‖F², sum = 8N
    wick = 8.0 * N
    residual_vs_wick = Q - wick  # can be negative for random B

    # Key inequality attempt: Q ≤ wick + |κ4 contraction|
    # κ4 contraction for maximizer
    return {
        "p": int(p),
        "Q": Q,
        "Q4": Q4,
        "wick": wick,
        "Q_minus_wick": Q - wick,
        "trB4": trB4,
        "Q4_over_trB4": Q4 / trB4 if abs(trB4) > 1e-12 else None,
        "thr_Q4": 10.0 * N,
        "ok": Q4 <= 10.0 * N + 1e-4,
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"q4_proof: workers={W}", flush=True)

    inv, alg, img, att = [], [], [], []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5, 7):
            futs[ex.submit(_worker_invariants, p)] = ("inv", p)
            futs[ex.submit(_worker_algebra_expand, p)] = ("alg", p)
            futs[ex.submit(_worker_proof_attempt_c4, p)] = ("att", p)
        for p in (3, 5):
            futs[ex.submit(_worker_disj_m4_image, p)] = ("img", p)

        for fut in as_completed(futs):
            kind, p = futs[fut]
            r = fut.result()
            if kind == "inv":
                inv.append(r)
                print(
                    f"  inv p={p}: qn={r['qn']:.6f} ({r['fracs']['qn_frac']}) "
                    f"q4n={r['q4n']:.6f} ({r['fracs']['q4n_frac']}) "
                    f"r16={r['ratio16']:.6f} ({r['fracs']['ratio16_frac']}) "
                    f"best_trial={r['best_q4n_trial']['name']} "
                    f"err={r['best_q4n_trial']['abs_err']:.2e}",
                    flush=True,
                )
            elif kind == "alg":
                alg.append(r)
                print(
                    f"  alg p={p}: resid={r['lstsq_max_abs_resid']:.3e} "
                    f"in_span={r['identity_in_span']} "
                    f"coef={r['lstsq_coef_Nbf2_NtrB4_Nb4_Nrow2sq_NtrBc2']} "
                    f"corr_Q4_trB4={r['Q4_trB4_corr']}",
                    flush=True,
                )
            elif kind == "img":
                img.append(r)
                print(
                    f"  img p={p}: Q4max={r.get('Q4_max_from_image')} "
                    f"thr={r.get('thr_10N')} ok={r.get('bound_ok')} "
                    f"ray_max={r.get('image_ray_max')} "
                    f"rays={r.get('image_ray_unique_r6')}",
                    flush=True,
                )
            else:
                att.append(r)
                print(
                    f"  att p={p}: Q4={r['Q4']:.4f} thr={r['thr_Q4']} ok={r['ok']} "
                    f"Q-wick={r['Q_minus_wick']:.4f} Q4/trB4={r['Q4_over_trB4']}",
                    flush=True,
                )

    inv.sort(key=lambda r: r["p"])
    alg.sort(key=lambda r: r["p"])
    img.sort(key=lambda r: r["p"])
    att.sort(key=lambda r: r["p"])

    # Analyze formula match across p
    formula_hits = {}
    for r in inv:
        for name, info in r["q4n_trials"].items():
            formula_hits.setdefault(name, []).append(
                {"p": r["p"], "err": info["abs_err"], "rel": info["rel_err"]}
            )
    uniform = {}
    for name, lst in formula_hits.items():
        max_err = max(x["err"] for x in lst)
        uniform[name] = {
            "max_abs_err": max_err,
            "all_p": lst,
            "plausible": bool(max_err < 1e-4),
        }

    out = {
        "invariants": inv,
        "algebra_expand": alg,
        "image_rayleigh": img,
        "wick_compare": att,
        "uniform_formula_screen": {
            k: v for k, v in sorted(uniform.items(), key=lambda kv: kv[1]["max_abs_err"])[:12]
        },
        "status": (
            "Q4 proof attack: screened closed forms + image Rayleigh + "
            "polynomial span. Uniform Q4≤10N still OPEN unless a formula locks."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_q4_proof.json"
    path.write_text(json.dumps(out, indent=2))
    print("uniform plausible:", [k for k, v in uniform.items() if v["plausible"]], flush=True)
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
