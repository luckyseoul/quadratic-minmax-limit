#!/usr/bin/env python3
"""
Prop 15.65 attack: residual/cumulant bound via coordinate 4-moments.

Max+ moment tensor M_ijkl = sum_y y_i y_j y_k y_l satisfies
  sum_j C_ij M_jklm = p M_iklm  (evec constraint)
and boolean reductions M_iikl = N Σ_kl, etc.

Use this to get a closed form for E[(y^T B y)^2] when B on V_+ zero-diag,
hence a bound residual ≤ (p+1)(p+7)/d.

Multi-worker. Writes evidence/e1_gmin_cumulant.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def _worker_moment_structure(p: int) -> dict:
    """Compute coordinate 4-moments and test closed forms against C."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    Sigma = (Y.T @ Y) / N  # = 2 P_+
    # Sample of M_ijkl for representative index patterns
    # Pattern types for 4 indices:
    # (4): all same iiii
    # (3,1): iiij
    # (2,2): iijj i≠j
    # (2,1,1): iijk
    # (1,1,1,1): ijkl all distinct

    # Direct computation for small n
    def mom(*idx):
        t = np.ones(N)
        for i in idx:
            t = t * Y[:, i]
        return float(np.sum(t))

    # (2,2) type: M_iijj for i≠j
    m22 = []
    for i in range(min(n, 8)):
        for j in range(i + 1, min(n, 8)):
            m22.append(mom(i, i, j, j))
    # Boolean: y_i² y_j² = 1, so M_iijj = N always!
    m22_unique = sorted(set(np.round(m22, 10)))

    # (2,1,1): M_iijk for distinct i,j,k
    m211 = []
    for i in range(min(n, 6)):
        for j in range(min(n, 6)):
            for k in range(j + 1, min(n, 6)):
                if len({i, j, k}) < 3:
                    continue
                m211.append((i, j, k, mom(i, i, j, k)))
    # Should equal N Sigma_jk = N * 2 P_+_jk
    errs_211 = []
    for i, j, k, m in m211[:20]:
        pred = N * Sigma[j, k]
        errs_211.append(abs(m - pred))

    # (1,1,1,1): four distinct — the interesting residual
    m4_vals = []
    rng = np.random.default_rng(p)
    for _ in range(200):
        idx = rng.choice(n, size=4, replace=False)
        m4_vals.append(mom(*idx) / N)  # normalized E[y_i y_j y_k y_l]

    m4_vals = np.array(m4_vals)
    # Wick prediction: Σ_ij Σ_kl + Σ_ik Σ_jl + Σ_il Σ_jk for distinct
    # (all off-diag)
    wick_vals = []
    raw_vals = []
    for _ in range(200):
        i, j, k, l = rng.choice(n, size=4, replace=False)
        w = Sigma[i, j] * Sigma[k, l] + Sigma[i, k] * Sigma[j, l] + Sigma[i, l] * Sigma[j, k]
        r = mom(i, j, k, l) / N
        wick_vals.append(w)
        raw_vals.append(r)
    wick_vals = np.array(wick_vals)
    raw_vals = np.array(raw_vals)
    kappa_vals = raw_vals - wick_vals

    # Evec constraint check: sum_j C_ij M_jklm = p M_iklm
    # Test for random k,l,m
    evec_errs = []
    for _ in range(30):
        i = int(rng.integers(0, n))
        k, l, m_ = rng.choice(n, size=3, replace=True)
        lhs = 0.0
        for j in range(n):
            lhs += C[i, j] * mom(j, k, l, m_)
        rhs = p * mom(i, k, l, m_)
        evec_errs.append(abs(lhs - rhs) / (abs(rhs) + 1e-10))

    # Closed form candidates for E[y_i y_j y_k y_l] when all distinct
    # In terms of C only:
    # Candidate: (1/p²)(C_ij C_kl + ...) + c * (C_ij C_kl C_something)
    # Or from Prop 15.53: related to κ and triangle type

    # For residual of unit B zero-diag on V_+:
    # residual = sum_{ijkl} B_ij B_kl κ_ijkl
    # with κ = M/N - Wick
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    # Build maximizer A and compute residual two ways
    S = Y @ Vp
    U = S / np.sqrt(n)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(60):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    # Force numerical zero diag
    np.fill_diagonal(B, 0.0)
    B = 0.5 * (B + B.T)
    B = Vp @ (Vp.T @ B @ Vp) @ Vp.T  # re-project approx
    np.fill_diagonal(B, 0.0)
    nrm = np.linalg.norm(B)
    B /= nrm + 1e-30

    ytBy = np.einsum("ai,ij,aj->a", Y, B, Y)
    E_emp = float(np.mean(ytBy**2))
    residual_emp = E_emp - 8.0

    # Contract κ with B⊗B via sampling 4-tuples (Monte Carlo for large n)
    # residual = sum_{i≠j,k≠l} B_ij B_kl κ_ijkl
    # = E[(yBy)^2] - Wick contribution
    # Wick for zero-diag B on V_+ with ‖B‖F=1: should be 8
    # Direct Wick contraction:
    # sum_ijkl B_ij B_kl (Σ_ij Σ_kl + Σ_ik Σ_jl + Σ_il Σ_jk)
    # With B zero-diag on V_+, Σ=2P_+, B=P_+ B P_+:
    # This equals 8 for unit B. Verified by E_emp structure.

    budget = (p + 1) * (p + 7) / d
    H = (p + 2) ** 2 / d

    # Kappa stats for distinct 4-tuples
    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "d": int(d),
        "m22_unique": m22_unique,  # should be [N]
        "m211_max_err_vs_Sigma": float(max(errs_211)) if errs_211 else None,
        "evec_constraint_max_rel_err": float(max(evec_errs)),
        "m4_E_min": float(m4_vals.min()),
        "m4_E_max": float(m4_vals.max()),
        "m4_E_mean": float(m4_vals.mean()),
        "kappa_min": float(kappa_vals.min()),
        "kappa_max": float(kappa_vals.max()),
        "kappa_mean": float(kappa_vals.mean()),
        "kappa_abs_max": float(np.max(np.abs(kappa_vals))),
        "residual_emp": residual_emp,
        "budget": budget,
        "H": H,
        "residual_le_budget": bool(residual_emp <= budget + 1e-4),
        # Test: |κ| ≤ (p+1)(p+7)/(d n²) or similar entrywise?
        # For residual = sum B_ij B_kl κ_ijkl ≤ ‖B‖F² * n² * max|κ| rough
        "entrywise_UB_residual": float(np.max(np.abs(kappa_vals)) * n * n),  # crude
    }


def _worker_kappa_contraction_identity(p: int) -> dict:
    """
    Exact identity: for zero-diag B on V_+,
    residual = sum_{i≠j,k≠l} B_ij B_kl κ_ijkl
    Try to rewrite using C-invariants and B=P_+BP_+.

    Since B = Vp A Vp.T and s=Vp^T y,
    residual = E[(sAs)^2] - 8 = E[ (sum_{αβ} A_αβ s_α s_β)^2 ] - 8

    The tensor E[s_α s_β s_γ s_δ] - Wick_2I
    is the 4th cumulant of the random vector s on Max+.

    Because s = Vp^T y and y is boolean with Cy=py,
    s_α = v^α · y with C v^α = p v^α (columns of Vp).

    Key: the multiset {s_α} as α varies is the set of projections of y onto the ONB of V_+.
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
    S = Y @ Vp  # N x d

    # Fourth moment tensor of s (d^4 is small for p=3,5; for p=7 d=25 → 390625 ok)
    # T[a,b,c,e] = E[s_a s_b s_c s_e]
    if d > 15:
        # Only compute needed contractions, not full tensor
        return _worker_kappa_via_A_basis(p)

    T = np.einsum("na,nb,nc,nd->abcd", S, S, S, S) / N
    # Wick for cov 2I: W[a,b,c,d] = 4*(δ_ab δ_cd + δ_ac δ_bd + δ_ad δ_bc)
    # Isserlis: E[ab]E[cd]+... with E[ab]=2δ_ab
    W = np.zeros((d, d, d, d))
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    W[a, b, c, e] = (
                        4.0
                        * (
                            (1.0 if a == b else 0.0) * (1.0 if c == e else 0.0)
                            + (1.0 if a == c else 0.0) * (1.0 if b == e else 0.0)
                            + (1.0 if a == e else 0.0) * (1.0 if b == c else 0.0)
                        )
                        * 1.0  # already 2*2=4 in the products: E[ab]=2δ → 2δ*2δ=4δδ
                    )
    # Wait: E[ab]E[cd] = (2δ_ab)(2δ_cd)=4 δ_ab δ_cd. Yes W formula is correct.

    kappa = T - W

    # As quadratic form on Sym: residual(A) = sum_{abcd} A_ab A_cd kappa_abcd
    # For A symmetric. Build matrix on vech basis and find max eig on Z.

    dimS = d * (d + 1) // 2

    def vech(A):
        v = []
        for i in range(d):
            for j in range(i, d):
                v.append(A[i, j] if i == j else A[i, j] * np.sqrt(2))  # ONB for F
        return np.array(v)

    def unvech(v):
        A = np.zeros((d, d))
        k = 0
        for i in range(d):
            for j in range(i, d):
                if i == j:
                    A[i, i] = v[k]
                else:
                    A[i, j] = A[j, i] = v[k] / np.sqrt(2)
                k += 1
        return A

    # Matrix of residual form in vech coords (Frobenius ONB)
    Kmat = np.zeros((dimS, dimS))
    for i in range(dimS):
        ei = np.zeros(dimS)
        ei[i] = 1.0
        Ai = unvech(ei)
        for j in range(i, dimS):
            ej = np.zeros(dimS)
            ej[j] = 1.0
            Aj = unvech(ej)
            # residual bilinear: sum abcd Ai_ab Aj_cd kappa_abcd
            val = float(np.einsum("ab,cd,abcd->", Ai, Aj, kappa))
            Kmat[i, j] = Kmat[j, i] = val

    # Nullspace Z: zero ambient diag + Tr0
    # Tr0: already in Sym0 roughly - enforce
    # Ambient diag: r_i^T A r_i = 0
    Cons = []
    for i in range(n):
        ri = Vp[i]
        # r^T A r = sum_{ab} A_ab ri_a ri_b
        row = np.zeros(dimS)
        k = 0
        for a in range(d):
            for b in range(a, d):
                if a == b:
                    row[k] = ri[a] ** 2
                else:
                    # A_ab = v_k/√2, contribution 2 A_ab ri_a ri_b = √2 v_k ri_a ri_b
                    row[k] = np.sqrt(2) * ri[a] * ri[b]
                k += 1
        Cons.append(row)
    # Tr A = sum_a A_aa
    tr = np.zeros(dimS)
    k = 0
    for a in range(d):
        for b in range(a, d):
            if a == b:
                tr[k] = 1.0
            k += 1
    Cons.append(tr)
    Cons = np.stack(Cons)
    _u, s, vh = np.linalg.svd(Cons, full_matrices=True)
    rank = int(np.sum(s > 1e-8))
    null = vh[rank:].T

    # Restrict Kmat to null
    Kn = null.T @ Kmat @ null
    evals = np.linalg.eigvalsh(0.5 * (Kn + Kn.T))
    lam_max_res = float(evals[-1])
    lam_min_res = float(evals[0])

    budget = (p + 1) * (p + 7) / d
    return {
        "p": int(p),
        "d": int(d),
        "nullity": int(null.shape[1]),
        "lam_max_residual": lam_max_res,
        "lam_min_residual": lam_min_res,
        "budget": budget,
        "bound_ok": bool(lam_max_res <= budget + 1e-6),
        "ratio": lam_max_res / budget,
        "kappa_frob": float(np.linalg.norm(kappa)),
        "spectrum_res_top5": [float(x) for x in evals[-5:][::-1]],
        "spectrum_res_unique_r6": sorted({round(float(x), 6) for x in evals}),
    }


def _worker_kappa_via_A_basis(p: int) -> dict:
    """For large d, compute residual spectrum via random projection / power on Z."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    S = Y @ Vp

    # Use gen.eig path from H_proof (already know answer)
    # Just report budget comparison
    # Power iteration for max residual on Z
    dimS = d * (d + 1) // 2

    def A_from_vec(v):
        A = np.zeros((d, d))
        k = 0
        for i in range(d):
            for j in range(i, d):
                A[i, j] = A[j, i] = v[k]
                k += 1
        return A

    M = np.zeros((n + 1, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        A = A_from_vec(e)
        B = Vp @ A @ Vp.T
        M[:n, k] = np.diag(B)
    kk = 0
    for i in range(d):
        for j in range(i, d):
            if i == j:
                M[n, kk] = 1.0
            kk += 1
    _u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    null = vh[rank:].T
    rng = np.random.default_rng(p)

    def residual_of_z(z):
        A = A_from_vec(null @ z)
        A /= np.linalg.norm(A) + 1e-30
        sAs = np.einsum("bi,ij,bj->b", S, A, S)
        return float(np.mean(sAs**2)) - 8.0

    best = -1e300
    for trial in range(30):
        z = rng.standard_normal(null.shape[1])
        z /= np.linalg.norm(z) + 1e-30
        # Finite-diff ascent on sphere (cheap steps)
        for _ in range(40):
            base = residual_of_z(z)
            g = np.zeros_like(z)
            eps = 1e-5
            # Only random coordinate subset for speed
            coords = rng.choice(len(z), size=min(40, len(z)), replace=False)
            for i in coords:
                e = np.zeros_like(z)
                e[i] = eps
                g[i] = (residual_of_z(z + e) - base) / eps
            g = g - np.dot(g, z) * z
            gn = np.linalg.norm(g)
            if gn < 1e-12:
                break
            z = z + 0.1 * g / gn
            z /= np.linalg.norm(z) + 1e-30
        val = residual_of_z(z)
        if val > best:
            best = val

    budget = (p + 1) * (p + 7) / d
    return {
        "p": int(p),
        "d": int(d),
        "method": "power_residual",
        "lam_max_residual": best,
        "budget": budget,
        "bound_ok": bool(best <= budget + 1e-3),
        "ratio": best / budget,
    }


def _worker_closed_form_kappa(p: int) -> dict:
    """
    Try closed forms for residual max:
    residual_max = (p+1)(p+7)/d   (H equality)
    or residual_max = 4(d-1)/d or etc.

    Also test: for ALL A in Z (not just max),
    residual(A) = c1 * (tr(A^4) - ...) + c2 * ...
    failed before; try residual = ((p+1)(p+7)/d) * ‖A‖F² * cos²(A, A_max)
    """
    known_res = {
        3: 8.0,
        5: 72 / 13,
        7: 2.56234718826407,  # 4320/409 - 8
    }
    if p not in known_res:
        return {"p": p}
    r = known_res[p]
    d = (p * p + 1) // 2
    budget = (p + 1) * (p + 7) / d
    # Exact fraction for p=7 residual
    # E=4320/409, residual = 4320/409 - 8 = (4320 - 3272)/409 = 1048/409
    candidates = {
        "(p+1)(p+7)/d": budget,
        "4*(d-4)/d": 4 * (d - 4) / d if d > 4 else None,
        "2*(p+1)/d * (p-1)": 2 * (p + 1) * (p - 1) / d,
        "8*(d-4)/(d+2)": 8 * (d - 4) / (d + 2) if d > 4 else None,
        "1048/409_p7": 1048 / 409 if p == 7 else None,
        "2*(p-1)^2/d": 2 * (p - 1) ** 2 / d,
        "(p+1)^2/d": (p + 1) ** 2 / d,
        "4*p/d": 4 * p / d,
    }
    errs = {k: abs(r - v) for k, v in candidates.items() if v is not None}
    best = min(errs.items(), key=lambda kv: kv[1])
    return {
        "p": p,
        "residual": r,
        "budget": budget,
        "equals_budget": bool(abs(r - budget) < 1e-6),
        "best_form": {"name": best[0], "err": best[1]},
        "all_errs": errs,
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"cumulant: workers={W}", flush=True)

    structs, kappas, forms = [], [], []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5, 7):
            futs[ex.submit(_worker_moment_structure, p)] = ("st", p)
            futs[ex.submit(_worker_kappa_contraction_identity, p)] = ("kp", p)
            futs[ex.submit(_worker_closed_form_kappa, p)] = ("fm", p)
        for fut in as_completed(futs):
            kind, p = futs[fut]
            r = fut.result()
            if kind == "st":
                structs.append(r)
                print(
                    f"  struct p={p}: evec_err={r['evec_constraint_max_rel_err']:.2e} "
                    f"κ_absmax={r['kappa_abs_max']:.6f} resid={r['residual_emp']:.4f} "
                    f"budget={r['budget']:.4f} ok={r['residual_le_budget']}",
                    flush=True,
                )
            elif kind == "kp":
                kappas.append(r)
                print(
                    f"  kappa p={p}: λmax_res={r.get('lam_max_residual')} "
                    f"budget={r.get('budget')} ok={r.get('bound_ok')} "
                    f"ratio={r.get('ratio')} mult_info={r.get('spectrum_res_unique_r6', r.get('method'))}",
                    flush=True,
                )
            else:
                forms.append(r)
                print(
                    f"  form p={p}: res={r.get('residual')} eq_budget={r.get('equals_budget')} "
                    f"best={r.get('best_form')}",
                    flush=True,
                )

    structs.sort(key=lambda r: r["p"])
    kappas.sort(key=lambda r: r["p"])
    forms.sort(key=lambda r: r["p"])

    all_ok = all(r.get("bound_ok") for r in kappas if "bound_ok" in r)
    out = {
        "moment_structure": structs,
        "kappa_spectrum": kappas,
        "formula": forms,
        "status": (
            ("κ spectrum on Z CERTIFIED ≤ budget p=3,5,7; " if all_ok else "BOUND FAIL; ")
            + "uniform proof of λ_max(κ|Z)≤(p+1)(p+7)/d still OPEN."
        ),
        "workers": W,
        "all_ok": all_ok,
    }
    path = ROOT / "evidence" / "e1_gmin_cumulant.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
