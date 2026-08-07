#!/usr/bin/env python3
"""
Prop 15.207 — ker(Gsum)=scheme⊕cross reduces to G+≻0 on ++ zero-diag
forms; scheme ++ formula; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.202–206)
  Dual-eq empty if free-e max κ_e < 2−α on ker∩box.  sc:=scheme⊕cross ⊆ ker
  with dim n−1+d² (15.202).  ker=sc at p=3,5,7 by dim match; free-e then
  equals free-e over sc, which is <2−α at p=5,7 (15.202–203).  General
  ker=sc still open.  PSD form (15.204): κ∈ker ⇔ yᵀ(κ⊙C)y≡0 on Max±.

PROVED Max+-free
  A. Scheme ++ formula.  Let P₊=(I+C/p)/2, Vp isometry R^d→V₊=range(P₊),
     d=n/2.  For diagonal Df (f∈R^n) and A=Df C+C Df:
        Vpᵀ A Vp = 2p Vpᵀ Df Vp.
     Proof: C=p(P₊−P₋), so P₊ A P₊ = p P₊ Df P₊ + p P₊ Df P₊ = 2p P₊ Df P₊.  ∎

  B. Scheme ⊆ ker Q₊.  For y∈Max+ (Cy=py, y_i=±1) and sum f=0:
        u=Vpᵀ y,  M=2p Vpᵀ Df Vp,
        uᵀ M u = 2p ∑_i f_i y_i² = 2p ∑ f_i = 0.
     Hence scheme-image ⊆ {M∈Sym(d): uᵀMu=0 ∀Max+}=:ker Q₊.  ∎

  C. Orthogonality ⇔ constant diagonal.  With B=Vp M Vpᵀ,
        ⟨M, Vpᵀ Df Vp⟩_F = ∑_i f_i B_ii.
     So M ⊥ scheme-image (all sum-zero f) ⇔ B has constant diagonal.
     If also Tr M=0 then ∑ B_ii=Tr B=0 ⇒ B_ii≡0 (zero diagonal).  ∎

  D. Wick identity on zero-diag ++ forms.  Let B=P₊ B P₊ with zero diagonal,
     S=2P₊=I+C/p.  Then SB=2B, Tr(SB)=2 Tr B=0, Tr(SBSB)=4‖B‖_F², so the
     Gaussian/Wick second moment of q=yᵀBy is
        E_Wick[q²] = 2 Tr(SBSB)+(Tr SB)² = 8 ‖B‖_F².  ∎
     (Uses only S=2P₊ from π, Max+-free.)

  E. Reduction of ker=sc.  Write zero-diag A with yᵀAy=0 on Max± as
     A=A₊₊+A₋₋+A×.  Cross A× vanishes on V±.  A₊₊↔M∈ker Q₊, A₋₋↔ker Q₋.
     scheme provides matching f on both sides; cross free.
     ker Q₊=scheme-image (dim n−1) and ker Q₋=scheme-image with the same f
     ⇔ ker=sc.  By (B)(C): ker Q₊=scheme-image ⇔ the only zero-diag
     B=P₊BP₊ with yᵀBy=0 on Max+ is B=0
     ⇔ G₊ is positive definite on the subspace
        𝒲₊₊⁰ := { w∈R^{E} : W⊙C = P₊(W⊙C)P₊, diag(W⊙C)=0 }
     (because yᵀBy = 2∑_{i<j} B_ij y_i y_j and with W=B⊙C one has
     E₊[q²]=4 wᵀ G₊ w on that subspace).  ∎

CERTIFIED (not general)
  F. dim ker Q₊ = n−1 = rank scheme-image at p=3,5; residual outside scheme
     <1e−14 (15.204/ this module).
  G. G₊|𝒲₊₊⁰ has min eigenvalue 8 (p=3) and 40/13 (p=5), both >0; so
     G₊≻0 on 𝒲₊₊⁰ at p=3,5 and ker Q₊=scheme-image there.
  H. Same dim match ker=sc from Gsum nullity at p=3,5,7 (15.202).

OPEN (blocks residual i / predicate flip)
  I. Prove G₊ ≻ 0 on 𝒲₊₊⁰ for all primes p≥5 (or ker Q₊=scheme-image
     by another route).  Then ker=sc.  Combined with free-e over sc <2−α
     (e.g. S*_sc ≤ n−2, cost_D<2−α, or free-e ≤ r(p)·scheme with r≤3/2
     from 15.205), residual-(i) dual-eq closes for all p≥5.
  J. Alternatively still: |μ₄|≤M_cand, K₄≤Wick_hi.

Until I or J: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15207.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)


def theorem_scheme_plusplus_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For A=Df C+C Df: Vpᵀ A Vp = 2p Vpᵀ Df Vp.  "
            "Proof: C=p(P₊−P₋) ⇒ P₊ A P₊ = 2p P₊ Df P₊."
        ),
    }


def theorem_scheme_in_ker_Qplus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "sum f=0 ⇒ uᵀ(2p Vpᵀ Df Vp)u = 2p ∑ f_i y_i² = 0 on Max+.  "
            "scheme-image ⊆ ker Q₊."
        ),
    }


def theorem_orthogonality_constant_diagonal() -> dict:
    return {
        "proved": True,
        "theorem": (
            "⟨M, Vpᵀ Df Vp⟩_F = ∑ f_i B_ii for B=Vp M Vpᵀ.  "
            "M ⊥ scheme-image (all ∑f=0) ⇔ B has constant diagonal; "
            "with Tr M=0 ⇒ zero diagonal."
        ),
    }


def theorem_wick_identity_zero_diag_plusplus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "B=P₊BP₊, zero diagonal, S=2P₊: SB=2B, Tr(SB)=0, "
            "Tr(SBSB)=4‖B‖_F² ⇒ E_Wick[(yᵀBy)²]=8‖B‖_F²."
        ),
    }


def theorem_ker_sc_reduction() -> dict:
    return {
        "proved": True,
        "theorem": (
            "ker=sc for all p≥5 once G₊ ≻ 0 on 𝒲₊₊⁰ (and the Max− twin), "
            "by the decomposition A=A₊₊+A₋₋+A× and (A)–(E).  "
            "G₊≻0 on 𝒲₊₊⁰ ⇔ ker Q₊ = scheme-image."
        ),
        "open_step": "G₊ ≻ 0 on 𝒲₊₊⁰ for all primes p≥5",
    }


def certify_Gplus_pd_on_subspace(primes: list[int] | None = None) -> dict:
    """Numerical: min eig of G+ on 𝒲₊₊⁰ at small p."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power
    from e1_bitight_infeas import boolean_evecs

    if primes is None:
        primes = [3, 5]
    rows = {}
    all_ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(float)
        n = int(C.shape[0])
        d = n // 2
        P = (np.eye(n) + C / p) / 2
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        evals, evecs = np.linalg.eigh(C)
        order = np.argsort(evals)
        Vp = evecs[:, order[d:]]
        ii, jj = np.triu_indices(n, 1)
        F = (Mp[:, ii] * Mp[:, jj]) * C[ii, jj]
        Gplus = (F.T @ F) / len(Mp)
        basis_M = []
        for i in range(d):
            for j in range(i, d):
                M = np.zeros((d, d))
                if i == j:
                    M[i, i] = 1.0
                else:
                    M[i, j] = M[j, i] = 1.0
                basis_M.append(M)
        K = P**2
        cols = []
        for M in basis_M:
            B = Vp @ M @ Vp.T
            f = np.linalg.lstsq(K, np.diag(B), rcond=None)[0]
            B = B - P @ np.diag(f) @ P
            B = P @ B @ P
            np.fill_diagonal(B, 0.0)
            W = B * C
            cols.append(W[ii, jj])
        L = np.column_stack(cols)
        U, s, _ = np.linalg.svd(L, full_matrices=False)
        r = int(np.sum(s > 1e-8 * s[0]))
        U = U[:, :r]
        ev = np.linalg.eigvalsh(U.T @ Gplus @ U)
        min_ev = float(ev[0])
        rows[str(p)] = {
            "dim_subspace": r,
            "min_eig_Gplus": min_ev,
            "positive_definite": min_ev > 1e-8,
            "expected_dim": d * (d + 1) // 2 - n,
        }
        if min_ev <= 1e-8:
            all_ok = False
    return {
        "certified": all_ok,
        "proved_general": False,
        "by_p": rows,
        "note": (
            "G₊≻0 on 𝒲₊₊⁰ at p=3,5 (min eig 8 and 40/13).  "
            "General PD is the open step for ker=sc."
        ),
    }


def certify_ker_Qplus_eq_scheme(primes: list[int] | None = None) -> dict:
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power
    from e1_bitight_infeas import boolean_evecs

    if primes is None:
        primes = [3, 5]
    rows = {}
    all_ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(float)
        n = int(C.shape[0])
        d = n // 2
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        evals, evecs = np.linalg.eigh(C)
        order = np.argsort(evals)
        Vp = evecs[:, order[d:]]
        U = Mp @ Vp
        idx = [(i, j) for i in range(d) for j in range(i, d)]
        nsym = len(idx)
        A = np.zeros((len(U), nsym))
        for t, (i, j) in enumerate(idx):
            A[:, t] = U[:, i] ** 2 if i == j else 2 * U[:, i] * U[:, j]
        s = np.linalg.svd(A, compute_uv=False)
        rank = int(np.sum(s > 1e-8 * s[0]))
        ker_dim = nsym - rank
        # scheme rank
        Ms = []
        for v in range(n - 1):
            f = np.zeros(n)
            f[v] = 1.0
            f[n - 1] = -1.0
            M = 2 * p * (Vp.T @ np.diag(f) @ Vp)
            Ms.append(np.array([M[i, j] for i, j in idx]))
        Ms = np.column_stack(Ms)
        sM = np.linalg.svd(Ms, compute_uv=False)
        rank_sc = int(np.sum(sM > 1e-8 * sM[0]))
        # residual
        _, _, Vt = np.linalg.svd(A, full_matrices=True)
        null = Vt[rank:].T
        coef, _, _, _ = np.linalg.lstsq(Ms, null, rcond=None)
        resid = float(np.linalg.norm(null - Ms @ coef))
        ok = ker_dim == n - 1 and rank_sc == n - 1 and resid < 1e-10
        rows[str(p)] = {
            "ker_Qplus_dim": ker_dim,
            "scheme_rank": rank_sc,
            "n_minus_1": n - 1,
            "residual_outside_scheme": resid,
            "match": ok,
        }
        if not ok:
            all_ok = False
    return {
        "certified": all_ok,
        "proved_general": False,
        "by_p": rows,
    }


def free_e_bound_proved_general() -> bool:
    return False


def ker_sc_proved_general() -> bool:
    return False


def residual_i_closed_via_207() -> bool:
    return free_e_bound_proved_general() and ker_sc_proved_general()


def hinge_status_207() -> dict:
    return {
        "scheme_plusplus_formula": True,
        "wick_identity_zero_diag_plusplus": True,
        "ker_sc_reduction_to_Gplus_pd": True,
        "Gplus_pd_on_Wpp0_general": False,
        "ker_sc_general": False,
        "residual_i_closed": residual_i_closed_via_207(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove G₊≻0 on 𝒲₊₊⁰ for all p≥5 (⇒ ker=sc), then free-e/sc bound; "
            "or |μ|≤M_cand / K₄≤Wick_hi."
        ),
    }


def main() -> dict:
    a = theorem_scheme_plusplus_formula()
    b = theorem_scheme_in_ker_Qplus()
    c = theorem_orthogonality_constant_diagonal()
    d = theorem_wick_identity_zero_diag_plusplus()
    e = theorem_ker_sc_reduction()
    cert_g = certify_Gplus_pd_on_subspace()
    cert_k = certify_ker_Qplus_eq_scheme()
    out = {
        "title": (
            "Prop 15.207 ker=sc reduces to G+≻0 on ++ zero-diag; "
            "residual (i) OPEN"
        ),
        "proved": {
            "scheme_plusplus_formula": a["proved"],
            "scheme_in_ker_Qplus": b["proved"],
            "orthogonality_constant_diagonal": c["proved"],
            "wick_identity": d["proved"],
            "ker_sc_reduction": e["proved"],
            "Gplus_pd_general": False,
            "ker_sc_general": False,
            "residual_i_closed": residual_i_closed_via_207(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "scheme_plusplus": a,
            "scheme_in_ker": b,
            "ortho_diag": c,
            "wick": d,
            "reduction": e,
        },
        "certified": {
            "Gplus_pd_subspace": cert_g,
            "ker_Qplus_eq_scheme": cert_k,
        },
        "hinge": hinge_status_207(),
        "F3": (
            "Reduction of ker=sc to G₊≻0 on 𝒲₊₊⁰ is Max+-free algebra; "
            "G₊≻0 general remains open.  Do not flip residual/E1/L."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15207.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.207 ker=sc reduces to G+≻0 on ++ zero-diag forms")
    print(f"  scheme ++ formula: {a['proved']}")
    print(f"  Wick identity: {d['proved']}")
    print(f"  ker=sc reduction: {e['proved']}")
    print(f"  G+ PD cert p=3,5: {cert_g['certified']}")
    for pk, v in cert_g["by_p"].items():
        print(f"    p={pk}: dim={v['dim_subspace']} min_eig={v['min_eig_Gplus']:.6f}")
    print(f"  ker Q+=scheme cert: {cert_k['certified']}")
    print(f"  residual_i closed: {residual_i_closed_via_207()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
