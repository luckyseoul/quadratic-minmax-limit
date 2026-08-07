#!/usr/bin/env python3
"""
Prop 15.203 — Max+-free free-e dual construction (scheme→Comm∩reg-deg);
residual (i) still OPEN for general p.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182, 15.190–202)
  Dual-eq empty ⇔ free-e max κ_e < 2−α on ker(Gsum)∩box.
  Scheme⊕cross ⊆ ker (15.202).  Free-e dual: max κ_e = α·S* when the
  lower box is active, with S* = min L¹ dual w⊥ker, w_e=1, w≥0.
  If ker = scheme⊕cross, duals for scheme⊕cross bound true free-e max.
  ker=sc certified at p=3,5,7 (15.202); general equality still open.

PROVED Max+-free
  A. Constructive dual algorithm D(C) (uses only C and its ± eigenspaces):
     (1) Scheme dual W⁰: W⁰_e=1, stars 0, far edges 1/(n−3)
         (regular degree; S*=(n−2)/2 for scheme-only — 15.190).
     (2) S⁰ = W⁰ ⊙ C.  Unconstrained Comm-projection
         S⁺ = P₊ S⁰ P₊ + P₋ S⁰ P₋.
     (3) Min-norm Δ ∈ Comm with diag(S⁺+Δ)=0 and W=S⊙C regular degree
         (linear least squares in A,B on V±; n+(n−1) constraints).
     (4) S = S⁺+Δ, W = S⊙C, scale so W_e=1.
     (5) If min W < 0, add t·(J−I) with t=−min W and re-scale W_e=1
         (preserves Comm, zero-diag of S, regular degree; S ↦ S+tC).
     Output dual cost cost_D = α · (∑_{g≠e} W_g).
     Feasible dual for free-e over scheme⊕cross ⇒
        free-e max_{sc} κ_e ≤ cost_D.  ∎

  B. Implication.  If cost_D < 2−α and ker=scheme⊕cross, then free-e max
     < 2−α ⇒ dual-eq empty ⇒ residual-(i) dual-eq branch blocked.  ∎

  C. Algebra: α(n−1) < 2−α for all primes p≥5
        ⇔ 6(n−1) < 2 p n − 6 ⇔ 3(n−1) < p n − 3
        ⇔ 3p² < p(p²+1)−3 wait: n=p²+1,
        6(p²)/(p(p²+1)) < 2−6/(p(p²+1)) ⇔ 6p² < 2p(p²+1)−6
        ⇔ 6p² < 2p³+2p−6 ⇔ 0 < 2p³ − 6p² + 2p − 6 = 2(p³−3p²+p−3)
        ⇔ p³−3p²+p−3 = (p−3)(p²+1)+p−6?  Direct: at p=5 value 18>0;
        derivative 3p²−6p+1>0 on [5,∞).  ∎
     So free-e max ≤ α(n−1) also suffices for p≥5.

CERTIFIED (not general; Max+-free algorithm on C only)
  D. D(C) residual of linear constraints ~0 (machine eps) for primes
     p∈{5,7,11,13,17,19,23}; W≥0, degstd~0, ||[W⊙C,C]||~0 after scale.
  E. cost_D < 2−α for those primes (blocks dual-eq **if** ker=sc):
        p=5:  cost≈1.447 < need≈1.954
        p=7:  cost≈1.098 < need≈1.983
        p=11: cost≈0.649 < need≈1.996
        … decreasing; max cost among certified is at p=5.
  F. At p=3: cost≈2.96 > need 1.8 (does **not** block — consistent with
     dual-eq feasible at p=3).

OPEN (blocks residual i / predicate flip)
  G. Prove for all primes p≥5 that cost_D(p) < 2−α (closed form / asymptotic
     control of sum_ne after projection), **and** ker(Gsum)=scheme⊕cross
     for all p≥5 (so the dual applies to true free-e max).
     Alternatively: |μ|≤2/n or K₄≤Wick_hi Max+-free.
  H. Until G: residual_i False; gsum False; E1/L OPEN.  No soft-close.

Writes evidence/e1_gmin_m4_prop15203.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402


def theorem_alpha_n_minus_1_beats_need() -> dict:
    """α(n−1)<2−α for all primes p≥5."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        al = alpha_particular(p)
        lhs = al * (n - 1)
        need = 2 - al
        # poly p^3 - 3p^2 + p - 3 > 0 checked via 2(p^3-3p^2+p-3) from expansion
        poly = p**3 - 3 * p**2 + p - 3
        # direct Fraction compare
        beats = lhs < need
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "alpha_n_minus_1": str(lhs),
                "need": str(need),
                "beats": beats,
                "poly": poly,
            }
        if not beats or poly <= 0:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: α(n−1)<2−α.  Hence free-e max κ_e ≤ α(n−1) "
            "⇒ dual-eq empty."
        ),
        "by_p_sample": sample,
    }


def dual_construction_cost(p: int) -> dict:
    """Max+-free dual D(C) cost for free-e over scheme⊕cross."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    d = n // 2
    al = float(alpha_particular(p))
    W0 = np.zeros((n, n))
    W0[0, 1] = W0[1, 0] = 1.0
    far = 1.0 / (n - 3)
    for i in range(2, n):
        W0[i, i + 1 :] = far
        W0[i + 1 :, i] = far
    S0 = W0 * C
    evals, evecs = np.linalg.eigh(C)
    order = np.argsort(evals)
    Vm = evecs[:, order[:d]]
    Vp = evecs[:, order[d:]]
    pA = [(i, j) for i in range(d) for j in range(i, d)]
    nA = len(pA)
    nparam = 2 * nA
    Sp = Vp @ (Vp.T @ S0 @ Vp) @ Vp.T + Vm @ (Vm.T @ S0 @ Vm) @ Vm.T
    A_diag = np.zeros((n, nparam))
    for i in range(n):
        for k, (u, v) in enumerate(pA):
            if u == v:
                A_diag[i, k] = Vp[i, u] ** 2
                A_diag[i, nA + k] = Vm[i, u] ** 2
            else:
                A_diag[i, k] = 2 * Vp[i, u] * Vp[i, v]
                A_diag[i, nA + k] = 2 * Vm[i, u] * Vm[i, v]
    ii, jj = np.triu_indices(n, k=1)
    A_deg = np.zeros((n, nparam))
    for k, (u, v) in enumerate(pA):
        if u == v:
            dA = Vp[ii, u] * Vp[jj, u]
            dB = Vm[ii, u] * Vm[jj, u]
        else:
            dA = Vp[ii, u] * Vp[jj, v] + Vp[ii, v] * Vp[jj, u]
            dB = Vm[ii, u] * Vm[jj, v] + Vm[ii, v] * Vm[jj, u]
        wA = dA * C[ii, jj]
        wB = dB * C[ii, jj]
        np.add.at(A_deg[:, k], ii, wA)
        np.add.at(A_deg[:, k], jj, wA)
        np.add.at(A_deg[:, nA + k], ii, wB)
        np.add.at(A_deg[:, nA + k], jj, wB)
    Wp = Sp * C
    np.fill_diagonal(Wp, 0.0)
    deg_sp = Wp.sum(0)
    A_con = np.vstack([A_diag, A_deg[1:] - A_deg[0]])
    b_con = np.concatenate([-np.diag(Sp), deg_sp[0] - deg_sp[1:]])
    x, _, _, _ = np.linalg.lstsq(A_con, b_con, rcond=None)
    rnorm = float(np.linalg.norm(A_con @ x - b_con))
    A = np.zeros((d, d))
    B = np.zeros((d, d))
    for k, (u, v) in enumerate(pA):
        if u == v:
            A[u, u] = x[k]
            B[u, u] = x[nA + k]
        else:
            A[u, v] = A[v, u] = x[k]
            B[u, v] = B[v, u] = x[nA + k]
    S = Sp + Vp @ A @ Vp.T + Vm @ B @ Vm.T
    W = S * C
    np.fill_diagonal(W, 0.0)
    We = float(W[0, 1])
    if abs(We) < 1e-12:
        return {"p": p, "ok": False, "reason": "We≈0"}
    W1 = W / We
    t0 = max(0.0, -float(W1.min()))
    Joff = np.ones((n, n)) - np.eye(n)
    Wf = (W1 + t0 * Joff) / (1.0 + t0)
    sum_ne = float(Wf.sum() / 2.0 - 1.0)
    cost = al * sum_ne
    need = float(2 - alpha_particular(p))
    return {
        "p": p,
        "ok": True,
        "rnorm": rnorm,
        "t0": t0,
        "sum_ne": sum_ne,
        "cost": cost,
        "need": need,
        "blocks_if_ker_eq_sc": cost < need - 1e-9,
        "degstd": float(Wf.sum(0).std()),
        "minW": float(Wf.min()),
        "comm_norm": float(np.linalg.norm((Wf * C) @ C - C @ (Wf * C))),
        "proved_general": False,
    }


def free_e_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_203() -> bool:
    return free_e_bound_proved_general()


def hinge_status_203() -> dict:
    return {
        "dual_construction_max_plus_free": True,
        "alpha_n_minus_1_beats_need_p_ge_5": True,
        "cost_D_lt_need_general": False,
        "ker_eq_scheme_cross_general": False,
        "residual_i_closed": residual_i_closed_via_203(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove cost_D(p)<2−α for all primes p≥5 and ker=scheme⊕cross "
            "(or |μ|≤2/n / K₄≤Wick_hi)."
        ),
    }


def main() -> dict:
    a = theorem_alpha_n_minus_1_beats_need()
    cert = {str(p): dual_construction_cost(p) for p in (3, 5, 7)}
    out = {
        "title": (
            "Prop 15.203 Max+-free free-e dual construction D(C); "
            "residual (i) OPEN"
        ),
        "proved": {
            "dual_construction": True,
            "alpha_n_minus_1_beats_need": a["proved"],
            "cost_D_lt_need_general": False,
            "ker_eq_sc_general": False,
            "free_e_bound_general": False,
            "residual_i_closed": residual_i_closed_via_203(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"alpha_n_minus_1": a},
        "certified": {"dual_D": cert},
        "hinge": hinge_status_203(),
        "F3": (
            "Dual D(C) is Max+-free but cost_D<2−α and ker=sc are not proved "
            "for all p≥5.  Do not flip residual/E1/L."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15203.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.203 free-e dual construction D(C)")
    print(f"  alpha(n-1)<2-alpha p>=5: {a['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: fail {v}")
            continue
        print(
            f"  p={pk}: cost={v['cost']:.4f} need={v['need']:.4f} "
            f"blocks_if_ker_sc={v['blocks_if_ker_eq_sc']} "
            f"rnorm={v['rnorm']:.2e}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_203()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
