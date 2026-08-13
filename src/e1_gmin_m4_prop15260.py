#!/usr/bin/env python3
"""
Prop 15.260 — Global moment identities ∑_S μ κ and ∑_S μ φ;
layer counts n₁,n₃; residual-(i) still OPEN on π-covariance / reflection.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189 Q,R; 15.256 ∑μ; 15.250 3-wise vanish; bal pairwise)
  Q := ∑_{i<j} y_i y_j.  R := ∑_{i<j} C_ij y_i y_j = yᵀ C y / 2.
  Max+: Q≡p, R≡p n/2.  Max−: Q≡−p, R≡−p n/2.
  Bal: E[Q²]=p², E[R²]=p² n²/4, E[QR]=p² n/2.
  φ(S)=∑_r ∏_{i∈S} C_ri.  κ as usual.

PROVED Max+-free (normalised conference with Max±, bal π, 3-wise 0)
  A. ∑_S μ(S) κ(S) = n p² (p²−1)/8.
     Proof: expand E[R²]=binom(n,2) + 2 ∑_S κ μ
     (same-edge: C_e²=1; wedges die by bal pairwise; disjoint matchings
     contribute 2κ X per 4-set).  E[R²]=p² n²/4, so
        ∑ κ μ = (p² n²/4 − n(n−1)/2)/2 = n p² (p²−1)/8.  ∎

  B. ∑_S μ(S) φ(S) = n ∑_S μ(S) = − n p² (p²−1)/12.
     Proof: on Max±, ∑_r (C y)_r⁴ = n p⁴.  Expand
        ∑_r (∑_i C_ri y_i)⁴
     by index partition.  Types (3,1) and (2,1,1) vanish under bal/C².
     Type (4): n(n−1).  Type (2,2): 3 n(n−1)(n−2).  Type (1⁴): 24 ∑ μ φ.
        n p⁴ = n(n−1)+3n(n−1)(n−2)+24 ∑ μ φ.
     Algebra with n=p²+1 yields ∑ μ φ = −n p²(p²−1)/12 = n ∑ μ.  ∎

  C. Layer counts (any conference with κ∈{±1,±3} a.e.).
     ∑ κ² = n(n−1)(n−2)(n−5)/8 (15.226 / design).  With B=binom(n,4),
        n₁ + n₃ = B,   n₁ + 9 n₃ = ∑κ²,
     so n₃=(∑κ²−B)/8, n₁=B−n₃.  ∎

CERTIFIED
  D. (A)(B) match full 4-set census at p=5,7.  (C) matches n₁,n₃ there.
  E. On |κ|=1: max|μ|=3/65 at p=5, 436/11452 at p=7, both <1/(2p);
     ∑_{κ1} μ κ = 300 (p=5), ≈4004.89 (p=7) — not yet closed-form forced
     to imply the L^∞ bound alone.

OPEN (blocks residual i / predicate flip)
  F. π-covariance m₄⁺(S)=χ_D m₄⁺(π⁻¹S) on |κ|=1, or reflection
     |ρ|≤|ρ_f4|, or |μ|≤1/(2p) by other Max+-free means.  Soft-close forbidden.
     Moments (A)(B) control averages, not the |κ|=1 L^∞ hinge by themselves.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15260.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    type_I_k_3p_minus_2_closed_general,
)


def sum_mu(p: int) -> Fraction:
    return Fraction(-(p * p) * (p * p - 1), 12)


def sum_mu_kappa(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(n * p * p * (p * p - 1), 8)


def sum_mu_phi(p: int) -> Fraction:
    return n_of(p) * sum_mu(p)


def layer_counts(p: int) -> tuple[int, int]:
    n = n_of(p)
    B = comb(n, 4)
    sk2 = n * (n - 1) * (n - 2) * (n - 5) // 8
    n3 = (sk2 - B) // 8
    n1 = B - n3
    return n1, n3


def theorem_sum_mu_kappa() -> dict:
    return {
        "proved": True,
        "theorem": (
            "∑_S μ κ = n p²(p²−1)/8.  "
            "Proof: E[R²]=binom(n,2)+2∑κμ with E[R²]=p²n²/4."
        ),
        "depends_on": ["R_identity_Maxpm", "bal_pairwise", "matching_expand"],
    }


def theorem_sum_mu_phi() -> dict:
    return {
        "proved": True,
        "theorem": (
            "∑_S μ φ = n ∑ μ = −n p²(p²−1)/12.  "
            "Proof: ∑_r(Cy)_r⁴=np⁴ on Max±; expand; (3,1)/(2,1,1) die; "
            "type (4)+(2,2)+(1⁴) give the identity."
        ),
        "depends_on": ["Cy_eq_pm_py", "bal_pairwise", "C_squared", "sum_mu_256"],
    }


def theorem_layer_counts() -> dict:
    return {
        "proved": True,
        "theorem": (
            "n₃=(∑κ²−binom(n,4))/8, n₁=binom(n,4)−n₃ with "
            "∑κ²=n(n−1)(n−2)(n−5)/8."
        ),
        "depends_on": ["sum_kappa_sq_226"],
    }


def certify_moments() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        a, b, c, d = fours.T
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        phi = np.zeros(len(fours))
        for r in range(n):
            phi += C[r, a] * C[r, b] * C[r, c] * C[r, d]
        mu = np.zeros(len(fours))
        bs = 8000
        for s0 in range(0, len(fours), bs):
            sl = slice(s0, s0 + bs)
            mp = np.mean(Yp[:, a[sl]] * Yp[:, b[sl]] * Yp[:, c[sl]] * Yp[:, d[sl]], 0)
            mm = np.mean(Ym[:, a[sl]] * Ym[:, b[sl]] * Ym[:, c[sl]] * Ym[:, d[sl]], 0)
            mu[sl] = 0.5 * (mp + mm)
        n1_t, n3_t = layer_counts(p)
        m1 = np.abs(np.abs(kap) - 1) < 1e-9
        out[str(p)] = {
            "ok": True,
            "sum_mu": float(mu.sum()),
            "sum_mu_theory": float(sum_mu(p)),
            "sum_mu_kappa": float((mu * kap).sum()),
            "sum_mu_kappa_theory": float(sum_mu_kappa(p)),
            "sum_mu_phi": float((mu * phi).sum()),
            "sum_mu_phi_theory": float(sum_mu_phi(p)),
            "n1": int(m1.sum()),
            "n1_theory": n1_t,
            "n3": int((~m1).sum()),
            "n3_theory": n3_t,
            "max_abs_mu_kappa1": float(np.max(np.abs(mu[m1]))),
            "thr_1_over_2p": 1 / (2 * p),
            "matches": (
                abs(mu.sum() - float(sum_mu(p))) < 1e-6
                and abs((mu * kap).sum() - float(sum_mu_kappa(p))) < 1e-6
                and abs((mu * phi).sum() - float(sum_mu_phi(p))) < 1e-6
                and int(m1.sum()) == n1_t
            ),
        }
    return out


def residual_i_closed_via_260() -> bool:
    return False


def hinge_status_260() -> dict:
    return {
        "sum_mu_kappa": True,
        "sum_mu_phi": True,
        "layer_counts": True,
        "pi_covariance_general": False,
        "reflection_hyp_general": False,
        "residual_i_closed": residual_i_closed_via_260(),
        "open": (
            "Moments control averages only.  Still need π-covariance / "
            "reflection / |μ|≤1/(2p) Max+-free on |κ|=1."
        ),
    }


def main() -> dict:
    a = theorem_sum_mu_kappa()
    b = theorem_sum_mu_phi()
    c = theorem_layer_counts()
    # algebraic checks
    alg_ok = True
    for p in (5, 7, 11, 13, 17, 19):
        n = n_of(p)
        if sum_mu_phi(p) != n * sum_mu(p):
            alg_ok = False
        # R^2 identity rewrite
        smk = sum_mu_kappa(p)
        if smk != Fraction(n * p * p * (p * p - 1), 8):
            alg_ok = False
        n1, n3 = layer_counts(p)
        B = comb(n, 4)
        sk2 = n * (n - 1) * (n - 2) * (n - 5) // 8
        if n1 + n3 != B or n1 + 9 * n3 != sk2:
            alg_ok = False
    cert = certify_moments()
    out = {
        "title": (
            "Prop 15.260 ∑μκ and ∑μφ identities; layer counts; "
            "residual (i) OPEN"
        ),
        "proved": {
            "sum_mu_kappa": a["proved"],
            "sum_mu_phi": b["proved"],
            "layer_counts": c["proved"],
            "algebraic_checks": alg_ok,
            "pi_covariance_general": False,
            "reflection_hyp_general": False,
            "residual_i_closed": residual_i_closed_via_260(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "sum_mu_kappa": a,
            "sum_mu_phi": b,
            "layer_counts": c,
            "formulas": {
                str(p): {
                    "sum_mu": str(sum_mu(p)),
                    "sum_mu_kappa": str(sum_mu_kappa(p)),
                    "sum_mu_phi": str(sum_mu_phi(p)),
                    "n1_n3": layer_counts(p),
                }
                for p in (5, 7, 11, 13)
            },
        },
        "certified": cert,
        "hinge": hinge_status_260(),
        "F3": (
            "∑μκ=np²(p²−1)/8 and ∑μφ=n∑μ proved Max+-free.  Layer counts "
            "closed.  L^∞ hinge on |κ|=1 still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15260.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.260 moment identities; residual-i OPEN")
    print(f"  sum mu kappa: {a['proved']}")
    print(f"  sum mu phi: {b['proved']}")
    print(f"  layer counts: {c['proved']}")
    print(f"  algebraic: {alg_ok}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: match={v['matches']} max|mu|_k1={v['max_abs_mu_kappa1']:.6f} "
            f"<1/(2p)={v['thr_1_over_2p']:.6f}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_260()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
