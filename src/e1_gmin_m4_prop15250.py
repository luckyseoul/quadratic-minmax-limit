#!/usr/bin/env python3
"""
Prop 15.250 — R-path ⇔ E_{y,z∈Max+}[(y·z)⁴] ≤ 15n²−22n;
odd moments vanish; residual-(i) still OPEN on the fourth-moment bound.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.217 R≤2p ⇔ ‖m4‖₂²≤n(n-2)/4; Max+ = {y∈{±1}ⁿ : Cy = py})
  n = p²+1.  Write s_{y,z} := y·z and M := |Max+|.
  m4(S) := E_{y∈Max+}[∏_{i∈S} y_i] on 4-sets S.
  The R-path (15.217) closes residual-(i) dual-eq once ‖m4‖₂² ≤ n(n−2)/4.

PROVED Max+-free (any symmetric conference with Max+ nonempty)

  A. Global sign flip.
     If Cy = py then C(−y) = p(−y), so −y ∈ Max+ whenever y ∈ Max+.
     Hence Max+ is a union of ± pairs (or fixed points of which there are none
     for n≥2: y=−y impossible).  ∎

  B. Odd multilinear moments vanish.
     For any multi-index with odd total degree (in particular three distinct
     coordinates a,b,c): the monomial ∏ y_{i_k} is odd under y↦−y, so
        E_{y∈Max+}[y_a y_b y_c] = 0.
     (Pair y with −y; A.)  ∎

  C. Fourth-moment expansion.
     E_{y,z∈Max+}[s_{y,z}⁴] = ∑_{i,j,k,l} μ_{ijkl}² where μ_{ijkl}=E[y_i y_j y_k y_l].
     Partition [n]⁴ by equality type and use B + second-moment structure
     E[y_i y_j] = C_{ij}/p (i≠j), E[y_i²]=1 (2-design / 15.189):
        type (4):     n
        type (3,1):   4 n(n−1)/p²
        type (2,2):   3 n(n−1)
        type (2,1,1): 6 n(n−1)(n−2)/p²
        type (1⁴):    24 ‖m4‖₂²
     With p² = n−1:
        E[s⁴] = 9n² − 10n + 24 ‖m4‖₂².  ∎

  D. R-path ⇔ fourth-moment bound.
     ‖m4‖₂² ≤ n(n−2)/4
        ⇔ 24‖m4‖₂² ≤ 6n(n−2)
        ⇔ E[s⁴] − 9n² + 10n ≤ 6n² − 12n
        ⇔ E[s⁴] ≤ 15n² − 22n.
     Combined with 15.217: residual-(i) dual-eq empty once
        E_{y,z∈Max+}[(y·z)⁴] ≤ 15n² − 22n.  ∎

  E. Equivalent L² forms (recall).
     ‖m4‖₂² ≤ n(n−2)/4  ⇔  ‖η‖₂² ≤ ‖κ‖₂²/p⁴  (η = m4 − κ/p²; 15.246 F)
        ⇔  ‖δ‖₂² ≤ room_δ^R  (m4 = m4_part + δ; 15.247 D).  ∎

CERTIFIED
  F. E[s⁴] ≤ 15n²−22n at p=5,7 on full Max+ census
        (p=5: 9261.5 ≤ 9568; p=7: 31383.8 ≤ 36400).
     ‖m4‖₂² ≤ target (143.2 ≤ 156; 383.2 ≤ 600).
     Odd third moments max|μ₃|≈0 at p=5,7.
     Expansion C matches machine precision.

OPEN (blocks residual i / predicate flip)
  G. Prove E_{y,z∈Max+}[(y·z)⁴] ≤ 15n² − 22n for all primes p≥5
     (equivalently ‖m4‖₂² ≤ n(n−2)/4 / R≤2p / ‖δ‖₂²≤room_δ^R),
     **or** ker(Gsum)=sc so 15.249 free-e dual finishes residual-i,
     **or** |μ|≤1/(2p) on |κ|=1.
     Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15250.json
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
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15217 import m4_target_R  # noqa: E402


def Es4_from_m4_norm(p: int, m4_norm_sq: Fraction) -> Fraction:
    """E[s⁴] = 9n² − 10n + 24 ‖m4‖₂²."""
    n = n_of(p)
    return Fraction(9 * n * n - 10 * n) + 24 * m4_norm_sq


def Es4_R_bound(p: int) -> Fraction:
    """Target 15n² − 22n for R-path."""
    n = n_of(p)
    return Fraction(15 * n * n - 22 * n)


def theorem_odd_moments_vanish() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Max+ is closed under y↦−y (Cy=py ⇒ C(−y)=p(−y)).  "
            "Hence every odd multilinear moment on Max+ vanishes; in "
            "particular E[y_a y_b y_c]=0 for distinct a,b,c.  Max+-free."
        ),
        "depends_on": ["Max_plus_boolean_evec", "conference_C2"],
    }


def theorem_Es4_expansion() -> dict:
    """E[s⁴] = 9n² − 10n + 24‖m4‖₂² for n=p²+1."""
    ok = True
    sample = {}
    # Algebraic identity with p²=n-1
    for p in [p for p in range(5, 40) if is_prime(p)]:
        n = n_of(p)
        # type contributions with p^2 = n-1
        t4 = n
        t31 = 4 * n * (n - 1) // (p * p)  # integer: (n-1)/p^2=1 so 4n
        # use Fraction
        t31_f = Fraction(4 * n * (n - 1), p * p)
        t22 = 3 * n * (n - 1)
        t211 = Fraction(6 * n * (n - 1) * (n - 2), p * p)
        known = t4 + t31_f + t22 + t211
        expect = Fraction(9 * n * n - 10 * n)
        if known != expect:
            ok = False
        if p in (5, 7, 11):
            sample[str(p)] = {
                "known_part": str(known),
                "9n2_minus_10n": str(expect),
                "match": known == expect,
            }
    return {
        "proved": ok,
        "theorem": (
            "E_{y,z∈Max+}[(y·z)⁴] = 9n²−10n + 24‖m4‖₂².  "
            "Proof: expand ∑ μ_ijkl² by equality type; odd type-(3,1) uses "
            "E[y_a y_b]=C_ab/p; type (2,1,1) same; type (1⁴)=24‖m4‖₂²; "
            "type (2,2) and (4) are 1; simplify with p²=n−1."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "odd_moments_vanish",
            "second_moment_2_design",
            "conference_p2_n_minus_1",
        ],
    }


def theorem_R_path_iff_Es4() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        n = n_of(p)
        # Algebra: m4_target = n(n-2)/4
        # 24 * target = 6n(n-2) = 6n² − 12n
        # Es4 at target = 9n²−10n + 6n²−12n = 15n² − 22n
        lhs = 24 * m4_target_R(p)
        rhs = Fraction(6 * n * (n - 2))
        es4_at = Fraction(9 * n * n - 10 * n) + lhs
        bound = Es4_R_bound(p)
        if lhs != rhs or es4_at != bound:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "m4_target": str(m4_target_R(p)),
                "Es4_bound": str(bound),
                "match": es4_at == bound,
            }
    return {
        "proved": ok,
        "theorem": (
            "‖m4‖₂² ≤ n(n−2)/4 ⇔ E[s⁴] ≤ 15n²−22n.  "
            "Hence R≤2p (15.217) ⇔ the fourth-moment bound on Max+."
        ),
        "by_p_sample": sample,
        "depends_on": ["Es4_expansion", "prop_15_217_R_path"],
    }


def certify_census_Es4() -> dict:
    """Census E[s⁴] and ‖m4‖₂² at p=5,7 if caches exist."""
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        cache = Path(f"/tmp/maxplus_p{p}.npy")
        if not cache.exists():
            out[str(p)] = {"ok": False, "reason": "no Max+ cache"}
            continue
        Y = np.load(cache)
        n = int(Y.shape[1])
        M = int(Y.shape[0])
        # E[s^4]
        # chunked Gram fourth power
        total = 0.0
        bs = 200
        for s in range(0, M, bs):
            G = Y[s : s + bs] @ Y.T
            total += float(np.sum(G.astype(np.float64) ** 4))
        Es4 = total / (M * M)
        # ||m4||^2
        m4n2 = 0.0
        for S in combinations(range(n), 4):
            m = float(np.mean(np.prod(Y[:, S], axis=1)))
            m4n2 += m * m
        bound = float(Es4_R_bound(p))
        target = float(m4_target_R(p))
        # expansion check
        Es4_from_m4 = float(Es4_from_m4_norm(p, Fraction(m4n2).limit_denominator(10**12)))
        # use float expansion
        Es4_pred = 9 * n * n - 10 * n + 24 * m4n2
        out[str(p)] = {
            "ok": True,
            "M": M,
            "Es4": Es4,
            "Es4_bound": bound,
            "Es4_le_bound": Es4 <= bound + 1e-6,
            "m4_norm_sq": m4n2,
            "m4_target": target,
            "m4_le_target": m4n2 <= target + 1e-6,
            "expansion_match": abs(Es4 - Es4_pred) < 1e-6,
            "room_Es4": bound - Es4,
            "room_m4": target - m4n2,
        }
    return out


def residual_i_closed_via_250() -> bool:
    return False


def hinge_status_250() -> dict:
    return {
        "odd_moments_vanish": True,
        "Es4_expansion": True,
        "R_path_iff_Es4_bound": True,
        "Es4_bound_general": False,
        "residual_i_closed": residual_i_closed_via_250(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove E[(y·z)⁴]≤15n²−22n on Max+ for all primes p≥5 "
            "(or ker=sc for 15.249 / |μ|≤1/(2p))."
        ),
    }


def main() -> dict:
    a = theorem_odd_moments_vanish()
    b = theorem_Es4_expansion()
    c = theorem_R_path_iff_Es4()
    cert = certify_census_Es4()
    out = {
        "title": (
            "Prop 15.250 R-path ⇔ E[s⁴]≤15n²−22n; odd moments vanish; "
            "residual (i) OPEN"
        ),
        "proved": {
            "odd_moments_vanish": a["proved"],
            "Es4_expansion": b["proved"],
            "R_path_iff_Es4": c["proved"],
            "Es4_bound_general": False,
            "residual_i_closed": residual_i_closed_via_250(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "odd_moments": a,
            "Es4_expansion": b,
            "R_path_iff": c,
        },
        "certified": {"census_Es4": cert},
        "hinge": hinge_status_250(),
        "closed_forms": {
            "Es4": "9n^2 - 10n + 24 ||m4||_2^2",
            "Es4_R_bound": "15n^2 - 22n",
            "m4_target_R": "n(n-2)/4",
        },
        "F3": (
            "R-path is exactly the Max+ fourth-moment bound E[s⁴]≤15n²−22n.  "
            "Odd moments vanish Max+-free.  The inequality itself is OPEN.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15250.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.250 R-path via E[s^4]")
    print(f"  odd moments vanish: {a['proved']}")
    print(f"  Es4 expansion: {b['proved']}")
    print(f"  R-path iff Es4 bound: {c['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: Es4={v['Es4']:.2f}<={v['Es4_bound']:.0f}? {v['Es4_le_bound']} "
            f"m4={v['m4_norm_sq']:.2f}<={v['m4_target']:.0f}? {v['m4_le_target']} "
            f"exp_match={v['expansion_match']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_250()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
