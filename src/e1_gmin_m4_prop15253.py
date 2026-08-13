#!/usr/bin/env python3
"""
Prop 15.253 — Wick residual of f4 beats L_abs; reflection criterion
|ρ|≤|ρ_f4| ⇔ μ∈[f4, 2κ/p²−f4] ⇒ residual-(i); census p=5,7; OPEN on
the reflection hypothesis.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176 Farkas, 15.205 L_abs, 15.191 f4, 15.252 ρ-path)
  ρ := μ − κ/p² (Wick residual of actual μ on 4-sets).
  f4 := (4κ − φ)/(p n),   ρ_f4 := f4 − κ/p².
  L_abs := (p−2)/(2p²).  On |κ|=1: |μ|≤1/p²+|ρ|, so |ρ|≤L_abs ⇒ |μ|≤1/(2p)
  (15.252 B).  Reflection of f4 over Wick: f4^♯ := 2κ/p² − f4.

PROVED Max+-free
  A. Closed form of ρ_f4 on |κ|=1.
        ρ_f4 = κ(4p − n)/(p² n) − φ/(p n)
             = κ(4p − p² − 1)/(p² n) − φ/(p n).  ∎

  B. |ρ_f4| ≤ L_abs for all primes p≥5 under |κ|=1, |φ|≤2(p−2).
     |4p−p²−1| = p²−4p+1 for p≥5, so
        |ρ_f4| ≤ (p²−4p+1)/(p² n) + 2(p−2)/(p n)
                = (3p² − 8p + 1)/(p² n).
     Need (3p²−8p+1)/(p²(p²+1)) ≤ (p−2)/(2p²)
     ⇔ 2(3p²−8p+1) ≤ (p−2)(p²+1)
     ⇔ 0 ≤ g(p) := p³ − 8p² + 17p − 4.
     g(5)=6>0; g'(p)=3p²−16p+17, g'(5)=12>0, and the parabola 3p²−16p+17
     opens upward with vertex at p=8/3<5, so g'>0 on [5,∞) and g>0.  ∎

  C. Reflection criterion (sufficient for residual-(i)).
     On |κ|=1 the following are equivalent:
        (i)  |ρ| ≤ |ρ_f4|
        (ii) μ lies in the closed interval with endpoints f4 and f4^♯
             (i.e. between f4 and its reflection across κ/p²).
     If either holds for every |κ|=1 four-set and all primes p≥5, then
     |ρ|≤|ρ_f4|≤L_abs (by B) ⇒ |μ|≤1/(2p) ⇒ Gsum≥−1/p ⇒ residual-(i)
     Farkas (15.176).  ∎

  D. Hull / triple majorant (structure).
     If μ ∈ conv{0, μ_part, f4} (hull, 15.241), then
        ρ ∈ conv{−κ/p², ρ_part, ρ_f4},
     so |ρ| ≤ max(1/p², |ρ_part|, |ρ_f4|).
     Each of 1/p², maj|ρ_part|, maj|ρ_f4| is ≤ L_abs for p≥5
     (1/p²≤L_abs ⇔ p≥4; |ρ_part| maj uses the same cubic as maj≤1/(2p);
     |ρ_f4| by B).  Hence hull ⇒ residual-(i) via the ρ-path as well.  ∎

CERTIFIED (Max± census; not general)
  E. p=5: μ=f4 on all |κ|=1 ⇒ |ρ|=|ρ_f4| (reflection equality at the
     f4 endpoint); max|ρ|=18/325 < 3/50 = L_abs.
  F. p=7: |ρ|≤|ρ_f4| on all 176400 |κ|=1 sets (0 viol); μ ∈ [f4,f4^♯]
     (0 viol); Ext ∈ [Ext_f4, Ext_f4^♯] (0 viol); max|ρ|≈0.0177 < 0.0261
     = max|ρ_f4| < L_abs=5/98.

OPEN (blocks residual i / predicate flip)
  G. Prove the reflection hypothesis |ρ|≤|ρ_f4| (equivalently μ∈[f4,f4^♯])
     on every |κ|=1 four-set for all primes p≥5 — Max+-free —
     **or** hull, **or** |Ext|≤2p²−4p+6, **or** Es4 / ker=sc.
     Soft-close forbidden.  Note: pure |μ|≤|f4| is strictly stronger and
     **false** at p=7; the Wick-centred reflection interval is the right
     target (census-stable, implies L_abs).

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15253.json
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
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import theorem_phi_bound_general  # noqa: E402
from e1_gmin_m4_prop15205 import L_abs  # noqa: E402


def rho_f4_formula(p: int, kappa: int, phi: int) -> Fraction:
    """ρ_f4 = f4 − κ/p² on a labelled 4-set."""
    n = n_of(p)
    return Fraction(4 * kappa - phi, p * n) - Fraction(kappa, p * p)


def f4_sharp(p: int, kappa: int, phi: int) -> Fraction:
    """Reflection of f4 across κ/p²."""
    return Fraction(2 * kappa, p * p) - Fraction(4 * kappa - phi, p * n_of(p))


def rho_f4_majorant(p: int) -> Fraction:
    """(3p² − 8p + 1)/(p² n) under |κ|=1, |φ|≤2(p−2), p≥5."""
    n = n_of(p)
    return Fraction(3 * p * p - 8 * p + 1, p * p * n)


def g_cubic(p: int) -> int:
    """p³ − 8p² + 17p − 4 (sign of L_abs − rho_f4 maj)."""
    return p**3 - 8 * p * p + 17 * p - 4


def theorem_rho_f4_closed_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1: ρ_f4 = κ(4p−p²−1)/(p² n) − φ/(p n).  "
            "Algebraic rearrangement of f4−κ/p²."
        ),
        "depends_on": ["f4_definition", "n_eq_p2_plus_1"],
    }


def theorem_rho_f4_le_L_abs() -> dict:
    primes = [p for p in range(5, 120) if is_prime(p)]
    ok = True
    sample = {}
    # g(5)>0 and g'>0 on [5,∞)
    g5 = g_cubic(5) > 0
    # g'(p)=3p²-16p+17; vertex 8/3; value at 5
    gp5 = 3 * 25 - 16 * 5 + 17  # 12
    for p in primes:
        maj = rho_f4_majorant(p)
        L = L_abs(p)
        beats = maj <= L
        # also check explicit max over discrete phi
        bound_phi = 2 * (p - 2)
        phis = list(range(-bound_phi, bound_phi + 1, 4))
        max_ex = max(
            abs(rho_f4_formula(p, k, ph)) for k in (-1, 1) for ph in phis
        )
        if not beats or max_ex > maj or g_cubic(p) <= 0:
            ok = False
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "rho_f4_maj": str(maj),
                "L_abs": str(L),
                "maj_le_L": beats,
                "max_explicit": str(max_ex),
                "explicit_le_maj": max_ex <= maj,
                "g_cubic": g_cubic(p),
            }
    return {
        "proved": ok and g5 and gp5 > 0 and theorem_phi_bound_general()["proved"],
        "theorem": (
            "On |κ|=1 with |φ|≤2(p−2): |ρ_f4|≤(3p²−8p+1)/(p²n)≤L_abs "
            "for all primes p≥5 (g(p)=p³−8p²+17p−4>0; g(5)=6, g'>0 on [5,∞))."
        ),
        "by_p_sample": sample,
        "g5": g5,
        "gp5": gp5,
        "depends_on": ["phi_bound_15_186", "L_abs_15_205"],
    }


def theorem_reflection_criterion() -> dict:
    # equivalence algebra: |mu - c| <= |f4 - c| iff mu between f4 and 2c-f4
    # check on samples
    ok = True
    for p in (5, 7, 11):
        for k in (-1, 1):
            for ph in range(-2 * (p - 2), 2 * (p - 2) + 1, 4):
                c = Fraction(k, p * p)
                f = Fraction(4 * k - ph, p * n_of(p))
                sharp = 2 * c - f
                lo, hi = min(f, sharp), max(f, sharp)
                # test midpoints and endpoints
                for t in (0, Fraction(1, 2), 1, Fraction(1, 3)):
                    mu = (1 - t) * f + t * sharp
                    rho = mu - c
                    rho_f = f - c
                    if abs(rho) > abs(rho_f) + Fraction(1, 10**9):
                        ok = False
                # outside should fail
                outside = hi + Fraction(1, 1000)
                if abs(outside - c) <= abs(f - c):
                    # may not be outside ball; skip
                    pass
    return {
        "proved": ok,  # equivalence is pure algebra; "proved" = criterion sound
        "theorem": (
            "On |κ|=1: |ρ|≤|ρ_f4| ⇔ μ ∈ [min(f4,f4^♯), max(f4,f4^♯)] with "
            "f4^♯=2κ/p²−f4.  Combined with (B): the reflection hypothesis "
            "⇒ |ρ|≤L_abs ⇒ |μ|≤1/(2p) ⇒ residual-(i) Farkas."
        ),
        "hypothesis_proved_general": False,
        "depends_on": ["theorem_rho_f4_le_L_abs", "prop_15_252_rho_path"],
    }


def theorem_hull_implies_Labs() -> dict:
    # 1/p^2 <= L_abs for p>=4; rho_part maj; rho_f4
    ok = True
    sample = {}
    for p in [p for p in range(5, 60) if is_prime(p)]:
        one = Fraction(1, p * p)
        L = L_abs(p)
        # |rho_part| <= 4(p-1)/(p^2 (p^2-5))
        rp = Fraction(4 * (p - 1), p * p * (p * p - 5))
        rf = rho_f4_majorant(p)
        triple = max(one, rp, rf)
        if not (one <= L and rp <= L and rf <= L and triple <= L):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "one_over_p2": str(one),
                "rho_part_maj": str(rp),
                "rho_f4_maj": str(rf),
                "triple_le_L": triple <= L,
            }
    return {
        "proved": ok,
        "theorem": (
            "Hull μ∈conv{0,μ_part,f4} on |κ|=1 ⇒ |ρ|≤max(1/p²,|ρ_part|,|ρ_f4|) "
            "≤L_abs for all primes p≥5 (each term ≤L_abs Max+-free).  "
            "Hence hull ⇒ residual-(i) via the ρ-path."
        ),
        "by_p_sample": sample,
        "depends_on": ["hull_15_241", "rho_part_maj", "theorem_rho_f4_le_L_abs"],
    }


def certify_reflection_census() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing Max± cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        a, b, c, d = fours.T
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        ids = np.where(np.abs(np.abs(kap) - 1) < 1e-9)[0]
        aa, bb, cc, dd = a[ids], b[ids], c[ids], d[ids]
        phi = np.zeros(len(ids))
        for r in range(n):
            phi += C[r, aa] * C[r, bb] * C[r, cc] * C[r, dd]
        kaps = kap[ids]
        mu = np.empty(len(ids))
        chunk = 4000
        for s0 in range(0, len(ids), chunk):
            sl = slice(s0, s0 + chunk)
            mp = np.mean(
                Yp[:, aa[sl]] * Yp[:, bb[sl]] * Yp[:, cc[sl]] * Yp[:, dd[sl]],
                axis=0,
            )
            mm = np.mean(
                Ym[:, aa[sl]] * Ym[:, bb[sl]] * Ym[:, cc[sl]] * Ym[:, dd[sl]],
                axis=0,
            )
            mu[sl] = 0.5 * (mp + mm)
        f4 = (4 * kaps - phi) / (p * n)
        sharp = 2 * kaps / (p * p) - f4
        rho = mu - kaps / (p * p)
        rho_f = f4 - kaps / (p * p)
        lo, hi = np.minimum(f4, sharp), np.maximum(f4, sharp)
        viol = int(np.sum((mu < lo - 1e-12) | (mu > hi + 1e-12)))
        viol_rho = int(np.sum(np.abs(rho) > np.abs(rho_f) + 1e-12))
        L = float(L_abs(p))
        out[str(p)] = {
            "ok": True,
            "n_kappa1": int(len(ids)),
            "reflection_viol": viol,
            "rho_le_rho_f4_viol": viol_rho,
            "max_abs_rho": float(np.max(np.abs(rho))),
            "max_abs_rho_f4": float(np.max(np.abs(rho_f))),
            "L_abs": L,
            "rho_le_L_abs": float(np.max(np.abs(rho))) <= L + 1e-12,
            "rho_f4_le_L_abs": float(np.max(np.abs(rho_f))) <= L + 1e-12,
            "reflection_ok": viol == 0,
        }
    return out


def residual_i_closed_via_253() -> bool:
    return False


def hinge_status_253() -> dict:
    return {
        "rho_f4_le_L_abs": True,
        "reflection_criterion": True,
        "reflection_hyp_general": False,
        "hull_implies_Labs": True,
        "residual_i_closed": residual_i_closed_via_253(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove |ρ|≤|ρ_f4| (μ∈[f4,f4^♯]) on |κ|=1 for all primes p≥5, "
            "or hull / Ext maj / Es4 / ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_rho_f4_closed_form()
    b = theorem_rho_f4_le_L_abs()
    c = theorem_reflection_criterion()
    d = theorem_hull_implies_Labs()
    cert = certify_reflection_census()
    out = {
        "title": (
            "Prop 15.253 |ρ_f4|≤L_abs; Wick-reflection criterion; "
            "residual (i) OPEN on reflection hyp"
        ),
        "proved": {
            "rho_f4_closed_form": a["proved"],
            "rho_f4_le_L_abs": b["proved"],
            "reflection_criterion": c["proved"],
            "reflection_hyp_general": False,
            "hull_implies_Labs": d["proved"],
            "residual_i_closed": residual_i_closed_via_253(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "rho_f4_form": a,
            "rho_f4_le_L_abs": b,
            "reflection_criterion": c,
            "hull_implies_Labs": d,
        },
        "certified": {"reflection_census": cert},
        "hinge": hinge_status_253(),
        "F3": (
            "|ρ_f4|≤L_abs proved Max+-free for p≥5.  Reflection hyp "
            "|ρ|≤|ρ_f4| (⇔ μ∈[f4,f4^♯]) would close residual-i; census "
            "holds at p=5,7 with 0 viol.  Pure |μ|≤|f4| is false at p=7.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15253.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.253 |rho_f4|<=L_abs; reflection criterion; residual-i OPEN")
    print(f"  rho_f4 form: {a['proved']}")
    print(f"  rho_f4 le L_abs: {b['proved']}")
    print(f"  reflection criterion: {c['proved']}")
    print(f"  hull => L_abs: {d['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: refl_viol={v['reflection_viol']} "
            f"rho_le_rho_f4_viol={v['rho_le_rho_f4_viol']} "
            f"max|rho|={v['max_abs_rho']:.6f} max|rho_f4|={v['max_abs_rho_f4']:.6f} "
            f"L={v['L_abs']:.6f}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_253()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
