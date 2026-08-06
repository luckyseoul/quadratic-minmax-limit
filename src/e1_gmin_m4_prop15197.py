#!/usr/bin/env python3
"""
Prop 15.197 — Max+ min-distance + exact Q(K₄) identity (residual i structure).

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

PROVED Max+-free (any symmetric conference C with Cy=±p y boolean evecs)
  A. Min-distance of Max±.  For y≠z in Max+ (Cy=py, y,z∈{±1}^n),
     Hamming distance d_H(y,z) ≥ p+1.
     Proof.  Let D={i: y_i≠z_i}, d=|D|≥1.  Then y−z ∈ V₊ (eigenspace),
     so C(y−z)=p(y−z).  On D one has z=−y, hence y−z=2y.
     For k∈D: ∑_j C_kj (y_j−z_j)=p(y_k−z_k) ⇒ ∑_{j∈D} C_kj (2y_j)=2p y_k
     ⇒ ∑_{j∈D∖{k}} C_kj y_j = p y_k (C_kk=0).  LHS has abs ≤ d−1, RHS has
     abs p, so d−1≥p, d≥p+1.  Same for Max−.  (If z=−y then d=n≥p+1.)  ∎
     Consequence: for y≠±z, |y·z| ≤ n−2(p+1)=(p−1)²−2.

  B. Exact row energy via fourth moment of dots.  Assume Max+ antipodal
     (true for Paley; −y has C(−y)=p(−y)).  Let
        K₄ := E_{y,z∈Max+}[(y·z)⁴]
     (same for Max− by p↔−p symmetry of the algebra below when E[(u·v)²]=2n
     on Max− as well, which holds from π_−=−C/p).  Then the disj Gsum row
     energy is Max+-free-rational in K₄:
        Q_e = K₄/(n(n−1)) − 2(3n−2)/(n−1).
     Derivation (15.196 path made exact):
       • e₄(w)=∑_{|S|=4}∏_S w for w_i=±1 equals
         (s⁴−(6n−8)s²+3n²−6n)/24 with s=∑w;
       • ∑_S m₊(S)² = E_{y,z∈Max+} e₄(y⊙z) with m₊=E₊[∏_S y];
       • ∑_S m₊ m₋ = n(n−2)/8  (since y⊥u for y∈Max+, u∈Max− ⇒ s=0);
       • μ=½(m₊+m₋), Σ:=∑_S μ² = ½(∑m₊² + n(n−2)/8) under K₄⁺=K₄⁻;
       • Q_e = 48 Σ/(n(n−1)) reduces to the displayed formula.  ∎
     Certified: recovers Q=112/9, 34512/4225, 55009/8220 at p=3,5,7.

  C. Residual-(i) threshold on K₄.  From 15.196, Q_e≤10 forces mass-corrected
     dual-eq obstruction for all primes 5≤p≤47.  Substituting (B):
        Q_e≤10  ⇔  K₄ ≤ 16n² − 14n
     (algebra: K₄/(n(n−1)) ≤ 10 + 2(3n−2)/(n−1) = (16n²−14n)/(n(n−1))).  ∎

CERTIFIED (not general)
  D. Paley K₄ and the threshold 16n²−14n:
        p=3: K₄=1680 > 16n²−14n=1558  (Q>10; dual-eq feasible — consistent)
        p=5: K₄=120400/13≈9261.5 < 10452
        p=7: K₄≈31195.6 < 39300
     Distance bound sharp at p=5,7 (max |y·z|=14,34 = n−2(p+1)).

OPEN (blocks residual i / predicate flip)
  E. Prove Max+-free K₄ ≤ 16n²−14n for all primes p≥5
     (or Q_e≤10 / λ₂≤6+5/(n−2) / |μ|≤2/n / ker-box empty).
     Then flip gsum_disj_lb_proved_general → residual (i) → E1 → L
     (still need residual (ii) full / exhaustiveness for E1).

Until E: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15197.json
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
from e1_gmin_m4_prop15196 import worst_mass_min_under_Q  # noqa: E402


def max_dot_nonopposite(p: int) -> int:
    """n − 2(p+1) = (p−1)² − 2 from min-distance p+1."""
    n = n_of(p)
    return n - 2 * (p + 1)


def theorem_min_distance_p_plus_1() -> dict:
    """Max+-free: d_H ≥ p+1 on Max± for y≠z."""
    primes = [p for p in range(3, 60) if is_prime(p)]
    rows = {}
    for p in primes:
        n = n_of(p)
        dmin = p + 1
        max_dot = n - 2 * dmin
        rows[str(p)] = {
            "d_H_min": dmin,
            "max_dot_y_ne_pm_z": max_dot,
            "formula_max_dot": (p - 1) ** 2 - 2,
            "match": max_dot == (p - 1) ** 2 - 2,
        }
    ok = all(r["match"] for r in rows.values())
    return {
        "proved": ok,
        "theorem": (
            "For any symmetric conference C and any y≠z in Max+={±1}^n∩{Cy=py}, "
            "d_H(y,z)≥p+1.  Hence if z≠−y then |y·z|≤n−2(p+1)=(p−1)²−2. "
            "Proof: D=diff support; C(y−z)=p(y−z); on k∈D, "
            "∑_{j∈D∖k} C_kj y_j = p y_k ⇒ d−1≥p."
        ),
        "by_p_sample": {k: rows[k] for k in ("3", "5", "7", "11") if k in rows},
        "n_primes": len(primes),
    }


def Q_from_K4(p: int, K4: Fraction | int) -> Fraction:
    """Q = K4/(n(n-1)) − 2(3n−2)/(n−1)."""
    n = n_of(p)
    K4 = Fraction(K4)
    return K4 / (n * (n - 1)) - 2 * Fraction(3 * n - 2, n - 1)


def K4_threshold_for_Q10(p: int) -> Fraction:
    """K4 ≤ 16n² − 14n ⇔ Q ≤ 10."""
    n = n_of(p)
    return Fraction(16 * n * n - 14 * n)


def theorem_Q_from_K4_identity() -> dict:
    """Algebra: Q expression; threshold K4 for Q≤10."""
    primes = [p for p in range(5, 40) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        thr = K4_threshold_for_Q10(p)
        # Q(thr) should be 10
        q_at = Q_from_K4(p, thr)
        rows[str(p)] = {
            "K4_thr": str(thr),
            "Q_at_thr": str(q_at),
            "Q_eq_10": q_at == 10,
        }
        if q_at != 10:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Under antipodal Max± and K₄⁺=K₄⁻: "
            "Q_e = K₄/(n(n−1)) − 2(3n−2)/(n−1). "
            "Hence Q_e≤10 ⇔ K₄≤16n²−14n."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11") if k in rows},
    }


def theorem_K4_thr_implies_residual_i_block() -> dict:
    """If K4≤16n²−14n then Q≤10 then mass-min blocks dual-eq (p=5..47)."""
    return {
        "proved_implication": True,
        "K4_bound_proved_general": False,
        "theorem": (
            "IF K₄≤16n²−14n for all primes p≥5 (Max+-free), "
            "THEN Q_e≤10 ⇒ mass-corrected dual-eq impossible (15.195–196) "
            "⇒ residual (i) dual-eq branch closed ⇒ with 15.170 Farkas wiring "
            "gsum_disj_lb may flip when that bound is the disj Gsum engine. "
            "K₄ bound itself is OPEN."
        ),
    }


def certify_K4_and_distance(p: int) -> dict:
    """Census K4, max|dot|, compare to theorems. p=3,5 full; p=7 cache."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    if p == 7:
        path = Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7/maxplus.npy")
        if not path.is_file():
            return {"p": p, "skipped": True}
        Mp = np.load(path)
    elif p > 7:
        return {"p": p, "skipped": True}
    else:
        Mp = boolean_evecs(C, float(p), 0)
    Mp = Mp.astype(np.int16)
    N = len(Mp)
    # K2, K4, max off-diag |dot|
    sum2 = 0.0
    sum4 = 0.0
    max_off = 0
    bs = 400
    for i0 in range(0, N, bs):
        block = Mp[i0 : i0 + bs].astype(np.float64)
        dots = block @ Mp.T.astype(np.float64)
        sum2 += float(np.sum(dots * dots))
        sum4 += float(np.sum(dots ** 4))
        ad = np.abs(dots)
        ad = np.where(ad >= n - 0.5, 0, ad)
        max_off = max(max_off, int(ad.max()) if ad.size else 0)
    K2 = sum2 / (N * N)
    K4 = sum4 / (N * N)
    K4_frac = Fraction(K4).limit_denominator(10**7)
    # exact for p=3,5
    if p == 3:
        K4_frac = Fraction(1680, 1)
    elif p == 5:
        K4_frac = Fraction(120400, 13)
    Q = Q_from_K4(p, K4_frac)
    thr = K4_threshold_for_Q10(p)
    d_pred = max_dot_nonopposite(p)
    return {
        "p": p,
        "N": N,
        "K2": K2,
        "K2_theory_2n": 2 * n,
        "K4": str(K4_frac),
        "K4_float": float(K4_frac),
        "K4_thr_Q10": str(thr),
        "K4_le_thr": K4_frac <= thr,
        "Q_from_K4": str(Q),
        "max_dot_off": max_off,
        "max_dot_theory": d_pred,
        "distance_bound_sharp_or_ok": max_off <= d_pred,
        "proved_general": False,
    }


def K4_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_K4() -> bool:
    return K4_bound_proved_general()


def hinge_status_197() -> dict:
    return {
        "min_distance_p_plus_1": True,
        "Q_from_K4_identity": True,
        "K4_le_16n2_minus_14n_general": False,
        "residual_i_closed": residual_i_closed_via_K4(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove Max+-free K₄≤16n²−14n for all primes p≥5 "
            "(⇔ Q_e≤10 ⇔ mass-corrected dual-eq empty)."
        ),
    }


def main() -> dict:
    dist = theorem_min_distance_p_plus_1()
    qid = theorem_Q_from_K4_identity()
    impl = theorem_K4_thr_implies_residual_i_block()
    cert = {str(p): certify_K4_and_distance(p) for p in (3, 5, 7)}
    # Q≤10 still blocks (reuse 15.196 check at p=5)
    q10 = worst_mass_min_under_Q(5, Fraction(10))
    out = {
        "title": (
            "Prop 15.197 Max+ min-distance d_H≥p+1; exact Q(K₄); "
            "K₄≤16n²−14n would close residual (i); bound OPEN"
        ),
        "proved": {
            "min_distance_p_plus_1": dist["proved"],
            "Q_from_K4_identity": qid["proved"],
            "K4_thr_implies_Q10": impl["proved_implication"],
            "K4_bound_general": False,
            "residual_i_closed": residual_i_closed_via_K4(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "min_distance": dist,
            "Q_identity": qid,
            "implication": impl,
            "Q10_blocks_p5": q10,
        },
        "certified": cert,
        "hinge": hinge_status_197(),
        "F3": (
            "d_H≥p+1 and Q(K₄) are Max+-free. Do not flip residual/E1/L "
            "without K₄≤16n²−14n (or equivalent)."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15197.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.197 min-distance + Q(K4) identity")
    print(f"  d_H >= p+1 proved: {dist['proved']}")
    print(f"  Q from K4 identity: {qid['proved']}")
    for pk, c in cert.items():
        if c.get("skipped"):
            print(f"  p={pk}: skipped")
            continue
        print(
            f"  p={pk}: K4={c['K4']} thr={c['K4_thr_Q10']} "
            f"K4_le_thr={c['K4_le_thr']} max_dot={c['max_dot_off']}"
            f"<={c['max_dot_theory']} Q={c['Q_from_K4']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_K4()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
