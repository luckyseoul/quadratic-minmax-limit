#!/usr/bin/env python3
"""
Prop 15.198 — Residual-(i) K₄ reductions: Wick_hi ≤ thr; exact K₄ partition;
census exact K₄ at p=7; hinge still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.195–197)
  K₄ := E_{y,z∈Max+}[(y·z)⁴],   n = p²+1,   S := E[yyᵀ] = I + C/p (15.189).
  Q_e = K₄/(n(n−1)) − 2(3n−2)/(n−1).
  Residual-(i) dual-eq blocked when Q_e ≤ 10 ⇔ K₄ ≤ 16n² − 14n (15.196–197).

PROVED Max+-free
  A. Wick envelope vs residual-(i) threshold.
        Wick_hi := 12n² + 48n.
     For every prime p≥5 (n=p²+1≥26):
        Wick_hi ≤ 16n² − 14n
     ⇔ 0 ≤ 4n² − 62n = 2n(2n − 31) ⇔ n ≥ 16, true on n≥26.  ∎
     At p=3 (n=10): Wick_hi = 1680 > thr = 1460 (dual-eq feasible — consistent).

  B. Wick_hi forces Q below 10.  Substituting K₄ ≤ Wick_hi into the Q identity:
        Q ≤ Q_Wick := (12n²+48n)/(n(n−1)) − 2(3n−2)/(n−1)
                    = (6n + 52)/(n−1) = 6 + 58/(n−1).
     For n≥26: Q_Wick ≤ 6 + 58/25 = 8.32 < 10.  ∎
     So K₄ ≤ Wick_hi ⇒ Q ≤ Q_Wick < 10 ⇒ mass-corrected dual-eq empty
     (15.195–196) for all primes p≥5 (and the p∈[5,47] Q≤10 check).

  C. Exact K₄ via boolean / Wick residual (any measure with E[yyᵀ]=S=I+C/p
     and y∈{±1}ⁿ, in particular Max+ with 15.189).
     Write m₄(ijkl)=E[y_i y_j y_k y_l], W_ijkl = S_ij S_kl+S_ik S_jl+S_il S_jk,
     η = m₄ − W.  Then (using yᵀ S y = 2n constant on Max+):
        ∑_{ijkl} W² = Wick_hi = 12n²+48n,
        ∑ m₄ W = 12 n²,
        K₄ = ∑ m₄² = 12n² − 48n + ‖η‖_F².
     Collision (repeated-index) contribution to ‖η‖_F² is Max+-free:
        C_diag = 4n(11n−14)/p² = 20n + 12n(2n−3)/(n−1)
     (same formula as 15.95; derived from S²=2S, E[y_i⁴]=1, E[y_i² y_j²]=1).
     The 4-distinct part is 24 ∑_{4-sets} η_S² with η_S = m₄(S) − T(S)/p²,
     T = sum of three matching products of C.  Hence
        K₄ = 12n² − 48n + C_diag + 24 ∑_S η_S²
            = 12n² − 4n − 12n/(n−1) + 24 ∑_S η_S².  ∎
     Equivalently K₄ ≤ Wick_hi ⇔ ∑_S η_S² ≤ n(13n−10)/(6(n−1)).

  D. Exact K₄ ↔ e₄ mass on Max+ (15.197 path).  With
        e₄(w)=(s⁴−(6n−8)s²+3n²−6n)/24,  s=∑ w_i,
     one has ∑_S m₊(S)² = E e₄(y⊙z) and
        K₄ = 24 ∑_S m₊(S)² + 9n² − 10n.  ∎

  E. Implication chain for residual (i).  IF for all primes p≥5 either
        (i)  K₄ ≤ 16n²−14n, or
        (ii) K₄ ≤ Wick_hi (⇒ (i) by A), or
        (iii) ∑_S η_S² ≤ n(13n−10)/(6(n−1)) (⇔ (ii) by C), or
        (iv) |μ₄|≤2/n on |κ|=1 (15.188), or ker-box empty,
     THEN dual-eq residual-(i) closes and gsum_disj_lb_proved_general may flip
     when the bound is wired as the disj Gsum engine.  ∎

CERTIFIED (not general)
  F. Exact K₄ fractions:
        p=3: 1680 = Wick_hi  (equality; thr fails)
        p=5: 120400/13 ≈ 9261.54 < Wick_hi=9360 < thr=10452
        p=7: 5218435600/167281 ≈ 31195.63 < Wick_hi=32400 < thr=39300
     Ratio K₄/Wick_hi: 1, 0.9895, 0.9628 (decreasing).
  G. ∑_S η_S² vs Wick budget n(13n−10)/(6(n−1)):
        p=3: equality 200/9;
        p=5: 17144/325 ≈ 52.75 < 56.85;
        p=7: strict.

OPEN (blocks residual i / predicate flip)
  H. Prove Max+-free K₄ ≤ Wick_hi (or ≤ 16n²−14n, or ∑ η_S² budget, or
     |μ|≤2/n / ker-box) for all primes p≥5.
     Note: K₄ ≤ Wick_hi is the same scalar as ∑M² ≤ Wick_hi in 15.95
     (Path C bi-tight companion); residual (i) thr is strictly weaker
     (room Wick_hi…16n²−14n) but still open.  Do **not** flip predicates.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15198.json
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
from e1_gmin_m4_prop15197 import (  # noqa: E402
    K4_threshold_for_Q10,
    Q_from_K4,
)


def wick_hi(p: int) -> int:
    n = n_of(p)
    return 12 * n * n + 48 * n


def wick_lo(p: int) -> int:
    n = n_of(p)
    return 12 * n * n - 48 * n


def C_diag(p: int) -> Fraction:
    """Repeated-index η-mass: 4n(11n−14)/p²."""
    n = n_of(p)
    return Fraction(4 * n * (11 * n - 14), p * p)


def eta_sumsq_budget_for_wick_hi(p: int) -> Fraction:
    """∑_S η_S² ≤ n(13n−10)/(6(n−1)) ⇔ K₄ ≤ Wick_hi."""
    n = n_of(p)
    return Fraction(n * (13 * n - 10), 6 * (n - 1))


def Q_from_wick_hi(p: int) -> Fraction:
    """Q when K₄ = Wick_hi: 6 + 58/(n−1)."""
    n = n_of(p)
    return Fraction(6 * n + 52, n - 1)


def theorem_wick_hi_le_K4_thr() -> dict:
    """Max+-free: Wick_hi ≤ 16n²−14n for all primes p≥5."""
    primes = [p for p in range(3, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        w = wick_hi(p)
        thr = int(K4_threshold_for_Q10(p))
        holds = w <= thr
        # algebra: thr - w = 4n² - 62n = 2n(2n-31)
        diff = 2 * n * (2 * n - 31)
        rows[str(p)] = {
            "n": n,
            "Wick_hi": w,
            "K4_thr": thr,
            "thr_minus_Wick": diff,
            "Wick_le_thr": holds,
            "expected": p >= 5,
            "match": holds == (p >= 5),
        }
        if holds != (p >= 5):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Wick_hi=12n²+48n ≤ 16n²−14n ⇔ 2n(2n−31)≥0 ⇔ n≥16. "
            "For primes p≥5, n=p²+1≥26, so Wick_hi ≤ K₄-threshold. "
            "At p=3, n=10, Wick_hi > thr (correct)."
        ),
        "by_p_sample": {k: rows[k] for k in ("3", "5", "7", "11") if k in rows},
    }


def theorem_Q_wick_lt_10() -> dict:
    """K₄≤Wick_hi ⇒ Q ≤ 6+58/(n−1) < 10 for p≥5."""
    primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        qw = Q_from_wick_hi(p)
        n = n_of(p)
        # verify equals formula from K4=Wick_hi
        q_alt = Q_from_K4(p, wick_hi(p))
        rows[str(p)] = {
            "Q_Wick": str(qw),
            "Q_from_K4_Wick": str(q_alt),
            "match": qw == q_alt,
            "Q_Wick_lt_10": qw < 10,
            "formula": "6+58/(n-1)",
            "check_formula": qw == 6 + Fraction(58, n - 1),
        }
        if qw != q_alt or qw >= 10 or qw != 6 + Fraction(58, n - 1):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Under K₄≤Wick_hi: Q≤(6n+52)/(n−1)=6+58/(n−1)<10 for all n≥26 "
            "(p≥5). At p=5, Q_Wick=8.32; at p=3, Q_Wick=112/9>10."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11") if k in rows},
    }


def theorem_K4_partition_identity() -> dict:
    return {
        "proved": True,
        "theorem": (
            "With S=I+C/p and boolean Max+: K₄=12n²−48n+C_diag+24∑_S η_S², "
            "C_diag=4n(11n−14)/p². Equivalent: K₄≤Wick_hi ⇔ "
            "∑_S η_S² ≤ n(13n−10)/(6(n−1)). Also K₄=24∑ m₊²+9n²−10n."
        ),
        "C_diag_formula": "4n(11n-14)/p^2",
        "eta_budget_wick": "n(13n-10)/(6(n-1))",
    }


def theorem_reduction_chain() -> dict:
    return {
        "proved_implications": True,
        "K4_bound_proved_general": False,
        "theorem": (
            "K₄≤Wick_hi ⇒ Q≤Q_Wick<10 ⇒ dual-eq empty (15.195–196) for p≥5. "
            "Wick_hi≤16n²−14n for p≥5 (A). Partition (C) equates K₄≤Wick_hi "
            "with an η_S² budget. OPEN: prove any of these for general p≥5."
        ),
    }


def certify_K4_exact(p: int) -> dict:
    """Exact K₄ from known fractions / p=7 integer sum."""
    n = n_of(p)
    thr = K4_threshold_for_Q10(p)
    w = wick_hi(p)
    if p == 3:
        K4 = Fraction(1680)
    elif p == 5:
        K4 = Fraction(120400, 13)
    elif p == 7:
        # 5218435600/167281 from full Max+ census (sum of fourth powers / M²)
        K4 = Fraction(5218435600, 167281)
    else:
        return {"p": p, "skipped": True}
    Q = Q_from_K4(p, K4)
    coll = C_diag(p)
    eta_tot = K4 - wick_lo(p)
    eta_dist = eta_tot - coll  # 24 ∑ η_S²
    sum_eta = eta_dist / 24
    bud = eta_sumsq_budget_for_wick_hi(p)
    return {
        "p": p,
        "n": n,
        "K4": str(K4),
        "K4_float": float(K4),
        "Wick_hi": w,
        "K4_thr": str(thr),
        "K4_le_Wick": K4 <= w,
        "K4_le_thr": K4 <= thr,
        "K4_over_Wick": float(K4 / w),
        "Q": str(Q),
        "Q_Wick": str(Q_from_wick_hi(p)),
        "C_diag": str(coll),
        "sum_eta_S2": str(sum_eta),
        "eta_budget_wick": str(bud),
        "eta_le_budget": sum_eta <= bud,
        "proved_general": False,
    }


def K4_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_K4_wick() -> bool:
    return K4_bound_proved_general()


def hinge_status_198() -> dict:
    return {
        "wick_hi_le_thr_p_ge_5": True,
        "Q_wick_lt_10_p_ge_5": True,
        "K4_partition_identity": True,
        "K4_le_Wick_hi_general": False,
        "K4_le_16n2_minus_14n_general": False,
        "residual_i_closed": residual_i_closed_via_K4_wick(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove Max+-free K₄≤Wick_hi=12n²+48n (or ≤16n²−14n / η-budget / "
            "|μ|≤2/n) for all primes p≥5. Then residual (i) dual-eq closes."
        ),
    }


def main() -> dict:
    a = theorem_wick_hi_le_K4_thr()
    b = theorem_Q_wick_lt_10()
    c = theorem_K4_partition_identity()
    d = theorem_reduction_chain()
    cert = {str(p): certify_K4_exact(p) for p in (3, 5, 7)}
    out = {
        "title": (
            "Prop 15.198 Wick_hi≤K₄-thr; Q_Wick<10; K₄ partition; "
            "K₄≤Wick_hi OPEN (residual i)"
        ),
        "proved": {
            "wick_hi_le_thr_p_ge_5": a["proved"],
            "Q_wick_lt_10": b["proved"],
            "K4_partition_identity": c["proved"],
            "reduction_chain": d["proved_implications"],
            "K4_bound_general": False,
            "residual_i_closed": residual_i_closed_via_K4_wick(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "wick_le_thr": a,
            "Q_wick": b,
            "partition": c,
            "reduction": d,
        },
        "certified": cert,
        "hinge": hinge_status_198(),
        "F3": (
            "Wick_hi≤thr and Q_Wick<10 are Max+-free. Do not flip residual/E1/L "
            "without K₄≤Wick_hi (or equivalent)."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15198.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.198 Wick_hi residual-(i) reductions")
    print(f"  Wick_hi <= thr for p>=5: {a['proved']}")
    print(f"  Q_Wick < 10 for p>=5: {b['proved']}")
    print(f"  K4 partition identity: {c['proved']}")
    for pk, crow in cert.items():
        if crow.get("skipped"):
            print(f"  p={pk}: skipped")
            continue
        print(
            f"  p={pk}: K4={crow['K4']} Wick={crow['Wick_hi']} "
            f"K4/Wick={crow['K4_over_Wick']:.4f} "
            f"K4_le_thr={crow['K4_le_thr']} eta_ok={crow['eta_le_budget']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_K4_wick()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
