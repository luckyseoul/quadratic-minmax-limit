#!/usr/bin/env python3
"""
Prop 15.190 — Scheme-ker of Gsum; dual-eq ker-box obstruction analysis;
3-point vanishing (certified); residual (i) still OPEN for general p.

Does **not** flip gsum_disj_lb_proved_general without a full Max+-free
ker-box emptiness proof for all p≥5. Soft-close forbidden.

PROVED Max+-free
  A. Scheme-ker ⊆ ker(F+) ∩ ker(F−) ∩ ker(Gsum).
     For f∈R^n with ∑_i f_i = 0 set κ_{ij} := f_i + f_j (i<j).
     Then for any y with Cy = ±p y:
        ∑_{i<j} κ_{ij} C_{ij} y_i y_j
          = (1/2) ∑_{i≠j} (f_i+f_j) C_{ij} y_i y_j
          = (1/2) ∑_i f_i y_i ∑_{j≠i} C_{ij} y_j
            + (1/2) ∑_j f_j y_j ∑_{i≠j} C_{ij} y_i
          = (1/2) ∑_i f_i y_i (p y_i) + (1/2) ∑_j f_j y_j (p y_j)
          = p ∑_i f_i = 0.  ∎
     Hence F±κ=0 and Gsum κ=0. Also 1ᵀκ = (n−1)∑f = 0.  ∎

  B. Scheme-ker dual-eq bound.
     Maximise f_a+f_b subject to ∑f=0 and −α ≤ f_i+f_j ≤ 1−α for all
     {i,j}≠{a,b}, with α=6/(pn).
     Two-level ansatz f_a=f_b=c/2, f_other=−c/(n−2) is optimal and gives
        max(f_a+f_b) = α(n−2)/2 = 3(n−2)/(p n).
     Proof of optimality: 2 min f ≥ −α ⇒ min f ≥ −α/2; with ∑f=0 the
     max of f_a+f_b is at most α(n−2)/2 by standard equalising.
     Check: 3(n−2)/(pn) < 2−α for all primes p≥3, since
        3(n−2) < 2pn − 6  ⇔  3n < 2pn  ⇔  3 < 2p.  ∎
     Thus pure scheme-ker never reaches dual-eq need κ_e=2−α.

  C. Adjacent Gsum=0 and π_ij=C_ij/p already in 15.189.

CERTIFIED (not general)
  D. 3-point moments E_±[y_i y_j y_k]=0 for all distinct triples (p=3,5 full;
     sample p=7). Consistent with (but stronger than) a pure 2-design.
  E. Dual-eq ker-box (Gsum κ=0, 1ᵀκ=0, κ_e=2−α, box on g≠e):
        p=3: FEASIBLE (max κ_e≈2.80 ≥ 1.80)
        p=5: INFEASIBLE (max κ_e≈0.811 < 1.954)
        p=7: INFEASIBLE (max κ_e≈0.593 < 1.983)
     So residual-(i) dual-eq is blocked at p=5,7 by ker-box; at p=3 it is
     not (extra ker beyond scheme-ker lifts max κ_e past the need).
  F. Boolean 4-cube moment LP with π and 3-pt=0 recovers only |μ|≤1−2/p
     (same as G+ PSD); not enough for residual (i).

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free that for all primes p≥5
        max{ κ_e : Gsum κ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e) } < 2−α.
     Equivalent: dual-eq ker-box empty for all p≥5.
     Partial: scheme-ker alone always blocked (B); full ker blocked at p=5,7
     by census (E). Need general bound on the extra-ker contribution.
  H. Or prove |μ₄|≤1/(2p) (or ≤2/n) on |κ|=1 for all p≥5.

Until G or H: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15190.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
)
from e1_gmin_m4_prop15189 import theorem_pairwise_pi_equals_C_over_p  # noqa: E402


def alpha_particular(p: int) -> Fraction:
    return Fraction(6, p * n_of(p))


def scheme_ker_max_kappa_e(p: int) -> Fraction:
    """α(n−2)/2 = 3(n−2)/(p n)."""
    n = n_of(p)
    return alpha_particular(p) * (n - 2) / 2


def theorem_scheme_ker_in_gsum_ker() -> dict:
    return {
        "proved": True,
        "theorem": (
            "κ_{ij}=f_i+f_j with ∑f=0 lies in ker(F+)∩ker(F−)∩ker(Gsum) "
            "for any conference boolean evecs (algebra of Cy=±py)."
        ),
    }


def theorem_scheme_ker_dual_eq_blocked() -> dict:
    primes = [p for p in range(3, 40) if p == 2 or all(p % d for d in range(2, int(p**0.5) + 1))]
    primes = [p for p in primes if p >= 3 and all(p % d for d in range(2, int(p**0.5) + 1))]
    # simple prime list
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    rows = {}
    ok = True
    for p in primes:
        mx = scheme_ker_max_kappa_e(p)
        need = 2 - alpha_particular(p)
        rows[str(p)] = {
            "scheme_max_kappa_e": str(mx),
            "need_2_minus_alpha": str(need),
            "blocked": mx < need,
        }
        if not (mx < need):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Scheme-ker max κ_e = α(n−2)/2 = 3(n−2)/(pn) < 2−α for all "
            "primes p≥3 (⇔ 3<2p)."
        ),
        "by_p": rows,
    }


def certify_three_point_vanish(p: int = 5, samples: int = 2000) -> dict:
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    rng = np.random.default_rng(0)
    max_abs = 0.0
    for _ in range(samples):
        i, j, k = rng.choice(n, 3, replace=False)
        max_abs = max(
            max_abs,
            abs(float(np.mean(Mp[:, i] * Mp[:, j] * Mp[:, k]))),
            abs(float(np.mean(Mm[:, i] * Mm[:, j] * Mm[:, k]))),
        )
    # full for p<=5
    if p <= 5:
        max_abs = 0.0
        for i, j, k in combinations(range(n), 3):
            max_abs = max(
                max_abs,
                abs(float(np.mean(Mp[:, i] * Mp[:, j] * Mp[:, k]))),
                abs(float(np.mean(Mm[:, i] * Mm[:, j] * Mm[:, k]))),
            )
    return {
        "p": p,
        "max_abs_3pt": max_abs,
        "vanishes": max_abs < 1e-12,
        "proved_general": False,
    }


def certify_kerbox(p: int = 5) -> dict:
    """Census: dual-eq ker-box feasibility and max κ_e."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power
    from scipy.optimize import linprog

    if p > 5:
        # p=7 uses cache if present
        from pathlib import Path

        cache = Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7")
        mp_path, mm_path = cache / "maxplus.npy", cache / "maxminus.npy"
        if not (mp_path.is_file() and mm_path.is_file()):
            return {"p": p, "skipped": True, "reason": "no Max± cache"}
        Mp = np.load(mp_path).astype(float)
        Mm = np.load(mm_path).astype(float)
        C = paley_conference_prime_power(p).astype(float)
    else:
        C = paley_conference_prime_power(p).astype(float)
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        Mm = boolean_evecs(C, -float(p), 0).astype(float)
    n = C.shape[0]
    edges = list(combinations(range(n), 2))
    E = len(edges)
    alpha = float(alpha_particular(p))
    ke = 2 - alpha
    Fp = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mp])
    Fm = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mm])
    Gsum = (Fp.T @ Fp) / len(Mp) + (Fm.T @ Fm) / len(Mm)
    A = np.vstack([Gsum, np.ones((1, E))])
    b = np.zeros(E + 1)
    # max κ_e
    bounds = [(-alpha, 1 - alpha) if i else (None, None) for i in range(E)]
    c = np.zeros(E)
    c[0] = -1
    res = linprog(c, A_eq=A, b_eq=b, bounds=bounds, method="highs")
    max_ke = float(-res.fun) if res.success else None
    # feasibility at need
    bounds_fix = [(-alpha, 1 - alpha) for _ in range(E)]
    bounds_fix[0] = (ke, ke)
    resf = linprog(np.zeros(E), A_eq=A, b_eq=b, bounds=bounds_fix, method="highs")
    return {
        "p": p,
        "alpha": alpha,
        "need_kappa_e": ke,
        "max_kappa_e": max_ke,
        "scheme_max": float(scheme_ker_max_kappa_e(p)),
        "dual_eq_feasible": bool(resf.success),
        "blocked": max_ke is not None and max_ke < ke - 1e-9,
    }


def hinge_status_190() -> dict:
    return {
        "residual_ii_closed": False,  # full open (15.193); affine branch closed by 15.179
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "scheme_ker_in_gsum_ker": True,
        "scheme_ker_blocks_dual_eq": True,
        "full_ker_blocks_dual_eq_general": False,
        "open": (
            "Prove full dual-eq ker-box empty for all p≥5 (max κ_e<2−α), "
            "or |μ|≤1/(2p) on |κ|=1."
        ),
    }


def main() -> dict:
    sk = theorem_scheme_ker_in_gsum_ker()
    sk_block = theorem_scheme_ker_dual_eq_blocked()
    pi = theorem_pairwise_pi_equals_C_over_p()
    cert3 = {str(p): certify_three_point_vanish(p) for p in (3, 5)}
    cert_box = {str(p): certify_kerbox(p) for p in (3, 5)}
    # p=7 if cache available
    c7 = certify_kerbox(7)
    if not c7.get("skipped"):
        cert_box["7"] = c7
    out = {
        "title": (
            "Prop 15.190 scheme-ker ⊆ Gsum-ker; dual-eq ker-box analysis; "
            "residual (i) OPEN for general p"
        ),
        "proved": {
            "scheme_ker_in_gsum_ker": sk["proved"],
            "scheme_ker_dual_eq_blocked": sk_block["proved"],
            "pairwise_pi": pi["proved"],
            "full_ker_box_empty_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "scheme_ker": sk,
            "scheme_ker_bound": sk_block,
            "pairwise_pi": pi,
        },
        "certified": {
            "three_point_vanish": cert3,
            "kerbox": cert_box,
        },
        "hinge": hinge_status_190(),
        "F3": (
            "Scheme-ker always blocks dual-eq need, but extra ker can restore "
            "feasibility (p=3). Need Max+-free proof that extra ker cannot "
            "reach κ_e=2−α under box for all p≥5. Soft-close forbidden."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15190.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.190 scheme-ker; ker-box analysis; residual (i) OPEN")
    print(f"  scheme-ker in Gsum-ker proved={sk['proved']}")
    print(f"  scheme-ker blocks dual-eq proved={sk_block['proved']}")
    print(f"  scheme max κ_e p5={scheme_ker_max_kappa_e(5)} need={2-alpha_particular(5)}")
    for k, v in cert_box.items():
        if v.get("skipped"):
            continue
        print(
            f"  kerbox p={k}: max_ke={v.get('max_kappa_e')} "
            f"feasible={v.get('dual_eq_feasible')} blocked={v.get('blocked')}"
        )
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
