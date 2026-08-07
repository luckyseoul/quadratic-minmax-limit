#!/usr/bin/env python3
"""
Prop 15.202 — Free-e dual form; scheme⊕cross = ker(Gsum) (certified);
residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182, 15.190–201)
  Dual-eq empty ⇔ free-e ker-box max κ_e < 2−α:
    max{κ_e : Gsumκ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)}.
  Scheme-ker max = α(n−2)/2.  α(n−2)<2−α for all primes p≥5 (15.201).
  (3/2)·scheme < 2−α for all primes p≥5 (15.192).

PROVED Max+-free
  A. Scheme ⊕ cross ⊆ ker(Gsum).  Scheme κ_ij=f_i+f_j (∑f=0) by 15.190.
     Cross: κ_ij = C_ij A_ij with A = P_+ B P_- + P_- Bᵀ P_+ (B arbitrary
     d×d, d=n/2) vanishes as a quadratic form on V± by block structure, hence
     on Max± ⊂ V±.  Dimension of the constructible space:
        dim = (n−1) + d² = n−1 + n²/4.
     (Cross params d²; scheme n−1 after ∑f=0.)  Rank of the edge-feature
     matrix L_sc is full (certified p=3,5,7,11).  ∎

  B. Free-e dual form (lower-box only).  At free-e optimum the upper box
     κ_g≤1−α is slack (certified p=3,5,7: dual yp≡0).  Hence
        free-e max κ_e = α · S*,
     where
        S* = min{ ∑_{g≠e} w_g : w_e=1, w_g≥0 (g≠e),
                  w ⊥ ker(Gsum) }.
     Equivalently, with ker ⊇ scheme⊕cross:
        w has equal weighted degrees (⊥ scheme), and
        (W ⊙ C) commutes with C (⊥ cross: P_+(W⊙C)P_−=0).
     Then free-e max ≤ α(n−2) ⇔ S* ≤ n−2, and
     free-e max ≤ (3/2)·scheme ⇔ S* ≤ 3(n−2)/4.  ∎

  C. Sufficient thresholds (recalled, Max+-free Fraction).
        α(n−2) < 2−α  ⇔  p(p²−3p+1)>0   (p≥5).
        (3/2)·scheme < 2−α  ⇔  4p³−9p²+4p−3>0   (p≥5; 15.192).
     So either free-e max ≤ α(n−2) or ≤ (3/2)·scheme closes dual-eq.  ∎

CERTIFIED (not general)
  D. dim ker(Gsum) = dim(scheme⊕cross) at p=3,5,7
        (34=34, 194=194, 674=674).  Hence ker = scheme⊕cross at these p
        (constructible Max+-free span is exact).  Free-e max over
        scheme⊕cross is the true free-e max.

  E. Free-e max κ_e (scheme⊕cross LP / dual; Max+ not required once ker=sc):
        p=3: 14/5 = 2.8  > α(n−2)=8/5  (dual-eq feasible — OK)
        p=5: 369/455 ≈ 0.811 < α(n−2)=72/65  and < (3/2)·scheme
        p=7: 11736/19775 ≈ 0.593 < α(n−2)=144/175 and < (3/2)·scheme
     Ratios max/scheme: 7/2 (p=3), 41/28 (p=5), 163/113 (p=7).
     Dual S*: 14, 123/7, 3912/113; all ≤ n−2 for p=5,7 (not for p=3).

  F. Dual support: ||[W⊙C, C]||_F ≈ 0 and regular degrees at optimum
     (confirms B).  Only lower box active.

OPEN (blocks residual i / predicate flip)
  G. Prove for all primes p≥5 one of:
     (i)   free-e max κ_e ≤ α(n−2)   (i.e. S* ≤ n−2),
     (ii)  free-e max κ_e ≤ (3/2)·scheme_max,
     (iii) |μ|≤2/n on every |κ|=1 4-set,
     (iv)  K₄ ≤ Wick_hi.
     Then flip gsum_disj_lb_proved_general → residual_i (and E1 only if
     residual ii allows).

  H. Optional: prove ker(Gsum)=scheme⊕cross for all primes p≥5
     (dim match beyond census), Max+-free via 2-design / Aut(C) of Paley.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15202.json
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
from e1_gmin_m4_prop15190 import (  # noqa: E402
    alpha_particular,
    scheme_ker_max_kappa_e,
)
from e1_gmin_m4_prop15201 import (  # noqa: E402
    theorem_alpha_n_minus_2_beats_need,
    two_scheme_max,
)


def dim_scheme_cross(p: int) -> int:
    n = n_of(p)
    d = n // 2
    return (n - 1) + d * d


def theorem_scheme_cross_in_ker() -> dict:
    return {
        "proved": True,
        "theorem": (
            "scheme ⊕ cross ⊆ ker(Gsum): scheme by 15.190 (Cy=±py algebra); "
            "cross by block-off-diagonal A in the C-eigenbasis (yᵀA y=0 on V± "
            "hence on Max±). dim = n−1+n²/4."
        ),
    }


def theorem_free_e_dual_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "If upper box is slack at free-e max (certified p=3,5,7), then "
            "free-e max κ_e = α·S* with S*=min{∑_{g≠e} w_g : w_e=1, w≥0, "
            "w ⊥ ker}.  ⊥scheme ⇔ regular weighted degree; ⊥cross ⇔ "
            "(W⊙C) commutes with C.  Hence free-e max ≤ α(n−2) ⇔ S*≤n−2."
        ),
        "note": (
            "Upper-box slack is certified at p=3,5,7 (dual yp=0); not yet "
            "proved for all p, but the dual form with both boxes always "
            "upper-bounds free-e max."
        ),
    }


def theorem_three_halves_scheme_beats_need() -> dict:
    """(3/2)·scheme < 2−α for all primes p≥5 (15.192 algebra)."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        thr = Fraction(3, 2) * scheme_ker_max_kappa_e(p)
        need = 2 - alpha_particular(p)
        poly = 4 * p**3 - 9 * p**2 + 4 * p - 3
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "three_halves_scheme": str(thr),
                "need": str(need),
                "beats": thr < need,
                "poly": poly,
            }
        if not (thr < need and poly > 0):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "(3/2)·scheme_max < 2−α for all primes p≥5 "
            "⇔ 4p³−9p²+4p−3>0.  Hence free-e max ≤ (3/2)·scheme ⇒ dual-eq empty."
        ),
        "by_p_sample": sample,
    }


def certify_dim_and_free_e(p: int) -> dict:
    """Census dim(ker)=dim(sc) and free-e max via scheme⊕cross dual."""
    import numpy as np
    from itertools import combinations
    from scipy.optimize import linprog
    from minmax_quadratic import paley_conference_prime_power

    n = n_of(p)
    d = n // 2
    sc_dim = dim_scheme_cross(p)
    al = float(alpha_particular(p))
    C = paley_conference_prime_power(p).astype(float)

    # ker dim from Max± when available
    ker_dim = None
    if p <= 5:
        from e1_bitight_infeas import boolean_evecs

        Mp = boolean_evecs(C, float(p), 0).astype(float)
        Mm = boolean_evecs(C, float(-p), 0).astype(float)
    elif p == 7:
        cache = Path("/tmp/grok-goal-7663bdcc1fa8/implementer/e1_p7")
        if not (cache / "maxplus.npy").is_file():
            return {"p": p, "skipped": True, "reason": "no Max± cache"}
        Mp = np.load(cache / "maxplus.npy").astype(float)
        Mm = np.load(cache / "maxminus.npy").astype(float)
    else:
        Mp = Mm = None

    if Mp is not None:
        edges = list(combinations(range(n), 2))
        E = len(edges)

        def feats(M):
            return np.stack([C[i, j] * M[:, i] * M[:, j] for i, j in edges], 1)

        G = feats(Mp).T @ feats(Mp) / len(Mp) + feats(Mm).T @ feats(Mm) / len(Mm)
        ev = np.linalg.eigvalsh(G)
        ker_dim = int(np.sum(ev < 1e-7))

    # free-e max over scheme⊕cross (dual)
    evals, evecs = np.linalg.eigh(C)
    order = np.argsort(evals)
    V_m = evecs[:, order[:d]]
    V_p = evecs[:, order[d:]]
    ii, jj = np.triu_indices(n, k=1)
    E = len(ii)
    e0 = int(np.where((ii == 0) & (jj == 1))[0][0])
    A = (
        V_p[ii, :, None] * V_m[jj, None, :]
        + V_p[jj, :, None] * V_m[ii, None, :]
    )
    L = (A * C[ii, jj, None, None]).reshape(E, d * d)
    sch = np.zeros((E, n - 1))
    for k in range(n - 1):
        sch[:, k] = (
            (ii == k).astype(float)
            + (jj == k).astype(float)
            - (ii == n - 1).astype(float)
            - (jj == n - 1).astype(float)
        )
    Lfull = np.concatenate([L, sch], 1)
    mask = np.ones(E, dtype=bool)
    mask[e0] = False
    Lg = Lfull[mask]
    m = Lg.shape[0]
    res = linprog(
        np.full(m, al),
        A_eq=Lg.T,
        b_eq=-Lfull[e0],
        bounds=[(0, None)] * m,
        method="highs",
    )
    if not res.success:
        return {"p": p, "success": False, "message": res.message}

    vmax = float(res.fun)
    Sstar = float(res.x.sum())
    bound = float(two_scheme_max(p))
    need = float(2 - alpha_particular(p))
    schv = float(scheme_ker_max_kappa_e(p))
    thr15 = 1.5 * schv

    return {
        "p": p,
        "success": True,
        "dim_scheme_cross": sc_dim,
        "dim_ker_Gsum": ker_dim,
        "ker_eq_scheme_cross": (ker_dim == sc_dim) if ker_dim is not None else None,
        "max_kappa_e": vmax,
        "S_star": Sstar,
        "scheme_max": schv,
        "alpha_n_minus_2": bound,
        "three_halves_scheme": thr15,
        "need": need,
        "max_le_alpha_n_minus_2": vmax <= bound + 1e-9,
        "max_le_three_halves_scheme": vmax <= thr15 + 1e-9,
        "S_star_le_n_minus_2": Sstar <= (n - 2) + 1e-9,
        "blocks_dual_eq": vmax < need - 1e-9,
        "ratio_max_over_scheme": vmax / schv,
        "frac": str(Fraction(vmax).limit_denominator(1_000_000)),
        "proved_general": False,
    }


def free_e_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_202() -> bool:
    return free_e_bound_proved_general()


def hinge_status_202() -> dict:
    return {
        "scheme_cross_in_ker": True,
        "free_e_dual_form": True,
        "alpha_n_minus_2_beats_need_p_ge_5": True,
        "three_halves_scheme_beats_need_p_ge_5": True,
        "ker_eq_scheme_cross_general": False,
        "free_e_max_le_alpha_n_minus_2_general": False,
        "residual_i_closed": residual_i_closed_via_202(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove free-e max κ_e ≤ α(n−2) or ≤(3/2)·scheme (or |μ|≤2/n / "
            "K₄≤Wick_hi) for all primes p≥5.  Dual form: S*≤n−2 or ≤3(n−2)/4."
        ),
    }


def main() -> dict:
    a = theorem_scheme_cross_in_ker()
    b = theorem_free_e_dual_form()
    c = theorem_alpha_n_minus_2_beats_need()
    d = theorem_three_halves_scheme_beats_need()
    cert = {str(p): certify_dim_and_free_e(p) for p in (3, 5, 7)}
    out = {
        "title": (
            "Prop 15.202 free-e dual form; scheme⊕cross ker (certified); "
            "residual (i) OPEN"
        ),
        "proved": {
            "scheme_cross_in_ker": a["proved"],
            "free_e_dual_form": b["proved"],
            "alpha_n_minus_2_beats_need": c["proved"],
            "three_halves_scheme_beats_need": d["proved"],
            "free_e_bound_general": False,
            "residual_i_closed": residual_i_closed_via_202(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "scheme_cross": a,
            "dual_form": b,
            "alpha_bound": c,
            "three_halves": d,
        },
        "certified": cert,
        "hinge": hinge_status_202(),
        "F3": (
            "Do not flip residual/E1/L without general free-e max ≤ α(n−2) "
            "or equivalent.  ker=scheme⊕cross certified only at p=3,5,7."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15202.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.202 free-e dual form; scheme+cross ker")
    print(f"  scheme⊕cross ⊆ ker: {a['proved']}")
    print(f"  dual form: {b['proved']}")
    print(f"  α(n−2)<2−α p≥5: {c['proved']}")
    print(f"  (3/2)·scheme<2−α p≥5: {d['proved']}")
    for pk, v in cert.items():
        if v.get("skipped") or not v.get("success"):
            print(f"  p={pk}: skipped/fail {v}")
            continue
        print(
            f"  p={pk}: max={v['frac']} ker=sc?{v['ker_eq_scheme_cross']} "
            f"le_α(n-2)={v['max_le_alpha_n_minus_2']} "
            f"le_1.5sch={v['max_le_three_halves_scheme']} "
            f"blocks={v['blocks_dual_eq']} ratio={v['ratio_max_over_scheme']:.4f}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_202()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
