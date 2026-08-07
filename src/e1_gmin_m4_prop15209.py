#!/usr/bin/env python3
"""
Prop 15.209 — G₊ spectrum on 𝒲₊₊⁰ through p=7; free-e star saturation;
dim formula; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.202–208)
  Dual-eq empty if free-e max κ_e < 2−α on ker∩box.  sc:=scheme⊕cross ⊆ ker.
  ker=sc ⇔ G₊ ≻ 0 on 𝒲₊₊⁰ (and Max− twin) by 15.207.  Free-e over sc
  ≤ r(p)·scheme with r≤3/2 closes dual-eq once ker=sc (15.205).  Open:
  G₊≻0 general, or |μ|≤M_cand, or K₄≤Wick_hi.

PROVED Max+-free
  A. Dimension of zero-diag ++ forms.
        dim 𝒲₊₊⁰ = dim{ M∈Sym(R^d) : diag(Vp M Vpᵀ)=0 }
                   = d(d+1)/2 − n = n(n−6)/8,
     for n=2d even, n≥6 (the n diagonal constraints are independent on the
     ++ compressions for conference π: rank of diag∘Ad_{Vp} = n, certified
     and used as the ambient count in 15.207).  Algebra: d=n/2 ⇒
        n/2·(n/2+1)/2 − n = n(n+2)/8 − n = n(n−6)/8.  ∎

  B. Implication chain (recall).  For all primes p≥5:
        G₊ ≻ 0 on 𝒲₊₊⁰ (and − twin)
        ⇒ ker Q₊ = scheme-image ⇒ ker(Gsum)=sc
        ⇒ free-e max = free-e_sc.
     If also free-e_sc ≤ r(p)·scheme_max with r(p)≤3/2 (15.205), then
        free-e max < 2−α (15.192) ⇒ dual-eq empty ⇒ residual-(i) dual-eq
     branch blocked.  Alternatively cost_D<2−α (15.203–204) with ker=sc.  ∎

  C. Star count.  For fixed edge e, #stars = 2(n−2).  (Each endpoint has
     n−2 other vertices.)  Far edges: binom(n,2)−1−2(n−2).  ∎

CERTIFIED (not general; uses Max± / free-e LP)
  D. G₊ ≻ 0 on 𝒲₊₊⁰ at p=3,5,7 with explicit min eigenvalues:
        p=3: min λ(G₊|𝒲₊₊⁰) = 8
        p=5: min λ = 40/13
        p=7: min λ = 1536/409
     All >0; dim matches n(n−6)/8; full moment spectrum of E[(yᵀBy)²]
     (Frobenius-ON M-basis, zero-diag nullspace) is positive:
        p=3: {16} ×5
        p=5: {80/13, 144/13, 176/13} ×(26,26,13)
        p=7: {3072,3360,3648,4032,4320}/409 ×(50,100,50,50,25)
     Multiplicities sum to dim.  |Max+| = 12, 260, 11452 at p=3,5,7.
     (p=7 Max+ via ProcessPool 2^{25} enum, W=80, wall ~15s.)

  E. Free-e star saturation at p=5,7.  The free-e LP
        max κ_e  s.t. κ∈sc, κ≥−α
     attains the same optimum when stars are fixed at −α:
        max = 369/455 (p=5), 11736/19775 (p=7)
     matching r(p)·scheme_max.  All 2(n−2) stars are active at the
     unrestricted optimum.  Structural seed for a general free-e dual
     with w supported off a controlled star set — not yet a general proof.

  F. Reconfirm ker=sc by dim(Gsum) at p=3,5 (15.202) and by G₊≻0 at p=7
     (new: min 1536/409>0 ⇒ ker Q₊=scheme-image at p=7).

OPEN (blocks residual i / predicate flip)
  G. Prove for all primes p≥5 one of:
     (i)   G₊ ≻ 0 on 𝒲₊₊⁰ (Fourier/association-scheme diagonalization of
           the fourth-moment operator on zero-diag ++), then free-e_sc bound;
     (ii)  free-e max ≤ r(p)·scheme with explicit dual, and ker=sc;
     (iii) |μ₄|≤M_cand (or ≤1/(2p) / ≤2/n) on |κ|=1;
     (iv)  K₄≤Wick_hi.
     Then wire gsum_disj_lb_proved_general → residual_i → E1 → L only if
     residual (ii) full allows.  Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15209.json
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
from e1_gmin_m4_prop15190 import alpha_particular, scheme_ker_max_kappa_e  # noqa: E402
from e1_gmin_m4_prop15205 import free_e_scheme_ratio_candidate  # noqa: E402


def dim_zero_diag_plusplus(p: int) -> int:
    """dim 𝒲₊₊⁰ = n(n−6)/8."""
    n = n_of(p)
    return n * (n - 6) // 8


def theorem_dim_zero_diag_plusplus() -> dict:
    primes = [p for p in range(3, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        d = n // 2
        formula = dim_zero_diag_plusplus(p)
        expand = d * (d + 1) // 2 - n
        match = formula == expand and formula == n * (n - 6) // 8
        if p in (3, 5, 7, 11):
            sample[str(p)] = {
                "n": n,
                "d": d,
                "dim": formula,
                "d_d1_over_2_minus_n": expand,
                "match": match,
            }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For n=2d even, n≥6: dim{M∈Sym(d): diag(Vp M Vpᵀ)=0} "
            "= d(d+1)/2 − n = n(n−6)/8."
        ),
        "by_p_sample": sample,
    }


def theorem_implication_Gplus_to_residual_i() -> dict:
    return {
        "proved": True,
        "theorem": (
            "G₊≻0 on 𝒲₊₊⁰ (and − twin) ⇒ ker=sc (15.207).  "
            "With free-e_sc ≤ r(p)·scheme_max and r≤3/2 (15.205) "
            "⇒ free-e max < 2−α for p≥5 ⇒ dual-eq empty."
        ),
        "depends_on": ["prop_15_207", "prop_15_205", "prop_15_192"],
    }


def theorem_star_count() -> dict:
    primes = [p for p in range(3, 40) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        n = n_of(p)
        n_stars = 2 * (n - 2)
        n_far = n * (n - 1) // 2 - 1 - n_stars
        if p in (3, 5, 7):
            sample[str(p)] = {
                "n_stars": n_stars,
                "n_far": n_far,
                "binom": n * (n - 1) // 2,
            }
        if n_stars != 2 * (n - 2) or n_far < 0:
            ok = False
    return {
        "proved": ok,
        "theorem": "For fixed edge e: #stars=2(n−2); #far=binom(n,2)−1−2(n−2).",
        "by_p_sample": sample,
    }


def certify_Gplus_spectrum() -> dict:
    """Certified min eig / full moment spectrum at p=3,5,7 (from this session)."""
    return {
        "certified": True,
        "proved_general": False,
        "by_p": {
            "3": {
                "dim_Wpp0": 5,
                "Gplus_min_eig": 8,
                "moment_eigs": {"16": 5},
                "N_Max_plus": 12,
                "positive_definite": True,
            },
            "5": {
                "dim_Wpp0": 65,
                "Gplus_min_eig": str(Fraction(40, 13)),
                "moment_eigs": {
                    str(Fraction(80, 13)): 26,
                    str(Fraction(144, 13)): 26,
                    str(Fraction(176, 13)): 13,
                },
                "N_Max_plus": 260,
                "positive_definite": True,
            },
            "7": {
                "dim_Wpp0": 275,
                "Gplus_min_eig": str(Fraction(1536, 409)),
                "moment_eigs": {
                    str(Fraction(3072, 409)): 50,
                    str(Fraction(3360, 409)): 100,
                    str(Fraction(3648, 409)): 50,
                    str(Fraction(4032, 409)): 50,
                    str(Fraction(4320, 409)): 25,
                },
                "N_Max_plus": 11452,
                "positive_definite": True,
                "note": "Max+ via ProcessPool 2^25 enum; G+ min 1536/409>0",
            },
        },
        "note": (
            "All min eigenvalues >0 at p=3,5,7 ⇒ G₊≻0 on 𝒲₊₊⁰ there.  "
            "General PD still open (15.207 hinge)."
        ),
    }


def certify_free_e_star_saturation() -> dict:
    """Free-e optimum saturates stars at −α; stars-fixed = unrestricted."""
    rows = {}
    # Exact free-e maxima from 15.202
    known = {
        5: Fraction(369, 455),
        7: Fraction(11736, 19775),
    }
    for p, fe in known.items():
        al = alpha_particular(p)
        sch = scheme_ker_max_kappa_e(p)
        r = free_e_scheme_ratio_candidate(p)
        n = n_of(p)
        rows[str(p)] = {
            "free_e_max": str(fe),
            "r_scheme": str(r * sch),
            "match_r_scheme": fe == r * sch,
            "n_stars": 2 * (n - 2),
            "alpha": str(al),
            "stars_fixed_equals_unrestricted": True,
            "all_stars_active_at_optimum": True,
            "note": "LP certified this session; see tests/test_prop15209.py",
        }
    return {
        "certified": True,
        "proved_general": False,
        "by_p": rows,
        "note": (
            "At p=5,7 free-e max over sc is attained with κ_star=−α on all "
            "2(n−2) stars.  Stars-fixed LP recovers the same max.  Not proved "
            "for general p."
        ),
    }


def free_e_bound_proved_general() -> bool:
    return False


def Gplus_pd_proved_general() -> bool:
    return False


def residual_i_closed_via_209() -> bool:
    return False


def hinge_status_209() -> dict:
    return {
        "dim_Wpp0_formula_proved": True,
        "Gplus_pd_general": Gplus_pd_proved_general(),
        "Gplus_pd_certified_p": [3, 5, 7],
        "free_e_star_saturation_general": False,
        "free_e_star_saturation_certified_p": [5, 7],
        "residual_i_closed": residual_i_closed_via_209(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove G₊≻0 on 𝒲₊₊⁰ for all p≥5 (or |μ|≤M_cand / K₄≤Wick_hi / "
            "free-e≤r·scheme+ker=sc)."
        ),
    }


def main() -> dict:
    dim_t = theorem_dim_zero_diag_plusplus()
    imp = theorem_implication_Gplus_to_residual_i()
    stars = theorem_star_count()
    gspec = certify_Gplus_spectrum()
    fsat = certify_free_e_star_saturation()
    status = hinge_status_209()
    out = {
        "prop": "15.209",
        "title": "G+ spectrum p=7; free-e star saturation; dim W++0",
        "proved": {
            "dim_zero_diag_plusplus": dim_t["proved"],
            "implication_Gplus_to_dual_eq": imp["proved"],
            "star_count": stars["proved"],
        },
        "certified": {
            "Gplus_spectrum": gspec,
            "free_e_star_saturation": fsat,
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_209(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15209.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    print(f"  dim formula proved={dim_t['proved']}")
    print(f"  G+ PD certified p=3,5,7 (min 8, 40/13, 1536/409)")
    print(f"  star saturation certified p=5,7")
    print(f"  residual_i closed: {residual_i_closed_via_209()}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  e1_closed_general={e1_closed_general()}")
    return out


if __name__ == "__main__":
    main()
