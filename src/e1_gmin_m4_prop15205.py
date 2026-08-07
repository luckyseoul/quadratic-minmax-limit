#!/usr/bin/env python3
"""
Prop 15.205 — Residual-(i) sufficient targets: M_cand ≤ 1/(2p);
free-e / scheme ratio algebra; residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP
  Residual (i) Farkas (15.176): disj Gsum ≥ −1/p on all matching pairs
  ⇔ |μ₄| ≤ 1/(2p) on every |κ|=1 four-set.
  Free-e dual-eq empty if free-e max κ_e < 2−α (15.182/200).
  M_cand := (p−2)/(p(2p+3)) from 15.74 (Path-C companion).
  Scheme max = α(n−2)/2; (3/2)·scheme < 2−α for p≥5 (15.192).

PROVED Max+-free (Fraction algebra)
  A. M_cand ≤ 1/(2p) for every prime p≥3.
        (p−2)/(p(2p+3)) ≤ 1/(2p)
     ⇔ 2(p−2) ≤ 2p+3  ⇔  −4 ≤ 3, true.  ∎
     Hence |μ₄| ≤ M_cand on |κ|=1 for all p≥5
     ⇒ |μ₄| ≤ 1/(2p) ⇒ Gsum ≥ −1/p ⇒ residual-(i) Farkas
     ⇒ gsum_disj_lb_proved_general may flip when wired.  ∎

  B. M_cand ≤ M_mid ≤ L_abs for p≥5, with
        M_mid = (p−2)/(2p(p+1)),   L_abs = (p−2)/(2p²)
     (15.74 algebra; rechecked Fraction).  ∎

  C. Free-e / scheme ratio candidate
        r(p) := (7p² − 3p + 4)/(5p² − 3p + 2)
     satisfies r(p) ≤ 3/2 for all primes p≥5:
        2(7p²−3p+4) ≤ 3(5p²−3p+2)
     ⇔ 14p²−6p+8 ≤ 15p²−9p+6  ⇔  0 ≤ p² − 3p − 2.
     At p=5: 8>0; derivative 2p−3>0 on [5,∞).  ∎
     Combined with (3/2)·scheme_max < 2−α (15.192):
        free-e max ≤ r(p)·scheme_max  ⇒ free-e max < 2−α
     for all primes p≥5 (no ker=sc needed if the free-e max is the
     true one; if only over sc, still need ker=sc).  ∎

  D. r(p) ≥ 1 for p≥2 (7p²−3p+4 ≥ 5p²−3p+2 ⇔ 2p²+2≥0).
     So the candidate is a genuine inflation of scheme.  ∎

  E. Asymptotic: r(p) → 7/5 = 1.4 < 3/2; M_cand ∼ 1/(2p).  ∎

CERTIFIED (not general)
  F. free-e max / scheme_max = r(p) exactly at p=5 (41/28) and p=7 (163/113)
     (15.202 fractions).  Matches the candidate closed form.
  G. max|μ₄| ≤ M_cand at p=5 (sharp 3/65) and p=7 (109/2863 < 5/119)
     by Max± census (15.74/188).
  H. M_cand ≤ 2/n at p=5 only; for p≥7 use 1/(2p) path (A), not 2/n.

OPEN (blocks residual i / predicate flip)
  I. Prove Max+-free for all primes p≥5 one of:
     (i)   |μ₄| ≤ M_cand on every |κ|=1 four-set
           (Path-C companion; star·S1≤0 + S3 budget, 15.74–78),
     (ii)  free-e max ≤ r(p)·scheme_max
           (or S*_sc ≤ (n−2) r(p)/2, or ≤ n−2), and ker=sc if dual is sc-only,
     (iii) K₄ ≤ Wick_hi, or cost_D < 2−α with ker=sc (15.203–204).
     Then wire gsum_disj_lb_proved_general → residual_i → E1 → L
     only if residual (ii) full allows.  Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15205.json
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
from e1_gmin_m4_prop15176 import theorem_minus_1_over_p_suffices  # noqa: E402
from e1_gmin_m4_prop15190 import alpha_particular, scheme_ker_max_kappa_e  # noqa: E402


def M_cand(p: int) -> Fraction:
    return Fraction(p - 2, p * (2 * p + 3))


def M_mid(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * (p + 1))


def L_abs(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * p)


def free_e_scheme_ratio_candidate(p: int) -> Fraction:
    """r(p) = (7p² − 3p + 4)/(5p² − 3p + 2)."""
    return Fraction(7 * p * p - 3 * p + 4, 5 * p * p - 3 * p + 2)


def theorem_Mcand_le_mu4_thr() -> dict:
    """M_cand ≤ 1/(2p) for all primes p≥3."""
    primes = [p for p in range(3, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        mc = M_cand(p)
        thr = Fraction(1, 2 * p)
        beats = mc <= thr
        # 2(p-2) ≤ 2p+3 ⇔ -4 ≤ 3
        if p in (3, 5, 7, 11, 13):
            sample[str(p)] = {
                "M_cand": str(mc),
                "mu4_thr_1_over_2p": str(thr),
                "M_cand_le_thr": beats,
                "two_over_n": str(Fraction(2, n_of(p))),
                "M_cand_le_two_over_n": mc <= Fraction(2, n_of(p)),
            }
        if not beats:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "M_cand=(p−2)/(p(2p+3)) ≤ 1/(2p) for all primes p≥3, "
            "since 2(p−2) ≤ 2p+3 ⇔ −4≤3.  Hence |μ₄|≤M_cand ⇒ "
            "|μ₄|≤1/(2p) ⇒ Gsum≥−1/p ⇒ residual-(i) Farkas (15.176)."
        ),
        "by_p_sample": sample,
        "depends_on": ["farkas_15_176"],
    }


def theorem_Mcand_le_Mmid_le_Labs() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        c, m, L = M_cand(p), M_mid(p), L_abs(p)
        chain = c <= m <= L
        if p in (5, 7, 11):
            sample[str(p)] = {
                "M_cand": str(c),
                "M_mid": str(m),
                "L_abs": str(L),
                "chain": chain,
            }
        if not chain:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For primes p≥5: M_cand ≤ M_mid ≤ L_abs with "
            "M_mid=(p−2)/(2p(p+1)), L_abs=(p−2)/(2p²)."
        ),
        "by_p_sample": sample,
    }


def theorem_ratio_le_three_halves() -> dict:
    """r(p) ≤ 3/2 for all primes p≥5."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        r = free_e_scheme_ratio_candidate(p)
        poly = p * p - 3 * p - 2
        beats = r <= Fraction(3, 2) and poly > 0
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "r": str(r),
                "r_float": float(r),
                "le_3_2": r <= Fraction(3, 2),
                "poly_p2_minus_3p_minus_2": poly,
            }
        if not beats:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "r(p)=(7p²−3p+4)/(5p²−3p+2) ≤ 3/2 for all primes p≥5 "
            "⇔ p²−3p−2 ≥ 0 (value 8 at p=5; increasing on [5,∞)).  "
            "Hence free-e ≤ r·scheme and (3/2)·scheme < 2−α (15.192) "
            "⇒ free-e < 2−α for p≥5."
        ),
        "by_p_sample": sample,
    }


def theorem_ratio_ge_one() -> dict:
    primes = [p for p in range(2, 40) if is_prime(p)]
    ok = all(free_e_scheme_ratio_candidate(p) >= 1 for p in primes)
    return {
        "proved": ok,
        "theorem": "r(p) ≥ 1 for primes p≥2 (⇔ 2p²+2 ≥ 0).",
    }


def certify_ratio_at_small_p() -> dict:
    """Census: free-e/scheme = r(p) at p=5,7 (known fractions from 15.202)."""
    known_free = {
        5: Fraction(369, 455),
        7: Fraction(11736, 19775),
    }
    rows = {}
    for p, fe in known_free.items():
        sch = scheme_ker_max_kappa_e(p)
        ratio = fe / sch
        r = free_e_scheme_ratio_candidate(p)
        rows[str(p)] = {
            "free_e_max": str(fe),
            "scheme_max": str(sch),
            "ratio": str(ratio),
            "r_candidate": str(r),
            "match": ratio == r,
        }
    return {
        "certified": all(v["match"] for v in rows.values()),
        "proved_general": False,
        "by_p": rows,
        "note": "Exact match at p=5,7; not a general proof.",
    }


def certify_mu_le_Mcand_small() -> dict:
    """Census max|μ| ≤ M_cand at p=5,7 from 15.188 fractions."""
    # p=5: 3/65 = M_cand; p=7: 109/2863 < 5/119 = M_cand
    rows = {
        "5": {
            "max_abs_mu": "3/65",
            "M_cand": str(M_cand(5)),
            "le": Fraction(3, 65) <= M_cand(5),
            "sharp": Fraction(3, 65) == M_cand(5),
        },
        "7": {
            "max_abs_mu": "109/2863",
            "M_cand": str(M_cand(7)),
            "le": Fraction(109, 2863) <= M_cand(7),
            "sharp": False,
        },
    }
    return {
        "certified": all(v["le"] for v in rows.values()),
        "proved_general": False,
        "by_p": rows,
    }


def free_e_bound_proved_general() -> bool:
    return False


def mu_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_205() -> bool:
    return free_e_bound_proved_general() or mu_bound_proved_general()


def hinge_status_205() -> dict:
    return {
        "Mcand_le_1_over_2p": True,
        "ratio_le_3_2": True,
        "Mcand_bound_general": False,
        "ratio_free_e_bound_general": False,
        "residual_i_closed": residual_i_closed_via_205(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove |μ₄|≤M_cand on |κ|=1 for all p≥5, or free-e≤r(p)·scheme "
            "(+ker=sc if sc-dual), or K₄≤Wick_hi / cost_D+ker=sc."
        ),
    }


def main() -> dict:
    mc = theorem_Mcand_le_mu4_thr()
    chain = theorem_Mcand_le_Mmid_le_Labs()
    ratio = theorem_ratio_le_three_halves()
    rge1 = theorem_ratio_ge_one()
    cert_r = certify_ratio_at_small_p()
    cert_m = certify_mu_le_Mcand_small()
    farkas_ok = bool(theorem_minus_1_over_p_suffices())
    out = {
        "title": (
            "Prop 15.205 M_cand≤1/(2p); free-e/scheme ratio algebra; "
            "residual (i) OPEN"
        ),
        "proved": {
            "Mcand_le_1_over_2p": mc["proved"],
            "Mcand_le_Mmid_le_Labs": chain["proved"],
            "ratio_le_3_2": ratio["proved"],
            "ratio_ge_1": rge1["proved"],
            "minus_1_over_p_farkas": farkas_ok,
            "Mcand_mu_bound_general": False,
            "ratio_free_e_bound_general": False,
            "residual_i_closed": residual_i_closed_via_205(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "Mcand_le_thr": mc,
            "Mcand_chain": chain,
            "ratio_le_3_2": ratio,
            "ratio_ge_1": rge1,
        },
        "certified": {
            "ratio_match_p5_p7": cert_r,
            "mu_le_Mcand_p5_p7": cert_m,
        },
        "hinge": hinge_status_205(),
        "F3": (
            "Threshold algebra is Max+-free, but |μ|≤M_cand and free-e≤r·scheme "
            "are not proved for all p≥5.  Do not flip residual/E1/L."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15205.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.205 M_cand≤1/(2p); free-e/scheme ratio algebra")
    print(f"  M_cand ≤ 1/(2p) p≥3: {mc['proved']}")
    print(f"  M_cand ≤ M_mid ≤ L: {chain['proved']}")
    print(f"  r(p) ≤ 3/2 p≥5: {ratio['proved']}")
    print(f"  r matches free-e/scheme p=5,7: {cert_r['certified']}")
    print(f"  |μ|≤M_cand p=5,7: {cert_m['certified']}")
    print(f"  residual_i closed: {residual_i_closed_via_205()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
