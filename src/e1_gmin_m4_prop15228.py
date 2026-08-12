#!/usr/bin/env python3
"""
Prop 15.228 — Boolean diamond |μ|+|ν|≤1; ν=Tμ/(4p); Tμ_part=c·star;
on |κ|=1: star=0⇒ν=δ_{+}−δ_{-}; Aut-SOS with diamond; residual OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.177, 15.215, 15.225–227)
  m₄⁺=μ+ν, m₄⁻=μ−ν with ν=½(m₄⁺−m₄⁻), μ=½(m₄⁺+m₄⁻).
  Masters: (4pI−T)m₄⁺=4κ/p, (4pI+T)m₄⁻=4κ/p (15.177 A–B).
  μ=μ_part+δ, δ∈E_{±4p}^G, δ=δ_{+}+δ_{-}, Tδ=4p(δ_{+}−δ_{-}).
  star=0 on |κ|=1, |star|=4 on |κ|=3 (15.215); Tκ=−6 star; Tφ=(n+2)star.

PROVED Max+-free
  A. Coupling identities from the two masters.
        Tμ = 4p ν,     ν = Tμ/(4p),
        μ = κ/p² + Tν/(4p),
     and (16p²I−T²)μ=16κ (master, 15.177 C).  ∎

  B. Boolean diamond.
     |m₄⁺_S|≤1 and |m₄⁻_S|≤1 for every 4-set S
     ⇔ |μ_S+ν_S|≤1 and |μ_S−ν_S|≤1
     ⇔ |μ_S| + |ν_S| ≤ 1.
     Proof: for real scalars max(|a+b|,|a−b|)=|a|+|b|.  ∎

  C. Closed form of Tμ_part.
     μ_part=aκ+bφ with a=(p²−1)/den, b=−2/den, den=p²(p²−5).
        Tμ_part = a Tκ + b Tφ = (−6a + b(n+2)) star
                = −8/(p²−5) · star.
     Proof: −6a+b(n+2)=[−6(p²−1)−2(n+2)]/den with n=p²+1 equals
     −8p²/den=−8/(p²−5).  ∎

  D. On |κ|=1: star=0 ⇒ (Tμ_part)_S=0 ⇒
        ν_S = (δ_{+} − δ_{-})_S,
        μ_S = μ_part_S + (δ_{+} + δ_{-})_S,
     and the diamond becomes
        |μ_part_S + δ_S| + |(δ_{+} − δ_{-})_S| ≤ 1.
     (Global vectors: δ_{+}∈E_{+4p}, δ_{-}∈E_{−4p}, so δ_{+}=δ_{-} as vectors
     forces δ=0; pointwise equality on |κ|=1 alone does not.)  ∎

  E. L² orthogonality upgrade.
     ⟨Tμ_part, Tδ⟩ = ⟨T²μ_part, δ⟩=0 because T²μ_part∈span{κ,φ} (from
     T²κ, T²φ) and δ⊥span{κ,φ} (15.225 H).  Hence
        ‖Tμ‖₂² = ‖Tμ_part‖₂² + ‖Tδ‖₂²
                = ‖Tμ_part‖₂² + 16 p² ‖δ‖₂²,
        ‖ν‖₂² = ‖Tμ_part‖₂²/(16 p²) + ‖δ‖₂²,
     with ‖Tμ_part‖₂² = 64/(p²−5)² · ‖star‖₂² and
     ‖star‖₂²=n(n−1)(n−2)(n−6)/6.  ∎

  F. Diamond L² consequence (weak).
     μ_S²+ν_S² ≤ (|μ_S|+|ν_S|)² ≤ 1 ⇒ ∑(μ²+ν²) ≤ binom(n,4).
     Combined with (E) and ‖μ‖₂²=‖μ_part‖₂²+‖δ‖₂² (15.226):
        2‖δ‖₂² ≤ binom(n,4) − ‖μ_part‖₂² − ‖Tμ_part‖₂²/(16p²).
     Positive but far above R-room / room_hyp/24 (not sufficient alone).  ∎

  G. Aut-SOS hinge with diamond (structure).
     Residual-(i) closes if a parameter-dependent certificate on
     (δ_{+},δ_{-})∈E_{+4p}^G × E_{−4p}^G proves
        |μ_part + δ_{+} + δ_{-}| ≤ 1/(2p)   on every |κ|=1 orbit,
     subject to the diamond (D), linear constraints ⟨δ,κ⟩=⟨δ,φ⟩=0,
     Tδ=4p(δ_{+}−δ_{-}), and Tν=4p(μ_part−κ/p²+δ) from (A).
     The diamond is a new quadratic constraint for SOS/Schur relative to
     15.225–227.  Bound still OPEN.  ∎

CERTIFIED
  H. C matches Fraction for all primes 5≤p≤100; diamond L² slack checked.
  I. Census p=5: m₄⁺=m₄⁻ on |κ|=1 (15.177 F) ⇒ ν≡0 there ⇒ diamond is
     |μ|≤1; actual max|μ|=3/65≤1/(2p).

OPEN (blocks residual i / predicate flip)
  J. Prove |μ|≤1/(2p) on |κ|=1 for all primes p≥5 using Aut-SOS+diamond
     on E_{±4p}^G, **or** Path C / K₄ / free-e+ker=sc / m4₂.  Soft-close forbidden.

Until J: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15228.json
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
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15217 import delta_room_for_R  # noqa: E402
from e1_gmin_m4_prop15226 import mu_part_norm_sq  # noqa: E402


def T_mupart_star_coeff(p: int) -> Fraction:
    """Tμ_part = c·star with c = −8/(p²−5)."""
    return Fraction(-8, p * p - 5)


def T_mupart_star_coeff_via_ab(p: int) -> Fraction:
    den = p * p * (p * p - 5)
    a = Fraction(p * p - 1, den)
    b = Fraction(-2, den)
    n = n_of(p)
    return -6 * a + b * (n + 2)


def star_norm_sq(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2) * (n - 6), 6)


def T_mupart_norm_sq(p: int) -> Fraction:
    c = T_mupart_star_coeff(p)
    return c * c * star_norm_sq(p)


def diamond_L2_delta_ub(p: int) -> Fraction:
    """2‖δ‖₂² ≤ binom − ‖μ_part‖₂² − ‖Tμ_part‖₂²/(16p²)."""
    B = comb(n_of(p), 4)
    return (
        Fraction(B)
        - mu_part_norm_sq(p)
        - T_mupart_norm_sq(p) / (16 * p * p)
    ) / 2


def theorem_coupling_masters() -> dict:
    return {
        "proved": True,
        "theorem": (
            "From (4pI−T)m₄⁺=4κ/p and (4pI+T)m₄⁻=4κ/p with "
            "μ=½(m₄⁺+m₄⁻), ν=½(m₄⁺−m₄⁻): Tμ=4p ν, μ=κ/p²+Tν/(4p), "
            "and (16p²I−T²)μ=16κ."
        ),
        "depends_on": ["prop_15_177_masters"],
    }


def theorem_boolean_diamond() -> dict:
    return {
        "proved": True,
        "theorem": (
            "|m₄⁺_S|≤1 and |m₄⁻_S|≤1 ⇔ |μ_S|+|ν_S|≤1 for every 4-set S "
            "(real scalar identity max(|a+b|,|a−b|)=|a|+|b|)."
        ),
    }


def theorem_T_mupart_closed() -> dict:
    ok = True
    sample = {}
    for p in range(5, 201):
        if not is_prime(p):
            continue
        c1 = T_mupart_star_coeff(p)
        c2 = T_mupart_star_coeff_via_ab(p)
        if c1 != c2:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31):
            sample[str(p)] = {
                "c": str(c1),
                "T_mupart_norm_sq": str(T_mupart_norm_sq(p)),
            }
    return {
        "proved": ok,
        "theorem": (
            "Tμ_part=−8/(p²−5)·star for all primes p≥5.  "
            "On |κ|=1 (star=0): Tμ_part=0."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "prop_15_215_Tkappa_star",
            "prop_15_215_Tphi_star",
            "prop_15_224_mu_part",
        ],
    }


def theorem_kappa1_nu_delta() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On every |κ|=1 four-set: ν_S=(δ_{+}−δ_{-})_S and "
            "|μ_part_S+δ_S|+|ν_S|≤1 (diamond + Tμ_part=0)."
        ),
        "depends_on": [
            "theorem_T_mupart_closed",
            "theorem_boolean_diamond",
            "prop_15_215_star_on_kappa1",
            "prop_15_225_delta_split",
        ],
    }


def theorem_L2_T_orthogonality() -> dict:
    return {
        "proved": True,
        "theorem": (
            "⟨Tμ_part,Tδ⟩=0 ⇒ ‖Tμ‖₂²=‖Tμ_part‖₂²+16p²‖δ‖₂² and "
            "‖ν‖₂²=‖Tμ_part‖₂²/(16p²)+‖δ‖₂²."
        ),
        "depends_on": [
            "theorem_T_mupart_closed",
            "prop_15_225_L2_orth",
            "prop_15_224_T2_phi",
        ],
    }


def theorem_diamond_L2_delta_ub() -> dict:
    ok = True
    sample = {}
    for p in range(5, 80):
        if not is_prime(p):
            continue
        ub = diamond_L2_delta_ub(p)
        if ub <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17):
            sample[str(p)] = {
                "delta2_ub": str(ub),
                "R_room": str(delta_room_for_R(p)),
                "ub_gt_R_room": ub > delta_room_for_R(p),
                "sufficient_for_R": False,
            }
    return {
        "proved": ok,
        "sufficient_for_residual_i": False,
        "theorem": (
            "2‖δ‖₂² ≤ binom(n,4)−‖μ_part‖₂²−‖Tμ_part‖₂²/(16p²) from the "
            "diamond L² inequality.  Positive for p≥5 but ≫ R-room "
            "(not a residual-(i) hinge by itself)."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_boolean_diamond",
            "theorem_L2_T_orthogonality",
            "prop_15_226_mu_part_norm_sq",
        ],
    }


def theorem_Aut_SOS_diamond_hinge() -> dict:
    return {
        "proved_structure": True,
        "bound_proved_general": False,
        "theorem": (
            "Aut-SOS/Schur on (δ_{+},δ_{-})∈E_{+4p}^G×E_{−4p}^G with diamond "
            "|μ_part+δ_{+}+δ_{-}|+|δ_{+}−δ_{-}|≤1 on |κ|=1, linear orth "
            "⟨δ,κ⟩=⟨δ,φ⟩=0, and Tδ=4p(δ_{+}−δ_{-}), would prove "
            "|μ|≤1/(2p) on |κ|=1.  Certificate still OPEN."
        ),
        "depends_on": [
            "theorem_kappa1_nu_delta",
            "prop_15_226_delta_G_invariant",
            "prop_15_227_onesided_form",
        ],
        "target": str(mu4_threshold(5)) + " at p=5 (general: 1/(2p))",
    }


def residual_i_closed_via_228() -> bool:
    return False


def hinge_status_228() -> dict:
    return {
        "coupling_masters": theorem_coupling_masters()["proved"],
        "boolean_diamond": theorem_boolean_diamond()["proved"],
        "T_mupart_closed": theorem_T_mupart_closed()["proved"],
        "kappa1_nu_delta": theorem_kappa1_nu_delta()["proved"],
        "L2_T_orth": theorem_L2_T_orthogonality()["proved"],
        "diamond_L2_ub": theorem_diamond_L2_delta_ub()["proved"],
        "Aut_SOS_structure": theorem_Aut_SOS_diamond_hinge()[
            "proved_structure"
        ],
        "Aut_SOS_bound": theorem_Aut_SOS_diamond_hinge()[
            "bound_proved_general"
        ],
        "residual_i_closed_via_228": residual_i_closed_via_228(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Aut-SOS+diamond on E_pm4p^G for |μ|≤1/(2p); "
            "or Path C / K4 / free-e+ker=sc."
        ),
    }


def main() -> dict:
    a = theorem_coupling_masters()
    b = theorem_boolean_diamond()
    c = theorem_T_mupart_closed()
    d = theorem_kappa1_nu_delta()
    e = theorem_L2_T_orthogonality()
    f = theorem_diamond_L2_delta_ub()
    g = theorem_Aut_SOS_diamond_hinge()
    status = hinge_status_228()
    out = {
        "prop": "15.228",
        "title": (
            "Boolean diamond |μ|+|ν|≤1; Tμ_part=c·star; |κ|=1⇒ν=δ_{+}−δ_{-}; "
            "Aut-SOS+diamond; residual OPEN"
        ),
        "proved": {
            "coupling_masters": a["proved"],
            "boolean_diamond": b["proved"],
            "T_mupart_closed": c["proved"],
            "kappa1_nu_delta": d["proved"],
            "L2_T_orth": e["proved"],
            "diamond_L2_ub": f["proved"],
            "Aut_SOS_structure": g["proved_structure"],
            "Aut_SOS_bound": g["bound_proved_general"],
            "residual_i_closed": residual_i_closed_via_228(),
        },
        "algebra": {
            "T_mupart": "-8/(p^2-5) * star",
            "nu": "T mu / (4p)",
            "diamond": "|mu| + |nu| <= 1",
            "on_kappa1": "nu = delta_+ - delta_-, T mu_part = 0",
            "c_p5": str(T_mupart_star_coeff(5)),
            "c_p7": str(T_mupart_star_coeff(7)),
            "diamond_L2_delta_ub_p5": str(diamond_L2_delta_ub(5)),
            "R_room_p5": str(delta_room_for_R(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_228(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Diamond + Tμ_part + |κ|=1 ν=δ_{+}−δ_{-} proved Max+-free.  "
            "Aut-SOS certificate still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15228.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.228  residual_i={residual_i_closed_via_228()}")
    print(
        f"  T_mupart={c['proved']} diamond={b['proved']} "
        f"kappa1={d['proved']} Aut_SOS_bound={g['bound_proved_general']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
