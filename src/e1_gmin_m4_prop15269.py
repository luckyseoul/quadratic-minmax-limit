#!/usr/bin/env python3
"""
Prop 15.269 — Additive Fourier support of Max+; Wick 3-point on ∞-sets;
residual-(i) reduces to a connected 3-point / ẑ-triple bound.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  Paley C on V=F_q ∪ {∞}, q=p², C_{∞,x}=1, C_{xy}=χ(x−y).
  Max+ = {y∈{±1}^n : Cy = p y}, n=q+1.  Write y_∞ and z=y|_{F_q}.
  Additive FT: ẑ(ξ)=∑_{x∈F_q} z(x) ψ(ξx), ψ(u)=exp(2πi Tr_{q/p}(u)/p).
  Quadratic Gauss sum τ=∑_x χ(x) ψ(x) ∈ {±p} (τ²=q, χ(−1)=1).
  Ω := {ξ∈F_q^× : χ̂(ξ)=p}, equivalently χ(ξ)=τ/p.  |Ω|=(q−1)/2.
  For an ∞-set S={∞,a,b,c} with a,b,c distinct in F_q:
        μ(S)=E_{y_∞=+1}[z(a)z(b)z(c)]
  (both sign classes of Max+ contribute equally; 15.268 ν=0 is not used
  for this identity).  Mean E_{+}[z(x)]=C_{∞x}/p=1/p (2-design / 15.189).

PROVED Max+-free (Paley)
  A. Convolution form of Cy=py.
     y_∞=+1 ⇒ ∑_x z(x)=p and (χ∗z)(a)=p z(a)−1 for every a∈F_q,
     where (χ∗z)(a)=∑_{x≠a} χ(a−x) z(x).  ∎

  B. Fourier support.
     χ̂(0)=0 and χ̂(ξ)=χ(ξ) τ for ξ≠0.  Taking FT of (A):
        ẑ(0)=p, and ẑ(ξ)(χ̂(ξ)−p)=0 for ξ≠0.
     Hence ẑ is supported on {0}∪Ω.  Parseval: ∑|ẑ|²=q∑|z|²=q²=p⁴.  ∎
     (Certified: every y_∞=+1 vector at p=5,7 has this support; τ=−p at
     p=5 so Ω=nonsquares; τ=+p at p=7 so Ω=squares.)

  C. Translation invariance of the y_∞=+1 slice.
     If z satisfies (A) then so does z(·−t) for every t∈F_q.  Uniform
     measure on that slice is translation-invariant.  ∎

  D. Two-point of ẑ (only the 2-design).
     E_{+}[z(x)z(w)]=C_{xw}/p (x≠w) and 1 (x=w).
     Hence for ξ,η≠0:
        E[ẑ(ξ)ẑ(η)] = q 1_{η=−ξ} (1 + χ̂(ξ)/p).
     On Ω one has χ̂=p, so E[ẑ(ξ)ẑ(−ξ)]=2q=2p²; off Ω the factor
     (1−1)=0, consistent with (B).  Also E[ẑ(ξ)]=0 for ξ≠0.  ∎

  E. Wick 3-point on ∞-sets.
     Write z=1/p+ζ.  Then E[ζ]=0 and
        E[z(a)z(b)z(c)] = Wick + κ₃,
     Wick = (C_{ab}+C_{ac}+C_{bc})/p² − 2/p³ = κ(S)/p² − 2/p³
     (κ(S)=C_{∞a}C_{bc}+C_{∞b}C_{ac}+C_{∞c}C_{ab}=C_{ab}+C_{ac}+C_{bc}).
     The connected remainder is
        κ₃=E[ζ(a)ζ(b)ζ(c)].  ∎

  F. |Wick|≤1/(2p) on |κ|=1, all primes p≥5.
     |κ|=1 ⇒ |Wick|≤1/p²+2/p³=(p+2)/p³.
     (p+2)/p³ ≤ 1/(2p) ⇔ p²−2p−4≥0, true on [5,∞).  ∎
     Room for the connected part:
        1/(2p)−(p+2)/p³ = (p²−2p−4)/(2p³) > 0.  ∎

  G. Frequency support of the 3-point of ẑ.
     By (C), E[ẑ(ξ)ẑ(η)ẑ(ζ)]=0 unless ξ+η+ζ=0.
     Hence (with ζ=z−1/p = q^{-1}∑_{ξ∈Ω} ẑ(ξ) e_{−ξ})
        κ₃(a,b,c)=q^{-3} ∑_{ξ,η,θ∈Ω, ξ+η+θ=0}
            E[ẑ(ξ)ẑ(η)ẑ(θ)] ψ(−ξa−ηb−θc).  ∎

  H. Collision types in Ω.
     (F_q^*)^2 acts regularly on Ω by multiplication and preserves
     ξ+η+θ=0.  Double collisions are (ξ,ξ,−2ξ) with −2Ω=Ω
     (χ(−2)=χ(2)=1 on F_{p²}).  Distinct triples are classified by
     r=η/ξ ∈ (F_q^*)^2 with 1+r a square, r∉{1,−2,−1/2}.  ∎

  I. Sufficient criterion.
     If |κ₃|≤(p²−2p−4)/(2p³) on every |κ|=1 ∞-set, for all primes p≥5,
     then |μ|≤1/(2p) on those sets.  PSL is 3-transitive, so every
     |κ|=1 four-set is Aut-equivalent to an ∞-set, and |μ| is Aut-
     invariant (15.268 E with ε=±1).  Hence residual-(i) via 15.176.  ∎

CERTIFIED (not general)
  J. p=5: 130 vectors y_∞=+1; two ẑ-profiles (30 halfspace, 100 full-Ω);
     E[ẑ(ξ)ẑ(η)ẑ(θ)] equals −500/13 on doubles and −2500/13 on distinct
     triples (one orbit).  max|μ|=3/65<1/10 on |κ|=1.
  K. p=7: 5726 vectors y_∞=+1; doubles give −31556/409; distinct triples
     take two values −48020/409 and −89180/409 (two r-classes).
     max|μ|=109/2863<1/14 on |κ|=1.
  L. The Ω-triple exponential sum
        S(u,v)=∑_{ξ,η∈Ω, −ξ−η∈Ω} ψ(ξu+ηv)
     for u,v,u−v≠0 takes values of size O(p) (max 4p at p=5, 26 at p=7).

OPEN (blocks residual i / predicate flip)
  M. Prove |κ₃|≤(p²−2p−4)/(2p³) for all primes p≥5, e.g. by a Hasse
     bound on S together with a Max+-free bound on
        A(ξ,η,θ):=E[ẑ(ξ)ẑ(η)ẑ(θ)]
     on Ω-triples (census: |A|≤2p³ at p=5,7), **or** envelope /
     reflection / Es4 / ker=sc.  Soft-close forbidden.

Until M: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15269.json
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
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15268 import nu_zero_on_kappa1_proved_general  # noqa: E402


def wick_infty(p: int, kappa: int) -> Fraction:
    """κ/p² − 2/p³ on an ∞-set with matching-sum κ."""
    return Fraction(kappa, p * p) - Fraction(2, p**3)


def wick_abs_majorant(p: int) -> Fraction:
    return Fraction(p + 2, p**3)


def kappa3_room(p: int) -> Fraction:
    """1/(2p) − (p+2)/p³ = (p²−2p−4)/(2p³)."""
    return Fraction(p * p - 2 * p - 4, 2 * p**3)


def theorem_wick_le_half_p() -> dict:
    primes = [p for p in range(5, 120) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        maj = wick_abs_majorant(p)
        half = Fraction(1, 2 * p)
        room = kappa3_room(p)
        if not (maj <= half and room > 0 and wick_infty(p, 1) + wick_infty(p, -1) == Fraction(-4, p**3)):
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "wick_maj": str(maj),
                "half_p": str(half),
                "room": str(room),
                "wick_plus": str(wick_infty(p, 1)),
                "wick_minus": str(wick_infty(p, -1)),
            }
    return {
        "proved": ok,
        "theorem": (
            "On |κ|=1 ∞-sets: |Wick|=|κ/p²−2/p³|≤(p+2)/p³≤1/(2p) for all "
            "primes p≥5 (p²−2p−4≥0).  Room for κ₃ is (p²−2p−4)/(2p³)."
        ),
        "by_p_sample": sample,
        "depends_on": ["2_design_15_189", "kappa_infty_matchings"],
    }


def theorem_fourier_support_structure() -> dict:
    return {
        "proved": True,
        "theorem": (
            "y∈Max+, y_∞=+1 ⇒ ẑ supported on {0}∪Ω with ẑ(0)=p and "
            "∑|ẑ|²=p⁴, where Ω={ξ: χ̂(ξ)=p}.  Proof: FT of χ∗z=pz−1."
        ),
        "depends_on": ["Cy_eq_py", "paley_kernel_convolution", "gauss_sum_eval"],
    }


def theorem_zhat_two_point() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E[ẑ(ξ)ẑ(η)]=2p² if η=−ξ∈Ω, else 0 (ξ,η≠0).  "
            "Proof: expand E[z(x)z(w)]=C_{xw}/p and evaluate Gauss sums."
        ),
        "depends_on": ["2_design_15_189", "theorem_fourier_support_structure"],
    }


def theorem_kappa3_criterion() -> dict:
    t = theorem_wick_le_half_p()
    return {
        "proved": t["proved"],
        "hypothesis_proved_general": False,
        "theorem": (
            "If |κ₃|≤(p²−2p−4)/(2p³) on every |κ|=1 ∞-set for all primes "
            "p≥5, then |μ|≤1/(2p) on all |κ|=1 four-sets (3-transitive Aut "
            "+ 15.176).  The κ₃ bound is OPEN."
        ),
        "depends_on": [
            "theorem_wick_le_half_p",
            "PSL_3_transitive",
            "prop_15_176_farkas",
        ],
    }


def certify_support_p5() -> dict:
    """Finite check of (B) at p=5 from the Max+ cache (not a general proof)."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    path = Path("/tmp/maxplus_p5.npy")
    if not path.exists():
        return {"certified": False, "reason": "no cache"}
    Y = np.load(path)
    C = paley_conference_prime_power(p)
    # Cy=py
    ev = np.max(np.abs(C @ Y.T - p * Y.T))
    Z = Y[Y[:, 0] > 0, 1:]
    # ẑ(0)=sum z
    zsum = Z.sum(axis=1)
    return {
        "certified": bool(ev < 1e-6 and np.allclose(zsum, p)),
        "n_yinf_plus": int(Z.shape[0]),
        "Cy_residual": float(ev),
        "zsum_ok": bool(np.allclose(zsum, p)),
    }


def residual_i_closed_via_269() -> bool:
    return False


def hinge_status_269() -> dict:
    return {
        "fourier_support": True,
        "zhat_two_point": True,
        "wick_le_half_p": theorem_wick_le_half_p()["proved"],
        "kappa3_criterion": True,
        "nu_zero": nu_zero_on_kappa1_proved_general(),
        "kappa3_bound": False,
        "residual_i_closed": residual_i_closed_via_269(),
        "open": (
            "Need |κ₃|≤(p²−2p−4)/(2p³) on |κ|=1 ∞-sets "
            "(bound E[ẑ(ξ)ẑ(η)ẑ(θ)] on Ω-triples, or another listed hinge)."
        ),
    }


def main() -> None:
    t_w = theorem_wick_le_half_p()
    t_f = theorem_fourier_support_structure()
    t_2 = theorem_zhat_two_point()
    t_c = theorem_kappa3_criterion()
    cert = certify_support_p5()
    h = hinge_status_269()
    out = {
        "prop": "15.269",
        "title": "Fourier support of Max+; Wick 3-point; κ₃ criterion",
        "theorems": {
            "fourier_support": t_f,
            "zhat_two_point": t_2,
            "wick_le_half_p": t_w,
            "kappa3_criterion": t_c,
        },
        "certified_p5_support": cert,
        "hinge": h,
        "predicates": {
            "residual_i_closed_via_269": residual_i_closed_via_269(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "e1_closed_general": e1_closed_general(),
        },
        "soft_close": False,
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15269.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}")
    print(f"  wick={t_w['proved']} support={t_f['proved']} 2pt={t_2['proved']}")
    print(f"  kappa3_hyp={t_c['hypothesis_proved_general']} residual_i={h['residual_i_closed']}")


if __name__ == "__main__":
    main()
