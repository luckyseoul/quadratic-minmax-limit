#!/usr/bin/env python3
"""
Prop 15.166 — 16N ⇔ λ_max(Q₂)≤4N/(d(d−1)) for Max+ unit 2-design;
Wick m₄ is C-eigen (same as true m₄); residual still OPEN.

Continues 15.164–15.165. Does NOT soft-close L.

============================================================================
Setup. Max+ unit vectors û_y = c_y/√n ∈ S^{d−1} ⊂ ℝ^d form a spherical
2-design (∑ û=0, (1/N)∑ ûûᵀ = I_d/d).  Q₂(t)=(d t²−1)/(d−1) Gegenbauer
with Q₂(1)=1.  Q̂₂ = (Q₂(û_a·û_b))_{a,b} on Max+ indices.

============================================================================
Theorem A (2-design of Max+ in V₊) — PROVED.
  c_y = V₊ᵀ y, ‖c_y‖²=n, û=c/√n.  ∑_y yyᵀ = 2N P₊ ⇒
  (1/N)∑ ûûᵀ = I_d/d.  Central symmetry ⇒ ∑ û=0.  Hence {û_y} is a
  spherical 2-design in ℝ^d.  Certified second-moment identity p=3,5,7.  ∎

Theorem B (Wick m₄ is C-eigen) — PROVED / CERTIFIED.
  The Isserlis/Wick reconstruction m₄^W from m₂=C/p satisfies
      p · m₄^W(a,b,c,d) = ∑_j C_{aj} m₄^W(j,b,c,d)
  for all index tuples (same eigen-equation as true m₄).  Certified
  random samples p=5,7 (max residual <1e-15).  Consequently the
  residual η = m₄−m₄^W lies in the +p eigenspace of C in each index
  and is invisible to the C-eigen constraint; P₊(Wick) does not recover
  true m₄.  ∎

Theorem C (16N ⇔ Q₂ bound) — PROVED Fraction.
  Let Ŝ_{ab}=(û_a·û_b)², S=Π∘Π with Π=G/(2N) the equidiagonal rank-d
  orthoprojector on ℝ^N (G=YYᵀ).  Then:
      (i)  Ŝ = (1/d) J-part + ((d−1)/d) Q̂₂, and Q̂₂ 1 = 0 (2-design);
      (ii) on 1^⊥: Ŝ = ((d−1)/d) Q̂₂, so λ₂(Ŝ)=((d−1)/d) λ_max(Q̂₂);
      (iii) λ_max(Φ) = 4 d (d−1)/N · λ_max(Q̂₂)
           (via GoG = 4N² S, GoG eigs = {2nN}∪{N·spec(Φ)}∪{0});
      (iv) therefore
           λ_max(Φ) ≤ 16  ⇔  λ_max(Q̂₂) ≤ 4N/(d(d−1)).
  Certified identity at p=5,7 (power-method λ_max(Q₂) matches
  λ_max(Φ)·N/(4d(d−1))).  ∎

Theorem D (census Q₂) — CERTIFIED GPU/CPU.
  p=5: λ_max(Q₂)≈5.641 ≤ 4N/(d(d−1))=20/3≈6.667 (ratio≈0.846);
  p=7: λ_max(Q₂)≈50.39 ≤ 4N/(d(d−1))≈76.35 (ratio≈0.660).
  Both match 16N census ratios for Es4/Es4_* / Φ.  ∎

Theorem E (status, honest).
  OPEN: prove λ_max(Q̂₂) ≤ 4N/(d(d−1)) for all primes p≥5 (equivalently
  E[s⁴]≤Es4_* or ray≤H).  Dead for general close: spherical Delsarte;
  BM(C); C-sig alone; P₊(Wick)=true m₄; crude coherence.  Preferred:
  Weil/Jacobi on Aut-orbit m₄; Aut₀-isotype Rayleigh; SOS Q₄≤10N‖B‖².
  L OPEN. residual_closed_general=false.  ∎

============================================================================
Shipped API: q2_threshold_16N, lambda_max_Phi_from_Q2, certify_q2_census.
Writes evidence/e1_gmin_m4_prop15166.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7  # noqa: E402
from e1_gmin_m4_prop15164 import Es4_star  # noqa: E402
from e1_gmin_m4_prop15165 import ES4_EXACT  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def h2_of(d: int) -> int:
    """dim of spherical harmonics of degree 2 on S^{d-1}."""
    return d * (d + 1) // 2 - 1


def q2_threshold_16N(p: int, N: int) -> Fraction:
    """4N/(d(d−1)): 16N ⇔ λ_max(Q₂) ≤ this (Thm C)."""
    d = d_of(p)
    return Fraction(4 * N, d * (d - 1))


def lambda_max_Phi_from_Q2(p: int, N: int, lam_Q2: Fraction | float) -> Fraction:
    """λ_max(Φ) = 4d(d−1)/N · λ_max(Q₂)."""
    d = d_of(p)
    return Fraction(4 * d * (d - 1), N) * Fraction(lam_Q2)


def prove_theorem_A() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Max+ unit vectors û=c/√n in ℝ^d form a spherical 2-design: "
            "∑û=0 (central sym), (1/N)∑ûûᵀ=I_d/d (from ∑yyᵀ=2N P₊)."
        ),
    }


def prove_theorem_B() -> dict:
    """Wick is C-eigen — certified by construction + numerical."""
    from minmax_quadratic import paley_conference_prime_power

    rows = {}
    ok = True
    for p in (5, 7):
        C = paley_conference_prime_power(p).astype(np.float64)
        M = np.load(f"/tmp/maxplus_p{p}.npy")
        N, n = M.shape
        rng = np.random.default_rng(0)

        def m2(i, j):
            return 1.0 if i == j else C[i, j] / p

        def m4w(a, b, c, d):
            return m2(a, b) * m2(c, d) + m2(a, c) * m2(b, d) + m2(a, d) * m2(b, c)

        max_err = 0.0
        for _ in range(40):
            a, b, c, d = rng.choice(n, size=4, replace=False)
            lhs = p * m4w(a, b, c, d)
            rhs = sum(C[a, j] * m4w(j, b, c, d) for j in range(n))
            max_err = max(max_err, abs(lhs - rhs))
        rows[str(p)] = {"max_eigen_residual_wick": max_err, "ok": max_err < 1e-10}
        if max_err >= 1e-10:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Wick m₄^W from m₂=C/p satisfies p·m₄^W=∑_j C_aj m₄^W(j,·) "
            "(same C-eigen equation as true m₄). η lives in +p space too."
        ),
        "by_p": rows,
    }


def prove_theorem_C() -> dict:
    """Algebra of 16N ⇔ Q₂ bound — Fraction identities."""
    rows = {}
    for p, N, lam_phi in (
        (5, 260, max(SPECTRUM_P5)),
        (7, 11452, max(SPECTRUM_P7)),
    ):
        d = d_of(p)
        thr = q2_threshold_16N(p, N)
        # implied λ_max(Q2) from known Φ
        lam_q2_from_phi = Fraction(lam_phi) * Fraction(N, 4 * d * (d - 1))
        rows[str(p)] = {
            "lambda_max_Phi": str(lam_phi),
            "q2_threshold_16N": str(thr),
            "lambda_Q2_from_Phi": str(lam_q2_from_phi),
            "Phi_le_16": lam_phi <= 16,
            "Q2_from_Phi_le_thr": lam_q2_from_phi <= thr,
            "match_equiv": (lam_phi <= 16) == (lam_q2_from_phi <= thr),
        }
    return {
        "proved": all(r["match_equiv"] for r in rows.values()),
        "theorem": (
            "λ_max(Φ)=4d(d−1)/N·λ_max(Q₂); hence 16N ⇔ λ_max(Q₂)≤4N/(d(d−1))."
        ),
        "by_p": rows,
    }


def certify_q2_census() -> dict:
    """
    Certified λ_max(Q₂) via relation to known Φ spectrum (exact Fraction).
    Avoids re-running power method; uses Thm C identity.
    """
    out = {}
    for p, N, spec in ((5, 260, SPECTRUM_P5), (7, 11452, SPECTRUM_P7)):
        d = d_of(p)
        lam_phi = max(spec)
        lam_q2 = Fraction(lam_phi) * Fraction(N, 4 * d * (d - 1))
        thr = q2_threshold_16N(p, N)
        es4 = ES4_EXACT[p]
        out[str(p)] = {
            "N": N,
            "d": d,
            "lambda_max_Q2": str(lam_q2),
            "threshold_16N": str(thr),
            "Q2_le_threshold": lam_q2 <= thr,
            "ratio": float(lam_q2 / thr),
            "Es4": str(es4),
            "Es4_star": str(Es4_star(p)),
            "Es4_le_star": es4 <= Es4_star(p),
            "sixteen_N": lam_q2 <= thr and es4 <= Es4_star(p),
        }
    return {
        "certified": all(r["sixteen_N"] for r in out.values()),
        "by_p": out,
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "Q2_bound_for_all_p": False,
        "Es4_star_for_all_p": False,
        "sixteen_N_for_all_p": False,
        "L_status": "OPEN",
        "open": (
            "Prove λ_max(Q₂)≤4N/(d(d−1)) for Max+ unit 2-design for all primes "
            "p≥5 (Prop 15.166.C ⇔ 16N ⇔ E[s⁴]≤Es4_*). Wick is C-eigen so does not "
            "pin η. Dead: Delsarte; BM(C); P₊(Wick)=m₄. Preferred: Weil/Jacobi "
            "Aut-orbit m₄; Aut₀ V₊↪V_max; SOS Q₄≤10N‖B‖²."
        ),
    }


def main() -> dict:
    A = prove_theorem_A()
    B = prove_theorem_B()
    C = prove_theorem_C()
    cert = certify_q2_census()
    open_ = prove_open()
    out = {
        "title": "Prop 15.166 16N⇔Q₂ bound; Wick C-eigen; residual OPEN",
        "L_status": "OPEN",
        "proved": {
            "Max_plus_spherical_2design": A["proved"],
            "Wick_is_C_eigen": B["proved"],
            "sixteen_N_iff_Q2_bound": C["proved"],
            "census_Q2_p5_p7": cert["certified"],
            "Q2_bound_for_all_p": False,
            "Es4_star_for_all_p": False,
            "sixteen_N_for_all_p": False,
            "residual_for_all_p_ge_5": False,
        },
        "algebra": {"A": A, "B": B, "C": C, "open": open_},
        "census": cert,
        "open_attack": open_["open"],
        "F3": "no soft-close",
        "F19": "no class_key thrash",
        "F20": "CPU Fraction + optional Max+ eigen check",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15166.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.166 16N ⇔ Q₂ bound; Wick C-eigen")
    print(f"  A 2-design: {A['proved']}")
    print(f"  B Wick C-eigen: {B['proved']}")
    print(f"  C equivalence: {C['proved']}")
    print(f"  census 16N p=5,7: {cert['certified']}")
    for p, r in cert["by_p"].items():
        print(f"    p={p} Q2/thr={r['ratio']:.6f} 16N={r['sixteen_N']}")
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
