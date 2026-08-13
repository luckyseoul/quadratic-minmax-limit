#!/usr/bin/env python3
"""
Prop 15.245 — Z-frame of Max+ in 𝒲₊₊⁰; Op = E[Z⊗Z]; scheme kernel;
average Rayleigh; residual-(i) still OPEN on λ_* floor.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189 S=2P₊, 15.207–213, 15.242–244)
  Max+ = boolean +p-eigenvectors of the conference C (Cy=py, y∈{±1}ⁿ).
  E[yyᵀ] = 2 P₊ = S (spherical 2-design in V₊).  n=p²+1, d=n/2.
  𝒲₊₊⁰ = {B = P₊ B P₊ : diag B = 0}, dim D = n(n−6)/8.
  λ_* = 8(n−6)/n.  Op(B) := E[(yᵀ B y) Z_y] with Rayleigh E[(yᵀBy)²]/‖B‖².

PROVED Max+-free
  A. Zero-diag ++ residual of yyᵀ.
     Set Z_y := yyᵀ − 2 P₊.  Then:
       (i)  Z_y is ++: range in V₊ (y∈V₊ and P₊ projects to V₊);
       (ii) diag(Z_y)_i = y_i² − 2(P₊)_{ii} = 1 − 1 = 0
            ((P₊)_{ii}=1/2 from C_ii=0);
       (iii) hence Z_y ∈ 𝒲₊₊⁰ for every y∈Max+.  ∎

  B. Constant Frobenius norm.
     ‖yyᵀ‖_F² = n², ⟨yyᵀ, 2P₊⟩_F = 2 yᵀ P₊ y = 2n,
     ‖2P₊‖_F² = 4 Tr(P₊) = 4d = 2n.
     So ‖Z_y‖_F² = n² − 4n + 2n = n(n−2).  ∎

  C. Pairing formula.
     ⟨Z_y, Z_z⟩_F = (y·z)² − 2n.  ∎

  D. Quadratic evaluation.
     For B∈𝒲₊₊⁰: ⟨B, 2P₊⟩_F = 2 Tr(B)=0, so
        yᵀ B y = ⟨B, yyᵀ⟩_F = ⟨B, Z_y⟩_F.
     Hence E[(yᵀ B y)²] = E[⟨B, Z_y⟩²] = ⟨B, Op B⟩_F
     with Op = E[Z ⊗ Z] (covariance of Z on 𝒲₊₊⁰).  ∎

  E. Scheme-image is ker Op (on traceless ++).
     For sum-zero f and A = Df C + C Df (scheme): yᵀ A y = 2p ∑ f_i = 0
     on Max+ (15.207), so Op vanishes on the scheme-image (dim n−1).
     Frobenius: scheme-image ⊥ 𝒲₊₊⁰ because ⟨P₊ Df P₊, B⟩ = ∑ f_i B_ii = 0
     for diag B=0.  Dimension check on traceless ++:
        dim{traceless ++} − (n−1) = d(d+1)/2 − 1 − (n−1) = D.
     So traceless ++ = scheme-image ⊕ 𝒲₊₊⁰ (orthogonal), Op|scheme=0,
     and the residual-(i) PD question is Op|𝒲₊₊⁰ ≻ 0.  ∎

  F. Trace and average Rayleigh (Max+-free).
     Tr(Op|𝒲₊₊⁰) = E[‖Z‖_F²] = n(n−2).
     Average eigenvalue = n(n−2)/D = 8(n−2)/(n−6).
     For all n≥26: λ_* = 8(n−6)/n ≤ 8(n−2)/(n−6)
        ⇔ (n−6)² ≤ n(n−2) ⇔ 36 ≤ 10n, true.  ∎
     (Average strictly above λ_*; the floor is a spectral gap question.)

  G. λ_* residual recall (15.243–244).
     E[q²] ≥ λ_* ⇔ ∑_S ρ κ_B ≥ −6/n − 1/(2p²) on unit 𝒲₊₊⁰.  ∎

CERTIFIED
  H. At p=3,5,7: Tr match; spectrum of Op on 𝒲₊₊⁰ has min ≥ λ_*
     (sharp mult n at p=5); average matches 8(n−2)/(n−6).

OPEN (blocks residual i / predicate flip)
  I. Prove Op ≽ λ_* I on 𝒲₊₊⁰ for all primes p≥5
     (equivalently ∑ρ κ_B ≥ −6/n−1/(2p²)), or hull/|μ|/K₄≤Wick_hi.
     Soft-close forbidden.  Note: CS on (ρ,κ_B) is dead (too weak);
     4-design shortcut dead (15.213); unrestricted traceless min~0
     (scheme kernel) so zero-diag is essential.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15245.json
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
)
from e1_gmin_m4_prop15211 import lambda_star  # noqa: E402
from e1_gmin_m4_prop15212 import dim_Wpp0  # noqa: E402
from e1_gmin_m4_prop15243 import lambda_star_rho_budget  # noqa: E402


def Z_frobenius_sq(n: int) -> int:
    """‖Z_y‖_F² = n(n−2)."""
    return n * (n - 2)


def average_rayleigh(n: int) -> Fraction:
    """Tr(Op)/D = 8(n−2)/(n−6)."""
    return Fraction(8 * (n - 2), n - 6)


def theorem_Z_in_Wpp0() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For y∈Max+ (Cy=py, y∈{±1}ⁿ): Z_y=yyᵀ−2P₊ lies in 𝒲₊₊⁰, "
            "‖Z_y‖_F²=n(n−2), and ⟨Z_y,Z_z⟩=(y·z)²−2n.  Max+-free "
            "(uses only Cy=py, C_ii=0, E not required for A–C)."
        ),
        "depends_on": ["Max_plus_boolean_evec", "Pplus_diag_half", "conference_C2"],
    }


def theorem_Op_trace_average() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(3, 80) if p == 3 or is_prime(p)]:
        n = n_of(p)
        D = dim_Wpp0(p) if p >= 5 or p == 3 else n * (n - 6) // 8
        if p == 3:
            D = n * (n - 6) // 8
        Tr = Z_frobenius_sq(n)
        if D <= 0:
            continue
        avg = Fraction(Tr, D)
        if avg != average_rayleigh(n):
            ok = False
        ls = lambda_star(n) if n >= 26 else Fraction(0)
        if n >= 26 and avg < ls:
            ok = False
        if p in (3, 5, 7, 11, 13, 17):
            sample[str(p)] = {
                "n": n,
                "D": D,
                "Tr": Tr,
                "average": str(avg),
                "lambda_star": str(ls) if n >= 26 else "n/a",
                "avg_ge_lstar": bool(n < 26 or avg >= ls),
            }
    return {
        "proved": ok,
        "theorem": (
            "Tr(Op|𝒲₊₊⁰)=n(n−2); average Rayleigh=8(n−2)/(n−6) ≥ λ_*=8(n−6)/n "
            "for all primes p≥5.  Scheme-image is ker Op on traceless ++; "
            "traceless ++ = scheme ⊕ 𝒲₊₊⁰ orthogonally."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_Z_in_Wpp0",
            "prop_15_207_scheme_ker_Q",
            "dim_Wpp0",
        ],
    }


def theorem_lambda_star_vs_average() -> dict:
    ok = True
    for p in [p for p in range(5, 200) if is_prime(p)]:
        n = n_of(p)
        # (n-6)^2 <= n(n-2) <=> 36 <= 10n
        if Fraction(n - 6, 1) ** 2 > Fraction(n * (n - 2), 1):
            ok = False
        if average_rayleigh(n) < lambda_star(n):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "λ_* ≤ average Rayleigh on 𝒲₊₊⁰ for all primes p≥5 "
            "(⇔ 36 ≤ 10n, true for n≥26)."
        ),
        "depends_on": ["theorem_Op_trace_average", "prop_15_211_lambda_star"],
    }


def residual_i_closed_via_245() -> bool:
    return False


def hinge_status_245() -> dict:
    return {
        "Z_frame_proved": theorem_Z_in_Wpp0()["proved"],
        "trace_average_proved": theorem_Op_trace_average()["proved"],
        "lstar_le_average": theorem_lambda_star_vs_average()["proved"],
        "Op_ge_lstar_general": False,
        "residual_i_closed_via_245": residual_i_closed_via_245(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Op ≽ λ_* I on W++0 for all p≥5 (∑ρ κ_B budget), "
            "or hull/|μ|≤1/(2p)/K4≤Wick_hi"
        ),
        "rho_budget_p5": str(lambda_star_rho_budget(5)),
    }


def main() -> dict:
    a = theorem_Z_in_Wpp0()
    b = theorem_Op_trace_average()
    c = theorem_lambda_star_vs_average()
    status = hinge_status_245()
    out = {
        "prop": "15.245",
        "title": (
            "Z-frame in W++0; Op=E[Z⊗Z]; scheme ker; average Rayleigh; "
            "λ_* floor OPEN"
        ),
        "proved": {
            "Z_in_Wpp0": a["proved"],
            "Op_trace_average": b["proved"],
            "lstar_le_average": c["proved"],
            "Op_ge_lstar_general": False,
            "residual_i_closed": residual_i_closed_via_245(),
        },
        "algebra": {
            "Z_norm_sq": "n(n-2)",
            "Z_pairing": "(y·z)^2 - 2n",
            "average_rayleigh": "8(n-2)/(n-6)",
            "lambda_star": "8(n-6)/n",
            "avg_p5": str(average_rayleigh(26)),
            "lstar_p5": str(lambda_star(26)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_245(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Z-frame + scheme kernel reduce residual-(i) PD to a spectral "
            "gap Op≽λ_* on W++0.  Average is strictly above λ_*.  Gap OPEN. "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15245.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.245  residual_i={residual_i_closed_via_245()}")
    print(f"  Z_frame={a['proved']} trace_avg={b['proved']} lstar_le_avg={c['proved']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
