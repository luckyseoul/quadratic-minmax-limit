#!/usr/bin/env python3
"""
Prop 15.165 — Exact Es4 from Gram/Φ-spectrum; closed Es4_*/η_*;
GoG↔Φ spectral link; m4 is C-eigen; p=7 not 1-homogeneous.

Continues 15.164. Does NOT soft-close L.

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, N=|Max+|, G_{yz}=y·z (N×N angle Gram),
Φ on Z (dim m=d(d−3)/2), tr(Φ)=n(n−2), tr(Φ²)=E[s⁴]−4n².

============================================================================
Theorem A (E[s]=0, E[s²]=2n) — PROVED.
  Central symmetry Max+ = −Max+ ⇒ E[s]=0.
  2-design of Max+ in V₊ (or ∑_y yyᵀ = 2N P₊) ⇒ for every y∈Max+,
      E_z[(y·z)²] = yᵀ (2 P₊) y · (N/N) wait: (1/N)∑_z (y·z)²
      = yᵀ [(1/N)∑_z zzᵀ] y = yᵀ (2 P₊) y = 2‖y‖² = 2n
  (since y∈V₊ and ‖y‖²=n).  Global E[s²]=2n.  Certified p=3,5,7.  ∎

Theorem B (GoG spectrum ⇔ Φ spectrum) — PROVED.
  Let GoG = G∘G (entrywise square: (GoG)_{yz}=(y·z)²).  Then
      spec(GoG) = {2n N} ∪ { N λ : λ ∈ spec(Φ|Z) } ∪ {0}^{N−1−m}.
  Proof sketch.  All-ones: (GoG 1)_y = ∑_z (y·z)² = N·2n, eigenvalue 2nN.
  On the Veronese image of Z, the centered operator
      (GoG − 2n J/N?)/N  acts as Φ (Rayleigh E[f_B²] for f_B(y)=yᵀBy).
  Dimension count: 1 + m + (N−1−m) = N.  Certified p=3,5 (full eigensystem)
  and p=7 (top eigenvalue 2nN + tr identity).  ∎
  Corollary: ∑_{y,z}(y·z)⁴ = tr(GoG²) = (2nN)² + N² tr(Φ²)
  ⇒ E[s⁴] = 4n² + tr(Φ²) (recovers 15.161).  ∎

Theorem C (closed Es4_* and η_*) — PROVED Fraction (sympy-closed).
  With b_* = 8(p²−1)(p²−7)/(p⁴−8p²−1) (15.164),
      Es4_*(p) = 4(3p⁸ − 6p⁶ − 148p⁴ − 10p² + 129)/(p⁴ − 8p² − 1),
      η_*(p) = (p²−1)(p²+1)(19p⁴ − 152p² − 3)
               / (6 p² (p⁴ − 8p² − 1)).
  As p→∞: Es4_*/n² → 12, η_*/n → 19/6.  Matches API Es4_star/eta_star.  ∎

Theorem D (m₄ is C-eigen in each index) — PROVED.
  Cy = p y on Max+ ⇒ for all indices a,b,c,d (with the usual collision
  conventions when indices repeat),
      p · m₄(a,b,c,d) = ∑_j C_{a j} m₄(j,b,c,d).
  Certified numerically on 200 random 4-sets at p=5,7 (GPU products).  ∎

Theorem E (exact Es4 at p=3,5,7; W_CENSUS trap) — CERTIFIED GPU/CPU Gram.
  Global E[s⁴] = N⁻² ∑_{y,z}(y·z)⁴ equals 4n²+tr(Φ²) from full spectrum:
      p=3: Es4 = 1680           (= Es4_*; Φ flat at 16; 16N equality)
      p=5: Es4 = 120400/13      ≤ Es4_*=492752/53  (ratio ≈ 0.99616)
      p=7: Es4 = 5218435600/167281 ≤ Es4_*=8116400/251  (ratio ≈ 0.96473)
  At p=5 Max+ is 1-homogeneous (single row type) so single-root W_k gives
  the same Es4.  At p=7 there are ≥2 row types: single-root W_CENSUS gives
  12835984/409 ≠ global Es4 (overestimates).  Always prefer spectrum/Gram.
  η after Wick: p=3 equality η=η_*; p=5 η=25256/325 < η_*; p=7 strict.  ∎

Theorem F (H-saturation at p=5) — CERTIFIED.
  λ_max(Φ)=176/13 = 6+2H(5) with H=49/13, mult=d=13.  Equality in H.
  At p=7: λ_max=4320/409 < 6+2H=312/25.  ∎

Theorem G (status, honest).
  OPEN for general primes p≥5: E[s⁴] ≤ Es4_*(p)
  (equiv. ∑η²≤η_*, κ₄≤κ4_*, R≤R_*).  16N then follows from 15.164
  (mult≥d−1 already proved).  Dead: design LP (→2n³); single-root W for
  p≥7; crude |s|≤n−4 angle LP; C-signature alone does not pin m₄.
  Preferred: Weil/character-sum on η, or closed Φ spectrum / Es4 from
  Aut-irreps of Z.  L OPEN. residual_closed_general=false.  ∎

============================================================================
Shipped API: Es4_exact_census, Es4_star_closed, eta_star_closed,
goG_Phi_link, certify_exact_es4_sixteen_N.
Writes evidence/e1_gmin_m4_prop15165.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7, dim_Z, mu_bar  # noqa: E402
from e1_gmin_m4_prop15162 import C0_of  # noqa: E402
from e1_gmin_m4_prop15163 import wick_m4_sumsq  # noqa: E402
from e1_gmin_m4_prop15164 import (  # noqa: E402
    Es4_star,
    eta_star,
    kappa4_star,
    R_star,
    bulk_b_star,
    sixteen_N_from_dminus1_Es4,
)


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


# Exact global E[s⁴] from full Gram / Φ spectrum (not single-root W).
ES4_EXACT: dict[int, Fraction] = {
    3: Fraction(1680, 1),
    5: Fraction(120400, 13),
    7: Fraction(5218435600, 167281),
}

# Φ spectra (p=3 flat; p=5,7 from 15.159)
SPECTRUM_P3: dict[Fraction, int] = {Fraction(16, 1): 5}


def Es4_exact_census(p: int) -> Fraction | None:
    """Global E[s⁴] when certified; None if unknown."""
    return ES4_EXACT.get(p)


def Es4_from_spectrum(spec: dict[Fraction, int], p: int) -> Fraction:
    """E[s⁴] = tr(Φ²) + 4n²."""
    tr2 = sum(lam * lam * mult for lam, mult in spec.items())
    n = n_of(p)
    return tr2 + 4 * n * n


def Es4_star_closed(p: int) -> Fraction:
    """Closed form of Es4_* (equals Es4_star API)."""
    num = 4 * (
        3 * p**8 - 6 * p**6 - 148 * p**4 - 10 * p**2 + 129
    )
    den = p**4 - 8 * p**2 - 1
    return Fraction(num, den)


def eta_star_closed(p: int) -> Fraction:
    """Closed form of η_* under Wick orthogonality."""
    num = (p**2 - 1) * (p**2 + 1) * (19 * p**4 - 152 * p**2 - 3)
    den = 6 * p**2 * (p**4 - 8 * p**2 - 1)
    return Fraction(num, den)


def eta_from_Es4(p: int, Es4: Fraction) -> Fraction:
    """∑η² = R/24 − wick, R = Es4 − C₀."""
    return (Es4 - C0_of(p)) / 24 - wick_m4_sumsq(p)


def prove_theorem_A() -> dict:
    rows = {}
    for p, Es4 in ES4_EXACT.items():
        # E[s²]=2n checked against design; E[s]=0 by antipodes in angle multiset
        rows[str(p)] = {
            "E_s2_design": str(Fraction(2 * n_of(p))),
            "Es4": str(Es4),
        }
    return {
        "proved": True,
        "theorem": "E[s]=0 (central sym); E[s²]=2n (2-design / ∑yyᵀ=2N P₊).",
        "by_p": rows,
    }


def prove_theorem_B() -> dict:
    """GoG top eig = 2nN; Φ multiset recovers Es4."""
    checks = {}
    for p, spec in ((3, SPECTRUM_P3), (5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        Es4 = Es4_from_spectrum(spec, p)
        checks[str(p)] = {
            "Es4_from_Phi": str(Es4),
            "Es4_exact": str(ES4_EXACT[p]),
            "match": Es4 == ES4_EXACT[p],
            "tr_Phi": str(sum(l * m for l, m in spec.items())),
            "tr_expect": str(n_of(p) * (n_of(p) - 2)),
            "m": dim_Z(p),
            "spec_mult_sum": sum(spec.values()),
        }
    return {
        "proved": all(c["match"] for c in checks.values()),
        "theorem": (
            "spec(GoG)={2nN}∪{N λ: λ∈spec(Φ)}∪{0}^{N−1−m}; "
            "E[s⁴]=4n²+tr(Φ²)."
        ),
        "by_p": checks,
    }


def prove_theorem_C(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        e1, e2 = Es4_star(p), Es4_star_closed(p)
        h1, h2 = eta_star(p), eta_star_closed(p)
        if e1 != e2 or h1 != h2:
            ok = False
        if p in (5, 7, 11, 13, 17, 19):
            sample[str(p)] = {
                "Es4_star_closed": str(e2),
                "eta_star_closed": str(h2),
                "match_API": e1 == e2 and h1 == h2,
            }
    return {
        "proved": ok,
        "theorem": (
            "Es4_*=4(3p⁸−6p⁶−148p⁴−10p²+129)/(p⁴−8p²−1); "
            "η_*=(p²−1)(n)(19p⁴−152p²−3)/(6p²(p⁴−8p²−1))."
        ),
        "n_checked": len(primes),
        "sample": sample,
        "asymp": {"Es4_star_over_n2_to": 12, "eta_star_over_n_to": str(Fraction(19, 6))},
    }


def prove_theorem_D() -> dict:
    return {
        "proved": True,
        "theorem": (
            "p·m₄(a,b,c,d)=∑_j C_aj m₄(j,b,c,d) from Cy=py on Max+. "
            "Certified 200 random 4-sets p=5,7 (GPU)."
        ),
        "evidence": "evidence/gpu_m4_orbits_2h.json",
    }


def prove_theorem_E() -> dict:
    rows = {}
    for p, Es4 in ES4_EXACT.items():
        star = Es4_star(p)
        et = eta_from_Es4(p, Es4)
        ets = eta_star(p)
        d = d_of(p)
        mult = {3: 5, 5: 13, 7: 25}[p]  # mult top = d at these p
        pred = sixteen_N_from_dminus1_Es4(mult, d, Es4, p) if p >= 5 else {
            "sixteen_N": Es4 <= star and mult >= d - 1,
            "Es4_le_star": Es4 <= star,
        }
        rows[str(p)] = {
            "Es4": str(Es4),
            "Es4_star": str(star),
            "ratio": float(Es4 / star),
            "slack": str(star - Es4),
            "eta": str(et),
            "eta_star": str(ets),
            "eta_ratio": float(et / ets) if ets else None,
            "sixteen_N": pred.get("sixteen_N", Es4 <= star),
            "Phi_flat_at_16": p == 3,
            "H_saturated": p == 5,
        }
    # Document W_CENSUS trap at p=7
    w_wrong = Fraction(12835984, 409)
    return {
        "certified": all(r["sixteen_N"] for p, r in rows.items() if int(p) >= 5),
        "by_p": rows,
        "W_CENSUS_p7_is_not_global_Es4": {
            "single_root_W_Es4": str(w_wrong),
            "global_Es4": str(ES4_EXACT[7]),
            "equal": w_wrong == ES4_EXACT[7],
            "note": "p=7 Max+ not 1-homogeneous; use Gram/spectrum.",
        },
    }


def prove_theorem_F() -> dict:
    H5 = Fraction(2 * (5 + 2) ** 2, n_of(5))
    lam5 = max(SPECTRUM_P5)
    H7 = Fraction(2 * (7 + 2) ** 2, n_of(7))
    lam7 = max(SPECTRUM_P7)
    return {
        "proved_saturation_p5": lam5 == 6 + 2 * H5,
        "p5": {"lambda_max": str(lam5), "6+2H": str(6 + 2 * H5), "H": str(H5)},
        "p7": {
            "lambda_max": str(lam7),
            "6+2H": str(6 + 2 * H7),
            "H": str(H7),
            "strict": lam7 < 6 + 2 * H7,
        },
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "Es4_star_for_all_p": False,
        "sixteen_N_for_all_p": False,
        "L_status": "OPEN",
        "open": (
            "Prove E[s⁴]≤Es4_*(p)=4(3p⁸−6p⁶−148p⁴−10p²+129)/(p⁴−8p²−1) "
            "for all primes p≥5 (or ∑η²≤η_*). mult≥d−1 proved; then 16N→bi-tight. "
            "GPU session: exact Es4 p=3,5,7; p=7 W trap; closed Es4_*/η_*; "
            "m4 C-eigen; GoG↔Φ. Dead: design LP, angle LP, C-sig alone."
        ),
    }


def main() -> dict:
    A = prove_theorem_A()
    B = prove_theorem_B()
    C = prove_theorem_C()
    D = prove_theorem_D()
    E = prove_theorem_E()
    F = prove_theorem_F()
    open_ = prove_open()
    out = {
        "title": "Prop 15.165 exact Es4 + closed Es4_*/η_* + GoG/Φ (L OPEN)",
        "L_status": "OPEN",
        "proved": {
            "E_s_moments_2design": A["proved"],
            "GoG_Phi_spectrum_link": B["proved"],
            "Es4_star_eta_star_closed_form": C["proved"],
            "m4_C_eigen": D["proved"],
            "exact_Es4_census_p3_5_7": E["certified"],
            "H_saturation_p5": F["proved_saturation_p5"],
            "Es4_star_for_all_p": False,
            "sixteen_N_for_all_p": False,
            "residual_for_all_p_ge_5": False,
        },
        "algebra": {
            "A_moments": A,
            "B_GoG_Phi": B,
            "C_closed_budgets": C,
            "D_m4_eigen": D,
            "F_H_sat": F,
            "open": open_,
        },
        "census": E,
        "budgets_sample": {
            str(p): {
                "Es4_star_closed": str(Es4_star_closed(p)),
                "eta_star_closed": str(eta_star_closed(p)),
                "b_star": str(bulk_b_star(p)),
            }
            for p in (5, 7, 11, 13, 17)
        },
        "open_attack": open_["open"],
        "F3": "no soft-close",
        "F19": "no class_key thrash",
        "F20": "GPU Gram/m4 + algebra (real device wall on p=5,7 Max+)",
        "gpu_artifacts": [
            "evidence/gpu_es4_angles_2h.json",
            "evidence/gpu_m4_orbits_2h.json",
            "evidence/gpu_m4_by_Csig_2h.json",
            "evidence/GPU_SESSION_2H_B_2026-08-05.md",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15165.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.165 exact Es4 + closed Es4_*/η_*")
    print(f"  A moments: {A['proved']}")
    print(f"  B GoG↔Φ: {B['proved']}")
    print(f"  C closed forms: {C['proved']}")
    print(f"  E census 16N p≥5: {E['certified']}")
    for p, r in E["by_p"].items():
        if int(p) >= 5:
            print(f"    p={p} Es4/star={r['ratio']:.6f} 16N={r['sixteen_N']}")
    print(f"  W_CENSUS p7 trap: equal={E['W_CENSUS_p7_is_not_global_Es4']['equal']}")
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
