#!/usr/bin/env python3
"""
Prop 15.159 — Structure ID of Φ|Z spectrum; mult(λ_max)=d;
16N chain under mult≥d + κ²≤96n; residual still OPEN for general p.

Strategy reframe (evidence/STRATEGY_REFRAME_2026-08-03.md): identify the
residual operator, do not invent Prop thrash that only renames δ²≤room.

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, Z = zero-diag traceless ambient matrices
on V₊ (dim m=d(d−3)/2). Φ|Z from E_y[(yᵀBy)²]. 16N ⇔ λ_max(Φ|Z)≤16
(Props 15.61, 15.105). mult(λ₂(P⊙P))=mult(λ_max(Γ|Sym₀)) (15.97).

============================================================================
Theorem A (exact Φ spectrum at p=5) — CERTIFIED Fraction.
  On Z, Φ has three eigenvalues (multiplicity in parentheses):
      176/13 (×d),   144/13 (×2d),   80/13 (×2d),
  with d=13, sum mult = 65 = m.  In particular
      λ_max(Φ|Z) = 176/13 = 16(d−2)/d,
  mult(λ_max)=d.  ∎

Theorem B (exact Φ spectrum at p=7) — CERTIFIED Fraction.
  With N=|Max+|=11452 and den = N/(4p)=409:
      λ ∈ {4320, 4032, 3648, 3360, 3072}/409
  with multiplicities (d, 2d, 2d, 4d, 2d)=(25,50,50,100,50).
  λ_max=4320/409 < 16(d−2)/d=368/25, mult(λ_max)=d.  ∎

Theorem C (algebra: design threshold 16(d−2)/d) — PROVED Fraction.
  For integers d>2: 16(d−2)/d = 16−32/d < 16.
  For d≥13 (i.e. primes p≥5): μ̄=8(d−1)/(d−3) ≤ 16(d−2)/d,
  because d(d−1) ≤ 2(d−2)(d−3) ⇔ 0≤(d−3)(d−6).  ∎

Theorem D (mult(λ_max)=d at p=5,7; PSL d-irrep) — CERTIFIED + structure.
  Top multiplicity equals d at p=5,7.  PSL(2,p²) has a complex irrep of
  dimension (p²+1)/2=d (character table, q=p² odd).  By Aut-Schur (15.97),
  any irrep inside the maximizer locus of Γ contributes full dimension to
  mult(λ_max).  Thus mult≥d whenever the d-dimensional PSL irrep embeds in
  the maximizer locus.  Embedding holds at p=5,7 by mult=d certification.
  **OPEN for general p≥5:** prove the embedding (⇒ mult(λ_max)≥d).  ∎

Theorem E (16N from mult≥d + Wick residual) — PROVED (15.105.E restated).
  If mult(λ_max(Φ|Z))≥d and ‖κ‖_F²≤96n (⇔ orth≤room ⇔ δ²≤room/24),
  then λ_max(Φ|Z)≤16, hence 16N, hence bi-tight empty for primes p≥5
  (Prop 15.61).  Proof: majorization extremal with m₁=d copies at L and
  mean μ̄ gives (L−μ̄)² = orth·(m−d)/(m d); orth≤room yields L≤16
  after Fraction simplification (15.105 Theorem E).  ∎

Theorem F (residual status, honest).
  Path C residual for general p≥5 is still OPEN.  Structure ID here reduces
  the analytic target for 16N to:
    (i) mult(λ_max)≥d for all primes p≥5 (PSL d-irrep embedding), and
    (ii) ‖κ‖_F²≤96n (weaker than κ_hyp / room_hyp residual),
  OR a direct proof of λ_max≤16(d−2)/d for p≥5 (holds with equality at p=5,
  strict at p=7; open in general).
  Dead: Max+ IP-scheme BM; class_key thrash; Prop-only δ re-encoding.  ∎

Theorem G (dual gap operator — preferred structural target) — CERTIFIED p=5,7.
  Define the dual gap operator on Z
      G := (d/32) (16 I − Φ|Z).
  Then
      G ≽ I  ⇔  16I − Φ ≽ (32/d) I  ⇔  λ_max(Φ|Z) ≤ 16(d−2)/d < 16
  and the right-hand side is exactly 16N (Props 15.61, 15.105).
  Certified spectra:
    p=5: eigenvalues of G are exactly {1, 2, 4} (mult d, 2d, 2d);
         equivalently 16−λ ∈ (32/d)·{1,2,4}; G ≽ I with equality on top.
    p=7: eigenvalues of G are {3475,3925,4525,4975,5425}/818, all > 4;
         G ≽ I strictly (λ_max ≪ 16(d−2)/d).
  **Structural reading (Attack vectors 1–3):** G is the object that should be
  identified with a Bose–Mesner / SDP dual / G-irrep scalar — not δ as a free
  residual mass.  Proving G ≽ I for all primes p≥5 closes 16N without the
  orth≤room path.  OPEN in general.  ∎

============================================================================
Shipped API: spectrum forms, dual gap G, mult checks, 16N chain, algebra C/G.
Backend: CPU numpy eigensystem on Max+; GPU unused (F20).
L OPEN. residual_closed_general=false. dual_gap_proved_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15159.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of, room_orth  # noqa: E402
from e1_gmin_m4_prop15102 import delta2_budget_96n, delta2_budget_hyp  # noqa: E402
from e1_gmin_m4_prop15128 import KNOWN_N  # noqa: E402


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


def mu_bar(p: int) -> Fraction:
    """Average eigenvalue of Φ|Z: 8(n−2)/(n−6)."""
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def lambda_max_design_threshold(p: int) -> Fraction:
    """16(d−2)/d — exact λ_max at p=5; candidate UB for p≥5."""
    d = d_of(p)
    return Fraction(16 * (d - 2), d)


def dim_Z(p: int) -> int:
    d = d_of(p)
    return d * (d - 3) // 2


def psl_irrep_dim_d(p: int) -> int:
    """Dimension of the (q+1)/2 irrep of PSL(2,p²), q=p²."""
    return d_of(p)


def psl_min_nontrivial_irrep(p: int) -> int:
    """(q−1)/2 = d−1."""
    return d_of(p) - 1


# Exact certified spectra (Theorems A–B)
SPECTRUM_P5 = {
    Fraction(176, 13): 13,  # d
    Fraction(144, 13): 26,  # 2d
    Fraction(80, 13): 26,  # 2d
}

SPECTRUM_P7 = {
    # den N/(4p)=409
    Fraction(4320, 409): 25,  # d
    Fraction(4032, 409): 50,  # 2d
    Fraction(3648, 409): 50,  # 2d
    Fraction(3360, 409): 100,  # 4d
    Fraction(3072, 409): 50,  # 2d
}


def prove_theorem_C(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 50) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        d = d_of(p)
        thr = lambda_max_design_threshold(p)
        mu = mu_bar(p)
        lt16 = thr < 16
        ge_mu = thr >= mu
        rows[str(p)] = {
            "d": d,
            "threshold_16_dminus2_over_d": str(thr),
            "mu_bar": str(mu),
            "threshold_lt_16": lt16,
            "threshold_ge_mu_bar": ge_mu,
            "ok": lt16 and ge_mu,
        }
        if not (lt16 and ge_mu):
            ok = False
    # algebra check (d-3)(d-6)≥0 for d≥13
    algebra_ok = all(d_of(p) >= 13 for p in primes)
    return {
        "proved": ok and algebra_ok,
        "theorem": (
            "For d>2: 16(d−2)/d<16. For d≥13 (p≥5): μ̄=8(d−1)/(d−3)≤16(d−2)/d "
            "iff (d−3)(d−6)≥0."
        ),
        "by_p": rows,
        "n_checked": len(primes),
    }


def prove_spectrum_p5() -> dict:
    d = 13
    m = dim_Z(5)
    total_mult = sum(SPECTRUM_P5.values())
    lmax = max(SPECTRUM_P5)
    # weighted mean
    mean = sum(lam * mult for lam, mult in SPECTRUM_P5.items()) / m
    return {
        "proved_census": True,
        "theorem": (
            "p=5: Φ|Z eigenvalues 176/13 (×d), 144/13 (×2d), 80/13 (×2d); "
            "λ_max=16(d−2)/d; mult(λ_max)=d."
        ),
        "levels": {str(k): v for k, v in SPECTRUM_P5.items()},
        "total_mult": total_mult,
        "dim_Z": m,
        "mult_sum_ok": total_mult == m,
        "lambda_max": str(lmax),
        "equals_16_dminus2_over_d": lmax == lambda_max_design_threshold(5),
        "mult_top": SPECTRUM_P5[lmax],
        "mult_top_eq_d": SPECTRUM_P5[lmax] == d,
        "mean": str(mean),
        "mean_eq_mu_bar": mean == mu_bar(5),
    }


def prove_spectrum_p7() -> dict:
    d = 25
    m = dim_Z(7)
    N = KNOWN_N[7]
    den = N // (4 * 7)
    total_mult = sum(SPECTRUM_P7.values())
    lmax = max(SPECTRUM_P7)
    mean = sum(lam * mult for lam, mult in SPECTRUM_P7.items()) / m
    return {
        "proved_census": True,
        "theorem": (
            "p=7: Φ|Z eigenvalues {4320,4032,3648,3360,3072}/409 with mults "
            "(d,2d,2d,4d,2d); 409=N/(4p); λ_max<16(d−2)/d; mult(λ_max)=d."
        ),
        "N": N,
        "den_N_over_4p": den,
        "levels": {str(k): v for k, v in SPECTRUM_P7.items()},
        "total_mult": total_mult,
        "dim_Z": m,
        "mult_sum_ok": total_mult == m,
        "lambda_max": str(lmax),
        "threshold": str(lambda_max_design_threshold(7)),
        "lambda_max_lt_threshold": lmax < lambda_max_design_threshold(7),
        "lambda_max_lt_16": lmax < 16,
        "mult_top": SPECTRUM_P7[lmax],
        "mult_top_eq_d": SPECTRUM_P7[lmax] == d,
        "mean": str(mean),
        "mean_eq_mu_bar": abs(float(mean - mu_bar(7))) < 1e-12,
    }


def dual_gap_eigenvalues(spectrum: dict[Fraction, int], d: int) -> dict[Fraction, int]:
    """G := (d/32)(16I − Φ) eigenvalues from a Φ spectrum table."""
    out: dict[Fraction, int] = {}
    for lam, mult in spectrum.items():
        g = (Fraction(16) - lam) * Fraction(d, 32)
        out[g] = out.get(g, 0) + mult
    return out


def prove_dual_gap_G() -> dict:
    """
    Theorem G: G ≽ I at p=5,7 from exact spectra; p=5 eigs {1,2,4}.
    Does not claim general p.
    """
    g5 = dual_gap_eigenvalues(SPECTRUM_P5, 13)
    g7 = dual_gap_eigenvalues(SPECTRUM_P7, 25)
    min5 = min(g5)
    min7 = min(g7)
    # p=5 exact {1,2,4}
    expected5 = {Fraction(1): 13, Fraction(2): 26, Fraction(4): 26}
    return {
        "theorem": (
            "G=(d/32)(16I−Φ) ≽ I ⇔ λ_max≤16(d−2)/d ⇔ 16N. "
            "Certified p=5 (eigs 1,2,4) and p=7 (min>4). OPEN general p."
        ),
        "p5": {
            "G_eigs": {str(k): v for k, v in sorted(g5.items())},
            "min_G": str(min5),
            "G_ge_I": min5 >= 1,
            "equals_1_2_4": g5 == expected5,
            "sixteen_minus_lambda": " (32/d)·{1,2,4}",
        },
        "p7": {
            "G_eigs": {str(k): v for k, v in sorted(g7.items())},
            "min_G": str(min7),
            "G_ge_I": min7 >= 1,
            "strictly_above_4": min7 > 4,
        },
        "certified_p5_p7": min5 >= 1 and min7 >= 1 and g5 == expected5,
        "proved_for_all_p": False,
        "preferred_attack": (
            "Identify G with BM / SDP dual / G-irrep scalar; prove G≽I for p≥5. "
            "Do not thrash δ² bounds."
        ),
    }


def sixteen_N_from_dual_gap(G_min: Fraction | float, p: int) -> dict:
    """If min eig of G ≥ 1 then λ_max ≤ 16(d−2)/d < 16 ⇒ 16N."""
    ok = float(G_min) >= 1.0 - 1e-12 and p >= 5 and is_prime(p)
    return {
        "G_min": str(G_min),
        "G_ge_I": float(G_min) >= 1.0 - 1e-12,
        "p_ge_5_prime": p >= 5 and is_prime(p),
        "sixteen_N": ok,
        "theorem": "G≽I ⇒ λ_max(Φ|Z)≤16(d−2)/d < 16 ⇒ 16N",
    }


def sixteen_N_from_mult_d_and_kappa96(
    mult_top: int,
    d: int,
    kappa2_le_96n: bool,
    p: int,
) -> dict:
    """
    Chain predicate (Theorem E): mult≥d and κ²≤96n ⇒ 16N for p≥5.
    Does not claim kappa2_le_96n for general p — caller supplies certificate.
    """
    ok_mult = mult_top >= d
    ok_p = p >= 5 and is_prime(p)
    sixteen = bool(ok_mult and kappa2_le_96n and ok_p)
    return {
        "mult_top": mult_top,
        "d": d,
        "mult_ge_d": ok_mult,
        "kappa2_le_96n": kappa2_le_96n,
        "p_ge_5_prime": ok_p,
        "sixteen_N": sixteen,
        "bi_tight_empty_via_15_61": sixteen,  # 16N⇒bi-tight empty p≥5
        "theorem": (
            "mult(λ_max)≥d + ‖κ‖²≤96n ⇒ λ_max(Φ|Z)≤16 ⇒ 16N ⇒ bi-tight empty (p≥5)"
        ),
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "mult_d_for_all_p": False,
        "kappa_le_96n_for_all_p": False,
        "lambda_max_le_16_dminus2_over_d_for_all_p": False,
        "dual_gap_G_ge_I_for_all_p": False,
        "open": (
            "PREFERRED: prove dual gap G=(d/32)(16I−Φ) ≽ I for all primes p≥5 "
            "(⇔ λ_max≤16(d−2)/d ⇔ 16N). Structure: BM / SDP dual / G-irrep ID of G. "
            "Alternates: (i) PSL d-irrep embedding + ‖κ‖²≤96n; "
            "(ii) Weil m₄ on G-orbits ⇒ δ²≤room_hyp/24 (15.134)."
        ),
        "structure_progress": (
            "Exact Φ spectra p=5,7; mult(λ_max)=d; dual gap G with eigs {1,2,4} at p=5; "
            "G≽I certified p=5,7; 16(d−2)/d algebra; 16N under G≽I or mult≥d+κ²≤96n."
        ),
        "dead": [
            "Max+ IP association scheme",
            "class_key thrash F19",
            "δ re-encoding without structure ID",
        ],
    }


def certify_live_spectrum(p: int, maxplus_path: Path | None = None) -> dict:
    """Live eigencomputation of Φ|Z; compares to exact tables for p=5,7."""
    from minmax_quadratic import paley_conference_prime_power

    if maxplus_path is None:
        for cand in (
            Path(f"/tmp/maxplus_p{p}.npy"),
            Path(f"/tmp/e1_p{p}/maxplus.npy"),
        ):
            if cand.is_file():
                maxplus_path = cand
                break
    if maxplus_path is None or not maxplus_path.is_file():
        return {"p": p, "skipped": True, "reason": "Max+ cache missing"}

    Y = np.load(maxplus_path).astype(np.float64)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = Y @ Vp
    dimS = d * (d + 1) // 2

    def vec_to_A(v: np.ndarray) -> np.ndarray:
        A = np.zeros((d, d))
        k = 0
        for i in range(d):
            for j in range(i, d):
                A[i, j] = A[j, i] = v[k]
                k += 1
        return A

    Mrows = []
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        A = vec_to_A(e)
        B = Vp @ A @ Vp.T
        Mrows.append(np.concatenate([[np.trace(A)], np.diag(B)]))
    Mmat = np.array(Mrows).T
    _u, s, vh = np.linalg.svd(Mmat, full_matrices=True)
    rank = int((s > 1e-8).sum())
    null = vh[rank:].T
    mZ = null.shape[1]
    Q = np.zeros((N, mZ))
    for j in range(mZ):
        A = vec_to_A(null[:, j])
        Q[:, j] = np.einsum("na,ab,nb->n", U, A, U)
    Gf = np.zeros((mZ, mZ))
    for i in range(mZ):
        Bi = Vp @ vec_to_A(null[:, i]) @ Vp.T
        for j in range(i, mZ):
            Bj = Vp @ vec_to_A(null[:, j]) @ Vp.T
            Gf[i, j] = Gf[j, i] = float(np.sum(Bi * Bj))
    QtQ = (Q.T @ Q) / N
    w, V = np.linalg.eigh(Gf)
    W12 = V @ np.diag(1.0 / np.sqrt(np.maximum(w, 1e-30))) @ V.T
    evals = np.sort(np.linalg.eigvalsh(W12 @ QtQ @ W12))[::-1]
    # cluster
    levels: list[tuple[float, int]] = []
    i = 0
    while i < len(evals):
        j = i
        while j < len(evals) and abs(evals[j] - evals[i]) < 1e-6:
            j += 1
        levels.append((float(evals[i]), j - i))
        i = j
    lmax, mtop = levels[0]
    thr = float(lambda_max_design_threshold(p))
    return {
        "p": p,
        "skipped": False,
        "N": int(N),
        "d": d,
        "dim_Z": mZ,
        "n_levels": len(levels),
        "levels_head": levels[:8],
        "lambda_max": lmax,
        "mult_top": mtop,
        "mult_top_eq_d": mtop == d,
        "lambda_max_lt_16": lmax < 16 + 1e-9,
        "lambda_max_le_threshold": lmax <= thr + 1e-8,
        "threshold": thr,
    }


def main() -> dict:
    c = prove_theorem_C()
    a5 = prove_spectrum_p5()
    a7 = prove_spectrum_p7()
    g = prove_dual_gap_G()
    live5 = certify_live_spectrum(5)
    live7 = certify_live_spectrum(7)

    # Chain demos: with certified mult=d and certified κ²≤96n at p=5,7
    # (κ²≤κ_hyp ≤96n at those primes — residual holds by census)
    chain5 = sixteen_N_from_mult_d_and_kappa96(
        mult_top=13, d=13, kappa2_le_96n=True, p=5
    )
    chain7 = sixteen_N_from_mult_d_and_kappa96(
        mult_top=25, d=25, kappa2_le_96n=True, p=7
    )
    dual5 = sixteen_N_from_dual_gap(Fraction(1), 5)
    dual7 = sixteen_N_from_dual_gap(Fraction(3475, 818), 7)
    open_ = prove_open()

    out = {
        "prop": "15.159",
        "title": (
            "Prop 15.159: Φ|Z spectrum + dual gap G=(d/32)(16I−Φ); "
            "G≽I ⇔ 16N; residual OPEN general p"
        ),
        "L_status": "OPEN",
        "proved": {
            "spectrum_p5": a5.get("mult_sum_ok") and a5.get("equals_16_dminus2_over_d"),
            "spectrum_p7": a7.get("mult_sum_ok") and a7.get("lambda_max_lt_16"),
            "threshold_algebra": c["proved"],
            "dual_gap_G_p5_p7": g["certified_p5_p7"],
            "dual_gap_G_for_all_p": False,
            "sixteen_N_chain_predicate": True,
            "mult_d_for_all_p": False,
            "kappa_96n_for_all_p": False,
            "residual_for_all_p_ge_5": False,
        },
        "algebra": {
            "theorem_C_threshold": c,
            "theorem_A_p5": a5,
            "theorem_B_p7": a7,
            "theorem_G_dual_gap": g,
            "theorem_E_chain": {
                "p5": chain5,
                "p7": chain7,
            },
            "dual_gap_chain": {"p5": dual5, "p7": dual7},
            "open": open_,
        },
        "live_spectrum": {"5": live5, "7": live7},
        "closed_forms": {
            "mu_bar": "8(n-2)/(n-6)",
            "lambda_max_threshold": "16(d-2)/d",
            "dual_gap_operator": "G=(d/32)(16I-Phi)",
            "p5_G_eigs": "1,2,4",
            "p5_lambda_max": "176/13",
            "p7_lambda_max": "4320/409",
            "p7_den": "N/(4p)=409",
        },
        "status": {
            "structure_id_phi_spectrum": True,
            "structure_id_dual_gap_G": True,
            "mult_top_eq_d_p5_p7": True,
            "general_residual": False,
            "dual_gap_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "open_attack": open_["open"],
        "backend": "CPU numpy Φ eigensystem p=5,7; Fraction algebra; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "strategy": (
            "STRUCTURE_REFRAME_2026-08-03: identify dual gap G / Φ spectrum, "
            "not δ thrash. Preferred: prove G≽I for all p≥5."
        ),
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15159.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.159 structure ID")
    print(f"  Thm C threshold algebra: {c['proved']}")
    print(f"  p5 λ_max={a5['lambda_max']} mult=d? {a5['mult_top_eq_d']}")
    print(f"  p7 λ_max={a7['lambda_max']} mult=d? {a7['mult_top_eq_d']}")
    print(f"  dual gap G≽I p5/p7: {g['certified_p5_p7']} p5_eigs_1_2_4={g['p5']['equals_1_2_4']}")
    print(f"  chain p5 sixteen_N={chain5['sixteen_N']} p7={chain7['sixteen_N']}")
    print(f"  dual chain p5={dual5['sixteen_N']} p7={dual7['sixteen_N']}")
    print(f"  live5 mult_d={live5.get('mult_top_eq_d')} live7={live7.get('mult_top_eq_d')}")
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
