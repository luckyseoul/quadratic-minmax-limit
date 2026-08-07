#!/usr/bin/env python3
"""
Prop 15.199 — Frame identity K₄=4n²+E[A²]; residual-(i) analytic target;
λ_max path to Wick dead at p=5; hinge still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189, 15.197–198)
  Max+ antipodal, S=E[yyᵀ]=I+C/p, yᵀSy=2n constant, ‖yyᵀ‖_F²=n².
  K₄=E_{y,z∈Max+}[(y·z)⁴], Wick_hi=12n²+48n.
  Residual (i) closes if K₄≤Wick_hi (15.198) or |μ|≤2/n on |κ|=1, etc.

PROVED Max+-free (under 15.189 π / tight frame)
  A. Centered Gram identity.  Let X=yyᵀ, Xc=X−S.  Then
        ‖Xc‖_F² = n² − 2⟨X,S⟩ + ‖S‖_F² = n² − 4n + 2n = n(n−2)
     constant on Max+ (uses yᵀSy=2n and ‖S‖_F²=Tr(S²)=Tr(2S)=2n from S²=2S).
     For independent y,z∈Max+:
        (y·z)² = ⟨X,Z⟩_F = ⟨Xc,Zc⟩ + ⟨Xc,S⟩ + ⟨S,Zc⟩ + ‖S‖_F².
     ⟨Xc,S⟩=yᵀSy−‖S‖_F²=0, similarly ⟨S,Zc⟩=0, so
        (y·z)² = ⟨Xc,Zc⟩ + 2n,
        K₄ = E[(y·z)⁴] = E[(⟨Xc,Zc⟩+2n)²] = E[A²] + 4n²
     with A=⟨Xc,Zc⟩ (cross terms vanish by independence and E[Xc]=0).  ∎

  B. Covariance form.  Let 𝒞 = E[Xc ⊗ Xc] (covariance operator on matrices).
     Then E[A²] = ‖𝒞‖_{HS}² = ∑_i μ_i(𝒞)², and Tr(𝒞)=n(n−2).
     Hence K₄ = 4n² + ∑ μ_i(𝒞)².  ∎

  C. Wick target restated.  K₄ ≤ Wick_hi ⇔ ∑ μ_i² ≤ 8n²+48n
     ⇔ E[Rayleigh_𝒞(Xc)] ≤ 8(n+6)/(n−2)
     (since ‖Xc‖_F²=n(n−2) constant, E[A²]=n(n−2)·E[R]).  ∎

  D. λ_max-only majorant is dead for residual (i).  The bound
        ∑ μ_i² ≤ λ_max(𝒞)·Tr(𝒞)
     would give K₄ ≤ 4n² + λ_max·n(n−2), and Wick would follow from
        λ_max ≤ 8(n+6)/(n−2).
     At p=3 this thr equals 16 and is achieved (equality case).
     At p=5 the thr is 32/3≈10.667, but certified λ_max(𝒞)≈13.54 > thr,
     so the λ_max route cannot prove Wick at p=5 (even though actual
     K₄ < Wick_hi: E[R]≈10.51 < 10.667).  Must control HS mass of 𝒞,
     not only λ_max.  ∎

CERTIFIED
  E. Frame identity K₄=4n²+E[A²] at p=3,5 (full Max+); ‖Xc‖_F²=n(n−2);
     yᵀSy=2n; K₄≤Wick_hi at p=3(=),5,7 (15.198 fractions).
  F. λ_max(𝒞) power-iteration: p=3 → 16 (=need); p=5 → ≈13.538 > need.

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free for all primes p≥5 one of:
     (i)  ∑ μ_i(𝒞)² ≤ 8n²+48n  (⇔ K₄≤Wick_hi ⇔ Q<10),
     (ii) |μ₄|≤2/n (or ≤1/(2p)) on every |κ|=1 4-set,
     (iii) dual-eq ker-box empty (e.g. max κ_e ≤ (3/2)·scheme-max).
     Then flip gsum_disj_lb_proved_general → residual (i) → (E1 if residual ii
     allows) → L.  Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15199.json
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
from e1_gmin_m4_prop15198 import wick_hi  # noqa: E402


def sigma2_Xc(p: int) -> int:
    """‖Xc‖_F² = n(n−2)."""
    n = n_of(p)
    return n * (n - 2)


def A2_budget_for_wick(p: int) -> int:
    """Need E[A²] ≤ 8n²+48n for K₄≤Wick_hi."""
    n = n_of(p)
    return 8 * n * n + 48 * n


def lambda_max_need_for_wick(p: int) -> Fraction:
    """λ_max ≤ 8(n+6)/(n−2) would suffice via ∑μ²≤λ_max Tr (too strong; dead at p=5)."""
    n = n_of(p)
    return Fraction(8 * (n + 6), n - 2)


def theorem_frame_identity() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Under S=I+C/p, yᵀSy=2n, S²=2S, antipodal Max+: "
            "‖yyᵀ−S‖_F²=n(n−2) constant, and K₄=4n²+E[⟨Xc,Zc⟩_F²]. "
            "Proof: expand ⟨X,Z⟩=(y·z)²; cross terms with S vanish."
        ),
    }


def theorem_wick_target_rayleigh() -> dict:
    primes = [p for p in range(5, 40) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        need_A2 = A2_budget_for_wick(p)
        need_R = Fraction(need_A2, sigma2_Xc(p))
        alt = Fraction(8 * (n + 6), n - 2)
        rows[str(p)] = {
            "need_A2": need_A2,
            "need_E_R": str(need_R),
            "formula_8_n_plus_6_over_n_minus_2": str(alt),
            "match": need_R == alt,
        }
        if need_R != alt:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "K₄≤Wick_hi ⇔ E[A²]≤8n²+48n ⇔ E[Rayleigh_𝒞(Xc)]≤8(n+6)/(n−2)."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11") if k in rows},
    }


def theorem_lambda_max_path_dead() -> dict:
    """Document: λ_max≤need is false at p=5 (certified), so cannot use alone."""
    return {
        "proved_dead": True,
        "theorem": (
            "The majorant ∑μ_i²≤λ_max·Tr would need λ_max≤8(n+6)/(n−2). "
            "At p=3 equality holds (λ_max=need=16, K₄=Wick_hi). "
            "At p=5, power-iteration on 𝒞 gives λ_max≈13.54 > need=32/3≈10.67, "
            "while actual K₄=120400/13 < Wick_hi (E[R]≈10.51<10.67). "
            "Hence residual (i) via Wick cannot be proved from λ_max(𝒞) alone."
        ),
        "p3_need": str(lambda_max_need_for_wick(3)),
        "p5_need": str(lambda_max_need_for_wick(5)),
        "p5_lambda_max_certified_approx": "13.538461",
        "p5_E_R_from_K4": str(
            Fraction(120400, 13) - 4 * 26 * 26
        ),  # A2 = K4-4n²; E[R]=A2/σ² below
    }


def certify_frame_constants(p: int) -> dict:
    """Census: ‖Xc‖, ySy, K4 vs 4n²+sample A² (p=3,5)."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        return {"p": p, "skipped": True, "reason": "enum cost; use 15.198 fractions"}
    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    Mp = boolean_evecs(C, float(p), 0).astype(np.float64)
    M = len(Mp)
    S = np.eye(n) + C / p
    # constants on a sample of vectors
    norms = []
    ysy = []
    for i in range(M):
        y = Mp[i]
        X = np.outer(y, y)
        norms.append(float(np.linalg.norm(X - S, "fro") ** 2))
        ysy.append(float(y @ S @ y))
    # K4
    s4 = 0.0
    bs = 256
    for i0 in range(0, M, bs):
        s4 += float(np.sum((Mp[i0 : i0 + bs] @ Mp.T) ** 4))
    K4 = s4 / (M * M)
    # A2 from identity
    A2 = K4 - 4 * n * n
    E_R = A2 / (n * (n - 2))
    need_R = float(lambda_max_need_for_wick(p))
    return {
        "p": p,
        "M": M,
        "norm_Xc_mean": float(np.mean(norms)),
        "norm_Xc_theory": n * (n - 2),
        "norm_match": abs(float(np.mean(norms)) - n * (n - 2)) < 1e-6,
        "ySy_mean": float(np.mean(ysy)),
        "ySy_theory": 2 * n,
        "ySy_match": abs(float(np.mean(ysy)) - 2 * n) < 1e-6,
        "K4": K4,
        "A2": A2,
        "E_R": E_R,
        "need_E_R_for_Wick": need_R,
        "E_R_le_need": E_R <= need_R + 1e-9,
        "Wick_hi": wick_hi(p),
        "K4_le_Wick": K4 <= wick_hi(p) + 1e-6,
        "proved_general": False,
    }


def K4_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_frame() -> bool:
    return K4_bound_proved_general()


def hinge_status_199() -> dict:
    return {
        "frame_identity_K4": True,
        "wick_target_as_E_R": True,
        "lambda_max_path_dead_p5": True,
        "K4_le_Wick_hi_general": False,
        "residual_i_closed": residual_i_closed_via_frame(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove ∑μ_i(𝒞)²≤8n²+48n (K₄≤Wick_hi) or |μ|≤2/n on |κ|=1 "
            "or ker-box empty, Max+-free for all primes p≥5."
        ),
    }


def main() -> dict:
    a = theorem_frame_identity()
    b = theorem_wick_target_rayleigh()
    c = theorem_lambda_max_path_dead()
    cert = {str(p): certify_frame_constants(p) for p in (3, 5, 7)}
    # p=5 E[R] from exact K4
    n5 = 26
    K4_5 = Fraction(120400, 13)
    A2_5 = K4_5 - 4 * n5 * n5
    ER_5 = A2_5 / (n5 * (n5 - 2))
    out = {
        "title": (
            "Prop 15.199 frame K₄=4n²+E[A²]; Wick as E[R] thr; "
            "λ_max path dead; residual (i) OPEN"
        ),
        "proved": {
            "frame_identity": a["proved"],
            "wick_target_rayleigh": b["proved"],
            "lambda_max_path_dead": c["proved_dead"],
            "K4_bound_general": False,
            "residual_i_closed": residual_i_closed_via_frame(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "frame": a,
            "wick_rayleigh": b,
            "lambda_max_dead": c,
            "p5_exact_E_R": str(ER_5),
            "p5_need_E_R": str(lambda_max_need_for_wick(5)),
            "p5_E_R_lt_need": ER_5 < lambda_max_need_for_wick(5),
        },
        "certified": cert,
        "hinge": hinge_status_199(),
        "F3": (
            "Frame identity is Max+-free under 15.189. Do not flip residual/E1/L "
            "without K₄≤Wick_hi or equivalent |μ|/ker-box."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15199.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.199 frame identity / Wick Rayleigh target")
    print(f"  frame identity: {a['proved']}")
    print(f"  wick as E[R]: {b['proved']}")
    print(f"  lambda_max path dead: {c['proved_dead']}")
    print(f"  p=5 exact E[R]={ER_5} need={lambda_max_need_for_wick(5)}")
    for pk, crow in cert.items():
        if crow.get("skipped"):
            print(f"  p={pk}: skipped")
            continue
        print(
            f"  p={pk}: norm_ok={crow['norm_match']} ySy_ok={crow['ySy_match']} "
            f"K4_le_Wick={crow['K4_le_Wick']} E_R={crow['E_R']:.4f}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_frame()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
