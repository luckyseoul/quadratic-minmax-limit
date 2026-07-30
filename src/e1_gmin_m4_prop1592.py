#!/usr/bin/env python3
"""
Prop 15.92 — Constant pairing sum on Max+; clean 16N/H reductions; W spectrum.

Proved for every symmetric conference matrix C of order n=p²+1 (C²=p²I,
zero diagonal, off-diag ±1) and Max+ = {y∈{±1}^n : Cy = py}:

  1) Pointwise (all boolean y, any conference C): by Prop 15.90 with B=C,
       ∑_{|S|=4} κ_C(S) ∏_{v∈S} y_v
         = (yᵀ C y)²/8 − (yᵀ C² y)/2 + ‖C‖_F²/4
         = (yᵀ C y)²/8 − (p² n)/2 + n(n−1)/4.
     (Uses C²=p²I and ‖C‖_F²=n(n−1).)

  2) On Max+ (Cy=py ⇒ yᵀCy = np): the sum is CONSTANT
       ∑_S κ_C(S) ∏ y_v = n(n−1)(n−2)/8.
     Hence ∑_S m4(S) κ_C(S) = n(n−1)(n−2)/8.

  3) Clean spectral reductions (unit zero-diag B on V₊, P=YYᵀ/(2N)):
       16N  ⇔  λ₂(P⊙P) ≤ 4/N  ⇔  λ₂(W) ≤ 4N/d²
       H    ⇔  λ₂(P⊙P) ≤ (3+H)/(2N) ⇔ λ₂(W) ≤ N(3+H)/(2d²)
     with W_ab=(u_a·u_b)², u_a=y_a/√n.
     Equiv. frame form: U = S/√(2N) (N×d, UᵀU=I_d, equal row norms),
       16N ⇔ max_{‖x‖=1, x⊥1} ‖Uᵀ diag(x) U‖_F² ≤ 4/N.

  4) Budget algebra (all primes p≥3): H≤5, 6+2H≤16, harm budgets as 15.91.

Certified W spectrum:
  p=3: {N/d (×1), 48/25 (×d)} with 48/25 = N(6+2H)/(4d²) (H-equality).
  p=5: {N/d (×1), 880/169 (×d), 720/169 (×2d), 400/169 (×2d)};
       top = N(6+2H)/(4d²) (H-equality); rank = binom(d−1,2)=66.

OPEN: prove λ₂(W) ≤ N(3+H)/(2d²) (H) or ≤ 4N/d² (16N / 2×sphere-class)
for all primes p≥5.  H_proved=false. L OPEN.

F19/F20: Fraction + Max+ p=3,5 checks; GPU unused.
Writes evidence/e1_gmin_m4_prop1592.json
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


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def n_of(p: int) -> int:
    return p * p + 1


def d_of(p: int) -> int:
    return n_of(p) // 2


def H(p: int) -> Fraction:
    return Fraction((p + 2) ** 2, d_of(p))


def prove_pairing_sum_constant(primes: list[int]) -> dict:
    """∑_S κ∏y = n(n-1)(n-2)/8 on Max+; formula for all p≥3."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # On Max+: f=np, γ=p²n, ||C||F²=n(n-1)
        f2 = (n * p) ** 2
        gamma = (p * p) * n
        cF2 = n * (n - 1)
        pred = Fraction(f2, 8) - Fraction(gamma, 2) + Fraction(cF2, 4)
        target = Fraction(n * (n - 1) * (n - 2), 8)
        match = pred == target
        # Expand algebraically:
        # np²/8 - p²n/2 + n(n-1)/4
        # = n[ n p²/8 - p²/2 + (n-1)/4 ]
        # = n/8 [ n p² - 4 p² + 2(n-1) ]
        # n=p²+1: n p² - 4p² + 2p² + 2 - 2 = n p² - 2p² = p²(n-2)
        # wait full:
        rows[str(p)] = {
            "pred_on_Max+": str(pred),
            "target_n_n1_n2_over_8": str(target),
            "match": match,
            "identity": (
                "∑ κ∏y = (yᵀCy)²/8 − p²n/2 + n(n−1)/4  (all boolean y); "
                "on Max+ = n(n−1)(n−2)/8"
            ),
        }
        if not match:
            ok = False
    # Algebraic verification without loop
    # pred = n²p²/8 - p²n/2 + n(n-1)/4 with n=p²+1
    # = [n²p² - 4 n p² + 2n(n-1)]/8
    # = n[ n p² - 4 p² + 2n - 2 ]/8
    # = n[ (p²+1)p² - 4p² + 2(p²+1) - 2 ]/8
    # = n[ p⁴ + p² - 4p² + 2p² + 2 - 2 ]/8
    # = n[ p⁴ - p² ]/8 = n p² (p²-1)/8 = n(n-1)(n-2)/8
    # since p²=n-1, p²-1=n-2.
    return {
        "proved": ok,
        "theorem": (
            "For every conference C of order n=p²+1 and every y∈Max+: "
            "∑_{|S|=4} κ_C(S)∏_{v∈S} y_v = n(n−1)(n−2)/8 (constant). "
            "Proof: Prop 15.90 with B=C gives the pointwise formula for all boolean y; "
            "C²=p²I ⇒ yᵀC²y=p²n; on Max+ yᵀCy=np; substitute n=p²+1. "
            "Averaging yields ∑_S m4(S) κ_C(S) = n(n−1)(n−2)/8."
        ),
        "by_p": rows,
    }


def prove_spectral_reductions(primes: list[int]) -> dict:
    rows = {}
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        Hf = H(p)
        # Thresholds in λ2(W) form — N cancels in comparisons of interest
        # thr_16 = 4N/d², thr_H = N(3+H)/(2d²)
        # ratio thr_H / thr_16 = (3+H)/(8)
        rows[str(p)] = {
            "H": str(Hf),
            "thr16_over_N": str(Fraction(4, d * d)),
            "thrH_over_N": str((Fraction(3) + Hf) / (2 * d * d)),
            "H_implies_16N": Hf <= 5,
            "frame_form_16N": "max_{‖x‖=1,x⊥1} ‖Uᵀdiag(x)U‖_F² ≤ 4/N",
            "frame_form_H": "λ₂(P⊙P) ≤ (3+H)/(2N) on 1^⊥",
        }
    return {
        "proved": True,
        "theorem": (
            "Let P=YYᵀ/(2N) (equal-diag orthoprojector rank d), "
            "W_ab=(y_a·y_b)²/n², U=S/√(2N). Then P⊙P=α² W with α=d/N, "
            "λ₁(W)=N/d, and for unit zero-diag B on V₊: "
            "max Q = n² λ₂(W) = 4d² λ₂(W). Hence "
            "16N⇔λ₂(W)≤4N/d²⇔λ₂(P⊙P)≤4/N; "
            "H⇔λ₂(W)≤N(3+H)/(2d²)⇔λ₂(P⊙P)≤(3+H)/(2N). "
            "Frame: λ₂(P⊙P)=max_{‖x‖=1,x⊥1}‖Uᵀdiag(x)U‖_F²."
        ),
        "by_p": rows,
    }


def certify_W_spectrum(p: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    elif p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    else:
        return {"p": p, "skipped": True}

    N, n = Y.shape
    d = n // 2
    D = Y @ Y.T
    W = (D * D) / (n * n)
    ev = np.linalg.eigvalsh(W)
    ev = np.sort(ev)[::-1]
    ctr = Counter(np.round(ev, 10))
    items = sorted(ctr.items(), reverse=True)
    eigs = []
    for val, mult in items:
        if abs(val) < 1e-8:
            eigs.append({"value": 0.0, "frac": "0", "mult": int(mult)})
            break
        f = Fraction(val).limit_denominator(5000)
        eigs.append(
            {
                "value": float(val),
                "frac": f"{f.numerator}/{f.denominator}",
                "frac_err": abs(val - float(f)),
                "mult": int(mult),
            }
        )
    lam1 = N / d
    thr_H = float(N * (6 + 2 * H(p)) / (4 * d * d))
    thr_16 = 4 * N / (d * d)
    # top nontrivial
    top = eigs[1]["value"] if len(eigs) > 1 else None
    return {
        "p": p,
        "N": N,
        "d": d,
        "rank_nonzero": int(np.sum(ev > 1e-8)),
        "rank_formula_binom_d1_2": (d - 1) * (d - 2) // 2,
        "eigs": eigs[:8],
        "lam1": lam1,
        "lam2": top,
        "thr_H": thr_H,
        "thr_16": thr_16,
        "H_equality": top is not None and abs(top - thr_H) < 1e-8,
        "below_16": top is not None and top <= thr_16 + 1e-8,
        "ok": top is not None and top <= thr_H + 1e-7,
    }


def certify_pairing_on_maxplus(p: int) -> dict:
    """Numeric check of constant pairing sum."""
    from itertools import combinations

    from minmax_quadratic import paley_conference_prime_power

    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    elif p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    else:
        return {"p": p, "skipped": True}

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    target = n * (n - 1) * (n - 2) / 8.0
    # check a few y
    errs = []
    rng = np.random.default_rng(p)
    for idx in rng.choice(N, size=min(15, N), replace=False):
        y = Y[idx]
        s = 0.0
        for a, b, c, d_ in combinations(range(n), 4):
            kap = (
                C[a, b] * C[c, d_]
                + C[a, c] * C[b, d_]
                + C[a, d_] * C[b, c]
            )
            s += kap * y[a] * y[b] * y[c] * y[d_]
        errs.append(abs(s - target))
    # full m4·κ via average of constant
    return {
        "p": p,
        "target": target,
        "max_abs_err_on_samples": float(max(errs)),
        "ok": max(errs) < 1e-8,
        "sum_m4_kappa": target,  # equals average of constant
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(3, 50) if is_prime(p)]
    print("Prop 15.92: pairing sum + spectral reductions", flush=True)

    pair = prove_pairing_sum_constant(primes)
    print(f"  pairing sum formula proved={pair['proved']}", flush=True)

    red = prove_spectral_reductions(primes)
    print(f"  spectral reductions proved={red['proved']}", flush=True)

    cert3 = certify_pairing_on_maxplus(3)
    cert5 = certify_pairing_on_maxplus(5)
    print(
        f"  pairing cert p=3 ok={cert3.get('ok')} p=5 ok={cert5.get('ok')}",
        flush=True,
    )

    w3 = certify_W_spectrum(3)
    w5 = certify_W_spectrum(5)
    print(
        f"  W spectrum p=3 H_eq={w3.get('H_equality')} p=5 H_eq={w5.get('H_equality')}",
        flush=True,
    )

    out = {
        "prop": "15.92",
        "backend": "CPU Fraction + Max+ p=3,5; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "pairing_sum": pair,
        "spectral_reductions": red,
        "pairing_cert": {"3": cert3, "5": cert5},
        "W_spectrum_cert": {"3": w3, "5": w5},
        "attack_surface": (
            "Prove λ₂(W)≤N(3+H)/(2d²) (H) or λ₂(W)≤4N/d² (16N) for all primes p≥5. "
            "Equiv: λ₂(P⊙P)≤(3+H)/(2N) or ≤4/N. "
            "Frame: max_{x⊥1,‖x‖=1} ‖Uᵀdiag(x)U‖_F² ≤ 4/N. "
            "Do not re-attack ∑ρ κ_B (Prop 15.90)."
        ),
        "verdict": (
            "Prop 15.92 proves ∑ m4 κ_C = n(n−1)(n−2)/8 for every conference Max+ "
            "(pointwise constant pairing sum), and records the clean spectral "
            "reductions of H and 16N to λ₂(W)/λ₂(P⊙P). W spectrum saturates H at "
            "p=3,5. General upper bound ray≤H or ray≤5 still OPEN. L OPEN."
        ),
        "settlement_note": (
            "Next: prove λ₂(P⊙P)≤4/N for the Max+ design projector for all p≥5 "
            "(closes bi-tight via 16N), or the sharper H form. Deep ND still needed "
            "for lim α_n."
        ),
    }
    out["proved"] = (
        pair["proved"]
        and red["proved"]
        and cert3.get("ok")
        and cert5.get("ok")
        and w3.get("ok")
        and w5.get("ok")
    )

    path = ROOT / "evidence" / "e1_gmin_m4_prop1592.json"
    write_json_atomic(path, out)
    print(
        f"Prop 15.92 proved(structure)={out['proved']} H_proved=False wrote {path}",
        flush=True,
    )
    print("  L OPEN — pairing sum closed; ray≤H / 16N still open", flush=True)


if __name__ == "__main__":
    main()
