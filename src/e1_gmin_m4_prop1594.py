#!/usr/bin/env python3
"""
Prop 15.94 — P⊙P annihilates range(P); gap criterion via mult≥d.

Proved (all conference Max+, central symmetry):

  1) P⊙P vanishes on range(P).
     Proof: range(P)=colspace(Y). For z=Yv, z_a=y_a·v.
     (P⊙P z)_a = ∑_b P_ab² z_b = (1/(4N²)) ∑_b (y_a·y_b)² (y_b·v).
     Coordinate i: ∑_b (y_a·y_b)² y_{b i}
       = ∑_{jk} y_{a j} y_{a k} ∑_b y_{b j} y_{b k} y_{b i}
       = N ∑_{jk} y_{a j} y_{a k} E[y_j y_k y_i].
     Third moments E[y_j y_k y_i]=0 by central symmetry (Max+ = −Max+).
     Hence P⊙P z=0.

  2) Corollary: ∑_{a,b} P_ab³ = Tr((P⊙P)P)=0.

  3) Spectral support: λ₁(P⊙P)=α=d/N on span(1); range(P)⊂ker(P⊙P);
     remaining positive spectrum lives in (span(1)⊕range(P))⊥ and has
     rank(P⊙P)−1 = binom(d−1,2)−1 = d(d−3)/2 eigenvalues (when rank formula holds)
     summing to S=d(d−1)/N.

  4) Gap criterion (proved algebra): if mult(λ₂(P⊙P))≥d and
     Q:=Tr((P⊙P)²)−α² satisfies Q ≤ d³/(4N²), then
     λ₂ ≤ √(Q/d) ≤ d/(2N) (the spectral gap of Prop 15.56),
     hence bi-tight empty for all primes p≥5 (Prop 15.55).
     Equiv. Q-bound: ∑_{a,b} P_ab⁴ ≤ d²/N² + d³/(4N²) = d²(d+4)/(4N²)
     i.e. sum_M² := ∑_{ijkl} M_ijkl² ≤ 4 d²(d+4) with M_ijkl=E[y_i y_j y_k y_l].

  5) Identity: sum_M² = 12n² − 48n + ∑ κ²
     (Wick algebra for Σ=2P₊, Σ²=2Σ; sum_M_Wick=12n²; sum_Wick²=12n²+48n).

  6) Certified: mult(λ₂)=d at p=3,5; gap criterion holds at p=5
     (√(Q/d)≈0.0216 ≤ 0.025=d/(2N)); fails at p=3 as required
     (√(Q/d)=1/3 > 5/24=d/(2N)). 16N criterion √(Q/d)≤4/N fails at p=5
     (0.0216 > 0.0154) — need mult or sharper tools for 16N.

OPEN: prove mult(λ₂)≥d for all primes p≥5, and sum_M² ≤ 4d²(d+4)
(or λ_max(FFT|1⊥)≤8N / H). L OPEN.

F19/F20: algebra + Max+ p=3,5; GPU unused.
Writes evidence/e1_gmin_m4_prop1594.json
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


def prove_annihilator() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For every centrally symmetric Max+ of a conference matrix, "
            "the Schur square P⊙P of the design projector P=YYᵀ/(2N) "
            "annihilates range(P). Proof: range(P)=colspace(Y); for z=Yv one has "
            "(P⊙P z)_a ∝ ∑_b (y_a·y_b)² (y_b·v); the coefficient of each ambient "
            "coordinate is N ∑_{jk} y_{aj} y_{ak} E[y_j y_k y_i]=0 because all "
            "third moments vanish by Max+=−Max+."
        ),
        "corollary_sum_P3": "∑_{a,b} P_ab³ = Tr((P⊙P)P) = 0",
    }


def prove_gap_criterion(primes: list[int]) -> dict:
    """Algebra: mult≥d and Q≤d³/(4N²) ⇒ gap. N cancels in the sum_M² form."""
    rows = {}
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        # sum_M² bound for gap: 4 d² (d+4)
        thr_M2 = 4 * d * d * (d + 4)
        # Wick base of sum_M²
        wick_base = 12 * n * n - 48 * n  # 12n²-48n
        # room for sum_κ²
        room_kappa2 = thr_M2 - wick_base
        rows[str(p)] = {
            "d": d,
            "n": n,
            "thr_sum_M2": thr_M2,
            "wick_base_12n2_minus_48n": wick_base,
            "room_for_sum_kappa2": room_kappa2,
            "room_positive": room_kappa2 > 0,
            "gap_thr_form": "d/(2N)",
            "sixteen_N_thr_form": "4/N",
            "note": (
                "gap via mult≥d: λ₂≤√(Q/d); need Q≤d³/(4N²) ⇔ sum_M²≤4d²(d+4)"
            ),
        }
    return {
        "proved": True,
        "theorem": (
            "If mult(λ₂(P⊙P))≥d and ∑_{ijkl}M_ijkl² ≤ 4d²(d+4), then "
            "λ₂(P⊙P)≤d/(2N) (spectral gap), hence for every prime p≥5 "
            "bi-tight covers are empty (Prop 15.55–15.56). "
            "Proof: PSD rest eigs with mult(λ₂)≥d ⇒ λ₂≤√(Q/d); "
            "Q=Tr((P⊙P)²)−α²=sum_M²/(16N²)−d²/N²; "
            "√(Q/d)≤d/(2N) rearranges to sum_M²≤4d²(d+4). "
            "Identity sum_M²=12n²−48n+∑κ² from Wick algebra on Σ=2P₊."
        ),
        "by_p": rows,
    }


def C_diag_formula(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(4 * n * (11 * n - 14), p * p)


def certify_annihilator_and_gap(p: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    elif p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
        if not path.is_file():
            return {"p": p, "skipped": True, "reason": "no /tmp/e1_p7/maxplus.npy"}
        Y = np.load(path).astype(np.float64)
    else:
        return {"p": p, "skipped": True}

    N, n = Y.shape
    d = n // 2
    C = paley_conference_prime_power(p).astype(np.float64)
    P = (Y @ Y.T) / (2.0 * N)
    PopP = P * P
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    S = Y @ Vp
    U = S / np.sqrt(2.0 * N)  # ONB of range(P)

    # Annihilator check
    ann_err = float(np.linalg.norm(PopP @ U))
    # sum P^3
    sum_P3 = float(np.sum(P**3))

    # Spectrum: full eig for N small; iterative for p=7
    alpha = d / N
    if N <= 400:
        ev = np.sort(np.linalg.eigvalsh(PopP))[::-1]
        pos = ev[ev > 1e-10]
        ctr = Counter(np.round(pos, 10))
        items = sorted(ctr.items(), reverse=True)
        if len(items) >= 2:
            lam2_val, lam2_mult = items[1]
        else:
            lam2_val, lam2_mult = items[0][0], 0
        tr2 = float(np.sum(PopP * PopP))
        rank_pop = int(np.sum(ev > 1e-10))
    else:
        # Iterative: top eigs via subspace iteration + matvec
        scale = 1.0 / (4 * N * N)

        def matvec(z):
            Wmat = Y.T @ (z[:, None] * Y)
            return scale * np.sum((Y @ Wmat) * Y, axis=1)

        rng = np.random.default_rng(p)
        k = min(40, N - 2)
        X = rng.standard_normal((N, k))
        X -= X.mean(axis=0, keepdims=True)
        q, _ = np.linalg.qr(X)
        evals = np.zeros(k)
        for _it in range(50):
            Ym = np.column_stack([matvec(q[:, j]) for j in range(k)])
            Ym -= Ym.mean(axis=0, keepdims=True)
            G = q.T @ Ym
            evals, vecs = np.linalg.eigh(0.5 * (G + G.T))
            idx = np.argsort(evals)[::-1]
            evals = evals[idx]
            q = np.linalg.qr(Ym @ vecs[:, idx])[0]
        lam2_val = float(evals[0])
        lam2_mult = int(np.sum(evals >= lam2_val - 1e-7))
        # s4 for Tr2
        s4 = 0.0
        bs = 400
        for i0 in range(0, N, bs):
            s4 += float(np.sum((Y[i0 : i0 + bs] @ Y.T) ** 4))
        tr2 = s4 / (16.0 * N**4)
        rank_pop = -1  # not computed
        pos = evals

    # Q and criterion
    Q = tr2 - alpha**2
    sqrt_Q_over_d = float(np.sqrt(max(Q, 0.0) / d))
    gap_thr = d / (2.0 * N)
    thr16 = 4.0 / N
    Hf = (p + 2) ** 2 / d
    thrH = (3 + Hf) / (2.0 * N)

    # sum_M2 from s4
    if N <= 400:
        D = Y @ Y.T
        s4 = float(np.sum(D**4))
    sum_M2 = s4 / (N * N)
    thr_M2 = 4.0 * d * d * (d + 4)
    wick_base = 12 * n * n - 48 * n
    sum_kappa2_est = sum_M2 - wick_base
    C_diag = float(C_diag_formula(p))
    sum_rho2_est = (sum_kappa2_est - C_diag) / 24.0

    return {
        "p": p,
        "N": int(N),
        "d": int(d),
        "annihilator_err": ann_err,
        "annihilator_ok": ann_err < 1e-8,
        "sum_P3": sum_P3,
        "sum_P3_ok": abs(sum_P3) < 1e-8,
        "lam1_or_top": float(pos[0]) if N <= 400 else float(evals[0]),
        "alpha": alpha,
        "lam2": float(lam2_val),
        "lam2_mult": int(lam2_mult),
        "mult_ge_d": int(lam2_mult) >= d,
        "Q": Q,
        "sqrt_Q_over_d": sqrt_Q_over_d,
        "gap_thr": gap_thr,
        "thr16": thr16,
        "thrH": thrH,
        "gap_by_mult_d": bool(int(lam2_mult) >= d and sqrt_Q_over_d <= gap_thr + 1e-9),
        "sixteen_N_by_mult_d": bool(
            int(lam2_mult) >= d and sqrt_Q_over_d <= thr16 + 1e-9
        ),
        "actual_gap_holds": float(lam2_val) <= gap_thr + 1e-9,
        "actual_16N_holds": float(lam2_val) <= thr16 + 1e-9,
        "sum_M2": sum_M2,
        "thr_M2_gap": thr_M2,
        "sum_M2_le_thr": sum_M2 <= thr_M2 + 1e-4,
        "wick_base": wick_base,
        "sum_kappa2_est": sum_kappa2_est,
        "C_diag_formula": C_diag,
        "sum_rho2_est": sum_rho2_est,
        "rank_PopP": rank_pop,
        "rank_formula": (d - 1) * (d - 2) // 2,
        "ok": bool(ann_err < 1e-8 and abs(sum_P3) < 1e-8),
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(3, 40) if is_prime(p)]
    print("Prop 15.94: P⊙P annihilates range(P); gap via mult≥d", flush=True)

    ann = prove_annihilator()
    print(f"  annihilator proved={ann['proved']}", flush=True)

    gapc = prove_gap_criterion(primes)
    print(f"  gap criterion algebra proved={gapc['proved']}", flush=True)

    cert3 = certify_annihilator_and_gap(3)
    cert5 = certify_annihilator_and_gap(5)
    print(
        f"  p=3: ann={cert3['annihilator_ok']} mult={cert3['lam2_mult']} "
        f"gap_by_mult={cert3['gap_by_mult_d']} (expect False)",
        flush=True,
    )
    print(
        f"  p=5: ann={cert5['annihilator_ok']} mult={cert5['lam2_mult']}≥d={cert5['d']} "
        f"gap_by_mult={cert5['gap_by_mult_d']} 16N_by_mult={cert5['sixteen_N_by_mult_d']}",
        flush=True,
    )
    print(
        f"    sqrt(Q/d)={cert5['sqrt_Q_over_d']:.8f} gap_thr={cert5['gap_thr']:.8f} "
        f"16N_thr={cert5['thr16']:.8f}",
        flush=True,
    )

    out = {
        "prop": "15.94",
        "backend": "CPU algebra + Max+ p=3,5; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "gap_proved_for_all_p": False,
        "annihilator": ann,
        "gap_criterion": gapc,
        "cert": {"3": cert3, "5": cert5},
        "attack_surface": (
            "To close bi-tight for all p≥5, prove ONE of:\n"
            "  (A) mult(λ₂(P⊙P))≥d and ∑M²≤4d²(d+4)  [gap via two-moment];\n"
            "  (B) λ_max(FFT|1⊥)≤8N  [16N];\n"
            "  (C) λ_max(FFT|1⊥)≤N(3+H)  [H];\n"
            "  (D) max|m4|≤(p−2)/(2p²) on |κ|=1  [m4 path, Prop 15.66].\n"
            "New structure: P⊙P kills range(P) (central symmetry). "
            "At p=5, (A) holds with certified mult=d and exact Q. "
            "Need mult≥d for general p, and control of ∑M² (or ∑κ²)."
        ),
        "verdict": (
            "Prop 15.94 proves P⊙P annihilates range(P) by vanishing third moments "
            "(central symmetry), and records the gap criterion mult(λ₂)≥d + "
            "∑M²≤4d²(d+4). Certified at p=5: mult=d and criterion holds for gap "
            "(not for 16N). General mult≥d and ∑M² bound OPEN. "
            "16N/H OPEN. lim α_n OPEN."
        ),
        "settlement_note": (
            "Next: prove mult(λ₂(P⊙P))≥d for all primes p≥5 (representation / "
            "explicit d-dimensional test space), and ∑M_ijkl²≤4d²(d+4). "
            "Then gap⇒bi-tight empty. Or prove λ_max(FFT|1⊥)≤8N. Deep ND remains."
        ),
    }
    out["proved"] = (
        ann["proved"]
        and gapc["proved"]
        and cert3.get("ok")
        and cert5.get("ok")
        and cert5.get("gap_by_mult_d") is True
        and cert3.get("gap_by_mult_d") is False  # p=3 must fail gap
    )

    path = ROOT / "evidence" / "e1_gmin_m4_prop1594.json"
    write_json_atomic(path, out)
    print(
        f"Prop 15.94 proved(structure)={out['proved']} gap_all_p=False wrote {path}",
        flush=True,
    )
    print(
        "  L OPEN — annihilator proved; gap at p=5 via mult=d; general mult OPEN",
        flush=True,
    )


if __name__ == "__main__":
    main()
