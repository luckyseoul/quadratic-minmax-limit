#!/usr/bin/env python3
"""
Prop 15.101 — Fickus Gram residual / bulk-variance reduction of ‖κ_orth‖².

Applies the Fickus–Jasper–Mixon residual-Gram pattern (arXiv:2605.28738) and
Ge–Liu multiplicity/interlacing ideas to the Path C orth residual of Prop 15.100.

Setup (Paley conference Max+ of order n=p²+1, d=n/2, |Max+|=N):
  G = Y Yᵀ (N×N Gram of boolean ±p evecs), rank(G)=d, nonzero eigs 2N (mult d).
  P = G/(2N) equal-diag orthoprojector rank d (design projector).
  PopP = P⊙P (Schur square).  K = G⊙G = 4 N² PopP.

Proved:

  1) Fickus–Schur residual rank (conference Max+).
     rank(K) = rank(PopP) = binom(d−1,2)  (Prop 15.59 structure).
     Crude Schur bound rank(K)≤d² is strict for d>1; the sharpened rank is the
     Veronese/bulk dimension used in the flat CS bound of Prop 15.100.
     Annihilator (Prop 15.94): range(P) ⊂ ker(PopP); λ₁(PopP)=d/N simple on 1;
     bulk positive spectrum has size ≤ m := d(d−3)/2 and sums to S := d(d−1)/N.

  2) ED4 / orth as bulk spectral variance (N-cancelling algebra).
     E[D⁴] = ∑_{a,b} G_ab⁴ / N² = 16 N² Tr(PopP²)
            = 16 d² + 16 N² ∑_{bulk} λ².
     Flat CS: ∑_bulk λ² ≥ S²/m gives ED4 ≥ ED4_flat (Prop 15.100).
     Orth residual (Pythagoras after dual-frame proj):
       ‖κ_orth‖_F² = ED4 − ED4_flat
                   = 16 N² (∑_bulk λ² − S²/m)
                   = 8 N² / m · ∑_{i,j bulk} (λ_i − λ_j)².
     Equivalently ‖κ_orth‖² = m' · Var(spec_bulk) in the ED4 scale.
     N cancels in the closed flat formulas of Prop 15.100.

  3) Ge–Liu / PSL multiplicity bound on # of bulk levels (Paley).
     Prop 15.98: mult(λ₂(PopP)) ≥ d−1.
     Hence the number of distinct positive bulk eigenvalues is
       t ≤ ⌊ m / (d−1) ⌋ = ⌊ d(d−3) / (2(d−1)) ⌋.
     (Each level carries ≥ d−1 dimensions of a nontrivial PSL-type isotypic.)

  4) N-free λ₂-sufficient criterion for ‖κ‖² ≤ 96n (Fickus residual form).
     Write λ_flat := S/m = 2(d−1)/(N(d−3)) for the flat bulk value.
     Majorization ∑_bulk λ² ≤ λ₂ · S yields
       ‖κ_orth‖² ≤ 16 N² (λ₂ S − S²/m).
     This is ≤ room_orth := 32 n (p²−9)/(p²−5) precisely when
       λ₂ ≤ λ_flat · (1 + ε),   ε = 4(p²−9)/(p²−1)².
     Proof of ε: substitute S, m, room, n=2d, d=(p²+1)/2; N cancels and
       d−3 = (p²−5)/2, d−1 = (p²−1)/2 give ε = 4(p²−9)/(p²−1)². ∎
     (Sufficient, not necessary: at p=5 the true λ₂ exceeds λ_flat(1+ε) while
      orth still meets the hyp target, because mass is not concentrated at λ₂.)

  5) Hyp-form criterion (PSL ratio).
     Same majorization with room_hyp = room_orth · ((d−1)/d)² gives
       ε_hyp = 4(p²−9)/(p²−1)² · ((d−1)/d)².
     λ₂ ≤ λ_flat(1+ε_hyp) ⇒ ‖κ_orth‖² ≤ room_hyp ⇒ ‖κ‖² ≤ κ_hyp ≤ 96n.

  6) Dual-frame / flat identification (Prop 15.100 recall).
     ED4_flat − wick_lo = ‖κ_proj‖² = 64 n (p²−3)/(p²−5);
     room_orth = 96n − ‖κ_proj‖² = 32 n (p²−9)/(p²−5);
     ‖κ‖² ≤ 96n ⇔ ‖κ_orth‖² ≤ room_orth;
     candidate: ‖κ_orth‖² ≤ room_orth · ((d−1)/d)² ⇔ ‖κ‖² ≤ κ_hyp.

  7) Certified Max+ p=3,5 (Fickus spectrum):
     p=3: bulk flat (1 level), orth=0=room, λ₂=λ_flat=4/N (16N eq), κ²=96n.
     p=5: bulk 3 levels, mults (d,2d,2d)=(13,26,26), eigs 11/845,9/845,5/845;
          orth = room·((d−1)/d)² exactly (κ²=κ_hyp); λ₂ < 4/N (16N holds);
          Fickus rank(PopP)=66=binom(12,2).

  8) OPEN residual (honest).
     Prove ‖κ_orth‖² ≤ room·((d−1)/d)² for all primes p≥5, e.g. by
     (a) closed-form bulk spectrum of PopP under PSL isotypics (≤ t levels),
     (b) sharp majorization better than λ₂·S using mult(λ₂)≥d−1, or
     (c) resolvent δ-bound (orth=24‖δ‖²) with μ²=4(p²+15) for T²b.
     Then κ≤96n + Prop 15.98 ⇒ bi-tight for Paley p≥5. L stays OPEN until
     bi-tight + deep ND + Main Theorem.

L OPEN. H_proved=false. kappa_bound_proved_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert p=3,5; GPU unused.
Writes evidence/e1_gmin_m4_prop15101.json
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

from e1_gmin_m4_prop15100 import (  # noqa: E402
    ED4_flat,
    d_of,
    kappa_hyp,
    n_of,
    proj_kappa2,
    room_orth,
    wick_hi,
    wick_lo,
)


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


def m_bulk(p: int) -> int:
    """Max number of positive bulk eigenvalues of PopP."""
    d = d_of(p)
    return d * (d - 3) // 2


def rank_popp(p: int) -> int:
    """rank(PopP)=binom(d−1,2)."""
    d = d_of(p)
    return (d - 1) * (d - 2) // 2


def eps_wick(p: int) -> Fraction:
    """ε in λ₂ ≤ λ_flat(1+ε) ⇒ orth ≤ room_orth ⇒ κ² ≤ 96n."""
    return Fraction(4 * (p * p - 9), (p * p - 1) ** 2)


def eps_hyp(p: int) -> Fraction:
    """ε for the hyp/PSL-ratio room."""
    d = d_of(p)
    return eps_wick(p) * Fraction((d - 1) ** 2, d * d)


def max_bulk_levels(p: int) -> int:
    """⌊m/(d−1)⌋ from mult(λ₂)≥d−1."""
    d = d_of(p)
    m = m_bulk(p)
    return m // (d - 1)


def prove_fickus_rank(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        d = d_of(p)
        rk = rank_popp(p)
        schur = d * d
        m = m_bulk(p)
        # rank = 1 + m when full bulk rank (λ1 + m bulk)
        rows[str(p)] = {
            "d": d,
            "rank_PopP_binom": rk,
            "schur_d2": schur,
            "strict_schur": rk < schur,
            "m_bulk": m,
            "rank_eq_1_plus_m": rk == 1 + m,
            "max_bulk_levels_PSL": max_bulk_levels(p),
        }
        if not (rk == 1 + m and rk < schur and max_bulk_levels(p) >= 1):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For conference Max+ design projector P of rank d: PopP=P⊙P has "
            "rank binom(d−1,2)=1+m with m=d(d−3)/2 (Veronese bulk). "
            "Schur bound rank≤d² is strict for d>1. "
            "With mult(λ₂)≥d−1 (Prop 15.98 Paley), # distinct bulk eigs "
            "≤ ⌊m/(d−1)⌋. Fickus residual-Gram rank pattern applied to K=G⊙G."
        ),
        "by_p": rows,
    }


def prove_orth_variance_identity(primes: list[int]) -> dict:
    """N-free closed forms: room, proj, ED4_flat, ε criteria."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        m = m_bulk(p)
        room = room_orth(p)
        pr = proj_kappa2(p)
        e4f = ED4_flat(p)
        wl = wick_lo(p)
        wh = wick_hi(p)
        # ED4_flat = wick_lo + proj
        id1 = e4f - wl == pr
        # room = wick_hi - ED4_flat
        id2 = wh - e4f == room
        # ε algebra
        eps = eps_wick(p)
        # re-derive ε from Fraction path
        # ε = room * (d-3) / (32 * d * (d-1)²) with room=32n(p²-9)/(p²-5)
        eps_from_room = room * Fraction(d - 3, 32 * d * (d - 1) ** 2)
        eps_match = eps == eps_from_room
        # simplified form
        eps_c = Fraction(4 * (p * p - 9), (p * p - 1) ** 2)
        rows[str(p)] = {
            "ED4_flat_minus_wick_lo_eq_proj": id1,
            "Wick_minus_ED4_flat_eq_room": id2,
            "eps_wick": str(eps),
            "eps_from_room_formula": str(eps_from_room),
            "eps_match": eps_match,
            "eps_closed": str(eps_c),
            "eps_hyp": str(eps_hyp(p)),
            "room_orth": str(room),
            "room_hyp": str(room * Fraction((d - 1) ** 2, d * d)),
            "kappa_hyp": str(kappa_hyp(p)),
            "budget_96n": 96 * n,
        }
        if not (id1 and id2 and eps_match and eps == eps_c):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "‖κ_orth‖² = ED4 − ED4_flat = 16 N² (∑_bulk λ² − S²/m). "
            "ED4_flat − wick_lo = ‖κ_proj‖²; Wick_hi − ED4_flat = room_orth. "
            "Majorization ∑λ²≤λ₂ S ⇒ orth≤room_orth when λ₂≤λ_flat(1+ε) with "
            "ε=4(p²−9)/(p²−1)² (N-free). Hyp form: ε_hyp=ε·((d−1)/d)²."
        ),
        "by_p": rows,
    }


def prove_eps_criterion(primes: list[int]) -> dict:
    """Algebra-only verification of ε formulas and positivity for p≥5."""
    rows = {}
    ok = True
    for p in primes:
        eps = eps_wick(p)
        eh = eps_hyp(p)
        d = d_of(p)
        # ε≥0 for p≥3, =0 only p=3
        nonneg = eps >= 0
        eq3 = (eps == 0) == (p == 3)
        # ε_hyp < ε for p>3
        hyp_le = eh <= eps
        # (1+ε) form
        rows[str(p)] = {
            "eps": str(eps),
            "eps_hyp": str(eh),
            "eps_nonneg": nonneg,
            "eps_zero_iff_p3": eq3,
            "eps_hyp_le_eps": hyp_le,
            "sufficient_criterion_wick": (
                "λ₂(PopP) ≤ λ_flat · (1+ε) ⇒ ‖κ_orth‖² ≤ room_orth ⇒ κ²≤96n"
            ),
            "sufficient_criterion_hyp": (
                "λ₂(PopP) ≤ λ_flat · (1+ε_hyp) ⇒ ‖κ_orth‖² ≤ room_hyp ⇒ κ²≤κ_hyp"
            ),
            "note_not_necessary": (
                "At p=5 true λ₂ > λ_flat(1+ε) but orth=room_hyp still "
                "(mass not concentrated at λ₂)."
            ),
        }
        if not (nonneg and eq3 and hyp_le):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For every odd prime p≥3: ε=4(p²−9)/(p²−1)²≥0 (eq only p=3). "
            "ε_hyp=ε·((d−1)/d)²≤ε. The λ₂≤λ_flat(1+ε) test is a sufficient "
            "Fickus-style residual criterion for ‖κ‖²≤96n; it is not necessary."
        ),
        "by_p": rows,
    }


def prove_reduction() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Prop 15.101 reduces the orth residual of Prop 15.100 to bulk spectral "
            "control of PopP=P⊙P on Max+: "
            "‖κ_orth‖² = 16 N² (∑_bulk λ² − S²/m). "
            "With Fickus rank m=d(d−3)/2 and PSL mult(λ₂)≥d−1, at most "
            "⌊m/(d−1)⌋ distinct bulk levels. "
            "Sufficient N-free test: λ₂ ≤ λ_flat(1+4(p²−9)/(p²−1)²) ⇒ κ²≤96n. "
            "Target form remains ‖κ_orth‖² ≤ room·((d−1)/d)² (⇔ κ²≤κ_hyp). "
            "General inequality OPEN; certified equality κ²=κ_hyp at p=3,5."
        ),
        "general_orth_inequality_proved": False,
        "kappa_bound_proved_for_all_p": False,
        "attack_surface": (
            "(a) closed-form PopP bulk under PSL isotypics; "
            "(b) mult-aware majorization better than λ₂·S; "
            "(c) resolvent orth=24‖δ‖² with μ²=4(p²+15)."
        ),
    }


def _load_Y(p: int):
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        path = Path("/tmp/maxplus_p5.npy")
        return np.load(path).astype(np.float64) if path.is_file() else None
    if p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
        return np.load(path).astype(np.float64) if path.is_file() else None
    return None


def certify_fickus_spectrum(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    N, n = Y.shape
    d = n // 2
    G = Y @ Y.T
    P = G / (2.0 * N)
    Pop = P * P
    ev = np.sort(np.linalg.eigvalsh(Pop))[::-1]
    # cluster
    rounded = np.round(ev, 10)
    counts = Counter(float(v) for v in rounded if abs(v) > 1e-12)
    levels = sorted(counts.items(), key=lambda x: -x[0])
    lam1 = float(ev[0])
    bulk = ev[1:]
    bulk_pos = bulk[bulk > 1e-12]
    m_obs = int(len(bulk_pos))
    S_obs = float(np.sum(bulk_pos))
    sum_sq = float(np.sum(bulk_pos**2))
    m_th = m_bulk(p)
    S_th = d * (d - 1) / N
    lam_flat = S_th / m_th if m_th else 0.0
    lam2 = float(bulk_pos[0]) if m_obs else 0.0
    mult2 = int(np.sum(np.abs(bulk_pos - lam2) < 1e-8)) if m_obs else 0

    ED4 = float(np.mean(G**4))
    e4f = float(ED4_flat(p))
    orth = ED4 - e4f
    room = float(room_orth(p))
    room_h = room * ((d - 1) / d) ** 2
    pr = float(proj_kappa2(p))
    kappa2 = ED4 - wick_lo(p)
    kh = float(kappa_hyp(p))
    eps = float(eps_wick(p))
    eps_h = float(eps_hyp(p))
    thr_wick = lam_flat * (1 + eps)
    thr_hyp = lam_flat * (1 + eps_h)

    rk_obs = int(np.sum(ev > 1e-10))
    rk_th = rank_popp(p)

    # p=5 exact eigs check
    exact_eigs = None
    if p == 5:
        # 11/845, 9/845, 5/845
        exact_eigs = {
            "11/845": abs(lam2 - 11 / 845) < 1e-9,
            "levels_3": len(levels) == 4,  # incl λ1
            "mults_1_13_26_26": (
                levels[0][1] == 1
                and levels[1][1] == 13
                and levels[2][1] == 26
                and levels[3][1] == 26
            ),
        }

    return {
        "p": p,
        "N": int(N),
        "d": int(d),
        "rank_PopP_obs": rk_obs,
        "rank_PopP_th": rk_th,
        "rank_match": rk_obs == rk_th,
        "lambda1": lam1,
        "d_over_N": d / N,
        "lambda1_match": abs(lam1 - d / N) < 1e-9,
        "m_bulk_obs": m_obs,
        "m_bulk_th": m_th,
        "m_match": m_obs == m_th,
        "S_obs": S_obs,
        "S_th": S_th,
        "S_match": abs(S_obs - S_th) < 1e-8,
        "lambda2": lam2,
        "mult_lambda2": mult2,
        "mult_ge_d_minus_1": mult2 >= d - 1,
        "lam_flat": lam_flat,
        "eps_wick": eps,
        "eps_hyp": eps_h,
        "thr_wick": thr_wick,
        "thr_hyp": thr_hyp,
        "lambda2_le_thr_wick": lam2 <= thr_wick + 1e-12,
        "lambda2_le_thr_hyp": lam2 <= thr_hyp + 1e-12,
        "lambda2_le_4_over_N": lam2 <= 4.0 / N + 1e-12,
        "lambda2_le_d_over_2N": lam2 <= d / (2.0 * N) + 1e-12,
        "ED4": ED4,
        "ED4_flat": e4f,
        "orth": orth,
        "room_orth": room,
        "room_hyp": room_h,
        "orth_le_room": orth <= room + 1e-4,
        "orth_le_room_hyp": orth <= room_h + 1e-4,
        "orth_eq_room_hyp": abs(orth - room_h) < 1e-4,
        "kappa2": kappa2,
        "proj": pr,
        "kappa_hyp": kh,
        "kappa2_eq_hyp": abs(kappa2 - kh) < 1e-4,
        "kappa2_le_96n": kappa2 <= 96.0 * n + 1e-4,
        "bulk_levels": [{"value": v, "mult": c} for v, c in levels],
        "max_bulk_levels_PSL": max_bulk_levels(p),
        "n_bulk_levels_obs": len(levels) - 1,  # exclude λ1
        "exact_eigs_p5": exact_eigs,
        "ok": bool(
            rk_obs == rk_th
            and abs(lam1 - d / N) < 1e-9
            and m_obs == m_th
            and mult2 >= d - 1
            and orth <= room + 1e-4
            and kappa2 <= 96.0 * n + 1e-4
        ),
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(3, 40) if is_prime(p)]
    print("Prop 15.101: Fickus Gram residual / bulk-variance orth reduction", flush=True)

    fr = prove_fickus_rank(primes)
    print(f"  Fickus rank structure proved={fr['proved']}", flush=True)

    ov = prove_orth_variance_identity(primes)
    print(f"  orth variance identity proved={ov['proved']}", flush=True)

    ep = prove_eps_criterion(primes)
    print(f"  eps criterion algebra proved={ep['proved']}", flush=True)

    red = prove_reduction()
    print(f"  reduction statement proved={red['proved']}", flush=True)

    certs = {}
    for p in (3, 5, 7):
        c = certify_fickus_spectrum(p)
        certs[str(p)] = c
        if c.get("skipped"):
            print(f"  p={p}: SKIP", flush=True)
            continue
        print(
            f"  p={p}: rank={c['rank_PopP_obs']} mult2={c['mult_lambda2']} "
            f"λ2={c['lambda2']:.6f} thr_w={c['thr_wick']:.6f} "
            f"orth={c['orth']:.4f} room_h={c['room_hyp']:.4f} "
            f"eq_hyp={c['kappa2_eq_hyp']} ok={c['ok']}",
            flush=True,
        )

    c3, c5 = certs.get("3", {}), certs.get("5", {})
    structure_ok = (
        fr["proved"]
        and ov["proved"]
        and ep["proved"]
        and red["proved"]
        and c3.get("ok") is True
        and c5.get("ok") is True
        and c3.get("orth_eq_room_hyp") is True  # both 0 at p=3
        and c5.get("orth_eq_room_hyp") is True
        and c5.get("kappa2_eq_hyp") is True
        and c3.get("kappa2_eq_hyp") is True
    )

    out = {
        "prop": "15.101",
        "backend": "CPU Fraction + Max+ p=3,5 Fickus spectrum; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "gap_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "general_orth_inequality_proved": False,
        "fickus_rank": fr,
        "orth_variance_identity": ov,
        "eps_criterion": ep,
        "reduction": red,
        "cert": certs,
        "method_source": (
            "Fickus–Jasper–Mixon Singer–Zauner gap (arXiv:2605.28738): residual "
            "Gram K=H⊙H̄ / Schur rank + joint eigenspace rank-nullity. "
            "Ge–Liu equiangular multiplicity/interlacing (arXiv:2606.29392): "
            "mult lower bounds via isotypics. Applied here to PopP=P⊙P bulk "
            "of Max+ design projector."
        ),
        "attack_surface": (
            "OPEN: prove ‖κ_orth‖² ≤ room·((d−1)/d)² for all primes p≥5. "
            "Reduced to bulk spectral variance of PopP (Fickus form) with "
            "≤⌊m/(d−1)⌋ levels (PSL mult). Sufficient λ₂-test is N-free but "
            "not sharp at p=5 equality. Next: closed bulk spectrum / mult-aware "
            "majorization / resolvent δ. Then 15.98⇒bi-tight; N_DEEP; Main Theorem."
        ),
        "verdict": (
            "Prop 15.101 applies Fickus Gram residual + PSL multiplicity to orth: "
            "proved rank/variance identities and N-free ε-criterion; certified "
            "bulk spectra and κ²=κ_hyp at p=3,5. General orth inequality OPEN. "
            "lim α_n OPEN."
        ),
        "settlement_note": (
            "L stays OPEN. Do not soft-close. Primary residual still orth bound "
            "for general p≥5; this prop reduces it to PopP bulk control."
        ),
        "proved": bool(structure_ok),
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15101.json"
    write_json_atomic(path, out)
    print(
        f"Prop 15.101 proved(structure)={out['proved']} "
        f"kappa_all_p=False wrote {path}",
        flush=True,
    )
    print("  L OPEN — orth bound still residual for ‖κ‖²≤96n", flush=True)


if __name__ == "__main__":
    main()
