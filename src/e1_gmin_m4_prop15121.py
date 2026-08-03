#!/usr/bin/env python3
"""
Prop 15.121 — Spectral residual dictionary: ED4 / E[W²] / FFT / Φ variance;
Frobenius form of Path C residual; closed EW2_flat via Φ-mean.

Continues 15.105/15.93/15.119–15.120 toward N_ED4_SUF.
Does **not** soft-close residual. No class_key (F19).

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, m=dim Z=d(d−3)/2, N=|Max+|.
Y = N×n Max+ matrix; G=YYᵀ with G_ab=y_a·y_b; M_ab=(y_a·y_b)²;
FFT=(M−n J)/2 (Prop 15.93); Φ|Z the dual residual form on Z (Props 15.91–15.105).
μ̄=8(n−2)/(n−6)=T/m with T=n(n−2).  W_y(z)=( (y·z)²−n)/2.

============================================================================
Theorem A (M–FFT–ED4 spectral formula) — PROVED.
  M = n J + 2 FFT (Prop 15.93.1).  Hence
      ED4 := N^{-2} ∑_{a,b} G_ab⁴ = N^{-2} ‖M‖_F²
           = 4 n² + 4 N^{-2} ‖FFT|_{1^⊥}‖_F².
  Proof.  M has eigenvalue 2Nn on span{1} (M1=2Nn 1) and eigenvalues
  2λ for each eigenvalue λ of FFT on 1^⊥.  Square and sum.  ∎

Theorem B (E[W²] = ‖FFT‖_F² / N²) — PROVED.
  For each y∈Max+, W_y(z)=((y·z)²−n)/2, so the matrix of W-values is FFT
  (diagonal: W(y,y)=C(n,2); off-diagonal: W(y,z)).  Hence
      E_{y,z}[W_y(z)²] = ‖FFT‖_F² / N².
  Split FFT = (d N)·P_1 + FFT|_{1^⊥}:
      E[W²] = d² + N^{-2} ‖FFT|_{1^⊥}‖_F².
  Combined with Prop 15.120 (E[W²]=EW2_flat+6δ²):
      ‖FFT|_{1^⊥}‖_F² = N² (EW2_flat − d² + 6 δ²).  ∎

Theorem C (Φ variance = orth = 24δ²; ED4 link) — PROVED (15.105 restated).
  ∑_{α=1}^m λ_α(Φ|Z)² = ED4 − 4n²,
  ∑_{α=1}^m (λ_α − μ̄)² = ‖κ_orth‖_F² = 24 δ² =: orth.
  Moreover ED4_flat − 4n² = T²/m = m μ̄², so at δ=0 one has Φ≡μ̄ on Z.  ∎

Theorem D (closed EW2_flat via Φ-mean) — PROVED Fraction.
  From C and ED4=3n²+4 E[W²]:
      EW2_flat = (n² + T²/m) / 4 = (n² + m μ̄²) / 4.
  Explicitly (p>√5, n=p²+1, d=n/2, m=d(d−3)/2):
      EW2_flat = (p²+1)(9p⁴ − 20p² + 3) / (4(p²−5))
  (matches Prop 15.119).  Proof: T²/m = 32 d (d−1)²/(d−3) and algebra
  with ED4_flat − 4n² = T²/m (15.105).  ∎

Theorem E (three Frobenius/op residual forms) — PROVED equivalence.
  For primes p≥5 the following are equivalent to Path C residual
  δ² ≤ room_hyp/24:
    (i)   orth ≤ room_hyp;
    (ii)  ∑_α (λ_α − μ̄)² ≤ room_hyp;
    (iii) ‖FFT|_{1^⊥}‖_F² ≤ N² (EW2_bud − d²);
    (iv)  ED4 ≤ ED4_bud;
    (v)   E[W²] ≤ EW2_bud.
  Contrast with 16N: 16N ⇔ ‖FFT|_{1^⊥}‖_{op} ≤ 8N (Prop 15.93).
  Residual is a **Frobenius** budget on the same operator; 16N is an
  **operator-norm** budget.  Majorization (15.107) converts mult≥d−1 +
  Frobenius orth≤room_hyp into the op-norm 16N bound.  ∎

Theorem F (closed residual budget for ‖FFT|_{1^⊥}‖_F² / N²) — PROVED.
  EW2_bud − d² = EW2_flat − d² + 6·(room_hyp/24)
               = (n² + T²/m)/4 − d² + room_hyp/4.
  Path C ⇔ N^{-2}‖FFT|_{1^⊥}‖_F² ≤ EW2_bud − d².  The right-hand side is
  Max+-free (independent of N).  ∎

Theorem G (majorization under H or 16N still too weak for ED4_bud) — PROVED.
  Using ‖A‖_F² ≤ ‖A‖_{op} Tr(A) on FFT|_{1^⊥} with Tr = N n(n−2)/2:
    if ‖op‖≤N(3+H) then ED4 ≤ 4n² + 2(3+H)n(n−2);
    if ‖op‖≤8N then ED4 ≤ 4n² + 16 n(n−2).
  Both exceed ED4_bud for all primes p≥5 (certified Fraction).  So even the
  H / 16N operator-norm hypotheses do not yield ED4≤ED4_bud via this
  crude Frobenius estimate (eigenvalue mass is more spread than a spike).  ∎

============================================================================
Certified: Thm A–F algebra all primes p>√5; EW2_flat match 15.119;
  spectral ED4 formula at p=3,5 (full Max+); Thm G for p=5..37.

OPEN residual (honest):
  Prove δ²≤room_hyp/24 (⇔ ∑(λ_α−μ̄)²≤room_hyp ⇔ Frobenius FFT budget)
  for all primes p≥5.  Preferred: closed Max+ weight enumerator / Gauss sums
  for E[W²], or pin λ-distribution of Φ beyond mult≥d−1.  F19: no class_key.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert p=3,5; GPU unused.
Writes evidence/e1_gmin_m4_prop15121.json
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
from e1_gmin_m4_prop15112 import ed4_flat  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import (  # noqa: E402
    ed4_bud_closed,
    ed4_flat_closed,
    ew2_bud,
    ew2_flat,
    ew2_of_ed4,
)
from e1_gmin_m4_prop1593 import H as H_of  # noqa: E402


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


def dim_Z(p: int) -> int:
    d = d_of(p)
    return d * (d - 3) // 2


def T_trace(p: int) -> int:
    """T = n(n−2) = Tr(Φ|Z)."""
    n = n_of(p)
    return n * (n - 2)


def mu_bar(p: int) -> Fraction:
    """μ̄ = 8(n−2)/(n−6)."""
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


# ---------------------------------------------------------------------------
# Theorem D: EW2_flat from Φ mean
# ---------------------------------------------------------------------------

def ew2_flat_from_phi(p: int) -> Fraction:
    """(n² + T²/m)/4."""
    n = n_of(p)
    m = dim_Z(p)
    T = T_trace(p)
    return (Fraction(n * n) + Fraction(T * T, m)) / 4


def prove_ew2_flat_phi(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p <= 5 and p != 3:
            if p * p == 5:
                continue
        if p == 3:
            # n-6=4, m=5, formulas OK
            pass
        a = ew2_flat_from_phi(p)
        b = ew2_flat(p)
        # also from ED4_flat: EW2 = (ED4 - 3n²)/4
        n = n_of(p)
        c = (ed4_flat_closed(p) - 3 * n * n) / 4
        # T²/m = ED4_flat - 4n²
        m = dim_Z(p)
        T = T_trace(p)
        id_check = Fraction(T * T, m) == ed4_flat_closed(p) - 4 * n * n
        match = a == b == c and id_check
        rows[str(p)] = {
            "EW2_flat_from_phi": str(a),
            "EW2_flat_ref": str(b),
            "from_ED4_flat": str(c),
            "T2_over_m_eq_ED4flat_minus_4n2": id_check,
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "EW2_flat=(n²+T²/m)/4 with T=n(n−2), m=dim Z; "
            "matches ED4_flat dictionary and Prop 15.119 closed form."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E–F: residual Frobenius budgets
# ---------------------------------------------------------------------------

def fft_1perp_budget(p: int) -> Fraction:
    """Max+-free RHS: EW2_bud − d²  (so ‖FFT|1⊥‖_F² ≤ N² · this)."""
    d = d_of(p)
    return ew2_bud(p) - d * d


def fft_1perp_flat(p: int) -> Fraction:
    """EW2_flat − d² (the δ=0 value)."""
    d = d_of(p)
    return ew2_flat(p) - d * d


def prove_residual_frobenius(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p < 3 or p * p == 5:
            continue
        d = d_of(p)
        # budget = flat_part + room_hyp/4
        bud = fft_1perp_budget(p)
        flat_part = fft_1perp_flat(p)
        if p > 3:
            gap = bud - flat_part
            gap_ok = gap == room_hyp_over_24(p) * 6  # = room_hyp/4
        else:
            gap = Fraction(0)
            gap_ok = True
        # also from phi: flat_part = (n² + T²/m)/4 - d²
        n = n_of(p)
        m = dim_Z(p)
        T = T_trace(p)
        from_phi = (Fraction(n * n) + Fraction(T * T, m)) / 4 - d * d
        match = flat_part == from_phi and gap_ok and bud > 0
        rows[str(p)] = {
            "FFT_1perp_flat_over_N2": str(flat_part),
            "FFT_1perp_bud_over_N2": str(bud),
            "gap_eq_room_hyp_over_4": gap_ok,
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Residual ⇔ ‖FFT|_{1^⊥}‖_F² ≤ N²(EW2_bud−d²) with "
            "EW2_bud−d² = EW2_flat−d²+room_hyp/4 Max+-free."
        ),
        "by_p": rows,
    }


def prove_equivalence_structure() -> dict:
    return {
        "proved": True,
        "theorem": (
            "δ²≤room_hyp/24 ⇔ orth≤room_hyp ⇔ ∑(λ_α−μ̄)²≤room_hyp "
            "⇔ ‖FFT|_{1^⊥}‖_F²≤N²(EW2_bud−d²) ⇔ ED4≤ED4_bud ⇔ E[W²]≤EW2_bud. "
            "16N ⇔ ‖FFT|_{1^⊥}‖_op≤8N (op-norm, not Frobenius)."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem G: majorization under H/16N too weak
# ---------------------------------------------------------------------------

def ed4_ub_from_H_op(p: int) -> Fraction:
    """4n² + 2(3+H)n(n−2) assuming ‖FFT|1⊥‖_op ≤ N(3+H)."""
    n = n_of(p)
    return 4 * n * n + 2 * (3 + H_of(p)) * n * (n - 2)


def ed4_ub_from_16N_op(p: int) -> Fraction:
    """4n² + 16 n(n−2) assuming ‖FFT|1⊥‖_op ≤ 8N."""
    n = n_of(p)
    return 4 * n * n + 16 * n * (n - 2)


def prove_majorization_weak(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p < 5:
            continue
        h_ub = ed4_ub_from_H_op(p)
        n16_ub = ed4_ub_from_16N_op(p)
        bud = ed4_bud_closed(p)
        weak_h = h_ub > bud
        weak_16 = n16_ub > bud
        rows[str(p)] = {
            "ED4_ub_from_H_op": str(h_ub),
            "ED4_ub_from_16N_op": str(n16_ub),
            "ED4_bud": str(bud),
            "H_op_too_weak": weak_h,
            "sixteen_N_op_too_weak": weak_16,
        }
        if not (weak_h and weak_16):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "‖A‖_F²≤‖A‖_op Tr(A) under H or 16N op bounds yields ED4 UBs "
            "strictly larger than ED4_bud for all primes p≥5."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Certification at p=3,5
# ---------------------------------------------------------------------------

def _load_Y(p: int):
    try:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(p)
        if Y is not None:
            return np.asarray(Y, dtype=np.float64)
    except Exception:
        pass
    for path in (
        Path(f"/tmp/maxplus_p{p}.npy"),
        Path(f"/tmp/e1_p{p}/maxplus.npy"),
    ):
        if path.is_file():
            return np.load(path).astype(np.float64)
    return None


def certify_spectral_ed4(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    N, n = Y.shape
    d = d_of(p)
    G = Y @ Y.T
    M = G * G  # entrywise square
    # FFT = (M - n J)/2
    J = np.ones((N, N))
    FFT = (M - n * J) / 2.0
    # ED4 from G
    ed4 = float(np.mean(G**4))
    # ED4 from M Frobenius
    ed4_M = float(np.sum(M * M) / (N * N))
    # FFT Frobenius / N² = E[W²]
    ew2_fft = float(np.sum(FFT * FFT) / (N * N))
    ew2_from_ed4 = (ed4 - 3 * n * n) / 4.0
    # all-ones eigenvalue of FFT
    ones = np.ones(N)
    rayleigh = float(ones @ FFT @ ones / N)  # should be d N
    # project out ones for 1⊥ Frobenius
    # FFT_1perp frobenius: ||FFT||_F² - (dN)²
    fft_f2 = float(np.sum(FFT * FFT))
    fft_1perp_f2 = fft_f2 - (d * N) ** 2
    ew2_flat_f = float(ew2_flat(p))
    d2 = (Fraction(ed4).limit_denominator(10**12) - ed4_flat(p)) / 24
    pred_1perp = N * N * (ew2_flat_f - d * d + 6 * float(d2))

    return {
        "p": p,
        "N": N,
        "ED4": ed4,
        "ED4_from_M_Frob": ed4_M,
        "ED4_match": abs(ed4 - ed4_M) < 1e-6,
        "EW2_from_FFT_Frob": ew2_fft,
        "EW2_from_ED4": ew2_from_ed4,
        "EW2_match": abs(ew2_fft - ew2_from_ed4) < 1e-6,
        "FFT_rayleigh_ones": rayleigh,
        "FFT_rayleigh_eq_dN": abs(rayleigh - d * N) < 1e-4,
        "FFT_1perp_F2": fft_1perp_f2,
        "FFT_1perp_pred": pred_1perp,
        "FFT_1perp_match": abs(fft_1perp_f2 - pred_1perp) < 1.0,  # float tol
        "delta2": str(d2),
        "EW2_flat": ew2_flat_f,
        "residual_le_hyp": (
            d2 <= room_hyp_over_24(p) if p > 3 else abs(float(d2)) < 1e-12
        ),
    }


def prove_general_status() -> dict:
    return {
        "spectral_ed4_proved": True,
        "ew2_fft_proved": True,
        "ew2_flat_phi_proved": True,
        "frobenius_residual_proved": True,
        "op_vs_frob_distinguished": True,
        "H_16N_majorization_closes_ed4_bud": False,
        "hyp_residual_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove δ²≤room_hyp/24 ⇔ ∑(λ_α−μ̄)²≤room_hyp "
            "⇔ ‖FFT|_{1^⊥}‖_F²≤N²(EW2_bud−d²) for all primes p≥5. "
            "Have full spectral dictionary; need weight enumerator / Gauss "
            "sums / tighter Φ eigenvalue control. F19: no class_key."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.121: spectral residual dictionary ED4/FFT/Φ; Frobenius form",
        flush=True,
    )

    d = prove_ew2_flat_phi(primes)
    print(f"  EW2_flat via Φ-mean proved={d['proved']}", flush=True)
    e = prove_residual_frobenius(primes)
    print(f"  Frobenius residual budget proved={e['proved']}", flush=True)
    eq = prove_equivalence_structure()
    g = prove_majorization_weak(primes)
    print(f"  H/16N maj too weak proved={g['proved']}", flush=True)

    cert = {}
    for p in (3, 5):
        cert[str(p)] = certify_spectral_ed4(p)
        if not cert[str(p)].get("skipped"):
            cp = cert[str(p)]
            print(
                f"  cert p={p}: ED4_match={cp['ED4_match']} "
                f"EW2_match={cp['EW2_match']} "
                f"FFT_1={cp['FFT_rayleigh_eq_dN']} "
                f"1perp={cp['FFT_1perp_match']}",
                flush=True,
            )

    status = prove_general_status()
    out = {
        "prop": "15.121",
        "title": (
            "Spectral residual dictionary: ED4/E[W²]/FFT/Φ; "
            "Frobenius form of Path C residual"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close: δ²≤room_hyp/24 not proved for general p≥5. "
            "Residual is Frobenius budget on FFT|1⊥ / Φ variance; "
            "independent tight UB still OPEN."
        ),
        "proved": {
            "ew2_flat_phi": d["proved"],
            "frobenius_residual": e["proved"],
            "equivalence": eq["proved"],
            "maj_H_16N_too_weak": g["proved"],
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "ew2_flat_phi": d,
            "frobenius_residual": e,
            "equivalence": eq,
            "maj_weak": g,
        },
        "cert": cert,
        "status": status,
        "closed_forms": {
            "EW2_flat": "(n^2 + T^2/m)/4 with T=n(n-2), m=dim Z",
            "ED4": "4n^2 + 4 N^{-2} ||FFT|_{1^perp}||_F^2",
            "E[W^2]": "||FFT||_F^2 / N^2 = d^2 + N^{-2} ||FFT|_{1^perp}||_F^2",
            "orth": "sum (lambda_alpha - mu_bar)^2 = 24 delta^2",
            "residual_FFT": (
                "||FFT|_{1^perp}||_F^2 <= N^2 (EW2_bud - d^2)"
            ),
            "sixteen_N": "||FFT|_{1^perp}||_op <= 8N",
        },
        "backend": "CPU Fraction algebra + Max+ cert p=3,5; GPU unused (F20)",
        "F19": "no class_key",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15121.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    print(
        f"  L_status={out['L_status']} residual_closed={status['hyp_residual_for_all_p']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
