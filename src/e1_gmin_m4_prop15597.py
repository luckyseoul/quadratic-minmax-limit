#!/usr/bin/env python3
r"""
Prop 15.597 — Φ_part = λ̄·I exactly: the particular solution is spectrally
invisible, so leftovers 1 and 3 are statements about Φ_δ alone.

Does **not** flip phi_F_ge_6 / type_I / residual_ii / e1 / L. Soft-close forbidden.

SETUP (15.593 Φ; 15.594 δ; 15.247 A m₄_part)
  For W ∈ Z: diag(W)=0 and CW = pW, so tr(P₊W) = tr(W)/2 = 0, hence
        ⟨B̃_y, W⟩ = ⟨y yᵀ − 2P₊, W⟩ = yᵀ W y .
  Therefore Φ's quadratic form is  ⟨W, ΦW⟩ = E_y[(yᵀWy)²],  and Φ is
  LINEAR in the four-point moment tensor m₄.  With m₄ = m₄_part + δ:
        Φ = Φ_part + Φ_δ .

PROVED
  A. **Φ_part = λ̄ I on Z**, λ̄ = 8(n−2)/(n−6).
     Verified exactly at p=5 (dim Z = 65): the full spectrum of Φ_part,
     built DATA-FREE from κ, φ, star alone, is the single value 9.6 =
     λ̄(5) with spread 4.1e−14 and ‖Φ_part − λ̄I‖_F = 1.6e−14.
     Verified DIRECTLY at p=7 by the quadratic form on random W ∈ Z
     (`scripts/phi_part_scalar_check.py`, no dim-Z basis needed, no Max±
     data): ⟨W,Φ_part W⟩/‖W‖² = 8.727272727 = λ̄(7) on three independent
     random W, deviation 0.0e0 / 0.0e0 / 1.8e−15 (Z-residual ≤ 7e−15).
     Verified DIRECTLY at p=11 the same way (8.8M four-sets, chunked):
     ⟨W,Φ_part W⟩/‖W‖² = 8.275862069 = λ̄(11) on three independent random
     W, deviations 1.8e−15 / 0.0e0 / 1.8e−15 (Z-residual ≤ 1.6e−14).
     Verified DIRECTLY at p=13 the same way (33.6M four-sets):
     8.195121951 = λ̄(13), deviations 0.0 / 5.3e−15 / 0.0.
     So Theorem A is directly verified at p=5,7,11,13 — data-free at
     every prime, needing no Max± ensemble.

  A''. **PROVED FOR ALL p (2026-08-22).**  See Theorem A* below; the
      contraction identities are now closed in closed form, so
      Phi_part = lambda_bar I is a general-p THEOREM, not a census.

  A'. (superseded by A*, retained for the reduction it records)
      For W ∈ Z:
        (i)  W C = p W  (transpose of CW = pW, both symmetric), hence
             W²C = pW² and  tr(W²C)/p = ‖W‖²_F.
        (ii) **(C∘W)·1 = 0**: every row sum of the Hadamard product
             vanishes, since Σ_j C_ij W_ij = Σ_j C_ij W_ji = (CW)_ii
             = p·W_ii = 0 (diag W = 0).
      With (i) the quadratic form collapses to
        ⟨W,Φ_part W⟩ = 8·Σ_S m₄_part(S)·t(S) + 6‖W‖²,
        t(S) := W_ij W_kl + W_ik W_jl + W_il W_jk,
      so Theorem A is EQUIVALENT to the single contraction identity
        Σ_S m₄_part(S) t(S) = (n+10)/(4(n−6)) · ‖W‖²_F ,
      i.e. (n−2)A − 2B − 2p·E = (n−1)(n+10)/4 · ‖W‖² with
      A=Σ_S κt, B=Σ_S φt, E=Σ_S star·t (using p²=n−1, D=(n−1)(n−6)).
      The κ-part is already done: its same-pairing piece equals
      ‖W‖²/4 by (ii) (all row sums R_i vanish), and its cross-pairing
      piece is governed by tr(CWCW) = p²‖W‖² = (n−1)‖W‖².  The φ- and
      star-contractions remain; both are traces of products of C and W
      subject to CW = pW, so the route is closed-form.  COMPLETED in A*.

THEOREM A* — Phi_part = lambda_bar * I for EVERY prime p >= 5 (PROVED).
  Lemmas (all for W in Z: W symmetric, diag W = 0, CW = pW; C^2 = (n-1)I):
    L1.  WC = pW (transpose), so W^2 C = pW^2 and tr(W^2 C)/p = ||W||^2.
    L2.  (CW)_ii = p W_ii = 0, and since C,W are symmetric
         sum_j C_ij W_ij = sum_j C_ij W_ji = (CW)_ii = 0:
         **every row sum of the Hadamard product C.W vanishes.**
    L3.  CWC = p(WC) = p^2 W, hence for c := C e_r (any column of C)
         c^T W c = (CWC)_rr = p^2 W_rr = 0,  and  W c = p W_{.r}.
  Contractions (t(S) := W_ij W_kl + W_ik W_jl + W_il W_jk):
    A = sum_S kappa(S) t(S) = (n+1)/4 * ||W||^2.
        Same-pairing part = sum over DISJOINT edge pairs of u_e u_f,
        u := C.W.  By L2 all row sums R_i = 0, so
        sum_{e<f} u_e u_f = -||W||^2/4 and sum_adjacent = -||W||^2/2,
        giving disjoint part = ||W||^2/4.  Cross-pairing part: two
        distinct perfect matchings of a 4-set union to a 4-cycle, and
        each C/W-alternating labelled 4-cycle is hit by exactly 4 index
        tuples of tr(CWCW), so it equals tr_distinct(CWCW)/4.  By
        inclusion-exclusion the a=c and b=d terms vanish (each is
        sum_a (CW)_aa^2 = 0 by L2) while a=c AND b=d contributes
        +||W||^2, so tr_distinct = tr(CWCW) + ||W||^2 = p^2||W||^2 +
        ||W||^2 = n||W||^2.  Hence A = ||W||^2/4 + n||W||^2/4.  QED
    B = sum_S phi(S) t(S) = -n/4 * ||W||^2.
        Fix r and set v_ij := C_ri C_rj W_ij on V minus {r}.  Then
        sum_e v_e = (1/2) c^T W c = 0 by L3, the v-row sums are
        p C_ri W_ir by L3 (W c = p W_{.r}), and sum_e v_e^2 =
        ||W||^2/2 - s_r with s_r := sum_j W_rj^2.  The same
        disjoint-pair identity gives ||W||^2/4 - (n/2) s_r per r
        (using 1 + p^2 = n); summing over r with sum_r s_r = ||W||^2
        gives (n/4)||W||^2 - (n/2)||W||^2.  QED
    E = sum_S star(S) t(S) = -p * ||W||^2.
        Fix s and d := C_{.s}.  For the term pairing c with s, the inner
        sum over {a,b} disjoint from {s,c} telescopes:
        sum_{a<b != s,c} d_a d_b W_ab = (1/2) d^T W d - d_c (W d)_c
        = 0 - d_c p W_cs by L3.  Hence E = -p sum_{s,c} d_c^2 W_cs^2
        = -p ||W||^2.  QED
  Combination (p^2 = n-1, D = p^2(p^2-5) = (n-1)(n-6)):
    (n-2)A - 2B - 2p E = (n-2)(n+1)/4 + n/2 + 2(n-1) = (n^2+9n-10)/4
                       = (n-1)(n+10)/4 ,
    so sum_S m4_part(S) t(S) = (n+10)/(4(n-6)) * ||W||^2, and with L1
    <W, Phi_part W> = 8*(n+10)/(4(n-6))*||W||^2 + 6||W||^2
                    = (8n-16)/(n-6) * ||W||^2 = lambda_bar * ||W||^2.
  Since W in Z was arbitrary, Phi_part = lambda_bar I on Z.  QED
  B. Consequence:  **all spectral deviation of Φ comes from δ.**  The
     particular solution — the entire explicit, Max-free part of the
     moment tensor — contributes a perfectly flat spectrum and cannot
     help or hurt either leftover.
  C. Sharpened criteria.  Since Φ = λ̄I + Φ_δ and Φ_δ is scalar on each
     constituent (equivariance + multiplicity-freeness):
        λ_min(Φ) = λ̄ + λ_min(Φ_δ),   ‖Φ_δ‖²_F = 24‖δ‖² ,
     so, with the 15.589 C multiplicity floor n/2,
        **leftover 1  ⟺  Φ_δ ⪰ −(2n+20)/(n−6) · I  on Z** ,
     and ‖Φ_δ‖_op ≤ √(48‖δ‖²/n) recovers the 15.595 threshold ‖δ‖² ≤ n/12.

WHY THIS DOES NOT FLIP ANYTHING (stated explicitly)
  D. A is a statement about the EXPLICIT part only.  It removes m₄_part
     from the problem but supplies no bound on δ, which is precisely the
     open content (δ spans ker(4pI − T), invisible to the master
     equation by construction — 15.595 F).  phi_F_ge_6_proved_general,
     type_I_aut_e_3AB_positive_general and multilevel_ND_k_ge_4p_proved
     all remain False; L remains OPEN.

OPEN (unchanged, now in its tightest form)
  E. Φ_δ ⪰ −(2n+20)/(n−6)·I on Z, equivalently ‖δ‖² ≤ n/12 via the
     multiplicity floor.  Measured ‖δ‖²/n = 0.9089, 0.2085, 0.01941 at
     p=5,7,11 (three points; p=5,7 exceed the threshold and are handed to
     census, p=11 clears it with 4.3× margin).  NOT extrapolated.

Writes evidence/e1_gmin_m4_prop15597.json
"""
from __future__ import annotations

import itertools
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    return p * p + 1


def lambda_bar(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def leftover1_operator_criterion(p: int) -> Fraction:
    """leftover 1 <=> Phi_delta >= -(2n+20)/(n-6) I."""
    n = n_of(p)
    return -Fraction(2 * n + 20, n - 6)


def phi_part_spectrum(p: int):
    """Build Phi_part DATA-FREE (kappa, phi, star only) and return its spectrum."""
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    n, Ci = lab.n, lab.C
    C = Ci.astype(np.float64)
    pairs = list(itertools.combinations(range(n), 2))
    A = []
    for (i, j) in pairs:
        E = np.zeros((n, n))
        E[i, j] = E[j, i] = 1
        A.append((C @ E - p * E).reshape(-1))
    Amat = np.array(A).T
    s = np.linalg.svd(Amat, compute_uv=False)
    _, _, Vt = np.linalg.svd(Amat)
    null = Vt[np.sum(s > 1e-8):]
    Q = []
    for v in null:
        B = np.zeros((n, n))
        for idx, (i, j) in enumerate(pairs):
            B[i, j] = B[j, i] = v[idx]
        for W in Q:
            B = B - (B * W).sum() * W
        nr = np.linalg.norm(B)
        if nr > 1e-8:
            Q.append(B / nr)
    d = len(Q)
    S4 = lab.S4
    i_, j_, k_, l_ = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    kap = (Ci[i_, j_] * Ci[k_, l_] + Ci[i_, k_] * Ci[j_, l_]
           + Ci[i_, l_] * Ci[j_, k_]).astype(np.int64)
    phi = np.zeros(len(S4), dtype=np.int64)
    star = np.zeros(len(S4), dtype=np.int64)
    for si, (a, b, c, dd) in enumerate(S4):
        col = Ci[:, a] * Ci[:, b] * Ci[:, c] * Ci[:, dd]
        col = col.copy()
        col[[a, b, c, dd]] = 0
        phi[si] = col.sum()
        star[si] = (Ci[b, a] * Ci[c, a] * Ci[dd, a] + Ci[a, b] * Ci[c, b] * Ci[dd, b]
                    + Ci[a, c] * Ci[b, c] * Ci[dd, c] + Ci[a, dd] * Ci[b, dd] * Ci[c, dd])
    D = p * p * (p * p - 5)
    m4p = ((p * p - 1) * kap - 2 * phi - 2 * p * star) / D
    pairmom = C / p
    W = np.array(Q)
    Phi = np.zeros((d, d))
    a4, b4, c4, d4_ = i_, j_, k_, l_
    for x in range(d):
        Wa = W[x]
        for y in range(x, d):
            Wb = W[y]
            t = (Wa[a4, b4] * Wb[c4, d4_] + Wa[c4, d4_] * Wb[a4, b4]
                 + Wa[a4, c4] * Wb[b4, d4_] + Wa[b4, d4_] * Wb[a4, c4]
                 + Wa[a4, d4_] * Wb[b4, c4] + Wa[b4, c4] * Wb[a4, d4_])
            val = 4 * np.sum(m4p * t)
            sh = 4 * np.sum((Wa @ Wb) * pairmom)
            Phi[x, y] = Phi[y, x] = val + 2 * (Wa * Wb).sum() + sh
    return d, np.linalg.eigvalsh(Phi), Phi


def theorem_A_phi_part_scalar(p: int, tol: float = 1e-8) -> dict:
    d, ev, Phi = phi_part_spectrum(p)
    lb = float(lambda_bar(p))
    return {"p": p, "dimZ": d, "lambda_bar": lb,
            "spec_min": float(ev.min()), "spec_max": float(ev.max()),
            "spread": float(ev.max() - ev.min()),
            "frob_dist_to_scalar": float(np.linalg.norm(Phi - lb * np.eye(d))),
            "is_scalar": bool(np.allclose(Phi, lb * np.eye(d), atol=tol))}


def main():
    t0 = time.time()
    r5 = theorem_A_phi_part_scalar(5)
    out = {
        "prop": "15.597",
        "title": "Phi_part = lambda_bar I: the particular solution is spectrally invisible",
        "p5": r5,
        "p7_forced_by": "G-equivariance + multiplicity-free Z + V = 24||delta||^2 (15.594)",
        "leftover1_operator_criterion": {
            str(p): str(leftover1_operator_criterion(p)) for p in (5, 7, 11, 13, 17)},
        "measured_delta_sq_per_n": {5: 0.9089, 7: 0.2085, 11: 0.01941},
        "flips_anything": False,
        "flags_not_flipped": ["phi_F_ge_6", "type_I", "residual_ii", "e1", "L"],
        "L_status": "OPEN",
        "seconds": round(time.time() - t0, 1),
    }
    (ROOT / "evidence" / "e1_gmin_m4_prop15597.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.597  Phi_part = lambda_bar I (particular solution spectrally invisible)")
    print(f"  p=5: dimZ={r5['dimZ']} lbar={r5['lambda_bar']:.6f} "
          f"spread={r5['spread']:.2e} ||Phi_part-lbar I||_F={r5['frob_dist_to_scalar']:.2e} "
          f"scalar={r5['is_scalar']}")
    print("  => all spectral deviation of Phi comes from delta alone.")
    print("  leftover 1 <=> Phi_delta >= -(2n+20)/(n-6) I on Z.  Still OPEN; no flag flipped.")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()


# ----------------------------------------------------------- Theorem A* checks
def contraction_closed_forms(p: int) -> dict:
    """The three proved contraction values (Theorem A*)."""
    n = n_of(p)
    return {"A": Fraction(n + 1, 4), "B": Fraction(-n, 4), "E": Fraction(-p)}


def theorem_A_star_algebra(p: int) -> dict:
    """Verify the PROVED combination closes symbolically for this p."""
    n = n_of(p)
    c = contraction_closed_forms(p)
    comb = (n - 2) * c["A"] - 2 * c["B"] - 2 * p * c["E"]
    target = Fraction((n - 1) * (n + 10), 4)
    m4p_contraction = comb / Fraction((n - 1) * (n - 6))
    qf = 8 * m4p_contraction + 6
    return {"p": p, "n": n,
            "combination": comb, "target": target, "combination_ok": comb == target,
            "m4part_contraction": m4p_contraction,
            "m4part_target": Fraction(n + 10, 4 * (n - 6)),
            "qf": qf, "lambda_bar": lambda_bar(p), "qf_equals_lambda_bar": qf == lambda_bar(p)}
