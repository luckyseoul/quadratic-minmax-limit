#!/usr/bin/env python3
"""
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
