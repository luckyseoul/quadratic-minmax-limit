#!/usr/bin/env python3
"""
Prop 15.593 — Es4 = 4n² + tr(Φ²); the exact design floor; both leftovers
reduce to the Frobenius spectral variance V of Φ on Z.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SETUP  (15.589 Z-decomposition; 15.592 Es4 reduction; s := y·z on Max+)
  P₊ = (I + C/p)/2 projects onto V₊.  For y ∈ Max+ put
        B̃_y := y yᵀ − 2 P₊ .

PROVED (exact identities; verified at p=5,7 — 15.590's MuLab ensembles)
  A. B̃_y ∈ Z pointwise.  diag(y yᵀ)=I and diag(2P₊)=I ⇒ zero diagonal;
     C(y yᵀ) = (Cy)yᵀ = p·y yᵀ and C P₊ = p P₊ ⇒ C B̃_y = p B̃_y.
     Moreover ‖B̃_y‖²_F = n² − 2n for EVERY y (no averaging).
  B. s² = ⟨y yᵀ, z zᵀ⟩ = 2n + ⟨B̃_y, B̃_z⟩, hence with the Gram operator
     Φ := E_y |B̃_y⟩⟨B̃_y| on Z (the 15.589 Φ):
        Es4 := E_{y,z}[s⁴] = 4n² + tr(Φ²),   tr Φ = n(n−2),
        λ̄ = tr Φ / dim Z = 8(n−2)/(n−6)   (15.589 spectral mean).
  C. **Design floor (new, exact).**  Cauchy–Schwarz tr(Φ²) ≥ (trΦ)²/dim Z:
        Es4  ≥  12n² + 16n + 128n/(n−6),
     with equality iff Φ is scalar on Z ("Max+ is a perfect 4-design in
     the Z-sense").  Verified: floor excess per n = 16 + 128/(n−6)
     = 22.40, 18.91, 17.10 at p=5,7,11 vs TRUE 44.21, 23.91, 17.57.
  D. Writing V := ‖Φ − λ̄ I‖²_F = Es4 − floor ≥ 0, multiplicity-freeness
     of Z (15.589 B) gives the exact constituent decomposition
        V = (n/2)(λ_exc − λ̄)² + n · Σ_α (λ_α − λ̄)² ,
     the sum over the (q−9)/8 principal-series constituents.  (Verified
     at p=5,7 against the full Φ spectrum; note p=7 has two COINCIDENT
     principal eigenvalues — mult 2n — which must be counted twice.)

CONSEQUENCE — both leftovers are the same variance bound, with constants
  E. **Leftover 1.**  If λ_min < 6 then, with m its multiplicity,
     V ≥ m(λ̄ − λ_min)² > m(λ̄ − 6)².  By 15.589 C the only possible
     sub-n multiplicity is λ_exc with m = n/2, so
        V ≤ (n/2)(λ̄ − 6)²  ⟹  λ_min ≥ 6 ,
     and (n/2)(λ̄−6)²/n = (2n+20)²/(2(n−6)²) → 2.  Threshold c₁ ≈ 2.
  F. **Leftover 3.**  By 15.592, leftover 3 ⟸ census(5,7) + Es4 ≤
     12n² + x(p)n, i.e.
        V ≤ (x(p) − 16 − 128/(n−6))·n .   Threshold c₃(11) = 15.50,
     growing to ≈ 69.  So c₃ > c₁ at every prime: **leftover 1's
     variance bound strictly implies leftover 3's** — the shared blocker
     of fable.md, now with explicit constants.
  G. Measured V/n = 21.81 (p=5), 5.00 (p=7), 0.4659 (p=11): collapsing.
     At p=11, Es4 sits within 2.7% of the design floor, and BOTH
     thresholds hold with margin (leftover 1: 5.6×, leftover 3: 33.3×).
     Recorded, NOT extrapolated (three points).

WHY 12 IS FORCED (interpretation, not needed for the logic)
  H. For vectors of norm √n uniform in a d-dimensional space,
     E[(y·z)⁴] = 3n⁴/(d(d+2)); at d = dim V₊ = n/2 this is
     12n³/(n+4) = 12n² − 48n + O(1).  So the leading constant 12 is the
     Gaussian/4-design value in V₊, and "Es4 = 12n² + O(n)" is exactly
     the statement that Max+ is 4-design-like in V₊.  This is why
     15.592 F holds: any majorant with leading constant (12+ε)n² is
     structurally, not numerically, insufficient.

OPEN (the one hard core, now for BOTH leftovers 1 and 3)
  I. Prove V = O(n), i.e.  (½)(λ_exc − λ̄)² + Σ_α (λ_α − λ̄)² = O(1)
     over the ≈ n/8 constituents.  Equivalently: the energies
     ‖P_c B̃_y‖² equidistribute across the PSL(2,q)-constituents of Z
     up to O(1) in the above sense.  Delsarte/2-design input is NOT
     enough (15.590's kill).

Writes evidence/e1_gmin_m4_prop15593.json
"""
from __future__ import annotations

import itertools
import json
import os
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


def design_floor(p: int) -> Fraction:
    """Es4 >= 12n^2 + 16n + 128n/(n-6), exact."""
    n = n_of(p)
    return Fraction(12 * n * n) + Fraction(16 * n) + Fraction(128 * n, n - 6)


def floor_via_cauchy_schwarz(p: int) -> Fraction:
    n = n_of(p)
    dZ = n * (n - 6) // 8
    return Fraction(4 * n * n) + Fraction(n * (n - 2)) ** 2 / dZ


def threshold_leftover1(p: int) -> Fraction:
    """V <= (n/2)(lambda_bar - 6)^2 implies lambda_min >= 6."""
    n = n_of(p)
    return Fraction(n, 2) * (lambda_bar(p) - 6) ** 2


def threshold_leftover3(p: int, x_of_p: Fraction) -> Fraction:
    """V <= (x(p) - 16 - 128/(n-6)) n  gives Es4 <= 12n^2 + x(p) n."""
    n = n_of(p)
    return x_of_p * n - Fraction(16 * n) - Fraction(128 * n, n - 6)


def theorem_A_B_pointwise(p: int) -> dict:
    """B_y in Z with constant norm; and the s^2 identity."""
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    n, C = lab.n, lab.C.astype(np.float64)
    Y = lab.Yp.astype(np.float64)
    Pp = (np.eye(n) + C / p) / 2
    ok_diag = ok_eig = ok_norm = True
    for y in Y[: min(50, len(Y))]:
        B = np.outer(y, y) - 2 * Pp
        ok_diag &= abs(np.diag(B)).max() < 1e-10
        ok_eig &= abs(C @ B - p * B).max() < 1e-8
        ok_norm &= abs((B * B).sum() - (n * n - 2 * n)) < 1e-6
    S = (Y @ Y.T)
    lhs = S * S
    rhs = 2 * n + np.array([[float(((np.outer(y, y) - 2 * Pp) *
                                    (np.outer(z, z) - 2 * Pp)).sum())
                             for z in Y[:8]] for y in Y[:8]])
    ok_s2 = abs(lhs[:8, :8] - rhs).max() < 1e-6
    return {"p": p, "B_in_Z": bool(ok_diag and ok_eig),
            "constant_norm": bool(ok_norm), "s2_identity": bool(ok_s2),
            "proved": bool(ok_diag and ok_eig and ok_norm and ok_s2)}


def es4_and_V(p: int) -> dict:
    """Exact Es4, tr(Phi^2), V from the ensemble."""
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    n = lab.n
    Y = lab.Yp.astype(np.float64)
    Ny = len(Y)
    S = (Y @ Y.T).astype(np.int64).astype(object)
    Es4 = Fraction(int((S ** 4).sum()), Ny * Ny)
    trPhi2 = Fraction(int(((S ** 2 - 2 * n) ** 2).sum()), Ny * Ny)
    V = Es4 - design_floor(p)
    return {"p": p, "Es4": Es4, "trPhi2": trPhi2,
            "identity_Es4": Es4 == 4 * n * n + trPhi2,
            "floor_agrees": design_floor(p) == floor_via_cauchy_schwarz(p),
            "V": V, "V_per_n": float(V / n),
            "Es4_excess_per_n": float((Es4 - 12 * n * n) / n),
            "floor_excess_per_n": float((design_floor(p) - 12 * n * n) / n)}


def phi_spectrum(p: int):
    """Full Phi spectrum on an explicit orthonormal basis of Z."""
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    n, C = lab.n, lab.C.astype(np.float64)
    Y = lab.Yp.astype(np.float64)
    pairs = list(itertools.combinations(range(n), 2))
    A = []
    for (i, j) in pairs:
        E = np.zeros((n, n))
        E[i, j] = E[j, i] = 1
        A.append((C @ E - p * E).reshape(-1))
    s = np.linalg.svd(np.array(A).T, compute_uv=False)
    _, _, Vt = np.linalg.svd(np.array(A).T)
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
    qf = np.array([[float(y @ B @ y) for B in Q] for y in Y])
    Phi = (qf.T @ qf) / len(Y)
    return n, np.linalg.eigvalsh(Phi), Phi


def theorem_D_decomposition(p: int) -> dict:
    """V = (n/2)(l_exc - lbar)^2 + n sum_a (l_a - lbar)^2, counting
    coincident principal eigenvalues with their full multiplicity."""
    n, ev, Phi = phi_spectrum(p)
    lbar = float(lambda_bar(p))
    V_spec = float(((ev - lbar) ** 2).sum())
    # cluster, then split each cluster into constituents of size n (principal)
    cl = []
    for e in sorted(ev):
        if cl and abs(e - cl[-1][0]) < 1e-6:
            cl[-1][1] += 1
        else:
            cl.append([e, 1])
    exc = [(a, m) for a, m in cl if m == n // 2]
    prin = [(a, m // n) for a, m in cl if m % n == 0 and m >= n]
    n_prin = sum(k for _, k in prin)
    Vdec = ((n / 2) * sum((a - lbar) ** 2 for a, _ in exc)
            + n * sum(k * (a - lbar) ** 2 for a, k in prin))
    return {"p": p, "dimZ": n * (n - 6) // 8,
            "n_principal_constituents": n_prin,
            "expected_principal": (p * p - 9) // 8,
            "lambda_exc": float(exc[0][0]) if exc else None,
            "lambda_min": float(ev.min()), "lambda_bar": lbar,
            "V_spectral": V_spec, "V_decomposition": Vdec,
            "trPhi": float(np.trace(Phi)), "trPhi_exact": n * (n - 2),
            "proved": bool(abs(V_spec - Vdec) < 1e-6
                           and n_prin == (p * p - 9) // 8
                           and abs(np.trace(Phi) - n * (n - 2)) < 1e-6)}


def main():
    t0 = time.time()
    full = os.environ.get("PROP15593_FULL", "") == "1"
    out = {"prop": "15.593",
           "title": "Es4 = 4n^2 + tr(Phi^2); design floor; both leftovers reduce to V"}
    primes = (5, 7) if full else (5,)
    for p in primes:
        out[f"pointwise_p{p}"] = theorem_A_B_pointwise(p)["proved"]
        e = es4_and_V(p)
        out[f"es4_p{p}"] = {k: (str(v) if isinstance(v, Fraction) else v)
                            for k, v in e.items()}
        out[f"decomposition_p{p}"] = theorem_D_decomposition(p)
    # p=11 from the 15.592 exact value (census-free)
    n11 = 122
    Es4_11 = Fraction(json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15592.json").read_text()
    )["Es4_p11"]["Es4_exact"])
    V11 = Es4_11 - design_floor(11)
    x11 = Fraction(3260, 100)
    out["p11"] = {
        "Es4_exact": str(Es4_11), "V": str(V11), "V_per_n": float(V11 / n11),
        "floor_excess_per_n": float((design_floor(11) - 12 * n11 * n11) / n11),
        # meaningful metric: V as a fraction of the O(n) excess over 12n^2
        "pct_of_excess": float(100 * V11 / (Es4_11 - 12 * n11 * n11)),
        "threshold_leftover1_per_n": float(threshold_leftover1(11) / n11),
        "threshold_leftover3_per_n": float(threshold_leftover3(11, x11) / n11),
        "margin_leftover1": float(threshold_leftover1(11) / V11),
        "margin_leftover3": float(threshold_leftover3(11, x11) / V11),
        "leftover1_implies_leftover3": threshold_leftover3(11, x11) > threshold_leftover1(11),
    }
    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    (ROOT / "evidence" / "e1_gmin_m4_prop15593.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.593  Es4 = 4n^2 + tr(Phi^2); design floor; V is the shared core")
    for p in primes:
        e = out[f"es4_p{p}"]
        d = out[f"decomposition_p{p}"]
        print(f"  p={p}: identity={e['identity_Es4']} floor_ok={e['floor_agrees']} "
              f"V/n={e['V_per_n']:.4f}  decomposition={d['proved']} "
              f"({d['n_principal_constituents']} principal)")
    q = out["p11"]
    print(f"  p=11: V/n={q['V_per_n']:.4f}  (V is {q['pct_of_excess']:.1f}% of the O(n) excess); "
          f"thresholds c1={q['threshold_leftover1_per_n']:.2f} "
          f"c3={q['threshold_leftover3_per_n']:.2f}; "
          f"margins {q['margin_leftover1']:.1f}x / {q['margin_leftover3']:.1f}x")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()
