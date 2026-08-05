#!/usr/bin/env python3
"""
Residual / Hypothesis H close attempt — load-bearing algebra + GPU census.

Does NOT soft-close L. residual_closed_general is True only if a general-p
proof of H (or equivalent) is present.

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, Z as in 15.104.
E[(yᵀBy)²] = 6 + 8⟨m₄,κ_B⟩ for unit B∈Z (15.109).
ray(B) := (E[(yᵀBy)²] − 6)/2 = 4⟨m₄,κ_B⟩.
H(p) := 2(p+2)²/n = (p+2)²/d.
α(p) := 1 + 2/p² = 4 α_κ with α_κ=(p²+2)/(4p²) (15.111.B).

============================================================================
Theorem A (H form and H≤5) — PROVED Fraction.
  H(p)=2(p+2)²/n. For primes p≥3: H(p)≤5, equality only at p=3.
  Proof. 2(p+2)² ≤ 5(p²+1) ⇔ 3p²−8p−3≥0 for p≥3 (disc 64+36=100;
  root (8+10)/6=3).  ∎

Theorem B (H ⇒ 16N ⇒ bi-tight empty p≥5) — PROVED.
  ray≤H ≤5 ⇒ E[f²]≤6+2H≤16 ⇒ λ_max(Φ|Z)≤16 ⇒ 16N
  ⇒ bi-tight empty for p≥5 (15.61).  Also thr_ray−H=(3p+7)(p−5)/(2d)≥0
  for p≥5 so H⇒ dual gap G≽I (15.160).  ∎

Theorem C (Wick floor) — PROVED (15.111.B restated).
  ⟨κ/p², κ_B⟩ = α_κ ‖B‖² on Z, so the Wick contribution to ray is
  α(p)=1+2/p², constant for all unit B∈Z.  ∎

Theorem D (λ_min≥6) — PROVED.
  Q₄ = Beᵀ Gu_disj Be with Gu_disj a Gram compression (nonneg on the
  image of edge maps of the form Be_e=2B_ij).  Hence Q₄≥0, ray≥0,
  E[f²]≥6, λ_min(Φ|Z)≥6.  ∎

Theorem E (certified H at p=5,7) — CERTIFIED Max+ eigensystem (CPU/GPU).
  ray_max = H at p=5; ray_max < H at p=7; both ≤5; mult(λ_max)=d;
  G_min≥1.  ∎

Theorem F (residual status).
  H for all primes p≥5 remains OPEN (no Max+-free general proof).
  residual_closed_general = False.
  L OPEN until H (or orth≤room / E(1)) closes and Main chain ships.  ∎

Writes evidence/e1_residual_h_close.json
"""
from __future__ import annotations

import json
import sys
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
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def n_of(p: int) -> int:
    return p * p + 1


def d_of(p: int) -> int:
    return n_of(p) // 2


def H_of(p: int) -> Fraction:
    """Hypothesis H threshold: 2(p+2)²/n = (p+2)²/d."""
    return Fraction(2 * (p + 2) ** 2, n_of(p))


def alpha_wick(p: int) -> Fraction:
    """Wick floor of ray: 1 + 2/p²."""
    return Fraction(p * p + 2, p * p)


def mu_bar(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def thr_lambda(p: int) -> Fraction:
    """16(d−2)/d dual-gap design threshold."""
    d = d_of(p)
    return Fraction(16 * (d - 2), d)


def prove_theorem_A(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(3, 80) if is_prime(p)]
    ok = True
    rows = {}
    for p in primes:
        H = H_of(p)
        le5 = H <= 5
        eq3 = (p == 3 and H == 5) or (p > 3 and H < 5)
        # algebra 3p²−8p−3 ≥ 0 for p≥3
        poly = 3 * p * p - 8 * p - 3
        rows[str(p)] = {
            "H": str(H),
            "H_le_5": le5,
            "poly_3p2_minus_8p_minus_3": poly,
            "ok": le5 and poly >= 0,
        }
        if not (le5 and poly >= 0):
            ok = False
    return {
        "proved": ok and H_of(3) == 5,
        "theorem": "H=2(p+2)²/n ≤5 for primes p≥3 (eq only p=3).",
        "by_p": rows,
        "n_checked": len(primes),
    }


def prove_theorem_B() -> dict:
    return {
        "proved": True,
        "theorem": (
            "ray≤H≤5 ⇒ λ_max(Φ)≤16 ⇒ 16N ⇒ bi-tight empty for p≥5 (15.61). "
            "Also H⇒G≽I for p≥5 (15.160)."
        ),
    }


def prove_theorem_C() -> dict:
    return {
        "proved": True,
        "theorem": (
            "⟨κ/p²,κ_B⟩=α_κ‖B‖² on Z with α_κ=(p²+2)/(4p²) (15.111.B); "
            "Wick ray floor α=4α_κ=1+2/p², constant on unit Z."
        ),
        "alpha": {str(p): str(alpha_wick(p)) for p in (3, 5, 7, 11, 13)},
    }


def prove_theorem_D() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Q₄=BeᵀGu_disj Be is a Gram Rayleigh (Gu_disj from Max+ edge signs); "
            "Q₄≥0 ⇒ ray≥0 ⇒ λ_min(Φ|Z)≥6 for any conference Max+."
        ),
    }


def certify_live_H(p: int) -> dict:
    """Live Φ eigensystem: ray_max vs H, 16N, dual gap G≽I."""
    from minmax_quadratic import paley_conference_prime_power
    from e1_gmin_typeA_wedge import _zero_diag_nullspace, _A_from_vec

    path = Path(f"/tmp/maxplus_p{p}.npy")
    if not path.is_file():
        return {"p": p, "skipped": True, "reason": "Max+ cache missing"}
    Y = np.load(path).astype(np.float64)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 1e-8]
    null, mZ = _zero_diag_nullspace(Vp)
    U = Y @ Vp
    use_gpu = False
    try:
        import cupy as cp

        Ug = cp.asarray(U)
        Q = cp.zeros((N, mZ), dtype=cp.float64)
        Gf = np.zeros((mZ, mZ))
        for j in range(mZ):
            A = _A_from_vec(null[:, j], d)
            Q[:, j] = cp.einsum("na,ab,nb->n", Ug, cp.asarray(A), Ug)
            B = Vp @ A @ Vp.T
            for i in range(j + 1):
                Bi = Vp @ _A_from_vec(null[:, i], d) @ Vp.T
                Gf[i, j] = Gf[j, i] = float(np.sum(Bi * B))
        QtQ = cp.asnumpy(Q.T @ Q) / N
        use_gpu = True
    except Exception:
        Q = np.zeros((N, mZ))
        Gf = np.zeros((mZ, mZ))
        for j in range(mZ):
            A = _A_from_vec(null[:, j], d)
            Q[:, j] = np.einsum("na,ab,nb->n", U, A, U)
            B = Vp @ A @ Vp.T
            for i in range(j + 1):
                Bi = Vp @ _A_from_vec(null[:, i], d) @ Vp.T
                Gf[i, j] = Gf[j, i] = float(np.sum(Bi * B))
        QtQ = (Q.T @ Q) / N
    w, V = np.linalg.eigh(Gf)
    W12 = V @ np.diag(1.0 / np.sqrt(np.maximum(w, 1e-30))) @ V.T
    Phi = W12 @ QtQ @ W12
    evals = np.sort(np.linalg.eigvalsh(Phi))[::-1]
    lmax = float(evals[0])
    mult = int(np.sum(np.abs(evals - lmax) < 1e-6))
    ray = (lmax - 6.0) / 2.0
    H = float(H_of(p))
    thr = float(thr_lambda(p))
    Gmat = (d / 32.0) * (16.0 * np.eye(mZ) - Phi)
    gmin = float(np.min(np.linalg.eigvalsh(Gmat)))
    return {
        "p": p,
        "skipped": False,
        "N": int(N),
        "d": d,
        "dim_Z": mZ,
        "lambda_max": lmax,
        "mult_top": mult,
        "mult_eq_d": mult == d,
        "ray": ray,
        "H": H,
        "ray_le_H": ray <= H + 1e-8,
        "lambda_max_le_16": lmax <= 16.0 + 1e-9,
        "threshold": thr,
        "G_min": gmin,
        "G_ge_I": gmin >= 1.0 - 1e-8,
        "sixteen_N": lmax <= 16.0 + 1e-9,
        "backend": "CuPy GPU" if use_gpu else "CPU numpy",
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "H_for_all_p": False,
        "sixteen_N_for_all_p": False,
        "L_status": "OPEN",
        "open": (
            "Prove ray_max≤H(p)=2(p+2)²/n for all primes p≥5, OR the single "
            "16N inequality E[s⁴]≤Es4_*(p)=4(3p⁸−6p⁶−148p⁴−10p²+129)/(p⁴−8p²−1) "
            "(Props 15.164–15.165; mult≥d−1 proved). Equiv: ∑η²≤η_* with "
            "η_*=(p²−1)n(19p⁴−152p²−3)/(6p²(p⁴−8p²−1)). "
            "Census exact Es4: p=3 equality; p=5 ratio 0.996; p=7 ratio 0.965 "
            "(Gram/Φ, not single-root W — p=7 not 1-homogeneous). "
            "Dead this session: spherical Delsarte (→two-level 12097>9297 at p=5); "
            "BM(C) (maximizers ⊥ I,C,J,C²); C-sig/CR alone; Gershgorin frame. "
            "Preferred: Weil/Jacobi on Aut-orbit m4; Aut₀ V₊↪V_max Rayleigh; "
            "SOS on Q₄≤10N‖B‖². See evidence/RESIDUAL_ATTACK_2026-08-05_GPU_B.md."
        ),
    }


def main() -> dict:
    a = prove_theorem_A()
    b = prove_theorem_B()
    c = prove_theorem_C()
    d = prove_theorem_D()
    live = {str(p): certify_live_H(p) for p in (5, 7)}
    open_ = prove_open()
    cert_ok = all(
        (not r.get("skipped")) and r.get("ray_le_H") and r.get("sixteen_N")
        for r in live.values()
    )
    out = {
        "title": "Residual H algebra + GPU census p=5,7; general OPEN",
        "L_status": "OPEN",
        "proved": {
            "H_le_5": a["proved"],
            "H_implies_16N": b["proved"],
            "wick_floor_alpha": c["proved"],
            "lambda_min_ge_6": d["proved"],
            "H_certified_p5_p7": cert_ok,
            "H_for_all_p": False,
            "residual_for_all_p_ge_5": False,
        },
        "algebra": {
            "theorem_A": a,
            "theorem_B": b,
            "theorem_C": c,
            "theorem_D": d,
            "open": open_,
        },
        "live": live,
        "status": {
            "structure_ready_for_H": True,
            "general_residual": False,
            "sixteen_N_for_all_p": False,
        },
        "open_attack": open_["open"],
        "F3": "no soft-close",
        "F19": "no class_key thrash",
        "F20": "GPU used for Φ eigensystem when CuPy available",
    }
    path = ROOT / "evidence" / "e1_residual_h_close.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("residual H close module")
    print(f"  A H≤5: {a['proved']}")
    print(f"  B H⇒16N: {b['proved']}")
    print(f"  C wick floor: {c['proved']}")
    print(f"  D λ_min≥6: {d['proved']}")
    for p, r in live.items():
        if r.get("skipped"):
            print(f"  p={p} skipped: {r.get('reason')}")
        else:
            print(
                f"  p={p} ray≤H={r['ray_le_H']} 16N={r['sixteen_N']} "
                f"G≥I={r['G_ge_I']} backend={r['backend']}"
            )
    print(f"  residual_closed_general={open_['residual_closed_general']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
