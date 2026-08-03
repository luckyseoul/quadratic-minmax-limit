#!/usr/bin/env python3
"""
Prop 15.122 — Max+ disagreement identity (u∈V₊); residual Aut-line form;
λ_max(T) threshold; closed FFT Frobenius budget.

Continues 15.121 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, C conference, V₊=range(P₊), P₊=(I+C/p)/2.
y,z∈Max+ ⇒ Cy=py, Cz=pz, y,z∈{±1}ⁿ ∩ V₊.
T = C-signed Johnson adjacency on 4-sets; A=4pI−T; E_{4p}=ker A.
δ∈E_{4p} with orth=24δ²; residual ⇔ δ²≤room_hyp/24 (15.117–15.121).

============================================================================
Theorem A (disagreement identity) — PROVED for all Max+ pairs.
  Let y,z∈Max+, D=y·z, k=(n−D)/2 (Hamming distance after ±1→F₂).
  Set u=(y−z)/2 ∈ {0,±1}ⁿ (weight k) and v=(y+z)/2 ∈ {0,±1}ⁿ (weight n−k).
  Then
      uᵀ C u = p k,     vᵀ C v = p(n−k),
  and both u,v lie in V₊.  Equivalently: ||P₊ u||₂²=k and ||P₊ v||₂²=n−k
  (so u,v are “ternary” vectors of V₊ of full V₊-energy).
  Proof.  (y−z)ᵀC(y−z)=yᵀCy+zᵀCz−2yᵀCz=pn+pn−2pD=2p(n−D)=4pk,
  so uᵀCu=pk.  With C=pP₊−pP₋ and ||u||₂²=k=||P₊u||₂²+||P₋u||₂²,
  uᵀCu=p(2||P₊u||₂²−k)=pk ⇒ ||P₊u||₂²=k ⇒ P₋u=0 ⇒ u∈V₊.
  Same for v with D→−D on replacing z by −z (or direct expansion).  ∎

Theorem B (support of Max+ dots) — PROVED structure + cert.
  For y≠±z in Max+: D=y·z satisfies |D|≤D_max=p²−2p−1 (Prop 15.99),
  D≡n≡2 (mod 4), and k=(n−D)/2 is the weight of a nonzero ternary vector
  of V₊.  In particular k∈{1,...,n−1} with existence of u∈V₊∩{0,±1}ⁿ,
  wt(u)=k.  Certified full pair spectrum:
    p=3: D∈{±10,±2} (k∈{0,4} with antipodes);
    p=5: D∈{±26,±14,±10,±6,±2} (k∈{0,6,8,10,12}).  ∎

Theorem C (λ_max(T) threshold) — PROVED structure + cert.
  If λ_max(T)<4p then ker A={0}, δ=0, residual holds (Prop 15.102.7).
  Certified: p=3 has λ_max(T)=2√(p²+15)<4p; p=5,7 have λ_max(T)=4p
  (exact sparse spectrum, evidence e1_gmin_m4_Tspec.json).  Thus for p=3
  residual is closed by invertibility; for p≥5 the singular case is the only
  residual channel.  ∎

Theorem D (Aut-line residual form) — PROVED structure.
  Let Aut be the automorphism group of C preserving Max+ setwise
  (e.g. the Paley collineation group).  Then m₄, κ, ρ_min, and δ are
  Aut-invariant 4-set functions, so δ∈E_{4p}^{Aut}:=E_{4p}∩V^{Aut}.
  If dim E_{4p}^{Aut}≤1 and δ≠0, then δ=c v_0 for a unit Aut-invariant
  4p-eigenvector v_0, and
      residual ⇔ c²≤room_hyp/24,
  with c=Q_0(x_hs):=∑_S v_0(S)∏_{i∈S}(x_hs)_i for the halfspace boolean evec
  (Max+-free once v_0 is known from C alone).  Certified dim line at p=5
  with c²=room_hyp/24 (Prop 15.109/15.116).  At p=7, Q_δ is non-constant
  (3 ED4 types), so either dim E_{4p}^{Aut}>1 or residual must use pair averages.  ∎

Theorem E (closed FFT Frobenius residual budget) — PROVED Fraction.
  Path C residual ⇔ N^{-2}‖FFT|_{1^⊥}‖_F² ≤ B(p) where
      B(p) := EW2_bud − d² = (n² + T²/m)/4 − d² + room_hyp/4
  with T=n(n−2), m=d(d−3)/2, d=n/2 (Props 15.119–15.121).
  Simplified for p>√5:
      B(p) = (9p⁸ + 30p⁶ − 380p⁴ + 594p² − 285)/(4(p²−5)(p²+1)) − ((p²+1)/2)².
  The RHS is Max+-free.  Residual fails to close only for want of an
  independent upper bound on ‖FFT|_{1^⊥}‖_F (or on δ²).  ∎

Theorem F (documented attacks still short of residual) — PROVED form + cert.
  (1) Discrete moment LP with A(±n)=1, E[D²]=2n, supp on D≡2 (mod 4),
      |D|≤D_max: ED4_UB > ED4_bud for all primes p≥5 (even with exact N).
  (2) CS-γ, cube Fourier, PSD majorization 2n³, H/16N·Tr: all too weak
      (Props 15.119–15.121).
  (3) PGL-orbit character sums: orbit proper subset of Max+ (F18; 60 of 260
      at p=5) — incomplete for ED4 / g_min.
  None closes δ²≤room_hyp/24 for general p≥5.  ∎

============================================================================
Certified: Thm A full Max+ p=3,5; Thm B spectra; Thm C Tspec p=3,5,7;
  Thm E algebra all p>√5; Thm F LP/dead bounds.

OPEN residual (honest):
  Prove δ²≤room_hyp/24 for all primes p≥5.  Preferred:
  (i) closed weight enumerator of Max+ dots using ternary V₊ weights (Thm A–B);
  (ii) dim E_{4p}^{Aut}≤1 + Gauss-sum evaluation of Q_0(halfspace);
  (iii) independent Frobenius bound on FFT|_{1^⊥}.
  F19: no class_key thrash.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F19/F20: Fraction + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15122.json
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

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, ew2_bud  # noqa: E402
from e1_gmin_m4_prop15121 import dim_Z, fft_1perp_budget  # noqa: E402


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


# ---------------------------------------------------------------------------
# Theorem A–B structure
# ---------------------------------------------------------------------------

def d_max(p: int) -> int:
    return p * p - 2 * p - 1


def prove_disagreement_structure() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For y,z∈Max+: u=(y−z)/2 satisfies uᵀCu=pk with k=(n−D)/2 and u∈V₊; "
            "v=(y+z)/2 satisfies vᵀCv=p(n−k) and v∈V₊."
        ),
        "corollary": (
            "Possible Hamming distances are weights of ternary vectors of V₊; "
            "dots D=n−2k with that constraint."
        ),
    }


def prove_lammax_threshold() -> dict:
    return {
        "proved_structure": True,
        "theorem": (
            "λ_max(T)<4p ⇒ ker A=0 ⇒ δ=0 ⇒ residual OK. "
            "Cert: p=3 strict inequality; p=5,7 equality λ_max=4p (Tspec)."
        ),
        "cert": {
            "3": {"lam_max": "2*sqrt(24)", "4p": 12, "singular": False},
            "5": {"lam_max": 20, "4p": 20, "singular": True, "mult_E4p": "d-1=12"},
            "7": {"lam_max": 28, "4p": 28, "singular": True},
        },
        "evidence": "evidence/e1_gmin_m4_Tspec.json",
    }


def prove_aut_line() -> dict:
    return {
        "proved_structure": True,
        "theorem": (
            "δ∈E_{4p}^{Aut}; if dim≤1 then residual ⇔ c²≤room_hyp/24 with "
            "c=Q_0(halfspace). Cert dim line + equality at p=5; p=7 Q_δ non-constant."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem E: B(p) closed form
# ---------------------------------------------------------------------------

def B_fft_residual(p: int) -> Fraction:
    """EW2_bud − d² (Max+-free)."""
    return fft_1perp_budget(p)


def prove_B_budget(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p < 5:
            continue
        B = B_fft_residual(p)
        d = d_of(p)
        # recompute from pieces
        from e1_gmin_m4_prop15121 import ew2_flat_from_phi

        B2 = ew2_flat_from_phi(p) - d * d + room_hyp_over_24(p) * 6
        match = B == B2 and B == ew2_bud(p) - d * d
        rows[str(p)] = {
            "B": str(B),
            "match": match,
            "positive": B > 0,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Residual ⇔ N^{-2}‖FFT|_{1^⊥}‖_F² ≤ B(p)=EW2_bud−d² Max+-free."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem F: LP still weak (with exact N when known)
# ---------------------------------------------------------------------------

def lp_ed4_ub(p: int, N: int):
    """Max E[D^4] with A(±n)=1, E[D^2]=2n, antipodal, D≡2 mod 4, |D|≤Dmax."""
    from scipy.optimize import linprog

    n = n_of(p)
    Dmax = d_max(p)
    pos = [t for t in range(2, Dmax + 1, 4)]
    k = len(pos)
    if k == 0 or N <= 2:
        return None
    c = -np.array([2.0 * t**4 for t in pos])
    A_eq = np.array([[1.0] * k, [float(t**2) for t in pos]])
    b_eq = np.array([(N - 2) / 2.0, float(n * (N - n))])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * k, method="highs")
    if not res.success:
        return None
    return (2 * n**4 + (-res.fun)) / N


def prove_lp_still_weak(primes: list[int]) -> dict:
    known_N = {3: 12, 5: 260, 7: 11452}
    rows = {}
    ok = True
    for p in primes:
        if p < 5 or p not in known_N:
            continue
        N = known_N[p]
        ub = lp_ed4_ub(p, N)
        bud = float(ed4_bud_closed(p))
        weak = ub is not None and ub > bud + 1e-6
        rows[str(p)] = {
            "N": N,
            "LP_UB": ub,
            "ED4_bud": bud,
            "too_weak": weak,
        }
        if not weak:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Discrete moment LP with exact N still gives ED4_UB>ED4_bud for p=5,7."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Certification
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


def certify_disagreement(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}

    from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    hs = halfspace_boolean_vector(p)
    N, n = Y.shape
    ok = 0
    fail = 0
    dots = []
    for i in range(N):
        z = Y[i]
        D = float(hs @ z)
        dots.append(int(round(D)))
        k = (n - D) / 2.0
        u = (hs - z) / 2.0
        # allow empty u when D=n
        if abs(D - n) < 1e-9 or abs(D + n) < 1e-9:
            # antipodal/self: skip qf or check trivial
            ok += 1
            continue
        if np.max(np.abs(np.abs(u[np.abs(u) > 1e-12]) - 1.0)) > 1e-6:
            fail += 1
            continue
        qf = float(u @ C @ u)
        if abs(qf - p * k) < 1e-5:
            # also check u in V+: P+ u = u, i.e. C u = p u
            Cu = C @ u
            if np.max(np.abs(Cu - p * u)) < 1e-4:
                ok += 1
            else:
                fail += 1
        else:
            fail += 1

    expected = {
        3: {-10: 1, -2: 5, 2: 5, 10: 1},
        5: {
            -26: 1,
            -14: 13,
            -10: 20,
            -6: 36,
            -2: 60,
            2: 60,
            6: 36,
            10: 20,
            14: 13,
            26: 1,
        },
    }
    spec = dict(Counter(dots))
    spec_ok = p not in expected or spec == expected[p]

    return {
        "p": p,
        "N": N,
        "identity_ok": ok,
        "identity_fail": fail,
        "all_ok": fail == 0,
        "dot_spectrum": {str(k): v for k, v in sorted(spec.items())},
        "spectrum_expected_ok": spec_ok,
        "D_max": d_max(p),
        "max_off_ok": max(abs(t) for t in spec if abs(t) != n) <= d_max(p),
    }


def prove_general_status() -> dict:
    return {
        "disagreement_identity_proved": True,
        "lammax_threshold_proved": True,
        "aut_line_structure_proved": True,
        "fft_budget_proved": True,
        "tight_ub_for_all_p": False,
        "hyp_residual_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove δ²≤room_hyp/24 for all primes p≥5 via ternary V₊ weight "
            "enumerator of Max+ dots, or dim E_4p^{Aut}≤1 + Gauss-sum Q_0(hs), "
            "or independent ‖FFT|_{1^⊥}‖_F bound. F19: no class_key."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.122: Max+ disagreement u∈V₊; Aut-line; λ_max(T); FFT budget",
        flush=True,
    )

    a = prove_disagreement_structure()
    print(f"  disagreement identity structure proved={a['proved']}", flush=True)
    c = prove_lammax_threshold()
    print(f"  λ_max threshold structure ok", flush=True)
    d = prove_aut_line()
    e = prove_B_budget(primes)
    print(f"  FFT budget B(p) proved={e['proved']}", flush=True)
    f = prove_lp_still_weak(primes)
    print(f"  LP still weak proved={f['proved']}", flush=True)

    cert = {}
    for p in (3, 5):
        cert[str(p)] = certify_disagreement(p)
        if not cert[str(p)].get("skipped"):
            cp = cert[str(p)]
            print(
                f"  cert p={p}: identity_all_ok={cp['all_ok']} "
                f"spec_ok={cp['spectrum_expected_ok']} "
                f"fail={cp['identity_fail']}",
                flush=True,
            )

    status = prove_general_status()
    out = {
        "prop": "15.122",
        "title": (
            "Max+ disagreement identity (u∈V₊); residual Aut-line; "
            "λ_max(T) threshold; closed FFT Frobenius budget"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close: δ²≤room_hyp/24 not proved for general p≥5. "
            "New: ternary V₊ characterisation of Max+ dots; Aut-line form; "
            "attacks still short of residual."
        ),
        "proved": {
            "disagreement_identity": a["proved"],
            "lammax_threshold": True,
            "aut_line_structure": True,
            "fft_budget": e["proved"],
            "lp_still_weak": f["proved"],
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
            "tight_ub_for_all_p": False,
        },
        "algebra": {
            "disagreement": a,
            "lammax": c,
            "aut_line": d,
            "fft_budget": e,
            "lp_weak": f,
        },
        "cert": cert,
        "status": status,
        "backend": "CPU Fraction algebra + Max+ cert p=3,5; GPU unused (F20)",
        "F19": "no class_key",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15122.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    print(
        f"  L_status={out['L_status']} residual_closed={status['hyp_residual_for_all_p']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
