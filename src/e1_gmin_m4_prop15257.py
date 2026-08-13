#!/usr/bin/env python3
"""
Prop 15.257 — Seidel 4×4 spectrum dichotomy by κ; resolvent of pI−C_S on
|κ|=1; structure toward Schur attack on m₄⁺=m₄⁻.  Residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP
  For a 4-set S, C_S is the principal 4×4 Seidel matrix (0 diagonal, ±1 off).
  κ(S)=C_ab C_cd + C_ac C_bd + C_ad C_bc ∈ {±1,±3} for any such matrix.
  On Max±: (Cy)=±p y ⇒ (pI−C_S)y_S = e (Max+) or (−pI−C_S)y_S = e (Max−)
  with e_i=∑_{j∉S} C_ij y_j the external field (Schur complement of S).

PROVED Max+-free (all 4×4 zero-diag ±1 symmetric matrices)
  A. Spectrum dichotomy.
        |κ|=1  ⇒  σ(C_S) = {±1, ±√5},   χ_{C_S}(λ) = λ⁴ − 6λ² + 5
                                        = (λ²−1)(λ²−5);
        |κ|=3  ⇒  σ(C_S) ∈ {{3,−1³}, {−3,1³}},
                  χ(λ)=λ⁴−6λ²∓8λ−3.
     Proof: 2⁶=64 exhaustive enumeration of edge-sign patterns; equivalently
     tr C_S²=12 forces the listed eigenvalue multisets once κ fixes the
     matching-product sum (classical Seidel 4-vertex classification).  ∎

  B. Resolvent on |κ|=1.
     det(pI−C_S)=(p²−1)(p²−5) for every prime p≥3 (and every |κ|=1 C_S).
     With A:=C_S and χ(λ)=λ⁴−6λ²+5,
        (pI−A)^{-1}
          = [p³ I + p² A + p(A²−6I) + (A³−6A)] / [(p²−1)(p²−5)].
     Proof: standard resolvent polynomial for monic χ with c₃=c₁=0,c₂=−6,c₀=5;
     verified on representatives.  Denominator matches μ_part / master denoms.  ∎

  C. Same det for Max− local: det(−pI−C_S)=det(pI+C_S)=χ(p)=(p²−1)(p²−5)
     (χ even).  ∎

STRUCTURE (toward hinge; not a bound)
  D. Schur form of 4-products.
     On Max+: y_S=(pI−C_S)^{-1} e with e the external field.
     On Max−: y_S=(−pI−C_S)^{-1} e′.
     Hence ∏_{i∈S} y_i is a multilinear form in the external field with
     resolvent coefficients.  On |κ|=1 the resolvent is the explicit
     cubic (B).  Difference m₄⁺−m₄⁻ is the average of these multilins
     under Max+ / Max− external laws.  OPEN: show the two averages agree
     on |κ|=1 (⇔ ν=0 ⇔ m₄⁺=m₄⁻).  ∎

  E. Why |κ|=1 is special.
     Minimal polynomial (A²−I)(A²−5I)=0 is even and excludes eigenvalues
     ±p for p≥5 (p≠√5).  On |κ|=3, ±3 can interact with p=3 (known
     exceptional); for p≥5, ±3≠±p still invertible, but the odd part of χ
     (the ±8λ term) is the algebraic distinction.  ∎

CERTIFIED
  F. (A)(B) on all edge patterns and on all |κ|=1 principal minors of Paley
     C at p=5,7.  ν=0 on |κ|=1 and ν≠0 on |κ|=3 at p=3,5,7 (prior census).

OPEN
  G. Close the Schur average identity of (D) Max+-free, **or** reflection,
     **or** Es4/ker=sc/hull.  Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15257.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    type_I_k_3p_minus_2_closed_general,
)


def kappa_of_matrix(M) -> int:
    return int(
        round(
            M[0, 1] * M[2, 3]
            + M[0, 2] * M[1, 3]
            + M[0, 3] * M[1, 2]
        )
    )


def theorem_seidel4_spectrum() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Every 4×4 zero-diag ±1 symmetric matrix has "
            "|κ|=1 ⇒ σ={±1,±√5}, χ=λ⁴−6λ²+5; "
            "|κ|=3 ⇒ σ∈{{3,−1³},{-3,1³}}.  Exhaustive 2⁶ / tr A²=12."
        ),
    }


def theorem_resolvent_kappa1() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1: det(pI−C_S)=(p²−1)(p²−5); "
            "(pI−A)^{-1}=[p³I+p²A+p(A²−6I)+(A³−6A)]/[(p²−1)(p²−5)].  "
            "Same det for −pI−A."
        ),
        "depends_on": ["theorem_seidel4_spectrum"],
    }


def det_pI_minus_CS(p: int) -> int:
    return (p * p - 1) * (p * p - 5)


def certify_seidel_and_resolvent() -> dict:
    import numpy as np

    spectra: dict[int, set] = {}
    resolvent_ok = True
    for bits in product([-1.0, 1.0], repeat=6):
        M = np.zeros((4, 4))
        idx = 0
        for i in range(4):
            for j in range(i + 1, 4):
                M[i, j] = M[j, i] = bits[idx]
                idx += 1
        k = kappa_of_matrix(M)
        ev = tuple(np.round(np.sort(np.linalg.eigvalsh(M)), 6))
        spectra.setdefault(k, set()).add(ev)
        if abs(abs(k) - 1) < 1e-9:
            A = M
            for p in (5, 7, 11, 13):
                R = np.linalg.inv(p * np.eye(4) - A)
                num = (
                    p**3 * np.eye(4)
                    + p**2 * A
                    + p * (A @ A - 6 * np.eye(4))
                    + (A @ A @ A - 6 * A)
                )
                chi_p = float(det_pI_minus_CS(p))
                if np.max(np.abs(R - num / chi_p)) > 1e-9:
                    resolvent_ok = False
                if abs(np.linalg.det(p * np.eye(4) - A) - chi_p) > 1e-6:
                    resolvent_ok = False

    # Expected spectra
    expected_k1 = {(-2.236068, -1.0, 1.0, 2.236068)}
    # round check
    k1_ok = all(
        abs(ev[0] + 2.236068) < 1e-5
        and abs(ev[1] + 1) < 1e-9
        and abs(ev[2] - 1) < 1e-9
        and abs(ev[3] - 2.236068) < 1e-5
        for ev in spectra.get(1, set()) | spectra.get(-1, set())
    )
    k3_ok = spectra.get(3, set()) | spectra.get(-3, set())
    k3_expected = {
        (-3.0, 1.0, 1.0, 1.0),
        (-1.0, -1.0, -1.0, 3.0),
    }
    # normalize rounding
    k3_norm = set()
    for ev in k3_ok:
        k3_norm.add(tuple(np.round(ev, 6)))

    return {
        "n_patterns": 64,
        "spectra_by_kappa": {str(k): [list(ev) for ev in sorted(v)] for k, v in spectra.items()},
        "kappa1_spectrum_ok": k1_ok,
        "kappa3_spectra": [list(ev) for ev in sorted(k3_norm)],
        "kappa3_ok": k3_norm == k3_expected,
        "resolvent_ok": resolvent_ok,
        "det_formula": {
            str(p): det_pI_minus_CS(p) for p in (5, 7, 11, 13)
        },
    }


def residual_i_closed_via_257() -> bool:
    return False


def hinge_status_257() -> dict:
    return {
        "seidel4_spectrum": True,
        "resolvent_kappa1": True,
        "schur_average_identity": False,
        "residual_i_closed": residual_i_closed_via_257(),
        "open": (
            "Close Schur multilin average identity m₄⁺=m₄⁻ on |κ|=1, "
            "or reflection |ρ|≤|ρ_f4|, or Es4/ker=sc/hull."
        ),
    }


def main() -> dict:
    a = theorem_seidel4_spectrum()
    b = theorem_resolvent_kappa1()
    cert = certify_seidel_and_resolvent()
    out = {
        "title": (
            "Prop 15.257 Seidel 4×4 spectrum / resolvent on |κ|=1; "
            "Schur structure; residual (i) OPEN"
        ),
        "proved": {
            "seidel4_spectrum": a["proved"],
            "resolvent_kappa1": b["proved"],
            "schur_average_identity": False,
            "residual_i_closed": residual_i_closed_via_257(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"spectrum": a, "resolvent": b},
        "certified": cert,
        "hinge": hinge_status_257(),
        "F3": (
            "Seidel 4×4: |κ|=1⇒{±1,±√5}, resolvent (pI−A)^{-1} with "
            "det=(p²−1)(p²−5) proved.  Schur form of ∏y_S set up for "
            "m₄⁺=m₄⁻ attack.  Average identity OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15257.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.257 Seidel spectrum / resolvent; residual-i OPEN")
    print(f"  spectrum: {a['proved']}")
    print(f"  resolvent: {b['proved']}")
    print(
        f"  cert: k1={cert['kappa1_spectrum_ok']} k3={cert['kappa3_ok']} "
        f"resolvent={cert['resolvent_ok']}"
    )
    print(f"  residual_i closed: {residual_i_closed_via_257()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
