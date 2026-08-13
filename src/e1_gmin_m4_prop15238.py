#!/usr/bin/env python3
"""
Prop 15.238 — Per-eigenrelations on |κ|=1; Per μ_part Cy-FE identity;
residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.229–231, 15.224–225)
  On |κ|=1: (p⁴−1)μ + 2φ = R̄₄ = (Per μ)_S − μ_S with
  (Per f)_S := ∑_T per(C[S,T]) f(T) (15.231).
  μ = μ_part + δ, μ_part = aκ + bφ, a=(p²−1)/den, b=−2/den,
  den=p²(p²−5), δ∈E_{±4p}.

PROVED Max+-free
  A. Coefficient identity (all primes p with p²≠5).
     Assume on every |κ|=1 four-set:
        (Per κ)_S = p⁴ κ_S − 6 φ_S,
        (Per φ)_S = (2n+1) φ_S,   n=p²+1.
     Then
        (Per μ_part)_S = p⁴ (μ_part)_S + 2 φ_S.
     Proof: Per μ_part = a Per κ + b Per φ
          = a(p⁴ κ − 6φ) + b(2n+1)φ
          = a p⁴ κ + (−6a + b(2n+1)) φ.
     Target: p⁴ μ_part + 2φ = a p⁴ κ + (p⁴ b + 2) φ.
     Equate φ-coeffs:
        −6a + b(2n+1) ≟ p⁴ b + 2.
     With den=p²(p²−5), a=(p²−1)/den, b=−2/den, n=p²+1:
        LHS num = −6(p²−1) − 2(2n+1) = −10p²,
        RHS num = −2p⁴ + 2 den = −2p⁴ + 2p²(p²−5) = −10p².
     Equal.  ∎

  B. Consequence for the δ-channel.
     On |κ|=1 the Cy-FE is (Per μ)_S = p⁴ μ_S + 2 φ_S.
     If (A) holds and μ=μ_part+δ, then equivalently
        (Per δ)_S = p⁴ δ_S    on every |κ|=1 four-set.
     (The particular solution alone matches the Cy functional equation on
     |κ|=1; the residual δ is a joint |κ|=1-eigenfunction of Per with
     eigenvalue p⁴.)  ∎

CERTIFIED (Paley; not general)
  C. (Per κ)_S = p⁴ κ_S − 6 φ_S and (Per φ)_S = (2n+1) φ_S on all |κ|=1
     four-sets of the Paley conference of order n=p²+1 for p∈{3,5}
     (full type census: 2 types at p=3, 4 types at p=5).  Pattern matches
     the conditional hypotheses of (A).  General-p / all-conference proof
     OPEN (character-sum / C² collision calculus for Per on κ,φ).

  D. At p=5 Max+ m₄: (Per δ)_S = p⁴ δ_S on |κ|=1 samples (δ=m₄−μ_part),
     consistent with (B).  On |κ|=3 the ratio (Per δ)/δ is not p⁴
     (Per is not p⁴·Id on the full δ vector).

OPEN (blocks residual i / predicate flip)
  E. Prove (C) for all primes p≥5 (Max+-free), **and/or** bound |δ| or
     |μ| on |κ|=1 (envelope, |μ|≤2/n, Aut-SOS, Path C, K₄ thr).
     Soft-close forbidden.

Until E: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15238.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402


def mu_part_coeffs(p: int) -> tuple[Fraction, Fraction]:
    den = p * p * (p * p - 5)
    return Fraction(p * p - 1, den), Fraction(-2, den)


def per_kappa_form(p: int, kappa: int, phi: int) -> int:
    """p⁴ κ − 6 φ (predicted on |κ|=1)."""
    return p**4 * kappa - 6 * phi


def per_phi_form(p: int, phi: int) -> int:
    """(2n+1) φ with n=p²+1."""
    n = n_of(p)
    return (2 * n + 1) * phi


def per_mupart_via_eigen(p: int, kappa: int, phi: int) -> Fraction:
    """a·(p⁴κ−6φ) + b·(2n+1)φ."""
    a, b = mu_part_coeffs(p)
    return a * per_kappa_form(p, kappa, phi) + b * per_phi_form(p, phi)


def p4_mupart_plus_2phi(p: int, kappa: int, phi: int) -> Fraction:
    a, b = mu_part_coeffs(p)
    return (p**4) * (a * kappa + b * phi) + 2 * phi


def theorem_per_mupart_identity() -> dict:
    """(A): eigenrelations ⇒ Per μ_part = p⁴ μ_part + 2φ on |κ|=1."""
    ok = True
    sample = {}
    # algebraic identity for all primes 5≤p≤300
    for p in range(5, 301):
        if not is_prime(p):
            continue
        a, b = mu_part_coeffs(p)
        n = n_of(p)
        lhs = -6 * a + b * (2 * n + 1)
        rhs = (p**4) * b + 2
        if lhs != rhs:
            ok = False
            break
        # spot-check scalar FE at sample (κ,φ)
        for kappa, phi in ((1, -2), (1, 2 * (p - 2)), (-1, 2), (-1, -2 * (p - 2))):
            if abs(kappa) != 1:
                continue
            if per_mupart_via_eigen(p, kappa, phi) != p4_mupart_plus_2phi(
                p, kappa, phi
            ):
                ok = False
                break
        if not ok:
            break
        if p in (5, 7, 11, 13, 17, 19, 23):
            sample[str(p)] = {
                "phi_coeff_match": str(lhs),
                "a": str(a),
                "b": str(b),
            }
    return {
        "proved": ok,
        "theorem": (
            "For every prime p≥5: if on |κ|=1 one has Per κ = p⁴κ−6φ and "
            "Per φ = (2n+1)φ (n=p²+1), then Per μ_part = p⁴ μ_part + 2φ.  "
            "Proof: equate coefficients using a=(p²−1)/den, b=−2/den; both "
            "φ-coefficients equal −10p²/den."
        ),
        "by_p_sample": sample,
        "depends_on": ["mu_part_coeffs", "n_of"],
    }


def theorem_delta_Per_eigenvalue() -> dict:
    """(B): structure — δ is Per-eigenvalue p⁴ on |κ|=1 under (A)+Cy-FE."""
    t = theorem_per_mupart_identity()
    return {
        "proved_structure": t["proved"],
        "bound_proved_general": False,
        "theorem": (
            "On |κ|=1: Cy-FE is Per μ = p⁴μ + 2φ.  Under the eigenrelations "
            "of (A), this is equivalent to Per δ = p⁴ δ on |κ|=1 for "
            "μ=μ_part+δ.  Structure only — does not bound |δ| or |μ|."
        ),
        "depends_on": [
            "theorem_per_mupart_identity",
            "prop_15_229_functional_eq",
            "prop_15_231_permanent_form",
        ],
    }


def theorem_per_eigen_certified_small_p() -> dict:
    """(C): certified at p=3,5 only — not general."""
    return {
        "proved": False,
        "certified_p": [3, 5],
        "forms": {
            "Per_kappa": "p^4 * kappa - 6 * phi",
            "Per_phi": "(2*n+1) * phi",
        },
        "theorem_open": (
            "For every prime p≥5 and every |κ|=1 four-set of the Paley "
            "conference of order n=p²+1: (Per κ)_S = p⁴ κ_S − 6 φ_S and "
            "(Per φ)_S = (2n+1) φ_S.  Certified by full (κ,φ)-type census at "
            "p=3 (2 types) and p=5 (4 types).  General proof OPEN."
        ),
        "depends_on": ["theorem_permanent_size4_form"],
    }


def residual_i_closed_via_238() -> bool:
    return False


def hinge_status_238() -> dict:
    return {
        "per_mupart_identity": theorem_per_mupart_identity()["proved"],
        "delta_Per_structure": theorem_delta_Per_eigenvalue()[
            "proved_structure"
        ],
        "per_eigen_general": theorem_per_eigen_certified_small_p()["proved"],
        "residual_i_closed_via_238": residual_i_closed_via_238(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "Prove Per κ / Per φ eigenrelations for all p≥5, then bound "
            "|μ| or |δ| on |κ|=1 (envelope / 2/n / Aut-SOS / Path C / K4)."
        ),
    }


def main() -> dict:
    a = theorem_per_mupart_identity()
    b = theorem_delta_Per_eigenvalue()
    c = theorem_per_eigen_certified_small_p()
    status = hinge_status_238()
    out = {
        "prop": "15.238",
        "title": (
            "Per-eigenrelations ⇒ Per μ_part Cy-FE; δ is Per-eigen p⁴ on "
            "|κ|=1; residual OPEN"
        ),
        "proved": {
            "per_mupart_identity": a["proved"],
            "delta_Per_structure": b["proved_structure"],
            "per_eigen_general": c["proved"],
            "residual_i_closed": residual_i_closed_via_238(),
        },
        "algebra": {
            "Per_kappa_form": "p^4*kappa - 6*phi",
            "Per_phi_form": "(2*n+1)*phi",
            "identity": "Per mu_part = p^4 mu_part + 2 phi on |kappa|=1",
            "thr_p5": str(mu4_threshold(5)),
            "maj_p5": str(majorant_mu_part(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_238(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Conditional eigenrelations close the particular Cy channel "
            "exactly.  Residual-(i) still needs a pointwise |μ| or |δ| bound. "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15238.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.238  residual_i={residual_i_closed_via_238()}")
    print(
        f"  mupart_id={a['proved']} eigen_general={c['proved']} "
        f"structure={b['proved_structure']}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
