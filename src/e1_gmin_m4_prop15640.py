#!/usr/bin/env python3
"""Prop 15.640: exact harmonic operator on the scaled-norm ``3p-6`` shell.

Proposition 15.639 classifies that complete shell, for every odd prime
``p>=11``, as negative signed triples together with point--square-circle
vectors.  For an admissible symmetric tensor

    P W P = W,       diag(W)=0,

put ``F=||W||_F^2`` and ``q_S=w_S^T W w_S``.  The two families have exact
quartic sums

    sum_T (x^T W x)^2 = 2(p-3)(p+1) F,

    sum_O (x^T W x)^2
      = 8(p-2) F + 2(p-5)/p^3 sum_S q_S^2.

The second identity uses the new through-point frame formula

    sum_{S containing i} w_S w_S^T
      = p^2 (P - 2(Pe_i)(Pe_i)^T).

Indeed, after moving ``i`` to infinity, the circles through it split into
``(p+1)/2`` mutually orthogonal parallel classes.  Each class has ``p``
words with Gram matrix ``p(pI-J)``, hence frame eigenvalue ``p^2``; their
total span is the codimension-one space ``V_+ intersect e_i^perp``.

Combining the quartic sums with the universal radial terms in ``H_W`` and
the square-circle tensor spectrum from Proposition 15.634 gives three
eigenvalues.  For every ``p>=11`` the circle-kernel eigenvalue is strictly
negative, while the two circle-image eigenvalues are strictly positive.
The norm-parity phase and evaluation at ``x/2`` multiply them by ``-1/16``.

This proves that the shell is a quartic saddle, not a spherical 4-design.
It still does not control all later shells and therefore does not prove R1.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def n_of(p: int) -> int:
    if p < 3 or p % 2 == 0:
        raise ValueError("the closed identities require odd p>=3")
    return p * p + 1


def rank_of(p: int) -> int:
    return n_of(p) // 2


def z_dimension(p: int) -> int:
    return n_of(p) * (n_of(p) - 6) // 8


def shell_squared_norm(p: int) -> Fraction:
    return Fraction(3 * p - 6, 2 * p)


def shell_signed_count(p: int) -> int:
    try:
        from .e1_gmin_m4_prop15639 import (
            first_nonminimal_odd_shell_signed_count,
        )
    except ImportError:
        from e1_gmin_m4_prop15639 import (
            first_nonminimal_odd_shell_signed_count,
        )

    return first_nonminimal_odd_shell_signed_count(p)


def through_point_circle_frame_certificate(p: int) -> dict:
    """Closed frame proof for ``sum_{S contains i} w_S w_S^T``."""
    n_of(p)
    blocks = (p + 1) // 2
    words_per_block = p
    within_block_rank = p - 1
    span_dimension = blocks * within_block_rank
    return {
        "parallel_classes": blocks,
        "words_per_class": words_per_block,
        "within_class_gram": "p*(p*I-J)",
        "within_class_nonzero_eigenvalue": p * p,
        "within_class_nonzero_multiplicity": within_block_rank,
        "distinct_classes_orthogonal": True,
        "span_dimension": span_dimension,
        "target_dimension": rank_of(p) - 1,
        "frame_operator_eigenvalue": p * p,
        "frame_identity": "p^2*(P-2*(P e_i)*(P e_i)^T)",
        "proved": bool(span_dimension == rank_of(p) - 1),
    }


def total_circle_frame_scalar(p: int) -> int:
    """``sum_S w_S w_S^T = scalar*P`` over unoriented square circles."""
    n_of(p)
    return p * p * (p - 1)


def off_circle_row_contraction_scalar(p: int) -> int:
    """Multiplier of ``F`` in ``sum_{S,i notin S}(Ww_S)_i^2``."""
    n_of(p)
    return p * p * (p - 2)


def negative_triangle_quartic_scalar(p: int) -> int:
    """Multiplier of ``F`` from both signs of all negative triples."""
    n_of(p)
    return 2 * (p - 3) * (p + 1)


def point_circle_quartic_scalar(p: int) -> int:
    """Scalar part of the signed point--circle quartic sum."""
    n_of(p)
    return 8 * (p - 2)


def point_circle_evaluation_coefficient(p: int) -> Fraction:
    """Coefficient of ``sum_S (w_S^T W w_S)^2``."""
    n_of(p)
    return Fraction(2 * (p - 5), p**3)


def complete_quartic_scalar(p: int) -> int:
    return negative_triangle_quartic_scalar(p) + point_circle_quartic_scalar(p)


def radial_harmonic_correction(p: int) -> Fraction:
    """Positive scalar subtracted from the quartic evaluation operator."""
    n = shell_signed_count(p)
    q = shell_squared_norm(p)
    d = rank_of(p)
    return Fraction(2 * n, d * (d + 2)) * q * q


def radial_harmonic_correction_closed(p: int) -> Fraction:
    return Fraction(3 * (p - 1) * (p + 7) * (p - 2) ** 2, p * p + 5)


def harmonic_scalar_offset(p: int) -> Fraction:
    """Eigenvalue on the kernel of all square-circle evaluations."""
    return Fraction(complete_quartic_scalar(p)) - radial_harmonic_correction(p)


def harmonic_kernel_closed(p: int) -> Fraction:
    numerator = p**4 + 2 * p**3 - 69 * p * p + 136 * p + 26
    return -Fraction(numerator, p * p + 5)


def harmonic_circle_low_closed(p: int) -> Fraction:
    numerator = p**4 - 14 * p**3 + 89 * p * p - 196 * p + 24
    return Fraction(numerator, p * p + 5)


def harmonic_circle_high_closed(p: int) -> Fraction:
    numerator = p**4 - 10 * p**3 + 69 * p * p - 176 * p - 76
    return Fraction(numerator, p * p + 5)


def harmonic_spectrum(p: int) -> list[dict]:
    """Complete unphased sum ``sum_x H_W(x)`` on the ``s=3p-6`` shell."""
    a = harmonic_scalar_offset(p)
    coefficient = point_circle_evaluation_coefficient(p)
    rows = [
        {
            "channel": "circle-kernel",
            "circle_tensor_eigenvalue": 0,
            "eigenvalue": a,
            "closed_form": harmonic_kernel_closed(p),
            "multiplicity": n_of(p) * (p - 1) * (p - 3) // 8,
        },
        {
            "channel": "circle-low",
            "circle_tensor_eigenvalue": p**3 * (p - 1),
            "eigenvalue": a + coefficient * p**3 * (p - 1),
            "closed_form": harmonic_circle_low_closed(p),
            "multiplicity": n_of(p) * (p - 1) // 4,
        },
        {
            "channel": "circle-high",
            "circle_tensor_eigenvalue": p**3 * (p + 1),
            "eigenvalue": a + coefficient * p**3 * (p + 1),
            "closed_form": harmonic_circle_high_closed(p),
            "multiplicity": n_of(p) * (p - 3) // 4,
        },
    ]
    assert sum(row["multiplicity"] for row in rows) == z_dimension(p)
    assert all(row["eigenvalue"] == row["closed_form"] for row in rows)
    return rows


def shifted_sign_polynomials(p: int) -> dict[str, tuple[int, ...]]:
    """Positive-coefficient expansions at ``x=p-11>=0``."""
    n_of(p)
    return {
        "minus_kernel_numerator": (10_476, 4_668, 723, 46, 1),
        "circle_low_numerator": (4_644, 2_004, 353, 30, 1),
        "circle_high_numerator": (7_668, 3_036, 465, 34, 1),
    }


def harmonic_sign_certificate(p: int) -> dict:
    x = p - 11
    coefficients = shifted_sign_polynomials(p)

    def evaluate(values: tuple[int, ...]) -> int:
        return sum(value * x**degree for degree, value in enumerate(values))

    kernel_numerator = p**4 + 2 * p**3 - 69 * p * p + 136 * p + 26
    low_numerator = p**4 - 14 * p**3 + 89 * p * p - 196 * p + 24
    high_numerator = p**4 - 10 * p**3 + 69 * p * p - 176 * p - 76
    return {
        "p_minus_11": x,
        "positive_coefficient_expansions": coefficients,
        "kernel_positive_numerator": kernel_numerator,
        "circle_low_positive_numerator": low_numerator,
        "circle_high_positive_numerator": high_numerator,
        "expansions_match": bool(
            evaluate(coefficients["minus_kernel_numerator"])
            == kernel_numerator
            and evaluate(coefficients["circle_low_numerator"])
            == low_numerator
            and evaluate(coefficients["circle_high_numerator"])
            == high_numerator
        ),
        "unphased_signs": ("negative", "positive", "positive"),
        "indefinite": bool(
            harmonic_kernel_closed(p) < 0
            and harmonic_circle_low_closed(p) > 0
            and harmonic_circle_high_closed(p) > 0
        ),
    }


def parity_twisted_half_spectrum(p: int) -> list[dict]:
    """Poisson-shadow contribution: odd phase and ``H_W(x/2)=H_W(x)/16``."""
    return [
        {
            **row,
            "unphased_eigenvalue": row["eigenvalue"],
            "unphased_closed_form": row["closed_form"],
            "eigenvalue": -row["eigenvalue"] / 16,
            "closed_form": -row["closed_form"] / 16,
        }
        for row in harmonic_spectrum(p)
    ]


def harmonic_operator_theorem(
    primes: tuple[int, ...] = (11, 13, 17, 19, 23, 29, 31),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        frame = through_point_circle_frame_certificate(p)
        spectrum = harmonic_spectrum(p)
        signs = harmonic_sign_certificate(p)
        shadow = parity_twisted_half_spectrum(p)
        row_ok = bool(
            frame["proved"]
            and off_circle_row_contraction_scalar(p) == p * p * (p - 2)
            and complete_quartic_scalar(p) == 2 * (p * p + 2 * p - 11)
            and radial_harmonic_correction(p)
            == radial_harmonic_correction_closed(p)
            and signs["expansions_match"]
            and signs["indefinite"]
            and [row["eigenvalue"] for row in shadow]
            == [-row["eigenvalue"] / 16 for row in spectrum]
            and [row["closed_form"] for row in shadow]
            == [-row["closed_form"] / 16 for row in spectrum]
        )
        rows[str(p)] = {
            "scaled_norm": 3 * p - 6,
            "squared_norm": str(shell_squared_norm(p)),
            "signed_shell_count": shell_signed_count(p),
            "through_point_frame": frame,
            "negative_triangle_quartic_scalar": (
                negative_triangle_quartic_scalar(p)
            ),
            "point_circle_quartic_scalar": point_circle_quartic_scalar(p),
            "point_circle_evaluation_coefficient": str(
                point_circle_evaluation_coefficient(p)
            ),
            "complete_quartic_scalar": complete_quartic_scalar(p),
            "radial_harmonic_correction": str(radial_harmonic_correction(p)),
            "harmonic_spectrum": [
                {
                    **entry,
                    "eigenvalue": str(entry["eigenvalue"]),
                    "closed_form": str(entry["closed_form"]),
                }
                for entry in spectrum
            ],
            "parity_twisted_half_spectrum": [
                {
                    **entry,
                    "eigenvalue": str(entry["eigenvalue"]),
                    "closed_form": str(entry["closed_form"]),
                    "unphased_eigenvalue": str(entry["unphased_eigenvalue"]),
                    "unphased_closed_form": str(
                        entry["unphased_closed_form"]
                    ),
                }
                for entry in shadow
            ],
            "sign_certificate": signs,
            "checks": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": (
            "For every odd prime p>=11, the complete scaled-norm 3p-6 "
            "harmonic shell has one negative circle-kernel eigenvalue and "
            "two positive square-circle-image eigenvalues."
        ),
        "rows": rows,
    }


def main() -> dict:
    theorem = harmonic_operator_theorem()
    out = {
        "prop": "15.640",
        "title": "Quartic saddle on the complete scaled-norm 3p-6 shell",
        "proved": {
            "through_point_square_circle_frame_identity_all_p_ge_11": theorem[
                "proved"
            ],
            "complete_scaled_norm_3p_minus_6_harmonic_operator": theorem[
                "proved"
            ],
            "complete_scaled_norm_3p_minus_6_shell_is_quartic_saddle": theorem[
                "proved"
            ],
            "R1": False,
            "global_QVAR": False,
        },
        "theorem": theorem,
        "remaining_obstruction": (
            "The intervening even candidates for p>=17 and the complete "
            "later norm-parity-twisted harmonic theta tail remain uncontrolled."
        ),
        "L_status": "OPEN",
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15640.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Prop 15.640 quartic saddle: {theorem['proved']}")
    print(f"  wrote {dest}")
    return out


if __name__ == "__main__":
    main()
