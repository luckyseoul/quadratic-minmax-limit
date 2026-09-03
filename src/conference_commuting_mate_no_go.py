#!/usr/bin/env python3
r"""Exact obstruction to a commuting skew-signing mate.

Let ``A`` be symmetric with zero diagonal and signs off the diagonal, and
let ``R`` be skew-symmetric with the same support.  In even order,

    (AR-RA)_ii = -2 sum_(j != i) A_ij R_ij = 2 (mod 4).

Thus ``AR != RA``.  In particular, a real symmetric conference matrix can
never commute with any skew signing, whether or not the latter is itself a
conference matrix.

The module records three exact checks used by the accompanying evidence
note: the diagonal parity identity, the mod-two invertibility of every
even-order skew signing, and an exhaustive order-six Paley audit.  None of
these finite checks is used as a substitute for the displayed all-orders
proof.
"""
from __future__ import annotations

from itertools import combinations
import json

import numpy as np


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _is_prime(q: int) -> bool:
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    d = 3
    while d * d <= q:
        if q % d == 0:
            return False
        d += 2
    return True


def _validate_pair(A: np.ndarray, R: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A = np.asarray(A, dtype=np.int64)
    R = np.asarray(R, dtype=np.int64)
    if A.ndim != 2 or A.shape[0] != A.shape[1] or R.shape != A.shape:
        raise ValueError("A and R must be square matrices of the same order")
    n = A.shape[0]
    off_diagonal = ~np.eye(n, dtype=bool)
    if not np.array_equal(A, A.T) or np.any(np.diag(A)):
        raise ValueError("A must be symmetric with zero diagonal")
    if not np.array_equal(R, -R.T) or np.any(np.diag(R)):
        raise ValueError("R must be skew-symmetric with zero diagonal")
    if not np.all(np.abs(A[off_diagonal]) == 1):
        raise ValueError("A must have signs off the diagonal")
    if not np.all(np.abs(R[off_diagonal]) == 1):
        raise ValueError("R must have signs off the diagonal")
    return A, R


def commutator_parity_audit(A: np.ndarray, R: np.ndarray) -> dict[str, object]:
    """Check the exact diagonal identity and its even-order consequence."""
    A, R = _validate_pair(A, R)
    n = A.shape[0]
    commutator = A @ R - R @ A
    signed_row_sums = np.sum(A * R, axis=1, dtype=np.int64)
    predicted_diagonal = -2 * signed_row_sums
    identity_exact = np.array_equal(np.diag(commutator), predicted_diagonal)
    _require(identity_exact, "commutator diagonal identity failed")

    even_order = n % 2 == 0
    diagonal_nonzero = bool(np.all(np.diag(commutator) != 0)) if even_order else None
    diagonal_mod_four = [int(x % 4) for x in np.diag(commutator)]
    if even_order:
        _require(diagonal_nonzero, "even-order commutator acquired a zero diagonal entry")
        _require(
            all(x == 2 for x in diagonal_mod_four),
            "even-order commutator diagonal is not 2 modulo 4",
        )

    diagonal_frobenius_squared = int(np.dot(np.diag(commutator), np.diag(commutator)))
    frobenius_squared = int(np.sum(commutator * commutator, dtype=np.int64))
    return {
        "order": n,
        "even_order": even_order,
        "signed_row_sums": [int(x) for x in signed_row_sums],
        "commutator_diagonal": [int(x) for x in np.diag(commutator)],
        "diagonal_identity_exact": identity_exact,
        "diagonal_mod_four": diagonal_mod_four,
        "diagonal_nonzero": diagonal_nonzero,
        "diagonal_frobenius_squared": diagonal_frobenius_squared,
        "frobenius_squared": frobenius_squared,
        "even_order_lower_bound": 4 * n if even_order else None,
        "commutes": bool(np.all(commutator == 0)),
    }


def mod_two_skew_inverse_audit(n: int) -> dict[str, object]:
    """Verify ``(J+I)^2=I`` over F_2 for even ``n``.

    Every skew signing reduces modulo two to ``J-I=J+I``.  For even order,
    ``J^2=nJ=0`` over F_2, so this common reduction is its own inverse.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 2 or n % 2:
        raise ValueError("n must be an even integer at least two")
    reduction = (np.ones((n, n), dtype=np.int8) + np.eye(n, dtype=np.int8)) % 2
    square = (reduction.astype(np.int64) @ reduction.astype(np.int64)) % 2
    identity = np.eye(n, dtype=np.int64)
    involution_exact = np.array_equal(square, identity)
    _require(involution_exact, "the common mod-two skew reduction is not an involution")
    return {
        "order": n,
        "reduction": "J+I over F_2",
        "square_is_identity": involution_exact,
        "all_even_order_skew_signings_have_odd_determinant": True,
        "all_even_order_skew_signings_are_invertible_over_R": True,
    }


def symmetric_conference_commuting_no_go(n: int) -> dict[str, object]:
    """Return the two all-orders parity ledgers for conference order ``n``."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 2 or n % 4 != 2:
        raise ValueError("a real symmetric conference order must be 2 modulo 4")
    mod_two = mod_two_skew_inverse_audit(n)
    half_dimension = n // 2
    return {
        "order": n,
        "symmetric_conference_order_class": "2 mod 4",
        "commutator_diagonal_formula": "(AR-RA)_ii=-2*sum_(j!=i) A_ij R_ij",
        "number_of_diagonal_summands": n - 1,
        "number_of_diagonal_summands_is_odd": (n - 1) % 2 == 1,
        "commutator_diagonal_is_2_mod_4": True,
        "commutator_diagonal_frobenius_lower_bound": 4 * n,
        "conference_eigenspace_dimensions": [half_dimension, half_dimension],
        "conference_eigenspace_dimensions_are_odd": half_dimension % 2 == 1,
        "skew_signing_is_invertible_mod_two": mod_two["square_is_identity"],
        "commuting_skew_signing_exists": False,
        "proof_scope": "every skew zero-diagonal off-diagonal-sign matrix R",
    }


def paley_symmetric_conference(q: int) -> np.ndarray:
    """Construct the normalized real Paley conference matrix of order ``q+1``.

    This small verifier handles prime ``q=1 mod 4``; the theorem itself does
    not depend on the Paley construction.
    """
    if not isinstance(q, int) or isinstance(q, bool) or not _is_prime(q) or q % 4 != 1:
        raise ValueError("q must be a prime congruent to one modulo four")
    n = q + 1
    A = np.zeros((n, n), dtype=np.int64)
    A[0, 1:] = 1
    A[1:, 0] = 1
    for x in range(q):
        for y in range(q):
            if x == y:
                continue
            delta = (x - y) % q
            A[x + 1, y + 1] = 1 if pow(delta, (q - 1) // 2, q) == 1 else -1
    _require(np.array_equal(A, A.T), "Paley matrix is not symmetric")
    _require(np.array_equal(A @ A, q * np.eye(n, dtype=np.int64)), "Paley square failed")
    return A


def skew_signing_from_mask(n: int, mask: int) -> np.ndarray:
    """Decode an upper-triangular sign assignment from a nonnegative mask."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer at least two")
    edges = tuple(combinations(range(n), 2))
    if not isinstance(mask, int) or isinstance(mask, bool) or not 0 <= mask < 1 << len(edges):
        raise ValueError("mask is outside the upper-triangular signing range")
    R = np.zeros((n, n), dtype=np.int64)
    for bit, (i, j) in enumerate(edges):
        value = 1 if (mask >> bit) & 1 else -1
        R[i, j] = value
        R[j, i] = -value
    return R


def nomenclature_audit(A: np.ndarray, R: np.ndarray) -> dict[str, object]:
    """Translate OD amicability conventions for symmetric ``A`` and skew ``R``."""
    A, R = _validate_pair(A, R)
    amicable_defect = A @ R.T - R @ A.T
    anti_amicable_defect = A @ R.T + R @ A.T
    anticommutator = A @ R + R @ A
    commutator = A @ R - R @ A
    return {
        "amicable_defect_equals_negative_anticommutator": bool(
            np.array_equal(amicable_defect, -anticommutator)
        ),
        "anti_amicable_defect_equals_negative_commutator": bool(
            np.array_equal(anti_amicable_defect, -commutator)
        ),
        "amicable_means": "AR+RA=0",
        "anti_amicable_means": "AR-RA=0",
    }


def exhaustive_paley_order_six_audit() -> dict[str, object]:
    """Enumerate every skew signing against the Paley conference matrix at n=6."""
    A = paley_symmetric_conference(5)
    n = A.shape[0]
    number_of_signings = 1 << (n * (n - 1) // 2)
    commuting = 0
    minimum_diagonal_frobenius_squared: int | None = None
    minimum_frobenius_squared: int | None = None
    minimum_mask: int | None = None
    for mask in range(number_of_signings):
        R = skew_signing_from_mask(n, mask)
        commutator = A @ R - R @ A
        diagonal = np.diag(commutator)
        diagonal_square = int(np.dot(diagonal, diagonal))
        full_square = int(np.sum(commutator * commutator, dtype=np.int64))
        if np.all(commutator == 0):
            commuting += 1
        if minimum_frobenius_squared is None or full_square < minimum_frobenius_squared:
            minimum_frobenius_squared = full_square
            minimum_mask = mask
        if (
            minimum_diagonal_frobenius_squared is None
            or diagonal_square < minimum_diagonal_frobenius_squared
        ):
            minimum_diagonal_frobenius_squared = diagonal_square

    _require(commuting == 0, "an order-six commuting skew mate appeared")
    _require(
        minimum_diagonal_frobenius_squared is not None
        and minimum_diagonal_frobenius_squared >= 4 * n,
        "the all-orders diagonal lower bound failed at order six",
    )
    return {
        "order": n,
        "skew_signings_checked": number_of_signings,
        "commuting_signings": commuting,
        "minimum_diagonal_frobenius_squared": minimum_diagonal_frobenius_squared,
        "theoretical_diagonal_lower_bound": 4 * n,
        "minimum_full_commutator_frobenius_squared": minimum_frobenius_squared,
        "one_minimizing_mask": minimum_mask,
        "all_checks": True,
        "finite_audit_is_not_the_proof": True,
    }


def full_audit() -> dict[str, object]:
    A = paley_symmetric_conference(5)
    R = skew_signing_from_mask(6, 0x1234)
    return {
        "all_orders": [symmetric_conference_commuting_no_go(n) for n in (2, 6, 10, 14)],
        "representative_parity": commutator_parity_audit(A, R),
        "nomenclature": nomenclature_audit(A, R),
        "order_six_exhaustion": exhaustive_paley_order_six_audit(),
        "residual_ii_closed": False,
        "conference_commuting_route_closed": True,
    }


if __name__ == "__main__":
    print(json.dumps(full_audit(), indent=2, sort_keys=True))
