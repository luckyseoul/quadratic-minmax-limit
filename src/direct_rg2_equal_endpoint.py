#!/usr/bin/env python3
r"""Exact finite checks for the equal-endpoint skew RG2 construction.

Let ``A`` be a symmetric zero-diagonal signing and ``R`` a skew signing.
For ``i<j`` use the Hadamard cloud block

    [[A_ij, -R_ij],
     [R_ij,  A_ij]].

The two constant cloud states induce ``A`` at both endpoints.  If
``T={i:t_i=1}``, the mixed state induces

    Q_{C_t}(s) = I_A(T,s) - C_R(T,s),

where ``I_A`` is the energy internal to ``T`` and its complement and
``C_R=sum_{i in T,j not in T} R_ij s_i s_j``.  Complementing ``T`` changes
only the sign of ``C_R``.  Consequently the exact frame minimax is

    K(A,R) = max_{T,s} (|I_A(T,s)| + |C_R(T,s)|)
           = 1/2 max_{x,y} (|Q_A(x)+Q_A(y)| + |x^T R y|).

This module is deliberately finite and exact.  It verifies the algebra and
records two obstructions; it does not prove a multiplier-two estimate.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from typing import Iterator, Sequence


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _matrix(rows: Sequence[Sequence[int]]) -> Matrix:
    matrix = tuple(tuple(int(value) for value in row) for row in rows)
    if not matrix or any(len(row) != len(matrix) for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    return matrix


def _sign_vector(values: Sequence[int], n: int) -> Vector:
    vector = tuple(int(value) for value in values)
    if len(vector) != n or any(value not in (-1, 1) for value in vector):
        raise ValueError(f"expected a length-{n} sign vector")
    return vector


def _bits(values: Sequence[int], n: int) -> tuple[int, ...]:
    bits = tuple(int(value) for value in values)
    if len(bits) != n or any(value not in (0, 1) for value in bits):
        raise ValueError(f"expected a length-{n} bit vector")
    return bits


def require_symmetric_signing(rows: Sequence[Sequence[int]]) -> Matrix:
    """Return ``rows`` after exact validation as a Seidel signing."""
    matrix = _matrix(rows)
    n = len(matrix)
    if any(matrix[i][i] != 0 for i in range(n)):
        raise ValueError("symmetric signing must have zero diagonal")
    if any(
        matrix[i][j] not in (-1, 1) or matrix[i][j] != matrix[j][i]
        for i in range(n)
        for j in range(i + 1, n)
    ):
        raise ValueError("off-diagonal entries must be symmetric signs")
    return matrix


def require_skew_signing(rows: Sequence[Sequence[int]]) -> Matrix:
    """Return ``rows`` after exact validation as a skew signing."""
    matrix = _matrix(rows)
    n = len(matrix)
    if any(matrix[i][i] != 0 for i in range(n)):
        raise ValueError("skew signing must have zero diagonal")
    if any(
        matrix[i][j] not in (-1, 1) or matrix[j][i] != -matrix[i][j]
        for i in range(n)
        for j in range(i + 1, n)
    ):
        raise ValueError("off-diagonal entries must be skew signs")
    return matrix


@lru_cache(maxsize=None)
def sign_vectors(n: int) -> tuple[Vector, ...]:
    if n < 0:
        raise ValueError("dimension must be nonnegative")
    return tuple(product((-1, 1), repeat=n))


@lru_cache(maxsize=None)
def bit_vectors(n: int) -> tuple[tuple[int, ...], ...]:
    if n < 0:
        raise ValueError("dimension must be nonnegative")
    return tuple(product((0, 1), repeat=n))


def symmetric_signings(n: int) -> Iterator[Matrix]:
    """Generate all order-``n`` symmetric signings."""
    if n < 1:
        raise ValueError("order must be positive")
    edges = tuple((i, j) for i in range(n) for j in range(i + 1, n))
    for signs in product((-1, 1), repeat=len(edges)):
        matrix = [[0] * n for _ in range(n)]
        for (i, j), value in zip(edges, signs):
            matrix[i][j] = matrix[j][i] = value
        yield tuple(tuple(row) for row in matrix)


def skew_signings(n: int) -> Iterator[Matrix]:
    """Generate all order-``n`` skew signings."""
    if n < 1:
        raise ValueError("order must be positive")
    edges = tuple((i, j) for i in range(n) for j in range(i + 1, n))
    for signs in product((-1, 1), repeat=len(edges)):
        matrix = [[0] * n for _ in range(n)]
        for (i, j), value in zip(edges, signs):
            matrix[i][j] = value
            matrix[j][i] = -value
        yield tuple(tuple(row) for row in matrix)


def quadratic_energy(rows: Sequence[Sequence[int]], x: Sequence[int]) -> int:
    """Return ``Q_A(x)=sum_{i<j} A_ij x_i x_j``."""
    matrix = _matrix(rows)
    vector = _sign_vector(x, len(matrix))
    return sum(
        matrix[i][j] * vector[i] * vector[j]
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
    )


def phi(rows: Sequence[Sequence[int]]) -> int:
    matrix = require_symmetric_signing(rows)
    return max(abs(quadratic_energy(matrix, x)) for x in sign_vectors(len(matrix)))


def skew_bilinear(
    rows: Sequence[Sequence[int]], x: Sequence[int], y: Sequence[int]
) -> int:
    matrix = require_skew_signing(rows)
    left = _sign_vector(x, len(matrix))
    right = _sign_vector(y, len(matrix))
    return sum(
        left[i] * matrix[i][j] * right[j]
        for i in range(len(matrix))
        for j in range(len(matrix))
    )


def equal_endpoint_block(a: int, r: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """The oriented ``2 by 2`` Hadamard block for one edge."""
    if a not in (-1, 1) or r not in (-1, 1):
        raise ValueError("a and r must be signs")
    return ((a, -r), (r, a))


def induced_edge(a: int, r: int, t_i: int, t_j: int) -> int:
    """Return the signing induced by the two cloud-state bits."""
    if t_i not in (0, 1) or t_j not in (0, 1):
        raise ValueError("cloud states must be bits")
    if t_i == t_j:
        return a
    return r if (t_i, t_j) == (0, 1) else -r


def block_state_value(
    a: int,
    r: int,
    s_i: int,
    t_i: int,
    s_j: int,
    t_j: int,
) -> int:
    """Compute ``z_i^T B z_j / 2`` directly from the four block entries."""
    if s_i not in (-1, 1) or s_j not in (-1, 1):
        raise ValueError("cloud signs must be signs")
    block = equal_endpoint_block(a, r)
    left = (s_i, s_i * (-1) ** t_i)
    right = (s_j, s_j * (-1) ** t_j)
    numerator = sum(
        left[u] * block[u][v] * right[v] for u in range(2) for v in range(2)
    )
    if numerator % 2:
        raise ArithmeticError("Hadamard block value was not integral")
    return numerator // 2


def exhaustive_local_block_certificate() -> dict[str, object]:
    """Check every one-edge choice and every one-cloud state exactly."""
    checked = 0
    for a, r, s_i, s_j in product((-1, 1), repeat=4):
        for t_i, t_j in product((0, 1), repeat=2):
            actual = block_state_value(a, r, s_i, t_i, s_j, t_j)
            expected = s_i * s_j * induced_edge(a, r, t_i, t_j)
            if actual != expected:
                raise ArithmeticError("equal-endpoint block identity failed")
            checked += 1
    return {
        "cases_checked": checked,
        "both_endpoints_equal_A": True,
        "mixed_states_are_skew_cuts": True,
        "proved": True,
    }


def induced_signing(
    a_rows: Sequence[Sequence[int]],
    r_rows: Sequence[Sequence[int]],
    t: Sequence[int],
) -> Matrix:
    """Return the order-``n`` signing induced by a cloud-state vector."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    bits = _bits(t, len(a))
    result = [[0] * len(a) for _ in a]
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            value = induced_edge(a[i][j], r[i][j], bits[i], bits[j])
            result[i][j] = result[j][i] = value
    return tuple(tuple(row) for row in result)


def equal_endpoint_lift(
    a_rows: Sequence[Sequence[int]],
    r_rows: Sequence[Sequence[int]],
    internal_signs: Sequence[int],
) -> Matrix:
    """Build the full order-``2n`` signing, with vertices grouped by cloud."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    d = _sign_vector(internal_signs, len(a))
    result = [[0] * (2 * len(a)) for _ in range(2 * len(a))]
    for i in range(len(a)):
        result[2 * i][2 * i + 1] = result[2 * i + 1][2 * i] = d[i]
        for j in range(i + 1, len(a)):
            block = equal_endpoint_block(a[i][j], r[i][j])
            for u in range(2):
                for v in range(2):
                    value = block[u][v]
                    result[2 * i + u][2 * j + v] = value
                    result[2 * j + v][2 * i + u] = value
    return tuple(tuple(row) for row in result)


def cloud_vector(s: Sequence[int], t: Sequence[int]) -> Vector:
    signs = _sign_vector(s, len(s))
    bits = _bits(t, len(signs))
    return tuple(
        coordinate
        for sign, bit in zip(signs, bits)
        for coordinate in (sign, sign * (-1) ** bit)
    )


def cut_internal_energy(
    a_rows: Sequence[Sequence[int]], t: Sequence[int], s: Sequence[int]
) -> int:
    """Energy of ``A`` on pairs lying on the same side of the cut."""
    a = require_symmetric_signing(a_rows)
    bits = _bits(t, len(a))
    signs = _sign_vector(s, len(a))
    return sum(
        a[i][j] * signs[i] * signs[j]
        for i in range(len(a))
        for j in range(i + 1, len(a))
        if bits[i] == bits[j]
    )


def symmetric_cut_energy(
    a_rows: Sequence[Sequence[int]], t: Sequence[int], s: Sequence[int]
) -> int:
    """The once-counted ``A`` energy crossing from ``T`` to its complement."""
    a = require_symmetric_signing(a_rows)
    bits = _bits(t, len(a))
    signs = _sign_vector(s, len(a))
    return sum(
        a[i][j] * signs[i] * signs[j]
        for i in range(len(a))
        if bits[i] == 1
        for j in range(len(a))
        if bits[j] == 0
    )


def skew_cut_energy(
    r_rows: Sequence[Sequence[int]], t: Sequence[int], s: Sequence[int]
) -> int:
    """Return ``sum_{i in T,j outside T} R_ij s_i s_j``."""
    r = require_skew_signing(r_rows)
    bits = _bits(t, len(r))
    signs = _sign_vector(s, len(r))
    return sum(
        r[i][j] * signs[i] * signs[j]
        for i in range(len(r))
        if bits[i] == 1
        for j in range(len(r))
        if bits[j] == 0
    )


def equal_endpoint_k_by_frames(
    a_rows: Sequence[Sequence[int]], r_rows: Sequence[Sequence[int]]
) -> int:
    """Compute ``max_t Phi(C_t)`` from the induced signings."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    return max(phi(induced_signing(a, r, t)) for t in bit_vectors(len(a)))


def equal_endpoint_k_by_cuts(
    a_rows: Sequence[Sequence[int]], r_rows: Sequence[Sequence[int]]
) -> int:
    """Compute ``max_(T,s) |I_A|+|C_R|``."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    return max(
        abs(cut_internal_energy(a, t, s)) + abs(skew_cut_energy(r, t, s))
        for t in bit_vectors(len(a))
        for s in sign_vectors(len(a))
    )


def equal_endpoint_k_by_pairs(
    a_rows: Sequence[Sequence[int]], r_rows: Sequence[Sequence[int]]
) -> int:
    r"""Compute ``max (|Q_A(x)+Q_A(y)|+|x^T R y|)/2``."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    numerator = max(
        abs(quadratic_energy(a, x) + quadratic_energy(a, y))
        + abs(skew_bilinear(r, x, y))
        for x in sign_vectors(len(a))
        for y in sign_vectors(len(a))
    )
    if numerator % 2:
        raise ArithmeticError("pair formula did not have even numerator")
    return numerator // 2


def verify_equal_endpoint_algebra(
    a_rows: Sequence[Sequence[int]], r_rows: Sequence[Sequence[int]]
) -> dict[str, object]:
    """Exhaustively replay all equivalent RG2 formulas for one ``(A,R)``."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    n = len(a)
    zero = (0,) * n
    one = (1,) * n
    if induced_signing(a, r, zero) != a or induced_signing(a, r, one) != a:
        raise ArithmeticError("the two endpoint signings differ from A")

    cut_cases = 0
    for t in bit_vectors(n):
        complement = tuple(1 - bit for bit in t)
        for s in sign_vectors(n):
            internal = cut_internal_energy(a, t, s)
            cut = skew_cut_energy(r, t, s)
            paired = tuple(sign * (-1) ** bit for sign, bit in zip(s, t))
            induced = quadratic_energy(induced_signing(a, r, t), s)
            induced_complement = quadratic_energy(
                induced_signing(a, r, complement), s
            )
            if induced != internal - cut:
                raise ArithmeticError("skew cut identity failed")
            if induced_complement != internal + cut:
                raise ArithmeticError("complementary skew cut identity failed")
            if quadratic_energy(a, s) + quadratic_energy(a, paired) != 2 * internal:
                raise ArithmeticError("pair-to-internal identity failed")
            if skew_bilinear(r, s, paired) != 2 * cut:
                raise ArithmeticError("pair-to-skew-cut identity failed")
            cut_cases += 1

    by_frames = equal_endpoint_k_by_frames(a, r)
    by_cuts = equal_endpoint_k_by_cuts(a, r)
    by_pairs = equal_endpoint_k_by_pairs(a, r)
    if len({by_frames, by_cuts, by_pairs}) != 1:
        raise ArithmeticError("the three equal-endpoint minimax formulas disagree")

    # Include arbitrary internal signs in an exact full-lift replay.
    d = tuple(1 if i % 2 == 0 else -1 for i in range(n))
    lift = equal_endpoint_lift(a, r, d)
    lift_cases = 0
    for t in bit_vectors(n):
        induced = induced_signing(a, r, t)
        matching = sum(d[i] * (-1) ** t[i] for i in range(n))
        for s in sign_vectors(n):
            actual = quadratic_energy(lift, cloud_vector(s, t))
            expected = 2 * quadratic_energy(induced, s) + matching
            if actual != expected:
                raise ArithmeticError("full equal-endpoint lift identity failed")
            lift_cases += 1

    return {
        "n": n,
        "endpoint_zero_equals_A": True,
        "endpoint_one_equals_A": True,
        "cut_cases_checked": cut_cases,
        "lift_cases_checked": lift_cases,
        "K": by_frames,
        "K_by_frames": by_frames,
        "K_by_cuts": by_cuts,
        "K_by_pairs": by_pairs,
        "proved_for_input": True,
    }


def _principal_extrema(a: Matrix, indices: tuple[int, ...]) -> tuple[int, int]:
    values = [
        sum(
            a[indices[u]][indices[v]] * x[u] * x[v]
            for u in range(len(indices))
            for v in range(u + 1, len(indices))
        )
        for x in sign_vectors(len(indices))
    ]
    return max(values), -min(values)


def hereditary_endpoint_certificate(
    a_rows: Sequence[Sequence[int]],
) -> dict[str, object]:
    """Replay the endpoint hereditary bounds, which hold for every ``A``.

    For fixed spins on both sides write the full energies as ``D+X`` and
    ``D-X`` after flipping one side.  Both have magnitude at most ``Phi(A)``,
    hence ``|D|+|X|<=Phi(A)``.  Taking independent within-side maxima or
    minima gives both inequalities from Proposition 6.4 automatically when
    the two endpoints coincide.
    """
    a = require_symmetric_signing(a_rows)
    n = len(a)
    bound = phi(a)
    max_positive_sum = 0
    max_negative_sum = 0
    max_cut_l1 = 0
    for t in bit_vectors(n):
        inside = tuple(i for i, bit in enumerate(t) if bit)
        outside = tuple(i for i, bit in enumerate(t) if not bit)
        positive_inside, negative_inside = _principal_extrema(a, inside)
        positive_outside, negative_outside = _principal_extrema(a, outside)
        max_positive_sum = max(
            max_positive_sum, positive_inside + positive_outside
        )
        max_negative_sum = max(
            max_negative_sum, negative_inside + negative_outside
        )
        for s in sign_vectors(n):
            internal = cut_internal_energy(a, t, s)
            cross = symmetric_cut_energy(a, t, s)
            direct = quadratic_energy(a, s)
            flipped = tuple(-s[i] if t[i] else s[i] for i in range(n))
            if direct != internal + cross:
                raise ArithmeticError("symmetric cut decomposition failed")
            if quadratic_energy(a, flipped) != internal - cross:
                raise ArithmeticError("symmetric cut-flip decomposition failed")
            if max(abs(direct), abs(quadratic_energy(a, flipped))) != (
                abs(internal) + abs(cross)
            ):
                raise ArithmeticError("cut l1 identity failed")
            max_cut_l1 = max(max_cut_l1, abs(internal) + abs(cross))

    if max(max_positive_sum, max_negative_sum, max_cut_l1) > bound:
        raise ArithmeticError("automatic hereditary endpoint bound failed")
    return {
        "n": n,
        "Phi_A": bound,
        "max_P_A_T_plus_P_A_Tc": max_positive_sum,
        "max_N_A_T_plus_N_A_Tc": max_negative_sum,
        "max_cut_abs_internal_plus_abs_cross": max_cut_l1,
        "hereditary_endpoint_bounds_automatic": True,
        "proved_for_input": True,
    }


def simple_walk_absolute_mean(length: int) -> Fraction:
    """Exact ``E|epsilon_1+...+epsilon_length|`` for fair signs."""
    if length < 0:
        raise ValueError("walk length must be nonnegative")
    direct = Fraction(
        sum(comb(length, k) * abs(length - 2 * k) for k in range(length + 1)),
        2**length,
    )
    if length == 0:
        closed = Fraction(0)
    else:
        closed = Fraction(
            length * comb(length - 1, (length - 1) // 2), 2 ** (length - 1)
        )
    if direct != closed:
        raise ArithmeticError("simple-walk closed form failed")
    return closed


def skew_norm_floor_exact(n: int) -> Fraction:
    r"""Exact averaging floor ``n E|S_(n-1)|`` for ``||R||_(inf->1)``."""
    if n < 1:
        raise ValueError("order must be positive")
    return n * simple_walk_absolute_mean(n - 1)


def skew_infinity_to_one(r_rows: Sequence[Sequence[int]]) -> int:
    r"""Compute ``max_(x,y) |x^T R y| = max_y ||R y||_1`` exactly."""
    r = require_skew_signing(r_rows)
    return max(
        sum(
            abs(sum(r[i][j] * y[j] for j in range(len(r))))
            for i in range(len(r))
        )
        for y in sign_vectors(len(r))
    )


def skew_norm_floor_certificate(
    r_rows: Sequence[Sequence[int]],
) -> dict[str, object]:
    """Replay the exact random-sign average underlying the skew norm floor."""
    r = require_skew_signing(r_rows)
    n = len(r)
    row_l1_values = [
        sum(abs(sum(r[i][j] * y[j] for j in range(n))) for i in range(n))
        for y in sign_vectors(n)
    ]
    average = Fraction(sum(row_l1_values), len(row_l1_values))
    formula = skew_norm_floor_exact(n)
    norm = max(row_l1_values)
    if average != formula or norm < formula:
        raise ArithmeticError("skew infinity-to-one averaging floor failed")
    return {
        "n": n,
        "infinity_to_one_norm": norm,
        "exact_average_over_y": average,
        "exact_floor_formula": formula,
        "floor_verified": True,
        "proved_for_input": True,
    }


def rg2_disk_value_squared(
    a_rows: Sequence[Sequence[int]],
    r_rows: Sequence[Sequence[int]],
    x: Sequence[int],
    y: Sequence[int],
) -> int:
    r"""Return ``((Q_A(x)+Q_A(y))/2)^2 + ((x^T R y)/2)^2``."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    internal_numerator = quadratic_energy(a, x) + quadratic_energy(a, y)
    cross_numerator = skew_bilinear(r, x, y)
    if internal_numerator % 2 or cross_numerator % 2:
        raise ArithmeticError("RG2 disk coordinates were not integral")
    return (internal_numerator // 2) ** 2 + (cross_numerator // 2) ** 2


def zero_error_disk_holds(
    a_rows: Sequence[Sequence[int]], r_rows: Sequence[Sequence[int]]
) -> bool:
    """Whether the pointwise Euclidean disk has radius exactly ``Phi(A)``."""
    a = require_symmetric_signing(a_rows)
    r = require_skew_signing(r_rows)
    if len(a) != len(r):
        raise ValueError("A and R must have the same order")
    radius_squared = phi(a) ** 2
    return all(
        rg2_disk_value_squared(a, r, x, y) <= radius_squared
        for x in sign_vectors(len(a))
        for y in sign_vectors(len(a))
    )


def n5_cycle_chord_signing() -> Matrix:
    """The ``C5`` signing: negative cycle edges and positive chords."""
    n = 5
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = -1 if (j - i) % n in (1, n - 1) else 1
            result[i][j] = result[j][i] = value
    return tuple(tuple(row) for row in result)


def n5_maximizer_basis() -> Matrix:
    """Rows ``v_r`` have negative coordinates ``r,r+2 (mod 5)``."""
    rows: list[Vector] = []
    for r in range(5):
        row = [1] * 5
        row[r] = -1
        row[(r + 2) % 5] = -1
        rows.append(tuple(row))
    return tuple(rows)


def exact_determinant(rows: Sequence[Sequence[int]]) -> int:
    """Fraction-free Bareiss determinant over the integers."""
    matrix = [list(row) for row in _matrix(rows)]
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if matrix[k][k] == 0:
            pivot_row = next((i for i in range(k + 1, n) if matrix[i][k]), None)
            if pivot_row is None:
                return 0
            matrix[k], matrix[pivot_row] = matrix[pivot_row], matrix[k]
            sign = -sign
        pivot = matrix[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = matrix[i][j] * pivot - matrix[i][k] * matrix[k][j]
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                matrix[i][j] = numerator // previous
            matrix[i][k] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def n5_disk_counterexample() -> dict[str, object]:
    """Certify that the zero-error Euclidean disk ansatz already fails at 5.

    The five displayed positive maximizers form an invertible matrix ``V``.
    A radius-``Phi(A)`` disk would force ``v_r^T R v_s=0`` for every pair,
    hence ``V R V^T=0`` and therefore ``R=0``.  A skew signing cannot be
    zero.  The finite replay additionally checks all ``2^10`` skew signings.
    This says nothing against an ``o(n^3)`` disk error asymptotically.
    """
    a = n5_cycle_chord_signing()
    vectors = sign_vectors(5)
    energies = tuple(quadratic_energy(a, x) for x in vectors)
    moment_two = Fraction(sum(value * value for value in energies), len(energies))
    bound = phi(a)
    basis = n5_maximizer_basis()
    basis_energies = tuple(quadratic_energy(a, row) for row in basis)
    determinant = exact_determinant(basis)

    if moment_two != 10 or any(value % 2 for value in energies):
        raise ArithmeticError("n=5 moment/parity replay failed")
    if bound != 4 or basis_energies != (4,) * 5 or determinant == 0:
        raise ArithmeticError("n=5 maximizer basis replay failed")

    signing_count = 0
    zero_error_count = 0
    minimum_anchor_cross = None
    for r in skew_signings(5):
        signing_count += 1
        anchor_cross = max(
            abs(skew_bilinear(r, x, y)) for x in basis for y in basis
        )
        minimum_anchor_cross = (
            anchor_cross
            if minimum_anchor_cross is None
            else min(minimum_anchor_cross, anchor_cross)
        )
        if anchor_cross == 0:
            zero_error_count += 1

    if signing_count != 2**10 or zero_error_count != 0:
        raise ArithmeticError("n=5 skew-signing exclusion replay failed")
    if minimum_anchor_cross is None:
        raise ArithmeticError("no n=5 skew signings were generated")

    minimum_anchor_disk_value_squared = bound**2 + (minimum_anchor_cross // 2) ** 2
    return {
        "n": 5,
        "Phi_A": bound,
        "energy_alphabet": sorted(set(energies)),
        "all_energies_even": True,
        "E_Q_squared": moment_two,
        "moment_and_parity_force_Phi_at_least_4": True,
        "maximizer_basis_energies": list(basis_energies),
        "determinant_V": determinant,
        "V_is_invertible": True,
        "zero_error_disk_forces_V_R_V_transpose_zero": True,
        "V_R_V_transpose_zero_forces_R_zero": True,
        "skew_signings_checked": signing_count,
        "skew_signings_passing_zero_error_anchor_constraints": zero_error_count,
        "minimum_max_anchor_abs_x_R_y": minimum_anchor_cross,
        "disk_radius_squared": bound**2,
        "minimum_anchor_disk_value_squared": minimum_anchor_disk_value_squared,
        "zero_error_disk_impossible": True,
        "refutes_only_zero_error_disk": True,
        "does_not_refute_asymptotic_o_n_cubed_error": True,
        "result_status": "counterexample to the zero-error disk ansatz only",
    }


def exhaustive_small_order_certificate(max_order: int = 3) -> dict[str, object]:
    """Replay the global identities for every ``(A,R)`` through an order."""
    if max_order < 1 or max_order > 4:
        raise ValueError("small-order replay supports 1<=max_order<=4")
    local = exhaustive_local_block_certificate()
    pairs_checked = 0
    by_order: dict[int, int] = {}
    for n in range(1, max_order + 1):
        count = 0
        for a in symmetric_signings(n):
            hereditary_endpoint_certificate(a)
            for r in skew_signings(n):
                verify_equal_endpoint_algebra(a, r)
                skew_norm_floor_certificate(r)
                pairs_checked += 1
                count += 1
        by_order[n] = count
    return {
        "max_order": max_order,
        "local_block_cases": local["cases_checked"],
        "A_R_pairs_by_order": by_order,
        "A_R_pairs_checked": pairs_checked,
        "exact_integer_arithmetic": True,
        "proved_for_enumerated_orders": True,
    }


__all__ = [
    "bit_vectors",
    "block_state_value",
    "cloud_vector",
    "cut_internal_energy",
    "equal_endpoint_block",
    "equal_endpoint_k_by_cuts",
    "equal_endpoint_k_by_frames",
    "equal_endpoint_k_by_pairs",
    "equal_endpoint_lift",
    "exact_determinant",
    "exhaustive_local_block_certificate",
    "exhaustive_small_order_certificate",
    "hereditary_endpoint_certificate",
    "induced_edge",
    "induced_signing",
    "n5_cycle_chord_signing",
    "n5_disk_counterexample",
    "n5_maximizer_basis",
    "phi",
    "quadratic_energy",
    "require_skew_signing",
    "require_symmetric_signing",
    "rg2_disk_value_squared",
    "sign_vectors",
    "simple_walk_absolute_mean",
    "skew_bilinear",
    "skew_cut_energy",
    "skew_infinity_to_one",
    "skew_norm_floor_certificate",
    "skew_norm_floor_exact",
    "skew_signings",
    "symmetric_cut_energy",
    "symmetric_signings",
    "verify_equal_endpoint_algebra",
    "zero_error_disk_holds",
]
