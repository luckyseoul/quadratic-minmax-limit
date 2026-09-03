#!/usr/bin/env python3
r"""Exact audit for the Kiss--Somlai triangular all-bad three-spike family.

This is an unnumbered construction/barrier, not a residual-(ii) closure.
For every odd prime ``p>=5`` it constructs an integral ``+p`` eigenvector
of the Paley conference matrix on ``P^1(F_(p^2))`` whose coordinates have
absolute value one except for exactly three coordinates equal to ``+3``.
Replacing those three entries by ``+1`` gives a Boolean vector of conference
defect ``6p-12`` in the all-bad signed-triple case.

The geometric input is the triangular set of Kiss and Somlai,
``S={(a,b): 0<=b<a<=p-1}``, which has exactly three special directions.

References:
  https://arxiv.org/abs/2109.13992
  https://doi.org/10.1007/s10623-024-01404-y

The exact checks below deliberately keep three signs separate:

* ``chi(d)`` is the character of a *spatial* affine-line direction;
* the Fourier support of that line lies on its trace annihilator;
* the quadratic Gauss sign converts the annihilator character into the
  eigenvalue sign of the finite Paley convolution operator.

Thus a square spatial direction has ``+p`` eigensign even when its Fourier
annihilator direction is nonsquare.  The direct integer line identity is
also checked, so the construction does not depend on a floating-point
Fourier transform.
"""
from __future__ import annotations

import math
from itertools import combinations

import numpy as np

from e1_gmin_m4_prop15598 import field_ctx, legendre
from e1_gmin_m4_prop15721 import is_prime
from minmax_quadratic import paley_conference_prime_power


BASE_SPECIAL_DIRECTIONS = ((1, 0), (0, 1), (1, 1))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _validate_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 5
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("the triangular three-spike construction needs an odd prime p>=5")


def _encode(p: int, a: int, b: int) -> int:
    return (a % p) + (b % p) * p


def _coordinates(p: int, u: int) -> tuple[int, int]:
    return u % p, u // p


def _add(p: int, u: int, v: int) -> int:
    a, b = _coordinates(p, u)
    c, d = _coordinates(p, v)
    return _encode(p, a + c, b + d)


def _neg(p: int, u: int) -> int:
    a, b = _coordinates(p, u)
    return _encode(p, -a, -b)


def _sub(p: int, u: int, v: int) -> int:
    return _add(p, u, _neg(p, v))


def _scalar(p: int, t: int, u: int) -> int:
    a, b = _coordinates(p, u)
    return _encode(p, t * a, t * b)


def _field_power(mul, x: int, exponent: int) -> int:
    out = 1
    base = x
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        base = mul(base, base)
        exponent >>= 1
    return out


def projective_directions(p: int) -> tuple[tuple[int, int], ...]:
    """Canonical direction representatives ``(1,m)`` and ``(0,1)``."""
    _validate_prime(p)
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def triangular_indicator(p: int) -> np.ndarray:
    """Return ``1_S`` for ``S={(a,b): 0<=b<a<=p-1}``."""
    _validate_prime(p)
    out = np.zeros(p * p, dtype=np.int64)
    for b in range(p):
        for a in range(p):
            out[_encode(p, a, b)] = int(b < a)
    return out


def triangular_augmented_function(p: int) -> np.ndarray:
    r"""Return ``f0=1_S+1_{a=0}+1_{b=p-2}`` on ``F_p^2``."""
    out = triangular_indicator(p)
    for b in range(p):
        out[_encode(p, 0, b)] += 1
    for a in range(p):
        out[_encode(p, a, p - 2)] += 1
    return out


def parallel_line_sums(
    p: int, values: np.ndarray, direction: tuple[int, int]
) -> tuple[int, ...]:
    """Line sums in a geometric direction, indexed by a normal coordinate."""
    if values.shape != (p * p,):
        raise ValueError("values must be a function on F_p^2")
    dx, dy = direction
    if dx % p == 0 and dy % p == 0:
        raise ValueError("zero is not a projective direction")
    sums = [0] * p
    for b in range(p):
        for a in range(p):
            intercept = (-dy * a + dx * b) % p
            sums[intercept] += int(values[_encode(p, a, b)])
    return tuple(sums)


def kiss_somlai_direction_audit(p: int) -> dict[str, object]:
    """Representative-prime audit of the three special geometric directions."""
    indicator = triangular_indicator(p)
    expected_flat = (p - 1) // 2
    rows: dict[str, list[int]] = {}
    special: list[tuple[int, int]] = []
    for direction in projective_directions(p):
        sums = parallel_line_sums(p, indicator, direction)
        rows[str(direction)] = list(sums)
        if len(set(sums)) != 1:
            special.append(direction)
        else:
            _require(
                sums == (expected_flat,) * p,
                "a nonspecial Kiss--Somlai direction has the wrong line sum",
            )
    _require(
        set(special) == set(BASE_SPECIAL_DIRECTIONS),
        "the Kiss--Somlai special-direction set changed",
    )
    return {
        "special_spatial_directions": [list(row) for row in special],
        "nonspecial_line_sum": expected_flat,
        "line_sums": rows,
        "proved_for_this_prime": True,
    }


def find_square_triangle_linear_map(p: int) -> tuple[int, int]:
    r"""Deterministically find the images of ``e_a,e_b`` under ``T``.

    The lexicographically first pair ``(A,B)`` is used subject to
    ``A,B`` being independent and ``A,B,A+B`` all having quadratic
    character ``+1``.  Hence ``T`` sends the horizontal, vertical, and
    slope-one geometric directions to three distinct square directions.

    Existence for every ``p>=5`` is not an experimental assertion: the
    ``p+1`` projective directions split into ``(p+1)/2`` square directions.
    Choose any three distinct square directions and rescale representatives
    in ``F_p^*`` so that the third is the sum of the first two.
    """
    _validate_prime(p)
    q, _mul, _field_add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    for image_a in range(1, q):
        if chi(image_a) != 1:
            continue
        aa, ab = _coordinates(p, image_a)
        for image_b in range(1, q):
            if chi(image_b) != 1:
                continue
            ba, bb = _coordinates(p, image_b)
            if (aa * bb - ab * ba) % p == 0:
                continue
            if chi(_add(p, image_a, image_b)) == 1:
                return image_a, image_b
    raise ArithmeticError("failed to find the guaranteed square direction triangle")


def transform_function(
    p: int, values: np.ndarray, image_a: int, image_b: int
) -> np.ndarray:
    """Return ``values o T^{-1}``, where columns of ``T`` are supplied."""
    if values.shape != (p * p,):
        raise ValueError("values must be a function on F_p^2")
    aa, ab = _coordinates(p, image_a)
    ba, bb = _coordinates(p, image_b)
    _require((aa * bb - ab * ba) % p != 0, "T must be invertible")
    out = np.empty_like(values)
    seen: set[int] = set()
    for b in range(p):
        for a in range(p):
            target = _add(
                p,
                _scalar(p, a, image_a),
                _scalar(p, b, image_b),
            )
            out[target] = values[_encode(p, a, b)]
            seen.add(target)
    _require(len(seen) == p * p, "T did not permute the affine plane")
    return out


def finite_line_identity(
    p: int, Q: np.ndarray, direction: int
) -> dict[str, object]:
    r"""Audit ``Q 1_L=chi(d)(p 1_L-1)`` for ``L=d F_p``."""
    q, mul, _field_add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    line = np.zeros(q, dtype=np.int64)
    for t in range(p):
        line[mul(direction, t)] = 1
    spatial_character = int(chi(direction))
    lhs = Q @ line
    rhs = spatial_character * (p * line - np.ones(q, dtype=np.int64))
    exact = bool(np.array_equal(lhs, rhs))
    _require(exact, "the spatial-direction Paley line identity failed")
    return {
        "direction": list(_coordinates(p, direction)),
        "spatial_direction_character": spatial_character,
        "line_contrast_eigenvalue": spatial_character * p,
        "Q_1L_equals_chi_d_times_p1L_minus_1": exact,
    }


def fourier_direction_ledger(p: int, directions: tuple[int, ...]) -> dict[str, object]:
    r"""Keep spatial type, Fourier-annihilator type, and eigensign separate."""
    q, mul, field_add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    eta_minus_one = legendre(-1, p)
    gauss_sign = -eta_minus_one

    trace_zero = next(
        u for u in range(1, q) if field_add(u, frob(u)) == 0
    )
    trace_zero_character = int(chi(trace_zero))
    _require(
        trace_zero_character == -eta_minus_one,
        "the trace-zero direction character has the wrong sign",
    )

    rows = []
    for direction in directions:
        inverse = _field_power(mul, direction, q - 2)
        annihilator = mul(trace_zero, inverse)
        spatial_character = int(chi(direction))
        annihilator_character = int(chi(annihilator))
        trace_pair = field_add(mul(annihilator, direction), frob(mul(annihilator, direction)))
        multiplier_sign = gauss_sign * annihilator_character
        _require(trace_pair == 0, "the claimed Fourier direction is not the annihilator")
        _require(
            annihilator_character == trace_zero_character * spatial_character,
            "the Fourier-annihilator character twist failed",
        )
        _require(
            multiplier_sign == spatial_character,
            "the Fourier multiplier and spatial line eigensign disagree",
        )
        rows.append(
            {
                "spatial_direction": list(_coordinates(p, direction)),
                "spatial_direction_character": spatial_character,
                "annihilator_direction": list(_coordinates(p, annihilator)),
                "annihilator_direction_character": annihilator_character,
                "quadratic_gauss_sign": gauss_sign,
                "finite_Paley_multiplier_sign": multiplier_sign,
            }
        )
    return {
        "fourier_convention": "fhat(xi)=sum_z f(z) exp(-2*pi*i*Tr(xi*z)/p)",
        "gauss_sum_over_Fp2": f"{gauss_sign}*p",
        "trace_zero_direction_character": trace_zero_character,
        "rows": rows,
        "proved_for_this_prime": True,
    }


def signed_pgl_triangle_transitivity(p: int) -> dict[str, object]:
    r'''Audit the single signed-PSL orbit of positive conference triples.

    For canonical homogeneous representatives r_P and g in GL(2,q), write
    g r_P=lambda_P r_(gP). Taking determinant characters gives

        C_(gP,gQ)=chi(det g) chi(lambda_P) chi(lambda_Q) C_(P,Q).

    Thus square-determinant projectivities preserve conference-triangle
    sign after switching by s_P=chi(lambda_P). PGL is sharply
    three-transitive; the unique projectivity between two positive ordered
    triples has square determinant.

    This proves one orbit for the signed-triple shell datum. It supplies a
    transported triangular completion for every datum, but does not assert
    uniqueness of all Boolean/eigenvector completions of a fixed datum.
    '''
    _validate_prime(p)
    q, mul, field_add, chi, _frob, _norm, _ia, _ib = field_ctx(p)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)

    def fdiv(u: int, v: int) -> int:
        if v == 0:
            raise ZeroDivisionError('field division by zero')
        return mul(u, _field_power(mul, v, q - 2))

    def matrix_action_and_scale(
        point: int, matrix: tuple[int, int, int, int]
    ) -> tuple[int, int]:
        a, b, c, d = matrix
        if point == 0:
            X, Y = 1, 0
        else:
            X, Y = point - 1, 1
        new_X = field_add(mul(a, X), mul(b, Y))
        new_Y = field_add(mul(c, X), mul(d, Y))
        if new_Y == 0:
            _require(new_X != 0, 'an invertible projectivity mapped a point to zero')
            return 0, new_X
        return 1 + fdiv(new_X, new_Y), new_Y

    # A nontrivial square-determinant lift, valid in every characteristic p>=5.
    matrix = (1, 1, 1, 2)
    determinant = _sub(p, mul(matrix[0], matrix[3]), mul(matrix[1], matrix[2]))
    _require(determinant != 0 and chi(determinant) == 1, 'bad PSL audit matrix')
    images: list[int] = []
    switches: list[int] = []
    for point in range(q + 1):
        image, scale = matrix_action_and_scale(point, matrix)
        images.append(image)
        switches.append(int(chi(scale)))
    _require(len(set(images)) == q + 1, 'the projectivity is not a permutation')

    factor_exact = True
    determinant_character = int(chi(determinant))
    for left in range(q + 1):
        for right in range(left + 1, q + 1):
            expected = (
                determinant_character
                * switches[left]
                * switches[right]
                * int(C[left, right])
            )
            if int(C[images[left], images[right]]) != expected:
                factor_exact = False
                break
        if not factor_exact:
            break
    _require(factor_exact, 'the signed-PGL switching factor failed')

    conference_square_exact = bool(
        np.array_equal(C @ C, q * np.eye(q + 1, dtype=np.int64))
    )
    _require(conference_square_exact, 'C^2=qI failed in the orbit audit')
    total_unordered = math.comb(q + 1, 3)
    signed_triangle_sum = int(np.trace(C @ C @ C)) // 6
    positive_unordered = (total_unordered + signed_triangle_sum) // 2
    negative_unordered = total_unordered - positive_unordered
    psl_order = q * (q * q - 1) // 2
    ordered_positive = 6 * positive_unordered
    _require(signed_triangle_sum == 0, 'positive and negative triangle counts differ')
    _require(
        ordered_positive == psl_order,
        'the positive ordered triples do not have PSL orbit size',
    )

    return {
        'q': q,
        'switching_formula': (
            'C[gP,gQ]=chi(det(g))*chi(lambda_P)*chi(lambda_Q)*C[P,Q]'
        ),
        'signed_vector_action_for_square_determinant': (
            '(U_g v)[gP]=chi(lambda_P)*v[P]'
        ),
        'audit_matrix': [list(_coordinates(p, entry)) for entry in matrix],
        'audit_matrix_determinant_character': determinant_character,
        'switching_factor_exact_on_all_edges': factor_exact,
        'positive_unordered_triangles': positive_unordered,
        'negative_unordered_triangles': negative_unordered,
        'ordered_positive_triangles': ordered_positive,
        'PSL_2_q_order': psl_order,
        'positive_support_triangles_single_PSL_orbit': True,
        'positive_signed_triples_single_signed_PSL_orbit': True,
        'nonsquare_scalar_lift_changes_all_switches': True,
        'triangular_completion_exists_for_every_signed_triple_datum': True,
        'all_completions_of_fixed_datum_classified': False,
        'proved_for_this_prime': True,
    }


def triangular_three_spike_certificate(p: int) -> dict[str, object]:
    """Run the complete exact construction audit at one representative prime."""
    _validate_prime(p)
    q = p * p
    n = q + 1
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    Q = C[1:, 1:]
    image_a, image_b = find_square_triangle_linear_map(p)
    image_diagonal = _add(p, image_a, image_b)
    directions = (image_a, image_b, image_diagonal)

    f0 = triangular_augmented_function(p)
    f = transform_function(p, f0, image_a, image_b)
    _require(int(f.sum()) == p * (p + 3) // 2, "the augmented mass changed")
    values, counts = np.unique(f, return_counts=True)
    value_counts = {int(value): int(count) for value, count in zip(values, counts)}
    # Keep the simple all-prime multiplicities explicit and independently checked.
    expected_zero = (p * p - 3 * p + 4) // 2
    expected_one = p * (p + 3) // 2 - 4
    _require(
        value_counts == {0: expected_zero, 1: expected_one, 2: 2},
        "the exact augmented-value multiplicities changed",
    )

    constant = (p + 3) // 2
    Q_identity = bool(
        np.array_equal(Q @ f, p * f - constant * np.ones(q, dtype=np.int64))
    )
    _require(Q_identity, "Qf=p f-(p+3)/2 failed")

    y = np.concatenate((np.array([3], dtype=np.int64), 2 * f - 1))
    source_u = _encode(p, 0, p - 2)
    source_v = _encode(p, p - 1, p - 2)
    target_u = _add(
        p,
        _scalar(p, 0, image_a),
        _scalar(p, p - 2, image_b),
    )
    target_v = _add(
        p,
        _scalar(p, p - 1, image_a),
        _scalar(p, p - 2, image_b),
    )
    _require(f0[source_u] == f0[source_v] == 2, "the two source overlaps changed")
    spikes = (0, 1 + target_u, 1 + target_v)
    actual_spikes = tuple(int(i) for i in np.flatnonzero(np.abs(y) == 3))
    _require(set(actual_spikes) == set(spikes), "the three-spike support changed")
    _require(all(int(y[i]) == 3 for i in spikes), "all spikes must be +3")

    eig_exact = bool(np.array_equal(C @ y, p * y))
    _require(eig_exact, "the integral vector is not an exact +p eigenvector")
    spike_edge_values = [int(C[i, j]) for i, j in combinations(spikes, 2)]
    _require(spike_edge_values == [1, 1, 1], "the spike triangle is not positive")

    shadow = y.copy()
    shadow[list(spikes)] = 1
    _require(bool(np.all(np.abs(shadow) == 1)), "the spike shadow is not Boolean")
    signed_triple = np.zeros(n, dtype=np.int64)
    signed_triple[list(spikes)] = -1
    _require(
        np.array_equal(y, shadow - 2 * signed_triple),
        "the shell relation y=x-2z failed",
    )
    _require(
        all(signed_triple[i] == -shadow[i] for i in spikes),
        "the three spike signs are not all bad",
    )

    Phi = p * n // 2
    q_shadow = int(shadow @ C @ shadow) // 2
    defect = Phi - q_shadow
    _require(defect == 6 * p - 12, "the Boolean shadow has the wrong defect")

    line_identities = [finite_line_identity(p, Q, direction) for direction in directions]
    _require(
        all(row["spatial_direction_character"] == 1 for row in line_identities),
        "T did not map every special spatial direction to square type",
    )
    fourier = fourier_direction_ledger(p, directions)
    KS = kiss_somlai_direction_audit(p)
    orbit = signed_pgl_triangle_transitivity(p)

    return {
        "p": p,
        "dimension": n,
        "linear_map_columns": [
            list(_coordinates(p, image_a)),
            list(_coordinates(p, image_b)),
        ],
        "mapped_special_spatial_directions": [
            list(_coordinates(p, direction)) for direction in directions
        ],
        "kiss_somlai_direction_audit": KS,
        "augmented_value_counts": value_counts,
        "augmented_mass": int(f.sum()),
        "Qf_identity_constant": constant,
        "Qf_equals_pf_minus_constant": Q_identity,
        "spike_indices": list(spikes),
        "spike_coordinates": [
            "infinity",
            list(_coordinates(p, target_u)),
            list(_coordinates(p, target_v)),
        ],
        "spike_edge_values": spike_edge_values,
        "C_y_equals_p_y": eig_exact,
        "shadow_is_boolean": True,
        "signed_triple_values_on_support": [-1, -1, -1],
        "all_three_shell_signs_bad": True,
        "boolean_shadow_q_C": q_shadow,
        "boolean_shadow_defect": defect,
        "line_direction_identities": line_identities,
        "fourier_direction_ledger": fourier,
        "signed_PGL_triangle_orbit": orbit,
        "classification": "all-prime construction/barrier for odd primes p>=5",
        "residual_ii_closed": False,
        "proved_for_this_prime": True,
    }


def representative_audit(primes: tuple[int, ...] = (5, 7, 11, 13)) -> dict[str, object]:
    """Exact fail-when-wrong checks; the mathematical proof is all-prime."""
    rows = {str(p): triangular_three_spike_certificate(p) for p in primes}
    return {
        "scope": "representative exact audit of an all-prime p>=5 construction",
        "rows": rows,
        "all_checks": all(row["proved_for_this_prime"] for row in rows.values()),
        "residual_ii_closed": False,
    }


if __name__ == "__main__":
    out = representative_audit()
    print("Kiss--Somlai triangular all-bad three-spike audit: exact")
    for p, row in out["rows"].items():
        print(
            f"  p={p}: T={row['linear_map_columns']} "
            f"spikes={row['spike_indices']} defect={row['boolean_shadow_defect']}"
        )
    print("  residual (ii): OPEN")
