#!/usr/bin/env python3
r"""Symbolic higher-moment audit for the compact rays of Proposition 15.758.

This module does not search primes or graphs.  It records two exact facts.

* Antipodal three-label atoms annihilate every odd member of the complete
  Proposition 15.759 hierarchy below top degree.  In particular degree five
  cannot exclude either compact ray.
* Separately at degrees six and eight, the whole ``p=1 mod 4`` compact ray
  and the lower endpoint of the ``p=3 mod 4`` ray have compatible atom labels
  for every live ``r>=7``.  Degree eight uses Cauchy--Davenport plus one
  order-four/eight cyclotomic number to obtain ``c+2a=1`` in nonzero eighth
  powers.  These separate degree-eight labels are not asserted to preserve
  the degree-six sums.
* On the first interior ``p=3 mod 4`` layers, the centered antipodal
  degree-six construction really does fail: a degree-six binary form is
  forced to have at least seven roots and also a prescribed nonzero value.
  This is an obstruction to that atom labelling, not to arbitrary labels.
* An exact 450-orbit certificate excludes every arbitrary compact atom plus
  six all-equal atoms when that row is required to vanish in all odd channels
  and in degrees six and eight.  The same profile occurs in every balanced
  ``p=31`` branch-C allocation from ``t=69`` through ``t=99``.
* The odd/Radon support theorem extends uniformly to every branch-C row with
  ``b`` arbitrary compact atoms and ``r-1`` all-equal atoms when
  ``p=4r+3``, ``r>=7``, and ``3b<=r+2``.  If all odd forms vanish, its
  aggregate edge chain is centrally symmetric.  This covers an explicit
  initial band of every balanced branch-C ray.
* The unrestricted joint degree-six/eight maps of four arbitrary compact
  atoms and of four all-equal atoms are both dominant.  Consequently there
  is no universal algebraic seven-channel identity extending the centered
  root obstruction; finite-field rationality and the other moment blocks
  remain essential.

No simultaneous degree-six/eight, integral, or Boolean common edge lift is
asserted.  The joint even hierarchy remains open at both all-orders fronts.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd, isqrt

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15758 import p1_local_survivor, p3_local_survivor


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 7
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=7")


def moment_polynomial(p: int, s: int, t: int, d: int, k: int) -> int:
    """Evaluate ``(s-t)^2(st)^k(s+t)^(d-2-2k)`` modulo ``p``."""
    _check_prime(p)
    if not 2 <= d <= p - 1 or not 0 <= k < d // 2:
        raise ValueError("invalid moment index")
    return (
        pow((s - t) % p, 2, p)
        * pow((s * t) % p, k, p)
        * pow((s + t) % p, d - 2 - 2 * k, p)
    ) % p


def star_moment(p: int, centre: int, d: int, k: int) -> int:
    """Moment of the unit star at ``centre``."""
    return sum(moment_polynomial(p, centre, t, d, k) for t in range(p)) % p


def compact_moment(
    p: int, a: int, b: int, distinguished: int, d: int, k: int
) -> int:
    """Moment of ``+ab-a*distinguished-b*distinguished``."""
    return (
        moment_polynomial(p, a, b, d, k)
        - moment_polynomial(p, a, distinguished, d, k)
        - moment_polynomial(p, b, distinguished, d, k)
    ) % p


def all_equal_moment(p: int, a: int, b: int, c: int, d: int, k: int) -> int:
    """Moment of the all-positive triangle on ``a,b,c``."""
    return (
        moment_polynomial(p, a, b, d, k)
        + moment_polynomial(p, a, c, d, k)
        + moment_polynomial(p, b, c, d, k)
    ) % p


def _dual_constant(value: int, size: int) -> tuple[int, tuple[int, ...]]:
    return value, (0,) * size


def _dual_add(
    left: tuple[int, tuple[int, ...]],
    right: tuple[int, tuple[int, ...]],
) -> tuple[int, tuple[int, ...]]:
    return (
        left[0] + right[0],
        tuple(a + b for a, b in zip(left[1], right[1])),
    )


def _dual_scale(
    value: tuple[int, tuple[int, ...]], scalar: int
) -> tuple[int, tuple[int, ...]]:
    return scalar * value[0], tuple(scalar * item for item in value[1])


def _dual_multiply(
    left: tuple[int, tuple[int, ...]],
    right: tuple[int, tuple[int, ...]],
) -> tuple[int, tuple[int, ...]]:
    return (
        left[0] * right[0],
        tuple(
            left[0] * b + right[0] * a
            for a, b in zip(left[1], right[1])
        ),
    )


def _dual_power(
    value: tuple[int, tuple[int, ...]], exponent: int
) -> tuple[int, tuple[int, ...]]:
    if exponent < 0:
        raise ValueError("negative dual exponent")
    out = _dual_constant(1, len(value[1]))
    base = value
    power = exponent
    while power:
        if power & 1:
            out = _dual_multiply(out, base)
        base = _dual_multiply(base, base)
        power >>= 1
    return out


def _dual_moment_polynomial(
    s: tuple[int, tuple[int, ...]],
    t: tuple[int, tuple[int, ...]],
    degree: int,
    channel: int,
) -> tuple[int, tuple[int, ...]]:
    difference = _dual_add(s, _dual_scale(t, -1))
    product_value = _dual_multiply(s, t)
    total = _dual_add(s, t)
    return _dual_multiply(
        _dual_multiply(
            _dual_power(difference, 2),
            _dual_power(product_value, channel),
        ),
        _dual_power(total, degree - 2 - 2 * channel),
    )


def _bareiss_determinant(matrix: list[list[int]]) -> int:
    """Return an exact integer determinant by fraction-free elimination."""
    values = [row[:] for row in matrix]
    size = len(values)
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if values[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            sign *= -1
        pivot_value = values[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    values[row][index] * pivot_value
                    - values[row][column] * values[column][index]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss division stopped being exact")
                values[row][index] = numerator // previous
        previous = pivot_value
        for row in range(column + 1, size):
            values[row][column] = 0
    return sign * values[-1][-1]


def _joint_even_four_atom_jacobian(atom_type: str) -> list[list[int]]:
    """Evaluate the seven-channel integer Jacobian at four fixed triples."""
    variable_count = 8
    pairs = ((2, 1), (3, 2), (4, 3), (5, 4))
    atoms = []
    for index, (a, b) in enumerate(pairs):
        a_gradient = [0] * variable_count
        b_gradient = [0] * variable_count
        a_gradient[2 * index] = 1
        b_gradient[2 * index + 1] = 1
        atoms.append(
            (
                (a, tuple(a_gradient)),
                (b, tuple(b_gradient)),
                _dual_constant(0, variable_count),
            )
        )

    rows: list[list[int]] = []
    for degree in (6, 8):
        for channel in range(degree // 2):
            output = _dual_constant(0, variable_count)
            for a, b, c in atoms:
                ab = _dual_moment_polynomial(a, b, degree, channel)
                ac = _dual_moment_polynomial(a, c, degree, channel)
                bc = _dual_moment_polynomial(b, c, degree, channel)
                if atom_type == "compact":
                    atom = _dual_add(ab, _dual_scale(_dual_add(ac, bc), -1))
                elif atom_type == "all_equal":
                    atom = _dual_add(ab, _dual_add(ac, bc))
                else:
                    raise ValueError("unknown atom type")
                output = _dual_add(output, atom)
            rows.append(list(output[1]))
    return rows


def joint_six_eight_atom_map_dominance_certificate() -> dict[str, object]:
    r"""Rule out a universal algebraic relation among the seven even channels.

    Vary four atoms on ``(a_i,b_i,0)`` and sum the three degree-six and four
    degree-eight contractions.  At

    ``((2,1),(3,2),(4,3),(5,4))``

    the Jacobian minor in variables
    ``(b_1,a_2,b_2,a_3,b_3,a_4,b_4)`` has the two displayed nonzero integer
    determinants.  Thus both the compact and all-equal four-atom maps are
    dominant on the distinct-label locus in every characteristic at least
    eleven.

    Every balanced branch-C hard row for ``r>=7`` has at least four compact
    atoms, and every opposite row has at least six all-equal atoms.  Freezing
    the extra atoms preserves dominance.  Over the algebraic closure, the
    dense-open images of the finitely many signed row types have a common
    point.  Scaling a base configuration by a linear form with no zero on
    ``P^1(F_p)`` then gives nonzero common degree-six/eight binary forms.
    This refutes a purely algebraic projective root-count closure of this
    seven-channel subsystem.  It does not provide labels or form
    coefficients in ``F_p`` and says nothing about odd/higher moments or the
    Boolean lift.
    """
    compact_jacobian = _joint_even_four_atom_jacobian("compact")
    all_equal_jacobian = _joint_even_four_atom_jacobian("all_equal")
    # Delete the first column, corresponding to a_1.
    compact_determinant = _bareiss_determinant(
        [row[1:] for row in compact_jacobian]
    )
    all_equal_determinant = _bareiss_determinant(
        [row[1:] for row in all_equal_jacobian]
    )
    expected_compact = 2**28 * 3**9 * 5**3 * 7**3
    expected_all_equal = 2**26 * 3**7 * 5**4 * 7**4
    proved = bool(
        compact_determinant == expected_compact
        and all_equal_determinant == expected_all_equal
    )
    if not proved:
        raise ArithmeticError("the joint degree-six/eight Jacobian changed")
    return {
        "channel_order": [
            "(6,0)",
            "(6,1)",
            "(6,2)",
            "(8,0)",
            "(8,1)",
            "(8,2)",
            "(8,3)",
        ],
        "base_triples": [[2, 1, 0], [3, 2, 0], [4, 3, 0], [5, 4, 0]],
        "minor_variables": ["b_1", "a_2", "b_2", "a_3", "b_3", "a_4", "b_4"],
        "compact_jacobian_determinant": compact_determinant,
        "compact_determinant_factorization": "2^28*3^9*5^3*7^3",
        "all_equal_jacobian_determinant": all_equal_determinant,
        "all_equal_determinant_factorization": "2^26*3^7*5^4*7^4",
        "full_rank_for_every_characteristic_at_least": 11,
        "branch_C_minimum_hard_compact_atoms": 4,
        "branch_C_minimum_opposite_all_equal_atoms": 6,
        "universal_polynomial_relation_among_seven_channels": False,
        "algebraic_closure_common_form_obstruction": False,
        "F_p_rational_common_forms_constructed": False,
        "odd_or_higher_moment_compatibility_proved": False,
        "Boolean_lift_constructed": False,
        "proved": proved,
    }


def omitted_pair_moment(p: int, a: int, b: int, d: int, k: int) -> int:
    """Moment of the canonical offset-minus-one omitted-pair row.

    Its coefficient graph is ``-S_a-S_b+{a,b}``.
    """
    return (
        -star_moment(p, a, d, k)
        - star_moment(p, b, d, k)
        + moment_polynomial(p, a, b, d, k)
    ) % p


def odd_antipodal_atom_certificate(p: int, scale: int = 1) -> dict[str, object]:
    """Replay the rowwise annihilation of all odd degrees below ``p-1``."""
    _check_prime(p)
    scale %= p
    if scale == 0:
        raise ValueError("the antipodal scale must be nonzero")
    degrees: dict[str, dict[str, object]] = {}
    for d in range(3, p - 1, 2):
        compact = tuple(
            compact_moment(p, scale, -scale, 0, d, k)
            for k in range(d // 2)
        )
        all_equal = tuple(
            all_equal_moment(p, scale, -scale, 0, d, k)
            for k in range(d // 2)
        )
        omitted = tuple(
            omitted_pair_moment(p, scale, -scale, d, k)
            for k in range(d // 2)
        )
        stars = tuple(star_moment(p, scale, d, k) for k in range(d // 2))
        proved = not any(compact + all_equal + omitted + stars)
        if not proved:
            raise ArithmeticError("an antipodal odd-moment identity changed")
        degrees[str(d)] = {
            "compact": list(compact),
            "all_equal": list(all_equal),
            "omitted_pair": list(omitted),
            "star": list(stars),
            "proved": proved,
        }
    return {
        "p": p,
        "labels": [scale % p, (-scale) % p, 0],
        "odd_degrees_checked": list(range(3, p - 1, 2)),
        "degree_five_is_rowwise_zero": degrees.get("5", {}).get("proved", False),
        "degree_rows": degrees,
        "all_odd_moments_below_top_are_zero": True,
        "proved": True,
    }


def degree_six_antipodal_vectors(p: int, scale: int = 1) -> dict[str, object]:
    """Return the three degree-six contractions of the antipodal atoms."""
    _check_prime(p)
    scale %= p
    if scale == 0:
        raise ValueError("the antipodal scale must be nonzero")
    compact = tuple(
        compact_moment(p, scale, -scale, 0, 6, k) for k in range(3)
    )
    all_equal = tuple(
        all_equal_moment(p, scale, -scale, 0, 6, k) for k in range(3)
    )
    omitted = tuple(
        omitted_pair_moment(p, scale, -scale, 6, k) for k in range(3)
    )
    a6 = pow(scale, 6, p)
    expected = {
        "compact": ((-2 * a6) % p, 0, (4 * a6) % p),
        "all_equal": ((2 * a6) % p, 0, (4 * a6) % p),
        "omitted_pair": (0, 0, (4 * a6) % p),
    }
    proved = compact == expected["compact"] and all_equal == expected["all_equal"] and omitted == expected["omitted_pair"]
    if not proved:
        raise ArithmeticError("the antipodal degree-six vectors changed")
    return {
        "p": p,
        "scale_sixth_power": a6,
        "compact": list(compact),
        "all_equal": list(all_equal),
        "omitted_pair": list(omitted),
        "middle_channel_zero": True,
        "proved": proved,
    }


def degree_eight_antipodal_vectors(p: int, scale: int = 1) -> dict[str, object]:
    """Return the four degree-eight contractions of the antipodal atoms."""
    _check_prime(p)
    scale %= p
    if scale == 0:
        raise ValueError("the antipodal scale must be nonzero")
    compact = tuple(
        compact_moment(p, scale, -scale, 0, 8, k) for k in range(4)
    )
    all_equal = tuple(
        all_equal_moment(p, scale, -scale, 0, 8, k) for k in range(4)
    )
    omitted = tuple(
        omitted_pair_moment(p, scale, -scale, 8, k) for k in range(4)
    )
    a8 = pow(scale, 8, p)
    expected = {
        "compact": ((-2 * a8) % p, 0, 0, (-4 * a8) % p),
        "all_equal": ((2 * a8) % p, 0, 0, (-4 * a8) % p),
        "omitted_pair": (0, 0, 0, (-4 * a8) % p),
    }
    proved = (
        compact == expected["compact"]
        and all_equal == expected["all_equal"]
        and omitted == expected["omitted_pair"]
    )
    if not proved:
        raise ArithmeticError("the antipodal degree-eight vectors changed")
    return {
        "p": p,
        "scale_eighth_power": a8,
        "compact": list(compact),
        "all_equal": list(all_equal),
        "omitted_pair": list(omitted),
        "two_middle_channels_zero": True,
        "proved": proved,
    }


def sixth_power_sumset_bound(p: int) -> dict[str, object]:
    """Cauchy--Davenport threshold for exact sums of nonzero sixth powers."""
    _check_prime(p)
    divisor = gcd(6, p - 1)
    size = (p - 1) // divisor
    threshold = (p - 2 + (size - 1)) // (size - 1)
    lower_at_threshold = min(p, threshold * size - (threshold - 1))
    proved = lower_at_threshold == p
    if not proved:
        raise ArithmeticError("the Cauchy--Davenport threshold changed")
    return {
        "p": p,
        "gcd_6_p_minus_1": divisor,
        "nonzero_sixth_power_set_size": size,
        "cauchy_davenport_threshold": threshold,
        "lower_bound_at_threshold": lower_at_threshold,
        "every_field_element_is_an_exact_N_term_sum_for_every_N_at_least_threshold": True,
        "proved": proved,
    }


def eighth_power_sumset_bound(p: int) -> dict[str, object]:
    """Cauchy--Davenport threshold for exact sums of nonzero eighth powers."""
    _check_prime(p)
    divisor = gcd(8, p - 1)
    size = (p - 1) // divisor
    threshold = (p - 2 + (size - 1)) // (size - 1)
    lower_at_threshold = min(p, threshold * size - (threshold - 1))
    proved = lower_at_threshold == p
    if not proved:
        raise ArithmeticError("the Cauchy--Davenport threshold changed")
    return {
        "p": p,
        "gcd_8_p_minus_1": divisor,
        "nonzero_eighth_power_set_size": size,
        "cauchy_davenport_threshold": threshold,
        "lower_bound_at_threshold": lower_at_threshold,
        "every_field_element_is_an_exact_N_term_sum_for_every_N_at_least_threshold": True,
        "proved": proved,
    }


_P37_BASES: dict[tuple[int, int], tuple[int, ...]] = {
    (6, 0): (1, 1, 1, 36, 36, 36),
    (7, 0): (1, 1, 1, 1, 11, 11, 11),
    (6, 18): (1, 1, 1, 26, 27, 36),
    (7, 18): (1, 1, 1, 1, 26, 26, 36),
    (6, 19): (1, 1, 1, 1, 26, 26),
    (7, 19): (1, 1, 1, 1, 26, 27, 36),
}


def p37_sixth_power_witness(term_count: int, target: int) -> tuple[int, ...]:
    """Exact non-census witnesses for the sole Cauchy--Davenport exception."""
    if term_count < 6 or target % 37 not in (0, 18, 19):
        raise ValueError("the p=37 witness covers N>=6 and targets 0,18,19")
    parity_base = 6 if term_count % 2 == 0 else 7
    values = list(_P37_BASES[(parity_base, target % 37)])
    while len(values) < term_count:
        values.extend((1, 36))
    if len(values) != term_count:
        raise ArithmeticError("the p=37 parity padding changed")
    residues = {pow(value, 6, 37) for value in range(1, 37)}
    if not all(value in residues for value in values) or sum(values) % 37 != target % 37:
        raise ArithmeticError("a displayed p=37 sixth-power identity failed")
    return tuple(values)


def _p1_sixth_power_existence(p: int, counts: tuple[int, ...]) -> dict[str, object]:
    r = (p - 1) // 4
    if p == 37:
        for count in counts:
            p37_sixth_power_witness(count, 0)
            p37_sixth_power_witness(count, 18)
        p37_sixth_power_witness(r - 2, 18)
        p37_sixth_power_witness(r - 2, 19)
        return {"method": "displayed p=37 identities", "proved": True}
    bound = sixth_power_sumset_bound(p)
    threshold = int(bound["cauchy_davenport_threshold"])
    proved = min(counts) >= threshold and r - 2 >= threshold
    if not proved:
        raise ArithmeticError("the sixth-power sumset bound does not cover this ray")
    return {"method": "Cauchy-Davenport", "threshold": threshold, "proved": True}


def p1_degree_six_ray_certificate(p: int, t: int) -> dict[str, object]:
    """Degree-six compatibility for the full ``p=4r+1`` compact ray.

    A projective linear form ``s_L=L(v)`` is chosen with its single zero in
    a hard direction.  All atoms use antipodal labels.  If ``b_L=Q_L-r`` is
    the number of opposite compact atoms, choose their sixth-power sum to be
    zero for even ``b_L`` and ``s_L^6`` for odd ``b_L``.  The omitted pair
    has sixth power ``s_L^6`` or ``-s_L^6`` respectively.  The all-equal
    atoms have sum ``-c_L^6/2`` and the hard compact atoms have sum
    ``-s_L^6/2``.  Every direction then evaluates the same two binary forms
    ``F_0(L)=s_L^6`` and ``F_2(L)=-2s_L^6``; the middle form is zero.
    """
    row = p1_local_survivor(p, t)
    r = int(row["r"])
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    hard_counts = tuple(int(item["e"]) for item in row["hard_rows"])
    compact_counts = tuple(int(item["Q"]) - r for item in row["opposite_rows"])
    existence = _p1_sixth_power_existence(p, hard_counts)
    minus_one_is_sixth = pow(-1, (p - 1) // gcd(6, p - 1), p) == 1
    # The exponent check above is equivalent to -1 lying in the subgroup of
    # sixth powers.  It holds because p=1 mod 4.
    compact_parity_constructible = minus_one_is_sixth and all(count >= 0 for count in compact_counts)
    proved = bool(existence["proved"] and compact_parity_constructible)
    if not proved:
        raise ArithmeticError("the p=1 mod 4 degree-six construction failed")
    return {
        "p": p,
        "r": r,
        "t": t,
        "hard_compact_atom_counts": list(hard_counts),
        "opposite_all_equal_count_per_direction": r - 2,
        "opposite_compact_atom_counts": list(compact_counts),
        "one_omitted_pair_per_opposite_direction": True,
        "minus_one_is_a_nonzero_sixth_power": minus_one_is_sixth,
        "sixth_power_existence": existence,
        "global_degree_six_forms": {
            "k_0": "L(v)^6",
            "k_1": "0",
            "k_2": "-2*L(v)^6",
        },
        "all_degree_five_rows_zero_simultaneously": True,
        "all_degree_six_moment_relations_pass": True,
        "proved": proved,
    }


def _sum_of_two_squares_parameter(p: int) -> int:
    """Return the uniquely signed ``s=1 mod 4`` in ``p=s^2+t^2``."""
    for first in range(isqrt(p) + 1):
        second_square = p - first * first
        second = isqrt(second_square)
        if second * second != second_square:
            continue
        for candidate in (first, second):
            if candidate and candidate % 2:
                return candidate if candidate % 4 == 1 else -candidate
    raise ArithmeticError("the sum-of-two-squares parameter was not found")


def _x_squared_plus_four_y_squared_parameter(p: int) -> int:
    """Return the uniquely signed ``x=1 mod 4`` in ``p=x^2+4y^2``."""
    for y_value in range(isqrt(p // 4) + 1):
        x_square = p - 4 * y_value * y_value
        x_value = isqrt(x_square)
        if x_value * x_value == x_square and x_value:
            return x_value if x_value % 4 == 1 else -x_value
    raise ArithmeticError("the x^2+4y^2 parameter was not found")


def p1_eighth_power_affine_pair_lemma(p: int) -> dict[str, object]:
    r"""Prove that ``c+2a=1`` has ``a,c`` nonzero eighth powers.

    Put ``H=(F_p^*)^8``.  For ``p=1 mod 16`` the pair is simply
    ``(a,c)=(1,-1)``.  In the other two residue classes the number of
    solutions is a classical order-four or order-eight cyclotomic number.
    The formulas are Theorems 10.2 and 10.4 of Huczynska--Johnson,
    arXiv:2201.07553 (Theorems 7.2 and 7.4 in the arXiv version).

    A cyclotomic number ``(i,0)`` counts ``x in (-2)H`` for which
    ``1+x in H``.  Taking ``x=-2a`` gives the asserted pair.
    """
    _check_prime(p)
    if p % 4 != 1:
        raise ValueError("need p=1 mod 4")
    r = (p - 1) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")

    explicit_pair: tuple[int, int] | None = None
    if p % 16 == 1:
        # Here 8 divides (p-1)/2, so -1 itself is an eighth power.
        explicit_pair = (1, p - 1)
        count: int | None = None
        parameter: dict[str, object] = {}
        method = "explicit (a,c)=(1,-1)"
        positivity = True
    elif p % 8 == 5:
        # H is the quartic-residue subgroup.  Since 2 is a nonsquare and
        # -1 is a square, -2 belongs to order-four class 1 or 3.  For odd
        # f=(p-1)/4 both relevant cyclotomic numbers have this value.
        s_value = _sum_of_two_squares_parameter(p)
        count = (p - 3 - 2 * s_value) // 16
        parameter = {"s": s_value, "identity": f"{p}=({s_value})^2+t^2"}
        method = "order-four cyclotomic number (1,0)=(3,0)"
        positivity = count > 0 and p - 3 - 2 * isqrt(p) > 0
    else:
        # Now p=9 mod 16 and f=(p-1)/8 is odd.  The class of -1 is 4.
        # If 2 is quartic, -2 is in class 0 or 4; otherwise it is in class
        # 2 or 6.  The corresponding pairs of order-eight numbers coincide.
        x_value = _x_squared_plus_four_y_squared_parameter(p)
        two_is_quartic = pow(2, (p - 1) // 4, p) == 1
        if two_is_quartic:
            count = (p - 15 - 2 * x_value) // 64
            method = "order-eight cyclotomic number (0,0)=(4,0)"
            positivity = count > 0 and p - 15 - 2 * isqrt(p) > 0
        else:
            count = (p - 7 + 6 * x_value) // 64
            method = "order-eight cyclotomic number (2,0)=(6,0)"
            positivity = count > 0 and (
                p == 41 or p - 7 - 6 * isqrt(p) > 0
            )
        parameter = {
            "x": x_value,
            "identity": f"{p}=({x_value})^2+4*y^2",
            "two_is_a_quartic_residue": two_is_quartic,
        }

    proved = bool(positivity)
    if explicit_pair is not None:
        a_value, c_value = explicit_pair
        proved = proved and (c_value + 2 * a_value) % p == 1
        proved = proved and pow(a_value, (p - 1) // gcd(8, p - 1), p) == 1
        proved = proved and pow(c_value, (p - 1) // gcd(8, p - 1), p) == 1
    if not proved:
        raise ArithmeticError("the eighth-power affine-pair lemma failed")
    return {
        "p": p,
        "r": r,
        "equation": "c+2*a=1 with a,c in (F_p^*)^8",
        "method": method,
        "cyclotomic_solution_count": count,
        "representation_parameter": parameter,
        "explicit_pair_a_c": list(explicit_pair) if explicit_pair else None,
        "source": "Huczynska--Johnson, arXiv:2201.07553, cyclotomic-number appendix",
        "proved": proved,
    }


def eighth_power_geometric_allocation(
    p: int, compact_count: int, a_value: int, c_value: int
) -> dict[str, object]:
    r"""Lift one ``c+2a=1`` pair to every exact compact count.

    The identity

    ``1=c^b+2*a*(1+c+...+c^(b-1))``

    uses one eighth power for the omitted pair and exactly ``b`` eighth
    powers for the compact atoms.  Multiplication by any eighth power keeps
    every summand in the same subgroup.
    """
    _check_prime(p)
    if not isinstance(compact_count, int) or isinstance(compact_count, bool) or compact_count < 0:
        raise ValueError("compact_count must be a nonnegative integer")
    a_value %= p
    c_value %= p
    exponent = (p - 1) // gcd(8, p - 1)
    if (
        not a_value
        or not c_value
        or pow(a_value, exponent, p) != 1
        or pow(c_value, exponent, p) != 1
        or (c_value + 2 * a_value) % p != 1
    ):
        raise ValueError("need nonzero eighth powers satisfying c+2a=1")
    compact_terms = tuple(
        a_value * pow(c_value, index, p) % p for index in range(compact_count)
    )
    omitted_value = pow(c_value, compact_count, p)
    proved = (omitted_value + 2 * sum(compact_terms)) % p == 1
    if not proved:
        raise ArithmeticError("the geometric eighth-power identity failed")
    return {
        "p": p,
        "compact_count": compact_count,
        "compact_eighth_power_terms": list(compact_terms),
        "omitted_pair_eighth_power": omitted_value,
        "identity_value": (omitted_value + 2 * sum(compact_terms)) % p,
        "proved": proved,
    }


_P1_EIGHTH_MINUS_HALF_SPECIALS: dict[tuple[int, int], tuple[int, ...]] = {
    (29, 4): (1, 1, 16, 25),
    (41, 7): (1, 1, 1, 1, 10, 10, 37),
    (41, 8): (1, 1, 1, 1, 1, 1, 18, 37),
    (41, 9): (1, 1, 1, 1, 1, 10, 10, 18, 18),
}

_P1_EIGHTH_ZERO_SPECIALS: dict[tuple[int, int], tuple[int, ...]] = {
    (29, 4): (1, 1, 7, 20),
    (41, 7): (1, 1, 1, 1, 1, 18, 18),
    (41, 8): (1, 1, 1, 1, 1, 10, 10, 16),
    (41, 9): (1, 1, 1, 1, 1, 1, 1, 16, 18),
}


def _p1_eighth_power_existence(p: int, counts: tuple[int, ...]) -> dict[str, object]:
    """Check the exact ``0`` and ``-1/2`` sums needed in branch B."""
    minus_half = (-pow(2, -1, p)) % p
    bound = eighth_power_sumset_bound(p)
    threshold = int(bound["cauchy_davenport_threshold"])
    special_witnesses: dict[str, dict[str, list[int]]] = {"0": {}, "-1/2": {}}
    for count in sorted(set(counts)):
        if count >= threshold:
            continue
        residue_set = {pow(value, 8, p) for value in range(1, p)}
        for label, target, table in (
            ("0", 0, _P1_EIGHTH_ZERO_SPECIALS),
            ("-1/2", minus_half, _P1_EIGHTH_MINUS_HALF_SPECIALS),
        ):
            witness = table.get((p, count))
            if witness is None:
                raise ArithmeticError("an eighth-power exact-count exception is uncovered")
            if len(witness) != count or not all(value in residue_set for value in witness):
                raise ArithmeticError("a displayed eighth-power witness is not in H_8")
            if sum(witness) % p != target:
                raise ArithmeticError("a displayed eighth-power sum changed")
            special_witnesses[label][str(count)] = list(witness)
    return {
        "targets": {"zero_at_the_L_root": 0, "normalized_nonzero": minus_half},
        "cauchy_davenport_threshold": threshold,
        "displayed_short_witnesses": special_witnesses,
        "all_requested_exact_counts_covered": True,
        "proved": True,
    }


def p1_degree_eight_ray_certificate(p: int, t: int) -> dict[str, object]:
    r"""Degree-eight compatibility for the full ``p=4r+1`` compact ray."""
    row = p1_local_survivor(p, t)
    r = int(row["r"])
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    hard_counts = tuple(int(item["e"]) for item in row["hard_rows"])
    compact_counts = tuple(int(item["Q"]) - r for item in row["opposite_rows"])
    affine_pair = p1_eighth_power_affine_pair_lemma(p)
    exact_sums = _p1_eighth_power_existence(p, hard_counts + (r - 2,))
    proved = bool(
        affine_pair["proved"]
        and exact_sums["proved"]
        and all(count >= 0 for count in compact_counts)
    )
    if not proved:
        raise ArithmeticError("the p=1 mod 4 degree-eight construction failed")
    return {
        "p": p,
        "r": r,
        "t": t,
        "hard_compact_atom_counts": list(hard_counts),
        "opposite_all_equal_count_per_direction": r - 2,
        "opposite_compact_atom_counts": list(compact_counts),
        "one_omitted_pair_per_opposite_direction": True,
        "eighth_power_affine_pair": affine_pair,
        "eighth_power_exact_sum_existence": exact_sums,
        "opposite_geometric_identity": (
            "C=c^b*L(v)^8; B=sum_{j=0}^{b-1} a*c^j*L(v)^8; C+2B=L(v)^8"
        ),
        "global_degree_eight_forms": {
            "k_0": "L(v)^8",
            "k_1": "0",
            "k_2": "0",
            "k_3": "2*L(v)^8",
        },
        "standalone_degree_eight_only": True,
        "same_labels_as_degree_six_not_proved": True,
        "all_degree_eight_moment_relations_pass": True,
        "proved": proved,
    }


_P31_COMPACT_FOUR = (
    (0, 1, 8),
    (0, 30, 23),
    (1, 5, 14),
    (30, 26, 17),
)
_P31_COMPACT_FIVE_SCALES = (1, 2, 4, 8, 3)
_P31_ALL_EQUAL_SIX_SCALES = (1, 1, 1, 4, 8, 3)


def p31_exceptional_lower_endpoint_certificate() -> dict[str, object]:
    """Verify the explicit ``p=31`` degree-five/six zero identities."""
    p = 31
    compact_four = {
        d: [
            sum(compact_moment(p, a, b, c, d, k) for a, b, c in _P31_COMPACT_FOUR) % p
            for k in range(d // 2)
        ]
        for d in (5, 6)
    }
    compact_five = {
        d: [
            sum(compact_moment(p, scale, -scale, 0, d, k) for scale in _P31_COMPACT_FIVE_SCALES) % p
            for k in range(d // 2)
        ]
        for d in (5, 6)
    }
    all_equal_six = {
        d: [
            sum(all_equal_moment(p, scale, -scale, 0, d, k) for scale in _P31_ALL_EQUAL_SIX_SCALES) % p
            for k in range(d // 2)
        ]
        for d in (5, 6)
    }
    proved = not any(
        value
        for family in (compact_four, compact_five, all_equal_six)
        for values in family.values()
        for value in values
    )
    if not proved:
        raise ArithmeticError("a displayed p=31 atom identity failed")
    return {
        "p": p,
        "four_compact_triples": [list(row) for row in _P31_COMPACT_FOUR],
        "five_compact_antipodal_scales": list(_P31_COMPACT_FIVE_SCALES),
        "six_all_equal_antipodal_scales": list(_P31_ALL_EQUAL_SIX_SCALES),
        "compact_four_moments": compact_four,
        "compact_five_moments": compact_five,
        "all_equal_six_moments": all_equal_six,
        "proved": proved,
    }


def p3_lower_endpoint_degree_six_certificate(p: int) -> dict[str, object]:
    """Degree-five/six compatibility at the lower ``p=4r+3`` endpoint."""
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    t = 2 * r * r - 4 * r - 2
    row = p3_local_survivor(p, t)
    hard_counts = tuple(int(item["e"]) for item in row["hard_rows"])
    opposite_count = r - 1
    if p == 31:
        special = p31_exceptional_lower_endpoint_certificate()
        proved = set(hard_counts) == {4, 5} and opposite_count == 6 and special["proved"]
        existence: dict[str, object] = {
            "method": "displayed p=31 identities",
            "special": special,
            "proved": proved,
        }
    else:
        bound = sixth_power_sumset_bound(p)
        threshold = int(bound["cauchy_davenport_threshold"])
        proved = min(hard_counts) >= threshold and opposite_count >= threshold
        existence = {
            "method": "Cauchy-Davenport",
            "threshold": threshold,
            "proved": proved,
        }
    if not proved:
        raise ArithmeticError("the p=3 mod 4 lower endpoint escaped the sumset proof")
    return {
        "p": p,
        "r": r,
        "t_lower_endpoint": t,
        "hard_compact_atom_counts": list(hard_counts),
        "opposite_all_equal_count_per_direction": opposite_count,
        "opposite_compact_atom_count": 0,
        "sixth_power_existence": existence,
        "global_degree_six_forms": {"k_0": "0", "k_1": "0", "k_2": "0"},
        "all_degree_five_rows_zero_simultaneously": True,
        "all_degree_six_moment_relations_pass": True,
        "proved": proved,
    }


def p3_centered_degree_six_interior_gate(p: int, t: int) -> dict[str, object]:
    r"""Exact polynomial gate for centered antipodal labels in branch C.

    Let ``F_0,F_2`` be the global degree-six forms.  Hard compact atoms have
    ``F_2+2F_0=0``.  There are ``m=2r+2>6`` hard directions, hence the binary
    form ``H=F_2+2F_0`` is identically zero.  On an opposite row, if ``A`` is
    the all-equal sixth-power sum and ``B`` the compact sum, the outer sign
    gives

    ``F_0=-2A+2B, F_2=-4A-4B, H=-8A``.

    Thus ``A=0``.  A row with compact count ``b=0`` is a root of ``F_0``;
    a row with ``b=1`` has ``F_0=2a^6 != 0``.  Seven zero rows therefore
    contradict one unit row.  Arbitrary triangle labels need not obey this
    centered relation; ``p31_first_interior_odd_six_certificate`` supplies
    the explicit escape at the first parameter.
    """
    row = p3_local_survivor(p, t)
    r = int(row["r"])
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    m = 2 * r + 2
    lower = 2 * r * r - 4 * r - 2
    delta = t - lower
    compact_counts = tuple(
        int(item["Q"]) - r - 2 for item in row["opposite_rows"]
    )
    zero_rows = compact_counts.count(0)
    unit_rows = compact_counts.count(1)
    root_count_obstruction = zero_rows >= 7 and unit_rows >= 1
    return {
        "p": p,
        "r": r,
        "t": t,
        "delta_from_lower_endpoint": delta,
        "balanced_opposite_compact_counts": list(compact_counts),
        "degree_six_form": "H=F_{6,2}+2*F_{6,0}",
        "hard_identity": "H=0 in 2r+2>6 directions, hence H is identically zero",
        "opposite_identity": "H=-8*A, hence the all-equal sixth-power sum A is zero",
        "forced_F0_zero_rows": zero_rows,
        "forced_F0_nonzero_unit_rows": unit_rows,
        "degree_six_root_bound": 6,
        "centered_antipodal_construction_obstructed": root_count_obstruction,
        "arbitrary_triangle_labels_obstructed": False,
        "proved": True,
    }


def normalized_compact_degree_six_h(p: int, parameter: int) -> dict[str, object]:
    r"""Verify the noncentered compact correction to ``H=F_2+2F_0``.

    For the compact atom on ``(parameter,1,0)``, direct factorization gives

    ``H=t(t+1)^2(4t^2-9t+4)``.

    Thus the centered value ``t=-1`` is a special zero and its root-count
    relation cannot be imposed on arbitrary compact triangles.
    """
    _check_prime(p)
    parameter %= p
    vector = tuple(compact_moment(p, parameter, 1, 0, 6, k) for k in range(3))
    lhs = (vector[2] + 2 * vector[0]) % p
    rhs = (
        parameter
        * pow(parameter + 1, 2, p)
        * (4 * parameter * parameter - 9 * parameter + 4)
    ) % p
    proved = lhs == rhs
    if not proved:
        raise ArithmeticError("the normalized compact H factorization changed")
    return {
        "p": p,
        "parameter": parameter,
        "degree_six_vector": list(vector),
        "H_value": lhs,
        "factorization": "t*(t+1)^2*(4*t^2-9*t+4)",
        "proved": proved,
    }


_P31_FIRST_INTERIOR_COMPACT = (1, 30, 0)
_P31_FIRST_INTERIOR_ALL_EQUAL = (
    (0, 1, 2),
    (0, 30, 29),
    (0, 1, 24),
    (0, 30, 7),
    (0, 4, 23),
    (0, 27, 8),
)


def p31_first_interior_odd_six_certificate() -> dict[str, object]:
    """Escape the centered ``p=31,t=69`` gate with arbitrary AE labels.

    The compact atom is centered, so every odd contraction vanishes.  The
    six all-equal triangles occur in three pairs ``T,-T``, so their odd
    contractions cancel in every odd degree.  Their degree-six vector plus
    the compact vector is zero.  Repetition/overlap is permitted because the
    Proposition 15.758 row is a sum of nonnegative integral sharp atoms.
    """
    p = 31
    t = 69
    row = p3_local_survivor(p, t)
    lower_endpoint_blocks = p31_exceptional_lower_endpoint_certificate()
    compact = _P31_FIRST_INTERIOR_COMPACT
    all_equal = _P31_FIRST_INTERIOR_ALL_EQUAL
    compact_six = tuple(compact_moment(p, *compact, 6, k) for k in range(3))
    all_equal_six = tuple(
        tuple(all_equal_moment(p, *triangle, 6, k) for k in range(3))
        for triangle in all_equal
    )
    total_six = tuple(
        (compact_six[k] + sum(vector[k] for vector in all_equal_six)) % p
        for k in range(3)
    )
    odd_totals = {
        degree: tuple(
            (
                compact_moment(p, *compact, degree, k)
                + sum(
                    all_equal_moment(p, *triangle, degree, k)
                    for triangle in all_equal
                )
            )
            % p
            for k in range(degree // 2)
        )
        for degree in range(3, p - 1, 2)
    }
    degree_eight_total = tuple(
        (
            compact_moment(p, *compact, 8, k)
            + sum(all_equal_moment(p, *triangle, 8, k) for triangle in all_equal)
        )
        % p
        for k in range(4)
    )
    compact_counts = [
        int(item["Q"]) - int(row["r"]) - 2 for item in row["opposite_rows"]
    ]
    hard_counts = [int(item["e"]) for item in row["hard_rows"]]
    hard_four_odd_totals = {
        degree: tuple(
            sum(
                compact_moment(p, *triangle, degree, k)
                for triangle in _P31_COMPACT_FOUR
            )
            % p
            for k in range(degree // 2)
        )
        for degree in range(3, p - 1, 2)
    }
    hard_five_odd_totals = {
        degree: tuple(
            sum(
                compact_moment(p, scale, -scale, 0, degree, k)
                for scale in _P31_COMPACT_FIVE_SCALES
            )
            % p
            for k in range(degree // 2)
        )
        for degree in range(3, p - 1, 2)
    }
    opposite_zero_odd_totals = {
        degree: tuple(
            sum(
                all_equal_moment(p, scale, -scale, 0, degree, k)
                for scale in _P31_ALL_EQUAL_SIX_SCALES
            )
            % p
            for k in range(degree // 2)
        )
        for degree in range(3, p - 1, 2)
    }
    distinct_labels = all(len(set(triangle)) == 3 for triangle in all_equal)
    proved = bool(
        lower_endpoint_blocks["proved"]
        and compact_six == (29, 0, 4)
        and total_six == (0, 0, 0)
        and all(not any(vector) for vector in odd_totals.values())
        and all(not any(vector) for vector in hard_four_odd_totals.values())
        and all(not any(vector) for vector in hard_five_odd_totals.values())
        and all(not any(vector) for vector in opposite_zero_odd_totals.values())
        and sorted(hard_counts) == [4] * 10 + [5] * 6
        and sorted(compact_counts) == [0] * 15 + [1]
        and distinct_labels
    )
    if not proved:
        raise ArithmeticError("the p=31 first-interior odd/six block changed")
    return {
        "p": p,
        "r": 7,
        "t": t,
        "opposite_compact_count_multiset": {"0": 15, "1": 1},
        "hard_compact_count_multiset": {"4": 10, "5": 6},
        "compact_triangle": list(compact),
        "all_equal_triangles": [list(triangle) for triangle in all_equal],
        "compact_degree_six_vector": list(compact_six),
        "all_equal_degree_six_vectors": [list(vector) for vector in all_equal_six],
        "combined_degree_six_vector": list(total_six),
        "all_odd_degrees_below_top_combined_zero": True,
        "lower_endpoint_blocks_reused_on_hard_and_b0_rows": True,
        "whole_balanced_row_allocation_has_zero_odd_and_degree_six_forms": True,
        "degree_eight_combined_vector": list(degree_eight_total),
        "degree_eight_also_zero": not any(degree_eight_total),
        "same_block_passes_degree_eight": False,
        "atom_repetition_or_overlap_forbidden_by_prop15758": False,
        "proved": proved,
    }


def p3_lower_endpoint_degree_eight_standalone_certificate(p: int) -> dict[str, object]:
    """Standalone degree-eight compatibility at the branch-C lower endpoint.

    This deliberately does not claim the same scales satisfy degree six.
    For ``p=3 mod 4``, nonzero eighth powers equal nonzero squares, and
    Cauchy--Davenport gives ``3H_8=F_p``.  Hence each exact atom count can
    be assigned eighth-power sum zero, independently of its sixth-power sum.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    t = 2 * r * r - 4 * r - 2
    row = p3_local_survivor(p, t)
    hard_counts = tuple(int(item["e"]) for item in row["hard_rows"])
    opposite_count = r - 1
    bound = eighth_power_sumset_bound(p)
    threshold = int(bound["cauchy_davenport_threshold"])
    proved = bool(
        gcd(8, p - 1) == 2
        and threshold == 3
        and min(hard_counts) >= threshold
        and opposite_count >= threshold
    )
    if not proved:
        raise ArithmeticError("the standalone p=3 degree-eight proof failed")
    return {
        "p": p,
        "r": r,
        "t_lower_endpoint": t,
        "hard_compact_atom_counts": list(hard_counts),
        "opposite_all_equal_count_per_direction": opposite_count,
        "nonzero_eighth_powers_are_nonzero_squares": True,
        "cauchy_davenport_threshold": threshold,
        "global_degree_eight_forms": {"k_0": "0", "k_1": "0", "k_2": "0", "k_3": "0"},
        "standalone_degree_eight_only": True,
        "same_labels_as_degree_six_not_proved": True,
        "all_degree_eight_moment_relations_pass_separately": True,
        "proved": proved,
    }


_P31_EVEN_CHANNELS = tuple(
    [(6, k) for k in range(3)] + [(8, k) for k in range(4)]
)


@lru_cache(maxsize=1)
def _p31_ae_even_catalog() -> tuple[
    tuple[tuple[int, int, int], ...], tuple[tuple[int, ...], ...]
]:
    triples = tuple(combinations(range(31), 3))
    vectors = tuple(
        tuple(all_equal_moment(31, *triangle, degree, k) for degree, k in _P31_EVEN_CHANNELS)
        for triangle in triples
    )
    return triples, vectors


def _p31_vector_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((a + b) % 31 for a, b in zip(left, right))


def _p31_vector_sub(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((a - b) % 31 for a, b in zip(left, right))


def _p31_vector_scale(scalar: int, vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(scalar * value % 31 for value in vector)


def _p31_vector_encode(vector: tuple[int, ...]) -> int:
    return sum(value * 31**index for index, value in enumerate(vector))


def _p31_centered_compact_even_vector() -> tuple[int, ...]:
    return tuple(
        compact_moment(31, *_P31_FIRST_INTERIOR_COMPACT, degree, k)
        for degree, k in _P31_EVEN_CHANNELS
    )


def p31_centered_joint_six_eight_gate_certificate() -> dict[str, object]:
    """Exact one-field certificate that four centered atoms cannot kill d6,d8.

    A centered atom scale contributes the pair ``(a^6,a^8)``.  There are
    only 15 distinct pairs modulo 31 (the signs agree).  Exact four-fold
    sumset propagation excludes ``(0,0)``.  This is a targeted joint-moment
    calculation at the first live endpoint, not a prime census, and says
    nothing about noncentered compact triangles.
    """
    p = 31
    pair_to_scale: dict[tuple[int, int], int] = {}
    for scale in range(1, p):
        pair_to_scale.setdefault((pow(scale, 6, p), pow(scale, 8, p)), scale)
    reachable: dict[tuple[int, int], tuple[int, ...]] = {(0, 0): ()}
    cardinalities: dict[int, int] = {}
    zero_witnesses: dict[int, list[int] | None] = {}
    for count in range(1, 7):
        next_reachable: dict[tuple[int, int], tuple[int, ...]] = {}
        for state, witness in reachable.items():
            for pair, scale in pair_to_scale.items():
                target = ((state[0] + pair[0]) % p, (state[1] + pair[1]) % p)
                next_reachable.setdefault(target, witness + (scale,))
        reachable = next_reachable
        cardinalities[count] = len(reachable)
        zero = reachable.get((0, 0))
        zero_witnesses[count] = list(zero) if zero is not None else None
    proved = bool(
        len(pair_to_scale) == 15
        and zero_witnesses[4] is None
        and zero_witnesses[5] is not None
        and zero_witnesses[6] is not None
    )
    if not proved:
        raise ArithmeticError("the p=31 centered joint sumset changed")
    return {
        "p": p,
        "distinct_nonzero_scale_moment_pairs": len(pair_to_scale),
        "reachable_pair_counts": cardinalities,
        "zero_pair_witnesses": zero_witnesses,
        "four_centered_compact_atoms_can_sum_to_zero_in_degrees_six_and_eight": False,
        "noncentered_compact_atoms_ruled_out": False,
        "proved": proved,
    }


def p31_centered_compact_three_ae_pairs_joint_no_go() -> dict[str, object]:
    r"""Exclude every blockwise-negation-fixed six-AE extension.

    Fix the centered compact atom ``(1,-1,0)`` in the unique ``b=1`` row at
    ``p=31,t=69``.  Three nonfixed negation pairs would have to have
    representative sum

    ``(1,0,29 | 1,0,0,2)``

    across degrees six and eight.  This function exhausts the 4,495
    unordered distinct-label representatives and their 10,104,760 pair
    sums in this one field.  It also checks every other six-block multiset
    fixed blockwise by negation: two pairs plus two invariant triangles, one
    pair plus four invariant triangles, and six invariant triangles.  The
    later Radon/trade certificate is needed to pass from odd-zero aggregate
    edge symmetry to these blockwise configurations and the Pasch/volume-six
    alternatives.  No nonzero global forms are excluded here.
    """
    import numpy as np

    p = 31
    triples, vector_rows = _p31_ae_even_catalog()
    vectors = np.asarray(vector_rows, dtype=np.int16)
    compact_vector = np.asarray(_p31_centered_compact_even_vector(), dtype=np.int16)
    target = (-compact_vector * pow(2, -1, p)) % p
    expected_target = np.asarray((1, 0, 29, 1, 0, 0, 2), dtype=np.int16)
    powers = np.asarray([p**index for index in range(7)], dtype=np.int64)
    pair_count = len(triples) * (len(triples) + 1) // 2
    encoded_pair_sums = np.empty(pair_count, dtype=np.int64)
    offset = 0
    for index, vector in enumerate(vectors):
        row = ((vectors[index:] + vector) % p).astype(np.int64) @ powers
        encoded_pair_sums[offset : offset + len(row)] = row
        offset += len(row)
    distinct_pair_sums = np.unique(encoded_pair_sums)
    target_minus_third = ((target - vectors) % p).astype(np.int64) @ powers
    solution_exists = bool(np.isin(target_minus_third, distinct_pair_sums).any())
    full_target = (-compact_vector) % p
    raw_vector_codes = np.unique(vectors.astype(np.int64) @ powers)

    invariant_vectors = tuple(
        tuple(
            all_equal_moment(p, 0, scale, -scale, degree, k)
            for degree in (6, 8)
            for k in range(degree // 2)
        )
        for scale in range(1, 16)
    )

    def invariant_sumset(term_count: int) -> set[tuple[int, ...]]:
        states: set[tuple[int, ...]] = {(0,) * 7}
        for _ in range(term_count):
            states = {
                tuple((state[index] + vector[index]) % p for index in range(7))
                for state in states
                for vector in invariant_vectors
            }
        return states

    invariant_two = invariant_sumset(2)
    invariant_four = invariant_sumset(4)
    invariant_six = invariant_sumset(6)

    def encode(vector: tuple[int, ...]) -> int:
        return sum(value * int(powers[index]) for index, value in enumerate(vector))

    def sorted_contains(values: np.ndarray, value: int) -> bool:
        location = int(np.searchsorted(values, value))
        return location < len(values) and int(values[location]) == value

    inverse_two = pow(2, -1, p)
    two_invariant_two_pair_solution = any(
        sorted_contains(
            distinct_pair_sums,
            encode(
                tuple(
                    inverse_two * (int(full_target[index]) - invariant[index]) % p
                    for index in range(7)
                )
            ),
        )
        for invariant in invariant_two
    )
    four_invariant_one_pair_solution = any(
        sorted_contains(
            raw_vector_codes,
            encode(
                tuple(
                    inverse_two * (int(full_target[index]) - invariant[index]) % p
                    for index in range(7)
                )
            ),
        )
        for invariant in invariant_four
    )
    six_invariant_solution = tuple(map(int, full_target)) in invariant_six
    category_solutions = {
        "0_invariant_3_pairs": solution_exists,
        "2_invariant_2_pairs": two_invariant_two_pair_solution,
        "4_invariant_1_pair": four_invariant_one_pair_solution,
        "6_invariant": six_invariant_solution,
    }

    import hashlib
    import json

    audit_payload = {
        "p": p,
        "target": tuple(map(int, full_target)),
        "triple_count": len(triples),
        "distinct_raw_vectors": len(raw_vector_codes),
        "distinct_unordered_pair_sums": len(distinct_pair_sums),
        "invariant_triangle_count": len(invariant_vectors),
        "distinct_invariant_2_sums": len(invariant_two),
        "distinct_invariant_4_sums": len(invariant_four),
        "distinct_invariant_6_sums": len(invariant_six),
        "solutions": {
            name: True for name, exists in category_solutions.items() if exists
        },
    }
    payload = json.dumps(audit_payload, sort_keys=True, separators=(",", ":"))
    evidence_sha256 = hashlib.sha256(payload.encode()).hexdigest()
    proved = bool(
        len(triples) == 4495
        and pair_count == 10_104_760
        and len(distinct_pair_sums) == 2_543_460
        and len(raw_vector_codes) == 2_255
        and (len(invariant_two), len(invariant_four), len(invariant_six))
        == (120, 925, 961)
        and np.array_equal(target, expected_target)
        and not any(category_solutions.values())
        and evidence_sha256
        == "26bea31c9906b005ff4fc1dc0121d43eb07ef7f62369b90b902026ae0d293c95"
    )
    if not proved:
        raise ArithmeticError("the p=31 structured joint no-go changed")
    return {
        "p": p,
        "t": 69,
        "fixed_centered_compact": list(_P31_FIRST_INTERIOR_COMPACT),
        "required_three_representative_sum": target.tolist(),
        "distinct_label_all_equal_triples": len(triples),
        "unordered_pair_sums_with_repetition": pair_count,
        "distinct_pair_sum_vectors": len(distinct_pair_sums),
        "three_all_equal_representatives_exist": solution_exists,
        "centered_compact_plus_three_AE_negation_pairs_can_kill_d6_and_d8": False,
        "full_joint_target": full_target.tolist(),
        "distinct_raw_AE_vectors": len(raw_vector_codes),
        "invariant_triangle_count": len(invariant_vectors),
        "distinct_invariant_sum_vectors": {
            "2": len(invariant_two),
            "4": len(invariant_four),
            "6": len(invariant_six),
        },
        "blockwise_negation_fixed_category_solutions": category_solutions,
        "all_blockwise_negation_fixed_six_AE_configurations_excluded": True,
        "evidence_sha256": evidence_sha256,
        "arbitrary_odd_zero_atom_configurations_ruled_out": False,
        "proved": proved,
    }


def p3_full_balanced_maximal_line_exclusion_certificate(
    p: int, compact_count: int
) -> dict[str, object]:
    r"""Exclude every maximal-line odd/Radon word for ``0<=b<=r``.

    Let ``p=4r+3``, ``h=2r+1``, and let the signed edge chain contain
    ``b`` compact atoms and ``r-1`` positive all-equal triangle boundaries.
    Its integer orbit differences have total ``l1`` mass at most

    ``N=3(r+b-1)<=3h-6``

    and individual absolute value at most ``B=r-1+2b<p``.  Conditional on
    the odd/Radon word being supported on one of the maximal lines in
    ``H x H``, these bounds exclude all three line types for every
    ``0<=b<=r``.

    On a horizontal or diagonal maximal line the unique dual relation makes
    the ``h`` nonzero integers ``n_E`` represent all projective residue
    classes in ``F_p^*/{+1,-1}``.  Their ``l1`` mass is therefore at least
    ``1+...+h=h(h+1)/2>3h-6``.

    On a vertical line all ``n_E`` have one nonzero residue ``kappa`` modulo
    ``p``.  If ``a=min(|kappa|,p-|kappa|)``, then ``ha<=N<3h``, so
    ``a<=2``.  A coordinate using the other integer lift has total ``l1``
    at least ``(h-1)a+(p-a)``; already at ``a=1`` this is ``3h-1>N``.
    Hence all coordinates equal the same integer of absolute value one or
    two.

    For the unit case, reduce the chain modulo two and project vertices by
    ``x~-x``.  Signs disappear, every all-equal or compact atom becomes a
    projected triangle boundary, and the projected orbit-difference graph
    is Eulerian.  The ``h`` fixed-sum edges have degree two at every quotient
    vertex except ``[0]`` and ``[sigma/2]``, where they have degree one.
    Thus a constant odd coefficient is impossible.  A constant coefficient
    of absolute value two needs at least ``2h`` occurrences aligned with its
    sign, while the atomwise aligned capacity is only
    ``(r-1)+2b<=3r-1<2h``.

    This is only a maximal-line support exclusion.  It does not extend the
    line-isolation theorem past support ``2h-3``, exclude conic or cubic
    supports, or prove centrality outside the previously certified band.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    if (
        not isinstance(compact_count, int)
        or isinstance(compact_count, bool)
        or compact_count < 0
    ):
        raise ValueError("compact_count must be a nonnegative integer")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    b = compact_count
    if b > r:
        raise ValueError("the full balanced branch has compact_count<=r")

    h = 2 * r + 1
    signed_occurrences = 3 * (r + b - 1)
    full_balanced_occurrence_bound = 3 * h - 6
    orbit_difference_bound = r - 1 + 2 * b

    horizontal_diagonal_l1_floor = h * (h + 1) // 2
    horizontal_diagonal_l1_margin = (
        horizontal_diagonal_l1_floor - signed_occurrences
    )

    maximum_vertical_canonical_absolute_value = signed_occurrences // h
    unit_alternative_lift_l1_floor = (h - 1) + (p - 1)
    two_alternative_lift_l1_floor = 2 * (h - 1) + (p - 2)
    alternative_lift_l1_margin = (
        unit_alternative_lift_l1_floor - signed_occurrences
    )

    quotient_vertex_count = h + 1
    fixed_sum_projected_edge_count = h
    fixed_sum_projected_degree_one_vertices = ("[0]", "[sigma/2]")
    fixed_sum_projected_degree_two_vertex_count = h - 1
    double_vertical_aligned_demand = 2 * h
    atomwise_aligned_capacity = r - 1 + 2 * b
    double_vertical_aligned_deficit = (
        double_vertical_aligned_demand - atomwise_aligned_capacity
    )

    proved = bool(
        b <= r
        and h == (p - 1) // 2
        and signed_occurrences <= full_balanced_occurrence_bound
        and full_balanced_occurrence_bound < 3 * h
        and orbit_difference_bound < p
        and horizontal_diagonal_l1_floor > full_balanced_occurrence_bound
        and horizontal_diagonal_l1_margin > 0
        and maximum_vertical_canonical_absolute_value <= 2
        and unit_alternative_lift_l1_floor == 3 * h - 1
        and two_alternative_lift_l1_floor == 4 * h - 3
        and alternative_lift_l1_margin > 0
        and quotient_vertex_count == h + 1
        and fixed_sum_projected_edge_count == h
        and len(fixed_sum_projected_degree_one_vertices) == 2
        and fixed_sum_projected_degree_two_vertex_count == h - 1
        and 2 + 2 * fixed_sum_projected_degree_two_vertex_count
        == 2 * fixed_sum_projected_edge_count
        and atomwise_aligned_capacity == orbit_difference_bound
        and double_vertical_aligned_deficit > 0
    )
    if not proved:
        raise ArithmeticError("the full balanced maximal-line exclusion changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": b,
        "full_balanced_compact_count_hypothesis": "0<=b<=r",
        "square_parameter_count": h,
        "signed_edge_occurrence_bound": signed_occurrences,
        "full_balanced_signed_occurrence_bound": full_balanced_occurrence_bound,
        "total_orbit_difference_bound": orbit_difference_bound,
        "horizontal_diagonal_projective_class_count": h,
        "horizontal_diagonal_l1_floor": horizontal_diagonal_l1_floor,
        "horizontal_diagonal_l1_margin": horizontal_diagonal_l1_margin,
        "maximum_vertical_canonical_absolute_value": (
            maximum_vertical_canonical_absolute_value
        ),
        "unit_alternative_lift_l1_floor": unit_alternative_lift_l1_floor,
        "two_alternative_lift_l1_floor": two_alternative_lift_l1_floor,
        "alternative_lift_l1_margin": alternative_lift_l1_margin,
        "projected_vertex_set": "F_p/{x~-x}",
        "projected_atom_chains_are_Eulerian_mod_two": True,
        "projected_fixed_sum_edge_count": fixed_sum_projected_edge_count,
        "projected_fixed_sum_degree_one_vertices": list(
            fixed_sum_projected_degree_one_vertices
        ),
        "projected_fixed_sum_degree_two_vertex_count": (
            fixed_sum_projected_degree_two_vertex_count
        ),
        "unit_vertical_line_excluded_by_projected_parity": True,
        "double_vertical_aligned_occurrence_demand": (
            double_vertical_aligned_demand
        ),
        "atomwise_aligned_occurrence_capacity": atomwise_aligned_capacity,
        "double_vertical_aligned_occurrence_deficit": (
            double_vertical_aligned_deficit
        ),
        "double_vertical_line_excluded_by_aligned_incidence": True,
        "horizontal_diagonal_maximal_lines_excluded": True,
        "vertical_maximal_lines_excluded": True,
        "all_maximal_line_supports_excluded": True,
        "support_isolation_extended_past_2h_minus_3": False,
        "conic_supports_excluded": False,
        "cubic_supports_excluded": False,
        "aggregate_signed_edge_chain_is_centrally_symmetric": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_low_weight_line_peeling_certificate(p: int) -> dict[str, object]:
    r"""Reduce every low-weight line-containing support to two lines.

    Put ``h=(p-1)/2`` and ``m=h-2``.  Let a nonzero word orthogonal to all
    polynomials of total degree at most ``m`` on ``H x H`` have support
    ``S`` of size at most ``3m=3h-6``.  If ``S`` contains ``h`` collinear
    points on a line ``L``, then either ``S`` is that line or

    ``S subset L union L2`` and ``h-1 <= |S minus L| <= h``,

    where ``L2`` is another maximal line of ``H x H``.

    Indeed, if ``R=S minus L`` is nonempty, multiplication of every test
    polynomial by an equation of ``L`` makes the nonzero coefficients on
    ``R`` a degree-``m-1`` dual word.  Since ``|R|<=2h-6``, Couvreur's
    first two linked-configuration thresholds force ``h-1`` points of
    ``R`` onto a line ``L2``: without them a linked set has size at least
    ``2(m-1)+2=2h-4``.  If a point remained off ``L2``, multiplying once
    more would make the remainder degree-``m-2`` linked, of size at least
    ``m=h-2``.  This would give ``|S|>=3h-3``, a contradiction.  Finally,
    a nonmaximal line meets ``H x H`` in at most ``r+1<h-1`` points.

    The linked-set thresholds are Theorem 3.8 of Alain Couvreur, *The dual
    minimum distance of arbitrary-dimensional algebraic-geometric codes*,
    J. Algebra 350 (2012), 84--107.  Lemma 2.13 there makes them invariant
    under algebraic field extension, so they apply to these F_p-rational
    affine points after homogenization.

    This support theorem does not by itself exclude the resulting two-line
    coefficients.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")

    h = 2 * r + 1
    m = h - 2
    support_bound = 3 * m
    first_residual_degree = m - 1
    first_residual_support_bound = support_bound - h
    forced_second_line_points = first_residual_degree + 2
    no_second_line_linked_threshold = 2 * first_residual_degree + 2
    second_residual_degree = first_residual_degree - 1
    off_second_line_support_floor = second_residual_degree + 2
    three_piece_support_floor = (
        h + forced_second_line_points + off_second_line_support_floor
    )
    nonmaximal_line_intersection_bound = r + 1

    proved = bool(
        m == h - 2
        and support_bound == 3 * h - 6
        and first_residual_degree == h - 3
        and first_residual_support_bound == 2 * h - 6
        and forced_second_line_points == h - 1
        and no_second_line_linked_threshold == 2 * h - 4
        and first_residual_support_bound
        < no_second_line_linked_threshold
        and second_residual_degree == h - 4
        and off_second_line_support_floor == h - 2
        and three_piece_support_floor == 3 * h - 3
        and three_piece_support_floor > support_bound
        and nonmaximal_line_intersection_bound < forced_second_line_points
    )
    if not proved:
        raise ArithmeticError("the low-weight line-peeling theorem changed")
    return {
        "p": p,
        "r": r,
        "square_parameter_count": h,
        "dual_polynomial_total_degree": m,
        "support_hypothesis_bound": support_bound,
        "first_residual_dual_degree": first_residual_degree,
        "first_residual_support_bound": first_residual_support_bound,
        "forced_second_line_point_count": forced_second_line_points,
        "no_second_line_linked_support_threshold": (
            no_second_line_linked_threshold
        ),
        "second_residual_dual_degree": second_residual_degree,
        "off_second_line_support_floor": off_second_line_support_floor,
        "forbidden_three_piece_support_floor": three_piece_support_floor,
        "three_piece_support_margin": three_piece_support_floor - support_bound,
        "nonmaximal_line_intersection_bound": (
            nonmaximal_line_intersection_bound
        ),
        "support_if_not_one_line": (
            "L union L2 with L maximal, L2 maximal, and "
            "h-1<=|S\\L|<=h"
        ),
        "couvreur_reference": (
            "A. Couvreur, J. Algebra 350 (2012), Theorem 3.8 and Lemma 2.13"
        ),
        "finite_field_extension_hypothesis_checked": True,
        "line_containing_support_reduced_to_two_maximal_lines": True,
        "two_maximal_line_coefficients_excluded": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_full_balanced_two_maximal_line_exclusion_certificate(
    p: int, compact_count: int
) -> dict[str, object]:
    r"""Exclude the two-maximal-line supports produced by line peeling.

    This combines the affine-Cartesian dual identity

    ``C_{H x H}(h-2)^perp = diag(U D) C_{H x H}(h-1)``

    (up to the harmless scalar ``h^-2``) with exact ``l1`` and quotient
    parity bounds for a chain of ``b`` compact and ``r-1`` all-equal atoms.
    It applies throughout ``0<=b<=r``.

    Two maximal lines from different families meet in ``H x H``.  Their
    union has ``2h-1`` points, and the degree-``h-2`` dual space on the union
    is two-dimensional, spanned by the two individual line relations.  At
    least one component is nonvertical; on its ``h-1`` exclusive points its
    orbit differences occupy distinct projective classes.  This costs at
    least ``1+...+(h-1)>3h-6`` in integer ``l1`` mass.

    For two lines in the same family, a degree-``h-1`` Cartesian polynomial
    supported on them factors as the product of the other ``h-2`` line
    equations times one affine-linear factor.  On a full vertical component,
    a nonconstant factor makes ``n_E`` affine and injective in ``D``; the
    least ``l1`` mass of ``h=2r+1`` distinct nonzero residues is
    ``(r+1)^2>6r-3``.  A constant factor makes both components full and
    constant.  The total mass forces both constants to be actual units; the
    two projected fixed-sum graphs leave the distinct odd vertices
    ``[sigma_1/2]`` and ``[sigma_2/2]``.

    On a full horizontal component, projective equality of coefficients is
    governed by the cubic ``U(A+B U)^2``.  On a full diagonal component it
    is governed by ``U^-3(A+B U)^2``, again a cubic after clearing
    denominators.  Hence one projective class occurs at most three times.
    Filling the smallest classes three times gives an exact ``l1`` floor
    still larger than ``3h-6`` for ``h>=15``.
    """
    peeling = p3_low_weight_line_peeling_certificate(p)
    line_only = p3_full_balanced_maximal_line_exclusion_certificate(
        p, compact_count
    )
    r = int(peeling["r"])
    h = int(peeling["square_parameter_count"])
    b = compact_count

    signed_occurrences = 3 * (r + b - 1)
    full_balanced_occurrence_bound = 3 * h - 6

    different_family_support_size = 2 * h - 1
    different_family_evaluation_rank = 2 * h - 3
    different_family_dual_nullity = (
        different_family_support_size - different_family_evaluation_rank
    )
    different_family_l1_floor = h * (h - 1) // 2
    different_family_l1_margin = (
        different_family_l1_floor - signed_occurrences
    )

    same_vertical_distinct_residue_l1_floor = (r + 1) ** 2
    same_vertical_distinct_residue_l1_margin = (
        same_vertical_distinct_residue_l1_floor - signed_occurrences
    )
    two_vertical_canonical_absolute_sum_bound = signed_occurrences // h
    two_vertical_alternative_lift_l1_floor = 4 * h - 1
    same_vertical_projected_odd_vertices = (
        "[sigma_1/2]",
        "[sigma_2/2]",
    )

    projective_fibre_bound = 3
    quotient, remainder = divmod(h, projective_fibre_bound)
    three_to_one_projective_l1_floor = (
        projective_fibre_bound * quotient * (quotient + 1) // 2
        + remainder * (quotient + 1)
    )
    three_to_one_projective_l1_margin = (
        three_to_one_projective_l1_floor - signed_occurrences
    )

    proved = bool(
        peeling["proved"]
        and line_only["proved"]
        and b <= r
        and signed_occurrences <= full_balanced_occurrence_bound
        and different_family_dual_nullity == 2
        and different_family_l1_floor > full_balanced_occurrence_bound
        and different_family_l1_margin > 0
        and same_vertical_distinct_residue_l1_floor
        - full_balanced_occurrence_bound
        == (r - 2) ** 2
        and same_vertical_distinct_residue_l1_margin > 0
        and two_vertical_canonical_absolute_sum_bound <= 2
        and two_vertical_alternative_lift_l1_floor
        > full_balanced_occurrence_bound
        and len(same_vertical_projected_odd_vertices) == 2
        and three_to_one_projective_l1_floor
        > full_balanced_occurrence_bound
        and three_to_one_projective_l1_margin > 0
    )
    if not proved:
        raise ArithmeticError("the two-maximal-line exclusion changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": b,
        "full_balanced_compact_count_hypothesis": "0<=b<=r",
        "square_parameter_count": h,
        "signed_edge_occurrence_bound": signed_occurrences,
        "full_balanced_signed_occurrence_bound": full_balanced_occurrence_bound,
        "affine_Cartesian_dual_multiplier": "U*D/h^2",
        "same_family_factor_degree": h - 2,
        "same_family_residual_factor_degree": 1,
        "different_family_support_size": different_family_support_size,
        "different_family_evaluation_rank": different_family_evaluation_rank,
        "different_family_dual_nullity": different_family_dual_nullity,
        "different_family_dual_basis": "the two individual line relations",
        "different_family_l1_floor": different_family_l1_floor,
        "different_family_l1_margin": different_family_l1_margin,
        "same_vertical_nonconstant_coefficients_are_injective": True,
        "same_vertical_distinct_residue_l1_floor": (
            same_vertical_distinct_residue_l1_floor
        ),
        "same_vertical_distinct_residue_l1_margin": (
            same_vertical_distinct_residue_l1_margin
        ),
        "two_vertical_canonical_absolute_sum_bound": (
            two_vertical_canonical_absolute_sum_bound
        ),
        "two_vertical_alternative_lift_l1_floor": (
            two_vertical_alternative_lift_l1_floor
        ),
        "same_vertical_projected_degree_one_vertices": list(
            same_vertical_projected_odd_vertices
        ),
        "same_vertical_constant_case_excluded_by_projected_parity": True,
        "horizontal_projective_class_equation": "z=U*(A+B*U)^2",
        "diagonal_projective_class_equation": "z*U^3=(A+B*U)^2",
        "horizontal_diagonal_projective_fibre_bound": projective_fibre_bound,
        "three_to_one_projective_l1_floor": (
            three_to_one_projective_l1_floor
        ),
        "three_to_one_projective_l1_margin": (
            three_to_one_projective_l1_margin
        ),
        "all_two_maximal_line_supports_excluded": True,
        "all_supports_containing_h_collinear_points_excluded": True,
        "conic_without_h_collinear_points_excluded": False,
        "cubic_supports_excluded": False,
        "aggregate_signed_edge_chain_is_centrally_symmetric": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_boundary_cubic_unit_reduction_certificate(p: int) -> dict[str, object]:
    r"""Exclude every boundary cubic for ``p=4r+3>=31``.

    Let ``p=4r+3``, ``h=2r+1``, and ``m=h-2``.  In Couvreur's third
    linked-configuration alternative the support has size
    ``3m=3h-6`` and is the complete intersection of a cubic ``F`` and a
    degree-``m`` curve ``G``.  The balanced atom budget reaches this size
    only at ``b=r``.  Equality of support and occurrence counts forces every
    integer orbit difference to be ``+1`` or ``-1``.

    The complete intersection is reduced and transverse.  Its unique dual
    relation has Cayley--Bacharach residue weights

    ``W(P)=lambda/J(F,G)(P)``.

    Suppose first that the cubic is reducible and write ``F=L*Q``.  Its
    line component contains exactly ``m=h-2`` support points.  Since a
    nonmaximal line meets ``H x H`` in at most ``r+1<m``, ``L`` is maximal
    and its support omits two parameters ``a,b`` from ``H``.  On ``L``,
    ``G`` restricts, up to a scalar, to

    ``(z^h-1)/((z-a)(z-b))``.

    The reciprocal derivative at a support point is proportional to
    ``z(z-a)(z-b)``.  Combining this with ``W=n*D*sigma`` and ``n^2=1``
    gives a polynomial identity on ``h-2>5`` values.  A horizontal line
    would require ``Q_L(z)^2`` to be proportional to
    ``z(z-a)^2(z-b)^2``; a diagonal line would require
    ``(z-a)^2(z-b)^2`` to be proportional to ``z Q_L(z)^2``.  Both are
    impossible by the odd valuation at ``z=0``.  Thus ``L`` is vertical,
    and then the identity says

    ``Q_L(z) proportional to (z-a)(z-b)``.

    Hence the conic meets the vertical line at exactly its two omitted grid
    points.  If the conic split into two more lines, applying the same
    argument to all three components would make them vertical; the other
    two line factors are constant on the first component, contradicting the
    displayed quadratic restriction.  A reducible candidate is therefore a
    vertical line plus an absolutely irreducible conic containing at least
    ``2m+2=2h-2=p-3`` grid points.

    The high-intersection conic theorem in
    ``NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md`` puts that conic in the
    tangent normal form

    ``U=u*z^2, D=d*(z-1)^2``.

    Write the vertical line as ``U=u*c^2``.  Its two intersections with the
    conic are ``z=+/-c``.  The conic component has ``2m=p-5`` support
    points, so among the ``p-2`` grid parameters it omits ``+/-c`` and one
    further parameter ``e``.  With

    ``P_T(x)=product_(t != 0,1) (x-t)``,

    the restriction of the degree-``m`` curve is

    ``G|_Q = C*P_T(z)/((z^2-c^2)*(z-e))``.

    Since ``P_T'(z)=-1/(z*(z-1))`` at its roots and the Hamiltonian tangent
    to the displayed conic is a constant multiple of its parametrized
    tangent, the factor ``z^2-c^2`` from the line cancels in
    ``J(F,G)``.  Thus on the conic support

    ``W(z)=C*z*(z-1)*(z-e)`` and
    ``n(z)=W/(D*sigma)=C'*(z-e)/(z-1)``.

    Boundary saturation gives ``n(z)^2=1`` at ``p-5>2`` points.  Therefore
    the degree-two polynomial

    ``(z-e)^2-C''*(z-1)^2``

    vanishes identically, which forces ``e=1``.  This contradicts the fact
    that ``e`` is a grid parameter.  Hence every reducible boundary cubic
    is excluded; this does not exclude a general conic-supported word.

    A singular geometrically integral cubic has at most ``p+2`` rational
    points and is already too small.  For a smooth cubic ``C``, the
    coordinate ``U=X/Z`` is a nonconstant rational function of degree two
    or three.  It is geometrically nonsquare: a square of degree at most
    three would have a degree-one square root, impossible on a genus-one
    curve.  The connected double cover ``Y^2=U`` has at most six branch
    points and genus at most four.  Every one of the ``3h-6`` support points
    has two rational lifts, whereas Weil gives at most
    ``p+1+8*sqrt(p)`` points.  The contradiction

    ``3p-15 > p+1+8*sqrt(p)``

    is equivalent to ``p-8>4*sqrt(p)`` and already holds at ``p=31``.
    Thus the smooth case is excluded for every prime in scope as well.

    This theorem excludes the boundary cubic alternative, not the separate
    high-intersection conic alternative.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")

    h = 2 * r + 1
    m = h - 2
    boundary_support = 3 * m
    boundary_compact_count = r
    boundary_occurrences = 3 * (r + boundary_compact_count - 1)
    nonmaximal_line_intersection_bound = r + 1
    line_component_support = m
    line_component_omitted_grid_points = h - line_component_support
    conic_component_support = 2 * m
    reciprocal_derivative_degree = 3
    squared_identity_degree_bound = 5

    singular_integral_cubic_point_bound = p + 2
    singular_integral_cubic_margin = (
        boundary_support - singular_integral_cubic_point_bound
    )
    hasse_squared_margin = (p - 17) ** 2 - 16 * p
    smooth_irreducible_excluded_by_hasse = bool(
        p >= 47 and p > 17 and hasse_squared_margin > 0
    )
    direct_cubic_hasse_forces_reducible = smooth_irreducible_excluded_by_hasse

    certified_conic_grid_point_lower_bound = (
        conic_component_support + line_component_omitted_grid_points
    )
    couvreur_conic_threshold = 2 * m + 2
    tangent_conic_grid_point_count = p - 2
    conic_omitted_grid_parameters = (
        tangent_conic_grid_point_count - conic_component_support
    )
    line_conic_intersection_parameters = 2
    additional_omitted_conic_parameters = (
        conic_omitted_grid_parameters - line_conic_intersection_parameters
    )
    conic_restriction_degree = 2 * m
    conic_unit_identity_degree_bound = 2
    reducible_boundary_cubic_excluded = True
    coordinate_function_degree_lower_bound = 2
    coordinate_function_degree_upper_bound = 3
    double_cover_branch_point_bound = 6
    double_cover_genus_bound = 4
    double_cover_lifted_support_point_lower_bound = 2 * boundary_support
    double_cover_weil_squared_margin = (p - 8) ** 2 - 16 * p
    smooth_irreducible_excluded_by_double_cover = bool(
        p >= 31 and p > 8 and double_cover_weil_squared_margin > 0
    )
    all_boundary_cubic_supports_excluded = (
        reducible_boundary_cubic_excluded
        and smooth_irreducible_excluded_by_double_cover
    )

    proved = bool(
        m == h - 2
        and boundary_support == 3 * h - 6
        and boundary_occurrences == boundary_support
        and line_component_support == h - 2
        and line_component_omitted_grid_points == 2
        and conic_component_support == 2 * h - 4
        and nonmaximal_line_intersection_bound < line_component_support
        and reciprocal_derivative_degree == 3
        and line_component_support > squared_identity_degree_bound
        and singular_integral_cubic_margin > 0
        and smooth_irreducible_excluded_by_hasse == (p >= 47)
        and certified_conic_grid_point_lower_bound == 2 * h - 2
        and certified_conic_grid_point_lower_bound == couvreur_conic_threshold
        and tangent_conic_grid_point_count == p - 2
        and conic_component_support == p - 5
        and conic_omitted_grid_parameters == 3
        and line_conic_intersection_parameters == 2
        and additional_omitted_conic_parameters == 1
        and conic_restriction_degree == p - 5
        and conic_component_support > conic_unit_identity_degree_bound
        and reducible_boundary_cubic_excluded
        and coordinate_function_degree_lower_bound == 2
        and coordinate_function_degree_upper_bound == 3
        and double_cover_branch_point_bound == 6
        and double_cover_genus_bound == 4
        and double_cover_lifted_support_point_lower_bound == 3 * p - 15
        and double_cover_weil_squared_margin > 0
        and smooth_irreducible_excluded_by_double_cover
        and all_boundary_cubic_supports_excluded
    )
    if not proved:
        raise ArithmeticError("the boundary cubic unit exclusion changed")
    return {
        "p": p,
        "r": r,
        "square_parameter_count": h,
        "dual_polynomial_total_degree": m,
        "boundary_compact_atom_count": boundary_compact_count,
        "boundary_support_size": boundary_support,
        "boundary_signed_edge_occurrences": boundary_occurrences,
        "support_saturates_occurrence_budget": True,
        "all_nonzero_integer_orbit_differences_are_units": True,
        "cayley_bacharach_weight_formula": "W(P)=lambda/J(F,G)(P)",
        "complete_intersection_is_reduced_and_transverse": True,
        "line_component_support_size": line_component_support,
        "line_component_omitted_grid_points": (
            line_component_omitted_grid_points
        ),
        "nonmaximal_line_intersection_bound": (
            nonmaximal_line_intersection_bound
        ),
        "line_component_is_maximal": True,
        "line_restriction_of_degree_m_curve": (
            "(z^h-1)/((z-a)*(z-b))"
        ),
        "reciprocal_derivative_on_line": "z*(z-a)*(z-b)",
        "squared_line_identity_degree_bound": squared_identity_degree_bound,
        "horizontal_line_component_excluded_by_odd_zero_valuation": True,
        "diagonal_line_component_excluded_by_odd_zero_valuation": True,
        "surviving_line_component_type": "vertical U=u0",
        "surviving_conic_line_restriction": "Q(u0,z)=c*(z-a)*(z-b)",
        "three_line_cubic_excluded": True,
        "last_reducible_candidate": (
            "vertical line plus absolutely irreducible conic"
        ),
        "conic_component_support_size": conic_component_support,
        "conic_additional_omitted_line_points": (
            line_component_omitted_grid_points
        ),
        "surviving_conic_grid_point_count_lower_bound": (
            certified_conic_grid_point_lower_bound
        ),
        "couvreur_conic_configuration_threshold": couvreur_conic_threshold,
        "high_intersection_conic_normal_form": (
            "U=u*z^2, D=d*(z-1)^2"
        ),
        "tangent_conic_grid_point_count": tangent_conic_grid_point_count,
        "conic_omitted_grid_parameters": conic_omitted_grid_parameters,
        "line_conic_intersection_parameters": (
            line_conic_intersection_parameters
        ),
        "additional_omitted_conic_parameters": (
            additional_omitted_conic_parameters
        ),
        "degree_m_curve_restriction_on_conic": (
            "P_T(z)/((z^2-c^2)*(z-e))"
        ),
        "conic_restriction_degree": conic_restriction_degree,
        "reciprocal_jacobian_weight_on_conic": "z*(z-1)*(z-e)",
        "orbit_difference_on_conic": "C*(z-e)/(z-1)",
        "conic_unit_identity_degree_bound": conic_unit_identity_degree_bound,
        "unit_identity_forces_forbidden_e_equals_one": True,
        "reducible_boundary_cubic_excluded": (
            reducible_boundary_cubic_excluded
        ),
        "smooth_cubic_coordinate_function": "U=X/Z",
        "smooth_cubic_coordinate_function_degree_range": [
            coordinate_function_degree_lower_bound,
            coordinate_function_degree_upper_bound,
        ],
        "smooth_cubic_coordinate_function_geometrically_nonsquare": True,
        "double_cover_branch_point_bound": double_cover_branch_point_bound,
        "double_cover_genus_bound": double_cover_genus_bound,
        "double_cover_lifted_support_point_lower_bound": (
            double_cover_lifted_support_point_lower_bound
        ),
        "double_cover_weil_point_bound": "p+1+8*sqrt(p)",
        "double_cover_weil_squared_margin": double_cover_weil_squared_margin,
        "smooth_irreducible_cubic_excluded_by_double_cover": (
            smooth_irreducible_excluded_by_double_cover
        ),
        "singular_integral_cubic_point_bound": (
            singular_integral_cubic_point_bound
        ),
        "singular_integral_cubic_point_margin": singular_integral_cubic_margin,
        "hasse_squared_margin": hasse_squared_margin,
        "smooth_irreducible_cubic_excluded_by_hasse": (
            smooth_irreducible_excluded_by_hasse
        ),
        "direct_cubic_hasse_forces_reducible": (
            direct_cubic_hasse_forces_reducible
        ),
        "all_boundary_cubic_supports_excluded": (
            all_boundary_cubic_supports_excluded
        ),
        "all_p_at_least_31_boundary_cubic_supports_excluded": True,
        "smooth_irreducible_cubic_case_remains": False,
        "high_intersection_conic_excluded": False,
        "aggregate_signed_edge_chain_is_centrally_symmetric": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_bounded_compact_odd_radon_centrality_certificate(
    p: int, compact_count: int
) -> dict[str, object]:
    r"""Force central symmetry for a bounded-compact branch-C row.

    Let ``p=4r+3``, ``r>=7``, and let the signed edge chain consist of
    ``compact_count=b`` arbitrary compact atoms together with ``r-1``
    positive all-equal triangle boundaries.  If ``3*b <= r+2`` and every
    odd contraction through degree ``p-2`` vanishes, then the chain is
    centrally symmetric.

    There are at most ``N=3(r+b-1)`` signed edge occurrences.  Thus the
    degree-``2r-1`` isolation argument applies because ``N<=4r-1``.  A
    nonzero odd/Radon word must be supported on a maximal line in
    ``H x H``.  Its orbit differences satisfy
    ``|n_E| <= (r-1)+2b``.

    This is too few integer projective classes for a horizontal or diagonal
    maximal line.  A vertical line instead gives a constant nonzero integer
    ``k`` on the ``2r+1`` fixed-sum matching orbits.  The occurrence bound
    forces ``|k|=1``.  After fixing its sign, each positive all-equal
    triangle supplies at most one aligned fixed-sum edge and each compact
    triangle at most two: two all-equal pair sums with the same sign, or all
    three compact signed pair sums, would force repeated labels.  Hence only
    ``r-1+2b<2r+1`` aligned occurrences are available, a contradiction.

    The parity version is recorded as an independent weaker check.  At most
    one baseline edge reverses per compact atom because its two negative
    edges share the distinguished vertex.  A triangle contributes either
    zero or two odd vertices to the selected baseline, so strict
    ``3b<r+2`` also gives a direct Eulerian contradiction.  The aligned-edge
    argument covers the boundary case ``3b=r+2`` as well.

    This theorem assumes zero odd global forms.  It does not constrain
    nonzero odd forms, even moments, an integral common edge lift, or the
    signed Boolean box.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    if (
        not isinstance(compact_count, int)
        or isinstance(compact_count, bool)
        or compact_count < 0
    ):
        raise ValueError("compact_count must be a nonnegative integer")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")
    b = compact_count
    maximum_compact_count = (r + 2) // 3
    if 3 * b > r + 2:
        raise ValueError(
            f"the support theorem requires compact_count<={maximum_compact_count}"
        )

    half_square_count = 2 * r + 1
    dual_degree = 2 * r - 1
    support_isolation_bound = 4 * r - 1
    signed_occurrences = 3 * (r + b - 1)
    support_isolation_margin = support_isolation_bound - signed_occurrences
    orbit_difference_bound = r + 2 * b - 1
    required_projective_classes = half_square_count
    available_integer_projective_classes = orbit_difference_bound

    nonaxis_line_intersections = {
        "chi_a=1,chi_b=1": r,
        "chi_a=1,chi_b=-1": r,
        "chi_a=-1,chi_b=1": r,
        "chi_a=-1,chi_b=-1": r + 1,
    }
    maximal_line_count = 3 * half_square_count

    matching_edges = half_square_count
    matching_odd_vertices = p - 1
    negative_compact_occurrences = 2 * b
    maximum_reversed_matching_edges = b
    minimum_after_reversals = (
        matching_odd_vertices - 4 * maximum_reversed_matching_edges
    )
    atom_count = r + b - 1
    triangle_group_parity_support_bound = 2 * atom_count
    triangle_group_parity_margin = (
        minimum_after_reversals - triangle_group_parity_support_bound
    )
    parity_argument_is_strict = triangle_group_parity_margin > 0

    aligned_occurrence_capacity = r - 1 + 2 * b
    aligned_occurrence_deficit = matching_edges - aligned_occurrence_capacity
    two_unit_vertical_deficit = 2 * matching_edges - signed_occurrences
    integer_residue_separation_margin = p - 2 * orbit_difference_bound
    integer_zero_lift_margin = p - orbit_difference_bound

    proved = bool(
        max(range(3, p - 1, 2)) == p - 2
        and dual_degree == (p - 5) // 2
        and 3 * b <= r + 2
        and b <= maximum_compact_count
        and signed_occurrences <= support_isolation_bound
        and support_isolation_margin == r - 3 * b + 2
        and max(nonaxis_line_intersections.values()) == r + 1
        and maximal_line_count == 6 * r + 3
        and orbit_difference_bound < half_square_count
        and 2 * orbit_difference_bound < p
        and two_unit_vertical_deficit > 0
        and maximum_reversed_matching_edges == b
        and triangle_group_parity_margin == 2 * (r + 2 - 3 * b)
        and parity_argument_is_strict == (3 * b < r + 2)
        and aligned_occurrence_capacity == orbit_difference_bound
        and aligned_occurrence_deficit == r - 2 * b + 2
        and aligned_occurrence_deficit > 0
        and integer_zero_lift_margin > 0
    )
    if not proved:
        raise ArithmeticError("the bounded-compact odd/Radon proof changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": b,
        "all_equal_triangle_count": r - 1,
        "atom_profile": f"{b} compact plus {r - 1} all-equal triangles",
        "compact_count_hypothesis": "3*b<=r+2",
        "maximum_compact_count_covered": maximum_compact_count,
        "odd_degrees": list(range(3, p - 1, 2)),
        "dual_polynomial_total_degree": dual_degree,
        "noncollinear_isolation_support_bound": support_isolation_bound,
        "signed_edge_occurrence_bound": signed_occurrences,
        "support_isolation_margin": support_isolation_margin,
        "square_parameter_count": half_square_count,
        "nonaxis_line_intersections": nonaxis_line_intersections,
        "maximal_line_types": ["U=u", "D=d", "D=aU"],
        "maximal_line_count": maximal_line_count,
        "all_equal_atom_orbit_difference_bound": 1,
        "compact_atom_orbit_difference_bound": 2,
        "total_orbit_difference_bound": orbit_difference_bound,
        "required_horizontal_diagonal_projective_classes": (
            required_projective_classes
        ),
        "available_bounded_integer_projective_classes": (
            available_integer_projective_classes
        ),
        "vertical_matching_edges": matching_edges,
        "vertical_matching_odd_vertices": matching_odd_vertices,
        "vertical_constant_residue_is_one_integer": True,
        "vertical_two_unit_l1_deficit": two_unit_vertical_deficit,
        "negative_compact_occurrence_bound": negative_compact_occurrences,
        "maximum_reversed_matching_edges": maximum_reversed_matching_edges,
        "minimum_odd_vertices_after_reversals": minimum_after_reversals,
        "underlying_triangle_count": atom_count,
        "triangle_group_parity_support_bound": (
            triangle_group_parity_support_bound
        ),
        "triangle_group_parity_margin": triangle_group_parity_margin,
        "parity_argument_strictly_closes": parity_argument_is_strict,
        "aligned_fixed_sum_occurrence_capacity": aligned_occurrence_capacity,
        "aligned_fixed_sum_occurrence_deficit": aligned_occurrence_deficit,
        "vertical_line_excluded_by_aligned_incidence": True,
        "integer_residue_separation_margin": integer_residue_separation_margin,
        "integer_zero_lift_margin": integer_zero_lift_margin,
        "aggregate_signed_edge_chain_is_centrally_symmetric": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_odd_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "F_p_common_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_balanced_odd_radon_centrality_band_certificate(
    p: int,
) -> dict[str, object]:
    r"""Map the general-``b`` theorem to the balanced branch-C ``t`` band."""
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")

    direction_count = 2 * r + 2
    lower = 2 * r * r - 4 * r - 2
    upper = 4 * r * r - 2 * r - 5
    compact_cap = (r + 2) // 3
    delta_max = direction_count * compact_cap
    t_max = lower + delta_max
    endpoint = p3_local_survivor(p, t_max)
    endpoint_counts = tuple(
        int(row["Q"]) - r - 2 for row in endpoint["opposite_rows"]
    )
    next_profile = p3_local_survivor(p, t_max + 1)
    next_counts = tuple(
        int(row["Q"]) - r - 2 for row in next_profile["opposite_rows"]
    )
    boundary_theorem = p3_bounded_compact_odd_radon_centrality_certificate(
        p, compact_cap
    )

    proved = bool(
        t_max < upper
        and delta_max == direction_count * compact_cap
        and set(endpoint_counts) == {compact_cap}
        and next_counts.count(compact_cap + 1) == 1
        and next_counts.count(compact_cap) == direction_count - 1
        and boundary_theorem["proved"]
    )
    if not proved:
        raise ArithmeticError("the balanced branch-C odd/Radon band changed")
    return {
        "p": p,
        "r": r,
        "branch_C_t_interval": [lower, upper],
        "balanced_direction_count": direction_count,
        "delta_definition": "delta=t-(2r^2-4r-2)=sum_L b_L",
        "balanced_compact_counts": "floor(delta/m) or ceil(delta/m)",
        "maximum_compact_count_covered": compact_cap,
        "centrality_delta_interval": [0, delta_max],
        "centrality_t_interval": [lower, t_max],
        "endpoint_compact_counts": list(endpoint_counts),
        "first_uncovered_balanced_profile_t": t_max + 1,
        "first_uncovered_balanced_profile_compact_counts": list(next_counts),
        "all_balanced_opposite_rows_central_when_odd_forms_zero": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_odd_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "unbalanced_allocations_ruled_out": False,
        "F_p_common_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "boundary_general_b_certificate": boundary_theorem,
        "no_prime_census": True,
        "proved": proved,
    }


def p3_first_interior_odd_radon_centrality_certificate(
    p: int,
) -> dict[str, object]:
    r"""Compatibility wrapper for the first-interior branch-C row.

    Let ``p=4r+3``, ``r>=7``.  The unique exceptional opposite row at the
    balanced first interior point contains one arbitrary compact atom and
    ``r-1`` all-equal triangles, hence at most ``3r`` signed edge
    occurrences.  If every odd contraction through degree ``p-2`` vanishes,
    pairing an edge with its negative gives

    ``W(U,D)=n_E*D*(s+t)``

    orthogonal on ``H x H`` to all bivariate polynomials of total degree at
    most ``2r-1``.  A word of support at most ``4r-1`` is either zero or
    supported on one maximal line.  Character summation shows that the
    maximal lines are exactly ``U=u``, ``D=d``, and ``D=aU`` for nonzero
    squares ``u,d,a``.

    Horizontal and diagonal lines need all ``2r+1`` projective nonzero
    coefficient classes, while ``|n_E|<=r+1`` supplies at most ``r+1``.
    On a vertical line, l1 mass forces ``|n_E|=1``.  The fixed-sum matching
    has ``4r+2`` odd vertices.  The compact atom has only two negative
    occurrences, so at most two baseline edges reverse, leaving at least
    ``4r-6`` odd vertices.  Removing the aligned baseline leaves l1 mass at
    most ``r-1`` and can toggle at most ``2r-2`` vertices.  This leaves
    ``2r-4>0`` odd vertices, whereas the mod-two sum of the ``r`` triangle
    boundaries is Eulerian.  Therefore the aggregate signed edge chain is
    centrally symmetric.

    The general bounded-compact theorem supplies a stronger aligned-incidence
    proof.  The legacy reversal and residual-mass fields below are retained
    so the original p31 evidence remains replayable.  The same theorem at
    compact count zero covers every pure ``r-1``-triangle opposite row at
    the lower endpoint and the nonexceptional first-interior rows.  It
    assumes zero odd global forms and does not settle nonzero odd forms or
    the even/Boolean gates.
    """
    _check_prime(p)
    if p % 4 != 3:
        raise ValueError("need p=3 mod 4")
    r = (p - 3) // 4
    if r < 7:
        raise ValueError("the all-prime compact-ray audit starts at r=7")

    general = p3_bounded_compact_odd_radon_centrality_certificate(p, 1)
    pure_general = p3_bounded_compact_odd_radon_centrality_certificate(p, 0)

    half_square_count = 2 * r + 1
    dual_degree = 2 * r - 1
    support_isolation_bound = 4 * r - 1
    signed_occurrences = 3 * r
    orbit_difference_bound = r + 1
    required_projective_classes = half_square_count
    available_integer_projective_classes = r + 1

    nonaxis_line_intersections = {
        "chi_a=1,chi_b=1": r,
        "chi_a=1,chi_b=-1": r,
        "chi_a=-1,chi_b=1": r,
        "chi_a=-1,chi_b=-1": r + 1,
    }
    maximal_line_count = 3 * half_square_count

    matching_edges = half_square_count
    matching_odd_vertices = p - 1
    negative_compact_occurrences = 2
    minimum_after_reversals = matching_odd_vertices - 4 * negative_compact_occurrences
    residual_l1 = signed_occurrences - matching_edges
    residual_vertex_toggles = 2 * residual_l1
    minimum_uncorrected_odd_vertices = (
        minimum_after_reversals - residual_vertex_toggles
    )

    pure_occurrences = 3 * (r - 1)
    pure_residual_l1 = pure_occurrences - matching_edges
    pure_minimum_uncorrected_odd_vertices = (
        matching_odd_vertices - 2 * pure_residual_l1
    )

    lower = 2 * r * r - 4 * r - 2
    first_interior = p3_local_survivor(p, lower + 1)
    compact_counts = tuple(
        int(row["Q"]) - r - 2 for row in first_interior["opposite_rows"]
    )
    proved = bool(
        general["proved"]
        and pure_general["proved"]
        and max(range(3, p - 1, 2)) == p - 2
        and dual_degree == (p - 5) // 2
        and signed_occurrences <= support_isolation_bound
        and max(nonaxis_line_intersections.values()) == r + 1
        and r + 1 < half_square_count
        and maximal_line_count == 6 * r + 3
        and orbit_difference_bound == r + 1
        and available_integer_projective_classes < required_projective_classes
        and 2 * matching_edges > signed_occurrences
        and minimum_after_reversals == 4 * r - 6
        and residual_l1 == r - 1
        and residual_vertex_toggles == 2 * r - 2
        and minimum_uncorrected_odd_vertices == 2 * r - 4
        and minimum_uncorrected_odd_vertices > 0
        and pure_residual_l1 == r - 4
        and pure_minimum_uncorrected_odd_vertices == 2 * r + 10
        and compact_counts.count(1) == 1
        and compact_counts.count(0) == 2 * r + 1
    )
    if not proved:
        raise ArithmeticError("the all-prime branch-C odd/Radon proof changed")
    return {
        "p": p,
        "r": r,
        "stronger_general_b_certificate": general,
        "first_interior_t": lower + 1,
        "first_interior_opposite_compact_counts": list(compact_counts),
        "exceptional_atom_profile": f"one compact plus {r - 1} all-equal triangles",
        "odd_degrees": list(range(3, p - 1, 2)),
        "dual_polynomial_total_degree": dual_degree,
        "noncollinear_isolation_support_bound": support_isolation_bound,
        "signed_edge_occurrence_bound": signed_occurrences,
        "square_parameter_count": half_square_count,
        "nonaxis_line_intersections": nonaxis_line_intersections,
        "maximal_line_types": ["U=u", "D=d", "D=aU"],
        "maximal_line_count": maximal_line_count,
        "all_equal_atom_orbit_difference_bound": 1,
        "compact_atom_orbit_difference_bound": 2,
        "total_orbit_difference_bound": orbit_difference_bound,
        "required_horizontal_diagonal_projective_classes": required_projective_classes,
        "available_bounded_integer_projective_classes": available_integer_projective_classes,
        "vertical_matching_edges": matching_edges,
        "vertical_matching_odd_vertices": matching_odd_vertices,
        "negative_compact_occurrence_bound": negative_compact_occurrences,
        "minimum_odd_vertices_after_two_matching_reversals": minimum_after_reversals,
        "remaining_l1_after_aligned_baseline": residual_l1,
        "remaining_vertex_toggle_bound": residual_vertex_toggles,
        "minimum_uncorrected_odd_vertices": minimum_uncorrected_odd_vertices,
        "pure_all_equal_lower_endpoint_corollary": {
            "signed_edge_occurrence_bound": pure_occurrences,
            "total_orbit_difference_bound": r - 1,
            "negative_occurrence_bound": 0,
            "remaining_l1_after_aligned_baseline": pure_residual_l1,
            "minimum_uncorrected_odd_vertices": pure_minimum_uncorrected_odd_vertices,
            "centrally_symmetric": True,
        },
        "first_interior_all_opposite_rows_centrally_symmetric_when_odd_forms_zero": True,
        "aggregate_signed_edge_chain_is_centrally_symmetric": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_odd_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "F_p_common_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "no_prime_census": True,
        "proved": proved,
    }


def p31_odd_ae_radon_symmetry_certificate() -> dict[str, object]:
    r"""Compress all odd rows and force antipodal aggregate edge symmetry.

    For a nonantipodal edge ``{s,t}``, put

    ``U=(s+t)^2, D=(s-t)^2, V=(U-D)/4``.

    Pairing the edge with its negative turns all odd contractions through
    degree 29 into the assertion that ``W(U,D)=n*D*(s+t)`` is orthogonal to
    every bivariate polynomial of total degree at most 13 on
    ``Omega=H x H``, where ``H`` is the set of 15 nonzero squares and the
    integral edge-multiplicity difference satisfies ``|n|<=6``.

    A support of at most 18 for such a dual Reed--Muller word is collinear:
    otherwise, through any support point, at most 13 lines avoiding that
    point cover the other at most 17 points, and their product isolates the
    point.  Every affine line meets ``Omega`` in 0, 7, 8, or 15 points.  The
    15-point lines are precisely ``U=u``, ``D=d``, and ``D=aU`` for nonzero
    squares ``u,d,a``.  The unique orthogonal weights along such a line are
    proportional to its square parameter.

    Horizontal and diagonal support would make ``n`` traverse all 15 sign
    classes, contradicting ``|n|<=6``.  Vertical support makes ``n`` a
    nonzero constant.  Its absolute value must be one; the 15 required
    fixed-sum edges form a matching on 30 odd-degree vertices, and the three
    remaining edge occurrences cannot make every degree even.  Thus ``W=0``
    and the aggregate edge multigraph of the six triangles is centrally
    symmetric.  This does not yet say that the triangles pair blockwise.
    """
    p = 31
    squares = {pow(value, 2, p) for value in range(1, p)}
    omega = {(u_value, d_value) for u_value in squares for d_value in squares}
    intersection_distribution: dict[int, int] = {}
    maximal_lines: list[tuple[object, ...]] = []

    def record(size: int, line: tuple[object, ...]) -> None:
        intersection_distribution[size] = intersection_distribution.get(size, 0) + 1
        if size == 15:
            maximal_lines.append(line)

    for u_value in range(p):
        record(sum(u == u_value for u, _ in omega), ("U", u_value))
    for slope in range(p):
        for intercept in range(p):
            record(
                sum(d == (slope * u + intercept) % p for u, d in omega),
                ("D", slope, intercept),
            )

    expected_maximal = {
        *(("U", value) for value in squares),
        *(("D", 0, value) for value in squares),
        *(("D", value, 0) for value in squares),
    }
    proved = bool(
        len(squares) == 15
        and len(omega) == 225
        and intersection_distribution == {15: 45, 0: 47, 7: 675, 8: 225}
        and set(maximal_lines) == expected_maximal
        and 18 <= 2 * 13
        and 15 > 2 * 6
        and 30 - 2 * 3 > 0
    )
    if not proved:
        raise ArithmeticError("the p=31 odd/Radon support certificate changed")
    return {
        "p": p,
        "odd_degrees": list(range(3, 30, 2)),
        "invariant_coordinates": ["U=(s+t)^2", "D=(s-t)^2", "V=(U-D)/4"],
        "dual_polynomial_total_degree": 13,
        "dual_Reed_Muller_minimum_support": 15,
        "six_triangle_edge_occurrence_bound": 18,
        "omega_size": len(omega),
        "affine_line_intersection_distribution": intersection_distribution,
        "maximal_line_count": len(maximal_lines),
        "edge_difference_absolute_bound": 6,
        "vertical_matching_odd_vertices": 30,
        "remaining_edge_occurrences": 3,
        "aggregate_AE_edge_multigraph_is_centrally_symmetric": True,
        "trianglewise_negation_pairs_forced": False,
        "proved": proved,
    }


def p31_arbitrary_compact_odd_radon_symmetry_certificate() -> dict[str, object]:
    r"""Force central symmetry with one arbitrarily labelled compact atom.

    Add the signed compact edge chain

    ``K(a,b;c)={a,b}-{a,c}-{b,c}``

    to six all-positive triangle boundaries.  It has at most 21 signed edge
    occurrences.  After pairing every nonantipodal edge with its negative,
    the odd rows through degree 29 again say that
    ``W(U,D)=n_E D(s+t)`` is orthogonal to all bivariate polynomials of total
    degree at most 13 on ``Omega=H x H``.

    Any noncollinear support of size at most 27 can be isolated by a product
    of at most 13 affine lines.  Hence a nonzero ``W`` of support at most 21
    must occupy one of the 45 maximal 15-point lines.  A single all-positive
    triangle changes an edge-orbit difference by at most one; the compact
    chain changes it by at most two.  Thus ``|n_E|<=8``.  Horizontal and
    diagonal maximal lines require all 15 projective nonzero field classes,
    which cannot be represented by the eight integer classes ``1,...,8``.

    On a vertical line the unique weights give ``n_E(s+t)=lambda``.  Scale
    the labels so that the fixed edge sum is one.  The total occurrence
    bound forces ``|n_E|=1`` on all 15 points.  Choosing the corresponding
    fixed-sum matching gives 30 odd vertices.  Only the two
    negative compact occurrences can reverse matching edges, and two such
    reversals leave at least 22 odd vertices.  After subtracting this signed
    baseline, the residual chain has coefficient ``l1`` mass at most six,
    so its parity support toggles at most 12 vertices.  This contradicts the
    fact that the compact
    triangle and six all-positive triangles form seven triangle boundaries,
    hence have even degree at every vertex modulo two.  Therefore ``W=0``.

    The conclusion concerns the *total signed edge chain*.  It neither makes
    the compact atom centered nor supplies the joint degree-six/eight lift.
    """
    p = 31

    def edge(first: int, second: int) -> tuple[int, int]:
        return tuple(sorted((first % p, second % p)))

    def negative(edge_value: tuple[int, int]) -> tuple[int, int]:
        return edge(-edge_value[0], -edge_value[1])

    edge_representatives = tuple(combinations(range(p), 2))

    def orbit_differences(
        terms: tuple[tuple[tuple[int, int], int], ...],
    ) -> dict[tuple[int, int], int]:
        coefficients: dict[tuple[int, int], int] = {}
        for edge_value, coefficient in terms:
            coefficients[edge_value] = coefficients.get(edge_value, 0) + coefficient
        differences: dict[tuple[int, int], int] = {}
        for edge_value in edge_representatives:
            negative_edge = negative(edge_value)
            if edge_value < negative_edge:
                differences[edge_value] = (
                    coefficients.get(edge_value, 0)
                    - coefficients.get(negative_edge, 0)
                )
        return differences

    ae_bound = 0
    for triangle in combinations(range(p), 3):
        terms = tuple((edge(*pair), 1) for pair in combinations(triangle, 2))
        ae_bound = max(ae_bound, *(abs(value) for value in orbit_differences(terms).values()))

    compact_bound = 0
    compact_overlap_witness: tuple[int, int, int] | None = None
    for positive_pair in combinations(range(p), 2):
        for distinguished in range(p):
            if distinguished in positive_pair:
                continue
            first, second = positive_pair
            terms = (
                (edge(first, second), 1),
                (edge(first, distinguished), -1),
                (edge(second, distinguished), -1),
            )
            local_bound = max(abs(value) for value in orbit_differences(terms).values())
            if local_bound > compact_bound:
                compact_bound = local_bound
                compact_overlap_witness = (first, second, distinguished)

    squares = {pow(value, 2, p) for value in range(1, p)}
    omega = {(u_value, d_value) for u_value in squares for d_value in squares}
    intersection_distribution: dict[int, int] = {}
    maximal_lines: set[tuple[object, ...]] = set()

    def record_line(size: int, line: tuple[object, ...]) -> None:
        intersection_distribution[size] = intersection_distribution.get(size, 0) + 1
        if size == 15:
            maximal_lines.add(line)

    for u_value in range(p):
        record_line(sum(u == u_value for u, _ in omega), ("U", u_value))
    for slope in range(p):
        for intercept in range(p):
            record_line(
                sum(d == (slope * u + intercept) % p for u, d in omega),
                ("D", slope, intercept),
            )

    expected_maximal = {
        *(("U", value) for value in squares),
        *(("D", 0, value) for value in squares),
        *(("D", value, 0) for value in squares),
    }

    def matching(total: int) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(
                {
                    edge(first, total - first)
                    for first in range(p)
                    if first != (total - first) % p
                }
            )
        )

    baseline = matching(1)
    opposite = tuple(negative(edge_value) for edge_value in baseline)
    minimum_odd_vertices = p
    for switch_count in range(3):
        for switched in combinations(range(len(baseline)), switch_count):
            selected = list(baseline)
            for index in switched:
                selected[index] = opposite[index]
            parity = [0] * p
            for first, second in selected:
                parity[first] ^= 1
                parity[second] ^= 1
            minimum_odd_vertices = min(minimum_odd_vertices, sum(parity))

    integer_projective_classes = {
        min(value % p, (-value) % p) for value in range(-8, 9) if value
    }
    required_projective_counts = {
        len(
            {
                min((scale * root) % p, (-scale * root) % p)
                for root in range(1, 16)
            }
        )
        for scale in range(1, p)
    }

    # Centrality itself does not force the compact atom to be centered.
    witness_compact = (0, 1, 2)
    witness_ae = (
        (0, 2, 30),
        (0, 2, 30),
        (0, 1, 2),
        (0, 1, 29),
        (0, 29, 30),
        (1, 2, 29),
    )
    witness_terms = (
        (edge(witness_compact[0], witness_compact[1]), 1),
        (edge(witness_compact[0], witness_compact[2]), -1),
        (edge(witness_compact[1], witness_compact[2]), -1),
        *tuple(
            (edge(*pair), 1)
            for triangle in witness_ae
            for pair in combinations(triangle, 2)
        ),
    )
    witness_central = not any(orbit_differences(witness_terms).values())
    witness_odd_zero = all(
        (
            compact_moment(p, *witness_compact, degree, channel)
            + sum(all_equal_moment(p, *triangle, degree, channel) for triangle in witness_ae)
        )
        % p
        == 0
        for degree in range(3, 30, 2)
        for channel in range(degree // 2)
    )
    witness_degree_six = tuple(
        (
            compact_moment(p, *witness_compact, 6, channel)
            + sum(all_equal_moment(p, *triangle, 6, channel) for triangle in witness_ae)
        )
        % p
        for channel in range(3)
    )
    witness_degree_eight = tuple(
        (
            compact_moment(p, *witness_compact, 8, channel)
            + sum(all_equal_moment(p, *triangle, 8, channel) for triangle in witness_ae)
        )
        % p
        for channel in range(4)
    )

    remaining_occurrences = 21 - 15
    proved = bool(
        ae_bound == 1
        and compact_bound == 2
        and compact_overlap_witness == (0, 1, 30)
        and 6 * ae_bound + compact_bound == 8
        and 21 <= 27
        and intersection_distribution == {15: 45, 0: 47, 7: 675, 8: 225}
        and maximal_lines == expected_maximal
        and required_projective_counts == {15}
        and len(integer_projective_classes) == 8
        and len(baseline) == 15
        and minimum_odd_vertices == 22
        and remaining_occurrences == 6
        and minimum_odd_vertices - 2 * remaining_occurrences == 10
        and witness_central
        and witness_odd_zero
        and witness_degree_six == (26, 26, 5)
        and witness_degree_eight == (13, 4, 30, 6)
    )
    if not proved:
        raise ArithmeticError("the arbitrary-compact p=31 odd/Radon certificate changed")
    return {
        "p": p,
        "odd_degrees": list(range(3, 30, 2)),
        "dual_polynomial_total_degree": 13,
        "noncollinear_isolation_support_bound": 27,
        "signed_edge_occurrence_bound": 21,
        "all_equal_atom_orbit_difference_bound": ae_bound,
        "compact_atom_orbit_difference_bound": compact_bound,
        "total_orbit_difference_bound": 6 * ae_bound + compact_bound,
        "compact_overlap_witness": list(compact_overlap_witness),
        "affine_line_intersection_distribution": intersection_distribution,
        "maximal_line_count": len(maximal_lines),
        "required_horizontal_diagonal_projective_classes": 15,
        "available_bounded_integer_projective_classes": len(integer_projective_classes),
        "vertical_matching_edges": len(baseline),
        "negative_compact_occurrence_bound": 2,
        "minimum_odd_vertices_after_two_matching_reversals": minimum_odd_vertices,
        "remaining_edge_occurrences": remaining_occurrences,
        "minimum_uncorrected_odd_vertices": minimum_odd_vertices - 2 * remaining_occurrences,
        "total_signed_edge_chain_is_centrally_symmetric": True,
        "centrality_forces_centered_compact": False,
        "noncentered_centrality_witness": {
            "compact": list(witness_compact),
            "all_equal_triangles": [list(triangle) for triangle in witness_ae],
            "all_odd_rows_zero": witness_odd_zero,
            "degree_six_vector": list(witness_degree_six),
            "degree_eight_vector": list(witness_degree_eight),
        },
        "joint_degree_six_eight_zero_ruled_out": False,
        "proved": proved,
    }


def p31_antipodal_pasch_joint_no_go() -> dict[str, object]:
    r"""Exclude every antipodal Pasch-four plus two symmetric blocks.

    One explicit parametrization of the positive Pasch leg is

    ``(a,-a,b), (a,-b,c), (-a,-b,-c), (b,c,-c)``.

    Its negative is the other leg.  The scan includes all seven fixed-point-
    free involutions from the standard Pasch leg to its mate, allows quotient
    collisions and repeated positive blocks, and then tries both possible
    two-block symmetric remainders: one pair ``T,-T`` or two individually
    invariant triangles ``(0,u,-u)``.
    """
    import hashlib

    pasch = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))
    mate = {
        frozenset(triangle)
        for triangle in ((0, 1, 3), (0, 2, 4), (1, 2, 5), (3, 4, 5))
    }
    involutions: list[tuple[tuple[int, ...], tuple[tuple[int, int], ...]]] = []
    for permutation in permutations(range(6)):
        if not all(permutation[permutation[index]] == index for index in range(6)):
            continue
        image = {
            frozenset(permutation[index] for index in triangle)
            for triangle in pasch
        }
        if image != mate:
            continue
        cycles = tuple(
            (index, permutation[index])
            for index in range(6)
            if index < permutation[index]
        )
        if len(cycles) == 3:
            involutions.append((permutation, cycles))

    parameter_assignments = 0
    valid_parameter_assignments = 0
    configurations: set[tuple[tuple[int, int, int], ...]] = set()
    for _, cycles in involutions:
        for scales in product(range(31), repeat=3):
            parameter_assignments += 1
            labels = [0] * 6
            for (first, second), scale in zip(cycles, scales):
                labels[first], labels[second] = scale, (-scale) % 31
            concrete = tuple(
                sorted(tuple(sorted(labels[index] for index in triangle)) for triangle in pasch)
            )
            if not all(len(set(triangle)) == 3 for triangle in concrete):
                continue
            valid_parameter_assignments += 1
            configurations.add(concrete)

    triples, vectors = _p31_ae_even_catalog()
    vector_by_triangle = dict(zip(triples, vectors))
    pair_unit_codes = {
        _p31_vector_encode(_p31_vector_scale(2, vector)) for vector in vectors
    }
    invariant_vectors = tuple(
        vector_by_triangle[tuple(sorted((0, scale, -scale % 31)))]
        for scale in range(1, 16)
    )
    invariant_two_codes = {
        _p31_vector_encode(_p31_vector_add(invariant_vectors[first], invariant_vectors[second]))
        for first, second in combinations_with_replacement(range(15), 2)
    }
    target = _p31_vector_scale(-1, _p31_centered_compact_even_vector())
    transcript = hashlib.sha256()
    pair_solutions = 0
    invariant_solutions = 0
    for concrete in sorted(configurations):
        total = (0,) * 7
        for triangle in concrete:
            total = _p31_vector_add(total, vector_by_triangle[triangle])
        residual = _p31_vector_sub(target, total)
        transcript.update(bytes(value for triangle in concrete for value in triangle))
        transcript.update(_p31_vector_encode(residual).to_bytes(8, "little"))
        residual_code = _p31_vector_encode(residual)
        pair_solutions += residual_code in pair_unit_codes
        invariant_solutions += residual_code in invariant_two_codes

    reference_internal = 0
    reference_simple = 0
    reference_injective = 0
    for a_value, b_value, c_value in product(range(31), repeat=3):
        concrete = tuple(
            tuple(sorted(triangle))
            for triangle in (
                (a_value, -a_value % 31, b_value),
                (a_value, -b_value % 31, c_value),
                (-a_value % 31, -b_value % 31, -c_value % 31),
                (b_value, c_value, -c_value % 31),
            )
        )
        if not all(len(set(triangle)) == 3 for triangle in concrete):
            continue
        reference_internal += 1
        reference_simple += len(set(concrete)) == 4
        reference_injective += len(
            {a_value, -a_value % 31, b_value, -b_value % 31, c_value, -c_value % 31}
        ) == 6

    evidence_sha256 = transcript.hexdigest()
    proved = bool(
        len(involutions) == 7
        and parameter_assignments == 208_537
        and valid_parameter_assignments == 165_660
        and len(configurations) == 6_910
        and (reference_internal, reference_simple, reference_injective)
        == (23_550, 22_680, 21_840)
        and len(pair_unit_codes) == 2_255
        and len(invariant_two_codes) == 120
        and pair_solutions == invariant_solutions == 0
        and evidence_sha256
        == "40889ecbc7e92660d045e547a7f532b1aaa1dcf5519c9185ef02f0f3eea910ce"
    )
    if not proved:
        raise ArithmeticError("the p=31 antipodal Pasch certificate changed")
    return {
        "p": 31,
        "reference_parametrization": [
            ["a", "-a", "b"],
            ["a", "-b", "c"],
            ["-a", "-b", "-c"],
            ["b", "c", "-c"],
        ],
        "fixed_point_free_involutions": len(involutions),
        "parameter_assignments": parameter_assignments,
        "internally_valid_assignments": valid_parameter_assignments,
        "distinct_concrete_block_multisets": len(configurations),
        "reference_assignment_counts": {
            "internally_valid": reference_internal,
            "four_distinct_blocks": reference_simple,
            "six_distinct_labels": reference_injective,
        },
        "pasch_plus_negation_pair_solution_count": pair_solutions,
        "pasch_plus_two_invariant_triangles_solution_count": invariant_solutions,
        "residual_transcript_sha256": evidence_sha256,
        "proved": proved,
    }


_P31_VOLUME_SIX_TRADE_TYPES = {
    "6-cycle": (
        ("123", "145", "167", "834", "856", "872"),
        ("134", "156", "172", "823", "845", "867"),
    ),
    "semihead": (
        ("127", "136", "145", "235", "246", "347"),
        ("126", "135", "147", "237", "245", "346"),
    ),
    "trade-X": (
        ("123", "124", "156", "256", "345", "346"),
        ("125", "126", "134", "234", "356", "456"),
    ),
    "trade-Y": (
        ("124", "125", "136", "137", "267", "345"),
        ("126", "127", "134", "135", "245", "367"),
    ),
}


def p31_volume_six_trade_joint_no_go() -> dict[str, object]:
    r"""Exclude the four simple volume-six trade types and all quotients.

    For every template isomorphism ``pi:T+ -> T-``, this enumerates all
    solutions of ``label(pi(v))=-label(v)``.  Cross-vertex collisions and
    repeated concrete blocks are allowed, while every individual triangle
    must retain three distinct labels.  Thus the finite scan contains the
    injective involutive embeddings required by the abstract simple-trade
    classification and additionally checks all such template quotients.
    """
    import hashlib

    def abstract_blocks(raw: tuple[str, ...]) -> tuple[tuple[int, int, int], ...]:
        return tuple(tuple(sorted(int(character) - 1 for character in word)) for word in raw)

    def template_isomorphisms(
        left: tuple[tuple[int, int, int], ...],
        right: tuple[tuple[int, int, int], ...],
        point_count: int,
    ) -> tuple[tuple[int, ...], ...]:
        right_set = set(right)
        return tuple(
            permutation
            for permutation in permutations(range(point_count))
            if {
                tuple(sorted(permutation[index] for index in triangle))
                for triangle in left
            }
            == right_set
        )

    def permutation_cycles(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
        cycles: list[tuple[int, ...]] = []
        seen: set[int] = set()
        for start in range(len(permutation)):
            if start in seen:
                continue
            cycle: list[int] = []
            point = start
            while point not in seen:
                seen.add(point)
                cycle.append(point)
                point = permutation[point]
            cycles.append(tuple(cycle))
        return tuple(cycles)

    triples, vectors = _p31_ae_even_catalog()
    vector_by_triangle = dict(zip(triples, vectors))
    target = _p31_vector_scale(-1, _p31_centered_compact_even_vector())
    seen_configurations: set[bytes] = set()
    generated = 0
    valid = 0
    solutions = 0
    per_type: dict[str, dict[str, int]] = {}

    for name, (raw_left, raw_right) in _P31_VOLUME_SIX_TRADE_TYPES.items():
        left = abstract_blocks(raw_left)
        right = abstract_blocks(raw_right)
        point_count = 1 + max(index for triangle in left + right for index in triangle)
        isomorphisms = template_isomorphisms(left, right, point_count)
        type_generated = 0
        type_valid = 0
        type_new = 0
        for permutation in isomorphisms:
            cycles = permutation_cycles(permutation)
            even_cycles = tuple(cycle for cycle in cycles if len(cycle) % 2 == 0)
            odd_cycles = tuple(cycle for cycle in cycles if len(cycle) % 2 == 1)
            for parameters in product(range(31), repeat=len(even_cycles)):
                generated += 1
                type_generated += 1
                labels = [0] * point_count
                for cycle, parameter in zip(even_cycles, parameters):
                    for position, point in enumerate(cycle):
                        labels[point] = parameter if position % 2 == 0 else -parameter % 31
                for cycle in odd_cycles:
                    for point in cycle:
                        labels[point] = 0
                concrete = tuple(
                    sorted(tuple(sorted(labels[index] for index in triangle)) for triangle in left)
                )
                if not all(len(set(triangle)) == 3 for triangle in concrete):
                    continue
                valid += 1
                type_valid += 1
                key = bytes(value for triangle in concrete for value in triangle)
                if key in seen_configurations:
                    continue
                seen_configurations.add(key)
                type_new += 1
                total = tuple(
                    sum(vector_by_triangle[triangle][channel] for triangle in concrete) % 31
                    for channel in range(7)
                )
                solutions += total == target
        per_type[name] = {
            "template_isomorphisms": len(isomorphisms),
            "parameter_assignments": type_generated,
            "internally_valid_assignments": type_valid,
            "new_concrete_block_multisets": type_new,
        }

    transcript = hashlib.sha256()
    for key in sorted(seen_configurations):
        transcript.update(key)
    evidence_sha256 = transcript.hexdigest()
    expected_per_type = {
        "6-cycle": {
            "template_isomorphisms": 12,
            "parameter_assignments": 2_803_392,
            "internally_valid_assignments": 2_007_780,
            "new_concrete_block_multisets": 164_375,
        },
        "semihead": {
            "template_isomorphisms": 24,
            "parameter_assignments": 14_424,
            "internally_valid_assignments": 2_520,
            "new_concrete_block_multisets": 105,
        },
        "trade-X": {
            "template_isomorphisms": 24,
            "parameter_assignments": 190_464,
            "internally_valid_assignments": 131_040,
            "new_concrete_block_multisets": 5_460,
        },
        "trade-Y": {
            "template_isomorphisms": 8,
            "parameter_assignments": 33_728,
            "internally_valid_assignments": 23_520,
            "new_concrete_block_multisets": 0,
        },
    }
    proved = bool(
        per_type == expected_per_type
        and generated == 3_042_008
        and valid == 2_164_860
        and len(seen_configurations) == 169_940
        and solutions == 0
        and evidence_sha256
        == "78ee0fc05757a9d332a8d2da3605a921b28207aafca746c806bf17e043f26dd0"
    )
    if not proved:
        raise ArithmeticError("the p=31 volume-six trade certificate changed")
    return {
        "p": 31,
        "target": list(target),
        "trade_types": list(_P31_VOLUME_SIX_TRADE_TYPES),
        "per_type": per_type,
        "parameter_assignments": generated,
        "internally_valid_assignments": valid,
        "distinct_concrete_block_multisets": len(seen_configurations),
        "solution_count": solutions,
        "concrete_multiset_transcript_sha256": evidence_sha256,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p31_centered_compact_six_ae_odd_six_eight_no_go() -> dict[str, object]:
    r"""Exhaust the zero-global-form first-interior centered-compact block.

    Odd degrees 3 through 29 force the six-AE aggregate edge multigraph to
    be centrally symmetric.  Cancel common blocks between a six-triangle
    multiset and its negative.  A nonempty reduced 2-(v,3) trade has volume
    at least four; volume four is the Pasch trade, volume five does not
    exist, and the simple volume-six trades are the four standard types.

    No reduced trade of volume at most six can have a repeated block.  If a
    positive block ``T`` has multiplicity ``m>=2``, the negative leg needs
    at least ``m`` different blocks through each of the three pairs of T, so
    its volume is at least ``3m``.  Equality can only be ``m=2, volume=6``.
    Then all six negative blocks contain one pair of T and give the three
    vertices total incidence 12.  On the positive side, the two copies of T
    plus four blocks containing at most one T-vertex give incidence at most
    10, a contradiction.

    The remaining cases are therefore: a blockwise symmetric six-block core;
    a Pasch four-trade plus a two-block symmetric core; or one of the four
    simple volume-six trades.  The three exact finite certificates below
    exclude every case for the joint degree-six/eight target.  Scope remains
    fixed centered compact and zero rowwise global forms.
    """
    radon = p31_odd_ae_radon_symmetry_certificate()
    symmetric_core = p31_centered_compact_three_ae_pairs_joint_no_go()
    pasch = p31_antipodal_pasch_joint_no_go()
    volume_six = p31_volume_six_trade_joint_no_go()
    proved = bool(
        radon["proved"]
        and symmetric_core["proved"]
        and pasch["proved"]
        and volume_six["proved"]
    )
    if not proved:
        raise ArithmeticError("the exhaustive p=31 centered block certificate changed")
    return {
        "p": 31,
        "t": 69,
        "fixed_compact_triangle": list(_P31_FIRST_INTERIOR_COMPACT),
        "required_odd_degrees": list(range(3, 30, 2)),
        "required_even_degrees": [6, 8],
        "odd_Radon_compression": radon,
        "reduced_trade_classification": {
            "minimum_nonempty_volume": 4,
            "volume_4": "unique Pasch trade",
            "volume_5": "does not exist",
            "volume_6": list(_P31_VOLUME_SIX_TRADE_TYPES),
            "repeated_blocks_at_volume_at_most_6": False,
        },
        "blockwise_symmetric_core_certificate": symmetric_core,
        "pasch_remainder_certificate": pasch,
        "volume_six_certificate": volume_six,
        "centered_compact_plus_six_AE_can_have_zero_odd_d6_d8_rows": False,
        "noncentered_compact_ruled_out": False,
        "nonzero_global_forms_ruled_out": False,
        "coordinated_changes_on_other_rows_ruled_out": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p31_arbitrary_compact_six_ae_odd_six_eight_no_go() -> dict[str, object]:
    r"""Close the arbitrary-compact local zero-form block at ``p=31,t=69``.

    The odd/Radon lemma first forces the total signed edge chain of one
    arbitrary compact atom and six all-equal atoms to be centrally symmetric.
    Multiplication of every label by ``F_31^*`` leaves centrality and zero
    degree-six/eight moments invariant.  The 13,485 labelled compact atoms
    split into 450 scaling orbits: 449 free orbits of size 30 and the unique
    centered orbit of size 15.

    For each noncentered representative, the exact sparse DFS branches on a
    nonzero edge-orbit imbalance.  Its ``l1`` and maximum-coordinate bounds
    are necessary.  Once a partial chain becomes central, at most five
    triangles remain; the complete volume-at-most-five trade classification
    says that they are an invariant/negation-pair core, or a Pasch four-trade
    possibly plus one invariant triangle.  The seven even channels are tracked
    throughout.  With one or two triangles left, doubled-moment membership in
    the exact one-pair/two-pair catalogs is an additional necessary prune.

    The archived run exhausts all 450 representatives.  Every noncentered
    orbit is infeasible; index 435 is the centered compact ``(1,-1;0)`` and
    is handed to :func:`p31_centered_compact_six_ae_odd_six_eight_no_go`.
    This proves only the local rowwise-zero odd/degree-six/degree-eight gate.
    Nonzero global forms and coordinated changes in other directions remain.
    """
    import hashlib
    import json
    import re
    from pathlib import Path

    p = 31
    root = Path(__file__).resolve().parents[1]
    source_path = root / "evidence" / "p31_arbitrary_compact_fiber.cpp"
    log_path = root / "evidence" / "p31_arbitrary_compact_fiber_v2_000_450.log"
    merge_path = root / "evidence" / "p31_arbitrary_compact_fiber_v2_merge.json"
    comparison_path = (
        root / "evidence" / "p31_arbitrary_compact_fiber_v1_v2_overlap_comparison.json"
    )
    overlap_path = (
        root / "evidence" / "p31_arbitrary_compact_fiber_v1_v2_overlap_0_434.txt"
    )
    v1_log_root = root / "evidence" / "p31_arbitrary_compact_fiber_v1_overlap"
    source_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    log_sha256 = hashlib.sha256(log_path.read_bytes()).hexdigest()
    merge_file_sha256 = hashlib.sha256(merge_path.read_bytes()).hexdigest()

    representatives: set[tuple[int, int, int]] = set()
    labelled_compacts = 0
    for first, second in combinations(range(p), 2):
        for distinguished in range(p):
            if distinguished in (first, second):
                continue
            labelled_compacts += 1
            orbit = []
            for scale in range(1, p):
                scaled_first, scaled_second = sorted(
                    (first * scale % p, second * scale % p)
                )
                orbit.append(
                    (scaled_first, scaled_second, distinguished * scale % p)
                )
            representatives.add(min(orbit))
    ordered_representatives = tuple(sorted(representatives))

    header = (
        "compact_orbits=450 range=0:450 invariant=15 "
        "fixed_sums=1,15,120,535,925,961, pair_units=2255 "
        "pair_pairs=2543460 pasch=3725"
    )
    unsat_pattern = re.compile(
        r"UNSAT_INDEX=(\d+) compact=(\d+),(\d+),(\d+) cumulative_nodes=(\d+)"
    )
    centered_pattern = re.compile(
        r"CENTERED_CERTIFIED_INDEX=(\d+) compact=(\d+),(\d+),(\d+)"
    )
    lines = log_path.read_text().splitlines()
    statuses: dict[int, tuple[str, tuple[int, int, int]]] = {}
    cumulative_nodes: list[int] = []
    for line in lines[1:-1]:
        unsat_match = unsat_pattern.fullmatch(line)
        centered_match = centered_pattern.fullmatch(line)
        if unsat_match is not None:
            index, first, second, distinguished, node_count = map(
                int, unsat_match.groups()
            )
            statuses[index] = ("UNSAT", (first, second, distinguished))
            cumulative_nodes.append(node_count)
        elif centered_match is not None:
            index, first, second, distinguished = map(int, centered_match.groups())
            statuses[index] = ("CENTERED", (first, second, distinguished))
        else:
            raise ArithmeticError(f"unrecognized compact-fiber log line: {line}")
    normalized_statuses = "\n".join(
        f"{index}:{statuses[index][0]}:{','.join(map(str, statuses[index][1]))}"
        for index in range(450)
    ).encode()
    normalized_sha256 = hashlib.sha256(normalized_statuses).hexdigest()

    merge_record = json.loads(merge_path.read_text())
    recorded_payload_sha256 = merge_record.pop("sha256_without_hash_field")
    merge_payload = json.dumps(
        merge_record, sort_keys=True, separators=(",", ":")
    ).encode()
    computed_payload_sha256 = hashlib.sha256(merge_payload).hexdigest()

    comparison_file_sha256 = hashlib.sha256(comparison_path.read_bytes()).hexdigest()
    comparison_record = json.loads(comparison_path.read_text())
    comparison_payload_sha256 = comparison_record.pop("sha256_without_hash_field")
    computed_comparison_payload_sha256 = hashlib.sha256(
        json.dumps(comparison_record, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    overlap_file_sha256 = hashlib.sha256(overlap_path.read_bytes()).hexdigest()
    archived_overlap = overlap_path.read_text().rstrip("\n")
    archived_overlap_sha256 = hashlib.sha256(archived_overlap.encode()).hexdigest()
    v1_statuses: dict[int, tuple[str, tuple[int, int, int]]] = {}
    v1_log_hashes_match = True
    for item in comparison_record["v1_logs"]:
        archived_log = v1_log_root / Path(item["path"]).name
        v1_log_hashes_match &= (
            hashlib.sha256(archived_log.read_bytes()).hexdigest() == item["sha256"]
        )
        for line in archived_log.read_text().splitlines()[1:]:
            unsat_match = unsat_pattern.fullmatch(line)
            centered_match = centered_pattern.fullmatch(line)
            if unsat_match is not None:
                index, first, second, distinguished, _ = map(int, unsat_match.groups())
                status = ("UNSAT", (first, second, distinguished))
            elif centered_match is not None:
                index, first, second, distinguished = map(int, centered_match.groups())
                status = ("CENTERED", (first, second, distinguished))
            elif line.startswith("SUMMARY "):
                continue
            else:
                raise ArithmeticError(f"unrecognized archived v1 log line: {line}")
            if index in v1_statuses and v1_statuses[index] != status:
                raise ArithmeticError(f"inconsistent archived v1 status at index {index}")
            v1_statuses[index] = status
    normalized_v1_overlap = "\n".join(
        f"{index}:{v1_statuses[index][0]}:{','.join(map(str, v1_statuses[index][1]))}"
        for index in range(435)
    )
    normalized_v2_overlap = "\n".join(
        f"{index}:{statuses[index][0]}:{','.join(map(str, statuses[index][1]))}"
        for index in range(435)
    )

    odd_radon = p31_arbitrary_compact_odd_radon_symmetry_certificate()
    centered = p31_centered_compact_six_ae_odd_six_eight_no_go()
    expected_source_sha256 = (
        "1dcfce7b5765630655d049413c4d9138c544a6d05fe19e3308a9a20a2880d1f2"
    )
    expected_log_sha256 = (
        "f3f77607181287095aa69644649d14d7b9b5e3a8f24044477b667549ef0512e3"
    )
    expected_merge_file_sha256 = (
        "c7f5dea5811a8d2aa25d7bd3224b1fceae3fce73bb49fd4c8fe3f335e2e71c2f"
    )
    expected_normalized_sha256 = (
        "ad3bf3c97b378c9cdebb0b77d486cced544199750ad689060bd2a24f6a2210cb"
    )
    expected_payload_sha256 = (
        "efcab50a9f0c67bb00aa6e11a53959205f4213f266072837f1f50fe87ef86459"
    )
    independent_regression_verified = bool(
        comparison_file_sha256
        == "32b1d64679b239fbc001eaf0182cb1978183880547a84c81541cee63a8483ce7"
        and comparison_payload_sha256
        == "b0fbbba6f26ef2f2579b6e72f8656c963091a46b99e21e5c5e30f55276f01890"
        and computed_comparison_payload_sha256 == comparison_payload_sha256
        and overlap_file_sha256
        == "86885254554eef4d616fc20aa937b1bde3fba35fb3cb2b8203d9454d5aad8d73"
        and archived_overlap_sha256
        == "8b6b6277cb63561f744865ecc6aa7012dacc20be7d062c6b53d8670cfd7d75fd"
        and v1_log_hashes_match
        and set(v1_statuses) == set(range(435))
        and normalized_v1_overlap == archived_overlap == normalized_v2_overlap
        and comparison_record["overlap_index_count"] == 435
        and comparison_record["disagreement_count"] == 0
        and comparison_record["all_v1_v2_overlap_verdicts_and_compacts_equal"]
        and comparison_record["normalized_overlap_status_sha256"]
        == archived_overlap_sha256
    )
    if not independent_regression_verified:
        raise ArithmeticError("the independent p=31 v1/v2 overlap audit changed")
    proved = bool(
        odd_radon["proved"]
        and centered["proved"]
        and labelled_compacts == 13_485
        and len(ordered_representatives) == 450
        and ordered_representatives[435] == (1, 30, 0)
        and lines[0] == header
        and lines[-1]
        == "SUMMARY range=0:450 compact_done=450 nodes=317916856 found=0"
        and len(statuses) == 450
        and all(statuses[index][1] == ordered_representatives[index] for index in range(450))
        and sum(status == "UNSAT" for status, _ in statuses.values()) == 449
        and statuses[435] == ("CENTERED", (1, 30, 0))
        and cumulative_nodes == sorted(set(cumulative_nodes))
        and cumulative_nodes[-1] == 317_916_856
        and source_sha256 == expected_source_sha256
        and log_sha256 == expected_log_sha256
        and merge_file_sha256 == expected_merge_file_sha256
        and normalized_sha256 == expected_normalized_sha256
        and recorded_payload_sha256 == expected_payload_sha256
        and computed_payload_sha256 == expected_payload_sha256
        and merge_record["solver_source_sha256"] == expected_source_sha256
        and merge_record["x86_64_binary_sha256"]
        == "4622dbcb2afbfdc4c0da3588e42cba86542a98763fb89429bed7ea2185915955"
        and merge_record["status_counts"]
        == {"UNSAT": 449, "CENTERED": 1, "SAT": 0}
        and merge_record["normalized_merged_status_sha256"]
        == expected_normalized_sha256
        and merge_record["executed_nodes_including_overlap"] == 317_916_856
        and merge_record["complete_coverage"]
        and merge_record["all_noncentered_fibers_UNSAT"]
        and merge_record["proved"]
    )
    if not proved:
        raise ArithmeticError("the p=31 arbitrary-compact finite certificate changed")
    return {
        "p": p,
        "t": 69,
        "local_atom_profile": "one arbitrary compact plus six all-equal triangles",
        "required_odd_degrees": list(range(3, 30, 2)),
        "required_even_degrees": [6, 8],
        "odd_Radon_compression": odd_radon,
        "labelled_compact_atoms": labelled_compacts,
        "compact_scaling_orbits": len(ordered_representatives),
        "noncentered_scaling_orbits_UNSAT": 449,
        "unique_centered_orbit": {
            "index": 435,
            "compact": [1, 30, 0],
            "handled_by_centered_certificate": True,
        },
        "central_remainder_catalog_sizes": {
            "invariant_triangles": 15,
            "fixed_sums_0_through_5": [1, 15, 120, 535, 925, 961],
            "pair_units": 2_255,
            "pair_pairs": 2_543_460,
            "Pasch_even_vectors": 3_725,
        },
        "executed_DFS_nodes": 317_916_856,
        "source_sha256": source_sha256,
        "binary_sha256": merge_record["x86_64_binary_sha256"],
        "raw_log_sha256": log_sha256,
        "normalized_status_sha256": normalized_sha256,
        "merge_file_sha256": merge_file_sha256,
        "merge_payload_sha256": recorded_payload_sha256,
        "independent_v1_v2_regression": {
            "proof_premise": False,
            "overlap_indices": [0, 434],
            "overlap_index_count": 435,
            "disagreement_count": 0,
            "normalized_overlap_status_sha256": archived_overlap_sha256,
            "normalized_overlap_file_sha256": overlap_file_sha256,
            "comparison_file_sha256": comparison_file_sha256,
            "comparison_payload_sha256": comparison_payload_sha256,
            "archived_v1_log_count": len(comparison_record["v1_logs"]),
            "proved": independent_regression_verified,
        },
        "arbitrary_compact_plus_six_AE_can_have_zero_odd_d6_d8_rows": False,
        "nonzero_global_forms_ruled_out": False,
        "coordinated_changes_on_other_rows_ruled_out": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p31_balanced_zero_form_band_certificate() -> dict[str, object]:
    r"""Reuse the one-compact local no-go on balanced ``69<=t<=99``.

    For ``p=31`` the lower branch-C endpoint is ``t_0=68`` and there are
    sixteen opposite directions.  Their deterministic balanced compact
    counts distribute ``delta=t-68`` as floor/ceiling values.  Exactly for
    ``1<=delta<=31`` at least one row has compact count one.  Every such row
    has the same atom profile as the archived exhaustive certificate: one
    arbitrary compact atom and six positive all-equal triangles.  Therefore
    simultaneous zero odd, degree-six, and degree-eight global forms are
    incompatible with each of these balanced local profiles.

    The ``t=69`` stored in the exhaustive wrapper records where the profile
    was first encountered; the solver and its hashes depend only on the
    one-compact/six-AE row.  This corollary does not cover unbalanced row
    allocations, nonzero global forms, or a Boolean common-edge lift.
    """
    p = 31
    r = 7
    lower = 68
    direction_count = 16
    first_t = 69
    last_t = 99
    local_no_go = p31_arbitrary_compact_six_ae_odd_six_eight_no_go()

    profiles: dict[str, dict[str, object]] = {}
    for t in range(first_t, last_t + 1):
        row = p3_local_survivor(p, t)
        compact_counts = tuple(
            int(item["Q"]) - r - 2 for item in row["opposite_rows"]
        )
        profiles[str(t)] = {
            "delta": t - lower,
            "compact_counts": list(compact_counts),
            "one_compact_row_count": compact_counts.count(1),
        }

    before_counts = tuple(
        int(item["Q"]) - r - 2
        for item in p3_local_survivor(p, first_t - 1)["opposite_rows"]
    )
    after_counts = tuple(
        int(item["Q"]) - r - 2
        for item in p3_local_survivor(p, last_t + 1)["opposite_rows"]
    )
    proved = bool(
        local_no_go["proved"]
        and local_no_go["local_atom_profile"]
        == "one arbitrary compact plus six all-equal triangles"
        and local_no_go["required_odd_degrees"] == list(range(3, 30, 2))
        and local_no_go["required_even_degrees"] == [6, 8]
        and len(profiles) == 31
        and all(
            int(profile["one_compact_row_count"]) >= 1
            for profile in profiles.values()
        )
        and set(before_counts) == {0}
        and set(after_counts) == {2}
        and direction_count == len(before_counts) == len(after_counts)
    )
    if not proved:
        raise ArithmeticError("the p31 balanced zero-form band changed")
    return {
        "p": p,
        "r": r,
        "lower_endpoint_t": lower,
        "balanced_direction_count": direction_count,
        "balanced_t_range_excluded_for_zero_odd_six_eight_forms": [
            first_t,
            last_t,
        ],
        "balanced_delta_range": [1, 31],
        "per_t_profiles": profiles,
        "profile_before_band": list(before_counts),
        "profile_after_band": list(after_counts),
        "reused_local_atom_profile": local_no_go["local_atom_profile"],
        "exhaustive_profile_certificate_t_is_provenance_only": True,
        "zero_global_odd_six_eight_forms_compatible_with_balanced_profile": False,
        "nonzero_global_forms_ruled_out": False,
        "unbalanced_allocations_ruled_out": False,
        "F_p_common_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    arbitrary_local_no_go = p31_arbitrary_compact_six_ae_odd_six_eight_no_go()
    all_prime_odd_radon = p3_first_interior_odd_radon_centrality_certificate(31)
    full_balanced_line_exclusion = (
        p3_full_balanced_maximal_line_exclusion_certificate(31, 7)
    )
    full_balanced_two_line_exclusion = (
        p3_full_balanced_two_maximal_line_exclusion_certificate(31, 7)
    )
    boundary_cubic_exclusion = p3_boundary_cubic_unit_reduction_certificate(31)
    general_b_odd_radon = p3_bounded_compact_odd_radon_centrality_certificate(
        31, 3
    )
    balanced_odd_radon_band = p3_balanced_odd_radon_centrality_band_certificate(31)
    p31_balanced_zero_form_band = p31_balanced_zero_form_band_certificate()
    joint_even_dominance = joint_six_eight_atom_map_dominance_certificate()
    return {
        "title": "Odd blindness, all-prime branch-C centrality, and the p31 local zero-form gate",
        "status": "PROVED METHOD BARRIER AND SHARPENED OPEN REDUCTION",
        "proved": {
            "both_full_rays_pass_every_odd_moment_below_top": True,
            "degree_five_can_exclude_either_ray": False,
            "p_1_mod_4_full_ray_passes_degree_six_for_r_at_least_7": True,
            "p_3_mod_4_lower_endpoint_passes_degree_six_for_r_at_least_7": True,
            "p_1_mod_4_full_ray_passes_degree_eight_separately": True,
            "p_3_mod_4_lower_endpoint_passes_degree_eight_separately": True,
            "p31_t69_passes_all_odd_moments_and_degree_six": True,
            "p31_t69_centered_compact_zero_form_odd_six_eight_block_exists": False,
            "p31_t69_arbitrary_compact_zero_form_odd_six_eight_block_exists": False,
            "all_prime_first_interior_opposite_rows_are_central_when_odd_forms_zero": True,
            "all_prime_bounded_compact_rows_are_central_when_odd_forms_zero": True,
            "all_prime_full_balanced_maximal_line_supports_are_excluded": True,
            "all_prime_full_balanced_line_containing_supports_are_excluded": True,
            "p_at_least_31_boundary_cubic_supports_are_excluded": True,
            "balanced_branch_C_initial_band_is_central_when_odd_forms_zero": True,
            "p31_t69_through_t99_balanced_zero_odd_six_eight_forms_compatible": False,
            "universal_algebraic_relation_among_joint_degree_six_eight_channels": False,
            "same_labels_pass_degrees_six_and_eight": False,
            "all_even_moments_pass": False,
            "signed_Boolean_affine_box_nonempty": False,
            "residual_ii_closed": False,
        },
        "first_live_moment_after_these_results": (
            "joint degree-6/degree-8 labels on the full p=1 mod 4 ray and the "
            "lower p=3 mod 4 endpoint; degree 6 remains open on general "
            "interior p=3 mod 4 layers; on the balanced p=31 profiles "
            "t=69..99 only nonzero/coupled global forms can evade the "
            "exhaustive one-compact/six-AE local zero-form no-go"
        ),
        "remaining_obstruction": (
            "Simultaneous even-degree atom moments for one fixed set of labels, "
            "followed by the signed Boolean affine-box intersection"
        ),
        "duplicate_work_guards": [
            "Do not retry degree five or any odd degree on the compact atom rays.",
            "Do not promote moment-compatible local labels to one integral or Boolean graph.",
            "Do not combine separate degree-six and degree-eight scale assignments.",
            "Do not promote the centered p=3 root obstruction to arbitrary triangle labels.",
            (
                "Do not promote the p=31 arbitrary-compact zero-form row no-go "
                "to nonzero or globally coupled forms."
            ),
            (
                "Do not seek a universal polynomial identity in the seven "
                "degree-six/eight channels: both unrestricted atom maps are dominant."
            ),
        ],
        "all_prime_first_interior_odd_Radon_centrality": all_prime_odd_radon,
        "all_prime_full_balanced_maximal_line_exclusion": (
            full_balanced_line_exclusion
        ),
        "all_prime_full_balanced_two_maximal_line_exclusion": (
            full_balanced_two_line_exclusion
        ),
        "p_at_least_31_boundary_cubic_unit_exclusion": (
            boundary_cubic_exclusion
        ),
        "all_prime_bounded_compact_odd_Radon_centrality": general_b_odd_radon,
        "balanced_branch_C_odd_Radon_centrality_band": balanced_odd_radon_band,
        "p31_balanced_zero_form_band": p31_balanced_zero_form_band,
        "joint_degree_six_eight_atom_map_dominance": joint_even_dominance,
        "p31_arbitrary_compact_local_no_go": arbitrary_local_no_go,
        "L_status": "OPEN",
    }
