#!/usr/bin/env python3
"""Exact even-syndrome barrier at the equianharmonic compact threshold.

This is a characteristic-zero symbolic calculation over
``Q(q)/(q^2+q+1)``.  It classifies the component-excess bookkeeping at
``b=(2r+7)/3`` and gives two genuine four-compact/two-cycle trade families.
Their mixed seven-channel syndrome map is dominant.  No finite-field
rational solution, common lift, or residual-(ii) closure is asserted.
"""

from __future__ import annotations

from functools import lru_cache
from fractions import Fraction

from e1_gmin_m4_compact_ray_moment_gate import (
    all_equal_moment,
    compact_moment,
)
from e1_gmin_m4_conic_odd_radon import (
    P31_COMPACT_ATOMS,
)
from e1_gmin_m4_equianharmonic_component_packing import (
    ONE,
    P43_COMPACT_ATOMS,
    Q,
    Q3,
    ZERO,
)


Poly = tuple[Q3, ...]
CHANNELS = ("F60", "F61", "F62", "F80", "F81", "F82", "F83")
JACOBIAN_POINTS = (
    ("U", -1),
    ("U", -2),
    ("U", 1),
    ("U", 0),
    ("V", -3),
    ("V", -2),
    ("V", -1),
)
JACOBIAN_DETERMINANT = 4_128_623_683_475_967_290_061_619_200
JACOBIAN_FACTORS = {2: 32, 3: 26, 5: 2, 7: 1, 2161: 1}
U_AFFINE_COEFFICIENTS = (
    Fraction(84461, 12096),
    Fraction(-187801, 6048),
    Fraction(619, 63),
    Fraction(-2509, 5184),
    Fraction(95099, 36288),
    Fraction(-3557, 1134),
    Fraction(1),
)
U_AFFINE_CONSTANT = Fraction(10176, 7)


def _trim(poly: Poly) -> Poly:
    values = list(poly)
    while len(values) > 1 and not values[-1]:
        values.pop()
    return tuple(values)


def _constant(value: int | Fraction | Q3) -> Poly:
    return (Q3.coerce(value),)


def _add(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return _trim(
        tuple(
            (left[index] if index < len(left) else ZERO)
            + (right[index] if index < len(right) else ZERO)
            for index in range(size)
        )
    )


def _neg(poly: Poly) -> Poly:
    return tuple(-value for value in poly)


def _sub(left: Poly, right: Poly) -> Poly:
    return _add(left, _neg(right))


def _mul(left: Poly, right: Poly) -> Poly:
    result = [ZERO] * (len(left) + len(right) - 1)
    for first, a in enumerate(left):
        for second, b in enumerate(right):
            result[first + second] = result[first + second] + a * b
    return _trim(tuple(result))


def _scale(poly: Poly, scalar: int | Fraction | Q3) -> Poly:
    scalar = Q3.coerce(scalar)
    return _trim(tuple(scalar * value for value in poly))


def _power(poly: Poly, exponent: int) -> Poly:
    result = _constant(1)
    base = poly
    while exponent:
        if exponent & 1:
            result = _mul(result, base)
        base = _mul(base, base)
        exponent >>= 1
    return result


def _derivative(poly: Poly) -> Poly:
    return _trim(
        tuple(index * poly[index] for index in range(1, len(poly)))
        or (ZERO,)
    )


def _evaluate(poly: Poly, value: int | Q3) -> Q3:
    value = Q3.coerce(value)
    result = ZERO
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def _phi(poly: Poly) -> Poly:
    return _add(_scale(poly, Q), _constant(1 - Q))


def _edge_moment(first: Poly, second: Poly, degree: int, channel: int) -> Poly:
    return _mul(
        _mul(
            _power(_sub(first, second), 2),
            _power(_mul(first, second), channel),
        ),
        _power(_add(first, second), degree - 2 - 2 * channel),
    )


def _compact(
    first: Poly, second: Poly, distinguished: Poly, degree: int, channel: int
) -> Poly:
    return _sub(
        _sub(
            _edge_moment(first, second, degree, channel),
            _edge_moment(first, distinguished, degree, channel),
        ),
        _edge_moment(second, distinguished, degree, channel),
    )


def _all_equal(first: Poly, second: Poly, third: Poly, degree: int, channel: int) -> Poly:
    return _add(
        _add(
            _edge_moment(first, second, degree, channel),
            _edge_moment(first, third, degree, channel),
        ),
        _edge_moment(second, third, degree, channel),
    )


def _trade_atoms_and_cycles(
    family: str,
) -> tuple[
    tuple[tuple[Poly, Poly, Poly], ...],
    tuple[tuple[Poly, Poly, Poly], ...],
]:
    """Return the symbolic compact atoms and two replaced AE cycles."""
    if family not in {"U", "V"}:
        raise ValueError("family must be U or V")
    x = (ZERO, ONE)
    a2 = x
    a0 = _phi(a2)
    a1 = _phi(a0)
    b2 = _neg(a2)
    b0 = _phi(b2)
    b1 = _phi(b0)
    if family == "U":
        atoms = (
            (b2, _neg(a1), _neg(a0)),
            (a1, a2, _neg(b1)),
            (a1, _neg(b1), _neg(b0)),
            (b0, b2, _neg(a1)),
        )
    else:
        atoms = (
            (a0, a2, _neg(b0)),
            (_neg(a1), b0, _neg(a0)),
            (_neg(b0), a1, _neg(b1)),
            (b2, b1, _neg(a1)),
        )
    return atoms, ((a0, a1, a2), (b0, b1, b2))


def _trade_orbit_chain_remainder(family: str) -> dict[object, int]:
    """Subtract the two AE chains from a U/V chain in a fixed orbit basis."""
    atoms, cycles = _trade_atoms_and_cycles(family)
    chain: dict[object, int] = {}

    def add_edge(first: Poly, second: Poly, coefficient: int) -> None:
        edge = tuple(sorted((first, second)))
        antipode = tuple(sorted((_neg(first), _neg(second))))
        if edge == antipode:
            raise ArithmeticError("generic U/V edge became self-antipodal")
        representative, orientation = (
            (edge, 1) if edge < antipode else (antipode, -1)
        )
        chain[representative] = chain.get(representative, 0) + (
            orientation * coefficient
        )

    for first, second, distinguished in atoms:
        add_edge(first, second, 1)
        add_edge(first, distinguished, -1)
        add_edge(second, distinguished, -1)
    for first, second, third in cycles:
        add_edge(first, second, -1)
        add_edge(first, third, -1)
        add_edge(second, third, -1)
    return {edge: value for edge, value in chain.items() if value}


@lru_cache(maxsize=2)
def trade_deviation_polynomials(family: str) -> tuple[Poly, ...]:
    """Return ``2^d`` times the seven even deviations for a U/V trade."""
    atoms, cycles = _trade_atoms_and_cycles(family)
    result = []
    for degree in (6, 8):
        for channel in range(degree // 2):
            compact_sum = _constant(0)
            for atom in atoms:
                compact_sum = _add(
                    compact_sum, _compact(*atom, degree, channel)
                )
            canonical = _add(
                *(
                    _all_equal(*cycle, degree, channel)
                    for cycle in cycles
                )
            )
            result.append(_sub(compact_sum, canonical))
    return tuple(result)


def _determinant(matrix: list[list[Q3]]) -> Q3:
    values = [row[:] for row in matrix]
    result = ONE
    for column in range(len(values)):
        pivot = next(
            (row for row in range(column, len(values)) if values[row][column]),
            None,
        )
        if pivot is None:
            return ZERO
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            result = -result
        pivot_value = values[column][column]
        result = result * pivot_value
        for row in range(column + 1, len(values)):
            factor = values[row][column] / pivot_value
            for index in range(column + 1, len(values)):
                values[row][index] = (
                    values[row][index] - factor * values[column][index]
                )
    return result


def _q3_mod_p(value: Q3, p: int, q_value: int) -> int:
    def reduce_fraction(number: Fraction) -> int:
        return number.numerator * pow(number.denominator, -1, p) % p

    return (
        reduce_fraction(value.a) + q_value * reduce_fraction(value.b)
    ) % p


def _symbolic_deviation_mod_p(
    family: str, p: int, q_value: int, x_value: int
) -> dict[int, tuple[int, ...]]:
    polynomials = trade_deviation_polynomials(family)
    result = {}
    offset = 0
    for degree in (6, 8):
        scale = pow(pow(2, degree, p), -1, p)
        result[degree] = tuple(
            _q3_mod_p(
                _evaluate(polynomials[offset + channel], x_value),
                p,
                q_value,
            )
            * scale
            % p
            for channel in range(degree // 2)
        )
        offset += degree // 2
    return result


def _actual_trade_deviation(
    p: int,
    k: int,
    x_value: int,
    compact_atoms: tuple[tuple[tuple[int, int, int], int], ...],
) -> dict[int, tuple[int, ...]]:
    q_value = (1 - k) * pow(1 + k, -1, p) % p

    def phi(value: int) -> int:
        return (q_value * value + 1 - q_value) % p

    a2 = x_value % p
    a0 = phi(a2)
    a1 = phi(a0)
    b2 = -a2 % p
    b0 = phi(b2)
    b1 = phi(b0)
    half = pow(2, -1, p)
    cycles = (
        tuple(value * half % p for value in (a0, a1, a2)),
        tuple(value * half % p for value in (b0, b1, b2)),
    )
    result = {}
    for degree in (6, 8):
        values = []
        for channel in range(degree // 2):
            compact_value = 0
            for triple, distinguished in compact_atoms:
                positive = tuple(
                    value for value in triple if value != distinguished
                )
                compact_value += compact_moment(
                    p,
                    positive[0],
                    positive[1],
                    distinguished,
                    degree,
                    channel,
                )
            canonical = sum(
                all_equal_moment(p, *cycle, degree, channel)
                for cycle in cycles
            )
            values.append((compact_value - canonical) % p)
        result[degree] = tuple(values)
    return result


def threshold_excess_assembly_certificate() -> dict[str, object]:
    """Classify all positive/negative excess totals when global excess is 1."""
    def tuples_with_excess(target: int) -> list[tuple[int, int, int, int]]:
        rows = []
        for compact_count in range(7):
            for ae_count in range(2):
                for cycle_rank in range(2):
                    for caps in range(4):
                        if compact_count + ae_count == 0:
                            continue
                        excess = (
                            4
                            - compact_count
                            - 4 * ae_count
                            - 4 * cycle_rank
                            - 2 * caps
                        )
                        if excess == target:
                            rows.append(
                                (compact_count, ae_count, cycle_rank, caps)
                            )
        return rows

    zero = tuples_with_excess(0)
    minus_one = tuples_with_excess(-1)
    minus_two = tuples_with_excess(-2)
    proved = bool(
        zero == [(0, 1, 0, 0), (2, 0, 0, 1), (4, 0, 0, 0)]
        and minus_one
        == [
            (1, 0, 0, 2),
            (1, 0, 1, 0),
            (1, 1, 0, 0),
            (3, 0, 0, 1),
            (5, 0, 0, 0),
        ]
        and minus_two
        == [
            (0, 1, 0, 1),
            (2, 0, 0, 2),
            (2, 0, 1, 0),
            (2, 1, 0, 0),
            (4, 0, 0, 1),
            (6, 0, 0, 0),
        ]
    )
    return {
        "threshold": "b=2L-1=(2r+7)/3",
        "global_deficit": "Delta=L-1",
        "global_component_excess": "b-2*Delta=1",
        "positive_mass_negative_mass_possibilities": [[1, 0], [2, 1], [3, 2]],
        "positive_block_assemblies": {
            "1": [["cap"], ["F"]],
            "2": [["HH"], ["F", "cap"]],
            "3": [["HH", "cap"]],
        },
        "zero_excess_tuples_K_AE_cycle_rank_caps": zero,
        "minus_one_tuples_K_AE_cycle_rank_caps": minus_one,
        "minus_two_tuples_K_AE_cycle_rank_caps": minus_two,
        "negative_mass_one": "one minus-one component",
        "negative_mass_two": (
            "one minus-two component or two minus-one components"
        ),
        "proved": proved,
    }


def equianharmonic_threshold_even_barrier_certificate() -> dict[str, object]:
    """Return exact U/V formulas, dominance, and witness replays."""
    assembly = threshold_excess_assembly_certificate()
    u = trade_deviation_polynomials("U")
    v = trade_deviation_polynomials("V")

    u_affine = _constant(0)
    for coefficient, polynomial in zip(U_AFFINE_COEFFICIENTS, u):
        u_affine = _add(u_affine, _scale(polynomial, coefficient))

    jacobian = [
        [
            _evaluate(
                _derivative(
                    u_polynomial if family == "U" else v_polynomial
                ),
                point,
            )
            for family, point in JACOBIAN_POINTS
        ]
        for u_polynomial, v_polynomial in zip(u, v)
    ]
    determinant = _determinant(jacobian)
    factor_product = 1
    for prime, exponent in JACOBIAN_FACTORS.items():
        factor_product *= prime**exponent

    witness_cases = {
        "p31_U_x18": (
            "U",
            31,
            11,
            18,
            tuple(P31_COMPACT_ATOMS[index] for index in (2, 3, 6, 5)),
        ),
        "p43_U_x38": (
            "U",
            43,
            13,
            38,
            tuple(P43_COMPACT_ATOMS[index] for index in (0, 5, 6, 7)),
        ),
        "p43_V_x7": (
            "V",
            43,
            13,
            7,
            tuple(P43_COMPACT_ATOMS[index] for index in (1, 3, 2, 4)),
        ),
    }
    replays = {}
    for name, (family, p, k, x_value, atoms) in witness_cases.items():
        q_value = (1 - k) * pow(1 + k, -1, p) % p
        symbolic = _symbolic_deviation_mod_p(family, p, q_value, x_value)
        actual = _actual_trade_deviation(p, k, x_value, atoms)
        replays[name] = {
            "family": family,
            "p": p,
            "k": k,
            "q": q_value,
            "x": x_value,
            "degree_six_deviation": list(actual[6]),
            "degree_eight_deviation": list(actual[8]),
            "symbolic_matches_atom_witness": symbolic == actual,
        }

    proved = bool(
        assembly["proved"]
        and not _trade_orbit_chain_remainder("U")
        and not _trade_orbit_chain_remainder("V")
        and u_affine == _constant(U_AFFINE_CONSTANT)
        and determinant == Q3(JACOBIAN_DETERMINANT)
        and factor_product == JACOBIAN_DETERMINANT
        and all(
            row["symbolic_matches_atom_witness"] for row in replays.values()
        )
    )
    if not proved:
        raise ArithmeticError("the threshold even-syndrome barrier changed")
    return {
        "coefficient_field": "Q(q)/(q^2+q+1)",
        "normalization": "X=2x",
        "channels": list(CHANNELS),
        "assembly": assembly,
        "trade_definition": (
            "A=(Phi(x),Phi^2(x),x), B=(Phi(-x),Phi^2(-x),-x); "
            "each U/V block is four compact atoms replacing AE(A)+AE(B)"
        ),
        "U_odd_edge_orbit_chain_exact": True,
        "V_odd_edge_orbit_chain_exact": True,
        "U_deviation_coefficients_low_to_high": [
            [coefficient.text() for coefficient in polynomial]
            for polynomial in u
        ],
        "V_deviation_coefficients_low_to_high": [
            [coefficient.text() for coefficient in polynomial]
            for polynomial in v
        ],
        "U_affine_syndrome_invariant": {
            "coefficients_in_channel_order": [
                str(value) for value in U_AFFINE_COEFFICIENTS
            ],
            "constant_per_trade": str(U_AFFINE_CONSTANT),
            "proved": True,
        },
        "mixed_jacobian": {
            "columns_family_and_x": [list(row) for row in JACOBIAN_POINTS],
            "determinant": JACOBIAN_DETERMINANT,
            "factorization": JACOBIAN_FACTORS,
            "nonzero_characteristic_zero": True,
            "dominant_outside_displayed_characteristics": True,
            "exceptional_characteristics": sorted(JACOBIAN_FACTORS),
        },
        "witness_replays": replays,
        "consequence": (
            "component excess alone supplies no common affine identity in "
            "the seven degree-six/eight channels once U and V trades mix"
        ),
        "finite_field_rational_zero_syndrome_trade_matching_constructed": False,
        "uniform_zero_degree_six_eight_exclusion_proved": False,
        "common_global_form_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(equianharmonic_threshold_even_barrier_certificate(), indent=2))
