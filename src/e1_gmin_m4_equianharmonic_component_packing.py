#!/usr/bin/env python3
"""Exact symbolic certificate for equianharmonic cancellation components.

Work over Q(q), q^2+q+1=0, in the normalized X=2x coordinate.  This
module classifies precisely the cancellation components with positive
``K-2*deficit`` in a constant tangent-conic atom fiber and proves the
resulting exact-fiber compact-count threshold.  It deliberately does not
assert a global common-form lift or residual-(ii).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import gcd


@dataclass(frozen=True, order=True)
class Q3:
    """Element ``a+b*q`` of Q[q]/(q^2+q+1)."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __init__(self, a: int | Fraction = 0, b: int | Fraction = 0):
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))

    @staticmethod
    def coerce(value: int | Fraction | Q3) -> Q3:
        return value if isinstance(value, Q3) else Q3(value)

    def __add__(self, other: int | Fraction | Q3) -> Q3:
        other = self.coerce(other)
        return Q3(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> Q3:
        return Q3(-self.a, -self.b)

    def __sub__(self, other: int | Fraction | Q3) -> Q3:
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Fraction | Q3) -> Q3:
        return self.coerce(other) - self

    def __mul__(self, other: int | Fraction | Q3) -> Q3:
        other = self.coerce(other)
        # q^2=-q-1.
        return Q3(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    __rmul__ = __mul__

    def inverse(self) -> Q3:
        norm = self.norm()
        if not norm:
            raise ZeroDivisionError("zero has no inverse")
        return Q3((self.a - self.b) / norm, -self.b / norm)

    def __truediv__(self, other: int | Fraction | Q3) -> Q3:
        return self * self.coerce(other).inverse()

    def __rtruediv__(self, other: int | Fraction | Q3) -> Q3:
        return self.coerce(other) * self.inverse()

    def norm(self) -> Fraction:
        return self.a * self.a - self.a * self.b + self.b * self.b

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    def text(self) -> str:
        if not self.b:
            return str(self.a)
        if not self.a:
            return "q" if self.b == 1 else "-q" if self.b == -1 else f"{self.b}q"
        sign = "+" if self.b > 0 else ""
        coefficient = "" if self.b == 1 else "-" if self.b == -1 else str(self.b)
        return f"{self.a}{sign}{coefficient}q"


ZERO = Q3()
ONE = Q3(1)
Q = Q3(0, 1)
Q2 = Q * Q
POWERS = (Q, Q2)

P43_AE_ATOMS = (
    (2, 31, 33),
    (3, 26, 37),
    (5, 6, 12),
    (9, 27, 30),
    (10, 20, 36),
    (11, 13, 42),
    (14, 17, 35),
    (16, 21, 29),
    (32, 38, 39),
)

P43_COMPACT_ATOMS = (
    ((0, 24, 39), 0),
    ((1, 25, 36), 36),
    ((2, 36, 40), 2),
    ((3, 7, 42), 42),
    ((3, 18, 41), 3),
    ((4, 9, 19), 9),
    ((4, 9, 35), 35),
    ((8, 24, 39), 39),
    ((15, 20, 23), 20),
)

Affine2 = tuple[Q3, Q3, Q3]
Affine1 = tuple[Q3, Q3]
Edge = frozenset[Q3]
AffineEdge = frozenset[Affine1]


def _phi(value: Q3, multiplier: Q3) -> Q3:
    return multiplier * value + 1 - multiplier


def _psi(value: Q3, multiplier: Q3) -> Q3:
    return multiplier * value + multiplier - 1


def _edge(a: Q3, b: Q3) -> Edge:
    if a == b:
        raise ValueError("loop")
    return frozenset((a, b))


def _negate_edge(edge: Edge) -> Edge:
    return frozenset(-value for value in edge)


def _affine(a: Q3 = ZERO, b: Q3 = ZERO, c: Q3 = ZERO) -> Affine2:
    return a, b, c


def _affine_neg(form: Affine2) -> Affine2:
    return tuple(-value for value in form)  # type: ignore[return-value]


def _affine_sub(left: Affine2, right: Affine2) -> Affine2:
    return tuple(left[i] - right[i] for i in range(3))  # type: ignore[return-value]


def _affine_scale(scale: Q3, form: Affine2) -> Affine2:
    return tuple(scale * value for value in form)  # type: ignore[return-value]


def _affine_phi(form: Affine2, multiplier: Q3) -> Affine2:
    result = _affine_scale(multiplier, form)
    return result[0], result[1], result[2] + 1 - multiplier


def _affine_psi(form: Affine2, multiplier: Q3) -> Affine2:
    result = _affine_scale(multiplier, form)
    return result[0], result[1], result[2] + multiplier - 1


def _evaluate(form: Affine2, first: Q3, second: Q3) -> Q3:
    return form[0] * first + form[1] * second + form[2]


def _high_defect_condition(
    first: Affine2,
    second: Affine2,
    phi_multiplier: Q3,
    psi_multiplier: Q3,
    swap: int,
) -> Affine2:
    """Linear condition that ``{first,second}`` is a high PN defect."""
    if swap:
        first, second = second, first
    ratio = psi_multiplier / phi_multiplier
    constant = 2 * psi_multiplier - ratio - 1
    result = _affine_sub(second, _affine_scale(ratio, first))
    return result[0], result[1], result[2] - constant


def _solve_rank_two(first: Affine2, second: Affine2) -> tuple[Q3, Q3] | None:
    a, b, c = first
    d, e, f = second
    determinant = a * e - b * d
    if not determinant:
        return None
    return (b * f - c * e) / determinant, (c * d - a * f) / determinant


def _system_class(first: Affine2, second: Affine2) -> str:
    a, b, c = first
    d, e, f = second
    if a * e - b * d:
        return "rank_two"
    minors = (a * f - d * c, b * f - e * c)
    if not a and not b and not d and not e:
        return "rank_zero" if not c and not f else "inconsistent"
    return "rank_one" if not any(minors) else "inconsistent"


def _rank_one_parameterization(first: Affine2, second: Affine2) -> tuple[Affine1, Affine1]:
    if _system_class(first, second) != "rank_one":
        raise ValueError("not a consistent rank-one system")
    equation = first if first[0] or first[1] else second
    a, b, c = equation
    if b:
        # first variable is the parameter t.
        return (ONE, ZERO), (-a / b, -c / b)
    return (-c / a, ZERO), (ONE, ZERO)


def _substitute(form: Affine2, first: Affine1, second: Affine1) -> Affine1:
    return (
        form[0] * first[0] + form[1] * second[0],
        form[0] * first[1] + form[1] * second[1] + form[2],
    )


def _high_data(
    defect_first: Affine2,
    defect_second: Affine2,
    phi_multiplier: Q3,
    psi_multiplier: Q3,
    swap: int,
    first: Q3,
    second: Q3,
) -> tuple[tuple[Q3, Q3, Q3], frozenset[Edge]]:
    atom = _high_atom(
        defect_first,
        defect_second,
        phi_multiplier,
        psi_multiplier,
        swap,
        first,
        second,
    )
    root, endpoint_b, endpoint_c = atom
    support = frozenset(
        (
            _edge(root, endpoint_b),
            _edge(-root, _phi(-root, psi_multiplier)),
        )
    )
    return atom, support


def _high_atom(
    defect_first: Affine2,
    defect_second: Affine2,
    phi_multiplier: Q3,
    psi_multiplier: Q3,
    swap: int,
    first: Q3,
    second: Q3,
) -> tuple[Q3, Q3, Q3]:
    endpoint_b = _evaluate(defect_first, first, second)
    endpoint_c = _evaluate(defect_second, first, second)
    if swap:
        endpoint_b, endpoint_c = endpoint_c, endpoint_b
    root = (endpoint_b - 1 + phi_multiplier) / phi_multiplier
    if _phi(root, phi_multiplier) != endpoint_b:
        raise AssertionError("bad Phi inversion")
    if _psi(root, psi_multiplier) != endpoint_c:
        raise AssertionError("bad Psi inversion")
    return root, endpoint_b, endpoint_c


def _high_target_affine_edges(
    defect_first: Affine2,
    defect_second: Affine2,
    phi_multiplier: Q3,
    psi_multiplier: Q3,
    swap: int,
    first: Affine1,
    second: Affine1,
) -> tuple[AffineEdge, AffineEdge]:
    if swap:
        defect_first, defect_second = defect_second, defect_first
    endpoint_b = _substitute(defect_first, first, second)
    root_form = _affine_scale(1 / phi_multiplier, defect_first)
    root_form = root_form[0], root_form[1], root_form[2] + 1 - 1 / phi_multiplier
    root = _substitute(root_form, first, second)
    negative_root = (-root[0], -root[1])
    phi_negative_root = (
        psi_multiplier * negative_root[0],
        psi_multiplier * negative_root[1] + 1 - psi_multiplier,
    )
    return (
        frozenset((root, endpoint_b)),
        frozenset((negative_root, phi_negative_root)),
    )


def _edge_kind(edge: Edge) -> str:
    if edge == _negate_edge(edge):
        return "self"
    x, y = tuple(edge)
    in_positive = any(
        y == _phi(x, multiplier) or x == _phi(y, multiplier)
        for multiplier in POWERS
    )
    in_negative = any(
        y == _psi(x, multiplier) or x == _psi(y, multiplier)
        for multiplier in POWERS
    )
    if in_positive and in_negative:
        raise AssertionError("positive and negative target graphs intersect")
    return "positive" if in_positive else "negative" if in_negative else "outside"


def _atom_alignment_values(atom: tuple[Q3, Q3, Q3]) -> tuple[int, int, int]:
    a, b, distinguished = atom
    occurrences = (
        (_edge(a, b), 1),
        (_edge(a, distinguished), -1),
        (_edge(b, distinguished), -1),
    )
    values = []
    for edge, coefficient in occurrences:
        kind = _edge_kind(edge)
        aligned_kind = "positive" if coefficient == 1 else "negative"
        reverse_kind = "negative" if coefficient == 1 else "positive"
        values.append(1 if kind == aligned_kind else -1 if kind == reverse_kind else 0)
    return tuple(values)  # type: ignore[return-value]


def _normalized_target_edge(edge: Edge) -> Edge:
    kind = _edge_kind(edge)
    if kind == "positive":
        return edge
    if kind == "negative":
        return _negate_edge(edge)
    raise ValueError("edge is not a target coordinate")


def _normalized_aligned_support(atom: tuple[Q3, Q3, Q3]) -> frozenset[Edge]:
    a, b, distinguished = atom
    occurrences = (
        (_edge(a, b), 1),
        (_edge(a, distinguished), -1),
        (_edge(b, distinguished), -1),
    )
    result = []
    for edge, coefficient in occurrences:
        kind = _edge_kind(edge)
        if (coefficient == 1 and kind == "positive") or (
            coefficient == -1 and kind == "negative"
        ):
            result.append(_normalized_target_edge(edge))
    return frozenset(result)


def _o2_free_supports() -> dict[str, object]:
    variable_u = _affine(ONE)
    branch_counts = {"rank_two": 0, "rank_one": 0, "inconsistent": 0}
    valid_assignments = 0
    rank_one_collisions = 0
    supports: set[frozenset[Edge]] = set()
    determinant_primes: set[int] = set()
    augmented_primes: set[int] = set()
    relation_change_primes: set[int] = set()

    for center_type, center_multiplier in product(("P", "N"), POWERS):
        if center_type == "P":
            variable_v = _affine_phi(variable_u, center_multiplier)
            variable_w = _affine(ZERO, ONE)
            desired = (
                (_affine_neg(variable_u), _affine_neg(variable_w)),
                (_affine_neg(variable_v), _affine_neg(variable_w)),
            )
            center_atom = lambda u, w: (u, _phi(u, center_multiplier), w)
            center_affine_edge = frozenset((variable_u, variable_v))
        else:
            variable_v = _affine(ZERO, ONE)
            variable_w = _affine_psi(variable_u, center_multiplier)
            desired = (
                (variable_u, variable_v),
                (_affine_neg(variable_v), _affine_neg(variable_w)),
            )
            center_atom = lambda u, v: (u, v, _psi(u, center_multiplier))
            center_affine_edge = frozenset(
                (_affine_neg(variable_u), _affine_neg(variable_w))
            )

        choices = product(POWERS, POWERS, (0, 1), POWERS, POWERS, (0, 1))
        for phi_one, psi_one, swap_one, phi_two, psi_two, swap_two in choices:
            equations = (
                _high_defect_condition(
                    *desired[0], phi_one, psi_one, swap_one
                ),
                _high_defect_condition(
                    *desired[1], phi_two, psi_two, swap_two
                ),
            )
            system_class = _system_class(*equations)
            if system_class == "rank_two":
                branch_counts["rank_two"] += 1
                determinant = (
                    equations[0][0] * equations[1][1]
                    - equations[0][1] * equations[1][0]
                )
                determinant_primes.update(_element_prime_factors(determinant))
                solution = _solve_rank_two(*equations)
                if solution is None:
                    raise AssertionError("rank-two system was not solved")
                first, second = solution
                center = center_atom(first, second)
                high_one_atom = _high_atom(
                    *desired[0], phi_one, psi_one, swap_one, first, second
                )
                high_two_atom = _high_atom(
                    *desired[1], phi_two, psi_two, swap_two, first, second
                )
                for atom in (center, high_one_atom, high_two_atom):
                    relation_change_primes.update(
                        _atom_relation_change_primes(atom)
                    )
                relation_change_primes.update(
                    _label_relation_change_primes(
                        center + high_one_atom + high_two_atom
                    )
                )
                try:
                    high_one, support_one = _high_data(
                        *desired[0], phi_one, psi_one, swap_one, first, second
                    )
                    high_two, support_two = _high_data(
                        *desired[1], phi_two, psi_two, swap_two, first, second
                    )
                except ValueError:
                    # A loop or repeated label is a degenerate atom, not a
                    # new free component.
                    continue
                if any(len(set(atom)) < 3 for atom in (center, high_one, high_two)):
                    continue
                if _atom_alignment_values(center).count(1) != 1:
                    continue
                if _atom_alignment_values(center).count(0) != 2:
                    continue
                if _atom_alignment_values(high_one).count(1) != 2:
                    continue
                if _atom_alignment_values(high_one).count(0) != 1:
                    continue
                if _atom_alignment_values(high_two).count(1) != 2:
                    continue
                if _atom_alignment_values(high_two).count(0) != 1:
                    continue
                support = frozenset(
                    set(_normalized_aligned_support(center))
                    | set(support_one)
                    | set(support_two)
                )
                if len(support) != 5:
                    continue
                valid_assignments += 1
                supports.add(support)
                continue

            if system_class == "rank_one":
                branch_counts["rank_one"] += 1
                first, second = _rank_one_parameterization(*equations)
                substituted_center = frozenset(
                    _substitute(endpoint, first, second)
                    for endpoint in center_affine_edge
                )
                affine_support = [substituted_center]
                affine_support.extend(
                    _high_target_affine_edges(
                        *desired[0], phi_one, psi_one, swap_one, first, second
                    )
                )
                affine_support.extend(
                    _high_target_affine_edges(
                        *desired[1], phi_two, psi_two, swap_two, first, second
                    )
                )
                if len(set(affine_support)) < 5:
                    rank_one_collisions += 1
                continue

            branch_counts["inconsistent"] += 1
            augmented = _augmented_witnesses(*equations)
            for witness in augmented:
                augmented_primes.update(_element_prime_factors(witness))

    return {
        "branch_counts": branch_counts,
        "valid_rank_two_assignments": valid_assignments,
        "rank_one_collision_count": rank_one_collisions,
        "supports": frozenset(supports),
        "rank_two_determinant_primes": sorted(determinant_primes),
        "inconsistent_augmented_minor_primes": sorted(augmented_primes),
        "rank_two_filter_relation_change_primes": sorted(
            relation_change_primes
        ),
    }


def _augmented_witnesses(first: Affine2, second: Affine2) -> tuple[Q3, ...]:
    a, b, c = first
    d, e, f = second
    if not a and not b and not d and not e:
        return tuple(value for value in (c, f) if value)
    return tuple(value for value in (a * f - d * c, b * f - e * c) if value)


def _prime_factors(value: int) -> set[int]:
    value = abs(value)
    factors: set[int] = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.add(value)
    return factors


def _element_prime_factors(value: Q3) -> set[int]:
    if not value:
        return set()
    norm = value.norm()
    return _prime_factors(norm.numerator) | _prime_factors(norm.denominator)


def _atom_relation_change_primes(atom: tuple[Q3, Q3, Q3]) -> set[int]:
    """Characteristics where a branch-filter equality could change.

    This includes label collisions, self-antipodal edges, and every ordered
    Phi/Psi incidence test used by ``_atom_alignment_values``.  Exact-zero
    identities remain identities in every characteristic and are skipped.
    """
    result: set[int] = set()
    for value in atom:
        result.update(_prime_factors(value.a.denominator))
        result.update(_prime_factors(value.b.denominator))
    result.update(_label_relation_change_primes(atom))
    return result


def _label_relation_change_primes(labels: tuple[Q3, ...]) -> set[int]:
    """Audit every equality used when filtering a whole candidate branch."""
    result: set[int] = set()
    for first, second in combinations(labels, 2):
        residuals = [first - second, first + second]
        for multiplier in POWERS:
            residuals.extend(
                (
                    second - _phi(first, multiplier),
                    first - _phi(second, multiplier),
                    second - _psi(first, multiplier),
                    first - _psi(second, multiplier),
                )
            )
        for residual in residuals:
            if residual:
                result.update(_element_prime_factors(residual))
    return result


def _hh_free_support() -> dict[str, object]:
    branch_counts = {"rank_two": 0, "rank_one": 0, "inconsistent": 0}
    valid_assignments = 0
    rank_one_collisions = 0
    supports: set[frozenset[Edge]] = set()
    representative_atoms: tuple[tuple[Q3, Q3, Q3], tuple[Q3, Q3, Q3]] | None = None
    determinant_primes: set[int] = set()
    augmented_primes: set[int] = set()
    relation_change_primes: set[int] = set()
    first_variable = _affine(ONE)
    second_variable = _affine(ZERO, ONE)

    for phi_one, psi_one, phi_two, psi_two, swap in product(
        POWERS, POWERS, POWERS, POWERS, (0, 1)
    ):
        first_defect = (
            _affine_phi(first_variable, phi_one),
            _affine_psi(first_variable, psi_one),
        )
        second_defect = (
            _affine_phi(second_variable, phi_two),
            _affine_psi(second_variable, psi_two),
        )
        desired = (
            _affine_neg(first_defect[swap]),
            _affine_neg(first_defect[1 - swap]),
        )
        equations = (
            _affine_sub(second_defect[0], desired[0]),
            _affine_sub(second_defect[1], desired[1]),
        )
        system_class = _system_class(*equations)
        if system_class == "rank_two":
            branch_counts["rank_two"] += 1
            determinant = (
                equations[0][0] * equations[1][1]
                - equations[0][1] * equations[1][0]
            )
            determinant_primes.update(_element_prime_factors(determinant))
            solution = _solve_rank_two(*equations)
            if solution is None:
                raise AssertionError("rank-two HH system was not solved")
            first, second = solution
            atom_one = (
                first,
                _phi(first, phi_one),
                _psi(first, psi_one),
            )
            atom_two = (
                second,
                _phi(second, phi_two),
                _psi(second, psi_two),
            )
            for atom in (atom_one, atom_two):
                relation_change_primes.update(
                    _atom_relation_change_primes(atom)
                )
            relation_change_primes.update(
                _label_relation_change_primes(atom_one + atom_two)
            )
            if any(len(set(atom)) < 3 for atom in (atom_one, atom_two)):
                continue
            if any(_atom_alignment_values(atom).count(1) != 2 for atom in (atom_one, atom_two)):
                continue
            support = frozenset(
                set(_normalized_aligned_support(atom_one))
                | set(_normalized_aligned_support(atom_two))
            )
            if len(support) != 4:
                continue
            valid_assignments += 1
            supports.add(support)
            representative_atoms = representative_atoms or (atom_one, atom_two)
            continue

        if system_class == "rank_one":
            branch_counts["rank_one"] += 1
            first, second = _rank_one_parameterization(*equations)
            first_edges = _high_target_affine_edges(
                *first_defect, phi_one, psi_one, 0, first, second
            )
            second_edges = _high_target_affine_edges(
                *second_defect, phi_two, psi_two, 0, first, second
            )
            if len(set(first_edges + second_edges)) < 4:
                rank_one_collisions += 1
            continue

        branch_counts["inconsistent"] += 1
        for witness in _augmented_witnesses(*equations):
            augmented_primes.update(_element_prime_factors(witness))

    if len(supports) != 1 or representative_atoms is None:
        raise AssertionError("HH classification did not collapse to one support")
    return {
        "branch_counts": branch_counts,
        "valid_rank_two_assignments": valid_assignments,
        "rank_one_collision_count": rank_one_collisions,
        "support": next(iter(supports)),
        "representative_atoms": representative_atoms,
        "rank_two_determinant_primes": sorted(determinant_primes),
        "inconsistent_augmented_minor_primes": sorted(augmented_primes),
        "rank_two_filter_relation_change_primes": sorted(
            relation_change_primes
        ),
    }


def _phi_triangle_edges(edge: Edge) -> frozenset[Edge]:
    first = next(iter(edge))
    vertices = (first, _phi(first, Q), _phi(_phi(first, Q), Q))
    if not set(edge) <= set(vertices):
        raise AssertionError("edge is not in its Phi triangle")
    return frozenset(_edge(vertices[i], vertices[j]) for i, j in combinations(range(3), 2))


def _r1_flip_supports(hh_support: frozenset[Edge]) -> frozenset[frozenset[Edge]]:
    result = set()
    for reverse_edge in hh_support:
        triangle = _phi_triangle_edges(reverse_edge)
        support = frozenset(
            (set(hh_support) - {reverse_edge})
            | (set(triangle) - {reverse_edge})
        )
        if len(support) != 5:
            raise AssertionError("R1 flip did not have five boundary coordinates")
        result.add(support)
    return frozenset(result)


def _ordered_antipodal_edge_lifts(edge: Edge) -> tuple[tuple[Q3, Q3], ...]:
    """Return both orders of both physical representatives of an edge orbit."""
    first, second = tuple(edge)
    return tuple(
        dict.fromkeys(
            (
                (first, second),
                (second, first),
                (-first, -second),
                (-second, -first),
            )
        )
    )


def _r1_exhaustive_supports(hh_support: frozenset[Edge]) -> dict[str, object]:
    """Enumerate all three placements of the unique reverse occurrence.

    An R1 low compact has one reverse target occurrence and two aligned
    occurrences.  The reverse occurrence may formally be the positive edge
    or either negative edge.  This audit tries both orders and both physical
    representatives of every HH target orbit, so the conclusion does not
    assume in advance that the positive edge is the reverse one.
    """
    expected_values = {
        "positive": (-1, 1, 1),
        "negative_ac": (1, -1, 1),
        "negative_bc": (1, 1, -1),
    }
    attempts = {name: 0 for name in expected_values}
    consistent = {name: 0 for name in expected_values}
    valid = {name: 0 for name in expected_values}
    inconsistent_primes: set[int] = set()
    relation_change_primes: set[int] = set()
    supports: set[frozenset[Edge]] = set()
    atoms: set[tuple[Q3, Q3, Q3]] = set()

    for reverse_edge in hh_support:
        for position in expected_values:
            for first, second in _ordered_antipodal_edge_lifts(reverse_edge):
                for multiplier_one, multiplier_two in product(POWERS, repeat=2):
                    attempts[position] += 1
                    if position == "positive":
                        a, b = first, second
                        c = _psi(a, multiplier_one)
                        residual = c - _psi(b, multiplier_two)
                        reverse_occurrence = lambda: _edge(a, b)
                    elif position == "negative_ac":
                        a, c = first, second
                        b = _phi(a, multiplier_one)
                        residual = c - _psi(b, multiplier_two)
                        reverse_occurrence = lambda: _edge(a, c)
                    else:
                        b, c = first, second
                        a = _phi(b, multiplier_one)
                        residual = c - _psi(a, multiplier_two)
                        reverse_occurrence = lambda: _edge(b, c)

                    if residual:
                        inconsistent_primes.update(
                            _element_prime_factors(residual)
                        )
                        continue
                    consistent[position] += 1
                    atom = (a, b, c)
                    relation_change_primes.update(
                        _atom_relation_change_primes(atom)
                    )
                    if len(set(atom)) < 3:
                        continue
                    if _atom_alignment_values(atom) != expected_values[position]:
                        continue
                    if _normalized_target_edge(reverse_occurrence()) != reverse_edge:
                        raise AssertionError("R1 reverse coordinate changed")
                    support = frozenset(
                        (set(hh_support) - {reverse_edge})
                        | set(_normalized_aligned_support(atom))
                    )
                    if len(support) != 5:
                        continue
                    valid[position] += 1
                    atoms.add(atom)
                    supports.add(support)

    return {
        "formal_branch_count": sum(attempts.values()),
        "attempts_by_reverse_position": attempts,
        "consistent_branches_by_reverse_position": consistent,
        "valid_branches_by_reverse_position": valid,
        "valid_atom_count": len(atoms),
        "supports": frozenset(supports),
        "inconsistent_relation_exception_primes": sorted(
            inconsistent_primes
        ),
        "valid_filter_relation_change_primes": sorted(
            relation_change_primes
        ),
    }


def _positive_outside_high_audit() -> dict[str, object]:
    """Classify the high type omitted by the PN parameterization.

    If both negative occurrences are aligned, write
    ``c=Psi_s(a)=Psi_t(b)``.  Equal multipliers force ``a=b``.  Unequal
    multipliers give ``b=Psi_(s/t)(a)``, so the positive edge is reverse,
    not outside, unless it is the deleted self-antipodal edge.  The two
    ordered self-antipodal solutions are the same capped singleton.
    """
    repeated_label_families = 0
    reverse_target_families = 0
    self_antipodal_candidates: list[tuple[Q3, Q3, Q3]] = []
    relation_change_primes: set[int] = set()
    for first_multiplier, second_multiplier in product(POWERS, repeat=2):
        if first_multiplier == second_multiplier:
            repeated_label_families += 1
            continue
        ratio = first_multiplier / second_multiplier
        # Equality of the two aligned negative endpoints gives
        # b=Psi_ratio(a), identically a reverse positive-edge occurrence.
        reverse_target_families += 1
        denominator = 1 + ratio
        relation_change_primes.update(
            _element_prime_factors(denominator)
        )
        a = (1 - ratio) / denominator
        b = -a
        c = _psi(a, first_multiplier)
        atom = (a, b, c)
        relation_change_primes.update(_atom_relation_change_primes(atom))
        if len(set(atom)) != 3 or _atom_alignment_values(atom) != (0, 1, 1):
            raise AssertionError("positive-outside high classification changed")
        self_antipodal_candidates.append(atom)

    supports = {
        _normalized_aligned_support(atom)
        for atom in self_antipodal_candidates
    }
    return {
        "formal_multiplier_pairs": 4,
        "equal_multiplier_repeated_label_families": repeated_label_families,
        "unequal_multiplier_reverse_target_families": reverse_target_families,
        "self_antipodal_ordered_candidates": len(self_antipodal_candidates),
        "distinct_self_antipodal_supports": len(supports),
        "support": next(iter(supports)),
        "relation_change_exception_primes": sorted(relation_change_primes),
    }


def _positive_excess_tuple_audit() -> dict[str, object]:
    """List every nonnegative integer tuple allowed by the excess identity."""
    rows = []
    for compact_count in range(1, 5):
        for all_equal_count in range(2):
            for cycle_rank in range(2):
                for cap_count in range(3):
                    excess = (
                        4
                        - compact_count
                        - 4 * all_equal_count
                        - 4 * cycle_rank
                        - 2 * cap_count
                    )
                    if excess > 0:
                        rows.append(
                            (
                                compact_count,
                                all_equal_count,
                                cycle_rank,
                                cap_count,
                                excess,
                            )
                        )
    expected = [
        (1, 0, 0, 0, 3),
        (1, 0, 0, 1, 1),
        (2, 0, 0, 0, 2),
        (3, 0, 0, 0, 1),
    ]
    return {
        "excess_identity": "K-2*delta=4-K-4*AE-4*cycle_rank-2*caps",
        "positive_integer_tuples_K_AE_cycle_rank_caps_excess": rows,
        "isolated_K1_cap0_requires_score_three_and_is_excluded": True,
        "proved": rows == expected,
    }


def _cap_support() -> dict[str, object]:
    candidates = []
    good_supports = []
    bad_duplicate_count = 0
    relation_change_primes: set[int] = set()
    for phi_multiplier, psi_multiplier in product(POWERS, POWERS):
        root = (phi_multiplier - psi_multiplier) / (
            phi_multiplier + psi_multiplier
        )
        atom = (
            root,
            _phi(root, phi_multiplier),
            _psi(root, psi_multiplier),
        )
        candidates.append(atom)
        relation_change_primes.update(_atom_relation_change_primes(atom))
        support = _normalized_aligned_support(atom)
        if len(support) == 1:
            bad_duplicate_count += 1
        elif len(support) == 2:
            good_supports.append(support)
        else:
            raise AssertionError("unexpected generic cap support")

    k = Q2 - Q
    exceptional_atom = (k, -k, Q3(-3))
    candidates.append(exceptional_atom)
    relation_change_primes.update(
        _atom_relation_change_primes(exceptional_atom)
    )
    exceptional_support = _normalized_aligned_support(exceptional_atom)
    if len(exceptional_support) != 2:
        raise AssertionError("exceptional NN cap is not valid")
    good_supports.append(exceptional_support)
    if len(set(good_supports)) != 1:
        raise AssertionError("valid caps do not share one support")
    return {
        "candidate_count": len(candidates),
        "bad_duplicate_count": bad_duplicate_count,
        "good_candidate_count": len(good_supports),
        "support": good_supports[0],
        "filter_relation_change_primes": sorted(relation_change_primes),
    }


def _edge_text(edge: Edge) -> tuple[str, str]:
    return tuple(sorted(value.text() for value in edge))  # type: ignore[return-value]


def _support_text(support: frozenset[Edge]) -> list[list[str]]:
    return [list(edge) for edge in sorted(_edge_text(edge) for edge in support)]


def _maximum_disjoint_weight(blocks: list[tuple[str, int, frozenset[Edge]]]) -> tuple[int, list[list[str]]]:
    best_weight = -1
    maximizers: list[list[str]] = []
    for mask in range(1 << len(blocks)):
        chosen = [blocks[index] for index in range(len(blocks)) if mask >> index & 1]
        supports = [block[2] for block in chosen]
        if any(left & right for left, right in combinations(supports, 2)):
            continue
        weight = sum(block[1] for block in chosen)
        names = [block[0] for block in chosen]
        if weight > best_weight:
            best_weight = weight
            maximizers = [names]
        elif weight == best_weight:
            maximizers.append(names)
    return best_weight, maximizers


def _score_three_compact_no_go() -> dict[str, object]:
    """Solve all eight formal score-three compact compatibility cases."""
    counts = {
        "degenerate_repeated_label": 0,
        "inconsistent": 0,
        "missing_self_antipodal_edge": 0,
        "actual_score_three": 0,
    }
    relation_change_primes: set[int] = set()
    for phi_multiplier, psi_first, psi_second in product(
        POWERS, POWERS, POWERS
    ):
        # psi_first(a)=psi_second(phi_multiplier(a)).
        coefficient = psi_first - psi_second * phi_multiplier
        constant = (psi_first - 1) - (
            psi_second * (1 - phi_multiplier) + psi_second - 1
        )
        if not coefficient:
            if constant:
                relation_change_primes.update(
                    _element_prime_factors(constant)
                )
                counts["inconsistent"] += 1
                continue
            raise AssertionError("unexpected one-parameter score-three family")
        root = -constant / coefficient
        atom = (
            root,
            _phi(root, phi_multiplier),
            _psi(root, psi_first),
        )
        relation_change_primes.update(_atom_relation_change_primes(atom))
        if len(set(atom)) < 3:
            counts["degenerate_repeated_label"] += 1
            continue
        values = _atom_alignment_values(atom)
        if values.count(1) == 2 and values.count(0) == 1:
            counts["missing_self_antipodal_edge"] += 1
        elif sum(values) == 3:
            counts["actual_score_three"] += 1
        else:
            raise AssertionError("unclassified score-three compatibility case")
    return {
        **counts,
        "formal_case_count": sum(counts.values()),
        "distinct_label_score_three_compact_exists": bool(
            counts["actual_score_three"]
        ),
        "filter_relation_change_primes": sorted(relation_change_primes),
        "proved": counts
        == {
            "degenerate_repeated_label": 4,
            "inconsistent": 2,
            "missing_self_antipodal_edge": 2,
            "actual_score_three": 0,
        },
    }


def _edge_collision_primes(first: Edge, second: Edge) -> set[int]:
    """Possible characteristics in which two unequal symbolic edges merge."""
    if first == second:
        return set()
    left = tuple(first)
    right = tuple(second)
    result: set[int] = set()
    for sign in (1, -1):
        signed_right = tuple(sign * value for value in right)
        for ordered_right in (signed_right, tuple(reversed(signed_right))):
            norms = [
                (left[index] - ordered_right[index]).norm()
                for index in range(2)
            ]
            numerator_gcd = gcd(
                abs(norms[0].numerator), abs(norms[1].numerator)
            )
            result.update(_prime_factors(numerator_gcd))
            result.update(_prime_factors(norms[0].denominator))
            result.update(_prime_factors(norms[1].denominator))
    return result


def _support_collision_primes(supports: list[frozenset[Edge]]) -> list[int]:
    result: set[int] = set()
    for support in supports:
        for edge in support:
            first, second = tuple(edge)
            result.update(_element_prime_factors(first - second))
            result.update(_element_prime_factors(first + second))
        for first, second in combinations(support, 2):
            result.update(_edge_collision_primes(first, second))
    for first_support, second_support in combinations(supports, 2):
        for first in first_support:
            for second in second_support:
                if first != second:
                    result.update(_edge_collision_primes(first, second))
    return sorted(result)


def equianharmonic_component_packing_certificate() -> dict[str, object]:
    """Return the exact positive-excess component and packing certificate."""
    o2 = _o2_free_supports()
    hh = _hh_free_support()
    cap = _cap_support()
    score_three = _score_three_compact_no_go()
    positive_outside_high = _positive_outside_high_audit()
    excess_tuples = _positive_excess_tuple_audit()
    o2_supports = o2["supports"]
    hh_support = hh["support"]
    cap_support = cap["support"]
    if not isinstance(o2_supports, frozenset):
        raise AssertionError("bad O2 support type")
    if not isinstance(hh_support, frozenset) or not isinstance(cap_support, frozenset):
        raise AssertionError("bad support type")
    r1_flips = _r1_flip_supports(hh_support)
    r1 = _r1_exhaustive_supports(hh_support)
    r1_supports = r1["supports"]
    if not isinstance(r1_supports, frozenset):
        raise AssertionError("bad R1 support type")

    o2_rows = sorted(o2_supports, key=lambda support: _support_text(support))
    pairwise_o2 = [len(left & right) for left, right in combinations(o2_rows, 2)]
    hh_o2 = [len(hh_support & support) for support in o2_rows]
    cap_o2 = [len(cap_support & support) for support in o2_rows]
    blocks = [("HH", 2, hh_support), ("cap", 1, cap_support)]
    blocks.extend((f"F{index + 1}", 1, support) for index, support in enumerate(o2_rows))
    maximum_weight, maximizers = _maximum_disjoint_weight(blocks)
    collision_primes = _support_collision_primes(
        o2_rows + [hh_support, cap_support]
    )
    relation_change_primes = sorted(
        set(o2["rank_two_filter_relation_change_primes"])
        | set(hh["rank_two_filter_relation_change_primes"])
        | set(cap["filter_relation_change_primes"])
        | set(score_three["filter_relation_change_primes"])
        | set(r1["inconsistent_relation_exception_primes"])
        | set(r1["valid_filter_relation_change_primes"])
        | set(positive_outside_high["relation_change_exception_primes"])
    )

    proved = (
        o2["branch_counts"]
        == {"rank_two": 176, "rank_one": 8, "inconsistent": 72}
        and o2["valid_rank_two_assignments"] == 48
        and o2["rank_one_collision_count"] == 8
        and len(o2_supports) == 4
        and hh["valid_rank_two_assignments"] == 8
        and hh["rank_one_collision_count"] == hh["branch_counts"]["rank_one"]
        and r1_flips == r1_supports == o2_supports
        and r1["formal_branch_count"] == 192
        and r1["attempts_by_reverse_position"]
        == {"positive": 64, "negative_ac": 64, "negative_bc": 64}
        and r1["consistent_branches_by_reverse_position"]
        == {"positive": 8, "negative_ac": 0, "negative_bc": 0}
        and r1["valid_branches_by_reverse_position"]
        == {"positive": 8, "negative_ac": 0, "negative_bc": 0}
        and r1["valid_atom_count"] == 8
        and cap["candidate_count"] == 5
        and cap["bad_duplicate_count"] == 2
        and cap["good_candidate_count"] == 3
        and pairwise_o2 == [2] * 6
        and hh_o2 == [3] * 4
        and not hh_support & cap_support
        and cap_o2 == [0] * 4
        and maximum_weight == 3
        and score_three["proved"]
        and excess_tuples["proved"]
        and positive_outside_high[
            "equal_multiplier_repeated_label_families"
        ]
        == 2
        and positive_outside_high[
            "unequal_multiplier_reverse_target_families"
        ]
        == 2
        and positive_outside_high["self_antipodal_ordered_candidates"] == 2
        and positive_outside_high["distinct_self_antipodal_supports"] == 1
        and positive_outside_high["support"] == cap_support
        and collision_primes == [2, 3, 7, 13]
        and relation_change_primes == [2, 3, 5, 7, 13, 19]
        and set(o2["rank_two_determinant_primes"]) <= {3}
        and set(o2["inconsistent_augmented_minor_primes"]) <= {2, 3}
        and set(hh["rank_two_determinant_primes"]) <= {3}
        and set(hh["inconsistent_augmented_minor_primes"]) <= {2, 3}
    )
    return {
        "coefficient_field": "Q(q)/(q^2+q+1)",
        "normalization": "X=2x",
        "valid_prime_scope": "p>=31, q^2+q+1=0 in F_p",
        "constant_target_sign_symmetry": (
            "negating every atom label replaces each physical edge by its "
            "antipode and negates the atom chain in a fixed orbit basis"
        ),
        "constant_target_signs_covered": [-1, 1],
        "component_identity": "delta=K+2*AE-2+2*cycle_rank+caps",
        "component_excess_identity": (
            "K-2*delta=4-K-4*AE-4*cycle_rank-2*caps"
        ),
        "distinct_component_support_rule": (
            "one aligned occurrence remains at each unit target coordinate, "
            "so different pairing components have disjoint unpaired supports"
        ),
        "nonexceptional_component_bound": "K<=2*delta",
        "positive_excess_types": {
            "HH": {"K": 2, "delta": 0, "weight": 2},
            "cap": {"K": 1, "delta": 0, "weight": 1},
            "O2_or_R1": {"K": 3, "delta": 1, "weight": 1},
        },
        "o2": {
            **{key: value for key, value in o2.items() if key != "supports"},
            "supports": [_support_text(support) for support in o2_rows],
        },
        "hh": {
            **{
                key: value
                for key, value in hh.items()
                if key not in ("support", "representative_atoms")
            },
            "support": _support_text(hh_support),
        },
        "cap": {
            **{key: value for key, value in cap.items() if key != "support"},
            "support": _support_text(cap_support),
        },
        "score_three_compact_no_go": score_three,
        "positive_excess_tuple_audit": excess_tuples,
        "positive_outside_high_audit": {
            **{
                key: value
                for key, value in positive_outside_high.items()
                if key != "support"
            },
            "support": _support_text(positive_outside_high["support"]),
        },
        "r1_exhaustive_audit": {
            **{key: value for key, value in r1.items() if key != "supports"},
            "supports": [
                _support_text(support)
                for support in sorted(
                    r1_supports, key=lambda support: _support_text(support)
                )
            ],
        },
        "r1_triangle_flips_equal_exhaustive_supports": (
            r1_flips == r1_supports
        ),
        "r1_supports_equal_o2_supports": r1_supports == o2_supports,
        "o2_pairwise_intersections": pairwise_o2,
        "hh_o2_intersections": hh_o2,
        "hh_cap_intersection": len(hh_support & cap_support),
        "cap_o2_intersections": cap_o2,
        "weighted_disjoint_packing_maximum": maximum_weight,
        "weighted_disjoint_packing_maximizers": maximizers,
        "orbit_support_collision_exception_primes": collision_primes,
        "discarded_branch_relation_change_exception_primes": (
            relation_change_primes
        ),
        "all_exceptional_characteristics_below_31": all(
            prime < 31
            for prime in (
                collision_primes
                + relation_change_primes
                + o2["rank_two_determinant_primes"]
                + o2["inconsistent_augmented_minor_primes"]
                + hh["rank_two_determinant_primes"]
                + hh["inconsistent_augmented_minor_primes"]
                + r1["inconsistent_relation_exception_primes"]
                + r1["valid_filter_relation_change_primes"]
                + positive_outside_high[
                    "relation_change_exception_primes"
                ]
            )
        ),
        "exact_fiber_threshold_dependency_ready": True,
        "global_common_form_threshold_asserted": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def equianharmonic_exact_fiber_threshold(
    p: int, compact_count: int
) -> dict[str, object]:
    """Prove the pure equianharmonic necessary compact-count threshold."""
    if not _is_prime(p) or p < 31 or p % 12 != 7:
        raise ValueError("need a prime p>=31 with p=7 mod 12")
    if (
        not isinstance(compact_count, int)
        or isinstance(compact_count, bool)
        or compact_count < 0
    ):
        raise ValueError("compact_count must be a nonnegative integer")
    r = (p - 3) // 4
    if compact_count > r:
        raise ValueError("need compact_count<=r")
    cycle_count = (p - 1) // 3
    all_equal_count = r - 1
    uncovered_cycle_count = cycle_count - all_equal_count
    minimum_compact_count = 2 * uncovered_cycle_count - 1
    total_deficit = 2 * compact_count - 3 * uncovered_cycle_count + 1
    required_positive_excess = compact_count - 2 * total_deficit
    packing = equianharmonic_component_packing_certificate()
    packing_maximum = int(packing["weighted_disjoint_packing_maximum"])
    below_threshold = compact_count < minimum_compact_count
    strict_packing_contradiction = bool(
        below_threshold and required_positive_excess > packing_maximum
    )
    proved = bool(
        packing["proved"]
        and r % 3 == 1
        and cycle_count == (4 * r + 2) // 3
        and uncovered_cycle_count == (r + 5) // 3
        and minimum_compact_count == (2 * r + 7) // 3
        and (
            not below_threshold
            or total_deficit < 0
            or strict_packing_contradiction
        )
    )
    if not proved:
        raise ArithmeticError("the equianharmonic threshold proof changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": compact_count,
        "all_equal_atom_count": all_equal_count,
        "phi_cycle_count": cycle_count,
        "uncovered_phi_cycle_count": uncovered_cycle_count,
        "target_score": p - 2,
        "total_deficit": total_deficit,
        "required_positive_excess": required_positive_excess,
        "component_packing_maximum": packing_maximum,
        "minimum_compact_atom_count": minimum_compact_count,
        "minimum_compact_atom_count_formula": "(2*r+7)/3",
        "constant_target_signs_covered": [-1, 1],
        "minus_constant_reduction": (
            "apply x->-x to every label; atom types are preserved and the "
            "signed chain is negated in the fixed edge-orbit coordinates"
        ),
        "below_threshold": below_threshold,
        "exact_fiber_excluded": below_threshold,
        "threshold_is_necessary_not_sufficient": True,
        "global_common_form_lift_asserted": False,
        "Boolean_lift_asserted": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def p3_odd_radon_centrality_component_upgrade(
    p: int, compact_count: int
) -> dict[str, object]:
    """Combine the support dichotomy with component packing centrality.

    If all odd contractions vanish, this proves centrality for
    ``3*b<=2*r+4``.  When ``p=11 mod 12`` the equianharmonic branch does
    not exist, so the same dependencies cover every ``b<=r``.
    """
    if not _is_prime(p) or p < 31 or p % 4 != 3:
        raise ValueError("need a prime p=3 mod 4 with p>=31")
    r = (p - 3) // 4
    if (
        not isinstance(compact_count, int)
        or isinstance(compact_count, bool)
        or not 0 <= compact_count <= r
    ):
        raise ValueError("need 0<=compact_count<=r")
    broad_no_equi_scope = p % 12 == 11
    bounded_equi_scope = 3 * compact_count <= 2 * r + 4
    if not broad_no_equi_scope and not bounded_equi_scope:
        raise ValueError("equianharmonic packing only covers 3*b<=2*r+4")

    from e1_gmin_m4_compact_ray_moment_gate import (
        p3_boundary_cubic_unit_reduction_certificate,
        p3_full_balanced_two_maximal_line_exclusion_certificate,
    )
    from e1_gmin_m4_conic_odd_radon import (
        conic_reduction_constants,
        theorem_record as conic_theorem_record,
    )

    two_line = p3_full_balanced_two_maximal_line_exclusion_certificate(
        p, compact_count
    )
    cubic = p3_boundary_cubic_unit_reduction_certificate(p)
    conic = conic_reduction_constants(p, compact_count)
    conic_theorem = conic_theorem_record()
    conic_flags = conic_theorem["proved"]
    required_conic_flags = {
        "conic_containing_low_weight_word_is_fully_conic_supported": True,
        "high_intersection_irreducible_conic_is_triangle_tangent": True,
        "nonconstant_affine_dual_weight_is_excluded_by_integer_l1": True,
        "only_constant_plus_or_minus_one_conic_word_survives": True,
        "star_constant_branch_is_excluded_by_quotient_parity": True,
        "nonequianharmonic_constant_branch_is_excluded": True,
        "constant_branch_forces_q_cubed_equals_one": True,
        "constant_branch_forces_p_congruent_7_mod_12": True,
    }
    conic_dichotomy_proved = all(
        conic_flags.get(name) is expected
        for name, expected in required_conic_flags.items()
    )
    equianharmonic = (
        None
        if broad_no_equi_scope
        else equianharmonic_exact_fiber_threshold(p, compact_count)
    )
    equianharmonic_branch_excluded = bool(
        broad_no_equi_scope
        or (
            equianharmonic is not None
            and equianharmonic["exact_fiber_excluded"]
        )
    )
    proved = bool(
        two_line["proved"]
        and two_line["all_supports_containing_h_collinear_points_excluded"]
        and cubic["proved"]
        and cubic["all_boundary_cubic_supports_excluded"]
        and conic["proved"]
        and conic_dichotomy_proved
        and equianharmonic_branch_excluded
    )
    if not proved:
        raise ArithmeticError("the upgraded odd-Radon centrality proof changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": compact_count,
        "compact_count_hypothesis": (
            "0<=b<=r" if broad_no_equi_scope else "3*b<=2*r+4"
        ),
        "p_mod_12": p % 12,
        "equianharmonic_branch_exists": not broad_no_equi_scope,
        "equianharmonic_branch_excluded": equianharmonic_branch_excluded,
        "equianharmonic_threshold": equianharmonic,
        "line_and_two_line_dependency_proved": two_line["proved"],
        "boundary_cubic_dependency_proved": cubic["proved"],
        "conic_peeling_dependency_proved": conic["proved"],
        "conic_dichotomy_dependency_proved": conic_dichotomy_proved,
        "required_conic_dichotomy_flags": required_conic_flags,
        "nonequianharmonic_conic_dependency": (
            "NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md, Section 4"
        ),
        "aggregate_signed_edge_chain_is_centrally_symmetric": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "global_common_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def _finite_edge(p: int, first: int, second: int) -> tuple[int, int]:
    first %= p
    second %= p
    if first == second:
        raise ValueError("loop")
    return tuple(sorted((first, second)))


def _finite_orbit_coordinate(
    p: int, edge: tuple[int, int]
) -> tuple[tuple[int, int] | None, int]:
    negative = _finite_edge(p, -edge[0], -edge[1])
    if edge == negative:
        return None, 0
    representative = min(edge, negative)
    return representative, 1 if edge == representative else -1


def _add_finite_orbit_edge(
    chain: dict[tuple[int, int], int],
    p: int,
    first: int,
    second: int,
    coefficient: int,
) -> None:
    representative, sign = _finite_orbit_coordinate(
        p, _finite_edge(p, first, second)
    )
    if representative is not None:
        chain[representative] = chain.get(representative, 0) + coefficient * sign


def p43_equianharmonic_threshold_witness_certificate() -> dict[str, object]:
    """Replay the exact p=43,b=9 fiber attaining the lower threshold."""
    from e1_gmin_m4_compact_ray_moment_gate import (
        all_equal_moment,
        compact_moment,
    )
    from e1_gmin_m4_conic_odd_radon import tangent_conic_target

    p = 43
    r = 10
    b = 9
    k = 13
    q = (1 - k) * pow(1 + k, -1, p) % p
    target = tangent_conic_target(p, k)
    chain: dict[tuple[int, int], int] = {}
    for atom in P43_AE_ATOMS:
        for first, second in combinations(atom, 2):
            _add_finite_orbit_edge(chain, p, first, second, 1)
    for atom, distinguished in P43_COMPACT_ATOMS:
        positive = [value for value in atom if value != distinguished]
        if len(positive) != 2:
            raise AssertionError("bad compact witness atom")
        _add_finite_orbit_edge(chain, p, positive[0], positive[1], 1)
        _add_finite_orbit_edge(chain, p, positive[0], distinguished, -1)
        _add_finite_orbit_edge(chain, p, positive[1], distinguished, -1)
    chain = {edge: value for edge, value in chain.items() if value}

    def moment_vector(degree: int) -> tuple[int, ...]:
        values = []
        for channel in range(degree // 2):
            value = sum(
                all_equal_moment(p, *atom, degree, channel)
                for atom in P43_AE_ATOMS
            )
            for atom, distinguished in P43_COMPACT_ATOMS:
                positive = [value for value in atom if value != distinguished]
                value += compact_moment(
                    p,
                    positive[0],
                    positive[1],
                    distinguished,
                    degree,
                    channel,
                )
            values.append(value % p)
        return tuple(values)

    odd_degrees = tuple(range(3, p - 1, 2))
    odd_vectors = {degree: moment_vector(degree) for degree in odd_degrees}
    degree_six = moment_vector(6)
    degree_eight = moment_vector(8)
    threshold = equianharmonic_exact_fiber_threshold(p, b)
    proved = bool(
        k * k % p == (-3) % p
        and pow(q, 3, p) == 1
        and q != 1
        and chain == target
        and len(chain) == p - 2
        and sum(abs(value) for value in chain.values()) == p - 2
        and all(not any(vector) for vector in odd_vectors.values())
        and degree_six == (37, 19, 8)
        and degree_eight == (18, 17, 10, 32)
        and b == threshold["minimum_compact_atom_count"]
    )
    if not proved:
        raise ArithmeticError("the p43 equianharmonic witness changed")
    return {
        "p": p,
        "r": r,
        "b": b,
        "k": k,
        "q": q,
        "ae_atoms": [list(atom) for atom in P43_AE_ATOMS],
        "compact_atoms": [
            {"triple": list(atom), "distinguished": distinguished}
            for atom, distinguished in P43_COMPACT_ATOMS
        ],
        "target_support": len(target),
        "target_l1": sum(abs(value) for value in target.values()),
        "edge_orbit_replay_exact": chain == target,
        "odd_degrees_checked": list(odd_degrees),
        "odd_channel_count": sum(degree // 2 for degree in odd_degrees),
        "all_odd_channels_zero": all(
            not any(vector) for vector in odd_vectors.values()
        ),
        "degree_six": list(degree_six),
        "degree_eight": list(degree_eight),
        "degree_six_and_eight_both_zero": False,
        "attains_equianharmonic_compact_threshold": True,
        "global_common_form_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(equianharmonic_component_packing_certificate(), indent=2))
