#!/usr/bin/env python3
"""Exact fractional obstruction to a fixed global semimetric separator.

At the p=31 top hard-fixed endpoint, average uniformly over graph edges
subject only to the exact parallel, fixed-edge, and no-double ledgers, and
average each nonzero hard center uniformly.  The resulting fractional graph
has an explicit rational compact-atom decomposition in every hard row and an
explicit six-triangle-plus-compact decomposition in every opposite row.

Consequently no *fixed sum* of the 32 row semimetric inequalities can be
positive on every graph/center allowed by this ledger relaxation.  This is a
method barrier, not an integral common graph and not a residual-(ii) witness.
Adaptive row separation and the integral/coupled atom fibre remain open.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations

from e1_gmin_m4_p31_top_j1_f3_case_split import (
    DIRECTION_SIGNS,
    HARD,
    TARGET_PROFILE,
)


P = 31
H = (P - 1) // 2
FIXED_DIRECTION = 1
EDGE_CATEGORIES = (
    "zero_incident",
    "nonzero_antipodal",
    "nonzero_nonantipodal",
)
CATEGORY_SIZES = (P - 1, H, 2 * H * (H - 1))
NONFIXED_EDGES_PER_DIRECTION = 14_400
NONFIXED_ORBITS_PER_DIRECTION = NONFIXED_EDGES_PER_DIRECTION // 2

# A type records the total coefficient of one atom on the three edge
# categories above.  Uniform averaging over all atoms of one type makes its
# coefficient on an individual category-i edge equal to type[i]/size[i].
K_GENERIC_NEGATIVE = (0, 0, -1)
K_ANTIPODAL_NEGATIVE = (0, -1, 0)
K_ZERO_NEGATIVE = (-2, 0, 1)
K_ANTIPODAL_POSITIVE = (0, 1, -2)
K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE = (-2, 1, 0)
T_GENERIC = (0, 0, 3)
T_ZERO = (2, 0, 1)
T_ANTIPODAL = (0, 1, 2)
T_ZERO_ANTIPODAL = (2, 1, 0)


DECOMPOSITION_WEIGHTS: dict[
    str, dict[str, tuple[tuple[tuple[int, int, int], Fraction], ...]]
] = {
    "hard_q15": {
        "compact": (
            (K_ANTIPODAL_NEGATIVE, Fraction(10_037, 1_920)),
            (K_ZERO_NEGATIVE, Fraction(853, 960)),
            (K_ANTIPODAL_POSITIVE, Fraction(11_297, 1_920)),
        ),
        "positive": (),
    },
    "hard_q14_fixed": {
        "compact": (
            (K_GENERIC_NEGATIVE, Fraction(9_499, 960)),
            (K_ANTIPODAL_NEGATIVE, Fraction(9, 32)),
            (K_ZERO_NEGATIVE, Fraction(791, 960)),
        ),
        "positive": (),
    },
    "hard_q14_other": {
        "compact": (
            (K_ANTIPODAL_NEGATIVE, Fraction(1_513, 320)),
            (K_ZERO_NEGATIVE, Fraction(137, 160)),
            (K_ANTIPODAL_POSITIVE, Fraction(1_733, 320)),
        ),
        "positive": (),
    },
    "opposite_q16": {
        "compact": (
            (K_ANTIPODAL_NEGATIVE, Fraction(1_201, 240)),
            (
                K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE,
                Fraction(479, 240),
            ),
        ),
        "positive": (
            (T_ZERO_ANTIPODAL, Fraction(143, 60)),
            (T_GENERIC, Fraction(217, 60)),
        ),
    },
    "opposite_q15": {
        "compact": (
            (K_ANTIPODAL_NEGATIVE, Fraction(12_493, 2_880)),
            (
                K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE,
                Fraction(4_787, 2_880),
            ),
        ),
        "positive": (
            (T_ZERO_ANTIPODAL, Fraction(1_499, 720)),
            (T_GENERIC, Fraction(2_821, 720)),
        ),
    },
}


def _edge_category(first: int, second: int) -> int:
    if not 0 <= first < second < P:
        raise ValueError("need an ordered pair of distinct p31 labels")
    if first == 0:
        return 0
    if (first + second) % P == 0:
        return 1
    return 2


def atom_type_counts() -> dict[str, dict[tuple[int, int, int], int]]:
    """Enumerate the invariant compact and positive-triangle type counts."""
    compact: Counter[tuple[int, int, int]] = Counter()
    positive: Counter[tuple[int, int, int]] = Counter()
    for first, second in combinations(range(P), 2):
        for distinguished in range(P):
            if distinguished in (first, second):
                continue
            vector = [0, 0, 0]
            for edge, coefficient in (
                ((first, second), 1),
                (tuple(sorted((first, distinguished))), -1),
                (tuple(sorted((second, distinguished))), -1),
            ):
                vector[_edge_category(*edge)] += coefficient
            compact[tuple(vector)] += 1
    for first, second, third in combinations(range(P), 3):
        vector = [0, 0, 0]
        for edge in ((first, second), (first, third), (second, third)):
            vector[_edge_category(*edge)] += 1
        positive[tuple(vector)] += 1
    return {
        "compact": dict(sorted(compact.items())),
        "positive": dict(sorted(positive.items())),
    }


def _row_class(direction_index: int) -> str:
    if not 0 <= direction_index < P + 1:
        raise ValueError("direction index is out of range")
    quota = TARGET_PROFILE[direction_index]
    if direction_index in HARD:
        if quota == 15:
            return "hard_q15"
        return (
            "hard_q14_fixed"
            if direction_index == FIXED_DIRECTION
            else "hard_q14_other"
        )
    return "opposite_q16" if quota == 16 else "opposite_q15"


def expected_row_category_coefficients(
    direction_index: int,
) -> tuple[Fraction, Fraction, Fraction]:
    """Return the exact three coefficients of the averaged transverse row.

    In spatial direction ``D`` there are 14,400 nonfixed actual edges.  For
    any different projection row, a generic label pair has 31 preimages and
    a nonzero antipodal pair has 30 nonfixed plus one fixed preimage.  The
    sole fixed edge is uniform over the 15 magnitudes in direction 1.
    """
    if not 0 <= direction_index < P + 1:
        raise ValueError("direction index is out of range")
    fixed = tuple(int(index == FIXED_DIRECTION) for index in range(P + 1))
    nonfixed_counts = tuple(
        TARGET_PROFILE[index] - fixed[index] for index in range(P + 1)
    )
    epsilon = DIRECTION_SIGNS[direction_index]
    signed_nonfixed = sum(
        DIRECTION_SIGNS[index] * nonfixed_counts[index]
        for index in range(P + 1)
        if index != direction_index
    )
    signed_fixed = sum(
        DIRECTION_SIGNS[index] * fixed[index]
        for index in range(P + 1)
        if index != direction_index
    )
    generic = Fraction(
        epsilon * P * signed_nonfixed,
        NONFIXED_EDGES_PER_DIRECTION,
    )
    antipodal = Fraction(
        epsilon * (P - 1) * signed_nonfixed,
        NONFIXED_EDGES_PER_DIRECTION,
    ) + Fraction(epsilon * signed_fixed, H)
    zero_incident = generic
    nonzero_generic = generic
    if direction_index in HARD:
        # The literal-star residual adds +S_j.  Averaging j uniformly over
        # F_p^* hits {0,a} once and {a,b}, a*b!=0, twice.
        zero_incident += Fraction(1, P - 1)
        antipodal += Fraction(2, P - 1)
        nonzero_generic += Fraction(2, P - 1)
    return zero_incident, antipodal, nonzero_generic


def _averaged_type_vector(
    atom_type: tuple[int, int, int],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        Fraction(value, size)
        for value, size in zip(atom_type, CATEGORY_SIZES, strict=True)
    )  # type: ignore[return-value]


def replay_fractional_atom_decomposition(direction_index: int) -> dict[str, object]:
    """Replay one of the five exact invariant atom decompositions."""
    name = _row_class(direction_index)
    weights = DECOMPOSITION_WEIGHTS[name]
    reconstructed = [Fraction(0) for _ in EDGE_CATEGORIES]
    for family in ("compact", "positive"):
        for atom_type, weight in weights[family]:
            averaged = _averaged_type_vector(atom_type)
            for index, value in enumerate(averaged):
                reconstructed[index] += weight * value
    target = expected_row_category_coefficients(direction_index)
    compact_weight = sum(weight for _type, weight in weights["compact"])
    positive_weight = sum(weight for _type, weight in weights["positive"])
    quota = TARGET_PROFILE[direction_index]
    hard = direction_index in HARD
    expected_compact = quota - (3 if hard else 9)
    expected_positive = 0 if hard else 6
    target_sum = sum(
        coefficient * size
        for coefficient, size in zip(target, CATEGORY_SIZES, strict=True)
    )
    expected_sum = -expected_compact if hard else 18 - expected_compact
    proved = bool(
        tuple(reconstructed) == target
        and compact_weight == expected_compact
        and positive_weight == expected_positive
        and target_sum == expected_sum
        and all(weight >= 0 for family in weights.values() for _type, weight in family)
    )
    if not proved:
        raise ArithmeticError("a fractional row decomposition changed")
    return {
        "direction_index": direction_index,
        "direction_sign": DIRECTION_SIGNS[direction_index],
        "parallel_quota": quota,
        "row_class": name,
        "category_coefficients": tuple(str(value) for value in target),
        "category_sizes": CATEGORY_SIZES,
        "compact_weight": str(compact_weight),
        "positive_triangle_weight": str(positive_weight),
        "target_edge_sum": str(target_sum),
        "compact_type_weights": tuple(
            {"type": atom_type, "weight": str(weight)}
            for atom_type, weight in weights["compact"]
        ),
        "positive_type_weights": tuple(
            {"type": atom_type, "weight": str(weight)}
            for atom_type, weight in weights["positive"]
        ),
        "fractional_atom_cone_membership": True,
        "proved": proved,
    }


def global_fractional_semimetric_obstruction() -> dict[str, object]:
    """Return the exact hard-fixed obstruction to any fixed summed bank."""
    if FIXED_DIRECTION not in HARD or TARGET_PROFILE[FIXED_DIRECTION] != 14:
        raise ArithmeticError("the frozen hard-fixed representative changed")
    fixed_counts = tuple(int(index == FIXED_DIRECTION) for index in range(P + 1))
    nonfixed_counts = tuple(
        TARGET_PROFILE[index] - fixed_counts[index] for index in range(P + 1)
    )
    pair_occupancies = tuple(
        Fraction(count, NONFIXED_ORBITS_PER_DIRECTION)
        for count in nonfixed_counts
    )
    rows = tuple(
        replay_fractional_atom_decomposition(index) for index in range(P + 1)
    )
    type_counts = atom_type_counts()
    expected_compact_counts = {
        K_GENERIC_NEGATIVE: 11_760,
        K_ANTIPODAL_NEGATIVE: 870,
        K_ZERO_NEGATIVE: 420,
        K_ANTIPODAL_POSITIVE: 420,
        K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE: 15,
    }
    expected_positive_counts = {
        T_GENERIC: 3_640,
        T_ZERO: 420,
        T_ANTIPODAL: 420,
        T_ZERO_ANTIPODAL: 15,
    }
    proved = bool(
        sum(TARGET_PROFILE) == 479
        and sum(fixed_counts) == 1
        and sum(nonfixed_counts) == 478
        and max(pair_occupancies) < 1
        and type_counts["compact"] == expected_compact_counts
        and type_counts["positive"] == expected_positive_counts
        and all(row["proved"] for row in rows)
    )
    if not proved:
        raise ArithmeticError("the global fractional obstruction changed")
    return {
        "p": P,
        "endpoint": "top t=177 hard-fixed j=0 ledger relaxation",
        "fixed_direction": FIXED_DIRECTION,
        "fixed_direction_type": "hard",
        "fixed_direction_quota": TARGET_PROFILE[FIXED_DIRECTION],
        "parallel_profile": TARGET_PROFILE,
        "graph_edge_count": sum(TARGET_PROFILE),
        "fixed_edge_count": sum(fixed_counts),
        "unused_double_orbit_count": 0,
        "nonfixed_edge_count": sum(nonfixed_counts),
        "fractional_nonfixed_orbit_pair_occupancy_max": str(
            max(pair_occupancies)
        ),
        "hard_center_average": "uniform on F_31^*",
        "fractional_graph_is_in_parallel_fixed_no_double_polytope": True,
        "fractional_graph_is_convex_combination_of_ledger_graphs": True,
        "row_count": len(rows),
        "row_class_counts": dict(Counter(str(row["row_class"]) for row in rows)),
        "atom_type_counts": {
            family: tuple(
                {"type": atom_type, "count": count}
                for atom_type, count in rows_by_type.items()
            )
            for family, rows_by_type in type_counts.items()
        },
        "rows": rows,
        "all_expected_rows_in_fractional_atom_cones": True,
        "fixed_sum_of_row_semimetrics_separates_every_ledger_graph": False,
        "reason": (
            "a fixed linear separator positive on every ledger graph and "
            "nonzero hard-center choice would be positive on their uniform "
            "average, but every averaged row satisfies its exact real atom cone"
        ),
        "adaptive_semimetric_oracle_excluded": False,
        "integral_atom_synchronization_solved": False,
        "common_graph_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(global_fractional_semimetric_obstruction(), indent=2, sort_keys=True))
