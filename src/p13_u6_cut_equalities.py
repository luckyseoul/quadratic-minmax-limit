#!/usr/bin/env python3
r"""Solver-free equality helpers for the open ``p=13,t=4,u=6`` row.

This module does not close the ``u=6`` residue.  It records two exact
consequences of the common-energy/translated-cut reduction:

* in the hard-excess partition ``1^5``, collision count ``C=2`` is
  impossible because the required two-unit energy deficit has no admissible
  row; and
* in the partition ``(2,2,1)``, collision count ``C=3`` is impossible by a
  quadratic-character product obstruction for a binary quartic with four
  projective roots.

The row catalogs here are deliberately solver-free.  They enumerate the
six-coordinate integer sphere with fixed sum and energy, and then apply all
74 pinned, one-sided translated-cut inequalities from Proposition 15.740.
The already established sharp upper energies are inputs to the equality
classification; this module does not independently prove that no row of
higher energy exists.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from itertools import combinations, permutations
from math import isqrt
from typing import Iterable

from e1_gmin_m4_prop15740 import translated_cut_vector_catalog


P = 13
DISTANCES = tuple(range(1, 7))
CUT_VECTOR_COUNT = 74
CUT_CATALOG_SHA256 = (
    "bfec2077a81acf1a6719caf93b066313445c55b4e2951c189d357731b437a265"
)
EMPTY_ROWS_SHA256 = hashlib.sha256(b"").hexdigest()


ROW_KINDS: dict[str, dict[str, object]] = {
    "H1": {
        "direction": "hard",
        "hard_excess": 1,
        "sum": 0,
        "l1_bound": 56,
        "cut_upper": 13,
        "sharp_energy": 28,
        "sharp_seed": (3, 2, 1, -1, -2, -3),
        "sharp_G_character": 1,
    },
    "H2": {
        "direction": "hard",
        "hard_excess": 2,
        "sum": -1,
        "l1_bound": 55,
        "cut_upper": 13,
        "sharp_energy": 63,
        "sharp_seed": (5, 3, 0, -2, -3, -4),
        "sharp_G_character": 1,
    },
    "O4": {
        "direction": "opposite",
        "opposite_parallel_count": 4,
        "sum": -9,
        "l1_bound": 57,
        "cut_upper": -52,
        "sharp_energy": 31,
        "sharp_seed": (1, 0, -1, -2, -3, -4),
        "sharp_G_character": -1,
    },
}


# These hashes use ``_rows_digest`` below.  The sphere-type hash is over the
# sorted coordinate multisets.  The admissible-row hash is over lexicographically
# sorted ordered rows after all 74 cuts.  The type/cut hash appends to each
# sorted type the minimum, over all its distinct permutations, of the maximum
# translated cut.  Making that last serialization explicit avoids relying on
# an opaque exploratory-report format.
EXPECTED_CATALOGS: dict[tuple[str, int], dict[str, object]] = {
    ("H1", 28): {
        "sphere_type_count": 5,
        "ordered_sphere_row_count": 1620,
        "sphere_type_sha256": (
            "b53e769083c366b37b80e717c25e2eb055f7a44e6b23d418c2f967a2fc7ddac7"
        ),
        "admissible_row_count": 6,
        "admissible_row_sha256": (
            "00da22b21560538fefac698faf41164ad4420b8a5f8612c24d85664f8ce3d074"
        ),
        "admissible_maximum_cuts": (12,),
    },
    ("H2", 63): {
        "sphere_type_count": 19,
        "ordered_sphere_row_count": 6360,
        "sphere_type_sha256": (
            "c6cfc7e5faad39a5ff6536ccc644a56ed7a62be6277f7b0f31da56e41d0e0f2a"
        ),
        "admissible_row_count": 6,
        "admissible_row_sha256": (
            "5cc45919fa60834028436a7f38230a31528bb8beac4d0e2b0bbd5faa08fec8a1"
        ),
        "admissible_maximum_cuts": (12,),
    },
    ("O4", 31): {
        "sphere_type_count": 3,
        "ordered_sphere_row_count": 1080,
        "sphere_type_sha256": (
            "7cec6c3503c5c93e36c789040d52c477efc6bcae58322cd0fd08b9f096582d03"
        ),
        "admissible_row_count": 6,
        "admissible_row_sha256": (
            "3652c03213140243465d4310070fdae42f846eaffae7c0f2e6b2b89e6ddd171e"
        ),
        "admissible_maximum_cuts": (-54,),
    },
    ("H1", 26): {
        "sphere_type_count": 7,
        "ordered_sphere_row_count": 1560,
        "sphere_type_sha256": (
            "ecf6aab14bbf202dcee24ce43c5e89d197d1925f0e50237d42b693f9c33920a1"
        ),
        "type_minimax_cut_sha256": (
            "db91e151d1414746bf45ff7db7a6804ed92e7c496cc617e6649a33d0751b0abc"
        ),
        "admissible_row_count": 0,
        "admissible_row_sha256": EMPTY_ROWS_SHA256,
        "expected_type_minimum_maximum_cuts": (
            ((-4, -1, 0, 1, 2, 2), 18),
            ((-4, 0, 0, 0, 1, 3), 16),
            ((-3, -2, -1, 2, 2, 2), 16),
            ((-3, -2, 0, 0, 2, 3), 14),
            ((-3, -1, 0, 0, 0, 4), 20),
            ((-2, -2, -2, 1, 2, 3), 18),
            ((-2, -2, -1, 0, 1, 4), 16),
        ),
    },
    ("O4", 29): {
        "sphere_type_count": 6,
        "ordered_sphere_row_count": 540,
        "sphere_type_sha256": (
            "3df77f1cbb288c81288a40c03ebfda7101303a4ac54bc6d5e8a543418fd7b65c"
        ),
        "type_minimax_cut_sha256": (
            "118851723cbd6fd01e1498a847560c86e5ca4d92ec017eac7da69d1a33c44547"
        ),
        "admissible_row_count": 0,
        "admissible_row_sha256": EMPTY_ROWS_SHA256,
        "expected_type_minimum_maximum_cuts": (
            ((-5, -1, -1, -1, -1, 0), -46),
            ((-4, -3, -2, 0, 0, 0), -50),
            ((-4, -3, -1, -1, -1, 1), -50),
            ((-4, -2, -2, -2, 0, 1), -48),
            ((-3, -3, -3, -1, 0, 1), -50),
            ((-3, -2, -2, -2, -2, 2), -46),
        ),
    },
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _rows_digest(rows: Iterable[tuple[int, ...]]) -> str:
    payload = ";".join(",".join(map(str, row)) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def translated_cut_vectors() -> tuple[tuple[int, ...], ...]:
    """Import and pin the canonical 74 translated-cut vectors."""
    catalog = translated_cut_vector_catalog()
    vectors = tuple(tuple(int(value) for value in row) for row in catalog["vectors"])
    digest = _rows_digest(vectors)
    proved = bool(
        catalog["proved"]
        and len(vectors) == CUT_VECTOR_COUNT
        and digest == CUT_CATALOG_SHA256
        and all(len(row) == 6 and sum(row) == 42 for row in vectors)
    )
    _require(proved, "the pinned p=13 translated-cut catalog changed")
    return vectors


def _sorted_sphere_types(total: int, energy: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate nondecreasing integral six-tuples of fixed sum and norm."""
    if energy < 0:
        raise ValueError("energy must be nonnegative")
    rows: list[tuple[int, ...]] = []

    def visit(
        prefix: tuple[int, ...],
        lower: int,
        remaining_sum: int,
        remaining_energy: int,
        slots: int,
    ) -> None:
        if slots == 0:
            if remaining_sum == 0 and remaining_energy == 0:
                rows.append(prefix)
            return
        radius = isqrt(remaining_energy)
        for value in range(max(lower, -radius), radius + 1):
            next_sum = remaining_sum - value
            next_energy = remaining_energy - value * value
            next_slots = slots - 1
            if next_slots == 0:
                if next_sum or next_energy:
                    continue
            else:
                # Later coordinates are at least ``value`` and Cauchy gives
                # the other exact, integer-only pruning inequality.
                if next_sum < next_slots * value:
                    continue
                if next_sum * next_sum > next_slots * next_energy:
                    continue
            visit(
                prefix + (value,),
                value,
                next_sum,
                next_energy,
                next_slots,
            )

    visit((), -isqrt(energy), total, energy, 6)
    return tuple(rows)


def _distinct_permutations(row: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(set(permutations(row))))


def translated_cut_values(row: Iterable[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in row)
    if len(values) != 6:
        raise ValueError("a p=13 distance row must have six coordinates")
    return tuple(
        sum(coefficient * value for coefficient, value in zip(cut, values))
        for cut in translated_cut_vectors()
    )


def maximum_translated_cut(row: Iterable[int]) -> int:
    """Return ``max_c c.q``; the inequalities are intentionally one-sided."""
    return max(translated_cut_values(row))


def moment_residues(row: Iterable[int]) -> tuple[int, int]:
    values = tuple(int(value) for value in row)
    if len(values) != 6:
        raise ValueError("a p=13 distance row must have six coordinates")
    return tuple(
        sum(pow(distance, degree, P) * value for distance, value in zip(DISTANCES, values))
        % P
        for degree in (2, 4)
    )  # type: ignore[return-value]


def sign_safe_quartic_value(kind: str, row: Iterable[int]) -> int:
    """Evaluate the global quartic in hard/opposite normalized coordinates."""
    if kind not in ROW_KINDS:
        raise ValueError(f"unknown row kind {kind!r}")
    n2, n4 = moment_residues(row)
    if ROW_KINDS[kind]["direction"] == "hard":
        return (n4 - n2 * n2) % P
    # If N_j=(-h)M_j, then h*M4-M2^2=-(N4+N2^2).
    return (-n4 - n2 * n2) % P


def quadratic_character(value: int, p: int = P) -> int:
    residue = value % p
    if residue == 0:
        return 0
    symbol = pow(residue, (p - 1) // 2, p)
    _require(symbol in (1, p - 1), "Euler's criterion returned a bad symbol")
    return 1 if symbol == 1 else -1


def _distance_bin(value: int) -> int:
    residue = value % P
    if residue == 0:
        raise ValueError("a distance multiplier must be nonzero")
    return min(residue, P - residue)


def multiply_distance_row(row: Iterable[int], multiplier: int) -> tuple[int, ...]:
    """Relabel a distance row by ``a -> multiplier*a`` in ``F_13^*/+-1``."""
    values = tuple(int(value) for value in row)
    if len(values) != 6:
        raise ValueError("a p=13 distance row must have six coordinates")
    if multiplier % P == 0:
        raise ValueError("the multiplier must be nonzero modulo 13")
    return tuple(values[_distance_bin(multiplier * distance) - 1] for distance in DISTANCES)


@lru_cache(maxsize=None)
def fixed_energy_row_catalog(kind: str, energy: int) -> dict[str, object]:
    """Classify one exact sphere under the pinned one-sided cut system."""
    if kind not in ROW_KINDS:
        raise ValueError(f"unknown row kind {kind!r}")
    if not isinstance(energy, int) or isinstance(energy, bool) or energy < 0:
        raise ValueError("energy must be a nonnegative integer")
    spec = ROW_KINDS[kind]
    total = int(spec["sum"])
    l1_bound = int(spec["l1_bound"])
    cut_upper = int(spec["cut_upper"])
    sphere_types = _sorted_sphere_types(total, energy)
    type_records: list[dict[str, object]] = []
    ordered_rows: list[tuple[int, ...]] = []
    admissible_rows: list[tuple[int, ...]] = []
    for sphere_type in sphere_types:
        rows = _distinct_permutations(sphere_type)
        row_cut_maxima = tuple(maximum_translated_cut(row) for row in rows)
        minimum_maximum_cut = min(row_cut_maxima)
        admitted = tuple(
            row
            for row, cut_maximum in zip(rows, row_cut_maxima)
            if sum(abs(value) for value in row) <= l1_bound
            and cut_maximum <= cut_upper
        )
        ordered_rows.extend(rows)
        admissible_rows.extend(admitted)
        type_records.append(
            {
                "sorted_type": list(sphere_type),
                "distinct_permutation_count": len(rows),
                "minimum_over_permutations_of_maximum_cut": minimum_maximum_cut,
                "permutations_attaining_minimum_count": row_cut_maxima.count(
                    minimum_maximum_cut
                ),
                "admissible_permutation_count": len(admitted),
            }
        )

    ordered = tuple(sorted(ordered_rows))
    admissible = tuple(sorted(admissible_rows))
    type_cut_rows = tuple(
        tuple(int(value) for value in record["sorted_type"])
        + (int(record["minimum_over_permutations_of_maximum_cut"]),)
        for record in type_records
    )
    maximum_cuts = tuple(sorted({maximum_translated_cut(row) for row in admissible}))
    result: dict[str, object] = {
        "kind": kind,
        "direction": spec["direction"],
        "sum": total,
        "energy": energy,
        "l1_bound": l1_bound,
        "translated_cut_inequality": f"c.q<={cut_upper}",
        "translated_cut_upper": cut_upper,
        "translated_cut_count": len(translated_cut_vectors()),
        "translated_cut_catalog_sha256": CUT_CATALOG_SHA256,
        "sphere_type_count": len(sphere_types),
        "ordered_sphere_row_count": len(ordered),
        "sphere_types": [list(row) for row in sphere_types],
        "sphere_type_sha256": _rows_digest(sphere_types),
        "type_cut_records": type_records,
        "type_minimax_cut_sha256": _rows_digest(type_cut_rows),
        "admissible_rows": [list(row) for row in admissible],
        "admissible_row_count": len(admissible),
        "admissible_row_sha256": _rows_digest(admissible),
        "admissible_maximum_cuts": list(maximum_cuts),
        "solver_backend": "none; exact recursive integer-sphere enumeration",
        "fixed_energy_only": True,
    }

    expected = EXPECTED_CATALOGS.get((kind, energy))
    if expected is not None:
        checks = [
            result["sphere_type_count"] == expected["sphere_type_count"],
            result["ordered_sphere_row_count"]
            == expected["ordered_sphere_row_count"],
            result["sphere_type_sha256"] == expected["sphere_type_sha256"],
            result["admissible_row_count"] == expected["admissible_row_count"],
            result["admissible_row_sha256"] == expected["admissible_row_sha256"],
        ]
        if "admissible_maximum_cuts" in expected:
            checks.append(
                tuple(result["admissible_maximum_cuts"])
                == expected["admissible_maximum_cuts"]
            )
        if "type_minimax_cut_sha256" in expected:
            checks.append(
                result["type_minimax_cut_sha256"]
                == expected["type_minimax_cut_sha256"]
            )
        if "expected_type_minimum_maximum_cuts" in expected:
            live_records = tuple(
                (
                    tuple(int(value) for value in record["sorted_type"]),
                    int(record["minimum_over_permutations_of_maximum_cut"]),
                )
                for record in type_records
            )
            checks.append(
                live_records == expected["expected_type_minimum_maximum_cuts"]
            )
        proved = all(checks)
        _require(proved, f"the pinned {kind} energy-{energy} catalog changed")
        result["proved"] = proved
    else:
        result["proved"] = True
    return result


def _sharp_row_coset_certificate(kind: str) -> dict[str, object]:
    spec = ROW_KINDS[kind]
    energy = int(spec["sharp_energy"])
    catalog = fixed_energy_row_catalog(kind, energy)
    rows = tuple(tuple(int(value) for value in row) for row in catalog["admissible_rows"])
    seed = tuple(int(value) for value in spec["sharp_seed"])
    orbit = tuple(sorted({multiply_distance_row(seed, multiplier) for multiplier in DISTANCES}))
    quartic_values = tuple(sorted(sign_safe_quartic_value(kind, row) for row in rows))
    characters = tuple(quadratic_character(value) for value in quartic_values)
    expected_values = (1, 1, 3, 3, 9, 9) if kind != "O4" else (2, 2, 5, 5, 6, 6)
    expected_character = int(spec["sharp_G_character"])
    proved = bool(
        catalog["proved"]
        and len(rows) == 6
        and orbit == rows
        and quartic_values == expected_values
        and characters == (expected_character,) * 6
        and all(
            sum(row) == int(spec["sum"])
            and sum(value * value for value in row) == energy
            and sum(abs(value) for value in row) <= int(spec["l1_bound"])
            and maximum_translated_cut(row) <= int(spec["cut_upper"])
            for row in rows
        )
    )
    _require(proved, f"the sharp {kind} row/coset certificate changed")
    return {
        "kind": kind,
        "sharp_energy": energy,
        "sharp_upper_bound_proved_here": False,
        "classification_scope": (
            "all equality rows conditional on the independently certified "
            "sharp row upper bound"
        ),
        "multiplicative_projective_coset_representatives": list(DISTANCES),
        "seed": list(seed),
        "orbit_is_exactly_the_six_admissible_rows": True,
        "ordered_row_sha256": catalog["admissible_row_sha256"],
        "sign_safe_global_quartic": (
            "G=N4-N2^2" if spec["direction"] == "hard" else "G=-(N4+N2^2)"
        ),
        "quartic_values_mod_13": list(quartic_values),
        "quadratic_characters": list(characters),
        "quartic_character_coset": expected_character,
        "proved": proved,
    }


def _projective_points(p: int) -> tuple[tuple[int, int], ...]:
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def _determinant(left: tuple[int, int], right: tuple[int, int], p: int) -> int:
    return (left[0] * right[1] - left[1] * right[0]) % p


@lru_cache(maxsize=1)
def quartic_character_product_certificate() -> dict[str, object]:
    r"""Prove and exhaustively audit the four-root character product at p=13.

    With normalized representatives ``(1,t)`` and ``(0,1)``, one has
    ``prod_(L!=r) det(L,r)=1``.  Pairing the ordered determinants inside a
    root set gives, for a squarefree degree-``d`` split form,

        prod_(L notin R) chi(F(L))
          = chi(c)^(p+1-d) chi(-1)^(d(d-1)/2).
    """
    points = _projective_points(P)
    determinant_products = []
    for root in points:
        value = 1
        for point in points:
            if point != root:
                value = value * _determinant(point, root, P) % P
        determinant_products.append(value)

    root_sets_checked = 0
    direct_products: list[int] = []
    for roots in combinations(points, 4):
        root_set = set(roots)
        for coefficient in (1, 2):  # representatives of chi(c)=+1 and -1
            character_product = 1
            for point in points:
                if point in root_set:
                    continue
                value = coefficient
                for root in roots:
                    value = value * _determinant(point, root, P) % P
                character_product *= quadratic_character(value)
            direct_products.append(character_product)
            root_sets_checked += 1

    degree = 4
    coefficient_exponent = P + 1 - degree
    minus_one_exponent = degree * (degree - 1) // 2
    formula_by_coefficient_character = {
        character: pow(character, coefficient_exponent)
        * pow(quadratic_character(-1), minus_one_exponent)
        for character in (-1, 1)
    }
    proved = bool(
        len(points) == 14
        and determinant_products == [1] * 14
        and _rows_digest(tuple((value,) for value in determinant_products))
        == "e60bc0aa0c7b044e9130b8808d3d3035597fe37c8b923afb04811e59ed5ba7b7"
        and coefficient_exponent == 10
        and minus_one_exponent == 6
        and quadratic_character(-1) == 1
        and formula_by_coefficient_character == {-1: 1, 1: 1}
        and root_sets_checked == 2 * 1001
        and set(direct_products) == {1}
    )
    _require(proved, "the universal quartic character-product audit changed")
    return {
        "p": P,
        "normalized_projective_points": [list(point) for point in points],
        "normalized_projective_point_count": len(points),
        "determinant_products": determinant_products,
        "determinant_product_sha256": _rows_digest(
            tuple((value,) for value in determinant_products)
        ),
        "general_formula": (
            "chi(c)^(p+1-d)*chi(-1)^(d*(d-1)/2)"
        ),
        "degree": degree,
        "coefficient_character_exponent": coefficient_exponent,
        "minus_one_character_exponent": minus_one_exponent,
        "formula_by_coefficient_character": {
            str(key): value for key, value in formula_by_coefficient_character.items()
        },
        "four_root_sets": 1001,
        "coefficient_character_representatives_per_root_set": 2,
        "direct_split_quartics_checked": root_sets_checked,
        "every_four_root_nonroot_character_product": 1,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def u6_energy_ledger_certificate() -> dict[str, object]:
    """Package the two new collision-boundary exclusions, not a u=6 close."""
    sharp = {kind: _sharp_row_coset_certificate(kind) for kind in ROW_KINDS}
    h1_deficit = fixed_energy_row_catalog("H1", 26)
    o4_deficit = fixed_energy_row_catalog("O4", 29)

    # Partition 1^5: two exact XNOR directions have already been subtracted
    # from the nonexact Parseval ledger.  At C=2 the exact total is two below
    # the sum of the twelve independent row maxima.
    one_five_upper = 5 * 28 + 7 * 31
    one_five_base = 303
    collision_increment = 26
    one_five_c2_target = one_five_base + 2 * collision_increment
    one_five_deficit = one_five_upper - one_five_c2_target
    parity_rule = all(
        sum(value * value for value in row) % 2 == sum(row) % 2
        for catalog in (
            fixed_energy_row_catalog("H1", 28),
            fixed_energy_row_catalog("O4", 31),
        )
        for row in catalog["admissible_rows"]
    )
    deficit_two_excluded = bool(
        one_five_deficit == 2
        and parity_rule
        and h1_deficit["admissible_row_count"] == 0
        and o4_deficit["admissible_row_count"] == 0
    )

    # Partition (2,2,1): four exact hard XNOR directions are the four roots
    # of the global binary quartic.  Equality at C=3 forces every other row
    # into its sharp catalog.
    two_two_one_upper = 2 * 63 + 28 + 7 * 31
    two_two_one_base = 293
    two_two_one_c3_target = two_two_one_base + 3 * collision_increment
    character_product = quartic_character_product_certificate()
    hard_nonroots = 3
    opposite_nonroots = 7
    equality_character_product = (
        int(sharp["H1"]["quartic_character_coset"])
        * int(sharp["H2"]["quartic_character_coset"]) ** 2
        * int(sharp["O4"]["quartic_character_coset"]) ** opposite_nonroots
    )
    quartic_equality_excluded = bool(
        two_two_one_upper == two_two_one_c3_target == 371
        and hard_nonroots + opposite_nonroots == 10
        and character_product["every_four_root_nonroot_character_product"] == 1
        and equality_character_product == -1
    )

    proved = bool(
        all(row["proved"] for row in sharp.values())
        and h1_deficit["proved"]
        and o4_deficit["proved"]
        and deficit_two_excluded
        and character_product["proved"]
        and quartic_equality_excluded
        and one_five_base + 3 * collision_increment > one_five_upper
        and two_two_one_base + 4 * collision_increment > two_two_one_upper
    )
    _require(proved, "the p13 u=6 equality ledger changed")
    return {
        "scope": "p=13,t=4,k=60,u=6 aggregate rows only",
        "result_status": "open reduction",
        "sharp_equality_catalogs": sharp,
        "sharp_row_upper_bounds_proved_here": False,
        "partition_1^5": {
            "hard_nonexact_rows": {"H1": 5},
            "opposite_nonexact_rows": {"O4": 7},
            "independent_row_energy_upper": one_five_upper,
            "parseval_base_at_C_zero": one_five_base,
            "collision_energy_increment": collision_increment,
            "C_equals_2_exact_energy": one_five_c2_target,
            "deficit_from_independent_maxima": one_five_deficit,
            "integer_square_parity_forces_each_row_deficit_even": parity_rule,
            "only_possible_deficit_two_rows": ["H1 energy 26", "O4 energy 29"],
            "H1_energy_26_catalog": h1_deficit,
            "O4_energy_29_catalog": o4_deficit,
            "C_equals_2_excluded": deficit_two_excluded,
            "C_at_least_3_excluded_by_energy": True,
            "remaining_collision_counts": [0, 1],
        },
        "partition_2_2_1": {
            "hard_nonexact_rows": {"H1": 1, "H2": 2},
            "opposite_nonexact_rows": {"O4": 7},
            "exact_XNOR_quartic_roots": 4,
            "independent_row_energy_upper": two_two_one_upper,
            "parseval_base_at_C_zero": two_two_one_base,
            "collision_energy_increment": collision_increment,
            "C_equals_3_exact_energy": two_two_one_c3_target,
            "equality_forces_every_row_sharp": True,
            "hard_nonroot_quartic_characters": [1] * hard_nonroots,
            "opposite_nonroot_quartic_characters": [-1] * opposite_nonroots,
            "forced_nonroot_character_product": equality_character_product,
            "universal_four_root_character_product": character_product,
            "C_equals_3_excluded": quartic_equality_excluded,
            "C_at_least_4_excluded_by_energy": True,
            "remaining_collision_counts": [0, 1, 2],
        },
        "p13_t4_u6_closed": False,
        "proved": proved,
    }


__all__ = [
    "CUT_CATALOG_SHA256",
    "CUT_VECTOR_COUNT",
    "EXPECTED_CATALOGS",
    "P",
    "ROW_KINDS",
    "fixed_energy_row_catalog",
    "maximum_translated_cut",
    "moment_residues",
    "multiply_distance_row",
    "quadratic_character",
    "quartic_character_product_certificate",
    "sign_safe_quartic_value",
    "translated_cut_values",
    "translated_cut_vectors",
    "u6_energy_ledger_certificate",
]
