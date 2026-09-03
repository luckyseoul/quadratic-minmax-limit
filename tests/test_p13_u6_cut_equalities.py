"""Fail-when-wrong tests for the p=13,u=6 cut-equality helpers."""
from __future__ import annotations

import pytest

import p13_u6_cut_equalities as M


def test_the_74_one_sided_cut_vectors_are_pinned() -> None:
    cuts = M.translated_cut_vectors()
    assert len(cuts) == 74
    assert M.CUT_CATALOG_SHA256 == (
        "bfec2077a81acf1a6719caf93b066313445c55b4e2951c189d357731b437a265"
    )
    assert all(len(row) == 6 and sum(row) == 42 for row in cuts)


@pytest.mark.parametrize(
    "kind,energy,type_count,type_hash,row_hash,cut_extremum",
    [
        (
            "H1",
            28,
            5,
            "b53e769083c366b37b80e717c25e2eb055f7a44e6b23d418c2f967a2fc7ddac7",
            "00da22b21560538fefac698faf41164ad4420b8a5f8612c24d85664f8ce3d074",
            12,
        ),
        (
            "H2",
            63,
            19,
            "c6cfc7e5faad39a5ff6536ccc644a56ed7a62be6277f7b0f31da56e41d0e0f2a",
            "5cc45919fa60834028436a7f38230a31528bb8beac4d0e2b0bbd5faa08fec8a1",
            12,
        ),
        (
            "O4",
            31,
            3,
            "7cec6c3503c5c93e36c789040d52c477efc6bcae58322cd0fd08b9f096582d03",
            "3652c03213140243465d4310070fdae42f846eaffae7c0f2e6b2b89e6ddd171e",
            -54,
        ),
    ],
)
def test_sharp_spheres_have_exactly_one_six_row_multiplicative_orbit(
    kind: str,
    energy: int,
    type_count: int,
    type_hash: str,
    row_hash: str,
    cut_extremum: int,
) -> None:
    row = M.fixed_energy_row_catalog(kind, energy)
    assert row["solver_backend"].startswith("none")
    assert row["sphere_type_count"] == type_count
    assert row["sphere_type_sha256"] == type_hash
    assert row["admissible_row_count"] == 6
    assert row["admissible_row_sha256"] == row_hash
    assert row["admissible_maximum_cuts"] == [cut_extremum]
    assert row["proved"] is True


def test_opposite_cut_is_not_silently_replaced_by_an_absolute_bound() -> None:
    row = (1, 0, -1, -2, -3, -4)
    values = M.translated_cut_values(row)
    assert min(values) == -98
    assert max(values) == -54 <= -52
    assert max(abs(value) for value in values) == 98 > 52
    assert row in tuple(
        tuple(values) for values in M.fixed_energy_row_catalog("O4", 31)["admissible_rows"]
    )


def test_all_seven_H1_deficit_two_types_are_present_and_fail_a_cut() -> None:
    expected = [
        ([-4, -1, 0, 1, 2, 2], 18),
        ([-4, 0, 0, 0, 1, 3], 16),
        ([-3, -2, -1, 2, 2, 2], 16),
        ([-3, -2, 0, 0, 2, 3], 14),
        ([-3, -1, 0, 0, 0, 4], 20),
        ([-2, -2, -2, 1, 2, 3], 18),
        ([-2, -2, -1, 0, 1, 4], 16),
    ]
    row = M.fixed_energy_row_catalog("H1", 26)
    live = [
        (
            record["sorted_type"],
            record["minimum_over_permutations_of_maximum_cut"],
        )
        for record in row["type_cut_records"]
    ]
    assert live == expected
    assert row["sphere_type_sha256"] == (
        "ecf6aab14bbf202dcee24ce43c5e89d197d1925f0e50237d42b693f9c33920a1"
    )
    assert row["type_minimax_cut_sha256"] == (
        "db91e151d1414746bf45ff7db7a6804ed92e7c496cc617e6649a33d0751b0abc"
    )
    assert all(minimum > 13 for _sphere_type, minimum in live)
    assert row["admissible_row_count"] == 0


def test_all_six_O4_deficit_two_types_fail_the_negative_cut_upper() -> None:
    expected = [
        ([-5, -1, -1, -1, -1, 0], -46),
        ([-4, -3, -2, 0, 0, 0], -50),
        ([-4, -3, -1, -1, -1, 1], -50),
        ([-4, -2, -2, -2, 0, 1], -48),
        ([-3, -3, -3, -1, 0, 1], -50),
        ([-3, -2, -2, -2, -2, 2], -46),
    ]
    row = M.fixed_energy_row_catalog("O4", 29)
    live = [
        (
            record["sorted_type"],
            record["minimum_over_permutations_of_maximum_cut"],
        )
        for record in row["type_cut_records"]
    ]
    assert live == expected
    assert row["sphere_type_sha256"] == (
        "3df77f1cbb288c81288a40c03ebfda7101303a4ac54bc6d5e8a543418fd7b65c"
    )
    assert row["type_minimax_cut_sha256"] == (
        "118851723cbd6fd01e1498a847560c86e5ca4d92ec017eac7da69d1a33c44547"
    )
    assert all(minimum > -52 for _sphere_type, minimum in live)
    assert row["admissible_row_count"] == 0


def test_sharp_G_values_land_in_the_required_character_cosets() -> None:
    ledger = M.u6_energy_ledger_certificate()
    sharp = ledger["sharp_equality_catalogs"]
    assert sharp["H1"]["quartic_values_mod_13"] == [1, 1, 3, 3, 9, 9]
    assert sharp["H2"]["quartic_values_mod_13"] == [1, 1, 3, 3, 9, 9]
    assert sharp["O4"]["quartic_values_mod_13"] == [2, 2, 5, 5, 6, 6]
    assert sharp["H1"]["quadratic_characters"] == [1] * 6
    assert sharp["H2"]["quadratic_characters"] == [1] * 6
    assert sharp["O4"]["quadratic_characters"] == [-1] * 6
    assert sharp["H1"]["sign_safe_global_quartic"] == "G=N4-N2^2"
    assert sharp["O4"]["sign_safe_global_quartic"] == "G=-(N4+N2^2)"
    assert all(
        record["orbit_is_exactly_the_six_admissible_rows"] is True
        for record in sharp.values()
    )


def test_deficit_two_argument_excludes_exactly_1_to_5_at_C_two_and_above() -> None:
    row = M.u6_energy_ledger_certificate()["partition_1^5"]
    assert row["independent_row_energy_upper"] == 357
    assert row["C_equals_2_exact_energy"] == 355
    assert row["deficit_from_independent_maxima"] == 2
    assert row["integer_square_parity_forces_each_row_deficit_even"] is True
    assert row["C_equals_2_excluded"] is True
    assert row["C_at_least_3_excluded_by_energy"] is True
    assert row["remaining_collision_counts"] == [0, 1]


def test_universal_four_root_character_product_kills_221_at_C_three() -> None:
    universal = M.quartic_character_product_certificate()
    assert universal["normalized_projective_point_count"] == 14
    assert universal["determinant_products"] == [1] * 14
    assert universal["determinant_product_sha256"] == (
        "e60bc0aa0c7b044e9130b8808d3d3035597fe37c8b923afb04811e59ed5ba7b7"
    )
    assert universal["general_formula"] == (
        "chi(c)^(p+1-d)*chi(-1)^(d*(d-1)/2)"
    )
    assert universal["direct_split_quartics_checked"] == 2002
    assert universal["every_four_root_nonroot_character_product"] == 1

    row = M.u6_energy_ledger_certificate()["partition_2_2_1"]
    assert row["independent_row_energy_upper"] == 371
    assert row["C_equals_3_exact_energy"] == 371
    assert row["forced_nonroot_character_product"] == -1
    assert row["C_equals_3_excluded"] is True
    assert row["C_at_least_4_excluded_by_energy"] is True
    assert row["remaining_collision_counts"] == [0, 1, 2]


def test_helper_reports_an_open_reduction_not_a_u6_close() -> None:
    row = M.u6_energy_ledger_certificate()
    assert row["result_status"] == "open reduction"
    assert row["sharp_row_upper_bounds_proved_here"] is False
    assert row["p13_t4_u6_closed"] is False
    assert row["proved"] is True


@pytest.mark.parametrize("kind", ["bad", "h1", "opposite"])
def test_unknown_row_kinds_fail_loudly(kind: str) -> None:
    with pytest.raises(ValueError):
        M.fixed_energy_row_catalog(kind, 28)
