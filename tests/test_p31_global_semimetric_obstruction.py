from fractions import Fraction

from e1_gmin_m4_p31_global_semimetric_obstruction import (
    CATEGORY_SIZES,
    DECOMPOSITION_WEIGHTS,
    FIXED_DIRECTION,
    K_ANTIPODAL_NEGATIVE,
    K_ANTIPODAL_POSITIVE,
    K_GENERIC_NEGATIVE,
    K_ZERO_NEGATIVE,
    K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE,
    T_ANTIPODAL,
    T_GENERIC,
    T_ZERO,
    T_ZERO_ANTIPODAL,
    atom_type_counts,
    expected_row_category_coefficients,
    global_fractional_semimetric_obstruction,
    replay_fractional_atom_decomposition,
)


def test_invariant_edge_and_atom_type_counts_are_exact():
    assert CATEGORY_SIZES == (30, 15, 420)
    rows = atom_type_counts()
    assert rows["compact"] == {
        K_GENERIC_NEGATIVE: 11_760,
        K_ANTIPODAL_NEGATIVE: 870,
        K_ZERO_NEGATIVE: 420,
        K_ANTIPODAL_POSITIVE: 420,
        K_ZERO_NEGATIVE_ANTIPODAL_POSITIVE: 15,
    }
    assert rows["positive"] == {
        T_GENERIC: 3_640,
        T_ZERO: 420,
        T_ANTIPODAL: 420,
        T_ZERO_ANTIPODAL: 15,
    }


def test_five_expected_row_types_have_exact_rational_coefficients():
    assert expected_row_category_coefficients(0) == (
        Fraction(-853, 14_400),
        Fraction(7, 160),
        Fraction(-373, 14_400),
    )
    assert expected_row_category_coefficients(FIXED_DIRECTION) == (
        Fraction(-791, 14_400),
        Fraction(-3, 160),
        Fraction(-311, 14_400),
    )
    assert expected_row_category_coefficients(2) == (
        Fraction(-137, 2_400),
        Fraction(11, 240),
        Fraction(-19, 800),
    )
    assert expected_row_category_coefficients(4) == (
        Fraction(31, 1_200),
        Fraction(-1, 24),
        Fraction(31, 1_200),
    )
    assert expected_row_category_coefficients(5) == (
        Fraction(403, 14_400),
        Fraction(-19, 480),
        Fraction(403, 14_400),
    )


def test_every_row_decomposition_matches_coefficients_and_atom_counts():
    classes = set()
    for direction in range(32):
        row = replay_fractional_atom_decomposition(direction)
        classes.add(row["row_class"])
        assert row["fractional_atom_cone_membership"]
        assert row["proved"]
    assert classes == set(DECOMPOSITION_WEIGHTS)


def test_uniform_fractional_graph_blocks_every_fixed_summed_semimetric_bank():
    row = global_fractional_semimetric_obstruction()
    assert row["graph_edge_count"] == 479
    assert row["fixed_edge_count"] == 1
    assert row["unused_double_orbit_count"] == 0
    assert row["nonfixed_edge_count"] == 478
    assert row["fractional_graph_is_convex_combination_of_ledger_graphs"]
    assert row["all_expected_rows_in_fractional_atom_cones"]
    assert not row["fixed_sum_of_row_semimetrics_separates_every_ledger_graph"]
    assert not row["adaptive_semimetric_oracle_excluded"]
    assert not row["common_graph_constructed"]
    assert not row["residual_ii_closed"]
    assert row["proved"]
