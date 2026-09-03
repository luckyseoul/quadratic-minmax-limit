import pytest

from e1_gmin_m4_prop15758 import (
    p1_local_survivor,
    p1_lower_endpoint_parseval,
    p3_local_survivor,
    p3_lower_endpoint_parseval,
    sharp_atom_l1_certificate,
    theorem_record,
)


def test_sharp_fixed_offset_l1_minimum_and_constant_cancellation():
    for a in range(7):
        for b in range(13):
            row = sharp_atom_l1_certificate(29, a, b)
            assert row["proved"]
            forced_sum = 3 * a - b
            expected = 3 if abs(forced_sum) == 1 else abs(forced_sum)
            assert row["forced_coefficient_sum"] == forced_sum
            assert row["exact_minimum_l1"] == expected

    cancelled = sharp_atom_l1_certificate(29, 4, 12)
    assert cancelled["scaled_mass"] == 16 * 29
    assert cancelled["exact_minimum_l1"] == 0
    assert cancelled["attaining_values_of_4B"] == [16]


def test_both_prime_classes_have_full_m_direction_local_survivors():
    for p in (17, 29, 37):
        r = (p - 1) // 4
        bounds = (2 * r * r - 5 * r, 4 * r * r - 6 * r - 3)
        for t in (bounds[0], sum(bounds) // 2, bounds[1]):
            row = p1_local_survivor(p, t)
            assert row["proved_local_aggregate"]
            assert len(row["hard_rows"]) == row["m"]
            assert len(row["opposite_rows"]) == row["m"]
            assert row["transported_infinity_degree_I"] == 0
            assert row["direction_count_per_type"] == row["m"]
            assert not row["missing_projective_direction"]
            assert row["hard_parallel_edge_total"] + row["opposite_parallel_edge_total"] == row["H_edge_count"]
            assert not row["one_common_simple_graph_constructed"]

    for p in (19, 31, 43):
        r = (p - 3) // 4
        bounds = (2 * r * r - 4 * r - 2, 4 * r * r - 2 * r - 5)
        for t in (bounds[0], sum(bounds) // 2, bounds[1]):
            row = p3_local_survivor(p, t)
            assert row["proved_local_aggregate"]
            assert len(row["hard_rows"]) == row["m"]
            assert len(row["opposite_rows"]) == row["m"]
            assert row["transported_infinity_degree_I"] == 0
            assert row["direction_count_per_type"] == row["m"]
            assert not row["missing_projective_direction"]
            assert row["hard_parallel_edge_total"] + row["opposite_parallel_edge_total"] == row["H_edge_count"]
            assert not row["one_common_simple_graph_constructed"]


def test_survivor_interval_endpoints_are_exact_isolated_chart_bounds():
    p = 37
    r = (p - 1) // 4
    lower = 2 * r * r - 5 * r
    upper = 4 * r * r - 6 * r - 3
    assert p1_local_survivor(p, upper)["isolated_vertex_gap"] == 4
    with pytest.raises(ValueError):
        p1_local_survivor(p, lower - 1)
    with pytest.raises(ValueError):
        p1_local_survivor(p, upper + 1)

    p = 31
    r = (p - 3) // 4
    lower = 2 * r * r - 4 * r - 2
    upper = 4 * r * r - 2 * r - 5
    assert p3_local_survivor(p, upper)["isolated_vertex_gap"] == 4
    with pytest.raises(ValueError):
        p3_local_survivor(p, lower - 1)
    with pytest.raises(ValueError):
        p3_local_survivor(p, upper + 1)


def test_scalar_parseval_has_no_upper_lower_gap_on_the_two_rays():
    for r in (7, 9, 12):
        row = p1_lower_endpoint_parseval(r)
        assert row["proved"]
        assert row["compact_max_minus_global_C0"] > 0
        assert not row["scalar_Parseval_upper_lower_gap_closes"]

        row = p3_lower_endpoint_parseval(r)
        assert row["proved"]
        assert row["compact_max_minus_global_C0"] > 0
        assert not row["scalar_Parseval_upper_lower_gap_closes"]


def test_theorem_record_preserves_the_common_graph_gate():
    out = theorem_record()
    assert out["prop"] == "15.758"
    assert out["proved"]["sharp_atom_l1_minimum"]
    assert out["proved"]["p_1_mod_4_local_survivor_interval"]
    assert out["proved"]["p_3_mod_4_local_survivor_interval"]
    assert not out["proved"]["scalar_Parseval_alone_excludes_survivors"]
    assert not out["proved"]["one_common_simple_graph_constructed"]
    assert not out["proved"]["residual_ii_closed"]
    assert out["L_status"] == "OPEN"
