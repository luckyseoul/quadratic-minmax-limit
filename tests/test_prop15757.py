import pytest

from e1_gmin_m4_prop15757 import (
    compact_survivor_parity_ledger,
    compact_triangle_atom,
    edge_radon_dimensions,
    exact_edge_radon_rank,
    sharp_atom_masses,
    theorem_record,
)


def test_binary_edge_radon_requires_an_odd_prime():
    for p in (2, 4, 9):
        with pytest.raises(ValueError, match="odd prime"):
            edge_radon_dimensions(p)


def test_exact_binary_edge_radon_ranks():
    for p in (3, 5, 7, 11):
        exact = exact_edge_radon_rank(p)
        dimensions = edge_radon_dimensions(p)
        assert exact["proved"]
        assert exact["exact_rank_over_F2"] == dimensions["image_rank"]
        assert dimensions["compatibility_codimension"] == p + 1


def test_compact_signed_triangle_is_the_p_plus_one_atom():
    for p in (5, 13, 17, 29):
        atom = compact_triangle_atom(p)
        assert atom["proved"]
        assert atom["scaled_mass_4pE"] == p + 1
        assert atom["coefficient_l1"] == 3
        assert sharp_atom_masses(p)["scaled_mass_4pE"] == p - 3


def test_compact_aggregate_survivor_passes_complete_f2_image_gate():
    for p, t in ((13, 5), (17, 20), (29, 70), (37, 117)):
        row = compact_survivor_parity_ledger(p, t)
        assert row["proved"]
        assert row["hard_direction_count"] == row["opposite_direction_count"]
        assert row["all_affine_direction_count"] == p + 1
        assert row["isolated_chart_infinity_parallel_count_I"] == 0
        assert row["all_direction_parallel_counts_sum_to_H_edges"]
        assert row["binary_edge_lift_exists_by_exact_image_theorem"]
        assert not row["integer_nonnegative_simple_edge_lift_proved"]


def test_theorem_record_keeps_global_gate_open():
    out = theorem_record()
    assert out["proved"]["boundary_and_total_parity_conditions_are_complete"]
    assert out["proved"]["aggregate_survivor_has_binary_edge_lift"]
    assert not out["proved"]["integer_nonnegative_simple_edge_lift"]
    assert not out["proved"]["residual_ii_closed"]
    assert out["L_status"] == "OPEN"
