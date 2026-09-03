from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_slice_exchange import (
    exact_small_slice_replay,
    symmetric_slice_kernel_theorem,
    theorem_record,
    used_orbit_deletion_bound,
)


def test_symbolic_one_slice_kernel_is_a_root_lattice():
    base = symmetric_slice_kernel_theorem(3)
    assert base["kernel_lattice"] == "A_0"
    assert base["kernel_rank"] == 0
    assert base["proved"] is True

    out = symmetric_slice_kernel_theorem(31)
    assert out["h"] == 15
    assert out["nonfixed_midpoint_orbits_per_difference_slice"] == 480
    assert out["zero_label_slab_size"] == 15
    assert out["nonzero_square_label_slab_count"] == 15
    assert out["nonzero_square_label_slab_size"] == 31
    assert out["kernel_lattice"] == "A_14"
    assert out["kernel_rank"] == 14
    assert out["circuit_positive_degree_in_pair_variables"] == 31
    assert out["circuit_support_in_pair_variables"] == 62
    assert out["physical_graph_edges_removed_and_added"] == 62
    assert out["weight_preserving"] is True
    assert out["all_directional_parallel_coordinates_preserved"] is True
    assert out["unsigned_kernel_valid_for_all_odd_primes"] is True
    assert out["paley_signed_specialization"].startswith(
        "only for p=3 mod 4"
    )
    assert out["full_unused_configuration_normality_proved"] is False
    assert out["proved"] is True


def test_mobius_deletion_cap_leaves_many_clean_slabs_in_one_slice():
    out = used_orbit_deletion_bound(31)
    assert out["difference_slice_count"] == 480
    assert out["mobius_used_orbit_cap"] == 480
    assert out["one_slice_used_orbits_at_most"] == 1
    assert out["clean_nonzero_slabs_in_that_slice_at_least"] == 14
    assert out["whole_slab_circuits_in_that_slice_at_least"] == 91
    assert out["surviving_unused_whole_slab_circuit_for_p_at_least_7"] is True
    assert out["global_connectivity_proved"] is False
    assert out["proved"] is True


def test_exact_p5_p7_replay_has_only_the_predicted_slice_kernel():
    expected = {5: (12, 11, 1), 7: (24, 22, 2)}
    for p, (columns, rank, nullity) in expected.items():
        out = exact_small_slice_replay(p)
        assert out["midpoint_orbit_columns"] == columns
        assert out["map_replayed"] == "unsigned reduced pair columns"
        assert out["exact_rational_column_rank"] == rank
        assert out["exact_integer_kernel_rank"] == nullity
        assert out["expected_integer_kernel_rank"] == nullity
        assert out["all_nonzero_slab_column_sums_equal"] is True
        assert out["proved"] is True


def test_theorem_record_stops_at_slice_connectivity():
    out = theorem_record(31)
    assert out["proved"] is True
    assert out["status"] == "PROVED SLICE KERNEL AND CONNECTIVITY; GLOBAL BOX OPEN"
    assert out["full_unused_configuration_normality_proved"] is False
    assert out["global_Boolean_fibre_nonempty_proved"] is False
    assert out["residual_ii_closed"] is False


def test_parameter_guards():
    with pytest.raises(ValueError):
        symmetric_slice_kernel_theorem(9)
    with pytest.raises(ValueError):
        exact_small_slice_replay(11)
    with pytest.raises(ValueError):
        used_orbit_deletion_bound(31, 481)
