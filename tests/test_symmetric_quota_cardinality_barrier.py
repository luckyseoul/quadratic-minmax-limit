from e1_gmin_m4_symmetric_quota_cardinality_barrier import (
    branch_c_uniform_information_barrier,
    quota_cardinality_barrier,
)


def test_uniform_branch_c_symbolic_gap_at_first_prime():
    theorem = branch_c_uniform_information_barrier(31)
    assert theorem["N"] == 480
    assert theorem["branch_C_H_max"] == 479
    assert theorem["R"] == 7680
    assert theorem["R_minus_d"] == 7648
    assert theorem["slice_log2_strict_upper_bound"] == 7185
    assert theorem["exponent_gap"] == 463
    assert theorem["actual_target_status"] == "OPEN"
    assert theorem["proved"]


def test_exact_product_slice_is_smaller_than_fixed_parity_fibre():
    p = 31
    d = p + 1
    h = (p - 1) // 2
    full_capacity = d * h * h
    # A deliberately arbitrary feasible vector: the theorem is independent
    # of the actual branch-C target and performs no Radon-matrix replay.
    quotas = (8,) * 15 + (7,) * 17
    capacities = (full_capacity,) * d
    theorem = quota_cardinality_barrier(p, quotas, capacities)
    assert theorem["quota_total_s"] == 239
    assert theorem["exact_quota_slice_size"] < theorem[
        "fixed_parallel_parity_fibre_size"
    ]
    assert theorem["onto_plus_quota_bounds_implies_slice_nonempty_for_every_target"] is False
    assert theorem["actual_branch_C_target_excluded"] is False
    assert theorem["proved"]


def test_p43_symbolic_not_a_prime_census():
    theorem = branch_c_uniform_information_barrier(43)
    assert theorem["N"] == (43 * 43 - 1) // 2
    assert theorem["exponent_gap"] > 0
    assert theorem["proved"]


def test_zero_quota_singleton_is_strictly_smaller_too():
    p = 31
    d = p + 1
    theorem = quota_cardinality_barrier(p, (0,) * d, (0,) * d)
    assert theorem["quota_total_s"] == 0
    assert theorem["exact_quota_slice_size"] == 1
    assert theorem["exact_quota_slice_size"] < theorem[
        "fixed_parallel_parity_fibre_size"
    ]
    assert theorem["proved"]


def test_rejects_wrong_scope():
    try:
        branch_c_uniform_information_barrier(29)
    except ValueError as error:
        assert "branch-C prime" in str(error)
    else:
        raise AssertionError("p=29 should be outside the theorem scope")

    p = 31
    d = p + 1
    try:
        quota_cardinality_barrier(p, (240,) + (0,) * (d - 1), (240,) * d)
    except ValueError as error:
        assert "total quota" in str(error)
    else:
        raise AssertionError("a quota above (N-1)/2 should be rejected")
