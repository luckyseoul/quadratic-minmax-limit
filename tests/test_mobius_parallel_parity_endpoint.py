from __future__ import annotations

import pytest

from e1_gmin_m4_mobius_parallel_parity_endpoint import (
    balanced_parallel_parity_profile,
    cancellation_offset_consequence,
    theorem_record,
)
from e1_gmin_m4_prop15758 import p3_local_survivor


def test_exact_quotas_are_the_canonical_balanced_branch_c_quotas():
    for t in (68, 79, 80, 98, 99, 177):
        out = balanced_parallel_parity_profile(31, t)
        canonical = p3_local_survivor(31, t)
        assert out["hard_parallel_quotas"] == [
            row["parallel_P"] for row in canonical["hard_rows"]
        ]
        assert out["opposite_parallel_quotas"] == [
            row["Q"] for row in canonical["opposite_rows"]
        ]
        assert sum(out["hard_compact_counts"]) == t + 1
        assert out["opposite_quota_total"] == t + 76
        assert out["opposite_minus_hard_total"] == 75
        assert out["proved"] is True


def test_closed_residue_formula_matches_direct_balanced_parities():
    # This is an algebra replay over one full period, not a prime census or a
    # source-support enumeration.
    for p, start in ((31, 68), (43, 158)):
        m = (p + 1) // 2
        for residue in range(2 * m):
            t = start + ((residue - (start + 1)) % (2 * m))
            out = balanced_parallel_parity_profile(p, t)
            assert out["residue_s_equals_E_mod_p_plus_1"] == residue
            assert out["base_direction_parity_weight_w0"] == sum(
                out["base_direction_parity"]
            )


def test_j0_and_j1_endpoint_bands_are_exact():
    # p=31 has m=16 and modulus p+1=32.
    fully_excluded = cancellation_offset_consequence(31, 68, 1)  # s=5
    assert fully_excluded["base_direction_parity_weight_w0"] == 21
    assert fully_excluded["minimum_forced_fixed_edge_weight"] == 5
    assert fully_excluded["minimum_extra_cancellations_beyond_size_floor"] == 2
    assert fully_excluded["excluded_by_parallel_parity"] is True

    fixed_only = cancellation_offset_consequence(31, 99, 1)  # s=4
    assert fixed_only["base_direction_parity_weight_w0"] == 19
    assert fixed_only["minimum_forced_fixed_edge_weight"] == 3
    assert fixed_only["minimum_extra_cancellations_beyond_size_floor"] == 1
    assert fixed_only["excluded_by_parallel_parity"] is False
    assert fixed_only["all_remaining_capacity_forced_to_fixed_edges_at_bound"] is True
    assert fixed_only["conditional_unused_double_orbits_at_bound"] == 0

    j0_boundary = cancellation_offset_consequence(31, 98, 0)  # s=3
    assert j0_boundary["base_direction_parity_weight_w0"] == 17
    assert j0_boundary["minimum_forced_fixed_edge_weight"] == 1
    assert j0_boundary["excluded_by_parallel_parity"] is False
    assert j0_boundary["symmetric_Boolean_completion_constructed"] is False

    upper_endpoint = cancellation_offset_consequence(31, 177, 0)  # s=18
    assert upper_endpoint["base_direction_parity_weight_w0"] == 17
    assert upper_endpoint["excluded_by_parallel_parity"] is False


def test_general_offset_bound_strengthens_the_support_floor_by_zero_one_or_two():
    cases = (
        (68, 2),   # s=5: w0=m+5
        (99, 1),   # s=4: w0=m+3
        (98, 0),   # s=3: w0=m+1
        (97, 0),   # s=2: w0=m-1
        (96, 0),   # s=1: w0=m-3
        (95, 0),   # s=0: w0=m-5
    )
    for t, expected_extra in cases:
        out = balanced_parallel_parity_profile(31, t)
        assert out["minimum_extra_cancellations_beyond_size_floor"] == expected_extra
        assert out["strengthened_cancellation_lower_bound"] == (
            out["size_floor_cancellations"] + expected_extra
        )


def test_theorem_record_keeps_the_integral_target_open():
    out = theorem_record(31, 68)
    assert out["proved"] is True
    assert out["j0_excluded_residues"] == "4<=s<=m+1"
    assert out["j1_fully_excluded_residues"] == "5<=s<=m"
    assert out["integral_transverse_fibre_solved"] is False
    assert out["residual_ii_closed"] is False


def test_parameter_validation():
    with pytest.raises(ValueError, match="p=3 mod 4"):
        balanced_parallel_parity_profile(29, 100)
    with pytest.raises(ValueError, match="68<=t<=177"):
        balanced_parallel_parity_profile(31, 67)
    with pytest.raises(ValueError, match="nonnegative"):
        cancellation_offset_consequence(31, 68, -1)
