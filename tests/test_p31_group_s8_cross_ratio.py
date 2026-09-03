from __future__ import annotations

import pytest

from e1_gmin_m4_p31_group_s8_cross_ratio import (
    anharmonic_orbit,
    fourth_silent_direction_reduction,
    p31_anharmonic_partition,
    quadruple_mitm_certificate_plan,
    theorem_record,
)


def test_anharmonic_orbits_partition_every_possible_fourth_direction():
    partition = p31_anharmonic_partition()
    assert partition == (
        (2, 16, 30),
        (3, 11, 15, 17, 21, 29),
        (4, 8, 10, 22, 24, 28),
        (5, 7, 9, 23, 25, 27),
        (6, 26),
        (12, 13, 14, 18, 19, 20),
    )
    assert set().union(*map(set, partition)) == set(range(2, 31))
    assert sum(map(len, partition)) == 29


def test_six_extra_silents_force_one_of_four_generic_cases():
    out = fourth_silent_direction_reduction()
    assert out["proved"] is True
    assert out["remaining_silent_directions"] == 6
    assert out["exceptional_union_size"] == 5
    assert out["harmonic_orbit"] == [2, 16, 30]
    assert out["equianharmonic_orbit"] == [6, 26]
    assert out["generic_orbit_representatives"] == [3, 4, 5, 12]
    assert out["large_search_run"] is False


def test_four_plus_four_plan_has_exact_coverage_and_storage_counts():
    out = quadruple_mitm_certificate_plan()
    assert out["proved"] is True
    assert out["fourth_direction_cases"] == [3, 4, 5, 12]
    assert out["signature"].startswith("the 60 parity bits")
    assert out["quadruples"] == 2_184_297_480
    assert out["quadruple_index_bits"] == 32
    assert out["unordered_4_plus_4_partitions_per_support"] == 35
    assert out["aligned_16_byte_record_gib"] < 33
    assert out["large_search_run"] is False


def test_record_does_not_promote_an_unrun_plan_to_a_certificate():
    out = theorem_record()
    assert out["proved"]["four_generic_cross_ratio_cases_are_complete"] is True
    assert out["proved"]["four_plus_four_signature_coverage"] is True
    assert out["proved"]["weight_eight_counterexample_excluded"] is False
    assert out["proved"]["group_support_lemma"] is False
    assert out["proved"]["row_code_minimum_distance"] is False
    assert out["proved"]["residual_ii_closed"] is False
    assert out["large_search_run"] is False


def test_cross_ratio_guard_rejects_a_member_of_the_fixed_triple():
    with pytest.raises(ValueError, match="fourth direction"):
        anharmonic_orbit(1)
