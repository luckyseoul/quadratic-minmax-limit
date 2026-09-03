from fractions import Fraction

import pytest

from e1_gmin_m4_mobius_endpoint_barrier import (
    actual_four_candidate_replay,
    centered_compact_atom_support,
    center_invisibility_theorem,
    fixed_family_center_coherence,
    opposite_swapped_locus,
    rational_clean_overlap_example,
    star_center_invisibility_replay,
    theorem_record,
)


def test_generic_opposite_swapped_locus_is_clean() -> None:
    result = opposite_swapped_locus(2, 3)
    assert result["A"] == Fraction(5, 6)
    assert result["B"] == Fraction(5, 3)
    assert result["candidate_verdicts"] == {
        "opposite_swapped": True,
        "same_orientation_swapped": False,
        "same_orientation_direct": False,
        "opposite_direct": False,
    }
    assert result["shared_orbit_count"] == 1
    assert result["desired_swapped_orbit_is_unique"]


@pytest.mark.parametrize(
    ("q", "r"),
    ((2, 3), (2, 2), (Fraction(3, 2), Fraction(5, 2))),
)
def test_all_four_actual_edge_candidates_are_compared(q: object, r: object) -> None:
    result = actual_four_candidate_replay(q, r)  # type: ignore[arg-type]
    assert result["candidate_verdicts"] == {
        "same_orientation_direct": False,
        "same_orientation_swapped": False,
        "opposite_direct": False,
        "opposite_swapped": True,
    }
    assert all("accepted" in row for row in result["candidates"].values())


def test_half_half_is_the_unique_double_overlap_point() -> None:
    result = opposite_swapped_locus(Fraction(1, 2), Fraction(1, 2))
    assert result["A"] == result["B"] == Fraction(3, 2)
    assert result["candidate_verdicts"]["opposite_direct"]
    assert result["shared_orbit_count"] == 2
    assert not result["desired_swapped_orbit_is_unique"]
    assert result["two_trade_sum_is_ternary"]


@pytest.mark.parametrize(
    ("q", "r"),
    ((0, 2), (2, 0), (1, 2), (2, 1), (2, Fraction(1, 2))),
)
def test_locus_rejects_degenerate_parameters(q: object, r: object) -> None:
    with pytest.raises(ValueError):
        opposite_swapped_locus(q, r)  # type: ignore[arg-type]


def test_q_equals_r_equals_two_rational_example() -> None:
    result = rational_clean_overlap_example()
    assert result["A"] == result["B"] == Fraction(3, 4)
    assert result["t"] == result["s"] == Fraction(-1, 2)
    assert result["opposite_edges"]
    assert result["M1_at_singleton"] == result["M2_at_singleton"] == 2
    assert result["other_three_candidates_rejected"]
    assert not result["global_complementary_family_constructed"]


def test_unit_star_moments_are_directly_center_independent_through_top_degree() -> None:
    replay = star_center_invisibility_replay(31, (1, 2, 30))
    assert replay["proved"]
    assert replay["top_degree_k0_value"] == 30
    expected = replay["expected_moments"]
    assert expected[30, 0] == 30
    assert all(value == 0 for key, value in expected.items() if key != (30, 0))
    assert len(set(tuple(row.items()) for row in replay["moments_by_center"].values())) == 1


def test_fixed_family_coherence_is_explicit_and_adaptive_choice_stays_open() -> None:
    coherent = fixed_family_center_coherence(
        (2, 3, 5), (1, Fraction(3, 2), Fraction(5, 2))
    )
    incoherent = fixed_family_center_coherence(
        (2, 3, 5), (2, Fraction(3, 2), Fraction(5, 2))
    )
    assert coherent["fixed_preassigned_family_coherent"]
    assert coherent["common_ratio"] == Fraction(1, 4)
    assert not incoherent["fixed_preassigned_family_coherent"]
    assert incoherent["common_ratio"] is None
    assert not coherent["adaptive_auxiliary_choice_analyzed"]


def test_centered_compact_atom_has_one_fixed_and_one_paired_source_orbit() -> None:
    result = centered_compact_atom_support(31)
    assert result["fixed_singleton_edge"] == (1, 30)
    assert result["fixed_singleton_support_size"] == 1
    assert result["paired_group_support_size"] == 1
    assert result["quotient_support_size"] == 2
    assert not any(result["canonical_odd_compact_moments"].values())

    # For v=(1,0), the p finite projective functionals (1,s) are nonzero
    # on v and the sole infinity functional (0,1) annihilates it. Hence the
    # singleton fixed word has exactly one silent affine-block group.
    v = (1, 0)
    directions = tuple((1, slope) for slope in range(31)) + ((0, 1),)
    projected = tuple((a * v[0] + b * v[1]) % 31 for a, b in directions)
    assert sum(value == 0 for value in projected) == 1
    assert sum(value != 0 for value in projected) == 31


def test_center_invisibility_scope_and_open_guard() -> None:
    result = center_invisibility_theorem()
    assert not result["unit_star_moments"]["depends_on_nonzero_center_j"]
    assert result["unit_star_moments"]["d=p-1 and k=0"] == "-1"
    assert not result["preassigned_fixed_auxiliary_family_automatically_coherent"]
    assert result["adaptive_auxiliary_choice_can_restore_coherence"] == "OPEN"
    assert not result["global_complementary_family_refuted"]
    assert "lambda^2!=1" in result["fixed_family_incoherence_witness"]


def test_theorem_record_keeps_global_problem_open() -> None:
    result = theorem_record()
    assert result["proved_all_claimed_statements"]
    assert not result["adaptive_global_complementary_choice_resolved"]
    assert not result["global_complementary_family_constructed"]
    assert not result["global_complementary_family_refuted"]
    assert not result["coherent_target_Boolean_completion_constructed"]
    assert not result["all_admissible_targets_excluded"]
    assert not result["residual_ii_closed"]
