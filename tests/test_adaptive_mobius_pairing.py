import pytest

from e1_gmin_m4_adaptive_mobius_pairing import (
    adaptive_center_pairing_certificate,
    color_count_theorem,
    complementary_pair_parameters,
    endpoint_color_profile,
    exact_magnitude_scope_guard,
    forced_affine_auxiliary_pair,
    paired_square_invariant_counterexample,
    prescribed_auxiliary_assignment_criterion,
    theorem_record,
)


@pytest.mark.parametrize(
    ("p", "m", "flexible", "mono_plus", "mono_minus"),
    ((31, 16, 14, 8, 8), (43, 22, 20, 10, 12)),
)
def test_endpoint_color_counts_include_both_mod_eight_classes(
    p: int, m: int, flexible: int, mono_plus: int, mono_minus: int
) -> None:
    result = color_count_theorem(p)
    assert result["m"] == m
    assert result["flexible_nonzero_evaluations"] == flexible
    assert result["monochrome_plus_evaluations"] == mono_plus
    assert result["monochrome_minus_evaluations"] == mono_minus
    assert endpoint_color_profile(31, 2)["valid_signs"] == (-1,)
    assert endpoint_color_profile(31, -2)["valid_signs"] == (1,)


def test_exact_pair_formulas_solve_both_singleton_equations() -> None:
    p = 31
    for a in range(1, p):
        first = endpoint_color_profile(p, a)
        for b in range(1, p):
            second = endpoint_color_profile(p, b)
            for eps, first_color in first["sign_colors"].items():
                for eps_prime, second_color in second["sign_colors"].items():
                    if first_color * second_color != -1:
                        continue
                    result = complementary_pair_parameters(
                        p, a, b, eps, eps_prime
                    )
                    assert result["first_auxiliary_evaluation"] == 2 * eps % p
                    assert result["second_auxiliary_evaluation"] == 2 * eps_prime % p
                    assert result["nu_character"] == -1


def test_arbitrary_nonzero_centers_get_a_target_perfect_matching() -> None:
    alphas = (1, 7, 4, 22, 9, 3, 29, 8, 14, 2, 25, 11, 6, 19, 13, 27)
    result = adaptive_center_pairing_certificate(31, alphas)
    assert result["target_perfect_matching_proved"]
    assert result["chosen_scale_monochrome_total"] <= 8
    assert result["sum_over_scales_of_monochrome_vertices"] == 16**2
    assert len(result["matched_target_pairs"]) == 8
    assert all(
        row["parameters"]["complementary_nonsquare"]
        for row in result["matched_target_pairs"]
    )
    assert not result["auxiliary_direction_SDR_proved"]
    assert not result["auxiliary_type_quota_proved"]


def test_forced_auxiliary_map_is_a_two_cycle_and_only_local() -> None:
    alphas = tuple(range(1, 17))
    matching = adaptive_center_pairing_certificate(31, alphas)
    row = matching["matched_target_pairs"][0]
    i = row["first_target"]
    k = row["second_target"]
    result = forced_affine_auxiliary_pair(
        31,
        i,
        k,
        alphas[i],
        alphas[k],
        matching["chosen_fixed_edge_scale_c"],
        row["first_sign"],
        row["second_sign"],
    )
    assert result["within_pair_auxiliaries_distinct"]
    assert result["auxiliaries_avoid_both_paired_targets"]
    assert result["nu_character"] == -1
    assert not result["cross_pair_auxiliary_distinctness_proved"]
    assert not result["auxiliary_type_quota_proved"]


def test_target_matching_can_have_cross_pair_auxiliary_collisions() -> None:
    """This fixed replay is a limitation witness, not a search/no-go claim."""
    p = 31
    alphas = tuple(range(1, 17))
    targets = tuple(range(16))
    matching = adaptive_center_pairing_certificate(p, alphas)
    auxiliaries = [0] * 16
    tau = [0] * 16
    signs = [0] * 16
    for row in matching["matched_target_pairs"]:
        i = row["first_target"]
        k = row["second_target"]
        local = forced_affine_auxiliary_pair(
            p,
            targets[i],
            targets[k],
            alphas[i],
            alphas[k],
            matching["chosen_fixed_edge_scale_c"],
            row["first_sign"],
            row["second_sign"],
        )
        auxiliaries[i] = local["first_auxiliary_coordinate_U"]
        auxiliaries[k] = local["second_auxiliary_coordinate_V"]
        tau[i] = k
        tau[k] = i
        signs[i] = row["first_sign"]
        signs[k] = row["second_sign"]

    criterion = prescribed_auxiliary_assignment_criterion(
        p,
        targets,
        alphas,
        auxiliaries,
        tau,
        tuple(range(16)),
        matching["chosen_fixed_edge_scale_c"],
        signs,
    )
    assert criterion["all_endpoint_equations_hold"]
    assert criterion["all_pair_parameters_nonsquare"]
    assert criterion["all_pair_square_invariants_hold"]
    assert len(set(auxiliaries)) == 9
    assert not criterion["auxiliary_occurrences_are_distinct_directions"]
    assert not criterion["pair_coherent_distinct_auxiliary_assignment"]
    assert not criterion["paley_type_quota_checked"]
    assert not criterion["full_endpoint_assignment_proved"]


def test_all_denominator_exceptions_are_explicit() -> None:
    with pytest.raises(ValueError, match=r"a=2\*eps"):
        complementary_pair_parameters(31, 2, 3, 1, 1)
    with pytest.raises(ValueError, match=r"b=2\*eps_prime"):
        complementary_pair_parameters(31, 3, 2, 1, 1)

    # w_i=2, w_k=3, c=5 makes nu=1 and the affine two-cycle singular.
    with pytest.raises(ValueError, match="nu=1"):
        forced_affine_auxiliary_pair(31, 0, 1, 1, 11, 5, 1, 1)


def test_bad_exact_magnitude_is_not_promoted_to_direction_obstruction() -> None:
    result = exact_magnitude_scope_guard(31)
    assert result["pair_nu_is_square"]
    assert not result["pairing_at_this_exact_magnitude"]
    assert not result["whole_fixed_direction_obstructed"]
    assert result["role"] == "scope guard, not an obstruction"


def test_paired_square_multiset_does_not_force_one_chord_ratio() -> None:
    result = paired_square_invariant_counterexample(31)
    assert result["square_fibre_polynomial"] == "(T-1)^4"
    assert result["perfect_matching_squared_chords"] == (
        (1, 4),
        (4, 9),
        (16, 1),
    )
    assert not result["one_common_chord_ratio_exists"]
    assert not result["paired_square_polynomial_condition_sufficient"]
    assert result["scope"] == "abstract invariant counterexample, not a branch-C target"


def test_theorem_record_preserves_exact_open_frontier() -> None:
    result = theorem_record()
    assert result["proved_all_claimed_statements"]
    assert result["target_perfect_matching_proved"]
    assert not result["cross_pair_auxiliary_distinctness_proved"]
    assert not result["auxiliary_type_quota_proved"]
    assert not result["full_parallel_quota_identity_constructed"]
    assert not result["fixed_word_singleton_constructed"]
    assert not result["common_graph_constructed"]
    assert not result["residual_ii_closed"]
