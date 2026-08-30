import pytest

from e1_gmin_m4_prop15726 import (
    SMALL_PRIMES,
    all_deletion_sizes_excluded,
    deletion_size_incidence_contradiction,
    linewise_slack_incidence_lemma,
    linear_tangent_envelope_cutoff,
    outside_linear_slack_exclusion,
    proposition_15726,
    small_prime_cutoff_table,
    tangent_envelope_dependency,
    universal_symbolic_ledger,
)


def test_exact_linear_cutoffs_close_all_former_R_four_cells():
    assert small_prime_cutoff_table() == {
        17: 4,
        19: 5,
        23: 6,
        29: 8,
        31: 9,
        37: 11,
        41: 12,
    }
    assert all(linear_tangent_envelope_cutoff(p) >= 4 for p in SMALL_PRIMES)


def test_line_slack_dominates_deleted_point_secant_incidence_universally():
    for u in range(257):
        row = linewise_slack_incidence_lemma(u)
        assert row["slack_contribution"] >= u
        assert row["gap"] >= 0
        assert row["formula_matches"] is True
        assert row["proved"] is True
    assert "a(a-1)" in linewise_slack_incidence_lemma(20)["gap_formula"]
    assert "a^2" in linewise_slack_incidence_lemma(21)["gap_formula"]


def test_every_minimal_deletion_size_has_the_strict_envelope_contradiction():
    for p in SMALL_PRIMES:
        cutoff = linear_tangent_envelope_cutoff(p)
        for R in range(1, cutoff + 1):
            whole_R = all_deletion_sizes_excluded(p, R)
            assert whole_R["possible_minimal_deletion_sizes"] == list(
                range(1, R + 1)
            )
            assert whole_R["every_deletion_size_contradictory"] is True
            for t, row in enumerate(whole_R["rows"], start=1):
                assert row == deletion_size_incidence_contradiction(p, R, t)
                assert row["stronger_strict_size_bound_met"] is True
                assert row["chosen_point_secant_index_lower"] == 1
                assert row["chosen_point_secant_index_upper"] >= 1
                assert row["tangent_count_minus_envelope_degree"] > 0
                assert row["dual_line_component_forced"] is True
                assert row["secant_evaluation_nonzero"] is True
                assert row["contradiction"] is True


def test_current_and_arxiv_tangent_envelope_numbering_is_explicit():
    theorem = tangent_envelope_dependency()
    assert theorem["current_manuscript_theorem"] == 13
    assert theorem["arxiv_v4_theorem"] == 11
    assert theorem["arxiv"] == "1705.10940v4"
    assert theorem["size_hypothesis_is_weak"] is True
    assert theorem["proof_meets_stronger_strict_size_bound"] is True
    assert "|A|>=2*tau+2" in theorem["odd_order_statement"]
    assert "degree-2*tau" in theorem["odd_order_statement"]
    assert theorem["proved"] is True


def test_symbolic_ledger_contains_the_complete_noncomputational_proof_chain():
    row = universal_symbolic_ledger()
    assert row["result_status"] == "proved theorem"
    assert "I=sum_(z in T) s_z<=R" == row["linewise_incidence_bound"]["conclusion"]
    assert "(R-t)(3-2/t)>=0" in row["algebraic_identity"]
    assert "p-2<p-1" in row["strict_tangent_chain"]
    assert "z*|Phi" in row["dual_root_count"]
    assert "f_P(z)^2!=0" in row["secant_evaluation"]
    assert row["finite_profile_search_used"] is False
    assert row["proved"] is True


def test_prime_parameter_validation_rejects_small_primes_and_composites():
    for invalid in (2, 3, 15, 16, 21, 25, 27, 35, 49):
        with pytest.raises(ValueError):
            linear_tangent_envelope_cutoff(invalid)
    with pytest.raises(ValueError):
        linewise_slack_incidence_lemma(-1)
    with pytest.raises(ValueError):
        all_deletion_sizes_excluded(17, 5)
    with pytest.raises(ValueError):
        deletion_size_incidence_contradiction(17, 4, 5)


def test_theorem_closes_the_linear_interval_but_not_the_shell_or_global_gates():
    for p in SMALL_PRIMES:
        row = outside_linear_slack_exclusion(p)
        cutoff = (p - 4) // 3
        assert row["excluded_positive_R_interval"] == [1, cutoff]
        assert row["first_possible_positive_R_at_least"] == cutoff + 1
        assert row["audited_integer_case_count"] == cutoff * (cutoff + 1) // 2
        assert row["all_integer_cases_verified"] is True
        assert row["minimum_strict_size_margin"] > 0
        assert row["minimum_tangent_degree_gap"] > 0
        assert row["former_R_four_cell_closed"] is True
        assert row["finite_profile_search_used"] is False
        assert row["p_plus_one_shell_closed"] is False
        assert row["non_walsh_residual_ii_closed"] is False
        assert row["multi_level_type_I_closed"] is False
        assert row["quadratic_minmax_limit_closed"] is False
        assert row["top_level_gates_changed"] is False
        assert row["result_status"] == "proved theorem"
        assert row["proved"] is True


def test_proposition_package_is_honest_about_the_remaining_frontier():
    result = proposition_15726()
    assert result["prop"] == "15.726"
    assert "floor((p-4)/3)" in result["statement"]
    assert result["R_four_closed_for_all_primes_p_ge_17"] is True
    assert "larger outside slack" in result["remaining_scope"]
    assert "Type I" in result["remaining_scope"]
    assert "L remain open" in result["remaining_scope"]
    assert result["top_level_gates_changed"] is False
    assert result["result_status"] == "proved theorem"
    assert result["proved"] is True
