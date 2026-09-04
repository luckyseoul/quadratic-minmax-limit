from e1_gmin_m4_p31_row_atom_incremental_bound import (
    AtomRowSpec,
    IncrementalAtomRowBound,
    frozen_graph_incremental_bound_certificate,
    strict_over_scalar_l1_witness,
    theorem_record,
)


def test_atom_row_spec_has_exact_occurrence_and_cut_budgets() -> None:
    hard = AtomRowSpec.hard(11)
    assert hard.target_edge_sum == -11
    assert hard.positive_occurrence_budget == 11
    assert hard.negative_occurrence_budget == 22
    assert hard.l1_budget == 33
    assert hard.cut_interval == (-22, 0)

    opposite = AtomRowSpec.opposite(6)
    assert opposite.target_edge_sum == 12
    assert opposite.positive_occurrence_budget == 24
    assert opposite.negative_occurrence_budget == 12
    assert opposite.l1_budget == 36
    assert opposite.cut_interval == (-12, 12)


def test_actual_atom_sum_has_zero_necessary_cost() -> None:
    # Six copies of +01+02+12 and six copies of K(3,4;5).
    coefficients = {
        (0, 1): 6,
        (0, 2): 6,
        (1, 2): 6,
        (3, 4): 6,
        (3, 5): -6,
        (4, 5): -6,
    }
    state = IncrementalAtomRowBound(coefficients, AtomRowSpec.opposite(6))
    row = state.summary()
    assert row["budget_defects"] == {
        "edge_sum": 0,
        "positive_mass": 0,
        "negative_mass": 0,
        "l1": 0,
    }
    assert row["odd_degree_count"] == 0
    assert row["degree_projection"]["feasible"] is True
    assert row["singleton_cut_violations"] == 0
    assert row["two_label_cut_violations"] == 0
    assert row["incremental_search_cost"] == 0
    assert row["coefficient_l1_edit_lower_bound"] == 0


def test_two_label_cut_is_strictly_stronger_than_scalar_l1_and_degrees() -> None:
    row = strict_over_scalar_l1_witness()
    assert row["proved"] is True
    assert row["compact_atoms"] == 11
    assert row["scalar_l1_and_mass_budgets_pass"] is True
    assert row["signed_degree_projection_passes"] is True
    assert row["violated_subset"] == [0, 1]
    assert row["violated_cut_value"] == 2
    assert row["required_cut_interval"] == [-22, 0]
    assert row["coefficient_l1_edit_lower_bound"] == 2


def test_cell_delta_update_exactly_matches_full_rebuild() -> None:
    spec = AtomRowSpec.hard(11)
    coefficients = {
        (0, 1): -3,
        (0, 2): 2,
        (3, 4): -1,
        (8, 13): 4,
        (13, 21): -2,
    }
    state = IncrementalAtomRowBound(coefficients, spec)
    for first, second, delta in (
        (0, 1, 2),
        (2, 17, -3),
        (8, 13, -5),
        (4, 30, 1),
    ):
        edge = tuple(sorted((first, second)))
        coefficients[edge] = coefficients.get(edge, 0) + delta
        if coefficients[edge] == 0:
            del coefficients[edge]
        state.apply_cell_delta(first, second, delta)
        rebuilt = IncrementalAtomRowBound(coefficients, spec)
        assert state.coefficients == rebuilt.coefficients
        assert state.degrees == rebuilt.degrees
        assert state.pair_cuts == rebuilt.pair_cuts
        assert state.summary() == rebuilt.summary()

    row = state.summary()
    assert row["cell_delta_affected_singletons"] == 2
    assert row["cell_delta_affected_two_label_cuts"] == 58
    assert row["cell_delta_update_complexity"] == "O(p)"


def test_exact_479_edge_api_is_scored_without_claiming_decomposition() -> None:
    row = frozen_graph_incremental_bound_certificate()
    assert row["proved"] is True
    assert row["graph_sha256"] == (
        "c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d"
    )
    assert row["row_count"] == 32
    assert row["hard_row_count"] == row["opposite_row_count"] == 16
    assert row["cut_bank_per_row"] == 496
    assert row["minimum_edit_lower_bound"] == 122
    assert row["maximum_edit_lower_bound"] == 194
    assert row["rows_with_two_label_cut_violations"] == 32
    assert row["full_atom_decomposition_certified"] is False
    assert row["residual_ii_closed"] is False


def test_theorem_record_keeps_the_integral_gate_open() -> None:
    row = theorem_record()
    assert row["proved"] is True
    assert row["incremental_cost_exactly_updated_in_O_p"] is True
    assert row["coefficient_edit_lower_bound_rigorous"] is True
    assert row["full_atom_transport_still_required_at_zero_cost"] is True
    assert row["residual_ii_closed"] is False
