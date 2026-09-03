from e1_gmin_m4_mobius_fixed_word_atom_coupling import (
    all_prime_collision_room_theorem,
    atom_fixed_incidence_realizable,
    branch_c_atom_quotas,
    coupled_atom_syndrome,
    fixed_word_atom_coupling_theorem,
    mobius_half_block_normal_form,
    p31_three_phi_block_counterexample,
    phi_collision_decomposition,
    phi_collision_ledger,
    theorem_record,
    two_half_phi_block_intersection_theorem,
)
from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_symmetric_fixed_edge_elimination import orbit_fixed_word


def _antipodal_class(p: int, point: tuple[int, int]) -> tuple[int, int]:
    point = point[0] % p, point[1] % p
    negative = (-point[0] % p, -point[1] % p)
    return min(point, negative)


def test_actual_branch_c_atom_quotas_at_both_endpoints() -> None:
    lower = branch_c_atom_quotas(31, 68)
    upper = branch_c_atom_quotas(31, 177)
    assert lower["hard_compact_counts"] == (5,) * 5 + (4,) * 11
    assert lower["opposite_all_equal_counts"] == (6,) * 16
    assert lower["opposite_compact_counts"] == (0,) * 16
    assert lower["total_atom_count"] == 165
    assert upper["hard_compact_counts"] == (12,) * 2 + (11,) * 14
    assert upper["opposite_all_equal_counts"] == (6,) * 16
    assert upper["opposite_compact_counts"] == (7,) * 13 + (6,) * 3
    assert upper["total_atom_count"] == 383


def test_coupled_syndrome_is_direct_binary_addition() -> None:
    # This is an abstract block-coordinate identity, deliberately not a
    # claimed Mobius U.  Desired atom rows have weights 2 and 1.
    result = coupled_atom_syndrome(
        block_parity_rows=((0, 0, 0), (1, 0, 0)),
        atom_counts=(2, 1),
        hard_literal_cells={0: 0},
        singleton_cells=(1, None),
    )
    assert result["required_antipodal_atom_incidence"] == (
        (1, 1, 0),
        (1, 0, 0),
    )
    assert result["required_weights"] == (2, 1)
    assert all(result["atom_count_parities_match"])
    assert all(result["atom_capacities_match"])
    assert result["fixed_word_layer_feasible"]


def test_centered_atoms_give_exact_weight_and_parity_criterion() -> None:
    assert atom_fixed_incidence_realizable(0, 0)
    assert atom_fixed_incidence_realizable(7, 1)
    assert atom_fixed_incidence_realizable(7, 5)
    assert atom_fixed_incidence_realizable(7, 7)
    assert not atom_fixed_incidence_realizable(7, 0)
    assert not atom_fixed_incidence_realizable(7, 6)
    assert not atom_fixed_incidence_realizable(7, 9)


def test_block_functional_normal_form_matches_actual_phi_blocks() -> None:
    p = 31
    center = 7
    edges = mobius_parameter_edges(p, (1, 0), (0, 1), center)
    functionals = []
    for z in range(2, p):
        coefficients = mobius_half_block_normal_form(p, z, center)
        functionals.append(coefficients)
        scaled_x = center * coefficients[0] % p
        scaled_y = center * coefficients[1] % p
        assert ((scaled_x + 2) * (scaled_x + scaled_y) + 1) % p == 0

        expected = {
            _antipodal_class(p, (x, y))
            for x in range(p)
            for y in range(p)
            if (x, y) != (0, 0)
            and (coefficients[0] * x + coefficients[1] * y) ** 2 % p == 1
        }
        record = orbit_fixed_word(p, edges[z - 1])
        actual = {tuple(point) for point in record["fixed_word_support"]}
        assert actual == expected

    assert len(set(functionals)) == p - 2
    assert not any(((-x) % p, (-y) % p) in functionals for x, y in functionals)


def test_exact_collision_identity_and_worst_case_room() -> None:
    top = phi_collision_ledger(31, 177)
    assert top["kappa_total"] == 1
    assert top["forced_Lambda_lower_bound"] == 17
    assert top["forced_sigma_lower_bound"] == 16
    assert top["used_distinct_nonzero_Phi_orbits"] == 462
    assert top["available_raw_upper_margin"] == 215
    assert not top["collision_bounds_contradict"]

    lower_worst = phi_collision_ledger(31, 68, literal_singleton_matches_q=16)
    assert lower_worst["kappa_total"] == 110
    assert lower_worst["forced_Lambda_lower_bound"] == 142
    assert lower_worst["forced_sigma_lower_bound"] == 32
    assert lower_worst["available_raw_upper_margin"] == 90
    room = all_prime_collision_room_theorem(31)
    assert room["minimum_raw_upper_margin"] == 90
    assert not room["scalar_collision_count_excludes_j0"]


def test_collision_decomposition_handles_triple_and_higher_overlaps() -> None:
    result = phi_collision_decomposition(
        {
            "A": (3, 1, 2),
            "B": (1, 1, 1),
            "C": (4, 5),
        }
    )
    assert result["raw_nonzero_Phi_occurrences"] == 18
    assert result["used_distinct_orbits_by_block"] == {"A": 2, "B": 3, "C": 1}
    assert result["used_distinct_nonzero_Phi_orbits"] == 6
    assert result["block_parity"] == {"A": 0, "B": 1, "C": 1}
    assert result["C"] == 2
    assert result["kappa_nonzero"] == 6
    assert result["sigma"] == 2
    assert result["Lambda"] == 8
    assert result["Lambda"] == result["kappa_nonzero"] + result["sigma"]


def test_pair_bound_one_is_false_but_conic_bound_eight_holds() -> None:
    theorem = two_half_phi_block_intersection_theorem(31)
    assert theorem["one_half_blocks_are_distinct"]
    assert theorem["per_orientation_intersection_bound"] == 4
    assert theorem["two_orientation_block_intersection_bound"] == 8
    assert not theorem["bound_one_valid"]

    witness = p31_three_phi_block_counterexample()
    assert witness["both_target_directions_Paley_hard"]
    assert witness["shared_Phi_block_parameter_pairs"] == (
        (9, 20),
        (19, 12),
        (25, 18),
    )
    assert witness["shared_Phi_block_count"] == 3
    assert witness["common_physical_edges"] == 0
    assert witness["sum_is_ternary"]


def test_theorem_keeps_the_actual_endpoint_open() -> None:
    coupling = fixed_word_atom_coupling_theorem(31, 177)
    assert coupling["block_basis_identity"] == "M^T*a_Y=ell+z"
    assert coupling["coupled_singleton_identity"] == "z=c_U+ell+s_x"
    assert not coupling["even_common_moments_and_nonfixed_target_cells_solved"]

    result = theorem_record()
    assert result["proved_all_claimed_statements"]
    assert not result["uniform_j0_exclusion_proved"]
    assert not result["j0_construction_proved"]
    assert not result["residual_ii_closed"]
