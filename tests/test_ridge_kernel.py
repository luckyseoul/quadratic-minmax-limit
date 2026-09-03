from e1_gmin_m4_ridge_kernel import (
    affine_radon_inversion,
    canonical_ridge_basis,
    edge_radon_image,
    half_difference_classes,
    parallel_ridge_move,
    projective_functionals,
    ridge_kernel_theorem_record,
    ridge_mod_p_dependency_certificate,
    ridge_psaturation_certificate,
    transverse_ridge_move,
)


def _value(functional, point, p):
    return (functional[0] * point[0] + functional[1] * point[1]) % p


def _simple_profile(p):
    return (1, -1) + (0,) * (p - 2)


def test_affine_radon_inversion_is_exact_not_numerical():
    for p in (3, 5):
        values = tuple(((7 * index + 3) % 9) - 4 for index in range(p * p - 1))
        values += (-sum(values),)
        row = affine_radon_inversion(p, values)
        assert row["proved"]
        assert row["reconstruction"] == [p * value for value in values]


def test_parallel_and_transverse_ridges_have_the_claimed_support_and_norm():
    for p in (3, 5):
        directions = projective_functionals(p)
        differences = half_difference_classes(p)
        functional = directions[0]
        parallel = [
            index
            for index, delta in enumerate(differences)
            if _value(functional, delta, p) == 0
        ]
        assert len(parallel) == (p - 1) // 2
        type_p = parallel_ridge_move(p, 0, parallel[0], _simple_profile(p))
        assert edge_radon_image(p, type_p) == {}
        assert sum(value != 0 for value in type_p) == 2 * p
        assert sum(value * value for value in type_p) == 2 * p

        by_square = {}
        for index, delta in enumerate(differences):
            projected = _value(functional, delta, p)
            if projected:
                by_square.setdefault(projected * projected % p, []).append(index)
        same_square = next(indices for indices in by_square.values() if len(indices) >= 2)
        type_k = transverse_ridge_move(
            p, 0, same_square[0], same_square[1], _simple_profile(p)
        )
        assert edge_radon_image(p, type_k) == {}
        assert sum(value != 0 for value in type_k) == 4 * p
        assert sum(value * value for value in type_k) == 4 * p


def test_p_saturation_reconstructs_a_nontrivial_mixed_kernel_vector():
    # These are fail-when-wrong identity checks for the all-prime proof, not a
    # finite-prime census or evidence for Boolean existence.
    for p in (3, 5):
        directions = projective_functionals(p)
        differences = half_difference_classes(p)
        profile = _simple_profile(p)
        pieces = []
        for direction_index, functional in enumerate(directions[:2]):
            parallel = next(
                index
                for index, delta in enumerate(differences)
                if _value(functional, delta, p) == 0
            )
            pieces.append(
                parallel_ridge_move(p, direction_index, parallel, profile)
            )
            by_square = {}
            for index, delta in enumerate(differences):
                projected = _value(functional, delta, p)
                if projected:
                    by_square.setdefault(projected * projected % p, []).append(index)
            same_square = next(indices for indices in by_square.values() if len(indices) >= 2)
            pieces.append(
                transverse_ridge_move(
                    p,
                    direction_index,
                    same_square[0],
                    same_square[1],
                    tuple(-value for value in profile),
                )
            )
        source = tuple(sum(piece[index] for piece in pieces) for index in range(len(pieces[0])))
        assert source != (0,) * len(source)
        assert edge_radon_image(p, source) == {}
        row = ridge_psaturation_certificate(p, source)
        assert row["proved"]
        assert row["ridge_decomposition_equals_p_times_source"]
        assert row["p_kernel_contained_in_ridge_lattice"]


def test_ridge_counts_equal_the_full_rational_kernel_dimension():
    for p in (3, 5, 7, 11):
        row = ridge_kernel_theorem_record(p)
        assert row["canonical_ridge_total"] == row["integer_kernel_rank"]
        assert row["proved"]["elementary_parallel_ridges_are_graver"]
        assert (
            row["proved"]["elementary_transverse_ridges_are_graver_for_p_at_least_5"]
            == (p >= 5)
        )
        assert row["proved"]["p_times_integer_kernel_is_in_ridge_lattice"]
        assert row["proved"]["ridge_moves_span_rational_kernel"]
        assert row["proved"]["kernel_mod_ridge_is_finite_elementary_p_torsion"]
        m = (p - 1) // 2
        expected_quotient_dimension = (
            (p + 1) * p * m * m + m * (m - 1) * (4 * m + 1) // 6
        )
        assert row["ridge_quotient_dimension_nu_p"] == expected_quotient_dimension
        assert row["proved"]["ridge_quotient_dimension_formula"]
        assert row["proved"]["ridge_lattice_is_proper"]
        assert row["proved"]["one_p_saturation_recovers_full_integer_kernel"]
        assert row["proved"]["mod_p_dependencies_parameterize_saturating_moves"]
        assert row["proved"]["displayed_ridges_are_not_the_complete_graver_basis"]
        assert (
            row["minimum_nonridge_graver_elements"]
            == 2 * expected_quotient_dimension
        )
        assert row["proved"]["complete_graver_basis_has_at_least_2nu_nonridge_elements"]
        assert not row["proved"]["residual_ii_closed"]


def test_canonical_ridge_matrix_has_the_proved_mod_p_dependency_nullity():
    # These are small fail-when-wrong checks of the uniform index proof, not
    # searches over residual targets or evidence by prime extrapolation.
    for p in (3, 5):
        columns = canonical_ridge_basis(p)
        row = ridge_mod_p_dependency_certificate(p)
        assert len(columns) == row["ridge_columns"]
        assert row["dependency_nullity_nu_p"] == row["closed_formula"]
        assert row["matches_uniform_quotient_dimension"]
        assert row["finite_identity_check_not_target_census"]
