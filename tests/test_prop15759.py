import pytest

from e1_gmin_m4_prop15759 import (
    exact_full_edge_radon_rank,
    moment_polynomial,
    moment_relation_rows,
    p_torsion_codimension,
    signed_residual_transport,
    theorem_record,
    top_degree_basis_complement,
    verify_moment_relations,
)


def test_moment_hierarchy_requires_an_odd_prime():
    for p in (2, 4, 9):
        with pytest.raises(ValueError, match="odd prime"):
            p_torsion_codimension(p)


def test_closed_p_primary_codimension_formula():
    expected = {3: 0, 5: 6, 7: 21, 11: 94, 13: 160}
    for p, codimension in expected.items():
        row = p_torsion_codimension(p)
        assert row["proved"]
        assert row["extra_p_primary_codimension_direct_sum"] == codimension
        assert row["extra_p_primary_codimension_closed_form"] == codimension
        assert row["equal_total_equation_count"] == p
        assert row["parallel_sum_equation_count"] == 1
        assert row["ordinary_boundary_total_codimension"] == p + 1
        assert row["rank_in_characteristic_p"] == row[
            "independent_bidegree_block_rank"
        ]
        assert len(moment_relation_rows(p)) == codimension


def test_top_degree_q_basis_omits_a_valid_ordinary_complement():
    for p in (3, 5, 7, 11, 13):
        row = top_degree_basis_complement(p)
        m = (p - 1) // 2
        assert row["proved"]
        assert len(row["expansion_coefficients_mod_p"]) == m
        assert all(row["expansion_coefficients_mod_p"])
        assert row["retained_new_relation_indices"] == list(range(m - 1))
        assert row["omitted_basis_index"] == m - 1
        assert row["omitted_coefficient_nonzero"]


def test_q_basis_is_symmetric_homogeneous_and_diagonal_vanishing():
    p = 13
    for d in range(2, p):
        for k in range(d // 2):
            for s, t, scale in ((2, 7, 3), (5, 5, 4), (0, 9, 6)):
                value = moment_polynomial(p, s, t, d, k)
                assert value == moment_polynomial(p, t, s, d, k)
                assert moment_polynomial(p, scale * s, scale * t, d, k) == (
                    pow(scale, d, p) * value
                ) % p
            assert moment_polynomial(p, 8, 8, d, k) == 0


def test_all_moment_relations_hold_on_common_graph_samples():
    for p in (3, 5, 7, 11):
        row = verify_moment_relations(p)
        assert row["proved"]
        assert row["relations_checked"] == p_torsion_codimension(p)[
            "extra_p_primary_codimension_closed_form"
        ]


def test_exact_full_edge_radon_ranks():
    for p in (3, 5, 7):
        for modulus in (2, p):
            row = exact_full_edge_radon_rank(p, modulus)
            assert row["proved"]
            assert row["exact_rank"] == row["formula_rank"]


def test_signed_normalization_is_unimodular_transport():
    row = signed_residual_transport()
    assert row["proved"]
    assert row["row_and_column_sign_changes_are_unimodular"]
    assert row["p_primary_codimension_unchanged"]
    assert row["projective_infinity_term_is_signed"]
    assert row["projective_infinity_coefficient"] == "epsilon_infinity"


def test_theorem_record_keeps_compact_survivor_open():
    out = theorem_record()
    assert out["prop"] == "15.759"
    assert out["proved"]["all_endpoint_moment_relations"]
    assert out["proved"]["relations_independent_and_exhaust_p_primary_linear_cokernel"]
    assert not out["proved"]["compact_aggregate_survivor_excluded"]
    assert not out["proved"]["residual_ii_closed"]
    assert out["L_status"] == "OPEN"
