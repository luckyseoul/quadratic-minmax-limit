from src.e1_gmin_m4_prop15644 import (
    baseline_parallel_bound,
    exceptional_parallel_bounds,
    negative_branch_normal_form,
    theorem_negative_branch_normal_form,
)


def test_baseline_and_exceptional_direction_bounds():
    for p in (31, 41, 101, 201):
        infinity_edges = 2 * p - 1
        assert baseline_parallel_bound(p, infinity_edges) == 3
        assert exceptional_parallel_bounds(p, infinity_edges) == (1, 4)


def test_unique_large_p_normal_form():
    for p in (31, 33, 41, 101, 201):
        result = negative_branch_normal_form(p)
        assert result["unique_normal_form"] is True
        normal = result["normal_form"]
        assert normal["residue"] == 2
        assert normal["k0"] == 4
        assert normal["infinity_edges"] == 2 * p - 1
        assert normal["finite_edges"] == 2 * p + 2
        assert normal["exceptional_pairs_positive_negative"] == [(1, 3), (3, 1)]


def test_theorem_is_a_normal_form_not_a_close():
    theorem = theorem_negative_branch_normal_form()
    assert theorem["proved"] is True
    assert theorem["all_sufficiently_large_odd_primes"] is True
    assert theorem["closes_negative_product_branch"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
