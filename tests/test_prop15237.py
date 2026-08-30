"""Prop 15.237 — residual (ii-a) ND closed; residual (i)/E1/L still open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed
from e1_gmin_m4_prop15193 import (
    residual_ii_bounded_even_k_le_4p_minus_2_closed,
    residual_ii_full_closed,
    residual_ii_open_subcases,
)
from e1_gmin_m4_prop15236 import residual_ii_a_ND_closed as iia_236
from e1_gmin_m4_prop15236 import residual_ii_b_ND_closed
from e1_gmin_m4_prop15237 import (
    even_type_pair_matrix,
    lemma_dual_bad_mass,
    lemma_even_types_span_R6,
    lemma_k4_rank_p5,
    lemma_pair_slice_mass_mismatch,
    lemma_three_equal_SH_values,
    lemma_triangle_gram_pd,
    lemma_two_edge_three_values,
    lemma_wedge_G_beats_thr,
    lemma_weight_hits_seven_types,
    main,
    residual_ii_a_ND_closed,
    theorem_A_mass,
    theorem_B_quadrants,
    theorem_D_U_empty,
    theorem_independence,
)


def test_mass_and_wedge():
    for p in (5, 7, 11, 13, 17, 19):
        assert lemma_wedge_G_beats_thr(p)["proved"] is True
        assert lemma_triangle_gram_pd(p)["proved"] is True
        assert lemma_weight_hits_seven_types(p)["proved"] is True
        from e1_gmin_m4_prop15236 import even_k_residual_ii

        for k in even_k_residual_ii(p):
            assert lemma_dual_bad_mass(p, k)["proved"] is True
            mm = lemma_pair_slice_mass_mismatch(p, k)
            assert mm["proved"] is True


def test_classification_algebra():
    assert lemma_even_types_span_R6()["proved"] is True
    M = even_type_pair_matrix()
    assert M.shape[1] == 6
    assert lemma_two_edge_three_values()["proved"] is True
    from e1_gmin_m4_prop15237 import (
        lemma_disjoint_product_not_pair_span,
        lemma_support_star_or_triangle,
    )

    assert lemma_disjoint_product_not_pair_span()["proved"] is True
    assert lemma_support_star_or_triangle()["proved"] is True
    sh = lemma_three_equal_SH_values()
    assert sh["proved"] is True
    assert -2.0 in sh["SH_on_3equal"]


def test_theorems_and_predicate():
    assert theorem_A_mass()["proved"] is True
    assert theorem_B_quadrants()["proved"] is True
    assert theorem_D_U_empty()["proved"] is True
    assert theorem_independence()["proved"] is True
    assert residual_ii_a_ND_closed() is True
    assert iia_236() is True  # lazy re-export
    assert residual_ii_b_ND_closed() is True
    assert residual_ii_bounded_even_k_le_4p_minus_2_closed() is True
    assert residual_ii_full_closed() is False
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False


def test_p5_k4_independence_anchor():
    """Recompute the p=5 rank-6 census used as an independence anchor."""
    r = lemma_k4_rank_p5()
    assert r["proved"] is True
    assert r["n_rank_lt_6"] == 0


def test_does_not_soft_close_e1():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    opens = residual_ii_open_subcases()
    assert not any("ii-a" in s for s in opens)
    assert not any("ii-b" in s for s in opens)
    assert any("k≥4p" in s for s in opens)
    out = main()
    assert out["proved"]["residual_ii_a_ND_closed"] is True
    assert out["proved"]["E1_closed"] is False
    assert gsum_disj_lb_proved_general() is False
    assert residual_ii_a_ND_closed.__module__.endswith("prop15237")
