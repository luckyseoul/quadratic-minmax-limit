"""Prop 15.236 — residual (ii-b) ND closed; E1/L still open (residual i)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed
from e1_gmin_m4_prop15193 import residual_ii_full_closed, residual_ii_open_subcases
from e1_gmin_m4_prop15236 import (
    B_minus_first_moment,
    B_plus_first_moment,
    binom_p_half,
    even_k_residual_ii,
    lemma_cleared_identity_mod_p,
    lemma_k_not_div_p,
    lemma_max_minus_ge_minus_2,
    lemma_N_mod_p2,
    lemma_Q_signs,
    lemma_slope_difference_identity,
    lemma_vp_binom,
    main,
    residual_ii_a_ND_closed,
    residual_ii_b_ND_closed,
    residue_N_mod_p2,
    theorem_A_all,
    theorem_C_all,
    theorem_D_dual_bad_empty,
    two_level_mass,
    wedge_slope_minus,
    wedge_slope_plus,
)


def test_even_k_range_and_mass():
    for p in (5, 7, 11, 13, 17, 19):
        ks = even_k_residual_ii(p)
        assert ks[0] >= 3 * p + 1
        assert ks[-1] <= 4 * p - 2
        assert all(k % 2 == 0 for k in ks)
        for k in ks:
            a = two_level_mass(p, k)
            assert 0 < a < 1
            assert lemma_max_minus_ge_minus_2(p, k)["proved"] is True
            assert lemma_k_not_div_p(p, k)["proved"] is True


def test_Q_signs_and_slopes():
    for p in (5, 7, 11, 13):
        assert lemma_Q_signs(p)["proved"] is True
        for k in even_k_residual_ii(p)[:3]:
            sl = lemma_slope_difference_identity(p, k)
            assert sl["proved"] is True
            assert B_plus_first_moment(p, k) == -B_minus_first_moment(p, k)
            wd = wedge_slope_plus(p, 1, 2) - wedge_slope_minus(p, 1, 2)
            from fractions import Fraction

            assert wd == Fraction(p, p * p - 1) * (4 - 2)


def test_translation_count_mod_p2():
    th = theorem_C_all([5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert th["proved"] is True
    assert th["census_congruence"] is True
    for p, N in ((5, 260), (7, 11452), (3, 12)):
        assert N % (p * p) == residue_N_mod_p2(p)
        assert lemma_vp_binom(p)["v_p"] == 1
        assert lemma_N_mod_p2(p)["p2_does_not_divide_N"] is True
        assert binom_p_half(p) % (p * p) != 0


def test_cleared_identity_rhs_not_div_p():
    for p in (5, 7, 11, 13, 17):
        for k in even_k_residual_ii(p):
            cl = lemma_cleared_identity_mod_p(p, k)
            assert cl["proved"] is True
            assert cl["contradiction"] is True
            assert cl["p_divides_rhs_coeff"] is False


def test_theorems_and_predicate():
    assert theorem_A_all()["proved"] is True
    assert theorem_D_dual_bad_empty()["proved"] is True
    assert residual_ii_b_ND_closed() is True
    assert residual_ii_a_ND_closed() is True


def test_wiring_does_not_soft_close_e1():
    """(ii-b) ND is imported into residual_ii_full; e1 still needs residual (i)."""
    assert residual_ii_b_ND_closed() is True
    assert residual_ii_full_closed() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    opens = residual_ii_open_subcases()
    assert not any("ii-a" in s for s in opens)
    assert not any("ii-b" in s for s in opens)
    out = main()
    assert out["proved"]["residual_ii_b_ND_closed"] is True
    assert out["proved"]["residual_ii_a_ND_closed"] is True
    assert out["proved"]["E1_closed"] is True
    assert gsum_disj_lb_proved_general() is False
    # hinge is the real function, not a stub True
    assert residual_ii_b_ND_closed.__module__.endswith("prop15236")
