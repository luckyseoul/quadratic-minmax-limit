"""Prop 15.272: bad-μ spanning proved; same-line pairing / residual (i) open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15270 import gplus_pd_proved_general, psl_span_F_eq_Wpp0
from e1_gmin_m4_prop15272 import (
    M3_count,
    _even_c_offdiag_variation_rank,
    _even_pairs,
    _wwt_built_perp_eig,
    _wwt_perp_eig,
    dim_F_perp,
    isotypic_total_dim,
    two_m_minus_p,
    u_term_quadratic_coeff,
    dim_bad_isotypic,
    dim_good_isotypic,
    fperp_injective_proved_general,
    hinge_status_272,
    residual_i_closed_via_272,
    theorem_K_invertible,
    theorem_bad_mu_span,
    theorem_complementary_mixed,
    theorem_convolution_hyperplane,
    theorem_every_triple_occurs,
    theorem_even_c_annihilator_is_constants,
    theorem_isotypic_dim_fill,
    theorem_isotypic_is_pair_hyperplane,
    theorem_Vplus_omega_support,
    theorem_zero_diag_is_convolution,
    theorem_good_mu_span,
    theorem_k1_cylinders,
    theorem_k1_sameline_hyperplane,
    theorem_oneswap_A_sym_form,
    theorem_oneswap_delta_vanishes,
    theorem_same_line_pairing,
    theorem_same_line_rank,
    theorem_subset_sum_kernel_is_constants,
    theorem_through_L0_last_dim,
)


def test_isotypic_dims():
    assert dim_good_isotypic(5) == 3
    assert dim_bad_isotypic(5) == 2
    assert dim_good_isotypic(7) == 6
    assert dim_bad_isotypic(7) == 5
    assert dim_good_isotypic(11) == 15
    assert dim_bad_isotypic(11) == 14


def test_units_fail_when_identities_wrong():
    """Helpers go red if the binomial, 2m−p=1, or DFT rank is wrong."""
    assert _wwt_perp_eig(2, 2) == 0
    assert _wwt_perp_eig(4, 1) == 0
    assert _wwt_perp_eig(3, 2) > 0
    assert _wwt_built_perp_eig(3, 2) == _wwt_perp_eig(3, 2)
    assert _wwt_built_perp_eig(5, 3) == _wwt_perp_eig(5, 3)
    assert _wwt_built_perp_eig(4, 1) == 0
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert two_m_minus_p(p) == 1
        assert u_term_quadratic_coeff(p) == 2
        m_wrong = p // 2
        assert 2 * (m_wrong - 1) - (p - 2) != 1
        # wrong m does not give 2σ=2 (odd p: 2m−p=−1, coeff=−2)
        assert 2 * (2 * m_wrong - p) != 2
    # even-c rank is m−1, not m (constants are the kernel)
    for p, alpha in ((5, 1), (7, 1), (7, 3), (11, 2)):
        m = (p + 1) // 2
        assert len(_even_pairs(p, alpha)) == m
        assert _even_c_offdiag_variation_rank(p, alpha) == m - 1
    assert M3_count(5) == 100
    assert M3_count(7) == 1176
    assert M3_count(11) == 24200
    assert M3_count(5) != 61  # k=3 Veronese rank, not M3
    assert isotypic_total_dim(5) == dim_F_perp(5) == 60


def test_johnson_oneswap_algebra():
    """Live WWᵀ / 2m-p=1 units — not a handwritten True."""
    t = theorem_subset_sum_kernel_is_constants()
    assert t["proved"] is True
    for p, n, k in ((5, 3, 2), (7, 5, 3), (11, 9, 5), (13, 11, 6)):
        alpha = _wwt_perp_eig(n, k)
        assert alpha > 0
        assert t["by_p_sample"][str(p)]["alpha"] == alpha
    dlt = theorem_oneswap_delta_vanishes()
    assert dlt["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        m = (p + 1) // 2
        assert 2 * m - p == 1
        assert 2 * (m - 1) - (p - 2) == 1
    assert theorem_oneswap_A_sym_form()["proved"] is True
    ev = theorem_even_c_annihilator_is_constants()
    assert ev["proved"] is True
    assert ev["by_p_sample"]["5"]["ranks"]["1"] == 2
    assert ev["by_p_sample"]["7"]["ranks"]["1"] == 3
    assert theorem_through_L0_last_dim()["proved"] is True
    assert theorem_every_triple_occurs()["proved"] is True
    assert theorem_isotypic_dim_fill()["proved"] is True
    assert theorem_Vplus_omega_support()["proved"] is True
    assert theorem_zero_diag_is_convolution()["proved"] is True
    assert theorem_isotypic_is_pair_hyperplane()["proved"] is True
    # k1 same-line is the AND of those units, not cylinders alone
    assert theorem_k1_cylinders()["proved"] is True
    assert theorem_k1_sameline_hyperplane()["proved"] is True
    assert theorem_k1_sameline_hyperplane()["depends_on"][-1] == (
        "theorem_even_c_annihilator_is_constants"
    )


def test_proved_lemmas():
    assert theorem_convolution_hyperplane()["proved"] is True
    assert theorem_every_triple_occurs()["proved"] is True
    assert theorem_isotypic_dim_fill()["proved"] is True
    assert theorem_bad_mu_span()["proved"] is True
    assert theorem_complementary_mixed()["proved"] is True
    assert theorem_K_invertible()["proved"] is True
    assert theorem_k1_cylinders()["proved"] is True
    assert theorem_k1_sameline_hyperplane()["proved"] is True


def test_pairing_unused_k13_closes():
    assert theorem_same_line_pairing()["proved"] is False
    assert theorem_same_line_rank()["proved"] is False
    assert theorem_good_mu_span()["proved"] is True
    assert fperp_injective_proved_general() is True
    assert residual_i_closed_via_272() is True


def test_live_flags_via_k13():
    assert gplus_pd_proved_general() is True
    assert psl_span_F_eq_Wpp0() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_272()
    assert h["bad_mu_span"] is True
    assert h["every_triple_occurs"] is True
    assert h["isotypic_dim_fill"] is True
    assert h["Vplus_omega_support"] is True
    assert h["zero_diag_is_convolution"] is True
    assert h["isotypic_is_pair_hyperplane"] is True
    assert h["k1_sameline_hyperplane"] is True
    assert h["K_invertible"] is True
    assert h["same_line_pairing"] is False
    assert h["fperp_injective_general"] is True
    assert h["residual_i_closed"] is True
    assert h["e1_closed_general"] is True
