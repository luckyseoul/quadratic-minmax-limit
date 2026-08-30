"""Prop 15.270: G₊≻0 via k=3 Singer Gram; live residual-(i) flags."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15207 import ker_sc_proved_general
from e1_gmin_m4_prop15249 import residual_i_closed_via_249
from e1_gmin_m4_prop15270 import (
    SMALL_DFT_PRIMES,
    certify_p5_k3_misses_Wpp0,
    certify_small_dft,
    gplus_pd_proved_general,
    psl_span_F_eq_Wpp0,
    singer_gram_dft_min,
    theorem_aut_schur,
    theorem_fejer_stretch,
    theorem_gershgorin_large_p,
    theorem_mu0_block,
    theorem_no_cuspidal,
    theorem_weil_T,
)


def test_algebra_lemmas():
    assert theorem_fejer_stretch()["proved"] is True
    assert theorem_mu0_block()["proved"] is True
    assert theorem_weil_T()["proved"] is True
    assert theorem_gershgorin_large_p()["proved"] is True
    assert theorem_no_cuspidal()["proved"] is True


def test_small_dft_live():
    c = certify_small_dft()
    assert c["certified"] is True
    for p in SMALL_DFT_PRIMES:
        assert c["by_p_min_abs_dft"][str(p)] > 1e-8
    assert singer_gram_dft_min(7) > 0.5


def test_failed_aut_schur_lift():
    """Jacquet ⇏ PSL-span(k=3 F)=W++0. Close is 15.272, not Aut-Schur."""
    w = certify_p5_k3_misses_Wpp0()
    assert w["dim_Wpp0"] == 65
    assert w["k3_veronese_rank"] == 61
    assert w["spans"] is False
    assert psl_span_F_eq_Wpp0() is False
    assert theorem_aut_schur()["proved"] is False
    assert gplus_pd_proved_general() is True
    assert ker_sc_proved_general() is True
    assert residual_i_closed_via_249() is True
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
