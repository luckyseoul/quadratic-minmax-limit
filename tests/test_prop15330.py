"""Tests for Prop 15.330 — 409 not GJ; shared |H+| den; sq-conf; K."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import (
    mu_part_majorant,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15276 import (
    lemma_D_2plane_amplitudes_proved,
    lemma_D_existence_written,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15330 import (
    HPLUS,
    LIVE_MU,
    LIVE_QPP,
    four_level_exists,
    main,
    named_gj_ints,
    prove_A,
    prove_B,
    prove_C_from_hunter,
    prove_D,
    prove_E,
    prove_open,
)


def test_409_not_named_gj():
    A = prove_A()
    assert A["proved"] is True
    assert A["S7_has_409"] is False
    assert A["q2_is_not_409"] is True
    assert 409 not in named_gj_ints(7)
    assert (7 ** 4) != 409
    assert (7 ** 2) != 409


def test_shared_Hplus_den_not_dimVp():
    B = prove_B()
    assert B["proved"] is True
    assert LIVE_QPP[7].denominator == HPLUS[7] // 14
    assert LIVE_MU[7].denominator == HPLUS[7] // 2
    assert LIVE_MU[7].denominator == 7 * LIVE_QPP[7].denominator
    assert LIVE_QPP[7].denominator != (7 * 7 + 1) // 2
    assert B["dimVp_is_not_den"] is True


def test_sqconf_does_not_recover_unique():
    C = prove_C_from_hunter()
    assert C["proved"] is True
    assert C["unique_recovery"] is False


def test_four_level_exists_does_not_empty_K():
    D = prove_D()
    assert D["proved"] is True
    assert four_level_exists(5) is True
    assert four_level_exists(7) is True
    assert residual_ii_k_eq_4p_empty() is False


def test_maj_does_not_close_type_I():
    E = prove_E()
    assert E["proved"] is True
    assert LIVE_MU[7] > mu_part_majorant(7)
    assert type_I_multilevel_bad_case_ND_closed() is False


def test_floor_and_imports_still_open():
    F = prove_open()
    assert F["proved"] is False
    assert F["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert F["lemma_D_existence"] is True
    assert F["lemma_D_2plane"] is True


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert lemma_D_existence_written() is True
    assert lemma_D_2plane_amplitudes_proved() is True


def test_main():
    out = main()
    assert out["proved"]["409_not_gauss_jacobi"] is True
    assert out["proved"]["shared_Hplus_den"] is True
    assert out["proved"]["sqconf_does_not_raise_rank"] is True
    assert out["proved"]["K_first_moment_does_not_empty"] is True
    assert out["proved"]["maj_does_not_close_type_I"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert Fraction(1544, 409).denominator == 409
