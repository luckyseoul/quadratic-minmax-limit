"""Tests for Prop 15.409 — D lattice; D_form off for p≥11."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15357 import D_live
from e1_gmin_m4_prop15371 import D_form
from e1_gmin_m4_prop15409 import (
    D_1d,
    D_1d_wrong_m,
    D_form_on_lattice_general,
    D_lattice,
    dform_minus_d1d_mod_p,
    live_u,
    main,
    n1d_two_forms,
    on_lattice,
    paley_region_counts,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_E,
    prove_open,
    u_window,
)


def test_congruence_live_and_form11():
    A = prove_A()
    assert A["proved"] is True
    assert (130 - n_1d(5)) % 25 == 0
    assert (5726 - n_1d(7)) % 49 == 0
    Nf11 = 2 * 11 * D_form(11)
    assert (Nf11 - n_1d(11)) % 11 == 0
    assert (Nf11 - n_1d(11)) % 121 != 0


def test_lattice_named():
    B = prove_B()
    assert B["proved"] is True
    assert D_1d(5) == 3
    assert D_1d(7) == 10
    assert D_1d(11) == 126
    assert D_1d_wrong_m(5) == 9
    assert D_1d_wrong_m(5) != D_1d(5)
    assert n1d_two_forms(5) and n1d_two_forms(7) and n1d_two_forms(11)
    assert D_live(5) == D_lattice(5, 2)
    assert D_live(7) == D_lattice(7, 57)
    assert live_u(5) == 2
    assert live_u(7) == 57
    assert 2 * 5 * D_1d(5) == n_1d(5)


def test_D_form_off_lattice():
    C = prove_C()
    assert C["proved"] is True
    assert on_lattice(5, D_form(5))
    assert on_lattice(7, D_form(7))
    assert not on_lattice(3, D_form(3))
    assert not on_lattice(11, D_form(11))
    assert (D_form(11) - D_1d(11)) % 11 == 35 % 11
    assert dform_minus_d1d_mod_p(5) == 30
    assert dform_minus_d1d_mod_p(7) == 35
    assert dform_minus_d1d_mod_p(13) == 30
    assert 30 % 13 != 0
    assert 35 % 11 != 0
    assert D_form_on_lattice_general() is False


def test_paley_region_misses_409():
    D = prove_D()
    assert D["proved"] is True
    c7 = paley_region_counts(7)
    assert c7["parallelogram"] == 144
    assert c7["pp_raw"] == 6
    assert 409 not in c7.values()
    assert 13 not in paley_region_counts(5).values()


def test_u_window_pins_p5():
    E = prove_E()
    assert E["proved"] is True
    assert u_window(5) == (2, 2)
    assert u_window(7) == (53, 64)
    assert 57 in range(53, 65)
    assert 0 not in range(2, 3)


def test_flags_untouched():
    F = prove_open()
    assert F["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_aut_e_3AB_positive_general() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    assert D_form_on_lattice_general() is False


def test_main():
    out = main()
    assert out["proved"]["Hplus_cong_n1d_mod_psq"] is True
    assert out["proved"]["D_lattice_D1d_plus_p_u"] is True
    assert out["proved"]["D_form_off_lattice_p_ge_11"] is True
    assert out["proved"]["paley_region_misses_409"] is True
    assert out["proved"]["u_window_pins_p5"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
    assert out["prop"] == "15.409"
