"""Tests for Prop 15.332 — NL ∩ Trans={id}; U⊆□."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15332 import (
    D_size,
    U_subseteq_squares,
    line_union_count,
    main,
    n_directions,
    nl_stab_N_ones,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_D_size_is_r_times_p():
    A = prove_A()
    assert A["proved"] is True
    assert A["D_always_0_mod_p"] is True
    for p in (5, 7, 11):
        assert D_size(p) % p == 0
        assert line_union_count(p) * p == D_size(p)


def test_nhs_equals_n1d_not_all_directions():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["n_hs"] == n_1d(5)
    assert B["by_p"]["7"]["n_hs"] == n_1d(7)
    assert B["by_p"]["7"]["n_hs_is_not_all_dir"] is True
    assert B["by_p"]["5"]["n_all_dir_choose"] != n_1d(5)


def test_U_subseteq_squares():
    C = prove_C()
    assert C["proved"] is True
    assert U_subseteq_squares(5) is True
    assert U_subseteq_squares(7) is True
    assert U_subseteq_squares(11) is True


def test_NL_stab_N_eq_1_embeds_U_rtimes_Frob():
    D = prove_D()
    assert D["proved"] is True
    assert D["embeds_in_U_rtimes_Frob"] is True
    assert D["order_divides_pplus1_general"] is False
    r5 = nl_stab_N_ones(5)
    assert r5["all_N_eq_1"] is True
    assert r5["n_frob"] > 0
    assert r5["n_N_not_1"] == 0


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["stab_embeds_in_U_rtimes_Frob"] is True
    assert E["stab_divides_pplus1_general"] is False
    assert E["Hplus_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["line_union"] is True
    assert out["proved"]["NL_trans_trivial"] is True
    assert out["proved"]["U_subseteq_squares"] is True
    assert out["proved"]["stab_embeds_in_U_rtimes_Frob"] is True
    assert out["proved"]["stab_divides_pplus1_general"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert n_directions(7) == 8
