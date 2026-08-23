"""Tests for Prop 15.611 — W cyclic F2[M]; ker2 dim 2 is a p-law."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15611 import (
    theorem_A_cyclic_regular,
    theorem_B_unique_hyperplane,
    theorem_C_ker2_dim_p_law,
    theorem_D_open,
    _v2,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_v2_N_at_least_2_odd_p():
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        N = (p * p - 1) // 2
        assert N % 4 == 0
        assert _v2(N) >= 2


def test_W_cyclic_regular():
    A = theorem_A_cyclic_regular()
    assert A["proved"] is True
    for p, rec in A["rows"].items():
        assert rec["dim_W"] == rec["N"]
        assert rec["D_N_is_I"] is True
        assert rec["W_cyclic_D"] is True
        assert rec["dim_WH"] == int(p) - 1
        assert rec["WH_cyclic_Fpstar"] is True
        assert rec["n_nsq"] == rec["index_M_over_Fpstar"]
        assert rec["WH_orbit_Cp"] == int(p) - 1
    assert A["rows"]["3"]["WH_simple_Cp"] is True
    assert A["rows"]["5"]["WH_simple_Cp"] is True
    assert A["rows"]["7"]["WH_simple_Cp"] is False


def test_unique_hyperplane_W0():
    B = theorem_B_unique_hyperplane()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["rank_im_DminusI"] == rec["N_minus_1"]
        assert rec["im_vanishes_at_0"] is True
        assert rec["dim_W0"] == rec["N_minus_1"]
        assert rec["extra_at_0"] == 0


def test_ker2_dim_p_law():
    C = theorem_C_ker2_dim_p_law()
    assert C["proved"] is True
    assert C["ker2_dim_is_p_law"] is True
    for rec in C["rows"].values():
        assert rec["v2_N"] >= 2
        assert rec["four_divides_N"] is True
        assert rec["dim_ker1_W0"] == 1
        assert rec["dim_ker2_W0"] == 2


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15611 import main

    out = main()
    assert out["proved"]["W_cyclic_F2M"] is True
    assert out["proved"]["W0_unique_D_hyperplane"] is True
    assert out["proved"]["ker2_dim_is_p_law"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
