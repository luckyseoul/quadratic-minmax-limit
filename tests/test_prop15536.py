"""Tests for Prop 15.536 — I_Ω = p S(r/(1+r)); doubles are CM."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15275 import (
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15523 import G_chi
from e1_gmin_m4_prop15536 import (
    S_minus1_named,
    ap2_cm,
    chi_Omega,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_I_Omega_equals_p_S():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        rec = A["rows"][str(p)]
        assert rec["I_eq_pS"] is True
        assert rec["n_ok"] == rec["n_tri"]
        assert chi_Omega(p) == (-1 if p % 4 == 1 else 1)
    # drop χ_Ω collides only when χ_Ω=1 (p≡3); fail-when-wrong is p=5
    assert A["rows"]["5"]["n_GS_collisions"] == 0
    # Orin / 15.529 doubles
    assert A["p5_I_dbl"] == [30]
    assert A["rows"]["7"]["I_dbl"] == [98]
    assert 98 == 2 * 7 * 7
    assert 30 != 2 * 5 * 5


def test_drop_chi_Omega_and_two_p2_fail():
    A = prove_A()
    assert A["p5_GS"] == -30
    assert A["p5_GS"] != 30
    assert A["p5_GS"] == G_chi(5) * A["rows"]["5"]["S_half"]
    assert 30 != 50


def test_doubles_are_CM_S_minus1():
    B = prove_B()
    assert B["proved"] is True
    assert S_minus1_named(5) == 6
    assert S_minus1_named(7) == 14
    assert ap2_cm(5) == 4
    assert ap2_cm(7) == 0
    assert B["by_p"]["5"]["I_dbl"] == [30]
    assert B["by_p"]["7"]["I_dbl"] == [98]
    assert B["by_p"]["5"]["is_2p2"] is False
    assert B["by_p"]["7"]["is_2p2"] is True
    assert B["by_p"]["11"]["is_2p2"] is True
    assert B["by_p"]["13"]["is_2p2"] is False
    # fail-when-wrong: I_dbl=2p² at all p
    assert B["by_p"]["5"]["S_m1"] != 2 * 5


def test_type_I_flags_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["type_I_multilevel_bad_case_ND_closed"] is True
    assert C["type_I_aut_e_3AB_positive_general"] is False
    assert C["Delta_conn_p_rational"] is False
    assert C["generic_S_is_p_law"] is False
    assert C["split_p11"] is True
    assert set(C["type111_p7"]) == {-2, 14}
    assert C["S_dist_p5"] == [-10]
    assert C["S_dist_p7"] == [-2]
    assert type_I_multilevel_bad_case_ND_closed() is True
    assert type_I_aut_e_3AB_positive_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["proved"]["I_Omega_eq_p_S"] is True
    assert out["proved"]["I_dbl_CM"] is True
    assert out["proved"]["type_I_multilevel_bad_case_ND_closed"] is True
    assert out["proved"]["type_I_aut_e_3AB_positive_general"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
