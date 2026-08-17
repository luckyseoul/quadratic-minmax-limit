"""Tests for Prop 15.308 — |G|, AP stab, orbit-stab sum."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15308 import (
    AP_orbit,
    AP_stab,
    G_order,
    G_order_no_frob,
    main,
    prove_AP_stab,
    prove_G_order,
    prove_open,
    prove_orbit_sum,
)


def test_G_order_needs_frob():
    A = prove_G_order()
    assert A["proved"] is True
    for p in (3, 5, 7):
        assert A["by_p"][str(p)]["G"] == G_order(p)
        assert A["by_p"][str(p)]["G_no_frob"] == G_order_no_frob(p)
        assert G_order_no_frob(p) * 2 == G_order(p)


def test_AP_stab_is_4p_not_2p():
    B = prove_AP_stab()
    assert B["proved"] is True
    for p in (3, 5, 7):
        assert B["by_p"][str(p)]["stab"] == AP_stab(p) == 4 * p
        assert B["by_p"][str(p)]["orbit"] == AP_orbit(p)
        assert B["by_p"][str(p)]["stab"] != 2 * p


def test_orbit_sum_needs_AP_split():
    C = prove_orbit_sum()
    assert C["proved"] is True
    assert C["by_p"]["5"]["sum_orbits"] == 130
    assert C["by_p"]["7"]["sum_orbits"] == 5726
    assert C["by_p"]["7"]["n_types"] >= 7
    # p=7 non-AP hs is the extra 56
    types7 = C["by_p"]["7"]["types"]
    na = next(t for t in types7 if t["type"] == "('hs-nonAP',)")
    assert na["orbit"] == 56
    assert na["stab"] == 42


def test_p7_NL_not_unique():
    C = prove_orbit_sum()
    nl = [t for t in C["by_p"]["7"]["types"] if not t["type"].startswith("('hs-")]
    assert len(nl) >= 6
    assert any(t["stab"] == 8 for t in nl)  # p+1
    assert any(t["stab"] == 2 for t in nl)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Hplus_named"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["lambda_ge_6"] is False
    assert D["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["G_order"] is True
    assert out["proved"]["AP_stab_4p"] is True
    assert out["proved"]["orbit_stab_sum"] is True
    assert out["proved"]["Hplus_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
