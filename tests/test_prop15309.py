"""Tests for Prop 15.309 — QR0 stab, p=7 dichotomy, NC not general."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15308 import AP_stab
from e1_gmin_m4_prop15309 import (
    QR0_orbit,
    QR0_stab,
    main,
    n_AGL_images_QR0,
    n_AP_subsets,
    prove_NC_not_general,
    prove_open,
    prove_p7_dichotomy,
    prove_QR0_stab,
    qr0,
)


def test_QR0_stab_is_p_pm1_not_4p_at_p7():
    A = prove_QR0_stab()
    assert A["proved"] is True
    assert A["by_p"]["5"]["orbit"] == 30
    assert A["by_p"]["7"]["stab"] == QR0_stab(7) == 42
    assert A["by_p"]["7"]["stab"] != AP_stab(7)
    assert A["by_p"]["7"]["orbit"] == 56 == QR0_orbit(7)
    assert A["by_p"]["11"]["orbit"] == 132
    assert A["by_p"]["5"]["S_is_AP"] is True
    assert A["by_p"]["7"]["S_is_AP"] is False


def test_p7_dichotomy_fails_p11():
    B = prove_p7_dichotomy()
    assert B["proved"] is True
    assert B["p7"]["n_other"] == 0
    assert B["p7"]["n_AP"] + B["p7"]["n_QR_images"] == B["p7"]["Cpm"] == 35
    assert B["p7"]["n_QR_images"] == n_AGL_images_QR0(7) == 14
    p11 = B["p11_not_dichotomy"]
    assert p11["n_AP"] + p11["n_QR_images"] < p11["Cpm"]
    assert p11["n_AP"] == n_AP_subsets(11) == 55
    assert p11["Cpm"] == comb(11, 6) == 462


def test_NC_empty_p13_and_hits_not_p_at_p5():
    C = prove_NC_not_general()
    assert C["proved"] is True
    assert C["by_p"][5]["hits"] == 10 == 2 * 5
    assert C["by_p"][5]["hits"] != 5
    assert C["by_p"][13]["hits"] == 0
    assert C["ystar_p13_found"] is False


def test_qr0_size():
    for p in (5, 7, 11, 13):
        assert len(qr0(p)) == (p + 1) // 2


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
    assert out["proved"]["QR0_stab_p_pm1"] is True
    assert out["proved"]["p7_AP_QR_dichotomy"] is True
    assert out["proved"]["NC_not_general"] is True
    assert out["proved"]["Hplus_named"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
