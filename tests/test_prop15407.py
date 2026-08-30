"""Tests for Prop 15.407 — Q_τ is not a GJ/Fermat integer ratio."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15407 import (
    Qtau_is_gj_ratio,
    catalog_ratio_equals_live,
    expanded_catalog,
    fermat_count,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_409_absent_fermat_not_409():
    A = prove_A()
    assert A["proved"] is True
    assert 409 not in expanded_catalog(7)
    assert fermat_count(7, 8) == 336
    assert fermat_count(7, 4) == 88
    assert fermat_count(7, 8) != 409
    assert 48 not in expanded_catalog(5)


def test_no_catalog_ratio_is_live_Q():
    B = prove_B()
    assert B["proved"] is True
    assert catalog_ratio_equals_live(5) is False
    assert catalog_ratio_equals_live(7) is False
    assert Qtau_is_gj_ratio() is False
    assert LIVE_QPP[5].denominator == 13
    assert LIVE_QPP[7].denominator == 409


def test_phi_F_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["409_absent_from_catalog"] is True
    assert out["proved"]["no_catalog_ratio_is_live_Q"] is True
    assert out["proved"]["Qtau_is_gj_ratio"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
