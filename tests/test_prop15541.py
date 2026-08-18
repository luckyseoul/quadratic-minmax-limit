"""Tests for Prop 15.541 — c_eq is the named Q_eq-Hoffman floor."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15457 import c_eq_named
from e1_gmin_m4_prop15541 import (
    P_CERT,
    c_eq_ceil,
    c_eq_drop_p2m9,
    c_eq_floor,
    main,
    prove_A,
    prove_B,
    prove_open,
    Q_hoffman,
)


def test_floor_equals_ceq_and_ceil_fails():
    A = prove_A()
    assert A["proved"] is True
    for p in P_CERT:
        assert c_eq_floor(p) == c_eq_named(p)
        assert c_eq_ceil(p) != c_eq_named(p)
    assert c_eq_drop_p2m9(7) != c_eq_named(7)
    assert c_eq_ceil(5) == 2
    assert c_eq_ceil(7) == 20
    assert c_eq_drop_p2m9(7) == 20


def test_Qhoff_hits_live_not_ceil():
    B = prove_B()
    assert B["proved"] is True
    assert Q_hoffman(5) == LIVE_QPP[5]
    assert Q_hoffman(7) == LIVE_QPP[7]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["live_c_eq_general"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
