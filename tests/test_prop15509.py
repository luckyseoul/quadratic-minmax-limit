"""Tests for Prop 15.509 — Burnside −id class; D / Q_τ still unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15509 import (
    class_size_named,
    class_size_wrong,
    fix_named,
    fix_wrong,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_class_size_is_q():
    A = prove_A()
    assert A["proved"] is True
    assert class_size_named(5) == 25
    assert class_size_named(7) == 49
    assert class_size_named(5) != class_size_wrong(5)
    assert class_size_named(7) != class_size_wrong(7)


def test_cycle_type_named():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["n_ok"] == 25
    assert B["rows"]["7"]["n_ok"] == 49


def test_fix_is_two_three_binom():
    C = prove_C()
    assert C["proved"] is True
    assert fix_named(3) == 2
    assert fix_named(5) == 6
    assert fix_named(7) == 54
    assert fix_wrong(5) == 6
    assert fix_wrong(7) == 20
    assert fix_named(7) != fix_wrong(7)
    assert C["rows"]["7"]["live_min"] == 54
    assert C["rows"]["5"]["wrong"] == 6


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["D_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert D["p11_predicted"] == 118098
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
