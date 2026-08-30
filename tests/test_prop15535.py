"""Tests for Prop 15.535 — twisted σ free at p=7, not p=5."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15532 import DIRS
from e1_gmin_m4_prop15535 import (
    a_pm,
    b_solutions,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_field_N_and_b_line():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        assert A["by_p"][str(p)]["N_plus"] == 1
        assert A["by_p"][str(p)]["N_minus"] == 1
        assert A["by_p"][str(p)]["n_b_plus"] == p
        assert A["by_p"][str(p)]["involution"] is True


def test_a_pm_norm_one_explicit():
    F = _kernel_field(5)
    va, vb = DIRS[5]["ns"]
    v = va + vb * 5
    ap, am = a_pm(F, v)
    assert b_solutions(F, ap)
    assert len(b_solutions(F, ap)) == 5


def test_free_at_p7_not_p5():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["n_free"] == 0
    assert B["by_p"]["5"]["n_maps"] == 10
    assert B["by_p"]["5"]["fix_set"] == [2, 30]
    assert B["by_p"]["5"]["inv_set"] == [2, 6]
    assert B["by_p"]["7"]["n_free"] == 14
    assert B["by_p"]["7"]["n_maps"] == 14
    assert B["by_p"]["7"]["fix_set"] == [0]
    assert B["by_p"]["7"]["inv_set"] == [0]


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["D_named_as_orbits"] is False
    assert C["free_pairing_p5"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
