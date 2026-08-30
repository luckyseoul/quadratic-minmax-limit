"""Tests for Prop 15.518 — p=11 affine leftover stab=2p; p=7 has no 2p."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15308 import AP_stab
from e1_gmin_m4_prop15309 import QR0_stab
from e1_gmin_m4_prop15518 import (
    leftover_stab_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_p7_no_twop_leftover():
    A = prove_A()
    assert A["proved"] is True
    assert leftover_stab_named(7) == 14
    assert QR0_stab(7) == 42
    assert leftover_stab_named(7) != QR0_stab(7)
    assert leftover_stab_named(7) not in A["fold"]["summary"]["QR0"]["stabs"]


def test_p11_leftover_has_twop_and_fourp():
    B = prove_B()
    assert B["proved"] is True
    other = B["fold"]["summary"]["other"]
    assert leftover_stab_named(11) in other["stabs"]
    assert AP_stab(11) in other["stabs"]
    assert other["stabs"] != [leftover_stab_named(11)]


def test_four_orbits_do_not_exhaust_n1d():
    C = prove_C()
    assert C["proved"] is True
    assert C["distinct_reps"] is True
    assert C["affine_sum_if_four"] == 1452
    assert C["n_1d"] == n_1d(11) == 2772
    assert C["affine_sum_if_four"] != n_1d(11)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
