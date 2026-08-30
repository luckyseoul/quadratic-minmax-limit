"""Tests for Prop 15.451 — Q_eq cut; p=5 half-net IFF; Type−."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15450 import c_window_ends
from e1_gmin_m4_prop15451 import (
    c_eq_ends,
    count_halfnets_p5,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    type_minus_halfnet_count,
)


def test_qeq_cut_tighter_than_qlo():
    cmin, cmax, nlo = c_window_ends(13)
    emin, emax, neq = c_eq_ends(13)
    assert emin == cmin == 551
    assert emax == 617
    assert cmax == 642
    assert neq == 67
    assert neq < nlo
    assert emax < cmax
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["n_eq"] == 1
    assert A["by_p"]["13"]["c_eq"] == 617


def test_p5_halfnet_count_is_hplus():
    n = count_halfnets_p5()
    assert n == HPLUS[5] == 130
    assert n != n_1d(5)
    B = prove_B()
    assert B["proved"] is True
    assert B["n_halfnet"] == 130


def test_type_minus_not_halfnet():
    assert type_minus_halfnet_count(5) == 0
    assert type_minus_halfnet_count(7) == 0
    C = prove_C()
    assert C["proved"] is True


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["n_eq_13"] == 67
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.451"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["c_eq_cut"] is True
    assert out["proved"]["p5_halfnet_iff"] is True
    assert out["proved"]["type_minus_not_halfnet"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["A"]["by_p"]["13"]["n_eq"] == 67
    assert out["algebra"]["B"]["n_halfnet"] == 130
