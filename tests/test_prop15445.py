"""Tests for Prop 15.445 — complete sums name the ++/N* χ-block."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15445 import (
    M_pn_named,
    M_pp_named,
    ST_named,
    TT_named,
    chsum,
    dC_N,
    dC_pp,
    L_ratio_need,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    sub_tor,
    sum_T_chi,
)
from fractions import Fraction


def test_hyperbola_and_T_sum():
    F = _kernel_field(5)
    assert sum_T_chi(F, 2) == -1
    assert sum_T_chi(F, 2) != 0
    A = prove_A()
    assert A["proved"] is True


def test_blocks_and_fails():
    assert ST_named(5) == -8
    assert ST_named(5) != 0
    assert TT_named(5) == 4
    assert M_pp_named(5) == -10
    assert M_pn_named(5) == 0
    F = _kernel_field(5)
    S, T = sub_tor(F, 5)
    assert chsum(F, S, T)[0] == -8
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["ST"] == -8


def test_dC_and_ratio():
    assert dC_pp(5) == Fraction(25, 3)
    assert dC_N(5) == Fraction(-15)
    assert L_ratio_need(5) == Fraction(9, 5)
    C = prove_C()
    assert C["proved"] is True


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.445"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["T_sum_minus_one"] is True
    assert out["proved"]["type_block_matrix"] is True
    assert out["proved"]["ns_drops_dC_named"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["B"]["by_p"]["5"]["ST"] == -8
