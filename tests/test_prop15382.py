"""Tests for Prop 15.382 — interior 4-level ≠ pure-pair double-star."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15382 import (
    E_fg_U_box,
    E_fg_U_boxtimes,
    main,
    n_box_muf0,
    n_boxtimes_muf0,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    pure_pair_S_mod4,
)


def test_interior_needs_both_mod4():
    A = prove_A()
    assert A["proved"] is True
    assert A["m4_mod"] == 0
    assert A["m6_mod"] == 2
    assert A["m4_mod"] != A["m6_mod"]


def test_pure_pair_not_interior():
    B = prove_B()
    assert B["proved"] is True
    assert pure_pair_S_mod4(5) == 2
    assert pure_pair_S_mod4(7) == 0
    assert pure_pair_S_mod4(5) != 0
    # one residue cannot be both 0 and 2
    assert not (pure_pair_S_mod4(5) == 0 and pure_pair_S_mod4(5) == 2)


def test_n_boxtimes_is_p_plus_1():
    C = prove_C()
    assert C["proved"] is True
    assert n_boxtimes_muf0(5) == 6
    assert n_box_muf0(5) == 14
    assert n_boxtimes_muf0(5) != n_box_muf0(5)
    assert E_fg_U_box(5) == 0
    assert E_fg_U_boxtimes(5) == Fraction(-2, 6)
    assert n_boxtimes_muf0(5) * E_fg_U_boxtimes(5) == -2


def test_residual_ii_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False
    assert residual_ii_k_eq_4p_empty() is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False


def test_main():
    out = main()
    assert out["proved"]["interior_needs_both_mod4"] is True
    assert out["proved"]["pure_pair_not_interior"] is True
    assert out["proved"]["n_boxtimes_eq_p_plus_1"] is True
    assert out["proved"]["residual_ii_k_eq_4p_empty"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
