"""Tests for Prop 15.459 — k≤3 undershoots; no p=13 slice in the window."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15459 import (
    main,
    prove_A,
    prove_B,
    prove_open,
    slice_count,
    sweep_p13,
)


def test_k3_undershoots():
    A = prove_A()
    assert A["proved"] is True
    assert slice_count(5, 1) == 1
    assert slice_count(7, 3) == 19
    assert slice_count(13, 3) == Fraction(3, 2)
    assert slice_count(13, 3) < c_eq_ends(13)[0]


def test_no_slice_in_window():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_in_win"] == 0
    assert Fraction(B["nearest_above"]["c"]) == 630
    assert Fraction(B["nearest_below"]["c"]) == Fraction(1011, 2)
    lo, hi, _ = c_eq_ends(13)
    assert not (lo <= 630 <= hi)
    # computed, not hard-coded: n0=2, max_k=6
    assert slice_count(13, 6, 2, 2) == 630
    rows = sweep_p13()
    assert all(not (lo <= Fraction(rw["c"]) <= hi) for rw in rows)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.459"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["k3_undershoots"] is True
    assert out["proved"]["no_slice_in_window"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
