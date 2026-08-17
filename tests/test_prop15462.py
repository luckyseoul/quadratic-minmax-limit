"""Tests for Prop 15.462 — no single mixed + C(r,3) in the window."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15460 import c_CRRR
from e1_gmin_m4_prop15462 import (
    extra_named,
    main,
    mixed_rows,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    singles_in_window,
)


def test_c1_gap():
    A = prove_A()
    assert A["proved"] is True
    c1s = set(A["rows"]["13"]["c1s"])
    assert 516 not in c1s and 582 not in c1s
    assert not any(516 <= c <= 582 for c in c1s)
    assert A["rows"]["13"]["nearest_below"] == 455
    assert A["rows"]["13"]["nearest_above"] == 665
    lo, hi, _ = c_eq_ends(13)
    assert not (lo <= 665 <= hi)
    assert not (lo <= 455 <= hi)
    assert c_CRRR(13) + 630 == 665
    assert c_CRRR(13) + 420 == 455
    assert A["rows"]["17"]["nearest_below"] == 5124
    assert A["rows"]["17"]["nearest_above"] == 7644
    lo17, hi17, _ = c_eq_ends(17)
    assert not (lo17 <= 5124 <= hi17)


def test_extras_half_integer_only():
    B = prove_B()
    assert B["proved"] is True
    assert B["rec"]["13"]["flat"]["n"] == 0
    assert B["rec"]["13"]["half"]["n"] == 0
    assert B["rec"]["13"]["three_halves"]["n"] == 0
    assert B["rec"]["13"]["p2"]["cs"] == ["1225/2"]
    assert Fraction(1225, 2).denominator != 1
    lo, hi, _ = c_eq_ends(13)
    assert lo <= Fraction(1225, 2) <= hi
    hits = singles_in_window(13, "p2")
    assert hits
    assert all(h["c1"] == 105 for h in hits)
    assert extra_named(13, {"J4": 2, "J6": 4, "J10": 1}, "p2") == Fraction(11, 2)
    assert B["rec"]["17"]["flat"]["n"] == 0
    assert B["rec"]["17"]["p2"]["n"] == 0


def test_not_live_p7():
    C = prove_C()
    assert C["proved"] is True
    assert 16 in C["vals"]
    assert 19 not in C["vals"]
    assert C["hits"]["flat"] == 0
    c1s = {rw["c1"] for rw in mixed_rows(7)}
    assert 12 in c1s
    assert c_CRRR(7) + 12 == 16 != 19


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.462"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["c1_gap"] is True
    assert out["proved"]["extras_no_integer_hit"] is True
    assert out["proved"]["not_live_at_p7"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
