"""Tests for Prop 15.463 — no two-occupancy support pair hits 1,19 and the window."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15463 import (
    P7_FIFTEEN,
    extra_flat,
    extra_p7,
    main,
    n_by_support,
    pair_c,
    prove_A,
    prove_B,
    prove_open,
)


def test_flat_no_19():
    A = prove_A()
    assert A["proved"] is True
    assert "19" not in A["c7_vals"]
    assert "20" in A["c7_vals"]
    assert "14" in A["c7_vals"]


def test_p7_five_miss_window():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_fifteen"] == 5
    assert all(r["c7"] == "19" for r in B["rows"])
    assert all(r["in_win"] is False for r in B["rows"])
    lo, hi, _ = c_eq_ends(13)
    c_bad = pair_c(13, frozenset({"J6"}), frozenset({"J2", "J4", "R"}), extra_p7)
    assert not (lo <= c_bad <= hi)
    assert pair_c(7, frozenset({"J6"}), frozenset({"J2", "J4", "R"}), extra_p7) == 19
    assert pair_c(5, frozenset({"J6"}), frozenset({"J2", "J4", "R"}), extra_p7) == 1
    got = {frozenset((frozenset(r["S"]), frozenset(r["T"]))) for r in B["rows"]}
    want = {frozenset(pt) for pt in P7_FIFTEEN}
    assert got == want


def test_flat_pair_not_fifteen():
    n7 = n_by_support(7, extra_flat)
    n5 = n_by_support(5, extra_flat)
    S7 = [S for S in n7 if S not in n5]
    assert all(n7[S] + n7[T] != 15 for i, S in enumerate(S7) for T in S7[i + 1 :])


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.463"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["flat_no_19"] is True
    assert out["proved"]["p7_five_miss_window"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
