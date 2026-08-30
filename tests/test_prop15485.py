"""Tests for Prop 15.485 — Aut_∞ orbits vs AG types at p=7."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15485 import fold_p7, main, prove_A, prove_B, prove_open


def test_orbit_table():
    A = prove_A()
    assert A["proved"] is True
    assert A["n_orbits"] == 8
    assert A["c_nl"] == 19
    f = fold_p7()
    ones = [r for r in f["recs"] if r["is_1d"]]
    assert sum(r["sz"] for r in ones) == n_1d(7) == 140
    assert sorted(r["sz"] for r in ones) == [56, 84]


def test_not_all_s2_are_TR():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_s2"] == 4
    assert B["n_s2_TR"] == 1
    assert B["n_s2_NF1"] == 3
    f = fold_p7()
    s2 = [r for r in f["recs"] if r["s"] == 2]
    assert not all(r["tags"].get("TR", 0) == r["sz"] for r in s2)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.485"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["orbit_table_p7"] is True
    assert out["proved"]["one_s2_is_TR"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
