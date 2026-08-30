"""Tests for Prop 15.492 — per-D Var(e_*) is not V_floor/κ."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15326 import V_floor, kappa_named
from e1_gmin_m4_prop15492 import main, prove_A, prove_B, prove_open


def test_p5_1d_exceeds_floor_nl_does_not():
    A = prove_A()
    assert A["proved"] is True
    row = A["row"]
    assert row["n_1d_flag"] == n_1d(5) == 30
    assert row["by"]["1d=True|Var=5.0"] == 30
    assert row["by"]["1d=False|Var=1.8"] == 100
    need = float(V_floor(5) / kappa_named(5))
    assert need == 18 / 7
    assert 5.0 > need
    assert 1.8 < need
    # fail-when-wrong: a single value, or every design under the cut
    assert len(row["by"]) == 2


def test_p7_1d_not_a_single_var():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_1d_flag"] == n_1d(7) == 140
    assert len(B["vars_1d"]) >= 2


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
