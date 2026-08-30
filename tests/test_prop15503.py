"""Tests for Prop 15.503 — per-y Q constant on sub, not on coarse ++."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15503 import main, prove_A, prove_B, prove_open


def test_per_y_constant_on_sub_not_plusplus():
    A = prove_A()
    assert A["proved"] is True
    r = A["by_p"][7]
    assert r["sub_max_spread"] < 1e-8
    assert r["pp_n_y_spread"] > 0
    assert r["pp_max_spread"] > 1.0
    assert set(r["sub"]) == {2, 3, 4, 5}


def test_ensemble_still_flat():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"][7]["++"]["ens_spread"] < 1e-8
    assert B["by_p"][5]["++"]["ens_spread"] < 1e-8


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
