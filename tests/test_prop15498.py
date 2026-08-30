"""Tests for Prop 15.498 — p≡3 J_{--} named; p=7 pairing identity."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15498 import (
    D_live,
    J_mm_named,
    J_mm_wrong_half,
    main,
    pair_p3_A2p,
    pair_p3_named,
    prove_A,
    prove_B,
    prove_open,
)


def test_J_mm_divisibility():
    A = prove_A()
    assert A["proved"] is True
    assert J_mm_named(5, 2) == 0
    assert J_mm_named(7, 8) == 4
    assert J_mm_named(7, 4) == -4
    assert J_mm_named(7, 2) == 0
    assert J_mm_wrong_half(7, 8) == 3
    assert A["by_p"][7]["n"] == 4
    assert A["by_p"][7]["values"] == [-4, 0, 4]


def test_pairing_identity_p7():
    B = prove_B()
    assert B["proved"] is True
    D = D_live(7)
    assert abs(float(pair_p3_named(7, 8, D)) - 200 / 409) < 1e-9
    assert abs(float(pair_p3_A2p(7, 8, D)) - float(pair_p3_named(7, 8, D))) > 1e-6


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
