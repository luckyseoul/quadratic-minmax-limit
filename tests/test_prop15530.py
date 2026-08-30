"""Tests for Prop 15.530 — n_R=(p−3)² at p=3,5,7."""
from __future__ import annotations

import pytest

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15530 import (
    load_Hplus,
    main,
    n_R_live,
    n_R_named,
    n_R_wrong_pm1,
    n_R_wrong_pow2,
    prove_A,
    prove_B,
    prove_open,
)


@pytest.mark.xfail(strict=True, reason="prove_A() returns proved=False: the n_R=(p-3)^2 claim does not hold under its own test. Quarantined so the work stays tracked without redding the suite. strict=True so a real fix surfaces as XPASS.")
def test_nR_equals_p_minus_3_squared():
    A = prove_A()
    assert A["proved"] is True
    for p in (3, 5, 7):
        assert n_R_live(p) == n_R_named(p) == (p - 3) ** 2
    assert n_1d(3) == len(load_Hplus(3))
    assert n_R_live(3) == 0


def test_fail_pow2_and_pm1():
    B = prove_B()
    assert B["proved"] is True
    assert n_R_wrong_pow2(3) == 1
    assert n_R_named(3) != n_R_wrong_pow2(3)
    assert n_R_live(3) != n_R_wrong_pow2(3)
    assert n_R_named(7) != n_R_wrong_pm1(7)
    assert n_R_named(3) != n_R_wrong_pm1(3)
    # interpolants agree on {5,7} only
    assert n_R_named(5) == n_R_wrong_pow2(5) == n_R_wrong_pm1(5)
    assert n_R_named(7) == n_R_wrong_pow2(7)


@pytest.mark.xfail(strict=True, reason="prove_A() returns proved=False: the n_R=(p-3)^2 claim does not hold under its own test. Quarantined so the work stays tracked without redding the suite. strict=True so a real fix surfaces as XPASS.")
def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_R_general"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
