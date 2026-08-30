"""Tests for Prop 15.487 — Paley masses and mixed type-indicator."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15487 import (
    LIVE_QNSTAR,
    fold_pm,
    main,
    n_N1_named,
    n_Nstar_half,
    n_Nstar_named,
    n_mm_named,
    n_pm_drop_aE,
    n_pm_named,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_paley_masses():
    assert n_pm_named(5) == 2
    assert n_pm_named(7) == 4
    assert n_mm_named(5) == 4
    assert n_mm_named(7) == 8
    assert n_pm_drop_aE(5) == Fraction(3, 2)
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["5"]["+-"] == 2
    assert A["rows"]["11"]["+-"] == 12


def test_norm_masses():
    assert n_N1_named(5) == 4
    assert n_N1_named(5) != 3
    assert n_Nstar_named(7) == 12
    assert n_Nstar_half(7) == 6
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["7"]["N*"] == 12


def test_mixed_pairing_not_Qpp():
    C = prove_C()
    assert C["proved"] is True
    rec5, rec7 = fold_pm(5), fold_pm(7)
    assert abs(rec5["pred"] - float(LIVE_QPP[5])) < 1e-9
    assert abs(rec7["pred"] - float(LIVE_QNSTAR[7])) < 1e-9
    assert abs(rec7["pred"] - float(LIVE_QPP[7])) > 0.01


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.487"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["paley_masses_CM"] is True
    assert out["proved"]["mixed_I_not_Qpp"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
