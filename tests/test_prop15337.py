"""Tests for Prop 15.337 — isolated orbit = harmonic c.s. family."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15337 import (
    K_HIST_ISOL,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    run_sweep,
)


def test_isol_khist():
    A = prove_A()
    assert A["proved"] is True
    assert A["n"] == 1176
    assert A["n_cs"] == 1176
    assert A["hists"][str(K_HIST_ISOL)] == 1176


def test_sweep_24_all_isolated():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_hits"] == 24
    assert B["n_isol_type"] == 24
    assert B["n_centered_times_q"] == 1176
    # fail-when-wrong: not the empty sweep
    assert B["n_hits"] != 0


def test_harmonic_opposite():
    hits = run_sweep()
    C = prove_C(hits)
    assert C["proved"] is True
    assert C["cross_ratio_0inf_34"] == 6
    assert C["n_full_slope1"] == 0
    assert C["n_opposite_FE"] == 24


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["isolated_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["isol_khist_3210"] is True
    assert out["proved"]["sweep_24_all_isol"] is True
    assert out["proved"]["harmonic_opposite_FE"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
