"""Tests for Prop 15.298 — Q is the half-Gauss image of L_χ."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general
from e1_gmin_m4_prop15290 import live_Q_on_T, type_key
from e1_gmin_m4_prop15298 import (
    coarse_Q_class,
    main,
    prove_gauss_identity,
    prove_open,
    prove_Q_coarsening,
)


def test_gauss_identity_holds():
    A = prove_gauss_identity([5, 7])
    assert A["proved"] is True
    for b in A["blocks"]:
        assert b["max_err_over_q2"] < 1e-9
        assert b["naive_ns_err_over_q2"] > 1e-3


def test_naive_nonsquare_mu_breaks():
    A = prove_gauss_identity([5])
    assert A["blocks"][0]["naive_ns_err_over_q2"] > 0.1


def test_Q_coarsening_p7_merges():
    B = prove_Q_coarsening([5, 7])
    assert B["proved"] is True
    assert B["merges"]["pp_merge"] is True
    assert B["merges"]["Nstar_merge"] is True


def test_fail_when_wrong_distinct_off_types():
    """15.290 off types are not all distinct at p=7 (the coarsening)."""
    F = _kernel_field(7)
    Q = live_Q_on_T(7)
    q2 = float(F["q"] ** 2)
    by = {}
    for r, Qr in Q.items():
        by.setdefault(type_key(F, r), []).append(Qr / q2)
    # ++ sub and ++ N1 must match (otherwise coarsening is false)
    sub = sum(by[(1, 1, "sub")]) / len(by[(1, 1, "sub")])
    n1 = sum(by[(1, 1, "N1")]) / len(by[(1, 1, "N1")])
    assert abs(sub - n1) < 1e-6
    # claiming they differ is the invert
    assert abs(sub - n1) < 0.01


def test_coarse_class_pm1():
    for p in (5, 7):
        F = _kernel_field(p)
        assert coarse_Q_class(F, 1) == "pm1"
        assert coarse_Q_class(F, F["neg1"]) == "pm1"


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["gauss_identity"] is True
    assert out["proved"]["Q_coarsening_certified"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
