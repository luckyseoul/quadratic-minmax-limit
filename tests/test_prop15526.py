"""Tests for Prop 15.526 — H_+ = n_1d + q n_free."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15526 import (
    D_drop_nfree,
    D_from_nfree,
    load_Hplus,
    main,
    n_free_wrong_pm1,
    n_free_wrong_zero,
    prove_A,
    prove_B,
    prove_open,
    translation_orbits,
)


def test_typeplus_1d_is_periodic_orbits():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        rec = A["by_p"][str(p)]
        H = load_Hplus(p)
        live = translation_orbits(H, p)
        assert live["nH"] == HPLUS[p]
        assert live["n_1d_in_H"] == n_1d(p)
        assert live["nH"] != n_1d(p)
        assert live["n_periodic_orbits"] * p == n_1d(p)
        assert live["n_free_orbits"] * p * p == HPLUS[p] - n_1d(p)
        assert rec["n_free_orbits"] == live["n_free_orbits"]


def test_D_formula_and_fail_nfree():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7):
        nf = B["by_p"][str(p)]["n_free"]
        assert D_from_nfree(p, nf) == D_named(p)
        assert nf != n_free_wrong_zero()
        assert D_drop_nfree(p) != D_named(p)
    assert B["by_p"]["7"]["n_free"] != n_free_wrong_pm1(7)
    assert B["by_p"]["5"]["n_free"] == n_free_wrong_pm1(5)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_free_named_in_p"] is False
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
