"""Tests for Prop 15.527 — n_free = n_R + n_X."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15526 import D_from_nfree
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15527 import (
    main,
    n_R_wrong_pm1,
    prove_A,
    prove_B,
    prove_open,
    split_free_orbits,
)


def test_split_matches_live_Hplus():
    A = prove_A()
    assert A["proved"] is True
    for p, wr, wx in ((5, 4, 0), (7, 16, 98)):
        rec = split_free_orbits(p)
        assert rec["n_R"] == wr
        assert rec["n_X"] == wx
        assert rec["n_R"] + rec["n_X"] == rec["n_free"]
        assert rec["n_reps"] == rec["n_free"]
        assert D_from_nfree(p, rec["n_free"]) == D_named(p)


def test_fail_nR_pm1_and_nX_zero_at_p7():
    B = prove_B()
    assert B["proved"] is True
    assert split_free_orbits(7)["n_R"] != n_R_wrong_pm1(7)
    assert split_free_orbits(7)["n_X"] != 0
    assert split_free_orbits(5)["n_R"] == n_R_wrong_pm1(5)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_free_named_in_p"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
