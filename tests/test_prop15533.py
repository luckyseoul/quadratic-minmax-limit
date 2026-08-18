"""Tests for Prop 15.533 — Aut_∞ involutions do not freely pair T_ns-orbits."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15533 import (
    main,
    n_inv_neg_named,
    n_inv_wrong_1d_orb,
    n_inv_wrong_zero,
    prove_A,
    prove_B,
    prove_open,
)


def test_n_inv_neg_is_1509_fix():
    assert n_inv_neg_named(5) == 6
    assert n_inv_neg_named(7) == 54
    assert n_inv_neg_named(5) != n_inv_wrong_zero()
    assert n_inv_neg_named(7) != n_inv_wrong_1d_orb(7)


def test_live_neg_matches_named():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["n_inv_neg"] == 6
    assert A["by_p"]["7"]["n_inv_neg"] == 54
    assert A["by_p"]["5"]["n_inv_g1"] == 6
    assert A["by_p"]["7"]["n_inv_g1"] == 54


def test_frob_not_free_pairing():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["frob"] == 2
    assert B["by_p"]["7"]["frob"] == 218
    assert B["by_p"]["5"]["frob_neg"] == 22
    assert B["by_p"]["7"]["frob_neg"] == 218
    assert B["by_p"]["5"]["frob"] != 0
    assert B["by_p"]["7"]["frob"] != 0


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["D_named_as_orbits"] is False
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
