"""Tests for Prop 15.511 — 1D involution Fix; D / Q_τ still unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15511 import (
    main,
    n_1d_inv_named,
    n_1d_inv_wrong_drop_m,
    n_1d_inv_wrong_n1dp,
    n_even_m_subsets_named,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_even_subsets_are_binom():
    A = prove_A()
    assert A["proved"] is True
    assert n_even_m_subsets_named(5) == 2
    assert n_even_m_subsets_named(7) == 3
    assert A["rows"]["13"]["live"] == 20


def test_named_product_and_fails():
    B = prove_B()
    assert B["proved"] is True
    assert n_1d_inv_named(3) == 2
    assert n_1d_inv_named(5) == 6
    assert n_1d_inv_named(7) == 12
    assert n_1d_inv_named(11) == 60
    assert n_1d_inv_wrong_n1dp(7) == 20
    assert n_1d_inv_named(7) != n_1d_inv_wrong_n1dp(7)
    assert n_1d_inv_wrong_drop_m(7) == 3
    assert n_1d_inv_named(7) != n_1d_inv_wrong_drop_m(7)


def test_cache_1d_matches_named():
    C = prove_C()
    assert C["proved"] is True
    assert C["rows"]["5"]["n_1d_live"] == 6
    assert C["rows"]["7"]["n_1d_live"] == 12
    assert C["rows"]["7"]["nl_left"] == 42
    assert C["rows"]["7"]["full_fix"] == 54


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert D["D_named_in_p"] is False
    assert D["nl_leftover"]["7"] == 42
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
