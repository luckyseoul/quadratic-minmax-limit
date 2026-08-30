"""Tests for Prop 15.542 — ns half-nets = H_+ at p=5,7."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15539 import n_free_from_c, n_free_live
from e1_gmin_m4_prop15542 import (
    HALFNET,
    count_halfnets,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_halfnet_equals_Hplus():
    A = prove_A()
    assert A["proved"] is True
    assert HALFNET[5] == 130 == HPLUS[5] == count_halfnets(5)
    assert HALFNET[7] == 5726 == HPLUS[7]
    assert A["by_p"]["7"]["n_H_halfnet"] == HPLUS[7]
    assert count_halfnets(5) != n_1d(5)
    assert HALFNET[7] != n_1d(7)
    assert HALFNET[7] != HPLUS[7] + 1


def test_nfree_from_hn_is_ceq_not_cmin():
    B = prove_B()
    assert B["proved"] is True
    assert n_free_live(5) == 4 == n_free_from_c(5, c_eq_named(5))
    assert n_free_live(7) == 114 == n_free_from_c(7, c_eq_named(7))
    assert n_free_from_c(7, c_min_named(7)) == 108 != 114


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_free_general"] is False
    assert C["halfnet_eq_H_general"] is False
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
