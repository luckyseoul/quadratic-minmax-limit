"""Tests for Prop 15.539 — n_free=c(p−1) ⇔ D=D_lattice(c)."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15409 import D_lattice
from e1_gmin_m4_prop15450 import u_from_c, u_from_c_wrong
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15539 import (
    D_from_nfree,
    main,
    n_free_from_c,
    n_free_live,
    prove_A,
    prove_B,
    prove_open,
)


def test_algebra_nfree_iff_Dlattice():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11, 13):
        ce = c_eq_named(p)
        nf = n_free_from_c(p, ce)
        assert abs(D_from_nfree(p, nf) - D_lattice(p, u_from_c(p, ce))) < 1e-12
        assert D_lattice(p, u_from_c_wrong(p, ce)) != D_lattice(p, u_from_c(p, ce))


def test_live_pin_and_fail_cmin():
    B = prove_B()
    assert B["proved"] is True
    assert n_free_live(5) == 4 == n_free_from_c(5, c_eq_named(5))
    assert n_free_live(7) == 114 == n_free_from_c(7, c_eq_named(7))
    assert n_free_from_c(7, c_min_named(7)) == 108
    assert n_free_from_c(7, c_min_named(7)) != 114
    assert D_named(7) == D_lattice(7, u_from_c(7, c_eq_named(7)))


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["n_free_general"] is False
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
