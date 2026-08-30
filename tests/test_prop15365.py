"""Tests for Prop 15.365 — GJ neighbours of Q++ die at p=7."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15365 import (
    Q_neigh_dim,
    Q_neigh_wick,
    dim_Vp,
    main,
    p5_u_atoms,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_dim_neighbour_hits_p5_dies_p7():
    A = prove_A()
    assert A["proved"] is True
    assert Q_neigh_dim(5) == LIVE_QPP[5] == Fraction(48, 13)
    assert Q_neigh_dim(7) == Fraction(48, 25)
    assert Q_neigh_dim(7) != LIVE_QPP[7]
    assert Q_neigh_dim(5) == Fraction(8 * n_pp_named(5), dim_Vp(5))
    assert dim_Vp(7) == 25


def test_wick_neighbour_hits_p5_dies_p7():
    B = prove_B()
    assert B["proved"] is True
    assert Q_neigh_wick(5) == Fraction(48, 13)
    assert Q_neigh_wick(7) == Fraction(96, 25)
    assert Q_neigh_wick(7) != LIVE_QPP[7]


def test_catalog_misses_p7():
    C = prove_C()
    assert C["proved"] is True
    assert C["S7_has_409"] is False
    assert C["gj_multiples_of_409"] == []
    assert 409 not in C["fermat"].values()
    assert C["weil_E_range"] == [36, 64]
    assert 2401 != 409
    assert dim_Vp(7) != 409


def test_p5_u_atom_fails_at_p7():
    D = prove_D()
    assert D["proved"] is True
    assert D["p5_n_atoms"] == 5
    assert D["p5_exact"] is True
    assert D["p7_outsider"] is not None
    assert D["p7_gap"] > 1.0
    atoms = p5_u_atoms()
    assert abs(atoms[1] - 10 * (5 - 5**0.5)) < 1e-9


def test_floor_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert E["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["neigh_dim_dies_p7"] is True
    assert out["proved"]["neigh_wick_dies_p7"] is True
    assert out["proved"]["catalog_misses_1544_409"] is True
    assert out["proved"]["p5_u_5atom_not_p7"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
