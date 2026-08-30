"""Tests for Prop 15.454 — complementary M-sum named; J2 not all O."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15454 import (
    J_r,
    M_all_named,
    M_comp_named,
    M_ns_named,
    initial_kset,
    is_Jr,
    main,
    one_d_halfnet,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_M_comp_identity():
    for p in (5, 7, 11, 13):
        assert M_comp_named(p) == M_all_named(p) - M_ns_named(p)
    assert M_comp_named(5) == 90
    assert M_comp_named(7) == 336
    A = prove_A()
    assert A["proved"] is True
    assert A["bad_p5"] != 90
    assert A["one_d_p5"] == 90
    assert len(one_d_halfnet(5)) == 10
    assert len(initial_kset(5)) == 10


def test_one_d_not_all():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_1d7"] == n_1d(7) == 140
    assert B["Hplus7"] == HPLUS[7] == 5726
    assert n_1d(7) != HPLUS[7]


def test_J2_not_residual():
    C = prove_C()
    assert C["proved"] is True
    assert J_r(7, (1, -1), (0, 1)) == (1, 1, 3, 3, 4, 4, 5)
    assert is_Jr(7, (0, 1, 1, 2, 5, 5, 7), rmax=4) is False
    assert is_Jr(7, (1, 1, 3, 3, 4, 4, 5), rmax=2) is True


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.454"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["M_comp_named"] is True
    assert out["proved"]["one_d_MH_tuple"] is True
    assert out["proved"]["J2_not_all_O"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
