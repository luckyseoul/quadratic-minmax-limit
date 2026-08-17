"""Tests for Prop 15.392 — first_nc t=0; c nsq; ib dies."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15392 import (
    chi_p,
    good_c_at_0,
    is_hoffman_c,
    least_nsq,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
    y_triangle_p1,
)


def test_t_eq_0():
    A = prove_A()
    assert A["proved"] is True
    assert all(A["rows"][str(p)]["t"] == 0 for p in (5, 13, 17, 29))


def test_c_nonsquare_min():
    B = prove_B()
    assert B["proved"] is True
    goods = good_c_at_0(5)
    assert set(goods) == {2, 3}
    assert all(chi_p(c, 5) == -1 for c in goods)
    assert not is_hoffman_c(5, y_triangle_p1(5), 1, t=0)


def test_ib_dies():
    C = prove_C()
    assert C["proved"] is True
    assert C["rows"]["13"]["min_c"] == 5
    assert C["rows"]["13"]["ib"] == 2
    assert least_nsq(13) == 2
    assert C["rows"]["13"]["eq"] is False


def test_catalog_empty():
    D = prove_D()
    assert D["proved"] is True
    assert D["eq_hits"] == []
    assert D["live"][13] == 5


def test_phi_F_still_open():
    E = prove_open()
    assert E["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["t_eq_0"] is True
    assert out["proved"]["c_nonsquare"] is True
    assert out["proved"]["ib_dies"] is True
    assert out["proved"]["catalog_empty"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
