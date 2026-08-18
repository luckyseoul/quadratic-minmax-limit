"""Tests for Prop 15.499 — Aut_∞ T-orbits free; n_T=c_eq(p−1); scale not p−1."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15499 import main, prove_A, prove_B, prove_C, prove_open


def test_translations_free_orbit_p2():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"][5]["nT"] == 4
    assert A["by_p"][5]["T_sizes"] == [25]
    assert A["by_p"][5]["n_fixed"] == 0
    assert A["by_p"][7]["nT"] == 114
    assert A["by_p"][7]["T_sizes"] == [49]
    assert A["by_p"][7]["n_fixed"] == 0
    # fail-when-wrong: T-orbit size is not p
    assert A["by_p"][5]["T_sizes"] != [5]
    assert A["by_p"][7]["T_sizes"] != [7]


def test_nT_is_ceq_not_cmin():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"][5]["nT"] == c_eq_named(5) * 4
    assert B["by_p"][7]["nT"] == c_eq_named(7) * 6
    assert B["by_p"][7]["nT"] != c_min_named(7) * 6
    assert c_eq_named(7) != c_min_named(7)


def test_scale_orbits_not_all_pm1():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"][5]["sizes"] == [4]
    assert C["by_p"][5]["nS"] == 1
    assert 3 in C["by_p"][7]["sizes"]
    assert 6 in C["by_p"][7]["sizes"]
    assert C["by_p"][7]["nS"] != c_eq_named(7)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
