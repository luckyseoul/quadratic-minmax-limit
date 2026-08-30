"""Tests for Prop 15.471 — lag-s quartic weight; order-4 pair halves."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_plus
from e1_gmin_m4_prop15467 import coarse_L
from e1_gmin_m4_prop15468 import chi_minus_named, chi_plus_named, k_named, S_chi_named
from e1_gmin_m4_prop15471 import (
    LIVE_LAM1_SQ,
    Lambda_tables,
    main,
    mean_on_chi,
    pair_n_mod4,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_nk_halves():
    A = prove_A()
    assert A["proved"] is True
    cnt = pair_n_mod4(5)
    assert abs(cnt[:, 0].mean() - float(chi_plus_named(5) / 2)) < 1e-8
    assert cnt[:, 1].max() - cnt[:, 1].min() >= 1
    cnt7 = pair_n_mod4(7)
    assert abs(cnt7[:, 1].max() - cnt7[:, 1].min()) < 1e-8
    assert abs(cnt7[:, 1].mean() - float(chi_minus_named(7) / 2)) < 1e-8


def test_Lam02_named():
    B = prove_B()
    assert B["proved"] is True
    F, _, Lam = Lambda_tables(5)
    k = k_named(5)
    assert abs(mean_on_chi(F, Lam[0], 1).real - float(mu_plus(5) * k * (k - 1))) < 1e-6
    assert abs(mean_on_chi(F, Lam[2], 1).real - float(mu_plus(5) * S_chi_named(5))) < 1e-6
    assert abs(mean_on_chi(F, Lam[2], 1).real) > 1.0


def test_Lam1_not_cov():
    C = prove_C()
    assert C["proved"] is True
    F, _, Lam = Lambda_tables(5)
    v = mean_on_chi(F, Lam[1], 1)
    assert abs(v - LIVE_LAM1_SQ[5]) < 1e-6
    cov = coarse_L(5)["++"][1] - mu_plus(5) ** 2
    assert abs(v.real - float(cov)) > 1.0


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.471"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["nk_halves"] is True
    assert out["proved"]["Lam1_not_cov"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
