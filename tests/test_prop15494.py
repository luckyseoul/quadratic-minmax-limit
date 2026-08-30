"""Tests for Prop 15.494 — A_r vanishes iff χ(r+1)=−1."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15290 import omega_set, paley_pair
from e1_gmin_m4_prop15494 import (
    main,
    prove_A,
    prove_B,
    prove_open,
    support_vanishes,
)


def test_support_chi_rp1_iff_outside_omega():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7, 11):
        row = A["by_p"][p]
        assert row["bad_m1"] == 0
        assert row["bad_p1"] == 0
        assert row["n_chi_rp1_m1"] > 0
        assert row["n_chi_rp1_p1"] > 0
        F = _kernel_field(p)
        Omega = omega_set(F)
        Oset = set(Omega)
        xi = Omega[0]
        # fail-when-wrong: invert to χ(r−1)
        inverted_ok = 0
        n = 0
        for r in F["squares"]:
            if r in (1, F["neg1"]):
                continue
            ch_rm1, ch_rp1 = paley_pair(F, r)
            n += 1
            if ch_rm1 == -1 and support_vanishes(F, r, xi, Oset):
                inverted_ok += 1
        # χ(r−1)=−1 is not the vanishing criterion
        assert inverted_ok < n


def test_live_A3_vanishes_only_on_chi_rp1_m1():
    B = prove_B()
    assert B["proved"] is True
    r5 = B["by_p"][5]
    assert abs(r5["(-1, -1, 'N*')"]["mean"]) < 1e-8
    assert abs(r5["(1, -1, 'N1')"]["mean"]) < 1e-8
    # fail: vanish on ++ 
    assert abs(r5["(1, 1, 'sub')"]["mean"] - 125 / 26) < 1e-6
    assert abs(r5["(-1, 1, 'N1')"]["mean"] - 625 / 26) < 1e-6
    r7 = B["by_p"][7]
    for tk, rec in r7.items():
        if rec["chi_rp1"] == -1:
            assert abs(rec["mean"]) < 1e-8
            assert rec["spread"] < 1e-8


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert D["Q_tau_named"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["phi_F_ge_6"] is False
