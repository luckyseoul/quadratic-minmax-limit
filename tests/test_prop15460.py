"""Tests for Prop 15.460 — k_R identity; CRRR = C(r,3) u; 35 misses window."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15455 import n_CRRR_named
from e1_gmin_m4_prop15456 import u_unit
from e1_gmin_m4_prop15460 import (
    c_CRRR,
    count_crrr_all_slots,
    k_R_from_sig,
    k_R_named,
    main,
    n_CRRR_binom,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    r_of,
)


def test_k_R_identity():
    A = prove_A()
    assert A["proved"] is True
    assert k_R_named(5) == 1
    assert k_R_named(7) == 2
    assert k_R_named(13) == 7
    assert k_R_named(7) != 1
    for p in (5, 7, 13, 17):
        assert k_R_named(p) == k_R_from_sig(p)
        assert 3 * k_R_named(p) == comb(r_of(p), 2)


def test_CRRR_binom_count():
    B = prove_B()
    assert B["proved"] is True
    assert n_CRRR_binom(5) == 100
    assert n_CRRR_binom(7) == 1176
    assert n_CRRR_binom(13) == 70980
    assert count_crrr_all_slots(5) == 100
    assert count_crrr_all_slots(7) == 1176
    assert n_CRRR_named(5) == 300
    assert n_CRRR_named(5) != n_CRRR_binom(5)
    assert comb(r_of(7), 4) * u_unit(7) != 1176


def test_c_misses_window():
    C = prove_C()
    assert C["proved"] is True
    assert c_CRRR(5) == 1
    assert c_CRRR(7) == 4
    assert c_CRRR(13) == 35
    lo, hi, _ = c_eq_ends(13)
    assert not (lo <= 35 <= hi)
    assert c_CRRR(7) != 19


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.460"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["k_R_identity"] is True
    assert out["proved"]["n_CRRR_binom"] is True
    assert out["proved"]["c_CRRR_misses_window"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
