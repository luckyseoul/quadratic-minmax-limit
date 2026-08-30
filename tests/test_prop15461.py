"""Tests for Prop 15.461 — p=7 remainder shapes do not extend."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15451 import c_eq_ends
from e1_gmin_m4_prop15460 import c_CRRR
from e1_gmin_m4_prop15461 import (
    affJ2_k_set,
    k_J4sq_plus_J2,
    k_J6_plus_J2,
    k_aff_needed,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_J_shapes_only_r4():
    A = prove_A()
    assert A["proved"] is True
    assert k_J6_plus_J2(4) == comb(4, 2) == 6
    assert k_J4sq_plus_J2(4) == 6
    assert k_J6_plus_J2(7) == 9
    assert k_J6_plus_J2(7) != 21
    assert k_J4sq_plus_J2(7) != comb(7, 2)


def test_affJ2_k15_absent():
    B = prove_B()
    assert B["proved"] is True
    assert k_aff_needed(7) == 15
    ks = affJ2_k_set(13)
    assert 15 not in ks
    assert ks == set(range(4, 13))
    assert 3 in affJ2_k_set(7)


def test_naive_50_misses():
    C = prove_C()
    assert C["proved"] is True
    assert c_CRRR(7) + 15 == 19
    assert c_CRRR(13) + 15 == 50
    lo, hi, _ = c_eq_ends(13)
    assert not (lo <= 50 <= hi)


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.461"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["J_shapes_r4_only"] is True
    assert out["proved"]["affJ2_k15_absent"] is True
    assert out["proved"]["naive_50_misses"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
