"""Tests for Prop 15.522 — Q_τ·D is a Z[√-2] norm; Z[i] fails."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15522 import (
    D_named,
    GAUSS_D,
    RING_D,
    a_pp_interpolant,
    all_off_types_are_norm,
    main,
    off_type_nums,
    prove_A,
    prove_B,
    prove_open,
    representations,
)


def test_off_types_are_zsqrt_minus2_norms():
    A = prove_A()
    assert A["proved"] is True
    for p in (5, 7):
        assert all_off_types_are_norm(p, RING_D) is True
        nums = off_type_nums(p)
        assert nums
        want_pp = int(LIVE_QPP[p] * D_named(p))
        assert want_pp in nums.values()
        for n in nums.values():
            assert n > 0
            assert representations(n, RING_D)


def test_gaussian_d1_fails_p5_plusplus():
    B = prove_B()
    assert B["proved"] is True
    assert all_off_types_are_norm(5, GAUSS_D) is False
    pp = int(LIVE_QPP[5] * D_named(5))
    assert representations(pp, GAUSS_D) == []
    assert B["p5_gauss_covers"] is False


def test_interpolant_not_imported_as_law():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["ab_named_in_p"] is False
    assert C["interpolant_imported"] is False
    assert C["phi_F_imported"] is False
    # 10p-46 hits a_++ at the two live primes and is still not a name
    a5 = representations(int(LIVE_QPP[5] * D_named(5)), RING_D)[0][0]
    a7 = representations(int(LIVE_QPP[7] * D_named(7)), RING_D)[0][0]
    assert a_pp_interpolant(5) == a5
    assert a_pp_interpolant(7) == a7
    assert C["interpolant_10p_46_hits_a_pp"] is True


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
