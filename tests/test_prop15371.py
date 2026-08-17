"""Tests for Prop 15.371 — D form at p=5,7; N4 split 2k^2."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15357 import D_live
from e1_gmin_m4_prop15367 import A_num, k_named
from e1_gmin_m4_prop15368 import dim_Vp
from e1_gmin_m4_prop15370 import N4_H, N4_O
from e1_gmin_m4_prop15371 import (
    D_form,
    D_form_drop_binom,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_D_form_hits_live_and_fails_when_wrong():
    A = prove_A()
    assert A["proved"] is True
    assert D_form(5) == D_live(5) == 13
    assert D_form(7) == D_live(7) == 409
    # formula is computed from named atoms, not a hardcoded table
    assert D_form(5) == 2 * A_num(5) - 3 * 3 + comb(5, 2)
    assert D_form(7) == 2 * A_num(7) - 3 * 4 + comb(7, 3)
    # fail-when-wrong
    assert D_form_drop_binom(5) == 3 != 13
    assert D_form_drop_binom(7) == 374 != 409
    assert D_form(3) == 3 != 1
    assert dim_Vp(7) != 409


def test_N4_split_is_two_k2():
    B = prove_B()
    assert B["proved"] is True
    assert N4_H - N4_O == 2 * k_named(7) ** 2
    assert N4_H - N4_O != k_named(7) ** 2
    assert B["p5_N4"] == 20
    assert B["p5_N4_constant"] is True


def test_D3_still_unnamed():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_hits_D3"] == 0


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False
    assert D["D_form_equals_live_p5_p7"] is True


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["D_form_p5_p7"] is True
    assert out["proved"]["N4_split_is_two_k2"] is True
    assert out["proved"]["D3_not_named"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
