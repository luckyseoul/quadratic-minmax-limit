"""Tests for Prop 15.363 — A008300 cut is full-spectrum; 13/204 dies."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15363 import (
    N7,
    T7,
    gen_T,
    main,
    occupancy,
    poly_alpha,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    wrap_diag_sums,
)


def test_fourier_identity_not_S0_over_p():
    A = prove_A()
    assert A["proved"] is True
    assert A["N_p5"] == 130
    assert A["T52"] == 2040
    assert A["fake_S0_over_p"] == 408.0


def test_affine_vanish_quad_cubic_occupancy():
    B = prove_B()
    assert B["proved"] is True
    assert B["n_T"] == 2040
    assert B["occ_x2"] == (440, 400, 400, 400, 400)
    assert B["occ_x3"] == (460, 395, 395, 395, 395)
    assert B["S_x2"] == 40
    mats = gen_T(5, 2)
    diags = [wrap_diag_sums(g, 1) for g in mats]
    occ_id = occupancy(diags, poly_alpha((0, 1), 5), 2)
    assert occ_id == (2040, 0, 0, 0, 0)


def test_p5_ratio_fails_at_p7():
    C = prove_C()
    assert C["proved"] is True
    assert C["remainder"] != 0
    assert T7 // 7 != N7
    assert T7 * 13 % 204 != 0


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["fourier_identity"] is True
    assert out["proved"]["S_on_deg_le_3"] is True
    assert out["proved"]["p5_ratio_fails_p7"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
