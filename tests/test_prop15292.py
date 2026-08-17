"""Tests for Prop 15.292 — Q_1D named in p; Q_NL / floor still open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import (
    E_hat_pm1_formula,
    E_hat_product,
    E_hat_sub_formula,
    Q_1d_off_Fp_over_q2,
    Q_1d_pm1_over_q2,
    Q_1d_pm1_over_q2_inverted,
    Q_1d_sub_over_q2,
    Q_1d_sub_over_q2_inverted,
    S_1d_over_q2_ntriv_on_Fp,
    live_Q_1d,
    main,
    n_1d,
    prove_closed_form,
    prove_open,
)


def test_n_1d_counts():
    assert n_1d(5) == 30
    assert n_1d(7) == 140


def test_closed_form_matches_johnson():
    B = prove_closed_form([5, 7, 11, 13])
    assert B["proved"] is True
    assert B["invert_hit"] < B["invert_tot"]


def test_E_formulas_at_known_p():
    assert E_hat_product(5, 1) == E_hat_pm1_formula(5) == Fraction(7, 2)
    assert E_hat_product(5, 2) == E_hat_sub_formula(5) == Fraction(1)
    assert E_hat_product(7, 1) == Fraction(34, 5)
    assert E_hat_product(7, 2) == Fraction(13, 5)


def test_Q_formulas():
    assert Q_1d_pm1_over_q2(5) == Fraction(56, 3)
    assert Q_1d_sub_over_q2(5) == Fraction(16, 3)
    assert Q_1d_pm1_over_q2(7) == Fraction(136, 5)
    assert Q_1d_sub_over_q2(7) == Fraction(52, 5)
    assert Q_1d_off_Fp_over_q2(7) == 0


def test_invert_denom_breaks():
    assert Q_1d_pm1_over_q2_inverted(5) != Q_1d_pm1_over_q2(5)
    assert Q_1d_sub_over_q2_inverted(7) != Q_1d_sub_over_q2(7)


def test_live_lifts_match_formula():
    for p in (5, 7):
        Qrel = live_Q_1d(p)
        F = _kernel_field(p)
        for r, val in Qrel.items():
            if r in (1, F["neg1"]):
                target = float(Q_1d_pm1_over_q2(p))
            elif r < p:
                target = float(Q_1d_sub_over_q2(p))
            else:
                target = 0.0
            assert abs(val - target) < 1e-8


def test_invert_sub_zero_fails():
    Qrel = live_Q_1d(5)
    F = _kernel_field(5)
    sub = [v for r, v in Qrel.items() if r not in (1, F["neg1"]) and r < 5]
    assert max(abs(v) for v in sub) > 1e-4


def test_S_1d_above_six_not_the_floor():
    assert S_1d_over_q2_ntriv_on_Fp(5) == Fraction(4 * 5 * 4, 3)
    assert float(S_1d_over_q2_ntriv_on_Fp(7)) > 6


def test_Q_NL_and_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Q_1d_zero_off_Fp"] is True
    assert out["proved"]["Q_1d_closed_form"] is True
    assert out["proved"]["live_lifts_match"] is True
    assert out["proved"]["Q_NL_named_in_p"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
