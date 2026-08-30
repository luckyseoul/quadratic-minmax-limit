"""Tests for Prop 15.353 — ensemble Q++ is a 2-/4-level mix."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15353 import main, mix_of, prove_A, prove_B, prove_open, twotype_Q


def test_p5_two_values_mix_to_live():
    A = prove_A()
    assert A["proved"] is True
    assert set(A["hist"]) == {"16/9", "64/15"}
    assert A["hist"]["16/9"] == n_1d(5)
    assert mix_of(A["hist"]) == LIVE_QPP[5]


def test_p7_four_values_and_2type_misses():
    B = prove_B()
    assert B["proved"] is True
    assert B["mix"] == "1544/409"
    assert twotype_Q(7) == Fraction(1796, 409)
    assert twotype_Q(7) != LIVE_QPP[7]


def test_twotype_hits_p5_only():
    assert twotype_Q(5) == LIVE_QPP[5]
    assert twotype_Q(7) != LIVE_QPP[7]


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["p5_two_level_mix"] is True
    assert out["proved"]["p7_four_level_mix_not_2type"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
