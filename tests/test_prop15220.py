"""Tests for Prop 15.220 — m4₂↔K₄ identity; residual-(i) still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15217 import m4_target_R  # noqa: E402
from e1_gmin_m4_prop15220 import (  # noqa: E402
    K4_threshold_for_m4_bound,
    e4_closed_form,
    e4_from_s,
    hinge_status_220,
    m4_norm_from_K4,
    residual_i_closed_via_220,
    theorem_e4_formula,
    theorem_m4_K4_identity,
    theorem_m4_bound_iff_K4,
)


def test_e4_formula_proved():
    assert theorem_e4_formula()["proved"] is True


def test_e4_matches_generating_function():
    for n in (10, 26, 50):
        for s in range(-n, n + 1, 2):
            assert e4_from_s(n, s) == e4_closed_form(n, s)


def test_m4_K4_identity_proved():
    assert theorem_m4_K4_identity()["proved"] is True


def test_m4_iff_K4_proved():
    assert theorem_m4_bound_iff_K4()["proved"] is True
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        thr_m4 = m4_target_R(p)
        thr_k4 = K4_threshold_for_m4_bound(p)
        assert thr_k4 == n * (15 * n - 22)
        # invert identity
        # if m4 = thr_m4 then K4 = 24*thr_m4 + 9n² - 10n
        assert 24 * thr_m4 + 9 * n * n - 10 * n == thr_k4


def test_p5_identity_corroboration():
    """Census: at p=5, m4₂=1862/13 and K4=120400/13 match the bridge."""
    n = 26
    K4 = Fraction(120400, 13)
    m42 = m4_norm_from_K4(K4, n)
    assert m42 == Fraction(1862, 13)
    assert m42 <= m4_target_R(5)
    assert K4 <= K4_threshold_for_m4_bound(5)


def test_residual_i_still_open():
    assert residual_i_closed_via_220() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_220()
    assert h["e4_formula"] is True
    assert h["m4_K4_identity"] is True
    assert h["m4_iff_K4"] is True
    assert h["residual_i_closed_via_220"] is False
