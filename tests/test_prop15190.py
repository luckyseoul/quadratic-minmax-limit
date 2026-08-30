"""Tests for Prop 15.190 — scheme-ker; ker-box analysis; residual i open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15190 import (  # noqa: E402
    alpha_particular,
    certify_kerbox,
    certify_three_point_vanish,
    scheme_ker_max_kappa_e,
    theorem_scheme_ker_dual_eq_blocked,
    theorem_scheme_ker_in_gsum_ker,
)


def test_scheme_ker_theorems():
    assert theorem_scheme_ker_in_gsum_ker()["proved"] is True
    r = theorem_scheme_ker_dual_eq_blocked()
    assert r["proved"] is True
    for p in (3, 5, 7, 11):
        mx = scheme_ker_max_kappa_e(p)
        need = 2 - alpha_particular(p)
        assert mx < need
        # closed form
        n = p * p + 1
        assert mx == Fraction(3 * (n - 2), p * n)


def test_three_point_vanish_certified():
    for p in (3, 5):
        r = certify_three_point_vanish(p)
        assert r["vanishes"] is True


def test_kerbox_census():
    r3 = certify_kerbox(3)
    assert r3["dual_eq_feasible"] is True  # p=3 not blocked by full ker
    assert r3["blocked"] is False
    r5 = certify_kerbox(5)
    assert r5["dual_eq_feasible"] is False
    assert r5["blocked"] is True
    assert r5["max_kappa_e"] < r5["need_kappa_e"]


def test_predicates_remain_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
