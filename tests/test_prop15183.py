"""Tests for Prop 15.183 — Max+⊥Max−, cross term; residual i still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15183 import (  # noqa: E402
    certify_cross_and_orth_p5,
    matching_nd2_psd_analytic,
    theorem_cross_gram,
    theorem_matching_nd2_psd,
    theorem_maxplus_orth_maxminus,
    theorem_ones_perp_ker,
    theorem_wedge_Gplus_sq,
)


def test_orth_proved():
    assert theorem_maxplus_orth_maxminus()["proved"] is True


def test_cross_gram_formula():
    r = theorem_cross_gram([p for p in range(3, 40) if is_prime(p)])
    assert r["proved"] is True
    for p in (3, 5, 7, 11):
        n = n_of(p)
        assert Fraction(n, 2 * p * p) == Fraction(n, 2 * p * p)


def test_ones_perp_ker():
    assert theorem_ones_perp_ker()["proved"] is True


def test_wedge_sq():
    r = theorem_wedge_Gplus_sq([3, 5, 7, 11])
    assert r["proved"] is True
    for p in (3, 5, 7):
        n = n_of(p)
        assert Fraction(2 * (n - 2), p * p) > 0


def test_matching_analytic():
    assert abs(matching_nd2_psd_analytic(0.0) + 2 * (2.0 ** 0.5)) < 1e-9
    assert abs(matching_nd2_psd_analytic(2.0) + 4.0) < 1e-9
    r = theorem_matching_nd2_psd()
    assert r["proved"] is True
    assert r["vs_need"]["7"]["wedge_kills"] is True
    assert r["vs_need"]["7"]["worst_matching_kills"] is False


def test_p5_census():
    r = certify_cross_and_orth_p5()
    assert bool(r["maxplus_orth_maxminus"]) is True
    assert bool(r["match"]) is True


def test_hinge_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
