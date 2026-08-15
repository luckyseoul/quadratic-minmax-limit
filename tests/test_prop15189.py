"""Tests for Prop 15.189 — π=C/p; G+ PSD bound; residual i still open."""
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
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15189 import (  # noqa: E402
    certify_CF_not_pF,
    certify_Gplus_eigenvalue_formula,
    certify_pairwise_pi,
    certify_sum_identity,
    gplus_psd_max_m,
    theorem_Gplus_psd_bound,
    theorem_adjacent_gsum_zero,
    theorem_pairwise_pi_equals_C_over_p,
    theorem_sum_identity,
)


def test_sum_identity_proved():
    assert theorem_sum_identity()["proved"] is True
    for p in (3, 5):
        r = certify_sum_identity(p)
        assert r["maxplus_sum_identity"] is True
        assert r["maxminus_sum_identity"] is True


def test_pairwise_pi_proved_and_certified():
    assert theorem_pairwise_pi_equals_C_over_p()["proved"] is True
    for p in (3, 5):
        r = certify_pairwise_pi(p)
        assert r["ok"] is True
        assert r["max_err_pi_plus"] < 1e-9


def test_CF_not_pF():
    for p in (3, 5, 7):
        assert certify_CF_not_pF(p)["nonzero"] is True


def test_adjacent_gsum_from_pi():
    assert theorem_adjacent_gsum_zero()["proved"] is True


def test_Gplus_psd_bound():
    assert theorem_Gplus_psd_bound()["proved"] is True
    assert theorem_Gplus_psd_bound()["sufficient_for_residual_i"] is False
    for p in (3, 5, 7):
        assert gplus_psd_max_m(p) == Fraction(p - 2, p)
        r = certify_Gplus_eigenvalue_formula(p)
        assert r["psd_up_to_bound"] is True
        # bound does not beat 1/(2p) for p>=5
        if p >= 5:
            assert gplus_psd_max_m(p) > mu4_threshold(p)


def test_predicates_remain_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
