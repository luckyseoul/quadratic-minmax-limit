"""Tests for Prop 15.181 — vertex-star, κ-counts; residual i still open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    dual_equality_correlation_need,
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15181 import (  # noqa: E402
    certify_ksparse_box_p5,
    gsum_star_sum_identity,
    n1_through_edge,
    n3_through_edge,
    sum_kappa_sq_through_edge,
    theorem_kappa_counts_through_edge,
    theorem_matching_nd2_psd_lb,
    vertex_star_feature_identity,
)


def test_vertex_star_identities_proved():
    assert vertex_star_feature_identity()["proved"] is True
    assert gsum_star_sum_identity()["proved"] is True


def test_kappa_counts_closed_forms():
    for p in [3, 5, 7, 11, 13, 17, 19, 23]:
        if not is_prime(p):
            continue
        n = n_of(p)
        D = Fraction((n - 2) * (n - 3), 2)
        n1 = n1_through_edge(p)
        n3 = n3_through_edge(p)
        assert n1 + n3 == D
        assert n1 + 9 * n3 == sum_kappa_sq_through_edge(p)
        assert n3 == Fraction((n - 2) * (n - 6), 8)
        assert n1 == Fraction(3 * (n - 2) * (n - 2), 8)


def test_kappa_counts_theorem_all_primes():
    r = theorem_kappa_counts_through_edge(
        [p for p in range(3, 60) if is_prime(p)]
    )
    assert r["proved"] is True


def test_matching_nd2_psd_floor():
    r = theorem_matching_nd2_psd_lb()
    assert r["proved"] is True
    assert r["min_g1_plus_g2"] == -4.0
    assert r["matching_weaker_than_wedge"] is True
    # does not kill p=5 dual eq by PSD alone
    need5 = float(dual_equality_correlation_need(5))
    # p=5: matching floor −4 and wedge −2√2 both sit below need −2.8 → no PSD kill
    assert -4.0 < need5
    assert need5 > r["wedge_lb"]


def test_p5_ksparse_box_census():
    r = certify_ksparse_box_p5()
    assert bool(r["dual_eq_linear_empty_p5"]) is True
    assert bool(r["need_lt_lp_min"]) is True
    assert bool(r["certified_minus_6_over_65"]) is True
    assert r["not_general_proof"] is True


def test_hinge_still_open_no_soft_close():
    """Predicates must stay false — 15.181 does not close residual i."""
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
