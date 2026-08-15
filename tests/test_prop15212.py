"""Tests for Prop 15.212 — Veronese spanning criterion; hinge open."""
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
from e1_gmin_m4_prop15212 import (  # noqa: E402
    Gplus_pd_proved_general,
    certify_veronese_full_rank,
    dim_Wpp0,
    dim_constant_diag_plusplus,
    free_e_bound_proved_general,
    hinge_status_212,
    residual_i_closed_via_212,
    theorem_dim_comparison,
    theorem_gram_spectrum_identity,
    theorem_implication_rank_to_residual_i,
    theorem_veronese_criterion,
    veronese_full_rank_proved_general,
)


def test_veronese_criterion_proved():
    assert theorem_veronese_criterion()["proved"]


def test_gram_spectrum_proved():
    assert theorem_gram_spectrum_identity()["proved"]


def test_dim_comparison():
    t = theorem_dim_comparison()
    assert t["proved"]
    assert dim_constant_diag_plusplus(5) == dim_Wpp0(5) + 1
    assert dim_Wpp0(5) == 65
    assert dim_Wpp0(7) == 275


def test_implication_recorded():
    t = theorem_implication_rank_to_residual_i()
    assert t["proved"]
    assert "rank(𝒱)" in t["open_hypothesis"] or "rank" in t["open_hypothesis"]


def test_veronese_rank_cert():
    c = certify_veronese_full_rank()
    assert c["certified"]
    assert c["proved_general"] is False
    for p in ("3", "5", "7"):
        row = c["by_p"][p]
        assert row["full_rank"] is True
        assert row["rank_veronese"] == row["dim_Wpp0"]
    assert c["by_p"]["5"]["sharp_lambda_star"] is True
    assert c["by_p"]["5"]["min_rayleigh"] == str(Fraction(80, 13))


def test_predicates_honest_false():
    assert free_e_bound_proved_general() is False
    assert veronese_full_rank_proved_general() is False
    assert Gplus_pd_proved_general() is False
    assert residual_i_closed_via_212() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_212()
    assert h["residual_i_closed"] is False
    assert h["veronese_criterion_proved"] is True
    assert h["veronese_full_rank_general"] is False
