"""Tests for Prop 15.192 — Gsum row algebra; Aut_e dual-eq; residual i open."""
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
from e1_gmin_m4_prop15192 import (  # noqa: E402
    certify_aut_e_kerbox,
    theorem_aut_e_averaging,
    theorem_avg_disj_gsum,
    theorem_gsum_diag_two,
    theorem_gsum_row_sum_n,
    theorem_three_halves_scheme_beats_need,
)


def test_gsum_row_algebra():
    assert theorem_gsum_diag_two()["proved"] is True
    assert theorem_gsum_row_sum_n()["proved"] is True
    r = theorem_avg_disj_gsum()
    assert r["proved"] is True
    assert r["by_p_sample"]["5"]["avg"] == str(Fraction(2, 23))


def test_aut_e_and_three_halves():
    assert theorem_aut_e_averaging()["proved"] is True
    r = theorem_three_halves_scheme_beats_need()
    assert r["proved"] is True
    assert r["by_p_sample"]["5"]["strict_lt"] is True
    assert r["by_p_sample"]["7"]["strict_lt"] is True


def test_aut_e_kerbox_census():
    r3 = certify_aut_e_kerbox(3)
    assert r3["blocked"] is False  # p=3 dual-eq Aut_e-feasible
    assert r3["max_kappa_e"] > r3["need"]
    r5 = certify_aut_e_kerbox(5)
    assert r5["blocked"] is True
    assert r5["max_kappa_e"] < r5["need"]
    assert r5["ratio_to_scheme"] < 1.5


def test_predicates_remain_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
