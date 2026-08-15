"""Tests for Prop 15.185 — |φ| path; residual i still open."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15185 import (  # noqa: E402
    certify_kappa_cr_dichotomy,
    certify_phi_bound_paley,
    certify_supersingular_only_11,
    theorem_H_mod_4,
    theorem_edge_only_pm_2p,
    theorem_hasse,
)


def test_H_mod_4_proved():
    assert theorem_H_mod_4()["proved"] is True


def test_hasse_proved():
    assert theorem_hasse()["proved"] is True


def test_edge_ladder_proved():
    assert theorem_edge_only_pm_2p()["proved"] is True


def test_supersingular_class_small_primes():
    r = certify_supersingular_only_11([3, 5, 7])
    assert r["all_ok"] is True
    assert r["proved_general_algebra"] is False
    assert r["by_p"]["5"]["max_abs_H_outside_11"] <= 6


def test_kappa_cr_dichotomy_p5():
    r = certify_kappa_cr_dichotomy(5, 800)
    assert r["dichotomy_holds"] is True


def test_phi_bound_census():
    r = certify_phi_bound_paley([3, 5])
    assert r["certified"] is True
    assert r["proved_general"] is False


def test_hinge_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
