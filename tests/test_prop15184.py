"""Tests for Prop 15.184 — T²κ identity; residual i still open."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15184 import (  # noqa: E402
    algebraic_T2_kappa_on_kappa1,
    certify_T2_identity_paley,
    certify_phi_bound_paley,
    theorem_T2_implies_abs_formula,
)


def test_T2_identity_algebra_proved():
    r = algebraic_T2_kappa_on_kappa1()
    assert r["proved"] is True
    assert r["n_kappa_abs_1_labelings"] == 48
    assert r["rest_equals_48_kappa"] is True


def test_abs_formula_statement():
    r = theorem_T2_implies_abs_formula()
    assert r["proved"] is True


def test_phi_bound_census_small():
    r = certify_phi_bound_paley([3, 5])
    assert r["certified_for_listed_primes"] is True
    assert r["proved_general"] is False
    assert r["by_p"]["5"]["max_abs_phi"] <= 6


def test_T2_paley_p5_census():
    r = certify_T2_identity_paley(5)
    assert bool(r["identity_holds"]) is True


def test_hinge_still_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
