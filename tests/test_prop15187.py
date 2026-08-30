"""Tests for Prop 15.187 — global T²κ; residual i still open."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15187 import (  # noqa: E402
    certify_T2_phi_formula,
    theorem_global_T2_kappa,
)


def test_global_T2_kappa_proved():
    r = theorem_global_T2_kappa()
    assert r["proved"] is True
    assert r["n_kappa1_labelings"] == 48
    assert r["n_kappa3_labelings"] == 16
    assert r["rest_equals_48_kappa_all"] is True


def test_T2_phi_census_p3_p5():
    r = certify_T2_phi_formula([3, 5])
    assert r["certified"] is True
    assert r["proved_general"] is False


def test_hinge_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
