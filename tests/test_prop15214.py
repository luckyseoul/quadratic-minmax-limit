"""Tests for Prop 15.214 — master/δ structure; residual (i) OPEN."""
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
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15214 import (  # noqa: E402
    certify_delta_structure_small_p,
    hinge_status_214,
    mu_bound_proved_general,
    residual_i_closed_via_214,
    theorem_T_symmetric,
    theorem_kappa_perp_Epm4p,
    theorem_matching_psd_half,
    theorem_master_spectral_form,
    theorem_parseval_maxplus,
    theorem_particular_majorant_recall,
    theorem_resolvent_when_invertible,
)


def test_T_symmetric_proved():
    assert theorem_T_symmetric()["proved"]


def test_master_spectral_form_proved():
    assert theorem_master_spectral_form()["proved"]


def test_kappa_perp_proved():
    assert theorem_kappa_perp_Epm4p()["proved"]


def test_resolvent_proved():
    assert theorem_resolvent_when_invertible()["proved"]


def test_particular_majorant_recall():
    t = theorem_particular_majorant_recall()
    assert t["proved"]
    assert t["majorant_beats_thr"]
    for p in (5, 7, 11):
        maj = majorant_mu_part(p)
        thr = mu4_threshold(p)
        assert maj <= thr


def test_matching_psd_half():
    t = theorem_matching_psd_half()
    assert t["proved"]
    assert t["sufficient_for_residual_i"] is False


def test_parseval_proved():
    assert theorem_parseval_maxplus()["proved"]


def test_census_delta_structure():
    c = certify_delta_structure_small_p()
    assert c["certified"]
    assert not c["proved_general"]
    assert c["by_p"]["5"]["delta_nonzero_on_kappa1"]
    assert c["by_p"]["5"]["beats_thr"]


def test_predicates_false():
    assert mu_bound_proved_general() is False
    assert residual_i_closed_via_214() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_open():
    h = hinge_status_214()
    assert h["T_symmetric"]
    assert h["kappa_perp_Epm4p"]
    assert h["delta_control_general"] is False
    assert h["residual_i_closed"] is False
