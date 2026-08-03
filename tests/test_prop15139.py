"""Tests for Prop 15.139 — affine halfspaces + double switch; p=7 size classes."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import (  # noqa: E402
    affine_halfspace,
    is_4_AP,
    main,
    prove_affine_halfspaces,
    prove_ap_dichotomy_p7,
    prove_double_switch,
    prove_open,
    prove_p7_all_size_classes,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_affine_halfspaces_evec():
    a = prove_affine_halfspaces([5, 7])
    assert a["proved"] is True
    assert a["by_p"]["5"]["all_evec"] is True
    assert a["by_p"]["7"]["all_evec"] is True
    C = paley_conference_prime_power(7).astype(float)
    y = affine_halfspace(7, (0, 1), {0, 1, 2, 4})
    assert np.allclose(C @ y, 7 * y, atol=1e-5)


def test_ap_dichotomy():
    assert is_4_AP({0, 1, 2, 3}, 7) is True
    assert is_4_AP({0, 1, 2, 4}, 7) is False
    b = prove_ap_dichotomy_p7()
    assert b["proved"] is True
    assert b["n_AP_size84"] == 21
    assert b["n_nonAP_size56"] == 14
    assert b["QR_half_is_nonAP"] is True


def test_double_switch_and_p7_sizes():
    assert prove_double_switch()["proved"] is True
    d = prove_p7_all_size_classes()
    assert d["proved"] is True
    assert d["orbit_sizes_verified"][56] == 56
    assert d["orbit_sizes_verified"][84] == 84
    assert d["orbit_sizes_verified"][294] == 294
    assert d["orbit_sizes_verified"][588] == 588
    assert d["orbit_sizes_verified"][1176] == 1176


def test_open_honest():
    e = prove_open()
    assert e["residual_closed_general"] is False
    assert e["p7_size_classes_maxplus_free"] is True
    assert e["p7_residual_maxplus_free"] is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["p7_all_size_classes_maxplus_free"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15139.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.139"
