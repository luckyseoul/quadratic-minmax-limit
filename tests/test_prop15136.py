"""Tests for Prop 15.136 — Max+-free flat on G; m4=flat+δ; residual OPEN."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_flat  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM, group_order  # noqa: E402
from e1_gmin_m4_prop15136 import (  # noqa: E402
    FLAT_CERT,
    P5_INF_M4,
    main,
    prove_flat_maxplus_free,
    prove_geom_insufficient,
    prove_m4_decomposition,
    prove_partial_inf_m4_p5,
)


def test_flat_maxplus_free():
    a = prove_flat_maxplus_free()
    assert a["proved"] is True
    for p in (3, 5, 7, 11, 13):
        assert "rho_min2" in a["by_p"][str(p)] or "construction_available" in a["by_p"][str(p)]
    # closed forms exist for all primes
    for p in (5, 7, 11, 13, 17, 19):
        assert rho_min2(p) > 0
        assert m4f_flat(p) > 0
    assert group_order(11) == 11**2 * (11**2 - 1)


def test_m4_decomposition_and_geom():
    b = prove_m4_decomposition()
    assert b["proved"] is True
    assert FLAT_CERT[5]["delta2"] == room_hyp_over_24(5) == Fraction(1536, 65)
    assert FLAT_CERT[7]["delta2"] == P7_DELTA2_PAIR
    assert FLAT_CERT[5]["nu"] == G_SPECTRUM[5]["nullity"] == 2
    assert FLAT_CERT[7]["nu"] == G_SPECTRUM[7]["nullity"] == 7
    c = prove_geom_insufficient()
    assert c["proved"] is True
    assert c["p5"]["types_with_unique_m4"] < c["p5"]["n_geom_types"]


def test_partial_p5_and_main_open():
    e = prove_partial_inf_m4_p5()
    assert e["proved_census"] is True
    assert e["proved_general"] is False
    assert P5_INF_M4["inf_k-3_fp1"] == Fraction(-1, 5)
    assert P5_INF_M4["inf_k-3_fp2"] == Fraction(-21, 65)

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["flat_maxplus_free_on_G"] is True
    assert out["proved"]["m4_equals_flat_plus_delta"] is True
    assert out["proved"]["gauss_cj_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["status"]["flat_for_general_p"] is True
    assert out["status"]["free_coefficients_remain"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15136.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.136"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert "Max+-free" in data["settlement_note"] or "flat" in data["settlement_note"]
