"""Tests for Prop 15.137 — c_j over Max+ G-orbits; p=5 two-type formula."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402
from e1_gmin_m4_prop15137 import (  # noqa: E402
    MAX_G_CENSUS,
    main,
    prove_G_equivariance,
    prove_Qhs_maxplus_free,
    prove_census,
    prove_hemisphere,
    prove_p5_two_type_formula,
)


def test_equivariance_and_hemisphere():
    assert prove_G_equivariance()["proved"] is True
    assert prove_hemisphere()["proved"] is True


def test_census_mass_and_weights():
    c = prove_census()
    assert c["proved"] is True
    # p=5 weights
    w = MAX_G_CENSUS[5]["weights_hemisphere"]
    assert w[0] == Fraction(3, 13)
    assert w[1] == Fraction(10, 13)
    assert sum(w) == 1
    assert 2 * (30 + 100) == 260
    # p=7 mass
    h = 56 + 84 + 294 + 588 + 4 * 1176
    assert 2 * h == 11452
    assert MAX_G_CENSUS[7]["N"] == 11452
    # p=3 single type
    assert MAX_G_CENSUS[3]["r_hemisphere_types"] == 1
    assert G_SPECTRUM[3]["nullity"] == 0


def test_p5_formula_and_Qhs():
    d = prove_p5_two_type_formula()
    assert d["proved"] is True
    assert d["match_room"] is True
    assert room_hyp_over_24(5) == Fraction(1536, 65)
    assert d["Q_hs_maxplus_free"] is True
    assert d["y_star_maxplus_free"] is False
    e = prove_Qhs_maxplus_free()
    assert e["proved"] is True
    assert e["sufficient_for_residual_alone"] is False


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["p5_two_type_formula"] is True
    assert out["proved"]["Qhs_maxplus_free"] is True
    assert out["proved"]["gauss_cj_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["status"]["hs_part_maxplus_free"] is True
    assert out["status"]["full_cj_general_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15137.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.137"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["closed_forms"]["p5"] == "c_j = (3/13) Q_j(hs) + (10/13) Q_j(y_*)"
