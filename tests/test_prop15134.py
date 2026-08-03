"""Tests for Prop 15.134 — strict Aut(C) Bose–Mesner; residual projection."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_flat  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import (  # noqa: E402
    G_SPECTRUM,
    group_order,
    main,
    prove_orbits_spectrum,
    prove_residual_projection,
    prove_strict_aut_group,
)


def test_group_order_and_spectrum():
    assert group_order(5) == 600
    assert group_order(7) == 2352
    assert group_order(3) == 72
    a = prove_strict_aut_group([3, 5, 7, 11])
    assert a["proved"] is True
    b = prove_orbits_spectrum()
    assert b["proved"] is True
    assert G_SPECTRUM[3]["nullity"] == 0
    assert G_SPECTRUM[5]["nullity"] == 2
    assert G_SPECTRUM[7]["nullity"] == 7
    assert G_SPECTRUM[5]["n_orbits"] == 42
    assert G_SPECTRUM[7]["n_orbits"] == 128
    assert G_SPECTRUM[5]["n_orbits"] > G_SPECTRUM[5]["class_key_orbits"]
    assert G_SPECTRUM[7]["n_orbits"] > G_SPECTRUM[7]["class_key_orbits"]
    # Aut-line dim≤1 fails for G at p≥5
    assert G_SPECTRUM[5]["nullity"] > 1
    assert G_SPECTRUM[7]["nullity"] > 1


def test_residual_projection_exact():
    c = prove_residual_projection()
    assert c["proved"] is True
    assert Fraction(c["p5"]["delta2"]) == room_hyp_over_24(5)
    assert Fraction(c["p5"]["delta2"]) == Fraction(1536, 65)
    assert Fraction(c["p7"]["delta2"]) == P7_DELTA2_PAIR
    assert Fraction(c["p7"]["delta2"]) <= room_hyp_over_24(7)
    # sum m4^2 = flat + delta2 at p=5
    assert Fraction(c["p5"]["sum_m4_sq"]) == m4f_flat(5) + room_hyp_over_24(5)


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["strict_aut_group"] is True
    assert out["proved"]["orbits_spectrum_p3_5_7"] is True
    assert out["proved"]["residual_projection_p5_p7"] is True
    assert out["proved"]["aut_line_dim_le_1_for_G"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["status"]["G_carries_residual_p7"] is True
    assert out["status"]["class_key_carries_residual_p7"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15134.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.134"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["closed_forms"]["nullity_E4p_G"]["5"] == 2
    assert data["closed_forms"]["nullity_E4p_G"]["7"] == 7
