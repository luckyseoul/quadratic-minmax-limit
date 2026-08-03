"""Tests for Prop 15.140 — weights, character-sum residual, geo coverage p=7."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15140 import (  # noqa: E402
    P7_DELTA2_FROM_Q,
    P7_G,
    P7_HPLUS,
    P7_N_FREE,
    P7_N_ORBITS,
    P7_ORBIT_TABLE,
    main,
    prove_character_sum_residual,
    prove_geo_coverage_live,
    prove_open,
    prove_weights,
)


def test_weights_stab_formula():
    a = prove_weights()
    assert a["proved"] is True
    assert a["G_order"] == P7_G == 2352
    assert a["H_plus"] == P7_HPLUS == 5726
    assert a["stab_formula_ok"] is True
    wsum = sum(Fraction(r["weight"]) for r in P7_ORBIT_TABLE)
    assert wsum == 1
    assert sum(r["size"] for r in P7_ORBIT_TABLE) == P7_HPLUS


def test_residual_le_room():
    b = prove_character_sum_residual()
    assert b["proved"] is True
    assert b["delta2_le_room"] is True
    assert b["matches_pair_avg"] is True
    assert P7_DELTA2_FROM_Q == P7_DELTA2_PAIR
    assert P7_DELTA2_PAIR <= room_hyp_over_24(7)


def test_geo_seven_of_eight():
    c = prove_geo_coverage_live()
    assert c["proved"] is True
    assert c["n_free"] == P7_N_FREE == 7
    assert c["n_orbits"] == P7_N_ORBITS == 8
    assert c["found_sizes"][56] == 56
    assert c["found_sizes"][84] == 84
    assert c["found_sizes"][294] == 294
    assert c["found_sizes"][588] == 588
    assert c["found_sizes"][1176] == 1176


def test_open_honest():
    d = prove_open()
    assert d["residual_closed_general"] is False
    assert d["p7_all_types_maxplus_free"] is False
    assert d["residual_form_p7_certified"] is True


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["character_sum_residual_p7"] is True
    assert out["proved"]["p7_all_eight_maxplus_free"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15140.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.140"
