"""Tests for Prop 15.145 — type-free residual package."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15112 import ed4_suf  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_flat  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15145 import (  # noqa: E402
    P5_DELTA2,
    P7_DELTA2,
    ed4_suf_closed,
    m4f_suf,
    main,
    prove_census,
    prove_dictionary,
    prove_gap_asymptotic,
    prove_open,
    prove_typefree_program,
    primes_ge,
)


def test_dictionary_closed_forms():
    a = prove_dictionary(primes_ge(5, 31))
    assert a["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert m4f_suf(p) == m4f_flat(p) + rho_min2(p)
        assert ed4_suf_closed(p) == ed4_suf(p)
        if p >= 7:
            assert rho_min2(p) < room_hyp_over_24(p)
        if p == 5:
            assert rho_min2(p) > room_hyp_over_24(p)


def test_gap_asymptotic():
    b = prove_gap_asymptotic(primes_ge(5, 97))
    assert b["proved"] is True
    assert b["asymptotic_ratio"] == "5/8"
    assert b["monotone_decreasing_to_5_8"] is True
    assert b["final_ratio"] < 0.63


def test_typefree_program_and_census():
    c = prove_typefree_program()
    assert c["proved_structure"] is True
    assert "sum_m4_sq_le_m4f_suf" in c["targets"]
    assert "type_list_residual_p_ge_11" in c["dead"]
    d = prove_census()
    assert d["proved"] is True
    assert d["p5_saturates_room"] is True
    assert P5_DELTA2 == room_hyp_over_24(5)
    assert P7_DELTA2 == P7_DELTA2_PAIR
    assert P5_DELTA2 <= rho_min2(5)
    assert P7_DELTA2 <= rho_min2(7)
    # equivalence ||ρ||² ≤ 2 ρ_min²
    for p, d2 in ((5, P5_DELTA2), (7, P7_DELTA2)):
        assert rho_min2(p) + d2 <= 2 * rho_min2(p)
        assert m4f_flat(p) + d2 <= m4f_suf(p)


def test_open_and_main():
    e = prove_open()
    assert e["residual_closed_general"] is False
    assert e["type_free_delta2_le_rho_min2_all_p"] is False
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["typefree_dictionary"] is True
    assert out["proved"]["gap_asymptotic_5_8"] is True
    assert out["proved"]["census_delta2_le_rho_min2_p5_p7"] is True
    assert out["proved"]["type_free_delta2_le_rho_min2_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15145.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.145"
    # sample closed forms parse as rationals
    for p in ("5", "7", "11"):
        Fraction(data["closed_forms"]["sample"][p]["m4f_suf"])
