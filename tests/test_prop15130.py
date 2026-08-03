"""Tests for Prop 15.130 — P m4=δ; gap algebra; census δ²≤ρ_min²."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15130 import (  # noqa: E402
    delta2_from_W,
    gap_room_minus_rho_min,
    is_prime,
    main,
    prove_P_m4_is_delta,
    prove_aut_line_program,
    prove_census_delta_le_rho_min,
    prove_gap_algebra,
    prove_jensen_refined,
)


def test_P_m4_and_jensen():
    assert prove_P_m4_is_delta()["proved"] is True
    assert prove_jensen_refined()["proved"] is True
    assert prove_aut_line_program()["proved_structure"] is True


def test_gap_algebra():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_gap_algebra(primes)
    assert r["proved"] is True
    assert r["all_primes_7_to_97_gap_positive"] is True
    # explicit checks
    assert gap_room_minus_rho_min(7) == room_hyp_over_24(7) - rho_min2(7)
    assert gap_room_minus_rho_min(7) > 0
    assert gap_room_minus_rho_min(5) < 0  # rho_min > room at p=5
    # closed form
    p = 11
    num = (p - 1) * (p + 1) * (3 * p**6 - 105 * p**4 + 37 * p * p - 15)
    den = 6 * p * p * (p * p - 5) * (p * p + 1)
    assert gap_room_minus_rho_min(p) == Fraction(num, den)


def test_census_delta_le_rho_min():
    r = prove_census_delta_le_rho_min()
    assert r["proved"] is True
    for p in (3, 5, 7):
        d2 = delta2_from_W(p)
        assert d2 <= rho_min2(p)
        assert d2 <= room_hyp_over_24(p)
    # p=7 strict ratios
    d2 = delta2_from_W(7)
    assert d2 == Fraction(82176, 4499)
    assert d2 < rho_min2(7)
    assert d2 < room_hyp_over_24(7)


def test_main_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["P_m4_equals_delta"] is True
    assert out["proved"]["gap_rho_min_lt_room_p_ge_7"] is True
    assert out["proved"]["census_delta_le_rho_min"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15130.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.130"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert "ρ_min" in data["settlement_note"] or "rho_min" in data["open_attack"]
