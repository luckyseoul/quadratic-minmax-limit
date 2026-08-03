"""Tests for Prop 15.107 — 16N from mult≥d−1 + room_hyp."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of, wick_hi  # noqa: E402
from e1_gmin_m4_prop15105 import room_hyp_of  # noqa: E402
from e1_gmin_m4_prop15107 import (  # noqa: E402
    L_majorization,
    certify_hyp_and_16n,
    is_prime,
    main,
    prove_four_design_le_wick,
    prove_gegenbauer_coeffs,
    prove_general_status,
    prove_theorem_A,
)


def test_theorem_A_all_primes_to_50():
    primes = [p for p in range(3, 50) if is_prime(p)]
    r = prove_theorem_A(primes)
    assert r["proved"] is True
    for p in primes:
        assert r["by_p"][str(p)]["L_le_16"] is True


def test_theorem_A_bound_below_16_strict_for_p_gt_3():
    for p in [5, 7, 11, 13, 17, 19]:
        L = L_majorization(p, room_hyp_of(p), d_of(p) - 1)
        assert L < 16.0


def test_gegenbauer_closed_forms():
    r = prove_gegenbauer_coeffs([3, 5, 7, 11, 13, 17, 19, 23])
    assert r["proved"] is True
    # spot-check α, β
    d = 13
    assert r["by_p"]["5"]["alpha"] == str(Fraction(-6, d + 4))
    assert r["by_p"]["5"]["beta"] == str(Fraction(3, (d + 2) * (d + 4)))


def test_four_design_le_wick():
    primes = [p for p in range(3, 40) if is_prime(p)]
    r = prove_four_design_le_wick(primes)
    assert r["proved"] is True
    for p in primes:
        d = d_of(p)
        e4 = Fraction(48 * d**3, d + 2)
        assert e4 <= wick_hi(p)


def test_general_status_honest():
    g = prove_general_status()
    assert g["theorem_A_proved"] is True
    assert g["sixteen_N_proved_for_all_p"] is False
    assert g["orth_le_room_hyp_for_all_p"] is False


def test_cert_p3_p5():
    for p in (3, 5):
        c = certify_hyp_and_16n(p)
        if c.get("skipped"):
            continue
        assert c["orth_le_room_hyp"] is True
        assert c["lambda_max_le_16"] is True
        assert c["thmA_le_16"] is True


def test_main_writes_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["theorem_A_16N_from_hyp"] is True
    assert out["proved"]["sixteen_N_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15107.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.107"
    assert data["theorem_A"]["proved"] is True
