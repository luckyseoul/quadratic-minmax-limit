"""Tests for Prop 15.120 — E[W²] factorization; Pythagoras; majorization."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, m4f_flat  # noqa: E402
from e1_gmin_m4_prop15120 import (  # noqa: E402
    F_min_mass,
    certify_factorization,
    cs_gamma_ub,
    crude_lp_ub,
    ed4_majorization_ub,
    is_prime,
    main,
    prove_dead_ubs,
    prove_majorization,
    prove_pythagoras,
)


def test_pythagoras_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_pythagoras(primes)
    assert r["proved"] is True
    assert F_min_mass(5) == m4f_flat(5) == Fraction(598, 5)
    assert F_min_mass(3) == Fraction(110, 3)


def test_majorization_too_weak():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_majorization(primes)
    assert r["proved"] is True
    for p in (5, 7, 11, 13):
        assert ed4_majorization_ub(p) > ed4_bud_closed(p)
    # 2n³ at p=5
    assert ed4_majorization_ub(5) == Fraction(2 * 26**3)


def test_dead_ubs_exceed_budget():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_dead_ubs(primes)
    assert r["proved"] is True
    assert cs_gamma_ub(5) > room_hyp_over_24(5)
    assert crude_lp_ub(5) > ed4_bud_closed(5)


def test_cert_factorization_and_main_open():
    for p in (3, 5):
        c = certify_factorization(p)
        if c.get("skipped"):
            continue
        assert c["mean_matches_factorization"] is True
        assert c["EW2_constant_on_Max"] is True
        assert c["residual_le_hyp"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["factorization"] is True
    assert out["proved"]["pythagoras"] is True
    assert out["proved"]["majorization_ub"] is True
    assert out["proved"]["majorization_closes_residual"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15120.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.120"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
