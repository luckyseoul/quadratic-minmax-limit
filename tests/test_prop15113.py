"""Tests for Prop 15.113 — ⟨f_y,Tκ⟩; ED4 via W; Q_δ criterion."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15112 import ed4_suf  # noqa: E402
from e1_gmin_m4_prop15113 import (  # noqa: E402
    certify_identities,
    fy_Tkappa_formula,
    is_prime,
    main,
    prove_fy_Tkappa_from_b,
)


def test_fy_Tkappa_algebra():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_fy_Tkappa_from_b(primes)
    assert r["proved"] is True
    assert fy_Tkappa_formula(3) == Fraction(480)
    assert fy_Tkappa_formula(5) == Fraction(6240)
    assert fy_Tkappa_formula(7) == Fraction(2 * 7 * (7**4 - 1))


def test_cert_p3_p5():
    for p in (3, 5):
        c = certify_identities(p)
        if c.get("skipped"):
            continue
        assert c["kappa_Tkappa_zero"] is True
        assert c["fy_Tkappa_constant"] is True
        assert c["fy_Tkappa_matches"] is True
        assert c["E_W_ok"] is True
        assert c["ED4_W_identity"] is True
        assert c["ED4_le_suf"] is True
        assert c["delta2_le_rho_min2"] is True


def test_cert_p7_if_available():
    c = certify_identities(7)
    if c.get("skipped"):
        return
    assert c["kappa_Tkappa_zero"] is True
    assert c["fy_Tkappa_matches"] is True
    assert c["E_W_ok"] is True
    assert c["ED4_le_suf"] is True


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["fy_Tkappa_formula"] is True
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15113.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.113"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
