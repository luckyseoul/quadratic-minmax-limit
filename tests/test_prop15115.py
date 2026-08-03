"""Tests for Prop 15.115 — Aρ=b; δ=P_E4p m4; spectral moments of f_y."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15115 import (  # noqa: E402
    certify_ed4_residual,
    certify_resolvent_and_moments,
    is_prime,
    main,
    prove_4p_not_pm_mu,
    prove_spectral_moments,
    spectral_m1,
    spectral_m2,
    spectral_var,
)


def test_spectral_moments_algebra():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_spectral_moments(primes)
    assert r["proved"] is True
    assert spectral_m1(5) == Fraction(4 * 5) - Fraction(12, 5)
    assert spectral_var(5) == Fraction(24 * (25 - 3) * (25 - 4), 25 * (25 - 2))
    # var form positive for p≥3
    for p in (3, 5, 7, 11, 13):
        assert spectral_var(p) > 0
        assert spectral_m2(p) > spectral_m1(p) ** 2


def test_4p_not_pm_mu():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_4p_not_pm_mu(primes)
    assert r["proved"] is True


def test_cert_resolvent_p3_p5():
    for p in (3, 5):
        c = certify_resolvent_and_moments(p)
        if c.get("skipped"):
            continue
        assert c["A_rho_eq_b"] is True
        assert c["E_gamma_f_eq_2k_over_p"] is True
        assert c["A_m4_eq_4k_over_p"] is True
        assert c["m1_ok"] is True
        assert c["m2_ok"] is True
        assert c["Q_delta_constant"] is True
        assert c["Q_eq_delta2"] is True
        assert c["delta2_le_rho_min2"] is True
        assert c["kappa_perp_delta"] is True


def test_ed4_cert_and_main_open():
    for p in (3, 5, 7):
        e = certify_ed4_residual(p)
        if e.get("skipped"):
            continue
        assert e["ED4_le_suf"] is True
        assert e["delta2_le_rho_min2"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["A_rho_eq_b"] is True
    assert out["proved"]["spectral_moments"] is True
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15115.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.115"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
