"""Tests for Prop 15.111 — pair Schur; closed α_κ, α_ρ; Φ residual = 8⟨δ,κ_B⟩."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15111 import (  # noqa: E402
    alpha_kappa,
    alpha_rho,
    beta_Tb,
    beta_b,
    certify_alpha_kappa_formula_on_Z,
    certify_schur_scalars,
    is_prime,
    main,
    mubar,
    pair_const,
    prove_alpha_kappa_on_Z,
    prove_channel_alpha_rho,
    prove_pair_and_alpha_rho,
    prove_zero_diag_pairing,
    sixteen_N_delta_budget,
)


def test_zero_diag_pairing_identity():
    r = prove_zero_diag_pairing([6, 7, 8], trials=15)
    assert r["proved"] is True
    assert r["max_abs_err"] < 1e-8


def test_alpha_kappa_algebra():
    primes = [p for p in range(3, 50) if is_prime(p)]
    r = prove_alpha_kappa_on_Z(primes)
    assert r["proved"] is True
    assert alpha_kappa(5) == Fraction(27, 100)
    assert alpha_kappa(3) == Fraction(11, 36)
    assert alpha_kappa(7) == Fraction(51, 196)


def test_pair_and_alpha_rho_closed():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_pair_and_alpha_rho(primes)
    assert r["proved"] is True
    # spot checks
    assert pair_const(5) == Fraction(9, 20)
    assert alpha_rho(5) == Fraction(9, 50)
    assert alpha_rho(3) == Fraction(17, 18)
    assert alpha_rho(7) == Fraction(87, 1078)
    assert pair_const(5) == alpha_kappa(5) + alpha_rho(5)
    assert pair_const(5) == (mubar(5) - 6) / 8


def test_channel_reconstructs_alpha_rho():
    primes = [p for p in range(3, 80) if is_prime(p)]
    r = prove_channel_alpha_rho(primes)
    assert r["proved"] is True
    assert beta_b(5) == Fraction(6, 5)
    assert beta_Tb(5) == Fraction(96, 5)


def test_sixteen_N_delta_budget():
    # (n-10)/(n-6); n=26 at p=5 → 16/20 = 4/5
    assert sixteen_N_delta_budget(5) == Fraction(4, 5)
    n7 = 50
    assert sixteen_N_delta_budget(7) == Fraction(n7 - 10, n7 - 6)


def test_schur_cert_p3_p5():
    for p in (3, 5):
        c = certify_schur_scalars(p, ntrials=5)
        assert c.get("skipped") is not True
        assert c["all_scalar"] is True
        assert c["all_match_closed"] is True


def test_alpha_kappa_collapse_on_Z():
    for p in (3, 5):
        c = certify_alpha_kappa_formula_on_Z(p, ntrials=4)
        assert c["holds"]


def test_main_honest_open():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["zero_diag_pairing"] is True
    assert out["proved"]["alpha_kappa_on_Z"] is True
    assert out["proved"]["pair_and_alpha_rho_forms"] is True
    assert out["proved"]["channel_alpha_rho"] is True
    assert out["proved"]["schur_certified_p3_5_7"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["sum_rho2_le_T_for_all_p"] is False
    assert out["proved"]["sixteen_N_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15111.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.111"
    assert data["L_status"] == "OPEN"
    assert "soft-close" not in data["settlement_note"].lower() or "No soft-close" in data["settlement_note"]
