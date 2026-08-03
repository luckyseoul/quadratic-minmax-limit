"""Tests for Prop 15.118 — pointwise bgf=0; T²κ; ρ_min·f closed."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15118 import (  # noqa: E402
    T2kappa_dot_m4_formula,
    certify_pointwise,
    is_prime,
    kappa_gamma_f_formula,
    main,
    prove_T2kappa_m4,
    prove_kappa_gamma_f,
    prove_rho_min_f_from_bgf0,
    rho_min_dot_f_if_bgf_zero,
)


def test_kappa_gamma_f_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_kappa_gamma_f(primes)
    assert r["proved"] is True
    assert kappa_gamma_f_formula(3) == Fraction(300)
    assert kappa_gamma_f_formula(5) == Fraction(16380)


def test_T2kappa_m4_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_T2kappa_m4(primes)
    assert r["proved"] is True
    assert T2kappa_dot_m4_formula(3) == Fraction(5760)
    assert T2kappa_dot_m4_formula(5) == Fraction(124800)


def test_rho_min_f_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_rho_min_f_from_bgf0(primes)
    assert r["proved"] is True
    assert rho_min_dot_f_if_bgf_zero(5) == Fraction(208, 5)
    assert rho_min_dot_f_if_bgf_zero(3) == Fraction(80, 3)


def test_cert_pointwise_and_main_open():
    for p in (3, 5):
        c = certify_pointwise(p)
        if c.get("skipped"):
            continue
        assert c["bgf_pointwise_zero"] is True
        assert c["T2kappa_f_constant_target"] is True
        assert c["rho_min_f_closed"] is True
        assert c["kappa_gf_closed"] is True
        assert c["b_f_closed"] is True
        assert c["delta2_le_hyp"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["pointwise_bgf_zero"] is True
    assert out["proved"]["kappa_gamma_f"] is True
    assert out["proved"]["T2kappa_m4"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15118.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.118"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
