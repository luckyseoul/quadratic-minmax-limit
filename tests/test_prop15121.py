"""Tests for Prop 15.121 — spectral residual dictionary; Frobenius form."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, ew2_flat  # noqa: E402
from e1_gmin_m4_prop15121 import (  # noqa: E402
    certify_spectral_ed4,
    dim_Z,
    ed4_ub_from_16N_op,
    ed4_ub_from_H_op,
    ew2_flat_from_phi,
    fft_1perp_budget,
    is_prime,
    main,
    mu_bar,
    prove_ew2_flat_phi,
    prove_majorization_weak,
    prove_residual_frobenius,
)


def test_ew2_flat_phi_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_ew2_flat_phi(primes)
    assert r["proved"] is True
    assert ew2_flat_from_phi(5) == ew2_flat(5) == Fraction(8333, 5)
    # T^2/m = m mu_bar^2
    p = 5
    m = dim_Z(p)
    T = n_of(p) * (n_of(p) - 2)
    assert m * mu_bar(p) ** 2 == Fraction(T * T, m)


def test_frobenius_residual_budget():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_residual_frobenius(primes)
    assert r["proved"] is True
    # gap = room_hyp/4
    d = d_of(5)
    gap = fft_1perp_budget(5) - (ew2_flat(5) - d * d)
    assert gap == room_hyp_over_24(5) * 6


def test_H_16N_maj_too_weak():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_majorization_weak(primes)
    assert r["proved"] is True
    assert ed4_ub_from_H_op(5) > ed4_bud_closed(5)
    assert ed4_ub_from_16N_op(5) > ed4_bud_closed(5)


def test_cert_spectral_and_main_open():
    for p in (3, 5):
        c = certify_spectral_ed4(p)
        if c.get("skipped"):
            continue
        assert c["ED4_match"] is True
        assert c["EW2_match"] is True
        assert c["FFT_rayleigh_eq_dN"] is True
        assert c["FFT_1perp_match"] is True
        assert c["residual_le_hyp"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["ew2_flat_phi"] is True
    assert out["proved"]["frobenius_residual"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15121.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.121"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert "Frobenius" in data["settlement_note"] or "Frobenius" in str(
        data["closed_forms"]
    )
