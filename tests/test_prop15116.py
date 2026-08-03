"""Tests for Prop 15.116 — e4↔ED4↔δ dictionary; coherent mass."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15110 import e4_newton  # noqa: E402
from e1_gmin_m4_prop15116 import (  # noqa: E402
    certify_dictionary,
    dmax_of,
    e4_poly,
    is_prime,
    kappa_dot_rho_min,
    main,
    prove_e4_poly,
    prove_flat_identity,
    prove_kappa_rho_min,
    prove_min_dist_envelope,
    sum_m4_sq_from_ed4,
)


def test_e4_poly():
    primes = [p for p in range(3, 40) if is_prime(p)]
    r = prove_e4_poly(primes)
    assert r["proved"] is True
    for n, s in [(10, 4), (26, 6), (50, 0), (50, 34)]:
        assert e4_poly(s, n) == e4_newton(s, n)


def test_kappa_rho_min_and_flat():
    primes = [p for p in range(3, 60) if is_prime(p)]
    k = prove_kappa_rho_min(primes)
    assert k["proved"] is True
    assert kappa_dot_rho_min(5) == Fraction(312)
    assert kappa_dot_rho_min(7) == Fraction(1200)
    f = prove_flat_identity(primes)
    assert f["proved"] is True


def test_min_dist_envelope_too_weak():
    r = prove_min_dist_envelope([3, 5, 7, 11, 13])
    assert r["proved_form"] is True
    assert r["closes_residual"] is False
    assert dmax_of(5) == 14
    assert dmax_of(7) == 34
    # p>=5: N=n envelope fails suf
    assert r["by_p_N_eq_n"]["5"]["ub_le_suf"] is False
    assert r["by_p_known_N"]["5"]["ub_le_suf"] is False


def test_cert_and_main_open():
    for p in (3, 5, 7):
        c = certify_dictionary(p)
        if c.get("skipped"):
            continue
        assert c["delta2_match"] is True
        assert c["delta2_le_rho_min2"] is True
        assert c["ED4_le_suf"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["e4_poly"] is True
    assert out["proved"]["kappa_rho_min"] is True
    assert out["proved"]["flat_identity"] is True
    assert out["proved"]["min_dist_closes_residual"] is False
    assert out["proved"]["ed4_le_suf_for_all_p"] is False
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15116.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.116"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
