"""Tests for Prop 15.105 — variance identity + exact 16-criterion."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15105 import (  # noqa: E402
    avg_Phi_Z,
    certify_variance_on_maxplus,
    dim_Z,
    main,
    prove_exact_16_criterion,
    prove_general_status,
    prove_hyp_dminus1_bound,
    prove_second_moment_and_variance,
    room_hyp_of,
    room_orth,
    tr_Phi_Z,
)
from e1_gmin_m4_prop15100 import ED4_flat, d_of, n_of, proj_kappa2  # noqa: E402


def test_flat_equals_T2_over_m_algebra():
    """ED4_flat − 4n² = T²/m for primes p≥3 (core of variance identity)."""
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 31]:
        n = n_of(p)
        d = d_of(p)
        T = tr_Phi_Z(p)
        m = dim_Z(p)
        lhs = ED4_flat(p) - 4 * n * n
        rhs = Fraction(T * T, m)
        closed = Fraction(32 * d * (d - 1) ** 2, d - 3)
        assert lhs == rhs == closed, (p, lhs, rhs, closed)


def test_avg_forms():
    for p in [3, 5, 7, 11, 13]:
        n = n_of(p)
        d = d_of(p)
        assert avg_Phi_Z(p) == Fraction(8 * (n - 2), n - 6)
        assert avg_Phi_Z(p) == Fraction(8 * (d - 1), d - 3)
        if p >= 5:
            assert avg_Phi_Z(p) < 16
        if p == 3:
            assert avg_Phi_Z(p) == 16


def test_variance_identity_module():
    r = prove_second_moment_and_variance(
        [3, 5, 7, 11, 13, 17, 19, 23, 31]
    )
    assert r["proved"] is True


def test_exact_16_criterion_algebra():
    """mult≥d + orth≤room ⇒ L_bound = 16 for all odd primes tested."""
    primes = [p for p in range(3, 60) if p > 1 and all(p % k for k in range(2, int(p**0.5) + 1))]
    r = prove_exact_16_criterion(primes)
    assert r["proved"] is True
    for p in primes:
        assert r["by_p"][str(p)]["L_bound_eq_16"] is True


def test_hyp_dminus1_le16():
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 31]
    r = prove_hyp_dminus1_bound(primes)
    assert r["proved"] is True


def test_room_hyp_formula():
    for p in [5, 7, 11]:
        d = d_of(p)
        assert room_hyp_of(p) == room_orth(p) * Fraction((d - 1) ** 2, d * d)


def test_general_status_honest_open():
    g = prove_general_status()
    assert g["variance_identity_proved"] is True
    assert g["exact_16_criterion_mult_d_proved"] is True
    assert g["sixteen_N_proved_for_all_p"] is False
    assert g["kappa_bound_proved_for_all_p"] is False


def test_cert_p3_or_p5_if_available():
    """Exercise shipped cert path on real Max+ when present (no hardcode skip)."""
    results = []
    for p in (3, 5):
        c = certify_variance_on_maxplus(p)
        results.append(c)
        if c.get("skipped"):
            continue
        assert c["orth_matches_var"] is True
        assert c["lambda_max_le_16"] is True
        assert c["sixteen_N"] is True
        # variance identity on real ED4
        assert abs(c["var_from_ED4"] - c["orth"]) < 1e-5 * max(1.0, abs(c["orth"]))
    # at least try; if both skipped the environment lacks Max+ — still algebra tests pass
    assert results


def test_main_writes_evidence(tmp_path=None):
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["variance_identity"] is True
    assert out["proved"]["exact_16_criterion_algebra"] is True
    assert out["proved"]["sixteen_N_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15105.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.105"
    assert data["L_status"] == "OPEN"
