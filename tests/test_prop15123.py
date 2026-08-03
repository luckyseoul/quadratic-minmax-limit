"""Tests for Prop 15.123 — switching; srg regular sets; dual residual."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_bud, m4f_flat  # noqa: E402
from e1_gmin_m4_prop15123 import (  # noqa: E402
    allowed_k,
    certify_switching_and_weights,
    is_prime,
    main,
    prove_allowed_k,
    prove_dual_residual,
    prove_srg_params,
    srg_degree,
    srg_lambda,
    srg_mu,
    srg_mult_tau,
    srg_mult_theta,
    srg_tau,
    srg_theta,
)


def test_srg_params_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_srg_params(primes)
    assert r["proved"] is True
    assert srg_degree(5) == 10
    assert srg_mu(5) == 4
    assert srg_lambda(5) == 3
    assert srg_theta(5) == Fraction(2)
    assert srg_tau(5) == Fraction(-3)
    assert srg_mult_theta(5) == d_of(5)
    assert srg_mult_tau(5) == d_of(5) - 1


def test_allowed_k():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_allowed_k(primes)
    assert r["proved"] is True
    assert allowed_k(3) == [0, 4, 6, 10]
    assert allowed_k(5) == [0, 6, 8, 10, 12, 14, 16, 18, 20, 26]
    # alpha - beta = tau
    k = 8
    p = 5
    alpha = Fraction(k - 1 - p, 2)
    beta = Fraction(k, 2)
    assert alpha - beta == srg_tau(p)


def test_dual_residual_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_dual_residual(primes)
    assert r["proved"] is True
    # m4f_bud = m4f_flat + room/24 is A'4 budget
    assert m4f_bud(5) == m4f_flat(5) + Fraction(1536, 65)


def test_cert_and_main_open():
    for p in (3, 5):
        c = certify_switching_and_weights(p)
        if c.get("skipped"):
            continue
        assert c["C2_ones_evec"]
        assert c["ones_in_Max"]
        assert c["w_evec_ok"]
        assert c["B_equals_W"]
        assert c["srg_params_ok"]
        assert c["weights_in_allowed_k"]
        assert c["ED4_le_bud"]

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["srg_params"] is True
    assert out["proved"]["dual_krawtchouk_residual"] is True
    assert out["proved"]["weight_enumerator_closed"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15123.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.123"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
