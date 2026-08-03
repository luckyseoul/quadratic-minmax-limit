"""Tests for Prop 15.124 — weight moments; E[k⁴] partition; residual as R₄."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, m4f_bud  # noqa: E402
from e1_gmin_m4_prop15124 import (  # noqa: E402
    CERT_W,
    Ek1,
    Ek2,
    Ek3,
    R4_bud,
    certify_moments_and_partition,
    ed4_from_exact3,
    exact_le3,
    is_prime,
    main,
    prove_delsarte_weak,
    prove_moments_le3,
    prove_partition,
    prove_residual_dictionary,
)


def test_moments_le3_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_moments_le3(primes)
    assert r["proved"] is True
    assert Ek1(5) == Fraction(13)
    assert Ek2(5) == Fraction(26 * 28, 4)
    assert Ek3(5) == Fraction(26 * 26 * 32, 8)  # n²(n+6)/8 = 2704
    assert Ek3(3) == 200


def test_partition_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_partition(primes)
    assert r["proved"] is True
    n = 10
    # exact_le3 at p=3
    ex = exact_le3(3)["exact_le3"]
    assert ex == Fraction(n, 2) + n * n + Fraction(3 * n * n, 4) + Fraction(
        3 * n * (n - 2) * (n + 2), 4
    )


def test_residual_dictionary():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_residual_dictionary(primes)
    assert r["proved"] is True
    # budget consistency
    p = 5
    assert ed4_from_exact3(p) + 16 * R4_bud(p) == ed4_bud_closed(p)


def test_cert_and_main_open():
    for p in (3, 5):
        c = certify_moments_and_partition(p)
        assert c["E_k_match"]
        assert c["E_k2_match"]
        assert c["E_k3_match"]
        assert c["ED4_via_R4_match"]
        assert c["ED4_le_bud"]
        assert c["R4_le_bud"]
        assert c["A4_eq_m4f_bud"]  # saturation at p=3,5
        assert c["W_p_plus_1_eq_d"]
        assert c["W_p_plus_1"] == d_of(p)

    d = prove_delsarte_weak([3, 5, 7])
    assert d["proved_weak_for_p_ge_5"] is True
    assert d["by_p"]["5"]["too_weak_for_residual"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["moments_le3_maxplus_free"] is True
    assert out["proved"]["exact_partition_k4"] is True
    assert out["proved"]["residual_dictionary_R4"] is True
    assert out["proved"]["weight_enumerator_closed"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15124.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.124"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    # CERT_W totals
    assert sum(CERT_W[3].values()) == 12
    assert sum(CERT_W[5].values()) == 260
