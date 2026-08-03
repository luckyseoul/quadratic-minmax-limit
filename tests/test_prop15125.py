"""Tests for Prop 15.125 — perfect 2-colorings; 4-design defect; R₄ budget."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed  # noqa: E402
from e1_gmin_m4_prop15124 import R4_bud, ed4_from_exact3  # noqa: E402
from e1_gmin_m4_prop15125 import (  # noqa: E402
    certify_hoffman_and_moments,
    ed4_from_exact3_factored,
    is_prime,
    main,
    prove_R4_budget,
    prove_antipodal_dual,
    prove_design_defect,
    prove_lp_weak,
    prove_perfect_2coloring,
    regular_alpha,
    regular_beta,
    sphere_Es4,
    srg_tau,
)


def test_perfect_2coloring():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_perfect_2coloring(primes)
    assert r["proved"] is True
    p, k = 5, 8
    assert regular_alpha(p, k) - regular_beta(p, k) == srg_tau(p)
    assert regular_alpha(5, 6) == 0  # Hoffman


def test_R4_budget_factored():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_R4_budget(primes)
    assert r["proved"] is True
    for p in (3, 5, 7, 11):
        assert ed4_from_exact3(p) == ed4_from_exact3_factored(p)
        assert ed4_from_exact3(p) + 16 * R4_bud(p) == ed4_bud_closed(p)


def test_design_defect_and_lp():
    primes = [p for p in range(3, 40) if is_prime(p)]
    b = prove_design_defect(primes)
    assert b["proved"] is True
    # p=5: empirical ED4 exceeds 4-design sphere value
    assert "5" in b["cert"]
    assert b["cert"]["5"]["exceeds_4design"] is True
    assert b["cert"]["5"]["ED4_eq_bud"] is True

    d = prove_lp_weak([3, 5, 7])
    assert d["proved_weak_for_p_ge_5"] is True
    assert d["by_p"]["5"]["too_weak"] is True
    assert d["by_p"]["3"]["A4_max_le_bud"] is True


def test_cert_and_main_open():
    e = prove_antipodal_dual([3, 5])
    assert e["proved"] is True

    for p in (3, 5):
        c = certify_hoffman_and_moments(p)
        assert c["moments_match"]
        assert c["R4_eq_bud"]
        assert c["W_p_plus_1_eq_d"]
        assert c["W_p_plus_1"] == d_of(p)

    # sphere formula positive
    assert sphere_Es4(5) == Fraction(3 * 26**4, 13 * 15)

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["perfect_2colorings"] is True
    assert out["proved"]["R4_budget_closed"] is True
    assert out["proved"]["delsarte_moments_weak_p_ge_5"] is True
    assert out["proved"]["weight_enumerator_closed"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15125.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.125"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
