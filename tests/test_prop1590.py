"""Tests for Prop 15.90: residual bound ≡ H; pointwise κ_B; orth form."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1590 import (  # noqa: E402
    H,
    certify_bound_from_known_ray,
    max_sum_rho_kappa_from_ray,
    prove_equivalence,
    prove_orth_reformulation,
    prove_pointwise_identity,
    residual_budget,
    wick_sum_m4_kappa,
)


def test_budget_matches_H_minus_wick():
    for p in (3, 5, 7, 11, 13, 17, 19):
        budget = residual_budget(p)
        lhs = H(p) / 4 - wick_sum_m4_kappa(p)
        assert budget == lhs
        assert budget > 0


def test_equivalence_proved():
    primes = [p for p in range(3, 40) if p == 2 or (p % 2 and all(p % d for d in range(3, int(p**0.5) + 1, 2)))]
    out = prove_equivalence(primes)
    assert out["proved"] is True
    for p in (3, 5, 7):
        assert out["by_p"][str(p)]["match"] is True
        assert out["by_p"][str(p)]["n_identity"] is True


def test_max_from_ray_equals_budget_at_H():
    for p in (3, 5, 7, 11):
        assert max_sum_rho_kappa_from_ray(p, H(p)) == residual_budget(p)


def test_pointwise_identity_structure():
    pt = prove_pointwise_identity()
    assert pt["proved"] is True
    assert "y" in pt["identity"] and "B" in pt["identity"]


def test_orth_reformulation():
    orth = prove_orth_reformulation([3, 5, 7, 11])
    assert orth["proved_equivalent"] is True
    # p=5: 2 - (6+2*49/13)/26 = 2 - (6+98/13)/26
    rhs5 = Fraction(orth["by_p"]["5"]["orth_lower_bound"])
    assert rhs5 == 2 - (6 + 2 * H(5)) / 26


def test_bound_cert_p357():
    bound = certify_bound_from_known_ray()
    assert bound["all_hold"] is True
    assert bound["by_p"]["3"]["equality_with_H"] is True
    assert bound["by_p"]["5"]["equality_with_H"] is True
    assert bound["by_p"]["7"]["holds"] is True
    assert bound["by_p"]["7"]["equality_with_H"] is False


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1590.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.90"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["residual_bound_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["equivalence"]["proved"] is True
    assert data["bound_certification"]["all_hold"] is True


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.90" in sol
    idx = sol.index("15.90")
    chunk = sol[idx : idx + 3500]
    assert "OPEN" in chunk
    assert "H" in chunk
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
