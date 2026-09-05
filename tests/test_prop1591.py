"""Tests for Prop 15.91: independent dual forms of H (orth / Φ / dim Z)."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1591 import (  # noqa: E402
    H,
    dim_Z,
    harm_budget,
    prove_budget_identities,
    prove_dim_Z,
    prove_orth_equivalence,
    prove_phi_equivalence,
    prove_two_sphere_implies_16N,
    qn_budget,
    residual_budget,
    sphere_baseline,
)


def _primes(lo: int = 3, hi: int = 40) -> list[int]:
    out = []
    for p in range(lo, hi):
        if p > 2 and p % 2 == 0:
            continue
        if all(p % d for d in range(3, int(p**0.5) + 1, 2)):
            out.append(p)
    return out


def test_dim_Z_formula():
    assert dim_Z(3) == 5
    assert dim_Z(5) == 65
    assert dim_Z(7) == 275
    out = prove_dim_Z(_primes())
    assert out["proved"] is True
    for p in (3, 5, 7, 11):
        assert out["by_p"][str(p)]["match"] is True
        assert out["by_p"][str(p)]["G_full_rank"] is True


def test_residual_budget_matches_qn_minus_wick():
    for p in _primes(3, 50):
        assert residual_budget(p) == qn_budget(p) - 8
        assert residual_budget(p) == Fraction((p + 1) * (p + 7), (p * p + 1) // 2)
        assert residual_budget(p) == 2 * (H(p) - 1)


def test_sphere_harm_budget():
    for p in _primes(3, 40):
        assert harm_budget(p) == qn_budget(p) - sphere_baseline(p)
        assert sphere_baseline(p) < 8
        assert harm_budget(p) > 0


def test_budget_chain_proved():
    out = prove_budget_identities(_primes(3, 55))
    assert out["proved"] is True
    # p=3 equality to 16; p>=5 strict
    assert Fraction(out["by_p"]["3"]["qn_budget"]) == 16
    for p in (5, 7, 11, 13):
        qn = Fraction(out["by_p"][str(p)]["qn_budget"])
        assert qn < 16
        assert out["by_p"][str(p)]["chain_sphere_wick_qn"] is True


def test_orth_and_phi_equivalence():
    primes = _primes()
    orth = prove_orth_equivalence(primes)
    phi = prove_phi_equivalence(primes)
    assert orth["proved_equivalent"] is True
    assert phi["proved_equivalent"] is True
    # p=5 orth LB = 2 - (176/13)/26 = 2 - 176/(13*26)
    n = 26
    lb = Fraction(orth["by_p"]["5"]["orth_lower_bound"])
    assert lb == 2 - qn_budget(5) / n


def test_two_sphere_implies_16N():
    out = prove_two_sphere_implies_16N(_primes())
    assert out["proved"] is True
    # 2×sphere does not imply H at p=5 (strictly weaker)
    assert out["by_p"]["5"]["2sphere_implies_H"] is False
    assert out["by_p"]["5"]["implies_qn_lt_16"] is True


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1591.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.91"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["H_proved_for_all_p_ge_5"] is False
    assert data["proved"] is True
    assert data["dim_Z"]["proved"] is True
    assert data["budget_identities"]["proved"] is True
    assert data["orth_cert"]["3"]["ok"] is True
    assert data["orth_cert"]["5"]["ok"] is True
    # attack surface must name independent targets, not only residual ∑ρ
    assert "orth" in data["attack_surface"].lower() or "Orth" in data["attack_surface"]
    assert "κ" in data["attack_surface"] or "kappa" in data["attack_surface"].lower() or "4th" in data["attack_surface"]


def test_solution_handoff():
    sol = (ROOT / "solution.md").read_text()
    assert "15.91" in sol
    idx = sol.index("15.91")
    chunk = sol[idx : idx + 4000]
    assert "OPEN" in chunk
    assert "H" in chunk
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
