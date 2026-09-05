"""Tests for Prop 15.97: Veronese mult ID + Ky Fan criterion."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop1597 import (  # noqa: E402
    prove_aut_schur,
    prove_ky_fan_criterion,
    prove_veronese_identification,
)


def test_veronese_proved():
    out = prove_veronese_identification()
    assert out["proved"] is True
    assert "Veronese" in out["theorem"] or "φ" in out["theorem"] or "Gamma" in out["theorem"] or "Γ" in out["theorem"]


def test_ky_fan_proved():
    out = prove_ky_fan_criterion()
    assert out["proved"] is True
    assert "Fan" in out["theorem"] or "orthonormal" in out["theorem"]


def test_aut_schur_proved():
    out = prove_aut_schur()
    assert out["proved"] is True
    assert "Schur" in out["theorem"] or "Aut" in out["theorem"]


def test_evidence_json():
    path = ROOT / "evidence" / "e1_gmin_m4_prop1597.json"
    assert path.is_file(), "run src/e1_gmin_m4_prop1597.py first"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.97"
    assert data["L_status"] == "OPEN"
    assert data["H_proved"] is False
    assert data["gap_proved_for_all_p"] is False
    assert data["mult_proved_for_all_p"] is False
    assert data["proved"] is True
    assert data["veronese_identification"]["proved"] is True
    assert data["ky_fan_criterion"]["proved"] is True

    for pk in ("3", "5"):
        c = data["cert"][pk]
        assert c.get("skipped") is not True
        assert c["mult_eq_d"] is True
        assert c["mult_Gamma_eq_d"] is True
        assert c["spectra_match"] is True
        assert c["ky_fan_equal_rays"] is True
        assert c["ok"] is True
        assert c["mult_PopP"] == c["d"]
        assert c["mult_Gamma_top"] == c["d"]


def test_solution_handoff_graph():
    sol = (ROOT / "solution.md").read_text()
    assert "15.97" in sol
    idx = sol.index("15.97")
    assert "OPEN" in sol[idx : idx + 4000]
    hand = (ROOT / "HANDOFF.md").read_text()
    assert "**Current mathematical status:** **L OPEN.**" in hand[:5000]
    graph = (ROOT / "evidence" / "P0_ENGINEERING_GRAPH.md").read_text()
    assert "15.97" in graph
