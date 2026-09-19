"""Smoke: exact m_5 cycle has nonempty residue (6.20) and Paley-R slack."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from residue_620_census import census_one  # noqa: E402


def test_cycle5_residue_nonempty_and_paley_fails():
    rec = census_one({"kind": "cycle5", "n": 5})
    assert rec["Phi"] == 4.0
    assert rec["n_residue"] > 0
    assert rec["max_slack_res"] > 0.0
