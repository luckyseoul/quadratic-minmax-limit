"""n=9 two-half CP-SAT: B<=34 infeasible, so zero-error diamond fails."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from original_mo_two_half_geometry import boolean_states, quadratic_values  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "two_half_n9_n10_20260919"


def test_a9_has_phi_12():
    a = np.load(EV / "A9_exact.npy")
    q = quadratic_values(a, boolean_states(9, fix_first=True))
    assert int(np.max(np.abs(q))) == 12


def test_n9_zero_error_diamond_infeasible_at_34():
    rec32 = json.loads((EV / "n9_feas32.json").read_text())
    rec34 = json.loads((EV / "n9_feas34.json").read_text())
    assert rec32["status"] == "INFEASIBLE"
    assert rec34["status"] == "INFEASIBLE"
    assert rec34["min_B"] == 36
