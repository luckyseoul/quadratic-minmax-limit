"""Recorded nuka GPU exact Phi(K_50)=169; smokes are n=10/26 in the log header."""
from __future__ import annotations

import json
from pathlib import Path

DIR = (
    Path(__file__).resolve().parents[1]
    / "evidence"
    / "k50_nuka_gpu_20260919"
)


def test_k50_receipt_undercuts_paley_by_six():
    rec = json.loads((DIR / "receipt.json").read_text())
    assert rec["phi_K"] == 169
    assert rec["paley_C50"] == 175
    assert rec["gap"] == 6
    assert rec["n"] == 50
    log = (DIR / "k50_gpu.log").read_text()
    assert '"phi":169.0' in log
    assert "hipblas-gfx1201" in log
