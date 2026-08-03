"""Tests for Prop 15.122 — Max+ disagreement; Aut-line; FFT budget."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of  # noqa: E402
from e1_gmin_m4_prop15119 import ew2_bud  # noqa: E402
from e1_gmin_m4_prop15122 import (  # noqa: E402
    B_fft_residual,
    certify_disagreement,
    d_max,
    is_prime,
    main,
    prove_B_budget,
    prove_lp_still_weak,
)


def test_B_budget_algebra():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_B_budget(primes)
    assert r["proved"] is True
    d = d_of(5)
    assert B_fft_residual(5) == ew2_bud(5) - d * d
    assert B_fft_residual(5) > 0


def test_lp_still_weak():
    primes = [p for p in range(3, 60) if is_prime(p)]
    r = prove_lp_still_weak(primes)
    assert r["proved"] is True
    assert r["by_p"]["5"]["too_weak"] is True


def test_dmax_formula():
    assert d_max(5) == 14
    assert d_max(7) == 34
    assert d_max(3) == 2


def test_cert_disagreement_and_main_open():
    for p in (3, 5):
        c = certify_disagreement(p)
        if c.get("skipped"):
            continue
        assert c["all_ok"] is True
        assert c["spectrum_expected_ok"] is True
        assert c["max_off_ok"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["disagreement_identity"] is True
    assert out["proved"]["fft_budget"] is True
    assert out["proved"]["hyp_residual_for_all_p"] is False
    assert out["proved"]["tight_ub_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15122.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.122"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
