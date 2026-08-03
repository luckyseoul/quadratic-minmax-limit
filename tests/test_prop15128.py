"""Tests for Prop 15.128 — full W_k census p=3,5,7; ED4 at p=7."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15124 import CERT_W, Ek1, Ek2, Ek3  # noqa: E402
from e1_gmin_m4_prop15127 import W_hoffman_int  # noqa: E402
from e1_gmin_m4_prop15128 import (  # noqa: E402
    ED4_from_W,
    KNOWN_N,
    W_CENSUS,
    load_p7_maxplus_verify,
    main,
    moments_from_W,
    prove_ED4_p7,
    prove_ED4_saturation_p3_p5,
    prove_census_consistency,
)


def test_census_consistency():
    r = prove_census_consistency()
    assert r["proved"] is True
    for p in (3, 5, 7):
        W = W_CENSUS[p]
        assert sum(W.values()) == KNOWN_N[p]
        assert W.get(p + 1) == W_hoffman_int(p)
        assert moments_from_W(W, 1) == Ek1(p)
        assert moments_from_W(W, 2) == Ek2(p)
        assert moments_from_W(W, 3) == Ek3(p)
        n = n_of(p)
        assert all(W.get(k, 0) == W.get(n - k, 0) for k in W)


def test_cert_W_match_p3_p5():
    assert W_CENSUS[3] == CERT_W[3]
    assert W_CENSUS[5] == CERT_W[5]


def test_ED4_values():
    # p=3,5 saturate
    s = prove_ED4_saturation_p3_p5()
    assert s["proved"] is True
    for p in (3, 5):
        assert ED4_from_W(W_CENSUS[p], n_of(p)) == ed4_bud_closed(p)

    # p=7 strict
    b = prove_ED4_p7()
    assert b["proved"] is True
    assert ED4_from_W(W_CENSUS[7], 50) == Fraction(12835984, 409)
    assert Fraction(12835984, 409) < ed4_bud_closed(7)
    assert Fraction(12835984, 409) > ed4_flat_closed(7)
    # structural zero
    assert W_CENSUS[7].get(10, 0) == 0
    assert W_CENSUS[7].get(40, 0) == 0


def test_live_and_main_open():
    live = load_p7_maxplus_verify()
    if not live.get("skipped"):
        assert live["W_match_census"] is True
        assert live["ED4_match"] is True

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["census_p3_p5_p7"] is True
    assert out["proved"]["ED4_p7_exact"] is True
    assert out["proved"]["weight_enumerator_closed_general"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15128.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.128"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
    assert data["census_W"]["7"]["8"] == 11
    assert data["census_W"]["7"]["24"] == 1834
