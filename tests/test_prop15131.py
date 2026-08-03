"""Tests for Prop 15.131 — pair-avg vs basepoint; p=7 Q_δ spectrum; pointwise."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_flat_closed  # noqa: E402
from e1_gmin_m4_prop15128 import ED4_from_W, W_CENSUS  # noqa: E402
from e1_gmin_m4_prop15131 import (  # noqa: E402
    P7_DELTA2_PAIR,
    P7_ED4_PAIR,
    P7_MAX_Q,
    P7_MIN_Q,
    P7_N,
    P7_SLACK_POINTWISE,
    P7_TYPES,
    P7_VAR_Q,
    main,
    prove_pair_vs_basepoint,
    prove_pointwise_p7,
    prove_spectrum_p7,
    prove_true_delta2_p7,
    prove_var_Q_p7,
)


def test_pair_vs_basepoint():
    r = prove_pair_vs_basepoint()
    assert r["proved"] is True
    assert r["p7"]["non_constant"] is True
    # W-based = middle type, not pair
    ED4_hs = ED4_from_W(W_CENSUS[7], 50)
    assert (ED4_hs - ed4_flat_closed(7)) / 24 == Fraction(82176, 4499)
    assert Fraction(82176, 4499) != P7_DELTA2_PAIR


def test_spectrum_three_types():
    r = prove_spectrum_p7()
    assert r["proved"] is True
    assert r["n_types"] == 3
    assert sum(t["count"] for t in P7_TYPES) == P7_N
    assert P7_MIN_Q < 0
    assert P7_MAX_Q > 0
    assert 4499 == 11 * 409
    assert P7_N == 28 * 409


def test_true_delta2_and_pointwise():
    c = prove_true_delta2_p7()
    assert c["proved"] is True
    assert P7_DELTA2_PAIR == Fraction(19180800, 1840091)
    assert P7_ED4_PAIR == Fraction(5218435600, 167281)
    assert P7_DELTA2_PAIR <= rho_min2(7)
    assert P7_DELTA2_PAIR <= room_hyp_over_24(7)
    # ED4_pair = flat + 24 d2
    assert ed4_flat_closed(7) + 24 * P7_DELTA2_PAIR == P7_ED4_PAIR

    d = prove_pointwise_p7()
    assert d["proved"] is True
    assert P7_MAX_Q <= rho_min2(7)
    assert rho_min2(7) - P7_MAX_Q == P7_SLACK_POINTWISE
    assert P7_MAX_Q / rho_min2(7) == Fraction(613872, 664625)


def test_var_and_main_open():
    e = prove_var_Q_p7()
    assert e["proved"] is True
    assert P7_VAR_Q > 0

    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["pair_vs_basepoint"] is True
    assert out["proved"]["spectrum_p7_three_types"] is True
    assert out["proved"]["true_delta2_p7"] is True
    assert out["proved"]["pointwise_maxQ_le_rho_min_p7"] is True
    assert out["proved"]["delta_le_rho_min_for_all_p"] is False
    assert out["proved"]["hyp_residual_for_all_p"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15131.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.131"
    assert data["L_status"] == "OPEN"
    assert "No soft-close" in data["settlement_note"]
