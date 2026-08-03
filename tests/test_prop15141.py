"""Tests for Prop 15.141 — size-12 partner; p=7 residual Max+-free."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15141 import (  # noqa: E402
    S_AFF_SIZE12,
    T_FIELD_SIZE12,
    construct_y_sharp,
    main,
    prove_all_eight_and_residual,
    prove_bitight_p5_p7,
    prove_open,
    prove_size12_partner,
)
from e1_gmin_m4_prop15139 import affine_halfspace, is_4_AP, seidel_switch  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_size12_construction():
    r = construct_y_sharp(7)
    assert r["found"] is True
    assert is_4_AP(S_AFF_SIZE12, 7)  # size-84 affine base
    assert len(T_FIELD_SIZE12) == 12
    C = paley_conference_prime_power(7).astype(float)
    y = r["y"]
    assert np.allclose(C @ y, 7 * y, atol=1e-5)
    assert abs(y[0] - 1.0) < 1e-9
    y0 = affine_halfspace(7, (0, 1), set(S_AFF_SIZE12))
    z = np.ones(50)
    for u in T_FIELD_SIZE12:
        z[1 + u] = -1.0
    C0 = seidel_switch(C, y0)
    assert np.allclose(C0 @ z, 7 * z, atol=1e-5)
    a = prove_size12_partner()
    assert a["proved"] is True


def test_all_eight_residual():
    b = prove_all_eight_and_residual()
    assert b["proved"] is True
    assert b["n_orbits"] == 8
    assert b["H_plus"] == 5726
    assert b["matches_true"] is True
    assert b["delta2_le_room"] is True
    assert b["sharp_orbit_size"] == 1176
    assert abs(b["delta2"] - float(P7_DELTA2_PAIR)) < 1e-8
    assert P7_DELTA2_PAIR <= room_hyp_over_24(7)


def test_bitight_and_open():
    c = prove_bitight_p5_p7()
    assert c["proved_form"] is True
    assert c["p7_residual_maxplus_free"] is True
    d = prove_open()
    assert d["residual_closed_general"] is False
    assert d["residual_closed_p7_maxplus_free"] is True


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["p7_residual_maxplus_free"] is True
    assert out["status"]["p7_all_eight_types_maxplus_free"] is True
    path = ROOT / "evidence" / "e1_gmin_m4_prop15141.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.141"
