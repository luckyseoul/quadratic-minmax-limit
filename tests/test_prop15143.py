"""Tests for Prop 15.143 — p=11 affine subtype law; double-switch LB."""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15143 import (  # noqa: E402
    P11_AFFINE_TYPES,
    P11_DBL_SIZE_MULT,
    P11_HPLUS_LB,
    main,
    max_gap,
    n_3ap,
    n_qr,
    phi_S,
    prove_affine_classifier_p11,
    prove_double_switch_lb,
    prove_open,
)
from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def test_phi_helpers_and_table():
    p = 11
    for row in P11_AFFINE_TYPES:
        S = row["sample_S"]
        assert phi_S(S, p) == tuple(row["phi"])
        assert n_qr(S, p) == row["n_qr"]
        assert n_3ap(S, p) == row["n3ap"]
        assert max_gap(S, p) == row["max_gap"]
        y = affine_halfspace(p, (0, 1), set(S))
        C = paley_conference_prime_power(p).astype(float)
        assert np.allclose(C @ y, p * y, atol=1e-5)
    # all phi distinct
    phis = [tuple(r["phi"]) for r in P11_AFFINE_TYPES]
    assert len(phis) == len(set(phis)) == 6


def test_affine_census_live():
    a = prove_affine_classifier_p11()
    assert a["proved"] is True
    assert a["n_affine_orbits"] == 6
    assert a["G_order"] == 14520
    assert a["size_mult"] == {"132": 1, "330": 2, "660": 3}
    for row in a["table_live"]:
        assert row["match"] is True


def test_double_switch_and_open():
    b = prove_double_switch_lb()
    assert b["proved"] is True
    assert b["ystar_orbit_size"] == 3630
    assert 7260 in b["dbl_sizes_from_ystar"]
    assert b["n_nc_on_size660"] == 0
    assert b["n_nc_on_size132"] > 0
    assert b["H_plus_lower_bound"] == P11_HPLUS_LB == 28182
    assert b["complete_type_list"] is False
    c = prove_open()
    assert c["residual_closed_general"] is False
    assert c["affine_types_complete_p11"] is True
    assert c["all_Hplus_types_complete_p11"] is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["affine_type_census_p11"] is True
    assert out["proved"]["p11_residual_maxplus_free"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15143.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.143"
    assert sum(P11_DBL_SIZE_MULT[s] * s for s in P11_DBL_SIZE_MULT) == 28182
