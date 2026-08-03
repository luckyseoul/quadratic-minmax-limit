"""Tests for Prop 15.144 — free orbits; type-enum dead; residual redirect."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15144 import (  # noqa: E402
    FREE_PATH,
    SIZE2420_PATH,
    main,
    prove_free_orbit,
    prove_open,
    prove_redirect,
    prove_size2420_chain,
    prove_type_enum_dead,
)


def test_free_orbit():
    a = prove_free_orbit()
    assert a["proved"] is True
    assert a["final_equals_G"] is True
    assert a["sizes_along_path"][0] == FREE_PATH["base_orbit_size"] == 132
    assert a["sizes_along_path"][-1] == 14520
    assert a["G_order"] == 14520


def test_size2420_chain():
    b = prove_size2420_chain()
    assert b["proved"] is True
    assert b["sizes_along_path"][-1] == SIZE2420_PATH["final_orbit_size"] == 2420
    # intermediate free orbit in chain
    assert 14520 in b["sizes_along_path"]


def test_type_enum_dead_and_redirect():
    c = prove_type_enum_dead()
    assert c["proved"] is True
    assert 14520 in c["nonaffine_sizes_seen"]
    assert 2420 in c["nonaffine_sizes_seen"]
    d = prove_redirect()
    assert d["proved_structure"] is True
    assert d["rho_lt_room_p11"] is True
    room = room_hyp_over_24(11)
    rho = Fraction(d["rho_min2_p11"])
    assert rho < room
    e = prove_open()
    assert e["residual_closed_general"] is False
    assert e["type_list_residual_viable_p11"] is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["residual_closed"] is False
    assert out["proved"]["free_orbit_p11"] is True
    assert out["proved"]["type_enum_dead_p_ge_11"] is True
    assert out["proved"]["p11_residual_maxplus_free"] is False
    path = ROOT / "evidence" / "e1_gmin_m4_prop15144.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.144"
