"""Tests for Prop 15.283 — n(S) Aut_∞-invariant; Paley type does not determine n.

Fail-when-wrong.  phi_F_ge_6 / e1 / L stay False / untouched.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15283 import (
    main,
    prove_aut_invariance,
    prove_dim1_ap,
    prove_open,
)


def test_aut_constant_and_paley_does_not_determine():
    a = prove_aut_invariance([5, 7])
    assert a["proved"] is True
    assert a["paley_type_determines_n"] is False
    for block in a["rows"]:
        for r in block.values():
            assert r["aut_constant"] is True
            assert r["aut_spanning"] == 0
            assert r["paley_type_constant"] is False
            assert r["paley_spanning_types"] > 0


def test_dim1_ap_names_n_dim2_splits():
    b = prove_dim1_ap([5, 7])
    assert b["proved"] is True
    assert b["dim2_paley_does_not_determine"] is True
    for r in b["rows"]:
        assert r["dim1_split"] == 0
        assert r["dim2_pal_split"] > 0


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["Nplus_named"] is False
    assert op["mixture_meets_floor"] is False
    assert op["n_S_named_formula"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["aut_invariance_k_le_4"] is True
    assert out["proved"]["paley_type_determines_n"] is False
    assert out["proved"]["dim1_ap_names_n"] is True
    assert out["proved"]["dim2_paley_does_not_determine"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15283.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.283"
    assert data["status"]["floor"] is False
