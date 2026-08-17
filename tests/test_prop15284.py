"""Tests for Prop 15.284 — dim-2 3-set type (Paley, minpoly).

Fail-when-wrong: pal+Tr and pal+N must NOT determine n.
Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15284 import main, prove_minpoly_type, prove_open


def test_minpoly_determines_n_partials_split():
    a = prove_minpoly_type([5, 7])
    assert a["proved"] is True
    assert a["pal_Tr_splits"] is True
    assert a["pal_N_splits"] is True
    for b in a["blocks"]:
        assert b["pal_Tr_N_split"] == 0
        assert b["n_dim2"] > 0
        if b["p"] == 7:
            assert b["pal_Tr_split"] > 0
            assert b["pal_N_split"] > 0


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["n_values_named_in_p"] is False
    assert op["Nplus_named"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["minpoly_determines_n"] is True
    assert out["proved"]["pal_Tr_is_identity"] is False
    assert out["proved"]["pal_N_is_identity"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15284.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.284"
    assert data["status"]["floor"] is False
