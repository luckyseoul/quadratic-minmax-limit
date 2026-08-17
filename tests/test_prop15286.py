"""Tests for Prop 15.286 — n_H 1D⊙nc; n_NL=n_H at p=5.

Fail-when-wrong: Type− parents, κ=1 at p=5, n_H=n_NL at p=7.
Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15286 import (
    main,
    prove_open,
    prove_overcount,
    prove_p5_nl_is_1dnc,
)


def test_p5_nl_equals_1dnc_and_swap_fails():
    A = prove_p5_nl_is_1dnc()
    assert A["proved"] is True
    assert A["all_orbits_nH_eq_nNL"] is True
    assert A["nH_sets"] == 100
    assert A["n1D"] == 30
    assert A["type_minus_swap_changes_nH"] is True
    assert A["type_minus_also_eq"] is False


def test_overcount_p5_is_three_not_one():
    B = prove_overcount([5])
    assert B["proved"] is True
    assert len(B["blocks"]) == 1
    assert B["blocks"][0]["kappa"] == 3.0
    assert B["blocks"][0]["constant"] is True


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["n_NL_named_in_p"] is False
    assert op["n_H_named_in_p"] is False
    assert op["Nplus_named"] is False
    assert op["kappa_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["p5_nNL_eq_nH"] is True
    assert out["proved"]["overcount_constant"] is True
    assert out["proved"]["p7_nH_eq_nNL"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    assert out["proved"]["n_NL_named_in_p"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15286.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["prop"] == "15.286"
    assert data["status"]["floor"] is False
    p7 = next(b for b in data["B"]["blocks"] if b["p"] == 7)
    assert p7["kappa"] == 1.0
    assert p7["p7_nH_eq_nNL"] is False
