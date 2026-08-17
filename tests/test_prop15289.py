"""Tests for Prop 15.289 — Paley-clique vanishing and n(ω).

Fail-when-wrong: inverted 3-fiber binomial; drop collapsing extra;
non-clique U with Hoffman mass; τ_2=0 on every Paley pair.
Floor / phi_F stay False.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15285 import n1d_set_count
from e1_gmin_m4_prop15288 import nU_formula
from e1_gmin_m4_prop15289 import (
    extra_collapsing,
    is_paley_clique,
    main,
    n_omega_formula,
    n_omega_formula_inverted,
    prove_nomega,
    prove_open,
    prove_vanish,
)


def test_nomega_sums_to_15288():
    for p, e in ((5, 0), (5, 1), (5, 2), (5, 3), (7, 0), (7, 2), (7, 3)):
        n0, n1, n2, n3 = nU_formula(p, e)
        # 1 empty + 3 singles + 3 pairs + 1 triple.  Special flags:
        # e Paley pairs, e opposite-to-vertex singles (one per Paley edge).
        got0 = n_omega_formula(p, e, 0, omega_is_special=False)
        got3 = n_omega_formula(p, e, 3, omega_is_special=False)
        # e special singles / pairs, 3-e ordinary
        got1 = e * n_omega_formula(p, e, 1, omega_is_special=True) + (
            3 - e
        ) * n_omega_formula(p, e, 1, omega_is_special=False)
        got2 = e * n_omega_formula(p, e, 2, omega_is_special=True) + (
            3 - e
        ) * n_omega_formula(p, e, 2, omega_is_special=False)
        assert got0 == n0
        assert got1 == n1
        assert got2 == n2
        assert got3 == n3
        assert got0 + got1 + got2 + got3 == n1d_set_count(p)


def test_nomega_live_and_invert():
    B = prove_nomega([5, 7])
    assert B["proved"] is True
    assert B["invert_hit"] < B["invert_tot"]
    assert B["invert_tot"] > 0
    for b in B["blocks"]:
        assert b["bad"] == 0
        assert b["checked"] > 0
    # e=3 k=0 collides; e=0 k=1 is the invert witness.
    inv = n_omega_formula_inverted(5, 0, 1, omega_is_special=False)
    assert inv != n_omega_formula(5, 0, 1, omega_is_special=False)
    assert n_omega_formula_inverted(5, 3, 0, omega_is_special=False) == n_omega_formula(
        5, 3, 0, omega_is_special=False
    )
    assert extra_collapsing(5, 3, 0, False) > 0


def test_clique_predicate():
    assert is_paley_clique(frozenset(), set()) is True
    assert is_paley_clique(frozenset([1]), set()) is True
    pairs = {frozenset((1, 2))}
    assert is_paley_clique(frozenset([1, 2]), pairs) is True
    assert is_paley_clique(frozenset([1, 3]), pairs) is False
    assert is_paley_clique(frozenset([1, 2, 3]), pairs) is False


def test_vanish_p5():
    A = prove_vanish([5])
    assert A["proved"] is True
    b = A["blocks"][0]
    assert b["nonclique_pos"] == 0
    assert b["nonclique_tot"] > 0
    assert b["paley_pair_pos"] > 0
    assert b["e3_u3_pos"] > 0
    assert b["e0_u23_pos"] == 0


def test_open_floor_untouched():
    op = prove_open()
    assert op["floor_closed"] is False
    assert op["tau_omega_named"] is False
    assert op["n_H_named_in_p"] is False
    assert op["n_NL_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_main_evidence():
    out = main()
    assert out["L_status"] == "OPEN"
    assert out["proved"]["paley_clique_vanishing"] is True
    assert out["proved"]["n_omega_formula"] is True
    assert out["proved"]["invert_is_identity"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["inequality_proved"] is False
    path = Path(__file__).resolve().parents[1] / "evidence" / "e1_gmin_m4_prop15289.json"
    data = json.loads(path.read_text())
    assert data["prop"] == "15.289"
    assert data["status"]["floor"] is False
