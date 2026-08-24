"""Tests for Prop 15.628 — affine circle completion closes W1/W2/Walsh."""
from __future__ import annotations

from e1_gmin_m4_prop15628 import (
    theorem_A_eligible_span,
    theorem_B_affine_completions,
    theorem_C_walsh_close,
)
from w2_affine_circle_close import run


def test_theorem_flags_and_open_scope():
    assert theorem_A_eligible_span()["proved"] is True
    assert theorem_B_affine_completions()["proved"] is True
    close = theorem_C_walsh_close()
    assert close["proved"] is True
    assert close["W1"] is True
    assert close["W2"] is True
    assert close["Walsh_15_406_E"] is True
    assert "L" in close["not_closed"]


def test_exact_small_and_p19_witnesses():
    out = run([5, 7, 19])
    assert out["all_checks"] is True
    assert out["p3_direct"]["direction_rank"] == 4
    assert out["p3_direct"]["full_slice"] is True
    assert out["rows"]["5"]["all_affine_subsets"] == {
        "tested": 10,
        "failures": 0,
    }
    assert out["rows"]["7"]["all_affine_subsets"] == {
        "tested": 35,
        "failures": 0,
    }
    p19 = out["rows"]["19"]
    assert p19["all_outside_pairs"]["tested"] == 58311
    assert p19["all_outside_pairs"]["selection_failures"] == 0
    assert p19["frobenius_witness"]["frobenius_pair"] == [2, 340]
    assert p19["frobenius_witness"]["chosen_subset"] == [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        18,
    ]
    assert p19["frobenius_witness"]["pair_in_U_before_and_after"] is True
