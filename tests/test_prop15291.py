"""Tests for Prop 15.291 — type cardinalities |τ| in p; Q_τ still open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15291 import (
    count_paley_off_formula,
    count_paley_off_formula_inverted,
    count_paley_off_live,
    live_type_counts,
    main,
    n_N1,
    n_Nstar,
    n_off,
    n_sub_off,
    prove_cardinalities,
    prove_open,
    prove_paley_slice,
)


def test_cardinalities_match_live():
    A = prove_cardinalities([3, 5, 7, 11, 13])
    assert A["proved"] is True
    assert A["invert_hit"] < A["invert_tot"]


def test_formulas_in_p():
    for p in (5, 7, 11, 13):
        live = live_type_counts(p)
        sub = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "sub")
        n1 = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "N1")
        ns = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "N*")
        assert sub == n_sub_off(p) == p - 3
        assert n1 == n_N1(p) == p - 1
        assert ns == n_Nstar(p) == (p - 1) * (p - 3) // 2
        assert sub + n1 + ns == n_off(p)


def test_invert_sub_off_is_wrong():
    # fail-when-wrong: |sub off| = p-1 includes ±1
    for p in (5, 7, 11):
        live = live_type_counts(p)
        sub = sum(v for k, v in live.items() if len(k) > 1 and k[-1] == "sub")
        assert sub != p - 1


def test_paley_slice_formula():
    B = prove_paley_slice([3, 5, 7, 11, 13])
    assert B["proved"] is True
    assert B["invert_hit"] < B["invert_tot"]


def test_paley_live_equals_formula():
    for p in (3, 5, 7, 11):
        assert count_paley_off_live(p) == count_paley_off_formula(p)


def test_paley_invert_breaks():
    for p in (5, 7):
        assert count_paley_off_formula_inverted(p) != count_paley_off_live(p)


def test_Q_tau_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["cardinalities"] is True
    assert out["proved"]["paley_slice"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
