"""Tests for Prop 15.290 — field-norm split of Q; floor is ⟨δ,ψ⟩≤2."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general
from e1_gmin_m4_prop15290 import (
    field_norm,
    in_subfield,
    live_Q_on_T,
    main,
    norm_type,
    paley_only_key,
    paley_pair,
    prove_deficit_bound_certified,
    prove_open,
    prove_type_constancy,
    type_key,
)


def test_norm_type_pm1():
    for p in (5, 7):
        F = _kernel_field(p)
        assert norm_type(F, 1) == "pm1"
        assert norm_type(F, F["neg1"]) == "pm1"


def test_subfield_squares_are_sub():
    F = _kernel_field(5)
    # F_5^× squares among 1,2,3,4: 1,4
    assert in_subfield(F, 1)
    assert in_subfield(F, 4)


def test_paley_only_splits_at_p7():
    F = _kernel_field(7)
    Q = live_Q_on_T(7)
    q2 = float(F["q"] ** 2)
    by = {}
    for r, Qr in Q.items():
        by.setdefault(paley_only_key(F, r), []).append(Qr)
    spreads = [(max(vs) - min(vs)) / q2 for vs in by.values()]
    assert max(spreads) > 1e-4


def test_norm_type_constant_Q():
    A = prove_type_constancy([5, 7])
    assert A["proved"] is True
    assert A["paley_class_splits"] >= 1


def test_N1_not_Nstar_on_minusminus():
    F = _kernel_field(7)
    Q = live_Q_on_T(7)
    q2 = float(F["q"] ** 2)
    vals = {}
    for r, Qr in Q.items():
        if paley_pair(F, r) != (-1, -1):
            continue
        vals.setdefault(norm_type(F, r), []).append(Qr / q2)
    assert "N1" in vals and "N*" in vals
    assert abs(sum(vals["N1"]) / len(vals["N1"]) - sum(vals["N*"]) / len(vals["N*"])) > 1e-4


def test_deficit_certified():
    B = prove_deficit_bound_certified()
    assert B["proved"] is True
    assert B["by_p"]["5"]["bound_2"] is True
    assert B["by_p"]["7"]["bound_2"] is True


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["type_constancy_p5_p7"] is True
    assert out["proved"]["paley_only_splits"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
