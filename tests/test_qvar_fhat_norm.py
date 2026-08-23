"""F̂ is not a Paley-field square or field-norm.  QVAR still OPEN."""
from __future__ import annotations

import inspect
from fractions import Fraction

import e1_gmin_global_qvar as G
import e1_gmin_leftover1_qvar_principal as L
import e1_gmin_qvar_fhat_norm as F
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general


def test_kronecker_and_orders():
    assert F.kronecker(5, 13) == -1
    assert F.kronecker(-7, 409) == -1
    assert F.mul_order_mod(13, 24) == 2
    assert F.mul_order_mod(25, 48) == 2
    assert F.quadratic_f(13, 5) == 2
    assert F.quadratic_f(409, -7) == 2
    assert F.cyclotomic_f(13, 24) == 2
    assert F.cyclotomic_f(409, 48) == 2


def test_fhat_fractions():
    assert F.fhat_from_lambda(5, Fraction(80, 13)) == Fraction(1250, 13)
    assert not F.is_q_square(Fraction(1250, 13))
    q7 = 49
    assert F.fhat_from_lambda(7, Fraction(3072, 409)) == (
        Fraction(3072, 409) - 6
    ) * q7 * q7


def test_kill_theorem():
    T = F.theorem_fhat_not_paley_field_norm()
    assert T["proved"] is True
    assert T["inequality_proved"] is False
    assert T["claim_Fhat_is_Q_square"] is False
    assert T["claim_Fhat_is_paley_field_norm"] is False
    assert T["p5_min_Fhat"] == "1250/13"
    for row in T["rows"]["5"] + T["rows"]["7"]:
        assert row["is_Q_square"] is False
        assert row["odd_val_blocks_norm"] is True
    assert 409 in T["rows"]["7"][0]["den_outside_paley_ramification"] or any(
        409 in r["den_outside_paley_ramification"] for r in T["rows"]["7"]
    )


def test_fail_when_claim_square_or_norm():
    T = F.theorem_fhat_not_paley_field_norm()
    min5 = next(r for r in T["rows"]["5"] if r["k"] == 4)
    min7 = next(r for r in T["rows"]["7"] if r["k"] == 16)
    assert min5["cyclotomic_f_at_den"]["13"] == 2
    assert min7["cyclotomic_f_at_den"]["409"] == 2
    assert F.valuation(Fraction(min5["Fhat"]), 13) == -1
    assert F.valuation(Fraction(min7["Fhat"]), 409) == -1


def test_does_not_flip_qvar_or_leftovers():
    assert G.global_qvar_proved_general() is False
    assert L.leftover1_qvar_and_principal_proved() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    src = inspect.getsource(G.global_qvar_proved_general)
    assert "return True" not in src
    src_t = inspect.getsource(F.theorem_fhat_not_paley_field_norm)
    assert "inequality_proved" in src_t
    assert "False" in src_t
