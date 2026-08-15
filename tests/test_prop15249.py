"""Prop 15.249: algebraic free-e dual D_alg; cost via Weil; residual_i open."""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15190 import alpha_particular
from e1_gmin_m4_prop15249 import (
    We_alg,
    cost_budget,
    free_e_sc_cost_bound_proved_general,
    hinge_status_249,
    residual_i_closed_via_249,
    star_weight_scaled,
    sum_ne_alg,
    sum_ne_final_ub,
    t_ub_from_weil,
    theorem_We_and_sum_ne_alg,
    theorem_cost_bound_via_weil,
    theorem_degree_and_diag_structure,
)


def test_We_and_sum_ne_closed():
    t = theorem_We_and_sum_ne_alg()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        we = We_alg(p)
        sne = sum_ne_alg(p)
        assert we == Fraction(p**4 - 4 * p * p + 1, 2 * p * p * (p * p - 2))
        assert sne == Fraction(n, 2) / we - 1
        assert 0 < we < Fraction(1, 2)
        assert sne < Fraction(3, 2) * (n - 1)
        assert star_weight_scaled(p) > 0


def test_cost_bound_via_weil():
    t = theorem_cost_bound_via_weil()
    assert t["proved"] is True
    assert free_e_sc_cost_bound_proved_general() is True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        Sf = sum_ne_final_ub(p)
        bud = cost_budget(p)
        al = alpha_particular(p)
        assert Sf < bud
        assert al * Sf < 2 - al
        assert t_ub_from_weil(p) > 0


def test_degree_structure_stated():
    t = theorem_degree_and_diag_structure()
    assert t["proved"] is True


def test_predicates_closed_via_gplus():
    assert residual_i_closed_via_249() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert e1_closed_general() is False
    h = hinge_status_249()
    assert h["residual_i_closed"] is False
    assert h["ker_eq_scheme_cross_general"] is False
    assert h["cost_D_lt_need_via_Weil_paley"] is True


def test_numeric_construction_p5_p7():
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    for p in (5, 7):
        C = paley_conference_prime_power(p).astype(float)
        n = int(C.shape[0])
        I = np.eye(n)
        Pp = (I + C / p) / 2
        Pm = (I - C / p) / 2

        def Comm(M):
            return Pp @ M @ Pp + Pm @ M @ Pm

        W0 = np.zeros((n, n))
        W0[0, 1] = W0[1, 0] = 1.0
        far = 1.0 / (n - 3)
        for i in range(2, n):
            W0[i, i + 1 :] = far
            W0[i + 1 :, i] = far
        Sp = Comm(W0 * C)
        # degree regular before repair
        Wp = Sp * C
        np.fill_diagonal(Wp, 0.0)
        assert float(Wp.sum(0).std()) < 1e-10
        d = np.diag(Sp)
        alpha = 2 * p * p / (p * p - 1)
        S = Sp - alpha * Comm(np.diag(d))
        assert float(np.max(np.abs(np.diag(S)))) < 1e-10
        W = S * C
        np.fill_diagonal(W, 0.0)
        We = float(W[0, 1])
        assert abs(We - float(We_alg(p))) < 1e-10
        Ws = W / We
        S1 = float(Ws.sum() / 2 - 1)
        assert abs(S1 - float(sum_ne_alg(p))) < 1e-8
        # m_p bound
        factor = p * p * (p * p - 2)
        m_p = min(
            float(S[i, j] * C[i, j] * factor)
            for i in range(n)
            for j in range(i + 1, n)
        )
        assert m_p >= 1 - 2 * p - 1e-6
        # stars positive constant
        stars = np.concatenate([Ws[0, 2:], Ws[1, 2:]])
        assert float(stars.min()) > 0
        assert float(stars.std()) < 1e-10
        assert abs(float(stars.mean()) - float(star_weight_scaled(p))) < 1e-8
        # t-repair under budget
        t0 = max(0.0, -float(Ws.min()))
        assert t0 <= float(t_ub_from_weil(p)) + 1e-12
        N = n * (n - 1) / 2
        Sf = (S1 + 1 + t0 * N) / (1 + t0) - 1
        al = float(alpha_particular(p))
        assert al * Sf < float(2 - alpha_particular(p))


def test_Sp_far_ternary_and_Q_weil_p5():
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    I = np.eye(n)
    Pp = (I + C / p) / 2
    Pm = (I - C / p) / 2

    def Comm(M):
        return Pp @ M @ Pp + Pm @ M @ Pm

    W0 = np.zeros((n, n))
    W0[0, 1] = W0[1, 0] = 1.0
    far = 1.0 / (n - 3)
    for i in range(2, n):
        W0[i, i + 1 :] = far
        W0[i + 1 :, i] = far
    Sp = Comm(W0 * C)
    factor = p * p * (p * p - 2)
    vals = {
        round(float(Sp[i, j] * C[i, j] * factor), 6)
        for i in range(2, n)
        for j in range(i + 1, n)
    }
    assert vals == {1.0, float(p * p), float(2 * p * p - 1)}
    M = C[:, 0] * C[:, 1]
    Q = C @ np.diag(M) @ C
    qmax = max(
        abs(float(Q[i, j])) for i in range(2, n) for j in range(i + 1, n)
    )
    assert qmax <= 2 * p + 1e-9
