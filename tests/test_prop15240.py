"""Tests for Prop 15.240 — maj≤2/n; envelope criterion; residual open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15240 import (  # noqa: E402
    hinge_status_240,
    maj_le_two_over_n_diff,
    residual_i_closed_via_240,
    theorem_envelope_criterion,
    theorem_free_e_sc_budget_recall,
    theorem_maj_le_two_over_n,
)


def test_maj_le_two_over_n_proved():
    assert theorem_maj_le_two_over_n()["proved"] is True


def test_envelope_criterion_chain_but_hyp_open():
    t = theorem_envelope_criterion()
    assert t["proved_criterion"] is True
    assert t["hypothesis_proved_general"] is False
    assert t["sufficient_for_residual_i"] is False


def test_free_e_structure_ker_open():
    t = theorem_free_e_sc_budget_recall()
    assert t["proved_structure"] is True
    assert t["ker_sc_proved_general"] is False


def test_residual_i_still_open():
    assert residual_i_closed_via_240() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_240()
    assert h["maj_le_two_over_n"] is True
    assert h["envelope_hypothesis"] is False
    assert h["residual_i_closed_via_240"] is False


def test_maj_le_two_n_values():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101, 103):
        assert is_prime(p)
        maj = majorant_mu_part(p)
        twn = two_over_n(p)
        assert maj <= twn
        assert maj_le_two_over_n_diff(p) == twn - maj
        assert maj_le_two_over_n_diff(p) > 0


def test_cubic_at_p5():
    assert (5**3 - 3 * 25 - 5 * 5 - 9) == 16
