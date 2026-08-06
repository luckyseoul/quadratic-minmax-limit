"""Drive shipped predicates: Gsum disj LB is the fatal hinge; L OPEN until proved."""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15167 import bitight_from_majorization
from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    gsum_disj_lb_status,
    gsum_farkas_poly,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed, main_L_from_e1

ROOT = Path(__file__).resolve().parents[1]


def test_maxplus_not_ip_scheme_blocks_scheme_gsum_lb():
    """15.158 evidence: Max+ is not an IP association scheme — scheme-min Gsum invalid."""
    data = json.loads((ROOT / "evidence" / "e1_gmin_m4_prop15158.json").read_text())
    assert data["proved"]["not_IP_association_scheme"] is True
    # Shipped status must not re-validate scheme_justification
    assert gsum_disj_lb_status()["scheme_justification_valid"] is False


def test_gsum_disj_lb_honest_false():
    assert gsum_disj_lb_proved_general() is False
    st = gsum_disj_lb_status()
    assert st["proved_adjacent_zero"] is True
    assert st["proved_disj_lb_general"] is False
    assert st["scheme_justification_valid"] is False
    assert st["certified_disj_at_p5"] is True


def test_farkas_poly_still_positive_when_lb_assumed():
    """Algebra of need < k·LB under candidate LB is real Fraction."""
    from e1_gmin_m4_prop15170 import dual_equality_farkas_algebra_if_lb, dual_equality_impossible_general

    for p in (5, 7, 11, 13, 17, 19, 23):
        assert gsum_farkas_poly(p) > 0
    assert dual_equality_farkas_algebra_if_lb() is True
    assert dual_equality_impossible_general() is False  # gated; no soft-close


def test_e1_and_L_open_until_hinge_closed():
    assert type_I_k_3p_minus_2_closed_general() is False
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is False
    assert e1_closed_general() is False
    bt = bitight_from_majorization(5)["bitight_empty"]
    assert bt is True
    w = main_L_from_e1(e1=False, bitight=bt)
    assert w["L_closed"] is False
    assert w["L_status"] == "OPEN"
