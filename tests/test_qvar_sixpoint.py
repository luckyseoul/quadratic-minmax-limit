"""6-point linear + SOS-4 cannot prove GLOBAL QVAR; flags stay False."""
from __future__ import annotations

import inspect

import e1_gmin_global_qvar as G
import e1_gmin_leftover1_qvar_principal as L
import e1_gmin_qvar_sixpoint as S
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

_T = None


def _thm():
    global _T
    if _T is None:
        _T = S.theorem_sixpoint_and_sos4_cannot_prove_qvar()
    return _T


def test_sixpoint_box_lp_min_is_negative():
    T = _thm()
    assert T["proved"]
    assert T["inequality_proved"] is False
    assert T["claim_sixpoint_min_ge_0"] is False
    assert T["p5_sixpoint"]["lp_ok"]
    assert T["p5_sixpoint"]["ker_dim"] >= 1
    assert T["p5_sixpoint"]["min_pairing"] < 0
    assert abs(T["p5_sixpoint"]["min_pairing"] - (-101 / 4)) < 1e-6
    assert T["p5_true"]["pairing"] > 0
    assert T["p5_true"]["pairing"] > T["p5_sixpoint"]["min_pairing"]


def test_sos4_edge_pairing_is_negative():
    T = _thm()
    assert T["claim_sos4_min_ge_0"] is False
    assert T["p5_sos4"]["min_edge_pairing"] < 0
    assert abs(T["p5_sos4"]["min_edge_pairing"] - (-45 / 4)) < 1e-6
    assert T["p5_sos4"]["true_pairing"] > 0
    assert T["p7_joint_deg6_recorded"]["imported_as_p_law"] is False
    assert T["p7_joint_deg6_recorded"]["kernel_dim"] == 4


def test_does_not_flip_qvar_or_leftovers():
    T = _thm()
    assert G.global_qvar_proved_general() is False
    assert L.leftover1_qvar_and_principal_proved() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    src = inspect.getsource(G.global_qvar_proved_general)
    assert "return True" not in src
    assert T["inequality_proved"] is False
