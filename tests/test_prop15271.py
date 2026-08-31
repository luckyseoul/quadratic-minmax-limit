"""Prop 15.271: k=3 Veronese ker⊆F reduction; residual (i) stays open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15270 import gplus_pd_proved_general, psl_span_F_eq_Wpp0
from e1_gmin_m4_prop15271 import (
    certify_fperp_injective,
    dim_F,
    dim_Wpp0,
    fperp_injective_proved_general,
    hinge_status_271,
    k3_spans_Wpp0_proved_general,
    residual_i_closed_via_271,
    theorem_fibre_dft,
    theorem_majority_identity,
    theorem_p5_k1_fills,
    theorem_phase_frequencies_distinct,
    theorem_phase_lock,
    theorem_singer_pd_on_F_p_ge_7,
    theorem_slice_sizes,
)


def test_dims():
    assert dim_F(5) == 5
    assert dim_F(7) == 11
    assert dim_F(11) == 29
    assert dim_Wpp0(5) == 65
    assert dim_Wpp0(7) == 275
    assert dim_Wpp0(11) == 1769
    assert dim_Wpp0(5) - dim_F(5) == 60


def test_structure_lemmas():
    assert theorem_majority_identity()["proved"] is True
    assert theorem_phase_lock()["proved"] is True
    assert theorem_fibre_dft()["proved"] is True
    assert theorem_slice_sizes()["proved"] is True
    assert theorem_phase_frequencies_distinct()["proved"] is True
    assert theorem_singer_pd_on_F_p_ge_7()["proved"] is True
    assert theorem_p5_k1_fills()["proved"] is True


def test_census_not_a_proof():
    c = certify_fperp_injective()
    assert c["certified"] is True
    assert c["proved_general"] is False
    assert c["by_p"]["5"]["rank_on_Fperp"] == 60
    assert c["by_p"]["7"]["rank_total"] == 275
    assert c["by_p"]["11"]["rank_on_Fperp"] == 1740
    assert fperp_injective_proved_general() is False
    assert k3_spans_Wpp0_proved_general() is False


def test_live_flags_open():
    assert residual_i_closed_via_271() is False
    assert gplus_pd_proved_general() is True
    assert psl_span_F_eq_Wpp0() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_271()
    assert h["fperp_injective_general"] is False
    assert h["residual_i_closed"] is False
    assert h["e1_closed_general"] is False
