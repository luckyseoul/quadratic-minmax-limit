"""Tests for Prop 15.456 — remaining 4410 = 4u+2u+9u; 15u dies at p=5."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15455 import n_CRRR_named, n_affJ2
from e1_gmin_m4_prop15456 import (
    main,
    n_J4_all_pairs_equal_wrong,
    n_J4_from_pairs,
    n_J4_p7,
    n_J4_pair_ordinary,
    n_J4_pair_special,
    n_J6_p7,
    n_NL_false_19u,
    n_affJ2_p7,
    n_rest_p7,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    u_unit,
)


def test_affJ2_not_J6():
    A = prove_A()
    assert A["proved"] is True
    u = u_unit(7)
    assert n_affJ2_p7(7) == 4 * u == 1176
    assert n_J6_p7(7) == 2 * u == 588
    assert n_affJ2_p7(7) != n_J6_p7(7)
    assert A["per_slot_aff"] == 294
    assert A["per_slot_J6"] == 147
    # residual still aff+J2, not a J2
    assert tuple(sorted(n_affJ2(7, 1, 0, 6, 0))) == (0, 1, 1, 2, 5, 5, 7)


def test_J4_pair_split():
    B = prove_B()
    assert B["proved"] is True
    assert n_J4_pair_special(7) == 735
    assert n_J4_pair_ordinary(7) == 294
    assert n_J4_from_pairs(7) == n_J4_p7(7) == 2646
    assert n_J4_all_pairs_equal_wrong(7) == 1764
    assert n_J4_all_pairs_equal_wrong(7) != 2646


def test_15u_dies_at_p5():
    C = prove_C()
    assert C["proved"] is True
    assert n_rest_p7(7) == 4410
    assert n_rest_p7(7) == HPLUS[7] - n_1d(7) - n_CRRR_named(7)
    assert n_rest_p7(5) == 1500
    assert n_rest_p7(5) != HPLUS[5] - n_1d(5)
    assert n_1d(5) + n_CRRR_named(5) + n_rest_p7(5) != HPLUS[5]
    assert n_1d(5) + n_NL_false_19u(5) != HPLUS[5]
    assert n_1d(7) + n_NL_false_19u(7) == HPLUS[7]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["Nr7"] == HPLUS[7]
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.456"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["affJ2_vs_J6"] is True
    assert out["proved"]["J4_pair_split"] is True
    assert out["proved"]["rest_15u_dies_at_p5"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
