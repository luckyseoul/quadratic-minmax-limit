from fractions import Fraction

from e1_gmin_m4_prop15755 import (
    affine_odd_parameter_ledger,
    affine_bad_sign_certificate,
    dangerous_spike_dichotomy,
    equality_sign_classification,
    hereditary_cut_redundancy,
    positive_triangle_endpoint,
    theorem_record,
)


def test_affine_bad_sign_family_exact():
    for p in (3, 5, 7, 11, 13):
        row = affine_bad_sign_certificate(p)
        assert row["proved"]
        assert row["delta_plus_x"] == 2 * p
        assert row["C_v_equals_p_v"]


def test_equality_sign_cases_are_exhaustive():
    rows = [
        equality_sign_classification(x_i, a)
        for x_i in (-1, 1)
        for a in (-1, 1)
    ]
    assert all(row["proved"] for row in rows)
    assert sum(row["one_bit_boolean_eigenvector_case"] for row in rows) == 2
    assert sum(row["one_coordinate_three_case"] for row in rows) == 2


def test_odd_shell_dichotomy_and_sharp_triangle_endpoint():
    for p in (11, 13, 17, 19):
        gap = dangerous_spike_dichotomy(p)
        triangle = positive_triangle_endpoint(p)
        assert gap["proved"] and triangle["proved"]
        assert triangle["boolean_defect_after_triangle_flip"] == 6 * p - 12


def test_shared_maximizer_cut_interval_is_redundant():
    for crosses in (False, True):
        for W_B in range(21):
            assert hereditary_cut_redundancy(crosses, W_B, 20)["proved"]


def test_affine_alias_summed_cut_bound():
    for p in (5, 7, 11, 13, 17):
        row = affine_odd_parameter_ledger(p, 1)
        assert row["proved"]
        assert row["H_edge_lower_bound"] == str(Fraction((p + 1) * (p + 3), 8))


def test_theorem_record_keeps_global_gate_open():
    out = theorem_record((3, 5, 7, 11, 13))
    assert out["proved"]["dangerous_defect_dichotomy_p_ge_11"]
    assert not out["proved"]["residual_ii_closed"]
    assert out["L_status"] == "OPEN"
