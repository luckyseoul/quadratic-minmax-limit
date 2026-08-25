from src.e1_gmin_m4_prop15645 import (
    additive_l1,
    classify_under_budget,
    positive_mass_lower_bound,
    theorem_baseline_fibre_profile,
)


def test_exact_ideal_and_transfer_costs():
    for p in range(7, 102, 2):
        ideal = (0,) * p
        transfer = (1, -1) + (0,) * (p - 2)
        assert additive_l1(ideal) == 0
        assert additive_l1(transfer) == 2 * p - 4
        assert classify_under_budget(p, ideal, 2 * p) == "ideal"
        assert classify_under_budget(p, transfer, 2 * p) == "one_transfer"


def test_two_units_of_positive_mass_exceed_budget():
    for p in range(7, 102, 2):
        assert positive_mass_lower_bound(p, 2) == 4 * p - 12
        assert positive_mass_lower_bound(p, 2) > 2 * p
        double = (1, 1, -1, -1) + (0,) * (p - 4)
        assert additive_l1(double) == 4 * p - 12
        assert classify_under_budget(p, double, 2 * p) == "over_budget"


def test_theorem_keeps_live_gates_open():
    theorem = theorem_baseline_fibre_profile()
    assert theorem["proved"] is True
    assert theorem["closes_negative_product_branch"] is False
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
