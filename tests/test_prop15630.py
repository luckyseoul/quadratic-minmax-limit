from src.e1_gmin_m4_prop15630 import (
    PARI_DUAL_MIN_CERT,
    balanced_square_sum,
    balancing_gap,
    balancing_gap_closed,
    direction_count,
    dual_kissing_number,
    dual_minimum,
    dual_minimum_theorem,
    pari_dual_minimum_certificate,
    zero_sum_energy_floor,
    zero_sum_mass_floor,
)


def test_balancing_formula_and_nonzero_gap():
    for p in (3, 5, 7, 11, 13, 17, 19):
        assert balanced_square_sum(p, 0) == 0
        for t in range(-3 * p, 3 * p + 1):
            assert balancing_gap(p, t) == balancing_gap_closed(p, t)
            if t:
                assert balancing_gap(p, t) >= p
        assert balancing_gap(p, 1) == p
        assert balancing_gap(p, p) == p


def test_zero_sum_mds_newton_mass_floor():
    for p in (3, 5, 7, 11, 13, 17, 19):
        r = direction_count(p)
        floors = [zero_sum_mass_floor(p, h) for h in range(1, r + 1)]
        assert min(floors) >= r - 1
        assert zero_sum_energy_floor(p) == 2 * (r - 1) == p - 1


def test_dual_minimum_shell_theorem_and_pari_audits():
    assert dual_minimum_theorem()["proved"] is True
    assert pari_dual_minimum_certificate()["certified"] is True
    for p, row in PARI_DUAL_MIN_CERT.items():
        assert dual_minimum(p) == (1, 2)
        assert row["scaled_min"] == p
        assert row["count"] == dual_kissing_number(p) == 2 * (p * p + 1)


def test_fail_when_wrong_perturbations():
    for p in (5, 7, 11):
        assert balancing_gap(p, 2) != p
        assert dual_kissing_number(p) != 2 * p * p
        assert zero_sum_energy_floor(p) != p
