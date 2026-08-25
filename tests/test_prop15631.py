from fractions import Fraction

from src.e1_gmin_m4_prop15631 import (
    coefficient_sum_phase,
    dual_first_norm,
    dual_scaled_first_norm,
    dual_scaled_second_norm_lower,
    dual_second_norm_lower,
    harmonic_min_shell_sum,
    phased_half_harmonic_min_shell_sum,
    radial_phase_from_scaled_norm,
    radial_shadow_theorem,
    rank_of,
    scaled_dual_norm,
    two_L_volume_squared,
)


def test_norm_parity_equals_coefficient_sum_parity():
    # z^T C z is even for every symmetric zero-diagonal integral C, while
    # ||z||^2 and sum z_i agree modulo two.
    samples = [
        ([1, 0, 0], 0),
        ([1, 1, 0], 2),
        ([2, -1, 3], -12),
        ([3, 4, 5, 6], 38),
    ]
    for p in (3, 5, 7, 11):
        for z, zCz in samples:
            q = scaled_dual_norm(p, sum(a * a for a in z), zCz)
            assert q % 2 == sum(z) % 2
            assert radial_phase_from_scaled_norm(q) == coefficient_sum_phase(sum(z))


def test_dual_first_gap_and_harmonic_shell():
    for p in (3, 5, 7, 11, 13, 17, 19):
        d = rank_of(p)
        assert dual_first_norm(p) == Fraction(1, 2)
        assert dual_second_norm_lower(p) == Fraction(p - 1, p)
        assert dual_scaled_first_norm(p) == p
        assert dual_scaled_second_norm_lower(p) == 2 * (p - 1)
        assert harmonic_min_shell_sum(p) == Fraction(-2, d + 2)
        assert phased_half_harmonic_min_shell_sum(p) == Fraction(1, 8 * (d + 2))


def test_radial_shadow_theorem_and_fail_when_wrong():
    assert radial_shadow_theorem()["proved"] is True
    for p in (5, 7, 11):
        assert two_L_volume_squared(p) > 1
        assert radial_phase_from_scaled_norm(p) == -1
        assert radial_phase_from_scaled_norm(2 * (p - 1)) == 1
        assert phased_half_harmonic_min_shell_sum(p) != Fraction(1, rank_of(p) + 2)
