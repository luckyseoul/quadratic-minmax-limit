from fractions import Fraction

import numpy as np

from src.e1_gmin_m4_prop15633 import (
    PARI_SECOND_SHELL_HALF,
    audit_half_conic_lemma,
    circle_shadow_psd_scale,
    circle_shadow_scalar_offset,
    half_conic_expected_counts,
    pair_half_count,
    pair_shadow_harmonic_coefficient,
    second_shell_half_count,
    second_shell_norm,
    second_shell_signed_count,
    second_shell_theorem,
    square_circle_half_count,
)
from src.minmax_quadratic import paley_conference_prime_power


def legendre(value: int, p: int) -> int:
    value %= p
    if value == 0:
        return 0
    return 1 if pow(value, (p - 1) // 2, p) == 1 else -1


def test_standard_signed_circle_complement_is_plus_eigenvector():
    for p in (3, 5, 7, 11):
        C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
        w = np.zeros(p * p + 1, dtype=np.int64)
        for x in range(p * p):
            w[1 + x] = legendre(x // p, p)
        assert np.array_equal(C @ w, p * w)
        assert np.count_nonzero(w == 0) == p + 1
        assert int(w @ w) == p * (p - 1)
        assert Fraction(int(w @ w), p * p) == second_shell_norm(p)


def test_point_pair_vectors_have_the_second_norm():
    for p in (5, 7, 11):
        C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
        for i, j in ((0, 1), (1, 2), (2, p + 3)):
            z = np.zeros(len(C), dtype=np.int64)
            z[i] = 1
            z[j] = -int(C[i, j])
            numerator = p * z + C @ z
            scaled_norm = int(numerator @ numerator) // (2 * p)
            assert scaled_norm == 2 * (p - 1)


def test_second_shell_counts_and_independent_pari_audits():
    expected_signed = {5: 780, 7: 2800, 11: 16104, 13: 30940}
    for p, signed in expected_signed.items():
        assert pair_half_count(p) == p * p * (p * p + 1) // 2
        assert square_circle_half_count(p) == p * (p * p + 1) // 2
        assert second_shell_signed_count(p) == signed
        assert second_shell_half_count(p) == signed // 2
    for p, half_count in PARI_SECOND_SHELL_HALF.items():
        assert second_shell_half_count(p) == half_count
    assert second_shell_signed_count(3) == 30


def test_half_conic_rigidity_counts_and_exhaustive_audits():
    for p in (5, 7, 11):
        expected = half_conic_expected_counts(p)
        assert expected == {
            "rank_one_squares": (p * p - 1) // 4,
            "anisotropic_norm_multiples": (p - 1) // 2,
            "total": (p - 1) * (p + 3) // 4,
        }
        for direction_type in (1, -1):
            audit = audit_half_conic_lemma(p, direction_type)
            assert audit["checks"] is True
            assert audit["observed"]["split_nondegenerate"] == 0


def test_exact_second_shell_harmonic_coefficients():
    # These two pair values independently match the exact-looking CUDA
    # shell sums obtained from arbitrary admissible W channels.
    assert pair_shadow_harmonic_coefficient(5) == Fraction(-1, 60)
    assert pair_shadow_harmonic_coefficient(7) == Fraction(-1, 12)
    assert pair_shadow_harmonic_coefficient(11) == Fraction(-37, 252)
    assert circle_shadow_scalar_offset(11) == Fraction(-25, 693)
    assert circle_shadow_psd_scale(11) == Fraction(1, 8 * 11**4)


def test_theorem_is_exact_but_keeps_r1_open():
    theorem = second_shell_theorem()
    assert theorem["proved"] is True
    assert theorem["p3_exception"]["candidate_orbits_coincide"] is True
    assert "point-pair" in theorem["classification"]
