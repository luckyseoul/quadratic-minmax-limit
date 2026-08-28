from fractions import Fraction

from e1_gmin_m4_prop15665 import (
    broad_channel_dimensions,
    constituent_partition,
    diagonal_gram_inverse,
    dimensions,
    p11_early_shell_audit,
    projected_rank_one_norm_coefficients,
    trace_harmonic_radial_correction,
    trace_harmonic_zonal_identity,
)


def test_constituent_partition_has_exact_Z_dimension():
    for p in (5, 7, 11, 13, 17, 19):
        rows = constituent_partition(p)
        assert sum(row["dimension"] * row["count"] for row in rows) == dimensions(p)[2]
        assert sum(broad_channel_dimensions(p).values()) == dimensions(p)[2]


def test_diagonal_gram_inverse_and_projection_coefficients():
    for p in (5, 7, 11, 13):
        n = p * p + 1
        a = Fraction(p * p - 1, 4 * p * p)
        b = Fraction(1, 4 * p * p)
        inverse = diagonal_gram_inverse(p)
        alpha = inverse["I_coefficient"]
        beta = inverse["J_coefficient"]
        assert a * alpha == 1
        assert a * beta + b * alpha + n * b * beta == 0
        projected = projected_rank_one_norm_coefficients(p)
        assert projected["radius_fourth_coefficient"] == 1 - beta
        assert projected["coordinate_fourth_coefficient"] == -alpha


def test_trace_radial_correction_and_p11_shell_conservation():
    _n, d, zdim = dimensions(11)
    correction = trace_harmonic_radial_correction(11)
    assert correction == -Fraction(4 * zdim, d * (d + 4)) + Fraction(
        2 * zdim, (d + 2) * (d + 4)
    )
    rows = p11_early_shell_audit()
    assert [row["scaled_norm"] for row in rows] == [11, 20, 24, 27]
    assert [row["raw_trace_mass"] for row in rows] == [
        "0",
        "89792/11",
        "7076",
        "538752",
    ]
    assert [row["harmonic_trace"] for row in rows] == [
        "-3538/63",
        "-85888/21",
        "-63684/7",
        "-527406/7",
    ]


def test_trace_polynomial_is_one_transitive_zonal_orbit_sum():
    for p in (3, 5, 7, 11, 13, 17, 19):
        n, d, _zdim = dimensions(p)
        identity = trace_harmonic_zonal_identity(p)
        coordinate = -Fraction(4 * p * p, p * p - 1)
        assert identity["coordinate_fourth_coefficient"] == coordinate
        assert identity["radius_fourth_coefficient"] == -coordinate * Fraction(
            3, 2 * (d + 2)
        )
        assert identity["one_coordinate_theta_factor"] == coordinate * n
