import itertools
from fractions import Fraction

import numpy as np
import pytest

from src.e1_gmin_m4_prop15634 import explicit_square_circles
from src.e1_gmin_m4_prop15640 import (
    complete_quartic_scalar,
    harmonic_circle_high_closed,
    harmonic_circle_low_closed,
    harmonic_kernel_closed,
    harmonic_operator_theorem,
    harmonic_scalar_offset,
    harmonic_sign_certificate,
    harmonic_spectrum,
    negative_triangle_quartic_scalar,
    off_circle_row_contraction_scalar,
    parity_twisted_half_spectrum,
    point_circle_evaluation_coefficient,
    point_circle_quartic_scalar,
    radial_harmonic_correction,
    radial_harmonic_correction_closed,
    rank_of,
    shell_signed_count,
    shell_squared_norm,
    through_point_circle_frame_certificate,
    total_circle_frame_scalar,
    z_dimension,
)
from src.minmax_quadratic import paley_conference_prime_power


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31])
def test_closed_quartic_and_harmonic_formulas(p):
    assert shell_squared_norm(p).numerator == 3 * p - 6
    assert shell_squared_norm(p).denominator == 2 * p
    assert negative_triangle_quartic_scalar(p) == 2 * (p - 3) * (p + 1)
    assert point_circle_quartic_scalar(p) == 8 * (p - 2)
    assert point_circle_evaluation_coefficient(p).numerator == 2 * (p - 5)
    assert point_circle_evaluation_coefficient(p).denominator == p**3
    assert complete_quartic_scalar(p) == 2 * (p * p + 2 * p - 11)
    assert radial_harmonic_correction(p) == radial_harmonic_correction_closed(p)
    assert harmonic_scalar_offset(p) == harmonic_kernel_closed(p)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23])
def test_through_point_frame_and_off_circle_contraction(p):
    cert = through_point_circle_frame_certificate(p)
    assert cert["proved"]
    assert cert["parallel_classes"] == (p + 1) // 2
    assert cert["within_class_nonzero_eigenvalue"] == p * p
    assert cert["span_dimension"] == rank_of(p) - 1
    assert total_circle_frame_scalar(p) == p * p * (p - 1)
    assert off_circle_row_contraction_scalar(p) == p * p * (p - 2)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29])
def test_complete_harmonic_spectrum_is_a_saddle(p):
    rows = harmonic_spectrum(p)
    assert [row["channel"] for row in rows] == [
        "circle-kernel",
        "circle-low",
        "circle-high",
    ]
    assert [row["eigenvalue"] for row in rows] == [
        harmonic_kernel_closed(p),
        harmonic_circle_low_closed(p),
        harmonic_circle_high_closed(p),
    ]
    assert [np.sign(row["eigenvalue"]) for row in rows] == [-1, 1, 1]
    assert sum(row["multiplicity"] for row in rows) == z_dimension(p)
    signs = harmonic_sign_certificate(p)
    assert signs["expansions_match"]
    assert signs["indefinite"]


def test_p11_exact_spectrum_and_multiplicities():
    rows = harmonic_spectrum(11)
    assert [row["eigenvalue"] for row in rows] == [
        -Fraction(582, 7),
        Fraction(258, 7),
        Fraction(426, 7),
    ]
    assert [row["multiplicity"] for row in rows] == [1220, 305, 244]
    assert shell_signed_count(11) == 442_860


@pytest.mark.parametrize("p", [11, 13, 17, 19])
def test_norm_parity_twist_reverses_every_spectrum_sign(p):
    raw = harmonic_spectrum(p)
    shadow = parity_twisted_half_spectrum(p)
    assert [row["eigenvalue"] for row in shadow] == [
        -row["eigenvalue"] / 16 for row in raw
    ]
    assert [row["closed_form"] for row in shadow] == [
        -row["closed_form"] / 16 for row in raw
    ]
    assert [row["unphased_closed_form"] for row in shadow] == [
        row["closed_form"] for row in raw
    ]
    assert [np.sign(row["eigenvalue"]) for row in shadow] == [1, -1, -1]


@pytest.mark.parametrize("p", [3, 5, 7, 11])
def test_explicit_through_point_circle_frame_identity(p):
    blocks, words = explicit_square_circles(p)
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    n = len(C)
    twice_p_times_P = p * np.eye(n, dtype=np.int64) + C
    for i in range(n):
        selected = [index for index, block in enumerate(blocks) if i in block]
        frame = words[selected].T @ words[selected]
        column = twice_p_times_P[:, i]
        twice_target = p * twice_p_times_P - np.outer(column, column)
        assert np.array_equal(2 * frame, twice_target)


def _random_admissible_matrix(C, p, seed=640):
    values, vectors = np.linalg.eigh(C.astype(np.float64))
    U = vectors[:, np.isclose(values, p, atol=1e-8)]
    d = U.shape[1]
    intrinsic_basis = []
    for a in range(d):
        matrix = np.zeros((d, d))
        matrix[a, a] = 1
        intrinsic_basis.append(matrix)
        for b in range(a + 1, d):
            matrix = np.zeros((d, d))
            matrix[a, b] = matrix[b, a] = 1 / np.sqrt(2)
            intrinsic_basis.append(matrix)
    intrinsic_basis = np.asarray(intrinsic_basis)
    diagonal_map = np.einsum(
        "ia,kab,ib->ik", U, intrinsic_basis, U, optimize=True
    )
    _left, singular, vh = np.linalg.svd(diagonal_map, full_matrices=True)
    tolerance = max(diagonal_map.shape) * singular[0] * np.finfo(float).eps
    rank = int(np.count_nonzero(singular > tolerance))
    rng = np.random.default_rng(seed)
    coefficients = rng.standard_normal(len(intrinsic_basis) - rank) @ vh[rank:]
    intrinsic = np.einsum("k,kab->ab", coefficients, intrinsic_basis)
    W = U @ intrinsic @ U.T
    W /= np.linalg.norm(W)
    assert np.max(np.abs(np.diag(W))) < 2e-12
    assert np.max(np.abs(((np.eye(len(C)) + C / p) / 2) @ W - W)) < 2e-12
    return W


@pytest.mark.parametrize("p", [5, 7, 11])
def test_direct_family_quartic_sums_match_closed_decomposition(p):
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    W = _random_admissible_matrix(C, p)
    F = float(np.sum(W * W))

    triangle_sum = 0.0
    for i, j, k in itertools.combinations(range(len(C)), 3):
        if C[i, j] * C[i, k] * C[j, k] != -1:
            continue
        value = C[i, j] * W[i, j] + C[i, k] * W[i, k] + C[j, k] * W[j, k]
        triangle_sum += 8 * value * value
    assert np.isclose(
        triangle_sum,
        negative_triangle_quartic_scalar(p) * F,
        atol=5e-10,
        rtol=5e-10,
    )

    blocks, words = explicit_square_circles(p)
    point_circle_sum = 0.0
    circle_evaluation_sum = 0.0
    for block, word in zip(blocks, words, strict=True):
        q = float(word @ W @ word)
        circle_evaluation_sum += q * q
        Ww = W @ word
        for i in range(len(C)):
            if i in block:
                continue
            value = q / (p * p) - 2 * word[i] * Ww[i] / p
            point_circle_sum += 2 * value * value
    predicted = (
        point_circle_quartic_scalar(p) * F
        + float(point_circle_evaluation_coefficient(p)) * circle_evaluation_sum
    )
    assert np.isclose(point_circle_sum, predicted, atol=5e-10, rtol=5e-10)


def test_theorem_keeps_r1_and_l_open():
    theorem = harmonic_operator_theorem()
    assert theorem["proved"]
    assert "scaled-norm 3p-6" in theorem["scope"]
    assert all(row["checks"] for row in theorem["rows"].values())
