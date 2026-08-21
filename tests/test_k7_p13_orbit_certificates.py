"""Consistency checks for the exact p=13 signed-orbit certificates."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "evidence"))

from e1_gmin_m4_prop15588 import _fit_poly_modp  # noqa: E402
from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k5_p29_coefficient_sieve import square_directions  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def _load(name: str) -> dict:
    return json.loads((ROOT / "evidence" / name).read_text())


def test_first_signed_orbit_is_free_and_balanced():
    report = _load("k7_p13_signed_psl_orbit.json")
    assert report["complete_orbit"] is True
    assert report["orbit_size"] == report["signed_psl_lift_order"] == 4_826_640
    assert report["projective_psl_order"] == 2_413_320
    assert report["signed_stabilizer_order"] == 1
    assert report["epsilon_plus_count"] == report["epsilon_minus_count"] == 2_413_320
    assert report["packed_sha256"] == (
        "7223169420a18477dbdf95f6c3685186fbc6a7a1916ac875d761b22800c01eb2"
    )


def test_orbit_and_k7_slice_clear_qvar_exactly():
    report = _load("k7_p13_orbit_quartic_xpu.json")
    threshold = Fraction(report["QVAR_threshold"])
    histogram = {int(k): int(v) for k, v in report["histogram"].items()}
    assert sum(histogram.values()) == report["evaluated_rows"] == 2_413_320
    assert sum(value * count for value, count in histogram.items()) == report[
        "sum_abs_Zpsi_sq"
    ]
    assert Fraction(report["E_abs_Zpsi_sq"]) == Fraction(806_468, 85) > threshold
    assert Fraction(report["induced_exceptional_scalar"]) == Fraction(19_088, 1_785) > 6

    by_k = report["by_activity"]
    assert {int(k): row["count"] for k, row in by_k.items()} == {
        5: 14_196,
        6: 28_392,
        7: 2_370_732,
    }
    assert Fraction(by_k["7"]["E_abs_Zpsi_sq"]) == Fraction(1_606_124, 167)
    assert Fraction(by_k["7"]["E_abs_Zpsi_sq"]) > threshold
    assert by_k["7"]["clears_QVAR"] is True


def test_nonzero_quintic_scalars_are_equidistributed_in_first_orbit():
    report = _load("k7_p13_extract_orbit_representatives.json")
    scalar_counts = {int(k): int(v) for k, v in report["top_scalar_histogram"].items()}
    assert scalar_counts == {scalar: 1_071 for scalar in range(1, 13)}
    assert report["depressed_k7_representatives"] == 12 * 1_071 == 12_852
    assert report["scalar7_representatives"] == 1_071
    assert report["representatives_sha256"] == (
        "341cbfec0e076e36df4f7b517d98b51689e5bd8d3bdc82478a12d92c1ba992eb"
    )


def test_second_orbit_seed_is_an_independent_k7_maxplus_vector():
    report = _load("k7_p13_second_orbit_seed.json")
    p = 13
    q = p * p
    negative_indices = np.asarray(report["negative_indices"], dtype=np.int64)
    y = np.ones(q + 1, dtype=np.int64)
    y[1 + negative_indices] = -1
    conference = paley_conference_prime_power(p).astype(np.int64)
    assert np.array_equal(conference @ y, p * y)

    coefficients = []
    for coordinate, _form in square_directions(p):
        line_sums = [int(y[1:][coordinate == value].sum()) for value in range(p)]
        rho = ((np.asarray(line_sums) + p - 2) // 2) % p
        coefficients.append(_fit_poly_modp(rho, p))
    assert all(max(i for i, value in enumerate(row) if value) == 5 for row in coefficients)
    assert all(row[4] == 0 for row in coefficients)

    kernel_real, kernel_imag = quartic_kernel(p)
    negative = (y[1:] < 0).astype(np.int64)
    real = int(negative @ kernel_real.astype(np.int64) @ negative)
    imag = int(negative @ kernel_imag.astype(np.int64) @ negative)
    assert (real, imag) == (-132, -198)
    assert real * real + imag * imag == 56_628


def test_second_signed_orbit_is_free_and_balanced():
    report = _load("k7_p13_second_signed_psl_orbit.json")
    assert report["complete_orbit"] is True
    assert report["orbit_size"] == report["signed_psl_lift_order"] == 4_826_640
    assert report["projective_psl_order"] == 2_413_320
    assert report["signed_stabilizer_order"] == 1
    assert report["epsilon_plus_count"] == report["epsilon_minus_count"] == 2_413_320
    assert report["packed_sha256"] == (
        "a3ce4e19e68770b41951b4ba28153fd5ed23884d0bcd912eeb43c421fa0e31c3"
    )


def test_second_orbit_and_k7_slice_clear_qvar_exactly():
    report = _load("k7_p13_second_orbit_quartic_xpu.json")
    threshold = Fraction(report["QVAR_threshold"])
    histogram = {int(k): int(v) for k, v in report["histogram"].items()}
    assert sum(histogram.values()) == report["evaluated_rows"] == 2_413_320
    assert sum(value * count for value, count in histogram.items()) == report[
        "sum_abs_Zpsi_sq"
    ]
    assert Fraction(report["E_abs_Zpsi_sq"]) == Fraction(806_468, 85) > threshold
    assert Fraction(report["induced_exceptional_scalar"]) == Fraction(19_088, 1_785) > 6

    by_k = report["by_activity"]
    assert {int(k): row["count"] for k, row in by_k.items()} == {
        6: 28_392,
        7: 2_384_928,
    }
    assert Fraction(by_k["7"]["E_abs_Zpsi_sq"]) == Fraction(198_692, 21)
    assert Fraction(by_k["7"]["E_abs_Zpsi_sq"]) > threshold
    assert by_k["7"]["clears_QVAR"] is True


def test_nonzero_quintic_scalars_are_equidistributed_in_second_orbit():
    report = _load("k7_p13_extract_second_orbit_representatives.json")
    scalar_counts = {int(k): int(v) for k, v in report["top_scalar_histogram"].items()}
    assert scalar_counts == {scalar: 1_071 for scalar in range(1, 13)}
    assert report["scalar7_representatives"] == 1_071
    assert report["representatives_sha256"] == (
        "49cacd5fe098de3aa1d2a88ecd94dba52debc532861026c3b798c429a553dbab"
    )


def test_first_and_second_scalar7_representatives_are_disjoint():
    report = _load("k7_p13_scalar7_union_orbit12.json")
    assert report["orbit1_count"] == report["orbit2_count"] == 1_071
    assert report["overlap"] == 0
    assert report["union_count"] == 2_142
    assert report["disjoint"] is True
    assert report["union_sha256"] == (
        "ac2615f788195cc0f446f5573d82040a72a2f28245207093cee561557fbaecf2"
    )


def test_independent_seeds_found_the_same_second_orbit_vector():
    seed17 = _load("k7_p13_orbit_completeness_seed17.json")
    seed41 = _load("k7_p13_orbit_completeness_seed41.json")
    assert seed17["outside_orbit_solution_found"] is True
    assert seed41["outside_orbit_solution_found"] is True
    assert seed17["outside_orbit_Zpsi"]["abs_sq"] == seed41["outside_orbit_Zpsi"]["abs_sq"] == 56_628
    assert seed17["outside_orbit_Zpsi"]["real"] == seed41["outside_orbit_Zpsi"]["real"] == -132
    assert seed17["outside_orbit_Zpsi"]["imag"] == seed41["outside_orbit_Zpsi"]["imag"] == -198
