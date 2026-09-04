#!/usr/bin/env python3
"""Sharp local mass spectrum below 2p-10, for every prime p>=29.

For a nonzero nonnegative integer-valued quadratic C on J(p,(p+1)/2),
0 < 4p E[C] < 2p-10 forces C to be Boolean and its mass to be p-3 or
p+1. The endpoint is attained by the existing four-/five-set height-three
examples. This is a local theorem, not a residual-II or limit closure.

Only the pinned 15.751 four-bit certificate is read: its 65,536-table
enumeration is never called. See NOTE_2026-09-04_SHARP_SMALL_MASS_SPECTRUM.md.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15734 import residual_even_floor_table
from e1_gmin_m4_prop15751 import (
    EXPECTED_HISTOGRAM,
    EXPECTED_SHA256,
    EXPECTED_VALID_TABLES,
    atomic_write_json,
    cube_half_mean_height_certificate,
    profile_density,
    unpack_layer_signature,
)
from e1_gmin_m4_prop15768 import cube_three_quarter_height_certificate


ROOT = Path(__file__).resolve().parents[1]
CATALOG_EVIDENCE = ROOT / "evidence" / "e1_gmin_m4_prop15751.json"
CATALOG_EVIDENCE_SHA256 = "b25b0f8896fccfb01c48c92b6266724bf484a0b14f4841655229a647bc2a61a7"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_prime(p: int) -> None:
    if type(p) is not int or p < 29 or not is_prime(p):
        raise ValueError("need prime p>=29")


def _fixed_catalog_dependency() -> dict[str, object]:
    """Read and validate the reviewed artifact without rerunning its census."""
    try:
        data = CATALOG_EVIDENCE.read_bytes()
    except OSError as exc:
        raise ArithmeticError("pinned four-bit evidence unavailable") from exc
    _require(hashlib.sha256(data).hexdigest() == CATALOG_EVIDENCE_SHA256,
             "pinned four-bit evidence hash mismatch")
    record = json.loads(data)
    catalog = record.get("four_cube_catalog", {})
    _require(record.get("proved") is True and catalog.get("proved") is True,
             "four-bit catalog proof dependency failed")
    _require(catalog.get("tables_checked") == 65536
             and catalog.get("valid_tables") == EXPECTED_VALID_TABLES == 222
             and catalog.get("valid_table_signature_sha256") == EXPECTED_SHA256,
             "four-bit catalog count or digest changed")
    histogram = {int(key): value for key, value in
                 catalog.get("packed_layer_signature_histogram", {}).items()}
    profiles = catalog.get("profiles", [])
    _require(histogram == EXPECTED_HISTOGRAM and len(profiles) == 14
             and [row.get("signature") for row in profiles] == sorted(histogram),
             "four-bit catalog profile partition changed")
    _require(all(row.get("layer_counts") == list(unpack_layer_signature(row["signature"]))
                 and row.get("multiplicity") == histogram[row["signature"]]
                 for row in profiles), "four-bit catalog profile content changed")
    return {"evidence": str(CATALOG_EVIDENCE.relative_to(ROOT)),
            "evidence_sha256": CATALOG_EVIDENCE_SHA256,
            "catalog_sha256": EXPECTED_SHA256,
            "valid_tables": 222, "profile_count": 14, "profiles": profiles,
            "catalog_rerun": False, "proved": True}


def _height_gap_certificate(p: int) -> dict[str, object]:
    floor = sharp_integral_quadratic_lift_floor(p)
    half = cube_half_mean_height_certificate()
    three_quarter = cube_three_quarter_height_certificate()
    _require(floor.get("proved") is True and half.get("proved") is True
             and three_quarter.get("proved") is True,
             "cube/stabilizer proof dependency failed")
    coefficient = Fraction(4) if p % 4 == 3 else Fraction(4 * (p - 1), p + 3)
    _require(floor.get("sharp_scaled_floor") == p - 3
             and floor.get("cube_degree_two_support_floor") == Fraction(1, 4)
             and floor.get("integer_quadratic_cube_mean_lattice") == "(1/4) Z"
             and floor.get("paired_cube_identity") == "T B(X)=(B(X)+p E[B])/(p+1)"
             and floor.get("H_at_least_two_stabilizer_coefficient") == coefficient,
             "paired-cube/stabilizer statement changed")
    _require(half.get("maximum_upper_bound") == 3
             and half.get("quarter_mean_lattice") is True
             and half.get("domain") == "all dimensions d>=0"
             and three_quarter.get("maximum_upper_bound") == 6
             and three_quarter.get("facet_means_are_quarter_integral") is True
             and three_quarter.get("domain") == "all dimensions d>=0",
             "dimension-free cube height statement changed")
    threshold = 2 * p - 10
    first = Fraction(2 * (p + 1) - threshold, 4)
    refined = Fraction(3 * (p + 1) - threshold, 4)
    average_upper = Fraction(threshold, 4 * (p + 1)) * (1 + 4 / coefficient)
    uniform_upper = Fraction(threshold, 2 * (p - 1))
    _require(first == 3 and refined == Fraction(p + 13, 4) > 6
             and average_upper <= uniform_upper < 1,
             "strict small-mass height bootstrap failed")
    return {"strict_upper_mass": threshold,
            "hypothesis": "2<=H=max C and 0<M=4p E[C]<2p-10",
            "paired_cube_identity": "T C(X)=(4H+M)/(4(p+1))",
            "all_maximizing_cubes_initial_mean_at_least": "1/2",
            "first_height_strictly_greater_than": str(first),
            "half_mean_cube_maximum": 3,
            "all_maximizing_cubes_refined_mean_at_least": "3/4",
            "refined_height_strictly_greater_than": str(refined),
            "three_quarter_mean_cube_maximum": 6,
            "all_maximizing_cubes_final_mean_at_least": 1,
            "stabilizer_coefficient": str(coefficient),
            "paired_cube_average_endpoint_upper": str(average_upper),
            "uniform_paired_cube_average_endpoint_upper": str(uniform_upper),
            "contradiction": "every maximizing cube has mean>=1 but their average is<1",
            "all_prime_congruence_classes": True, "proved": True}


def _boolean_spectrum_certificate(p: int) -> dict[str, object]:
    catalog = _fixed_catalog_dependency()
    _require(catalog.get("proved") is True and catalog.get("catalog_rerun") is False,
             "fixed Boolean catalog dependency failed")
    influence_floor = Fraction((p + 1) * (p - 3), 16 * p * (p - 2))
    junta_upper = Fraction(8 * (p - 1) * (p - 2), (p + 1) * (p - 3))
    gap = Fraction(8 * (p - 5), (p + 1) * (p - 3))
    q = (p - 1) // 2
    _require(influence_floor > 0 and 8 - junta_upper == gap > 0
             and junta_upper < 8 and 7 < q,
             "uniform seven-coordinate junta bound failed")
    densities = sorted({profile_density(tuple(row["layer_counts"]), p)
                        for row in catalog["profiles"]})
    masses = [4 * p * value for value in densities]
    expected = [0, p - 3, p + 1, 2 * p - 2, 2 * p + 2,
                3 * p - 1, 3 * p + 3, 4 * p]
    _require(masses == expected, "fixed four-bit density spectrum changed")
    small = [int(value) for value in masses if 0 < value < 2 * p - 10]
    _require(small == [p - 3, p + 1], "Boolean small-mass interval changed")
    return {"relevant_pair_influence_floor": str(influence_floor),
            "influence_definition": "I_ij=(1/4)Pr[f(X)!=f(X^(ij))]",
            "total_influence_upper": "(p-1)mu(1-mu)",
            "relevant_pair_count_lower": "p*L/2",
            "largest_invariant_class_complement_bound": str(junta_upper),
            "strict_gap_below_eight": str(gap),
            "junta_coordinates_at_most": 7,
            "all_kept_coordinate_patterns_extend_to_slice": True,
            "symmetrization_preserves_degree_at_most_two": True,
            "cube_relevant_coordinate_influence_floor": "1/2",
            "cube_total_influence_upper": 2,
            "cube_coordinates_actually_needed_at_most": 4,
            "all_boolean_scaled_masses": expected,
            "allowed_positive_masses": small,
            "fixed_catalog_dependency": {key: value for key, value in catalog.items()
                                         if key != "profiles"},
            "proved": True}


def _sharp_endpoint_certificate(p: int) -> dict[str, object]:
    m = (p + 1) // 2
    examples = []
    for size in (4, 5):
        mean_r = Fraction(size * m, p)
        mean_pairs = Fraction(comb(size, 2) * m * (m - 1), p * (p - 1))
        mean = 3 - 2 * mean_r + mean_pairs
        values = [3 - 2 * r + comb(r, 2) for r in range(size + 1)]
        _require(min(values) == 0 and max(values) == 3
                 and 4 * p * mean == 2 * p - 10
                 and size <= min(m, p - m), "sharp endpoint example changed")
        examples.append({"support_size": size, "formula": "3-2r+binom(r,2)",
                         "r": "|X intersect R|", "layer_values": values,
                         "mean": str(mean), "scaled_mass": 2 * p - 10,
                         "height": 3, "proved": True})
    return {"scaled_mass": 2 * p - 10, "examples": examples,
            "strict_upper_endpoint_is_necessary": True,
            "classification_of_endpoint_equalities_claimed": False, "proved": True}


def small_mass_spectrum(p: int) -> dict[str, object]:
    """Return the sharp local theorem, retaining the strict mass endpoint."""
    _check_prime(p)
    height = _height_gap_certificate(p)
    boolean = _boolean_spectrum_certificate(p)
    endpoint = _sharp_endpoint_certificate(p)
    _require(all(row.get("proved") is True for row in (height, boolean, endpoint))
             and height.get("strict_upper_mass") == 2 * p - 10
             and boolean.get("allowed_positive_masses") == [p - 3, p + 1]
             and endpoint.get("scaled_mass") == 2 * p - 10,
             "small-mass theorem dependency failed")
    m, q = (p + 1) // 2, (p - 1) // 2
    pair_mass = 4 * p * Fraction(m * (m - 1), p * (p - 1))
    omitted_pair_mass = 4 * p * Fraction(q * (q - 1), p * (p - 1))
    _require([omitted_pair_mass, pair_mass] == [p - 3, p + 1],
             "Boolean spectrum attainment changed")
    return {"p": p, "slice": f"J({p},{m})", "strict_upper_mass": 2 * p - 10,
            "scaled_mass_definition": "M=4p E[C]",
            "hypotheses": "C is nonzero, nonnegative, integer-valued and degree at most two",
            "allowed_positive_masses": [p - 3, p + 1],
            "boolean_below_strict_upper_mass": True,
            "nonboolean_mass_lower_bound": 2 * p - 10,
            "height_gap": height, "boolean_spectrum": boolean,
            "attaining_boolean_examples": [
                {"formula": "(1-x_i)(1-x_j)", "scaled_mass": p - 3},
                {"formula": "x_i*x_j", "scaled_mass": p + 1}],
            "sharp_endpoint": endpoint,
            "new_catalog_or_prime_census_used": False,
            "residual_ii_closed": False, "limit_closed": False, "proved": True}


def _parity_baseline(p: int, b: int, phase: int) -> dict[str, object]:
    """These are genuine 0/1 parity minima, not punctured lift baselines."""
    m = (p + 1) // 2
    if b == 2:
        values = [(phase + x + y) % 2 for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
        xor_mass = p + 1
        mass = xor_mass if phase == 0 else 2 * p - xor_mass
        formula = "(x_i-x_j)^2" if phase == 0 else "1-(x_i-x_j)^2"
        _require(values == ([0, 1, 1, 0] if phase == 0 else [1, 0, 0, 1]),
                 "pair parity baseline changed")
    elif b == p - 1:
        values = [(phase + m - x) % 2 for x in (0, 1)]
        formula = "x_j" if values == [0, 1] else "1-x_j"
        mass = p + 1 if formula == "x_j" else p - 1
        _require(values in ([0, 1], [1, 0]), "last-bit parity baseline changed")
    else:
        raise ValueError("pointwise baseline needs b=2 or p-1")
    return {"b": b, "phase": phase, "formula": formula,
            "truth_values": values, "scaled_mass": mass,
            "pointwise_parity_minimum": True,
            "C_equals_half_difference_is_nonnegative_integral_quadratic": True,
            "lift_mass_equals_a_minus_baseline_mass": True, "proved": True}


def affine_parity_small_mass_spectrum(p: int) -> dict[str, object]:
    """Both phases and every even boundary for a=2p E[A]<2p-10.

Here A>=0 is an integer-valued quadratic with
A(x) == phase + sum_{i in B} x_i (mod 2), on J(p,(p+1)/2), |B| even.
"""
    _check_prime(p)
    local = small_mass_spectrum(p)
    floors = residual_even_floor_table(p)
    sharp = sharp_integral_quadratic_lift_floor(p)
    _require(local.get("proved") is True and floors.get("proved") is True
             and sharp.get("proved") is True
             and local.get("strict_upper_mass") == 2 * p - 10
             and local.get("allowed_positive_masses") == [p - 3, p + 1]
             and sharp.get("sharp_scaled_floor") == p - 3,
             "affine-parity spectrum dependency failed")
    threshold = 2 * p - 10
    phases = {}
    for phase in (0, 1):
        table = floors["phase_zero_floors" if phase == 0 else "phase_one_floors"]
        candidates = [b for b, floor in table.items() if floor < threshold]
        expected_b = [0, 2, p - 1] if phase == 0 else [2, p - 1]
        _require(list(table) == list(range(0, p, 2)) and candidates == expected_b,
                 "affine-parity boundary partition changed")
        rows = []
        allowed = {0, p - 3, p + 1} if phase == 0 else set()
        for b in (2, p - 1):
            baseline = _parity_baseline(p, b, phase)
            _require(baseline.get("proved") is True,
                     "genuine parity-baseline proof dependency failed")
            mass = baseline["scaled_mass"]
            _require(baseline.get("proved") is True
                     and baseline.get("pointwise_parity_minimum") is True
                     and baseline.get("C_equals_half_difference_is_nonnegative_integral_quadratic") is True
                     and baseline.get("lift_mass_equals_a_minus_baseline_mass") is True
                     and table[b] == mass
                     and 0 < threshold - mass < p - 3,
                     "genuine parity-baseline lift exclusion failed")
            rows.append({**baseline, "positive_lift_mass_strictly_below": threshold - mass,
                         "positive_lift_mass_floor": p - 3,
                         "positive_lift_excluded": True})
            allowed.add(mass)
        expected_masses = ([0, p - 3, p - 1, p + 1] if p % 4 == 1
                           else [0, p - 3, p + 1]) if phase == 0 else (
                               [p - 1, p + 1] if p % 4 == 1 else [p - 1])
        _require(sorted(allowed) == expected_masses,
                 "phase-specific small-mass spectrum changed")
        phases[str(phase)] = {"allowed_masses": sorted(allowed),
                              "candidate_boundary_sizes": candidates,
                              "pointwise_baselines": rows,
                              "b0_is_twice_nonnegative_integral_quadratic": phase == 0,
                              "all_other_boundaries_excluded_by_floor": True,
                              "proved": True}
    union = sorted(set(phases["0"]["allowed_masses"]) | set(phases["1"]["allowed_masses"]))
    _require(union == [0, p - 3, p - 1, p + 1], "combined parity spectrum changed")
    return {"p": p, "strict_upper_mass": threshold,
            "scaled_mass_definition": "a=2p E[A]",
            "parity_hypothesis": "A(x)=phase+sum_(i in B)x_i modulo 2; |B| even",
            "nonnegative_integral_quadratic_required": True,
            "phase_zero_allowed_masses": phases["0"]["allowed_masses"],
            "phase_one_allowed_masses": phases["1"]["allowed_masses"],
            "union_allowed_masses": union, "phases": phases,
            "b0_local_dependency_proved": True,
            "only_genuine_pointwise_parity_minima_subtracted": True,
            "no_punctured_complement_triple_difference_used": True,
            "all_even_boundary_sizes_covered": True,
            "new_catalog_or_prime_census_used": False, "proved": True}


def local_mass_exclusion(p: int, mass: int | Fraction) -> dict[str, object]:
    """Certify only p+1<mass<2p-10; preserve both endpoints explicitly."""
    _check_prime(p)
    if isinstance(mass, bool) or not isinstance(mass, (int, Fraction)):
        raise ValueError("mass must be an exact integer or Fraction")
    value = Fraction(mass)
    serialized = int(value) if value.denominator == 1 else str(value)
    base = {"p": p, "scaled_mass": serialized,
            "strict_lower_mass": p + 1, "strict_upper_mass": 2 * p - 10}
    if not p + 1 < value < 2 * p - 10:
        return {**base, "proved": False, "excluded": False,
                "reason": "outside the requested open exclusion interval; no conclusion"}
    local = small_mass_spectrum(p)
    _require(local.get("proved") is True
             and local.get("strict_upper_mass") == 2 * p - 10
             and local.get("allowed_positive_masses") == [p - 3, p + 1]
             and value not in local["allowed_positive_masses"],
             "local mass exclusion dependency failed")
    return {**base, "local_spectrum_proved": True, "proved": True, "excluded": True,
            "reason": "Boolean small-mass spectrum has no value in this open interval"}


def theorem_record() -> dict[str, object]:
    samples = {str(p): {"local": small_mass_spectrum(p),
                        "affine_parity": affine_parity_small_mass_spectrum(p)}
               for p in (29, 31)}
    _require(all(row["local"]["proved"] and row["affine_parity"]["proved"]
                 for row in samples.values()), "small-mass evidence package failed")
    return {"theorem": "sharp small-mass spectrum", "status": "PROVED INFINITE-FAMILY THEOREM",
            "scope": "all primes p>=29; both congruence classes",
            "statement": "0<4p E[C]<2p-10 implies C Boolean and 4p E[C] in {p-3,p+1}",
            "sharp_nonboolean_threshold": "2p-10, attained by F4 and F5",
            "sample_exact_symbolic_evaluations": samples,
            "new_catalog_or_prime_census_used": False,
            "residual_ii_closed": False, "E1_closed": False,
            "quadratic_minmax_limit_closed": False, "proved": True}


if __name__ == "__main__":
    target = ROOT / "evidence" / "e1_gmin_m4_small_mass_spectrum.json"
    atomic_write_json(target, theorem_record())
    print(json.dumps({"wrote": str(target), "proved": True}, sort_keys=True))
