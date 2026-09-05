#!/usr/bin/env python3
"""Prop. 15.775: the whole first post-15.774 layer, for p>=259201.

The proof is analytic, not a prime, graph, slice, or coefficient census.
Bounded-mean integral slice quadratics have bounded height and a bounded
invariant-class complement.  Their full-cube extensions have quarter-integral
means, excluding the two local masses 2p+4 and 2p+6 eventually.  The exact
15.774 two-type capacities force these masses at |H|=5p+6.

The same argument also excludes a superlinear, explicitly bounded support
band.  Neither conclusion is an all-size or limit theorem.
See evidence/NOTE_2026-09-04_EVENTUAL_FIRST_LAYER_CLOSE.md.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_prop15642 import hypergeometric_moments, stabilizer_mass_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15774 import quotient_floor
from e1_gmin_m4_small_mass_spectrum import affine_parity_small_mass_spectrum
from io_atomic import write_json_atomic

ROOT = Path(__file__).resolve().parents[1]
PRIME_THRESHOLD = 259201
HEIGHT_BOUND = 1800
MEAN_UPPER = Fraction(9, 8)
JUNTA_STRICT_BOUND = 129600
CUBE_HEIGHT_CONSTANT = 324
SAMPLE_PRIMES = (524287, 6700417)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _odd(p: int, minimum: int = 29) -> None:
    if type(p) is not int or p < minimum or p % 2 == 0:
        raise ValueError(f"need odd p>={minimum}")


def _shift(s: int) -> None:
    if type(s) is not int or s not in (4, 6):
        raise ValueError("only the local shifts s=4,6 are proved here")


def cube_height_bound_certificate() -> dict[str, object]:
    """Check the exact constants in the self-contained noise proof."""
    rho_squared = Fraction(1, 3)
    lhs = [Fraction(1), 6 * rho_squared, rho_squared**2]
    rhs = [Fraction(1), Fraction(2), Fraction(1)]
    remainder = [b - a for a, b in zip(lhs, rhs)]
    inverse_noise = 1 / rho_squared
    l2_from_l1 = inverse_noise**2
    fourier_denominator = 4
    height_constant = fourier_denominator * l2_from_l1**2
    _require(lhs == [1, 2, Fraction(1, 9)]
             and remainder == [0, 0, Fraction(8, 9)]
             and all(value >= 0 for value in remainder)
             and inverse_noise == 3 and l2_from_l1 == 9
             and height_constant == CUBE_HEIGHT_CONSTANT == 324,
             "noise/Fourier height constants changed")
    return {
        "domain": "nonnegative integer-valued cube quadratics, every dimension",
        "rho_squared": str(rho_squared),
        "two_point_coefficients_a4_a2b2_b4": list(map(str, lhs)),
        "two_point_nonnegative_remainder": list(map(str, remainder)),
        "tensorization_uses_Cauchy_Schwarz": True,
        "degree_two_L4_over_L2": str(inverse_noise),
        "interpolation_exponents": ["1/3", "2/3"],
        "L2_over_L1": str(l2_from_l1),
        "Fourier_coefficient_lattice": "(1/4) Z",
        "height_over_mean_squared": int(height_constant),
        "proof_is_dimension_free": True, "proved": True,
    }


def paired_cube_basis_certificate(p: int) -> dict[str, object]:
    """Check T f(X)=(f(X)+p E[f])/(p+1) on every monomial orbit."""
    _odd(p)
    m = (p + 1) // 2
    rows = [
        ("constant", 1, Fraction(1), Fraction(1)),
        ("coordinate_in_X", 1, Fraction(m, p), Fraction(m + 1, 2 * m)),
        ("coordinate_outside_X", 0, Fraction(m, p), Fraction(1, 2)),
        ("pair_inside_X", 1, Fraction(m * (m - 1), p * (p - 1)),
         Fraction(m + 2, 4 * m)),
        ("pair_outside_X", 0, Fraction(m * (m - 1), p * (p - 1)), Fraction(1, 4)),
        ("pair_crossing_X", 0, Fraction(m * (m - 1), p * (p - 1)), Fraction(1, 4)),
    ]
    _require(all(value == (at_X + p * mean) / (p + 1)
                 for _, at_X, mean, value in rows), "paired-cube mean identity failed")
    return {"p": p, "identity": "T f(X)=(f(X)+p E[f])/(p+1)",
            "monomial_orbits": [dict(orbit=name, at_X=at_X, slice_mean=str(mean),
                                     paired_mean=str(value))
                                for name, at_X, mean, value in rows],
            "proved": True}


def _stabilizer_dependency(p: int) -> dict[str, object]:
    old = sharp_integral_quadratic_lift_floor(p)
    dual = stabilizer_mass_certificate(p)
    coefficient = Fraction(4) if p % 4 == 3 else Fraction(4 * (p - 1), p + 3)
    _require(old.get("proved") is True
             and old.get("paired_cube_identity") == "T B(X)=(B(X)+p E[B])/(p+1)"
             and old.get("integer_quadratic_cube_mean_lattice") == "(1/4) Z"
             and old.get("H_at_least_two_stabilizer_coefficient") == coefficient,
             "15.688 paired/stabilizer dependency failed")
    nodes, weights = dual["nodes"], dual["weights"]
    represented = tuple(sum(w * x**degree for x, w in zip(nodes, weights))
                        for degree in range(3))
    _require(represented == hypergeometric_moments(p)
             and all(w >= 0 for w in weights)
             and weights[-1] == dual["value"] == coefficient / (4 * p),
             "15.642 stabilizer moment certificate failed")
    return {"p": p, "endpoint_coefficient": str(dual["value"]),
            "nodes": list(nodes), "weights": list(map(str, weights)),
            "moments": list(map(str, represented)),
            "uniform_height_bound": "H<=p*mu*(p+3)/(p-1)", "proved": True}


def slice_height_bound_certificate(p: int, s: int) -> dict[str, object]:
    _odd(p)
    _shift(s)
    cube = cube_height_bound_certificate()
    paired, stabilizer = paired_cube_basis_certificate(p), _stabilizer_dependency(p)
    _require(all(row.get("proved") is True for row in (cube, paired, stabilizer))
             and cube.get("height_over_mean_squared") == 324
             and paired.get("identity") == "T f(X)=(f(X)+p E[f])/(p+1)",
             "height proof dependency failed")
    mu = Fraction(2 * p + s, 2 * p)
    preliminary = p * mu * Fraction(p + 3, p - 1)
    cube_mean = (preliminary + p * mu) / (p + 1)
    height = Fraction(cube["height_over_mean_squared"]) * cube_mean**2
    _require(mu <= Fraction(32, 29) < MEAN_UPPER == Fraction(9, 8)
             and cube_mean == Fraction(2 * p + s, p - 1) <= Fraction(16, 7)
             and height <= Fraction(82944, 49) < HEIGHT_BOUND == 1800,
             "uniform slice height bound failed")
    return {"p": p, "s": s, "scaled_mass": 2 * p + s, "mean": str(mu),
            "paired_cube_mean_upper": str(cube_mean), "height_upper_exact": str(height),
            "uniform_height_bound": HEIGHT_BOUND,
            "maximizing_point_belongs_to_every_selected_cube": True,
            "cube_dependency": cube, "paired_dependency": paired,
            "stabilizer_dependency": stabilizer, "proved": True}


def derivative_support_certificate(p: int) -> dict[str, object]:
    """Nonzero integer affine differences, not a Boolean-valued assumption."""
    _odd(p)
    m, N, K = (p + 1) // 2, p - 2, (p - 1) // 2
    one_of_pair = Fraction(2 * K * (N - K), N * (N - 1))
    support = one_of_pair / 2
    conditioning = Fraction(2 * m * (p - m), p * (p - 1))
    influence = conditioning * support / 4
    _require(one_of_pair == Fraction(p - 1, 2 * (p - 2))
             and support == Fraction(p - 1, 4 * (p - 2)) > Fraction(1, 4)
             and conditioning == Fraction(p + 1, 2 * p)
             and influence == Fraction(p * p - 1, 32 * p * (p - 2)) > Fraction(1, 32),
             "affine derivative support/influence bound failed")
    return {"p": p, "conditional_slice": [N, K],
            "exactly_one_of_unequal_coefficient_pair": str(one_of_pair),
            "conditional_nonzero_support_lower": str(support),
            "transposition_conditioning_probability": str(conditioning),
            "squared_influence_definition": "I_ij=(1/4)E[(A-A^(ij))^2]",
            "relevant_pair_influence_lower": str(influence),
            "uniform_strict_influence_lower": "1/32",
            "Johnson_degree_two_eigenvalue": 2 * (p - 1),
            "total_influence_variance_factor": p - 1,
            "constant_nonzero_affine_difference_has_full_support": True,
            "integer_nonzero_difference_has_square_at_least_one": True, "proved": True}


def bounded_junta_certificate(p: int, s: int) -> dict[str, object]:
    height, derivative = slice_height_bound_certificate(p, s), derivative_support_certificate(p)
    _require(height.get("proved") is True and derivative.get("proved") is True
             and height.get("uniform_height_bound") == HEIGHT_BOUND
             and derivative.get("uniform_strict_influence_lower") == "1/32"
             and derivative.get("total_influence_variance_factor") == p - 1,
             "bounded-junta height/influence dependency failed")
    mu = Fraction(height["mean"])
    influence = Fraction(derivative["relevant_pair_influence_lower"])
    L_upper = Fraction(2 * (p - 1) * HEIGHT_BOUND, p) * mu / influence
    _require(influence > Fraction(1, 32)
             and L_upper < 64 * HEIGHT_BOUND * mu
             and 64 * HEIGHT_BOUND * mu <= 64 * HEIGHT_BOUND * MEAN_UPPER
             == JUNTA_STRICT_BOUND == 129600,
             "strict invariant-class complement bound failed")
    q = (p - 1) // 2
    return {"p": p, "s": s, "height_dependency": height,
            "derivative_dependency": derivative,
            "relevant_pair_count_lower": "p*L/2",
            "largest_invariant_class_complement_upper": str(L_upper),
            "junta_size_strict_upper": JUNTA_STRICT_BOUND,
            "junta_coordinates_at_most": JUNTA_STRICT_BOUND - 1,
            "q": q, "all_kept_patterns_extend": q >= JUNTA_STRICT_BOUND - 1,
            "zero_influence_is_an_equivalence_relation": True,
            "symmetrization_then_sK_equals_m_minus_sJ_preserves_degree_two": True,
            "proved": True}


def middle_slice_mean_certificate(p: int) -> dict[str, object]:
    _odd(p)
    m = (p + 1) // 2
    slice_moments = [Fraction(1), Fraction(m, p), Fraction(m * (m - 1), p * (p - 1))]
    cube_moments = [Fraction(1), Fraction(1, 2), Fraction(1, 4)]
    at_zero = [1, 0, 0]
    _require(slice_moments == [1, Fraction(p + 1, 2 * p), Fraction(p + 1, 4 * p)]
             and all(p * a == (p + 1) * b - c
                     for a, b, c in zip(slice_moments, cube_moments, at_zero)),
             "exact middle-slice/cube mean identity failed")
    return {"p": p, "monomials": ["1", "x_i", "x_i*x_j"],
            "slice_moments": list(map(str, slice_moments)),
            "cube_moments": list(map(str, cube_moments)),
            "at_zero": at_zero, "identity": "p*mu=(p+1)*nu-f(0)",
            "full_cube_mean_lattice": "(1/4) Z", "proved": True}


def local_mass_exclusion(p: int, s: int) -> dict[str, object]:
    """A local theorem for every odd p>=259201; primality is unnecessary."""
    _require(PRIME_THRESHOLD == 259201, "proved order threshold changed")
    _odd(p, PRIME_THRESHOLD)
    _shift(s)
    junta, mean = bounded_junta_certificate(p, s), middle_slice_mean_certificate(p)
    _require(junta.get("proved") is True and mean.get("proved") is True
             and junta.get("all_kept_patterns_extend") is True
             and junta.get("junta_size_strict_upper") == JUNTA_STRICT_BOUND
             and mean.get("identity") == "p*mu=(p+1)*nu-f(0)"
             and mean.get("full_cube_mean_lattice") == "(1/4) Z",
             "local mean/junta/extension dependency failed")
    numerator_min, numerator_max = s // 2 - 1, HEIGHT_BOUND + s // 2 - 1
    gap_min, gap_max = Fraction(numerator_min, p + 1), Fraction(numerator_max, p + 1)
    _require(1 <= numerator_min <= numerator_max <= HEIGHT_BOUND + 2
             and p + 1 > 4 * (HEIGHT_BOUND + 2)
             and 0 < gap_min <= gap_max < Fraction(1, 4)
             and 4 < 4 * (1 + gap_min) <= 4 * (1 + gap_max) < 5,
             "strict quarter-grid contradiction failed")
    return {"p": p, "s": s, "scaled_mass": 2 * p + s,
            "hypotheses": "A>=0 integer-valued degree<=2 on J(p,(p+1)/2)",
            "affine_parity_hypothesis_needed": False,
            "cube_nu_minus_one_interval": [str(gap_min), str(gap_max)],
            "contradiction": "1<nu<5/4 but nu is quarter-integral",
            "junta_dependency": junta, "mean_dependency": mean,
            "excluded": True, "proved": True}


def _capacity_spectrum_dependency() -> dict[str, object]:
    """Read/replay the existing all-prime theorem at one fixed small receipt.

    Its proof is symbolic for p>=29.  Calling the full boundary-table API at
    the new enormous p is unnecessary: no O(p) boundary or residue scan is
    used here.  The quotient interval algebra is checked separately below.
    """
    old = affine_parity_small_mass_spectrum(37)
    _require(old.get("proved") is True
             and old.get("strict_upper_mass") == 64
             and old.get("union_allowed_masses") == [0, 34, 36, 38]
             and old.get("all_even_boundary_sizes_covered") is True
             and old.get("new_catalog_or_prime_census_used") is False,
             "15.774 all-prime spectrum dependency failed")
    return {"all_prime_theorem_minimum": 29, "fixed_receipt_prime": 37,
            "strict_upper_mass_formula": "2p-10",
            "union_mass_formula": "{0,p-3,p-1,p+1}",
            "new_catalog_or_prime_census_used": False, "proved": True}


def first_layer_quota_certificate(p: int) -> dict[str, object]:
    _odd(p, 37)
    if not is_prime(p):
        raise ValueError("two-type capacity needs prime p>=37")
    dependency = _capacity_spectrum_dependency()
    _require(dependency.get("proved") is True, "capacity spectrum dependency failed")
    m, t = (p + 1) // 2, p + 3
    endpoints = {0: 0, 1: 2, 2: 2, 3: 2, m - 7: 2,
                 m - 6: 1, m - 3: 1, m - 2: 0, m - 1: 0}
    _require(all(quotient_floor(p, u) == k for u, k in endpoints.items()),
             "15.774 quotient floor dependency failed")
    # Entire residue intervals are separated algebraically, not enumerated.
    _require(m >= 19 and t == 2 * m + 2
             and 2 * (m - 7) + 2 * m < 2 * p - 10
             and 2 * (m - 6) + 2 * m == 2 * p - 10
             and 2 + 2 * m <= t < 3 + 2 * m
             and m - 6 > 3 and m - 12 > 3,
             "first-layer low/high residue separation failed")
    lows = [0, 1, 2]
    collisions = [[u, v] for u in lows for v in lows if u + v == 3]
    u1_extra, u2_extra = t - 1 - 2 * m, t - 2 - 2 * m
    _require(collisions == [[1, 2], [2, 1]] and u1_extra == 1 and u2_extra == 0
             and 3 * p + 2 * t == 5 * p + 6
             and 2 + 2 * (p + 1) == 2 * p + 4
             and 4 + 2 * (p + 1) == 2 * p + 6,
             "first-layer exact quota/mass conclusion failed")
    return {"p": p, "m": m, "t_shell": t, "H": 5 * p + 6,
            "allowed_residue_intervals": [[0, 2], [m - 6, m - 1]],
            "only_ordered_residue_pairs": collisions,
            "u1_quotient_histogram": {"2": m - 1, "3": 1},
            "u2_quotient_histogram": {"2": m},
            "forced_low_masses": [2 * p + 4, 2 * p + 6],
            "P_and_Q_sum": 10,
            "signed_T_formula": "sigma*T=(p+1)*(P-5)+1, sigma is the u1 type",
            "P_equals_5_or_T_equals_pm1_assumed": False,
            "large_residue_or_boundary_scan_used": False,
            "spectrum_dependency": dependency, "proved": True}


def eventual_first_layer_exclusion(p: int) -> dict[str, object]:
    _odd(p, PRIME_THRESHOLD)
    if not is_prime(p):
        raise ValueError("residual layer needs prime p>=259201")
    quotas = first_layer_quota_certificate(p)
    local = [local_mass_exclusion(p, s) for s in (4, 6)]
    isolation_margin = p * p + 1 - 2 * (5 * p + 6)
    _require(quotas.get("proved") is True
             and quotas.get("H") == 5 * p + 6
             and quotas.get("only_ordered_residue_pairs") == [[1, 2], [2, 1]]
             and quotas.get("u1_quotient_histogram") == {"2": (p - 1) // 2, "3": 1}
             and quotas.get("u2_quotient_histogram") == {"2": (p + 1) // 2}
             and quotas.get("forced_low_masses") == [2 * p + 4, 2 * p + 6]
             and all(row.get("proved") is True and row.get("excluded") is True
                     and row.get("scaled_mass") == 2 * p + s
                     for row, s in zip(local, (4, 6)))
             and isolation_margin > 0,
             "eventual first-layer closure dependency failed")
    return {"p": p, "k": 5 * p + 5, "H": 5 * p + 6,
            "t_residual": (p - 1) // 2 + 3,
            "isolation_margin": isolation_margin,
            "both_signed_shell_floors": 3,
            "all_boundary_sizes_excluded": True,
            "whole_layer_excluded": True,
            "quota_dependency": quotas, "local_dependencies": local,
            "proved": True}


def generic_bounded_mean_certificate(p: int, mean_cap: Fraction) -> dict[str, object]:
    """Quantitative height/junta bounds, without a fixed local mass."""
    _odd(p)
    if type(mean_cap) not in (int, Fraction) or mean_cap <= 0:
        raise ValueError("need a positive exact mean cap")
    B = Fraction(mean_cap)
    cube, paired = cube_height_bound_certificate(), paired_cube_basis_certificate(p)
    stabilizer, derivative = _stabilizer_dependency(p), derivative_support_certificate(p)
    mean = middle_slice_mean_certificate(p)
    _require(all(row.get("proved") is True
                 for row in (cube, paired, stabilizer, derivative, mean))
             and cube.get("height_over_mean_squared") == 324
             and paired.get("identity") == "T f(X)=(f(X)+p E[f])/(p+1)"
             and derivative.get("uniform_strict_influence_lower") == "1/32"
             and derivative.get("total_influence_variance_factor") == p - 1
             and mean.get("identity") == "p*mu=(p+1)*nu-f(0)"
             and mean.get("full_cube_mean_lattice") == "(1/4) Z",
             "generic bounded-mean dependency failed")
    paired_upper = Fraction(2 * p, p - 1) * B
    height_upper, L_strict = 2916 * B**2, 186624 * B**3
    _require(paired_upper <= 3 * B
             and 324 * paired_upper**2 <= height_upper
             and 64 * height_upper * B == L_strict,
             "generic height/junta constants failed")
    return {"p": p, "mean_cap": str(B), "height_upper": str(height_upper),
            "junta_size_strict_upper": str(L_strict),
            "height_coefficient": 2916, "junta_coefficient": 186624,
            "cube_dependency": cube, "paired_dependency": paired,
            "stabilizer_dependency": stabilizer, "derivative_dependency": derivative,
            "mean_dependency": mean, "proved": True}


def power_band_exclusion(p: int, r: int, h: int) -> dict[str, object]:
    """Both signed shell floors r; preserve the indispensable h=r mod 2."""
    _odd(p)
    if not is_prime(p):
        raise ValueError("signed shell theorem needs prime p>=29")
    if type(r) is not int or r not in (3, 4, 5):
        raise ValueError("need shell floor r in {3,4,5}")
    if type(h) is not int or h < 0 or (h - r) % 2:
        raise ValueError("need nonnegative h congruent to r modulo two")
    if h < r * p:
        _require(Fraction(h, p) < r, "signed frame contradiction failed")
        return {"p": p, "r": r, "H": h, "case": "below the signed frame mean",
                "both_signed_shell_floors_required": True,
                "parity_required": "h=r mod 2", "excluded": True, "proved": True}
    if 46656 * h**3 > p**3 * (p - 1):
        raise ValueError("outside the proved cubic support band")
    B = Fraction(h, 2 * p)
    average = Fraction(h - r * p, 2 * p)
    dep = generic_bounded_mean_certificate(p, B)
    _require(dep.get("proved") is True
             and dep.get("height_coefficient") == 2916
             and dep.get("junta_coefficient") == 186624,
             "power-band bounded-mean dependency failed")
    height, L_strict = Fraction(dep["height_upper"]), Fraction(dep["junta_size_strict_upper"])
    q, M = (p - 1) // 2, p + 1
    left_upper = 4 * (2 * height + r)
    isolation = p * p + 1 - 2 * h
    _require(0 <= average <= B and B >= Fraction(r, 2) >= Fraction(3, 2)
             and height == 2916 * B**2 and L_strict == 186624 * B**3
             and 373248 * B**3 <= p - 1
             and L_strict <= q
             and left_upper == 4 * (5832 * B**2 + r)
             and left_upper <= 23336 * B**2 <= 23336 * B**3
             < 373248 * B**3 <= p - 1 < M
             and 4 * B <= p - 1 and isolation > 0,
             "power-band extension/divisibility contradiction failed")
    # For the two selected rows, eliminate T from
    # a_eps=M*P_eps-eps*T-r*p=2*M*nu_eps-2*alpha_eps.
    # The RHS bracket below is integral because each 4*nu_eps is integral.
    return {"p": p, "r": r, "H": h, "case": "quarter-grid divisibility",
            "parity_required": "h=r mod 2", "both_signed_shell_floors_required": True,
            "support_condition": "46656*h^3<=p^3*(p-1)",
            "type_average_A": str(average), "selected_mean_cap": str(B),
            "select_one_low_row_in_each_type": True,
            "isolation_margin": isolation, "junta_strict_upper": str(L_strict),
            "q": q, "full_cube_extensions_exist": True,
            "positive_multiple_lower": 4 * r,
            "positive_multiple_upper": str(left_upper), "modulus": M,
            "divisibility_identity": (
                "4*(alpha_plus+alpha_minus+r)="
                "(p+1)*(4*nu_plus+4*nu_minus-2*(P_plus+P_minus)+4*r)"),
            "neither_phase_omitted": True, "bounded_mean_dependency": dep,
            "excluded": True, "proved": True}


def proposition_15775() -> dict[str, object]:
    rows = [eventual_first_layer_exclusion(p) for p in SAMPLE_PRIMES]
    power = [power_band_exclusion(6700417, r, r * 6700417 + 2) for r in (3, 4, 5)]
    _require(len(rows) == 2 and all(row.get("proved") is True for row in rows)
             and len(power) == 3 and all(row.get("proved") is True for row in power)
             and {p % 4 for p in SAMPLE_PRIMES} == {1, 3},
             "first-layer representative proof checks failed")
    return {"prop": "15.775", "status": "PROVED_INFINITE_FAMILY",
            "proof_note": "evidence/NOTE_2026-09-04_EVENTUAL_FIRST_LAYER_CLOSE.md",
            "minimum_prime": PRIME_THRESHOLD,
            "theorem": "p>=259201 prime: signed shell floor 3 is impossible at |H|=5p+6",
            "residual_layer": "k=5p+5, t=(p-1)/2+3",
            "eventual_next_residual_frontier": "p>=259201: k>=5p+7",
            "smaller_prime_frontiers_changed": False,
            "uniform_height_bound": HEIGHT_BOUND,
            "junta_size_strict_upper": JUNTA_STRICT_BOUND,
            "records": rows, "records_are_identity_replays_not_a_prime_census": True,
            "superlinear_support_theorem": (
                "p>=29 prime,r in {3,4,5},h=r mod2: both signed floors r are "
                "impossible below rp or when 46656*h^3<=p^3*(p-1)"),
            "superlinear_support_records": power,
            "new_equality_catalog_used": False, "proved": True,
            "residual_ii_closed_general": False,
            "minimal_four_gap_bridge_closed_general": False,
            "eventual_E1_proved": False, "e1_closed_general": False,
            "original_MO_limit_closed": False}


if __name__ == "__main__":
    result = proposition_15775()
    write_json_atomic(ROOT / "evidence/e1_gmin_m4_prop15775.json", result)
    print(json.dumps({"prop": result["prop"], "proved": result["proved"],
                      "minimum_prime": result["minimum_prime"]}, sort_keys=True))
