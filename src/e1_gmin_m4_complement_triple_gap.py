#!/usr/bin/env python3
r"""Complement-triple punctured gap and two p=1 mod 4 local exclusions.

For odd p>=29, a nonnegative integral quadratic A on J(p,(p+1)/2)
with parity r=|X intersect C|, |C|=3, and scaled mean 2p-6+delta,
0<=delta<=4, is either the old (r-2)^2 baseline (delta=0), or one
of the three pair-plus-complement-literal forms (delta=4). In particular
delta=2 is impossible. The difference from the old baseline is NOT
globally nonnegative; the neighboring-slice section bound and explicit
small-side kernel reduction supply the missing implication.

For primes p=1 mod 4, p>=29, local lift masses p-1 and p+11 are
impossible. These statements use the previously proved half-mean cube
theorem and fixed four-bit catalog, not a new configuration census.

Proof: evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md.
This module deliberately does not depend on Proposition 15.770.
"""
from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15751 import (
    cube_half_mean_height_certificate,
    exact_four_cube_catalog,
    profile_density,
)
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
PROOF_NOTE = "evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _check_odd_order(p: int) -> None:
    if not isinstance(p, int) or isinstance(p, bool) or p < 29 or p % 2 == 0:
        raise ValueError("need an odd integer p>=29")


def _check_p1_prime(p: int) -> None:
    _check_odd_order(p)
    if p % 4 != 1 or not is_prime(p):
        raise ValueError("need a prime p>=29 congruent to 1 modulo 4")


def _evaluate_polynomial(coefficients: tuple[int, ...], value: int) -> int:
    result = 0
    for coefficient in coefficients:
        result = value * result + coefficient
    return result


def _positive_shift_certificate(coefficients: tuple[int, ...]) -> dict[str, object]:
    """Expand a polynomial identically at p=29+x, using integer arithmetic.

    Coefficients are in descending order. Positivity of the full translated
    coefficient vector proves the inequality on the entire ray x>=0;
    evaluating it at a few primes is not the proof.
    """
    shifted_ascending = [0] * len(coefficients)
    for power, coefficient in enumerate(reversed(coefficients)):
        for smaller_power in range(power + 1):
            shifted_ascending[smaller_power] += (
                coefficient * comb(power, smaller_power) * 29 ** (power - smaller_power)
            )
    shifted = tuple(reversed(shifted_ascending))
    _require(
        all(value >= 0 for value in shifted) and shifted[-1] > 0,
        "the generic positive polynomial certificate failed",
    )
    return {
        "coefficients_in_p_descending": list(coefficients),
        "translated_variable": "x=p-29>=0",
        "coefficients_in_x_descending": list(shifted),
        "translation_checked_coefficientwise": True,
        "strictly_positive_for_every_p_at_least_29": True,
        "proved": True,
    }


@lru_cache(maxsize=None)
def complement_triple_gap_certificate(p: int) -> dict[str, object]:
    """Exclude delta=2 and classify delta=4, for every odd order p>=29."""
    _check_odd_order(p)
    m = (p + 1) // 2
    outside = p - 3
    section_order = p - 4
    section_floor = sharp_integral_quadratic_lift_floor(section_order)
    floor_density = Fraction(section_order - 3, 4 * section_order)
    neighbor_density = Fraction(outside - 2, outside) * floor_density
    contact_contribution = outside * neighbor_density
    gap = _positive_shift_certificate((1, -28, 99))
    quadrature = parity_floor_certificate(p, 3, 0)
    nodes = tuple(quadrature["quadrature_nodes"])
    weights = tuple(quadrature["quadrature_weights"])
    moments = tuple(
        sum(weight * node**degree for node, weight in zip(nodes, weights))
        for degree in range(3)
    )
    forms = []
    for pair in combinations(range(3), 2):
        omitted = next(index for index in range(3) if index not in pair)
        values = []
        differences = []
        for mask in range(8):
            bits = tuple((mask >> index) & 1 for index in range(3))
            r = sum(bits)
            value = (1 - bits[pair[0]] - bits[pair[1]]) ** 2 + 1 - bits[omitted]
            difference = Fraction(value - (r - 2) ** 2, 2)
            signed = tuple(2 * bit - 1 for bit in bits)
            _require(
                value >= 0 and value % 2 == r % 2
                and difference.denominator == 1
                and 3 + 2 * value == 5 + signed[pair[0]] * signed[pair[1]] - signed[omitted],
                "a complement-triple delta-four form changed",
            )
            values.append(value)
            differences.append(int(difference))
        third_difference = sum(
            (-1) ** (3 - mask.bit_count()) * value
            for mask, value in enumerate(values)
        )
        _require(third_difference == 0, "a delta-four form ceased to be quadratic")
        forms.append({
            "pair": list(pair),
            "complement_literal_coordinate": omitted,
            "values_in_binary_mask_order": values,
            "half_difference_from_old_baseline": differences,
            "signed_target": f"5+z_{pair[0]}*z_{pair[1]}-z_{omitted}",
            "coefficient_offset": 4,
        })
    pair_literal_mean = Fraction(p - 1, p)
    proved = bool(
        section_floor["proved"]
        and section_floor["sharp_mass_floor"] == floor_density
        and contact_contribution == Fraction((p - 5) * (p - 7), 4 * (p - 4))
        and contact_contribution - 4
        == Fraction(_evaluate_polynomial((1, -28, 99), p), 4 * (p - 4))
        and contact_contribution > 4
        and gap["proved"]
        and 2 <= m - 3 <= outside - 2
        and 0 < m - 1 < outside
        and 0 <= m - 3 <= m <= outside
        and quadrature["exact_positive_quadrature_certificate"]
        and nodes == (1, 2, 3) and all(weight > 0 for weight in weights)
        and moments == (Fraction(1), Fraction(3 * (p + 1), 2 * p), Fraction(3 * (p + 1), p))
        and len(forms) == 3
        and 2 * p * pair_literal_mean == 2 * p - 2
    )
    _require(proved, "the complement-triple punctured gap certificate failed")
    return {
        "p": p,
        "scope": "all odd integers p>=29; large-boundary phase-one translation uses p=1 mod4",
        "small_side_size": 3,
        "small_side_parity": "r mod2",
        "old_baseline": "A0=(r-2)^2",
        "old_baseline_is_pointwise_parity_minimum_globally": False,
        "difference_L_is_nonnegative_on_contacts_only": [1, 2, 3],
        "neighboring_slice_bound": {
            "outside_order": outside,
            "conditional_weights": [m - 1, m - 3],
            "odd_section_order": section_order,
            "odd_section_need_not_be_prime": True,
            "identically_zero_coordinate_sections_at_most": 2,
            "zero_section_argument": "three forbidden coordinate sections give a three-swap cube zero away from its origin",
            "section_averaging_identity": "mean of all coordinate one-sections or all zero-sections equals E[f]",
            "nonzero_odd_section_mean_floor": str(floor_density),
            "nonzero_neighbor_mean_floor": str(neighbor_density),
            "one_nonzero_outer_contact_contribution_at_least": str(contact_contribution),
            "strict_gap_above_four": gap,
            "proved": True,
        },
        "quadrature_nodes": list(nodes),
        "quadrature_weights": [str(value) for value in weights],
        "conditional_means_are_symmetrized_over_small_patterns_of_equal_weight": True,
        "exact_excess_identity": "delta=(p-3)*(sum_singletons mu_S+mu_C)+4*sum_pairs mu_S",
        "delta_at_most_four_forces_singleton_and_triple_sections_zero": True,
        "globalization": {
            "complemented_bits": "w_i=1-x_i for the three small-side coordinates only",
            "outside_slice_identity": "sum y=m-3+sum w",
            "kernel_elimination_sign": "PLUS",
            "contacts_in_w": [0, 2],
            "all_small_side_patterns_extend": True,
            "conclusion": "L depends only on the three small-side coordinates",
        },
        "residual_integer_identity": "delta=4*sum(nonnegative integer pair values)",
        "allowed_excesses_in_zero_to_four": [0, 4],
        "excess_two_excluded": True,
        "excess_four_excluded": False,
        "excess_zero_form": "(r-2)^2",
        "excess_zero_coefficient_offset": 2,
        "excess_four_forms": forms,
        "excess_four_coefficient_offsets": [4],
        "proof_note": PROOF_NOTE,
        "new_graph_prime_slice_or_catalog_census_used": False,
        "proved": proved,
    }


def _boolean_mass_exclusion(p: int, offset: int) -> dict[str, object]:
    """Use the already proved influence reduction and fixed density formulas."""
    mass = p + offset
    density = Fraction(mass, 4 * p)
    denominator = p * p * (p + 1) * (p - 3)
    numerator = 2 * (p - 1) * (p - 2) * mass * (4 * p - mass)
    junta_bound = Fraction(numerator, denominator)
    if offset == -1:
        upper = 6
        gap_coefficients = (10, -40, 2, 4)
        bracket = (Fraction(p - 3, 4 * p), Fraction(p + 1, 4 * p))
    elif offset == 11:
        upper = 8
        gap_coefficients = (2, -42, 338, -814, 484)
        bracket = (Fraction(p + 1, 4 * p), Fraction(p - 1, 2 * p))
    else:  # pragma: no cover - private fixed-scope helper
        raise ValueError("unsupported local mass offset")
    positive_gap = _positive_shift_certificate(gap_coefficients)
    catalog = exact_four_cube_catalog()
    densities = sorted({
        profile_density(tuple(row["layer_counts"]), p)
        for row in catalog["profiles"]
    })
    expected = sorted({
        Fraction(0), Fraction(1), Fraction(p - 3, 4 * p),
        Fraction(p + 1, 4 * p), Fraction(p - 1, 2 * p),
        Fraction(p + 1, 2 * p), Fraction(3 * p - 1, 4 * p),
        Fraction(3 * (p + 1), 4 * p),
    })
    proved = bool(
        positive_gap["proved"]
        and upper * denominator - numerator == _evaluate_polynomial(gap_coefficients, p)
        and junta_bound < upper
        and upper - 1 < (p - 1) // 2
        and catalog["proved"] and densities == expected
        and bracket[0] < density < bracket[1]
        and density not in densities
    )
    _require(proved, "a local Boolean mass survived the fixed density catalog")
    return {
        "density": str(density),
        "relevant_pair_influence_floor": str(Fraction((p + 1) * (p - 3), 16 * p * (p - 2))),
        "total_influence_upper_bound": str((p - 1) * density * (1 - density)),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "strict_junta_upper_bound": upper,
        "generic_junta_gap_polynomial": positive_gap,
        "junta_coordinates_at_most": upper - 1,
        "all_junta_patterns_extend_to_both_slice_sides": True,
        "cube_active_coordinates_at_most": 4,
        "four_bit_density_values": [str(value) for value in densities],
        "density_formula_dependency": "Proposition15.751 fixed fourteen-profile symbolic density formulas",
        "target_strict_bracket": [str(value) for value in bracket],
        "target_density_absent": True,
        "fixed_catalog_sha256": catalog["valid_table_signature_sha256"],
        "fixed_four_bit_catalog_reused": True,
        "new_catalog_used": False,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def p1_p_minus_one_local_exclusion(p: int) -> dict[str, object]:
    """Exclude lift mass p-1, including the non-Boolean equality height."""
    _check_p1_prime(p)
    mass = p - 1
    sharp = sharp_integral_quadratic_lift_floor(p)
    half = cube_half_mean_height_certificate()
    height_lower = Fraction(2 * (p + 1) - mass, 4)
    height_upper = Fraction(mass * (p + 3), 4 * (p - 1))
    paired_mean = Fraction(height_lower + Fraction(mass, 4), p + 1)
    height_proved = bool(
        sharp["proved"] and half["proved"]
        and sharp["H_at_least_two_scaled_floor"] == mass
        and height_lower == height_upper == Fraction(p + 3, 4)
        and paired_mean == Fraction(1, 2)
        and height_lower >= 8 > int(half["maximum_upper_bound"])
    )
    boolean = _boolean_mass_exclusion(p, -1)
    proved = bool(height_proved and boolean["proved"])
    _require(proved, "the p=1 mod4 p-1 local mass survived")
    return {
        "p": p,
        "scaled_mass": mass,
        "scaled_mass_definition": "4p E[L]",
        "statement": "no nonzero nonnegative integral quadratic on J(p,(p+1)/2) has mass p-1",
        "height_at_least_two": {
            "paired_height_lower_bound": str(height_lower),
            "stabilizer_height_upper_bound": str(height_upper),
            "forced_height": int(height_lower),
            "every_maximizing_cube_mean_at_least": "1/2",
            "average_maximizing_cube_mean": str(paired_mean),
            "therefore_every_maximizing_cube_mean": "1/2",
            "half_mean_cube_maximum_upper_bound": int(half["maximum_upper_bound"]),
            "excluded": height_proved,
            "proved": height_proved,
        },
        "height_one_boolean": boolean,
        "proof_note": PROOF_NOTE,
        "new_graph_prime_slice_or_catalog_census_used": False,
        "excluded": proved,
        "proved": proved,
    }


@lru_cache(maxsize=None)
def p1_p_plus_eleven_local_exclusion(p: int) -> dict[str, object]:
    """Exclude lift mass p+11 by half-mean cubes and fixed Boolean densities."""
    _check_p1_prime(p)
    mass = p + 11
    sharp = sharp_integral_quadratic_lift_floor(p)
    half = cube_half_mean_height_certificate()
    height_lower = Fraction(p - 9, 4)
    height_upper = Fraction(mass * (p + 3), 4 * (p - 1))
    paired_average_upper = Fraction(mass, 2 * (p - 1))
    upper_gap = Fraction(3, 4) - paired_average_upper
    height_proved = bool(
        sharp["proved"] and half["proved"]
        and height_lower == Fraction(2 * (p + 1) - mass, 4)
        and height_lower >= 5 > int(half["maximum_upper_bound"])
        and Fraction(height_upper + Fraction(mass, 4), p + 1) == paired_average_upper
        and upper_gap == Fraction(p - 25, 4 * (p - 1)) > 0
    )
    boolean = _boolean_mass_exclusion(p, 11)
    proved = bool(height_proved and boolean["proved"])
    _require(proved, "the p=1 mod4 p+11 local mass survived")
    return {
        "p": p,
        "scaled_mass": mass,
        "scaled_mass_definition": "4p E[L]",
        "statement": "no nonzero nonnegative integral quadratic on J(p,(p+1)/2) has mass p+11",
        "height_at_least_two": {
            "paired_height_lower_bound": str(height_lower),
            "stabilizer_height_upper_bound": str(height_upper),
            "every_maximizing_cube_mean_at_least": "1/2",
            "cube_mean_lattice": "(1/4)Z",
            "average_maximizing_cube_mean_upper_bound": str(paired_average_upper),
            "strict_gap_below_three_quarters": str(upper_gap),
            "generic_strict_gap_identity": "3/4-(p+11)/(2(p-1))=(p-25)/(4(p-1))>0",
            "some_maximizing_cube_has_mean_exactly": "1/2",
            "half_mean_cube_maximum_upper_bound": int(half["maximum_upper_bound"]),
            "excluded": height_proved,
            "proved": height_proved,
        },
        "height_one_boolean": boolean,
        "proof_note": PROOF_NOTE,
        "new_graph_prime_slice_or_catalog_census_used": False,
        "excluded": proved,
        "proved": proved,
    }


def local_bridge_package() -> dict[str, object]:
    """Package generic proof certificates instantiated at their endpoint."""
    rows = {
        "complement_triple_punctured_gap": complement_triple_gap_certificate(29),
        "p1_p_minus_one": p1_p_minus_one_local_exclusion(29),
        "p1_p_plus_eleven": p1_p_plus_eleven_local_exclusion(29),
    }
    return {
        "classification": "proved generic local lemmas with exact symbolic inequalities",
        "proof_note": PROOF_NOTE,
        "endpoint_instantiations": rows,
        "all_parameter_inequalities_proved_by_exact_identities_not_endpoint_sampling": True,
        "complement_triple_excess_four_is_attained_not_excluded": True,
        "standalone_endpoint_or_global_closure_claimed": False,
        "proved": all(row["proved"] for row in rows.values()),
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = ROOT / "evidence" / "e1_gmin_m4_complement_triple_gap.json"
    write_json_atomic(path, local_bridge_package())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
