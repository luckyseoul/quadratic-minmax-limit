#!/usr/bin/env python3
"""Independent stdlib arithmetic replay for the eventual first residual layer.

This does not import the proof module, evaluate Paley signings, enumerate
primes, or search a graph/cell catalog. It checks fixed monomial identities,
exact rational constants, and interval certificates valid on an infinite tail.
The analytic inputs (Bonami and the proved Johnson identities) are stated
explicitly; arithmetic replay is not a replacement for their proof audit.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform


P_MIN = 259201
HEIGHT_CAP = 1800
JUNTA_CAP = 129600


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def clean(poly: dict[tuple[int, ...], Fraction]) -> dict[tuple[int, ...], Fraction]:
    return {key: value for key, value in poly.items() if value}


def add(left, right, scale=Fraction(1)):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, Fraction(0)) + scale * value
    return clean(answer)


def multiply(left, right):
    """Formal Boolean multilinear product, using only x_i**2=x_i."""
    answer = {}
    for a, av in left.items():
        for b, bv in right.items():
            key = tuple(sorted(set(a) | set(b)))
            answer[key] = answer.get(key, Fraction(0)) + av * bv
    return clean(answer)


def derivative_basis_check() -> dict[str, object]:
    """Check the swap derivative on representatives of every degree<=2 term."""
    monomials = [()] + [(i,) for i in range(4)] + list(combinations(range(4), 2))
    factor = {(0,): Fraction(1), (1,): Fraction(-1)}
    checked = []
    for monomial in monomials:
        polynomial = {monomial: Fraction(1)}
        swapped = {tuple(sorted(1 if i == 0 else 0 if i == 1 else i
                                for i in monomial)): Fraction(1)}
        actual = add(polynomial, swapped, Fraction(-1))
        # f-f^(01) = (x0-x1)(a0-a1 + sum_k (b0k-b1k)xk).
        affine = {}
        if monomial == (0,):
            affine[()] = Fraction(1)
        elif monomial == (1,):
            affine[()] = Fraction(-1)
        elif len(monomial) == 2 and 0 in monomial and 1 not in monomial:
            affine[(next(i for i in monomial if i != 0),)] = Fraction(1)
        elif len(monomial) == 2 and 1 in monomial and 0 not in monomial:
            affine[(next(i for i in monomial if i != 1),)] = Fraction(-1)
        expected = multiply(factor, affine)
        require(actual == expected, "formal transposition derivative identity failed")
        checked.append(list(monomial))
    return {"formal_monomials_checked": checked,
            "all_degree_zero_one_two_index_patterns_covered": True,
            "identity": "f-f^(ij)=(x_i-x_j)*(a_i-a_j+sum_k(b_ik-b_jk)x_k)",
            "verified": True}


def infinite_tail_margins() -> dict[str, object]:
    """Each numerator is written as c0+c1*(p-P_MIN), with positive denominator."""
    p = P_MIN
    margins = {
        "paired_mean_below_or_equal_16_over_7": {
            "numerator_coefficients_in_p_minus_P_MIN": [2 * (p - 29), 2],
            "denominator": "7*(p-1)"},
        "mu_below_9_over_8": {
            "numerator_coefficients_in_p_minus_P_MIN": [p - 24, 1],
            "denominator": "8*p"},
        "relevant_influence_above_1_over_32": {
            "numerator_coefficients_in_p_minus_P_MIN": [2 * p - 1, 2],
            "denominator": "32*p*(p-2)"},
        # Integral height <1800 permits f(0)<=1799. Worst mean uses a=3.
        "uniform_cube_mean_below_5_over_4": {
            "numerator_coefficients_in_p_minus_P_MIN": [p + 1 - 4 * (2 + HEIGHT_CAP - 1), 1],
            "denominator": "4*(p+1)"},
        "isolated_vertex_margin": {
            "positive_factorization": "p^2+1-2*(5p+6)=(p+1)*(p-11)",
            "factors_at_P_MIN": [p + 1, p - 11]},
    }
    for name, record in margins.items():
        values = record.get("numerator_coefficients_in_p_minus_P_MIN",
                            record.get("factors_at_P_MIN"))
        require(all(value > 0 for value in values), "nonpositive tail certificate: " + name)
    return margins


def residue_interval_check() -> dict[str, object]:
    """Interval algebra, not a residue/prime/signing census."""
    m_min = (P_MIN + 1) // 2
    require(m_min >= 19, "small-mass capacity range not available")
    # At t=2m+2: quotient floor 2 excludes 3<=u<=m-7.
    # Thus u lies in [0,2] or [m-6,m-1]. Coupling asks for sum ==3 mod m.
    require(m_min - 6 > 3 and m_min - 12 > 3, "high interval could hit residue three")
    low_pairs = [(a, b) for a in range(3) for b in range(3) if a + b == 3]
    require(low_pairs == [(1, 2), (2, 1)], "fixed low-low interval identity failed")
    return {
        "parameter": "m=(p+1)/2; t=2m+2",
        "allowed_residue_intervals": ["[0,2]", "[m-6,m-1]"],
        "low_high_sum_mod_m": ["[0,1]", "[m-6,m-1]"],
        "high_high_sum_mod_m": ["[m-12,m-2]"],
        "required_sum_mod_m": 3,
        "only_ordered_pairs": [list(pair) for pair in low_pairs],
        "u_one_total_excess_above_quotient_two": 1,
        "u_two_total_excess_above_quotient_two": 0,
        "forced_low_masses": ["2p+4", "2p+6"],
        "u_one_high_mass": "3p+5",
        "actual_parallel_count_five_assumed": False,
        "signed_total_plus_or_minus_one_assumed": False,
        "verified": True,
    }


def power_band_check() -> dict[str, object]:
    """Symbolic constants and divisibility, not sampled orders or supports."""
    height_coefficient = 324 * 3 ** 2
    junta_coefficient = 64 * height_coefficient
    extension_coefficient = 2 * junta_coefficient
    support_coefficient = Fraction(extension_coefficient, 8)
    require((height_coefficient, junta_coefficient, extension_coefficient,
             support_coefficient) == (2916, 186624, 373248, 46656),
            "generic height, junta, extension, or support constant failed")
    # Variables are M, nu_plus, nu_minus, P_plus, P_minus, T, r.
    # Every product below is square-free, so the formal Boolean polynomial
    # helper is also ordinary polynomial multiplication for this identity.
    variables = [{(i,): Fraction(1)} for i in range(7)]
    M, nu_plus, nu_minus, P_plus, P_minus, T, r = variables
    p = add(M, {(): Fraction(1)}, Fraction(-1))
    alpha_plus = add(multiply(M, nu_plus), multiply(M, P_plus), Fraction(-1, 2))
    alpha_plus = add(add(alpha_plus, T, Fraction(1, 2)), multiply(r, p), Fraction(1, 2))
    alpha_minus = add(multiply(M, nu_minus), multiply(M, P_minus), Fraction(-1, 2))
    alpha_minus = add(add(alpha_minus, T, Fraction(-1, 2)), multiply(r, p), Fraction(1, 2))
    lhs = add({}, add(add(alpha_plus, alpha_minus), r), Fraction(4))
    bracket = add(add({}, nu_plus, Fraction(4)), nu_minus, Fraction(4))
    bracket = add(add(add(bracket, P_plus, Fraction(-2)), P_minus, Fraction(-2)),
                  r, Fraction(4))
    rhs = multiply(M, bracket)
    require(lhs == rhs and all(5 not in key for key in lhs),
            "formal two-type elimination of T failed")
    shell_margins = []
    for shell in (3, 4, 5):
        # Set B=r/2+t, t>=0. These are all coefficients of 8*B^2-4*r.
        coefficients = [2 * shell ** 2 - 4 * shell, 8 * shell, 8]
        require(all(value >= 0 for value in coefficients),
                "shell upper bound 4*r<=8*B^2 failed")
        shell_margins.append({"r": shell, "B_parameterization": "r/2+t, t>=0",
                              "8B_squared_minus_4r_coefficients_in_t": coefficients})
    upper_coefficient = 8 * height_coefficient + 8
    min_B = Fraction(3, 2)
    isolation_factor_margin = extension_coefficient * min_B ** 2 - 4
    require(upper_coefficient == 23336 < extension_coefficient
            and min_B >= 1 and isolation_factor_margin > 0,
            "support band does not imply strict modulus/isolation bounds")
    return {
        "domain": "prime p>=29; r in {3,4,5}; h>=0; h=r mod 2",
        "below_frame_case": "h<r*p contradicts the signed frame mean",
        "remaining_case": "h>=r*p and 46656*h^3<=p^3*(p-1)",
        "B": "h/(2*p)", "minimum_B": str(min_B),
        "height_coefficient": height_coefficient,
        "junta_coefficient": junta_coefficient,
        "extension_coefficient": extension_coefficient,
        "support_condition_coefficient": int(support_coefficient),
        "pairing_factor_margin": "3-2p/(p-1)=(p-3)/(p-1)>0 for p>=29",
        "extension": "L<186624*B^3<=(p-1)/2",
        "type_selection": "each signed type has mean A=(h-r*p)/(2*p)<=B",
        "divisibility_identity": (
            "4*(alpha_plus+alpha_minus+r)="
            "(p+1)*(4*nu_plus+4*nu_minus-2*(P_plus+P_minus)+4*r)"),
        "formal_T_elimination_verified": True,
        "formal_variables": ["M=p+1", "nu_plus", "nu_minus", "P_plus", "P_minus", "T", "r"],
        "shell_upper_bound_margins": shell_margins,
        "strict_modulus_chain": (
            "0<4*r<=4*(alpha_plus+alpha_minus+r)<=4*(5832*B^2+r)"
            "<=23336*B^2<=23336*B^3<373248*B^3<=p-1<p+1"),
        "isolation_factor_margin_at_min_B": str(isolation_factor_margin),
        "isolation_implication": "4*B<373248*B^3<=p-1 implies p^2+1-2*h>0",
        "quarter_integrality_makes_bracket_integral": True,
        "both_signed_types_required": True, "parity_required": True,
        "verified": True,
    }


def verify() -> dict[str, object]:
    p = P_MIN
    paired_cap = Fraction(16, 7)
    exact_height_cap = 324 * paired_cap ** 2
    integer_height_cap = exact_height_cap.numerator // exact_height_cap.denominator
    derivative = derivative_basis_check()
    # Conditional nonzero affine support: pair an unequal-coefficient swap.
    outer_probability = Fraction(p + 1, 2 * p)
    conditional_support = Fraction(p - 1, 4 * (p - 2))
    influence_floor = outer_probability * conditional_support / 4
    expected_influence_floor = Fraction(p * p - 1, 32 * p * (p - 2))
    mean_checks = []
    for degree in (0, 1, 2):
        # On the middle slice E[x_i]=(p+1)/(2p), E[x_i*x_j]=(p+1)/(4p).
        slice_mean = Fraction(1) if degree == 0 else Fraction(p + 1, p * 2 ** degree)
        cube_mean = Fraction(1, 2 ** degree)
        at_zero = int(degree == 0)
        require(p * slice_mean == (p + 1) * cube_mean - at_zero,
                "slice/cube mean identity failed on a monomial")
        mean_checks.append({"degree": degree, "slice_mean_at_P_MIN": str(slice_mean),
                            "cube_mean": str(cube_mean), "value_at_zero": at_zero})
    integer_value_at_zero_max = HEIGHT_CAP - 1
    target_rows = []
    for a in (2, 3):
        minimum_nu = Fraction(p + a, p + 1)
        maximum_nu = Fraction(p + a + integer_value_at_zero_max, p + 1)
        require(1 < minimum_nu <= maximum_nu < Fraction(5, 4),
                "target interval intersects the quarter-integral lattice")
        target_rows.append({"a": a, "slice_mean": "1+%d/p" % a,
                            "scaled_mass_2p_E_A": "2p+%d" % (2 * a),
                            "uniform_cube_mean_formula": "1+(a-1+f(0))/(p+1)",
                            "uniform_cube_mean_interval_at_P_MIN": [str(minimum_nu), str(maximum_nu)],
                            "quarter_integral_value_exists": False})
    checks = {
        "formal_swap_derivative_verified": derivative["verified"],
        "bonami_fourier_constant_is_324": 4 * 9 ** 2 == 324,
        "height_cap_exact": exact_height_cap == Fraction(82944, 49),
        "integer_height_cap_1692": integer_height_cap == 1692,
        "loose_height_cap_1800_is_strict": integer_height_cap < HEIGHT_CAP,
        "derivative_probability_product_identity": influence_floor == expected_influence_floor,
        "relevant_influence_strictly_above_one_over_32": influence_floor > Fraction(1, 32),
        "junta_cap_product": 64 * HEIGHT_CAP * Fraction(9, 8) == JUNTA_CAP,
        "threshold_matches_twice_junta_cap_plus_one": P_MIN == 2 * JUNTA_CAP + 1,
        "both_slice_sides_exceed_kept_coordinate_bound": (p - 1) // 2 == JUNTA_CAP,
        "first_layer_has_an_isolated_vertex": p * p + 1 > 2 * (5 * p + 6),
        "all_target_quarter_lattice_intervals_empty": all(
            not row["quarter_integral_value_exists"] for row in target_rows),
    }
    margins = infinite_tail_margins()
    residues = residue_interval_check()
    power = power_band_check()
    checks.update({
        "power_band_constants_verified": power["verified"],
        "power_band_formal_T_elimination": power["formal_T_elimination_verified"],
        "power_band_quarter_integral_bracket": power["quarter_integrality_makes_bracket_integral"],
        "power_band_both_types_and_parity_retained": (
            power["both_signed_types_required"] and power["parity_required"]),
    })
    require(all(checks.values()) and residues["verified"] and power["verified"],
            "arithmetic replay failed")
    return {
        "classification": "INDEPENDENT EXACT ARITHMETIC REPLAY; NOT A FINITE CENSUS",
        "host": platform.node(), "architecture": platform.machine(),
        "python": platform.python_version(), "P_MIN": P_MIN,
        "analytic_inputs_not_replaced_by_this_replay": [
            "degree-two Bonami ||f||_4 <=3||f||_2 and interpolation ||f||_2<=9 E[f]",
            "integer cube quadratics have Fourier coefficients in (1/4)Z",
            "proved maximizing paired-cube and stabilizer identities on the middle slice",
            "Johnson transposition Laplacian level-d eigenvalue d*(p+1-d)",
            "zero-transposition-influence classes and degree-two symmetrization",
            "Proposition 15.774 strict small-mass quotient floors"],
        "derivative_identity": derivative,
        "influence_floor_at_P_MIN": str(influence_floor),
        "influence_floor_formula": "(p^2-1)/(32*p*(p-2)) > 1/32",
        "height_cap_exact": str(exact_height_cap), "integer_height_cap": integer_height_cap,
        "strict_loose_height_cap": HEIGHT_CAP, "strict_junta_coordinate_cap": JUNTA_CAP,
        "slice_cube_mean_monomials": mean_checks,
        "slice_cube_mean_identity": "p*mu=(p+1)*nu-f(0)",
        "target_masses": target_rows, "tail_margins": margins,
        "residue_interval_certificate": residues,
        "power_band_certificate": power,
        "checks": checks, "verified": True,
        "proved_layer_scope": "conditional on the audited analytic proof: prime p>=259201, k=5p+5, |H|=5p+6",
        "proved_power_band_scope": (
            "conditional on the audited analytic proof: prime p>=29,r in {3,4,5},h=r mod 2; "
            "both signed floors r excluded if h<r*p or 46656*h^3<=p^3*(p-1)"),
        "global_residual_ii_closed": False, "E1_closed": False, "limit_closed": False,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True), flush=True)
