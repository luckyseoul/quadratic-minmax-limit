#!/usr/bin/env python3
"""Prop. 15.725 -- parabola-plus-internal boundary family exclusion.

This proposition closes one explicit strict-deficit family in the first open
``|D|=p+1`` boundary shell.  In trace-zero coordinates

    F_(p^2) = F_p(omega),             omega^2 = nu, chi_p(nu) = -1,

the all-finite boundary is

    D = {x+x^2*omega : x in F_p} union {a*omega},
    chi_p(-a) = -1.

Orient its product sign so that the exceptional vertical direction is phase
one, and send ``z=a*omega`` to infinity by ``g(u)=1/(u-z)``.  The secant
involution proves that the transported product sign is ``+1``.  Explicitly,
for ``h(x)=x+(x^2-a)*omega`` and ``x!=0``,

    h(-a/x)=(-a/x^2)h(x).

The involution ``x -> -a/x`` has no fixed point because ``chi_p(-a)=-1``,
and every nonzero base-field scalar is a square in the quadratic extension.
The nonzero extension characters therefore multiply to one, while the
``x=0`` factor is ``chi_p(-nu)``.  This is also the source vertical type, so
the transported sign is one.  Because the transformed boundary contains
infinity and has ``p`` finite points, the two direction-type factors and the
two odd parity factors in the exact parity sign cancel.  Hence every
transformed direction has phase zero.  The finite points are

    (A_x,B_x) = (x/Q(x), -(x^2-a)/Q(x)),
    Q(x) = x^2-nu*(x^2-a)^2.

For a projective direction represented by the fibre label ``r*A+s*B``, its
actual Paley quadratic type is

    eps(r,s) = chi_p(s^2-nu*r^2).

The exact Proposition 15.669 phase-zero floor is summed separately over the
``(p+1)/2`` directions of each type and compared with the exact type budget
``(p+1)^2/2``.

There are two disjoint proof branches.

1. The symbolic character-curve bounds prove the family impossible for every
   prime ``p>=53``.  This module records the exact inequalities and their
   monotone integer threshold checks; it does not replace them by a finite
   prime scan.

2. The remaining primes

       17,19,23,29,31,37,41,43,47

   are exhausted with exact modular arithmetic.  Every nonsquare ``nu``, every
   admissible ``a``, and all ``p+1`` directions are evaluated directly.  A
   survivor raises ``ParabolaFamilySurvivor`` immediately with its complete
   typed profile.

The finite enumeration also checks the normalization

    nu' = lambda^2*nu  ==>  (nu',a) ~ (nu,lambda^2*a),
    [r:s] -> [lambda*r:s].

Indeed, the field isomorphism ``omega' -> lambda*omega``, followed by source
scaling by ``lambda``, sends the first boundary to the second; after inversion
the target coordinates are ``(A,B)->(A/lambda,B)``.  This preserves norm
characters and preserves fibre occupancies after the explicit scalar
permutation of their labels.  Direct enumeration is nevertheless retained for
every ``nu``; the reduction is an independent consistency check.

Only this explicit boundary family is closed.  The whole ``p+1`` shell,
residual (ii), Type I, and the limit remain open.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from typing import Iterable

from e1_gmin_m4_prop15669 import full_symbolic_floor


FINITE_PRIMES = (17, 19, 23, 29, 31, 37, 41, 43, 47)


class ParabolaFamilySurvivor(RuntimeError):
    """Raised immediately if an exact finite case meets both type budgets."""

    def __init__(self, record: dict[str, object]):
        self.record = record
        super().__init__(
            "parabola-plus-internal survivor: "
            + json.dumps(record, sort_keys=True, separators=(",", ":"))
        )


def legendre(value: int, p: int) -> int:
    """Return the quadratic character of ``value`` in ``F_p``."""
    value %= p
    if value == 0:
        return 0
    return 1 if pow(value, (p - 1) // 2, p) == 1 else -1


def nonsquares(p: int) -> tuple[int, ...]:
    """All nonzero nonsquares in deterministic increasing order."""
    return tuple(value for value in range(1, p) if legendre(value, p) == -1)


def admissible_internal_parameters(p: int) -> tuple[int, ...]:
    """All nonzero ``a`` for which ``(0,a)`` is internal to ``y=x^2``."""
    return tuple(value for value in range(1, p) if legendre(-value, p) == -1)


def projective_directions(p: int) -> tuple[tuple[int, int], ...]:
    """Canonical labels for ``P^1(F_p)``: ``(1,s)`` and ``(0,1)``."""
    return tuple((1, slope) for slope in range(p)) + ((0, 1),)


def normalize_direction(r: int, s: int, p: int) -> tuple[int, int]:
    """Normalize a nonzero projective pair to ``projective_directions`` form."""
    r %= p
    s %= p
    if r:
        return 1, s * pow(r, p - 2, p) % p
    if not s:
        raise ValueError("zero pair is not a projective direction")
    return 0, 1


def extension_character(real: int, imag: int, nu: int, p: int) -> int:
    """Quadratic character in ``F_p(omega)``, ``omega^2=nu``."""
    return legendre(real * real - nu * imag * imag, p)


def direction_type(p: int, nu: int, direction: tuple[int, int]) -> int:
    """Actual Paley type of the kernel of the fibre label ``r*A+s*B``."""
    r, s = direction
    eps = legendre(s * s - nu * r * r, p)
    if eps not in (-1, 1):
        raise ArithmeticError("anisotropic norm vanished on a direction")
    return eps


def transformed_points(
    p: int, nu: int, a: int
) -> tuple[tuple[tuple[int, int], ...], int, int]:
    """Return transformed points, boundary multiplier, and vertical type.

    The calculation uses only exact arithmetic in ``F_p``.  The returned
    multiplier is ``product_x chi_(p^2)(u_x-z)``.
    """
    if legendre(nu, p) != -1 or legendre(-a, p) != -1:
        raise ValueError("need nonsquare nu and chi(-a)=-1")
    points: list[tuple[int, int]] = []
    multiplier = 1
    for x in range(p):
        y = (x * x - a) % p
        denominator = (x * x - nu * y * y) % p
        if denominator == 0:
            raise ArithmeticError("internal point entered the parabola")
        inverse = pow(denominator, p - 2, p)
        points.append((x * inverse % p, -y * inverse % p))
        multiplier *= legendre(denominator, p)
    vertical_type = legendre(-nu, p)
    return tuple(points), multiplier, vertical_type


def direction_record(
    p: int,
    nu: int,
    points: Iterable[tuple[int, int]],
    direction: tuple[int, int],
) -> dict[str, object]:
    """Compute one exact occupancy, odd-fibre count, type, and floor."""
    r, s = direction
    occupancy = [0] * p
    for real, imag in points:
        occupancy[(r * real + s * imag) % p] += 1
    b = sum(value & 1 for value in occupancy)
    eps = direction_type(p, nu, direction)
    floor = full_symbolic_floor(p, b, 0)
    return {
        "direction": direction,
        "type": eps,
        "b": b,
        "floor": floor,
        "occupancy": tuple(occupancy),
    }


def _typed_profile(
    rows: Iterable[dict[str, object]],
) -> tuple[tuple[int, int, int], ...]:
    counts = Counter((int(row["type"]), int(row["b"])) for row in rows)
    return tuple((eps, b, count) for (eps, b), count in sorted(counts.items()))


def _profile_json(profile: tuple[tuple[int, int, int], ...]) -> list[list[int]]:
    return [[eps, b, count] for eps, b, count in profile]


def exact_case(p: int, nu: int, a: int) -> dict[str, object]:
    """Evaluate one ``(p,nu,a)`` case in all actual direction types."""
    points, multiplier, vertical_type = transformed_points(p, nu, a)
    rows = tuple(
        direction_record(p, nu, points, direction)
        for direction in projective_directions(p)
    )
    type_counts = Counter(int(row["type"]) for row in rows)
    m = (p + 1) // 2
    expected_type_counts = Counter({-1: m, 1: m})
    if type_counts != expected_type_counts:
        raise ArithmeticError(
            f"quadratic type split failed at {(p, nu, a)}: {type_counts}"
        )

    # Phase-one on the original exceptional vertical direction means
    # c_H=vertical_type.  The exact boundary multiplier must be the same sign.
    source_c_h = vertical_type
    transported_c_h = source_c_h * multiplier
    if multiplier != vertical_type or transported_c_h != 1:
        raise ArithmeticError(
            f"phase transport failed at {(p, nu, a, multiplier, vertical_type)}"
        )

    type_sums = {
        eps: sum(int(row["floor"]) for row in rows if row["type"] == eps)
        for eps in (-1, 1)
    }
    budget = m * (p + 1)
    profile = _typed_profile(rows)
    return {
        "p": p,
        "nu": nu,
        "a": a,
        "source_vertical_type": vertical_type,
        "boundary_multiplier": multiplier,
        "transported_c_H": transported_c_h,
        "common_phase": 0,
        "type_counts": {-1: type_counts[-1], 1: type_counts[1]},
        "type_sums": type_sums,
        "type_budget": budget,
        "typed_profile": profile,
        "b1_directions": sum(int(row["b"]) == 1 for row in rows),
        "excluded": type_sums[-1] > budget or type_sums[1] > budget,
        "both_types_strict": type_sums[-1] > budget and type_sums[1] > budget,
        "rows": rows,
    }


def normalization_to_reference(p: int, nu: int, a: int) -> dict[str, int]:
    """Map ``(nu,a)`` to the least nonsquare representative exactly."""
    reference_nu = nonsquares(p)[0]
    ratio = nu * pow(reference_nu, p - 2, p) % p
    lambdas = [value for value in range(1, p) if value * value % p == ratio]
    if len(lambdas) != 2:
        raise ArithmeticError("nonsquare ratio did not have two square roots")
    scale = min(lambdas)
    reference_a = scale * scale * a % p
    if legendre(-reference_a, p) != -1:
        raise ArithmeticError("normalization left the internal-point orbit")
    return {
        "reference_nu": reference_nu,
        "scale": scale,
        "reference_a": reference_a,
    }


def normalized_target_direction(
    p: int, scale: int, direction: tuple[int, int]
) -> tuple[int, int]:
    """Transport ``[r:s]`` under ``(A,B)->(A/scale,B)``."""
    r, s = direction
    return normalize_direction(scale * r, s, p)


def normalized_fibre_label_scale(
    p: int, scale: int, direction: tuple[int, int]
) -> int:
    """Return ``k`` with direct label equal to ``k`` times canonical label.

    In reference coordinates the unnormalized direction is
    ``(scale*r,s)``.  Canonical projective normalization divides its fibre
    label by its first nonzero coordinate, so the direct occupancy at label
    ``t`` equals the reference occupancy at ``t/k``.
    """
    r, s = direction
    first = scale * r % p
    label_scale = first if first else s % p
    if label_scale == 0:
        raise ValueError("zero pair is not a projective direction")
    return label_scale


def _stable_update(digest, record: object) -> None:
    digest.update(
        (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode()
    )


def _witness(case: dict[str, object]) -> dict[str, object]:
    sums = case["type_sums"]
    assert isinstance(sums, dict)
    return {
        "nu": int(case["nu"]),
        "a": int(case["a"]),
        "type_sum_minus": int(sums[-1]),
        "type_sum_plus": int(sums[1]),
        "typed_profile": _profile_json(case["typed_profile"]),  # type: ignore[arg-type]
    }


def symbolic_phase_transport() -> dict[str, object]:
    """Return the exact all-prime product-sign identity used by both branches."""
    return {
        "h": "h(x)=x+(x^2-a)*omega",
        "nonzero_involution": "iota(x)=-a/x",
        "involution_identity": "h(iota(x))=(-a/x^2)h(x)",
        "fixed_point_equation": "x^2=-a",
        "fixed_point_free_hypothesis": "chi_p(-a)=-1",
        "base_scalars_square_in_extension": True,
        "nonzero_factor_product": 1,
        "x_zero_norm": "N(-a*omega)=-nu*a^2",
        "boundary_multiplier": "chi_p(-nu)",
        "source_vertical_type": "chi_p(-nu)",
        "transported_c_H": 1,
        "directional_parity_sign": (
            "eps*(-1)^((|H|-3)/2)*c_H*eps^[infinity in D]*(-1)^b"
        ),
        "residual_edge_count": "|H|=4p+1",
        "residual_exponent": "(|H|-3)/2=2p-1 is odd",
        "infinity_in_transformed_boundary": True,
        "finite_boundary_size": "p, so every b is odd",
        "sign_reduction": "eps*(-1)*1*eps*(-1)=1",
        "phase_independent_of_direction_type": True,
        "common_phase": 0,
        "proved": True,
    }


def audit_prime(p: int) -> dict[str, object]:
    """Exhaust one finite prime and return compact exact certificate data."""
    if p not in FINITE_PRIMES:
        raise ValueError(f"finite certificate is restricted to {FINITE_PRIMES}")
    nu_values = nonsquares(p)
    a_values = admissible_internal_parameters(p)
    directions = projective_directions(p)
    m = (p + 1) // 2
    budget = m * (p + 1)

    direction_digest = sha256()
    case_digest = sha256()
    normalization_digest = sha256()
    canonical_cache: dict[int, dict[str, object]] = {}
    b_histogram: Counter[int] = Counter()
    b1_per_case: Counter[int] = Counter()
    floor_sum_pairs: Counter[tuple[int, int]] = Counter()
    typed_profiles: Counter[tuple[tuple[int, int, int], ...]] = Counter()
    both_types_strict = 0
    phase_checks = 0
    symmetry_checks = 0
    minimum_minus: tuple[int, dict[str, object]] | None = None
    minimum_plus: tuple[int, dict[str, object]] | None = None
    minimum_joint: tuple[int, dict[str, object]] | None = None

    for nu in nu_values:
        for a in a_values:
            case = exact_case(p, nu, a)
            phase_checks += 1
            sums = case["type_sums"]
            assert isinstance(sums, dict)
            gap_minus = int(sums[-1]) - budget
            gap_plus = int(sums[1]) - budget
            joint_gap = max(gap_minus, gap_plus)

            if gap_minus <= 0 and gap_plus <= 0:
                survivor = _witness(case)
                survivor["type_budget"] = budget
                survivor["all_directions"] = [
                    {
                        "direction": list(row["direction"]),
                        "type": int(row["type"]),
                        "b": int(row["b"]),
                        "floor": int(row["floor"]),
                        "occupancy": list(row["occupancy"]),
                    }
                    for row in case["rows"]  # type: ignore[union-attr]
                ]
                raise ParabolaFamilySurvivor(survivor)

            if bool(case["both_types_strict"]):
                both_types_strict += 1
            b1_per_case[int(case["b1_directions"])] += 1
            profile = case["typed_profile"]
            assert isinstance(profile, tuple)
            typed_profiles[profile] += 1
            floor_sum_pairs[(int(sums[-1]), int(sums[1]))] += 1

            witness = _witness(case)
            if minimum_minus is None or gap_minus < minimum_minus[0]:
                minimum_minus = (gap_minus, witness)
            if minimum_plus is None or gap_plus < minimum_plus[0]:
                minimum_plus = (gap_plus, witness)
            if minimum_joint is None or joint_gap < minimum_joint[0]:
                minimum_joint = (joint_gap, witness)

            normalization = normalization_to_reference(p, nu, a)
            reference_a = normalization["reference_a"]
            if reference_a not in canonical_cache:
                canonical_cache[reference_a] = exact_case(
                    p, normalization["reference_nu"], reference_a
                )
            reference = canonical_cache[reference_a]
            reference_rows = {
                tuple(row["direction"]): row
                for row in reference["rows"]  # type: ignore[union-attr]
            }

            _stable_update(
                normalization_digest,
                [p, nu, a, normalization["scale"], reference_a],
            )
            for row in case["rows"]:  # type: ignore[union-attr]
                direction = tuple(row["direction"])
                mapped = normalized_target_direction(
                    p, normalization["scale"], direction  # type: ignore[arg-type]
                )
                label_scale = normalized_fibre_label_scale(
                    p, normalization["scale"], direction  # type: ignore[arg-type]
                )
                reference_row = reference_rows[mapped]
                direct_signature = (
                    int(row["type"]),
                    int(row["b"]),
                    int(row["floor"]),
                    tuple(row["occupancy"]),
                )
                reference_signature = (
                    int(reference_row["type"]),
                    int(reference_row["b"]),
                    int(reference_row["floor"]),
                    tuple(
                        reference_row["occupancy"][
                            (
                                label
                                * pow(label_scale, p - 2, p)
                            )
                            % p
                        ]
                        for label in range(p)
                    ),
                )
                if direct_signature != reference_signature:
                    raise ArithmeticError(
                        "nonsquare normalization mismatch at "
                        f"{(p, nu, a, direction, mapped)}"
                    )
                symmetry_checks += 1
                _stable_update(
                    normalization_digest,
                    [
                        p,
                        nu,
                        a,
                        int(direction[0]),
                        int(direction[1]),
                        int(mapped[0]),
                        int(mapped[1]),
                        label_scale,
                    ],
                )
                b_histogram[int(row["b"])] += 1
                _stable_update(
                    direction_digest,
                    [
                        p,
                        nu,
                        a,
                        int(direction[0]),
                        int(direction[1]),
                        int(row["type"]),
                        int(row["b"]),
                        int(row["floor"]),
                        list(row["occupancy"]),
                    ],
                )

            _stable_update(
                case_digest,
                {
                    "p": p,
                    "nu": nu,
                    "a": a,
                    "source_vertical_type": case["source_vertical_type"],
                    "boundary_multiplier": case["boundary_multiplier"],
                    "transported_c_H": case["transported_c_H"],
                    "type_sum_minus": int(sums[-1]),
                    "type_sum_plus": int(sums[1]),
                    "typed_profile": _profile_json(profile),
                },
            )

    parameter_cases = len(nu_values) * len(a_values)
    direction_records = parameter_cases * len(directions)
    if minimum_minus is None or minimum_plus is None or minimum_joint is None:
        raise ArithmeticError("finite audit unexpectedly had no cases")
    if phase_checks != parameter_cases or symmetry_checks != direction_records:
        raise ArithmeticError("audit accounting mismatch")

    return {
        "p": p,
        "nonsquare_count": len(nu_values),
        "admissible_a_count": len(a_values),
        "parameter_case_count": parameter_cases,
        "direction_record_count": direction_records,
        "phase_transport_check_count": phase_checks,
        "normalization_direction_check_count": symmetry_checks,
        "canonical_nu": nu_values[0],
        "canonical_a_case_count": len(canonical_cache),
        "type_budget": budget,
        "both_types_strict_case_count": both_types_strict,
        "survivor_count": 0,
        "minimum_type_minus_gap": minimum_minus[0],
        "minimum_type_plus_gap": minimum_plus[0],
        "minimum_joint_exclusion_gap": minimum_joint[0],
        "minimum_type_minus_witness": minimum_minus[1],
        "minimum_type_plus_witness": minimum_plus[1],
        "minimum_joint_witness": minimum_joint[1],
        "distinct_type_sum_pairs": len(floor_sum_pairs),
        "distinct_typed_profiles": len(typed_profiles),
        "type_sum_pair_histogram": {
            f"{minus},{plus}": count
            for (minus, plus), count in sorted(floor_sum_pairs.items())
        },
        "typed_profile_count_histogram": {
            json.dumps(_profile_json(profile), separators=(",", ":")): count
            for profile, count in sorted(typed_profiles.items())
        },
        "b_histogram": {str(b): count for b, count in sorted(b_histogram.items())},
        "b1_directions_per_case": {
            str(count): cases for count, cases in sorted(b1_per_case.items())
        },
        "all_cases_have_one_b1_direction": b1_per_case == Counter({1: parameter_cases}),
        "all_cases_excluded_in_both_types": both_types_strict == parameter_cases,
        "direction_audit_sha256": direction_digest.hexdigest(),
        "case_audit_sha256": case_digest.hexdigest(),
        "normalization_audit_sha256": normalization_digest.hexdigest(),
        "proved": both_types_strict == parameter_cases,
    }


def large_prime_symbolic_branch() -> dict[str, object]:
    """Record the separate all-``p>=53`` character-bound proof.

    The four radical inequalities are converted to integer polynomial checks
    at the threshold.  Their forward differences are positive from 53 onward.
    """
    p = 53
    checks = {
        # p-6*sqrt(p)-8>0, after p-8>0 and squaring.
        "c_nonzero_lower_gt_1": (p - 8) ** 2 - 36 * p,
        # 2p-6*sqrt(p)-28>0, equivalently p-14>3*sqrt(p).
        "c_nonzero_upper_lt_p_minus_2": (p - 14) ** 2 - 9 * p,
        # p-4*sqrt(p)-5>0.
        "c_zero_lower_gt_1": (p - 5) ** 2 - 16 * p,
        # p-2*sqrt(p)-25>0.
        "c_zero_upper_lt_p_minus_2": (p - 25) ** 2 - 4 * p,
    }
    forward_differences_at_53 = {
        "c_nonzero_lower_gt_1": 2 * p - 51,
        "c_nonzero_upper_lt_p_minus_2": 2 * p - 36,
        "c_zero_lower_gt_1": 2 * p - 25,
        "c_zero_upper_lt_p_minus_2": 2 * p - 53,
    }
    squared_polynomials = {
        "c_nonzero_lower_gt_1": "(p-8)^2-36p=p^2-52p+64",
        "c_nonzero_upper_lt_p_minus_2": "(p-14)^2-9p=p^2-37p+196",
        "c_zero_lower_gt_1": "(p-5)^2-16p=p^2-26p+25",
        "c_zero_upper_lt_p_minus_2": "(p-25)^2-4p=p^2-54p+625",
    }
    forward_difference_formulas = {
        "c_nonzero_lower_gt_1": "2p-51",
        "c_nonzero_upper_lt_p_minus_2": "2p-36",
        "c_zero_lower_gt_1": "2p-25",
        "c_zero_upper_lt_p_minus_2": "2p-53",
    }
    proved = all(value > 0 for value in checks.values()) and all(
        value > 0 for value in forward_differences_at_53.values()
    )
    if not proved:
        raise ArithmeticError("large-prime threshold ledger failed")
    return {
        "scope": "every prime p>=53",
        "separate_from_finite_enumeration": True,
        "character_curve_bounds": {
            "c_nonzero": [
                "(p-6*sqrt(p)-5)/3 <= b_c",
                "b_c <= (p+22+6*sqrt(p))/3",
            ],
            "c_zero": [
                "(p-2-4*sqrt(p))/3 <= b_0",
                "b_0 <= (2p+2*sqrt(p)+19)/3",
            ],
        },
        "odd_integrality_consequence": "3 <= b_c <= p-4 for every finite c",
        "oddness_reason": "sum of the p fibre occupancies is odd",
        "exceptional_direction": {
            "direction": "[0:1]",
            "reason": "B(x) is even in x, leaving only x=0 unpaired",
            "b": 1,
            "type": 1,
        },
        "phase_transport": symbolic_phase_transport(),
        "phase": 0,
        "phase_zero_floor_bounds": {
            "b=1": "p+1",
            "3<=b<=p-4": "at least 2p-6",
        },
        "type_budget": "m(p+1), m=(p+1)/2",
        "type_minus_exact_lower_gap": "m(p-7)>0",
        "type_plus_exact_lower_gap": "(m-1)(p-7)>0",
        "positive_squaring_side_conditions": [
            "p-8>0",
            "p-14>0",
            "p-5>0",
            "p-25>0",
        ],
        "threshold_squared_polynomials": squared_polynomials,
        "threshold_squared_polynomial_values": checks,
        "forward_difference_formulas": forward_difference_formulas,
        "forward_differences_at_53": forward_differences_at_53,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def finite_exact_certificate() -> dict[str, object]:
    """Exhaust all nine remaining primes; cached for tests and theorem use."""
    rows = {str(p): audit_prime(p) for p in FINITE_PRIMES}
    parameter_cases = sum(int(row["parameter_case_count"]) for row in rows.values())
    direction_records = sum(
        int(row["direction_record_count"]) for row in rows.values()
    )
    symmetry_reduced_cases = sum(
        int(row["canonical_a_case_count"]) for row in rows.values()
    )
    proved = all(bool(row["proved"]) for row in rows.values())
    return {
        "scope": list(FINITE_PRIMES),
        "separate_from_large_prime_proof": True,
        "all_nonsquare_nu_enumerated_directly": True,
        "all_a_with_chi_minus_a_negative_enumerated": True,
        "all_projective_directions_enumerated": True,
        "actual_norm_character_types_used": True,
        "common_phase_zero_checked_per_case": True,
        "parameter_case_count": parameter_cases,
        "symmetry_reduced_case_count": symmetry_reduced_cases,
        "direction_record_count": direction_records,
        "survivor_count": 0,
        "rows": rows,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def theorem_parabola_internal_family() -> dict[str, object]:
    """Combine the disjoint large-prime proof and finite certificates."""
    large = large_prime_symbolic_branch()
    finite = finite_exact_certificate()
    proved = bool(large["proved"] and finite["proved"])
    return {
        "prop": "15.725",
        "title": "Parabola-plus-internal boundary family exclusion",
        "proved": proved,
        "scope": "every odd prime p>=17",
        "boundary_family": (
            "{x+x^2*omega:x in F_p} union {a*omega}, "
            "omega^2=nu nonsquare, chi(-a)=-1"
        ),
        "orientation": "exceptional original vertical direction has phase one",
        "normalization": "send a*omega to infinity by u -> 1/(u-a*omega)",
        "transported_common_phase": 0,
        "large_prime_branch": large,
        "finite_exact_branch": finite,
        "theorem": {
            "parabola_plus_internal_family": "EXCLUDED",
            "whole_p_plus_one_shell": "OPEN",
            "residual_ii": False,
            "type_I": False,
            "limit_exists": False,
        },
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    theorem = theorem_parabola_internal_family()
    if theorem["proved"] is not True:
        raise ArithmeticError("Proposition 15.725 audit failed")
    finite = theorem["finite_exact_branch"]
    assert isinstance(finite, dict)
    print("Prop 15.725 parabola-plus-internal family: excluded for p>=17")
    print(
        "  finite exact cases="
        f"{finite['parameter_case_count']} "
        f"directions={finite['direction_record_count']}"
    )
    for p, row in finite["rows"].items():  # type: ignore[union-attr]
        print(
            f"  p={p} cases={row['parameter_case_count']} "
            f"directions={row['direction_record_count']} "
            f"hash={row['direction_audit_sha256']}"
        )
    return theorem


if __name__ == "__main__":
    main()
