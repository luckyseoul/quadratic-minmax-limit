#!/usr/bin/env python3
"""Prop. 15.681 -- stronger quadratic lift and the p=29 endpoint.

The paired-cube operator introduced in Proposition 15.680 applies to every
nonzero nonnegative integer-valued quadratic, not only to Boolean ones.  If
``B(X)=h>0`` and ``c=4p E[B]``, cube distance gives

    c >= p+1-4h.

Combining this with Proposition 15.642's stabilizer inequality gives

    c >= (p+1)/2,  p=3 mod 4,
    c >= (p-1)/2,  p=1 mod 4.

At the second all-finite endpoint this removes every positive-residue row
for ``p=29,31,37,41``.  For ``p=29,s=24`` only residue zero remains.  Pair
slack modulo four leaves either a 24-arc with at least four undetermined
directions or a set with one 3-secant and six undetermined directions.

In the first case, adjoining two undetermined infinity points gives a
26-arc in ``PG(2,29)``.  In the second, deleting one point of the unique
collinear triple and adjoining two infinity points gives a 25-arc.
Coolsaet--Sticker's exhaustive classification has respectively five and
ten projective classes of 26- and 25-arcs.  Exact PGL(2,29) orbit counts for
four- and five-point complements on a conic are also five and ten, so every
such arc is conic-contained.  Using three undetermined infinity points
then forces one conic to contain three collinear points, a contradiction.

This closes only the ``p=29,s=24`` endpoint (and gives a shorter independent
exclusion of the already closed ``p=37`` endpoint).  At this boundary
``p=31,41`` retain only residue zero; ``p=17,19,23``, later all-finite
sizes, the infinity-present remainder, residual (ii), R1, global QVAR,
Type I, and the limit remain open.
"""
from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15642 import stabilizer_mass_certificate
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15675 import first_even_survivor


ROOT = Path(__file__).resolve().parents[1]
P = 29
M = 15
S = 24
PERIOD = 30
PAIR_DEFICIT_BUDGET = S * (S - 1)


def paired_cube_integral_quadratic_floor(p: int) -> dict[str, object]:
    """Mass floor for a nonzero nonnegative integral quadratic on the slice.

    Put ``p=2m-1`` and fix ``X`` with ``B(X)=h>=1``.  Leave one point of
    ``X`` fixed, pair the other points of ``X`` with its complement, and
    average over the resulting Boolean cube.  On degree at most two the
    transition is ``T B=rho B+(1-rho)E[B]``, where ``rho=1/(p+1)``.

    The restriction is nonzero at ``X``.  The elementary degree-two cube
    distance bound gives support at least ``1/4``; integrality and
    nonnegativity therefore give cube average at least ``1/4``.  This is
    the first displayed inequality in the module docstring.  The second
    inequality is Proposition 15.642 applied at the same point ``X``.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("need odd p>=5")
    m = (p + 1) // 2
    rho = Fraction(1, p + 1)
    coordinate_mean = Fraction(m, p)
    pair_mean = Fraction(m * (m - 1), p * (p - 1))
    if (1 - rho) * coordinate_mean != Fraction(1, 2):
        raise ArithmeticError("paired-cube coordinate transition changed")
    if (1 - rho) * pair_mean != Fraction(1, 4):
        raise ArithmeticError("paired-cube pair transition changed")

    stabilizer_weight = Fraction(stabilizer_mass_certificate(p)["value"])
    stabilizer_scaled_coefficient = 4 * p * stabilizer_weight
    if p % 4 == 3:
        r = (p - 3) // 4
        expected_coefficient = Fraction(4)
        universal_scaled_floor = Fraction(p + 1, 2)
        convex_weight_on_stabilizer = Fraction(1, 2)
    else:
        r = (p - 1) // 4
        expected_coefficient = Fraction(4 * r, r + 1)
        universal_scaled_floor = Fraction(p - 1, 2)
        convex_weight_on_stabilizer = Fraction(r + 1, 2 * r + 1)
    if stabilizer_scaled_coefficient != expected_coefficient:
        raise ArithmeticError("stabilizer coefficient changed")

    # The weighted average cancels h exactly between
    # c>=alpha*h and c>=p+1-4h.
    lam = convex_weight_on_stabilizer
    h_coefficient = lam * stabilizer_scaled_coefficient - 4 * (1 - lam)
    combined_constant = (1 - lam) * (p + 1)
    if h_coefficient or combined_constant != universal_scaled_floor:
        raise ArithmeticError("paired-cube/stabilizer combination changed")

    return {
        "p": p,
        "middle_weight": m,
        "rho": rho,
        "degree_two_operator": "T(B)=rho*B+(1-rho)*E[B]",
        "cube_support_floor": Fraction(1, 4),
        "cube_average_floor": Fraction(1, 4),
        "point_value": "h=B(X)>=1",
        "paired_cube_scaled_inequality": "4p*E[B] >= p+1-4h",
        "stabilizer_scaled_coefficient_of_h": stabilizer_scaled_coefficient,
        "convex_weight_on_stabilizer_inequality": lam,
        "universal_scaled_mass_floor": universal_scaled_floor,
        "universal_mass_floor": universal_scaled_floor / (4 * p),
        "closed_form": (
            "(p+1)/2 for p=3 mod 4; (p-1)/2 for p=1 mod 4"
        ),
        "proved": True,
    }


def second_even_boundary(p: int) -> int:
    """Second even integer strictly above ``3(p-1)/4``."""
    if p < 19 or p % 2 == 0:
        raise ValueError("need odd p>=19")
    return first_even_survivor(p) + 2


def exact_type_rows(p: int, phase: int) -> list[dict[str, object]]:
    """Exact quotient/floor minima for one type at the second boundary."""
    if phase not in (0, 1):
        raise ValueError("phase must be zero or one")
    s = second_even_boundary(p)
    m = (p + 1) // 2
    period = p + 1
    rows: list[dict[str, object]] = []
    for u in range(m):
        target = m - u
        best_by_quotient: dict[int, tuple[int, int]] = {}
        for b in range(0, s + 1, 2):
            floor_value = full_symbolic_floor(p, b, phase)
            for quotient in range(target + 1):
                excess = 2 * u + period * quotient - floor_value
                if excess < 0 or excess == 2:
                    continue
                candidate = (s - b, b)
                old = best_by_quotient.get(quotient)
                if old is None or candidate < old:
                    best_by_quotient[quotient] = candidate

        states: dict[int, tuple[int, tuple[int, ...]]] = {0: (0, ())}
        for _ in range(m):
            next_states: dict[int, tuple[int, tuple[int, ...]]] = {}
            for used, (deficit, profile) in states.items():
                for quotient, (added, b) in best_by_quotient.items():
                    new_used = used + quotient
                    if new_used > target:
                        continue
                    candidate = (deficit + added, profile + (b,))
                    old = next_states.get(new_used)
                    if old is None or candidate < old:
                        next_states[new_used] = candidate
            states = next_states
        if target in states:
            deficit, profile = states[target]
            rows.append(
                {
                    "u": u,
                    "quotient_sum": target,
                    "minimum_deficit": deficit,
                    "profile": dict(sorted(Counter(profile).items())),
                }
            )
    return rows


def endpoint_residue_ledger(p: int) -> dict[str, object]:
    """Apply the new lift floor to pair-surviving endpoint rows."""
    if p not in (29, 31, 37, 41):
        raise ValueError("the checked endpoint ledger is for p=29,31,37,41")
    s = second_even_boundary(p)
    m = (p + 1) // 2
    budget = s * (s - 1)
    phase_zero = exact_type_rows(p, 0)
    phase_one = exact_type_rows(p, 1)
    pair_rows = []
    for zero in phase_zero:
        for one in phase_one:
            required = int(zero["minimum_deficit"]) + int(
                one["minimum_deficit"]
            )
            if required <= budget:
                pair_rows.append(
                    {
                        "u0": int(zero["u"]),
                        "u1": int(one["u"]),
                        "phase_zero_minimum_deficit": int(
                            zero["minimum_deficit"]
                        ),
                        "phase_one_minimum_deficit": int(
                            one["minimum_deficit"]
                        ),
                        "pair_slack": budget - required,
                        "phase_zero_profile": zero["profile"],
                        "phase_one_profile": one["profile"],
                    }
                )

    lift = paired_cube_integral_quadratic_floor(p)
    scaled_floor = Fraction(lift["universal_scaled_mass_floor"])
    positive_rows = []
    for row in pair_rows:
        u = int(row["u0"])
        if not u:
            continue
        scaled_mean = 2 * u
        least_positive_b_floor = min(
            full_symbolic_floor(p, b, 0) for b in range(2, s + 1, 2)
        )
        positive_rows.append(
            {
                "u0": u,
                "quotient_sum": m - u,
                "forces_quotient_zero": m - u < m,
                "zero_quotient_scaled_mean": scaled_mean,
                "least_positive_b_floor": least_positive_b_floor,
                "therefore_b_zero": scaled_mean < least_positive_b_floor,
                "factorization": "A=2B with B nonzero nonnegative integral quadratic",
                "new_scaled_lift_floor": scaled_floor,
                "excluded": scaled_floor > scaled_mean,
            }
        )
    residue_zero = [row for row in pair_rows if int(row["u0"]) == 0]
    return {
        "p": p,
        "s": s,
        "pair_deficit_budget": budget,
        "phase_zero_rows": phase_zero,
        "phase_one_rows": phase_one,
        "pair_survivors": pair_rows,
        "positive_residue_rows": positive_rows,
        "positive_residues_all_excluded": all(
            bool(row["excluded"]) for row in positive_rows
        ),
        "residue_zero_rows": residue_zero,
        "residue_zero_remains": bool(residue_zero),
        "new_lift": lift,
    }


def _histogram(profile: tuple[int, ...]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(Counter(profile).items())}


@lru_cache(maxsize=None)
def _p29_profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate every exact p=29 profile within a tight deficit cap."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<15")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor_value = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor_value
            if excess >= 0 and excess != 2:
                options.append((quotient, S - b, b))

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for _ in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                if new_used <= target and new_deficit <= deficit_cap:
                    next_states.add(
                        (new_used, new_deficit, tuple(sorted(profile + (b,))))
                    )
        states = next_states
    return tuple(
        sorted(
            (deficit, profile)
            for used, deficit, profile in states
            if used == target
        )
    )


def pair_slack_divisibility() -> dict[str, object]:
    """Exact local contribution to pair-budget slack."""
    return {
        "line_contribution": "2*(C(n,2)-floor(n/2))",
        "n=2r": "4*r*(r-1)",
        "n=2r+1": "4*r*r",
        "global_slack_modulus": 4,
        "slack_zero_iff": "every affine line has occupancy at most two",
        "slack_four_iff": "exactly one affine line has occupancy three and all others at most two",
        "proved": True,
    }


def p29_residue_zero_profiles() -> dict[str, object]:
    """Classify all p=29 residue-zero profiles allowed by pair slack."""
    phase_zero = _p29_profile_rows(0, 0, PAIR_DEFICIT_BUDGET - 308)
    phase_one = _p29_profile_rows(1, 14, PAIR_DEFICIT_BUDGET - 240)
    candidates = []
    for deficit_zero, profile_zero in phase_zero:
        for deficit_one, profile_one in phase_one:
            total = deficit_zero + deficit_one
            slack = PAIR_DEFICIT_BUDGET - total
            if slack < 0 or slack % 4:
                continue
            secants = Counter((S - b) // 2 for b in profile_zero + profile_one)
            candidates.append(
                {
                    "phase_profiles_b": {
                        "0": _histogram(profile_zero),
                        "1": _histogram(profile_one),
                    },
                    "phase_deficits": {
                        "0": deficit_zero,
                        "1": deficit_one,
                    },
                    "total_deficit": total,
                    "pair_slack": slack,
                    "arc": slack == 0,
                    "global_secant_distribution": {
                        str(key): value for key, value in sorted(secants.items())
                    },
                    "undetermined_directions": secants[0],
                }
            )

    global_shapes = {
        tuple(sorted(row["global_secant_distribution"].items()))
        for row in candidates
    }
    expected_shapes = {
        tuple(sorted({"0": 6, "11": 14, "12": 10}.items())),
        tuple(sorted({"0": 5, "2": 1, "11": 14, "12": 10}.items())),
        tuple(sorted({"0": 4, "1": 2, "11": 14, "12": 10}.items())),
    }
    if len(candidates) != 5 or global_shapes != expected_shapes:
        raise ArithmeticError("p=29 endpoint profile classification changed")
    near_arc = [row for row in candidates if int(row["pair_slack"]) == 4]
    arcs = [row for row in candidates if int(row["pair_slack"]) == 0]
    if len(near_arc) != 1 or len(arcs) != 4:
        raise ArithmeticError("p=29 arc/near-arc split changed")
    if int(near_arc[0]["undetermined_directions"]) != 6:
        raise ArithmeticError("p=29 near-arc missing-direction count changed")
    if min(int(row["undetermined_directions"]) for row in arcs) != 4:
        raise ArithmeticError("p=29 arc missing-direction floor changed")
    return {
        "p": P,
        "s": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "pair_slack_divisibility": pair_slack_divisibility(),
        "profiles": candidates,
        "distinct_global_shapes": [
            dict(shape) for shape in sorted(global_shapes)
        ],
        "near_arc_profile_count": len(near_arc),
        "near_arc_geometry": (
            "one 3-secant, all other lines of occupancy at most two, "
            "and six undetermined directions"
        ),
        "arc_profile_count": len(arcs),
        "arc_minimum_undetermined_directions": min(
            int(row["undetermined_directions"]) for row in arcs
        ),
        "proved": True,
    }


def _p1_generators(p: int) -> tuple[tuple[int, ...], ...]:
    """Translation, inversion, and nonsquare scaling on P^1(F_p)."""
    infinity = p
    primitive = next(
        value
        for value in range(2, p)
        if len({pow(value, exponent, p) for exponent in range(1, p)}) == p - 1
    )
    translation = tuple(
        infinity if x == infinity else (x + 1) % p for x in range(p + 1)
    )
    inversion = tuple(
        0 if x == infinity else infinity if x == 0 else pow(x, -1, p)
        for x in range(p + 1)
    )
    scaling = tuple(
        infinity if x == infinity else primitive * x % p
        for x in range(p + 1)
    )
    return translation, inversion, scaling


@lru_cache(maxsize=None)
def pgl2_group_order(p: int) -> int:
    """Generate PGL(2,p) exactly from three standard transformations."""
    generators = _p1_generators(p)
    identity = tuple(range(p + 1))
    seen = {identity}
    stack = [identity]
    while stack:
        permutation = stack.pop()
        for generator in generators:
            composite = tuple(generator[permutation[x]] for x in range(p + 1))
            if composite not in seen:
                seen.add(composite)
                stack.append(composite)
    return len(seen)


@lru_cache(maxsize=None)
def pgl2_subset_orbit_audit(p: int, subset_size: int) -> dict[str, object]:
    """Count PGL(2,p)-orbits of subsets of the projective line exactly."""
    if p != 29 or subset_size not in (4, 5):
        raise ValueError("this finite audit is scoped to p=29 and sizes 4,5")
    generators = _p1_generators(p)
    remaining = set(combinations(range(p + 1), subset_size))
    orbit_sizes = []
    while remaining:
        seed = next(iter(remaining))
        seen = {seed}
        stack = [seed]
        while stack:
            subset = stack.pop()
            for generator in generators:
                image = tuple(sorted(generator[x] for x in subset))
                if image not in seen:
                    seen.add(image)
                    stack.append(image)
        remaining.difference_update(seen)
        orbit_sizes.append(len(seen))
    expected_group_order = p * (p * p - 1)
    group_order = pgl2_group_order(p)
    if group_order != expected_group_order:
        raise ArithmeticError("standard generators did not produce PGL(2,p)")
    return {
        "p": p,
        "projective_line_size": p + 1,
        "subset_size": subset_size,
        "PGL2_order": group_order,
        "subset_count": comb(p + 1, subset_size),
        "orbit_count": len(orbit_sizes),
        "orbit_sizes": sorted(orbit_sizes),
        "proved": True,
    }


def p29_arc_classification_ledger() -> dict[str, object]:
    """Match all classified 25/26-arcs with conic-subset classes."""
    complements = {
        25: pgl2_subset_orbit_audit(P, 5),
        26: pgl2_subset_orbit_audit(P, 4),
    }
    classified_total = {25: 10, 26: 5}
    matches = {
        size: int(complements[size]["orbit_count"]) == total
        for size, total in classified_total.items()
    }
    if not all(matches.values()):
        raise ArithmeticError("conic-subset orbit counts no longer exhaust arcs")
    return {
        "external_dependency": True,
        "source": (
            "K. Coolsaet and H. Sticker, The complete k-arcs of "
            "PG(2,27) and PG(2,29), J. Combin. Des. 19 (2011), 111-130"
        ),
        "doi": "10.1002/jcd.20261",
        "open_pdf": (
            "https://backoffice.biblio.ugent.be/download/1247338/1247417"
        ),
        "location": "Table 5 (arcs not necessarily complete)",
        "PGL_equals_PGammaL": "true because 29 is prime",
        "classified_projective_arc_classes": classified_total,
        "conic_complement_PGL2_orbits": {
            size: complements[size] for size in (25, 26)
        },
        "uniqueness_of_containing_conic": (
            "five arc points determine a unique nondegenerate conic; hence "
            "PGL3 classes of conic subsets equal PGL2 complement orbits"
        ),
        "all_25_and_26_arcs_conic_contained": all(matches.values()),
        "proved_conditional_on_external_classification": True,
    }


def p29_geometric_exclusion() -> dict[str, object]:
    """Exclude every residue-zero p=29 profile by arc extension."""
    profiles = p29_residue_zero_profiles()
    classification = p29_arc_classification_ledger()
    proved = bool(
        profiles["proved"]
        and int(profiles["arc_minimum_undetermined_directions"]) >= 3
        and int(profiles["near_arc_profile_count"]) == 1
        and classification["all_25_and_26_arcs_conic_contained"] is True
    )
    if not proved:
        raise ArithmeticError("p=29 geometric exclusion audit failed")
    return {
        "arc_case": {
            "starting_size": 24,
            "minimum_undetermined_infinity_points": 4,
            "adjoin_two_size": 26,
            "classification_consequence": "contained in a conic",
        },
        "near_arc_case": {
            "starting_size": 24,
            "unique_3_secant": True,
            "delete_one_point_from_triple_size": 23,
            "result_is_arc": True,
            "undetermined_infinity_points_preserved": 6,
            "adjoin_two_size": 25,
            "classification_consequence": "contained in a conic",
        },
        "common_contradiction": (
            "choose three undetermined infinity points; the conics through "
            "the two corresponding extensions share at least 23 arc points "
            "and hence coincide, forcing one nondegenerate conic to contain "
            "three collinear points"
        ),
        "classification": classification,
        "excluded": True,
    }


def theorem_record() -> dict[str, object]:
    endpoint_ledgers = {
        str(p): endpoint_residue_ledger(p) for p in (29, 31, 37, 41)
    }
    expected_survivors = {
        "29": [0, 2, 3, 4, 5],
        "31": [0, 2, 3, 4, 5, 6],
        "37": [2, 3, 4, 5],
        "41": [0, 2, 3, 4, 5, 6, 7],
    }
    observed_survivors = {
        key: [int(row["u0"]) for row in value["pair_survivors"]]
        for key, value in endpoint_ledgers.items()
    }
    if observed_survivors != expected_survivors:
        raise ArithmeticError("small endpoint residue ledger changed")
    p29_geometry = p29_geometric_exclusion()
    proved = bool(
        all(
            row["positive_residues_all_excluded"] is True
            for row in endpoint_ledgers.values()
        )
        and endpoint_ledgers["29"]["residue_zero_remains"] is True
        and p29_geometry["excluded"] is True
    )
    return {
        "prop": "15.681",
        "title": "Paired-cube integral lift and the p=29 next endpoint",
        "proved": proved,
        "theorem": {
            "all_odd_p_at_least_5": (
                "every nonzero nonnegative integer-valued quadratic B on "
                "J(p,(p+1)/2) has 4p E[B] >= (p+1)/2 for p=3 mod 4 "
                "and >=(p-1)/2 for p=1 mod 4"
            ),
            "p29_s24_next_all_finite_endpoint": "EXCLUDED",
            "p31_p41_same_boundary": "ONLY_RESIDUE_ZERO_REMAINS",
            "p37_same_boundary": "ALREADY_CLOSED_BY_15.680_AND_REPROVED",
            "remaining_smaller_endpoints": [17, 19, 23, 31, 41],
            "later_all_finite_sizes": "OPEN",
            "infinity_present_remainder": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "paired_cube_integral_lift_samples": {
            str(p): paired_cube_integral_quadratic_floor(p)
            for p in (5, 7, 11, 17, 29, 31, 37, 41, 101)
        },
        "endpoint_ledgers": endpoint_ledgers,
        "p29_profiles": p29_residue_zero_profiles(),
        "p29_geometric_exclusion": p29_geometry,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.681 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15681.json"
    destination.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
    print("Prop 15.681 p=29,s=24 next all-finite endpoint: excluded")
    print("  p=31,41: only residue zero remains")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
