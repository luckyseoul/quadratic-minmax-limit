#!/usr/bin/env python3
r"""Prop. 15.746 -- support-330 equality catalog and the ``u=4`` reduction.

This proposition has two deliberately separate parts.

First, positive quadrature makes the phase-one ``b=2`` equality baseline
pointwise on every intersection layer.  The residual lift is therefore a
globally nonnegative integral quadratic of scaled mass ten.  Proposition
15.688 makes it Boolean, of support 330 on ``J(13,7)``.  The exact classifier
in ``scripts/p13_support330_boolean_classifier.py`` proves that the 78
omitted-pair and 286 all-equal-triple supports exhaust this equality case.

Second, the two catalog families have coefficient offsets three and five.
The common p13 t4 ledger therefore has uniform hard parallel count ``P=3``
or ``P=5``.  At least two opposite cells have mean twelve.  In the ``P=3``
branch they are ``Q=5,b=0`` mass-twelve lifts, and the seven hard directions
force a common binary sextic to vanish.  This is an open reduction, not a
closure of ``u=4`` or ``u=6``.
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.p13_support330_boolean_classifier import (  # noqa: E402
    build_classifier_model,
    catalog_arithmetic_certificate,
)

from e1_gmin_m4_prop15652 import parity_floor_certificate  # noqa: E402
from e1_gmin_m4_prop15688 import (  # noqa: E402
    sharp_integral_quadratic_lift_floor,
)
from e1_gmin_m4_prop15734 import (  # noqa: E402
    BRANCH_B2,
    baseline_coefficient_rules,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15744 import t4_all_residue_sieve  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402


P = 13
Q_MODULUS = 6
DOMAIN_SIZE = 1716
H_EDGE_COUNT = 61
HARD_DIRECTION_COUNT = 7
SUPPORT330 = 330
SUPPORT396 = 396
EXPECTED_MODEL_SHA256 = (
    "7ace4efab04ada945381ad7826a8659dbd9c6d2292e543c6dced0894750da066"
)
EXPECTED_IDENTITY_SHA256 = (
    "ee92d6662f0f14523dc4c6620f89b407a66048dd4a6c0962dd9b058800136083"
)
EXPECTED_CATALOG_SHA256 = (
    "4edf1fe1b9c73f05598b667dba121f064807c68421a4df2c8db7090a3e3ff35f"
)
EXPECTED_ANCHORED_CATALOG_SHA256 = (
    "84ce6099dcca66f7cc2792dc60bcbb378672f2e9cac2b19e02812f2f20563c7a"
)
CLASSIFIER_EVIDENCE = ROOT / "evidence" / "p13_support330_boolean_classifier.json"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _mod_rank(rows: list[tuple[int, ...]], modulus: int) -> int:
    work = [[value % modulus for value in row] for row in rows]
    if not work:
        return 0
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(work))
                if work[row][column] % modulus
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column] % modulus, -1, modulus)
        work[rank] = [(value * inverse) % modulus for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column] % modulus:
                continue
            scalar = work[row][column] % modulus
            work[row] = [
                (left - scalar * right) % modulus
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == len(work[0]):
            break
    return rank


@lru_cache(maxsize=1)
def mass10_boolean_lift_bridge() -> dict[str, object]:
    """Certify the pointwise XNOR-to-Boolean bridge at scaled mean 22."""
    quadrature = parity_floor_certificate(P, 2, 1)
    rules = baseline_coefficient_rules(P)
    lift = sharp_integral_quadratic_lift_floor(P)

    nodes = tuple(int(node) for node in quadrature["quadrature_nodes"])
    weights = tuple(quadrature["quadrature_weights"])
    all_intersection_layers = tuple(range(3))
    baseline_values = tuple((1 - intersection) ** 2 for intersection in nodes)
    parity_values = tuple((intersection + 1) % 2 for intersection in nodes)
    baseline_scaled_mean = int(quadrature["scaled_floor"])
    hard_scaled_mean = 22
    lift_scaled_mass = hard_scaled_mean - baseline_scaled_mean
    support_density_numerator = lift_scaled_mass
    support_density_denominator = 4 * P
    support_size = (
        DOMAIN_SIZE * support_density_numerator // support_density_denominator
    )

    pointwise = bool(
        quadrature["exact_positive_quadrature_certificate"]
        and nodes == all_intersection_layers
        and all(weight > 0 for weight in weights)
        and baseline_values == parity_values == (1, 0, 1)
        and rules["proved"]
        and rules["b2_phase_one_equality_is_pointwise_XNOR"]
        and rules[BRANCH_B2]["baseline"] == "A=(1-x_i-x_j)^2"
    )
    proved = bool(
        pointwise
        and baseline_scaled_mean == 12
        and lift_scaled_mass == 10
        and lift["proved"]
        and lift["sharp_scaled_floor"] == 10
        and lift["H_at_least_two_scaled_floor"] == 12
        and support_size == SUPPORT330
        and DOMAIN_SIZE * support_density_numerator
        % support_density_denominator
        == 0
    )
    _require(proved, "the mass-ten pointwise Boolean bridge changed")
    return {
        "p": P,
        "slice": "J(13,7)",
        "hard_phase": 1,
        "hard_b": 2,
        "hard_scaled_mean": hard_scaled_mean,
        "quadrature_nodes": list(nodes),
        "quadrature_weights": [str(weight) for weight in weights],
        "every_intersection_layer_has_positive_weight": True,
        "baseline_values_on_layers": list(baseline_values),
        "parity_values_on_layers": list(parity_values),
        "pointwise_baseline": "A_0=(1-x_i-x_j)^2",
        "difference_lift": "B=(A-A_0)/2",
        "difference_is_integral_by_parity": True,
        "difference_is_globally_nonnegative": True,
        "baseline_scaled_mean": baseline_scaled_mean,
        "lift_scaled_mass": "4p*E[B]=10",
        "prop_15688_sharp_scaled_floor": int(lift["sharp_scaled_floor"]),
        "prop_15688_H_ge_2_scaled_floor": int(
            lift["H_at_least_two_scaled_floor"]
        ),
        "height_forced": 1,
        "lift_is_boolean": True,
        "support_density": "10/52=5/26",
        "support_size": support_size,
        "pointwise_bridge_precedes_nonnegative_lift_theorem": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def exact_mass10_boolean_classification() -> dict[str, object]:
    """Validate the exact support-330 model and its terminal evidence."""
    bridge = mass10_boolean_lift_bridge()
    catalog = catalog_arithmetic_certificate()
    model, _values, metadata = build_classifier_model()
    raw = CLASSIFIER_EVIDENCE.read_bytes()
    evidence = json.loads(raw)
    solver = evidence.get("solver", {})
    classification = evidence.get("classification", {})
    evidence_model = evidence.get("model", {})
    evidence_catalog = evidence.get("candidate_catalog", {})

    proved = bool(
        bridge["proved"]
        and model.Validate() == ""
        and metadata["boolean_variable_count"] == DOMAIN_SIZE
        and metadata["constraint_count"] == 1710
        and metadata["third_difference_equality_count"] == 1638
        and metadata["anchored_nogood_count"] == 70
        and metadata["shard_fixing_count"] == 0
        and metadata["model_textproto_sha256"] == EXPECTED_MODEL_SHA256
        and metadata["third_difference_identity_sha256"]
        == EXPECTED_IDENTITY_SHA256
        and metadata["candidate_catalog_sha256"] == EXPECTED_CATALOG_SHA256
        and metadata["anchored_candidate_catalog_sha256"]
        == EXPECTED_ANCHORED_CATALOG_SHA256
        and catalog["candidate_catalog_verified"]
        and catalog["full_candidate_count"] == 364
        and catalog["distinct_candidate_count"] == 364
        and catalog["family_counts"]
        == {"omitted_pair": 78, "all_equal_triple": 286}
        and catalog["anchored_family_counts"]
        == {"omitted_pair": 15, "all_equal_triple": 55}
        and catalog["support_point_anchor_is_wlog_by_S13_transitivity"]
        and catalog["candidate_catalog_is_S13_invariant"]
        and catalog["S13_symmetry_certificate"]["catalog_images_checked"]
        == 4368
        and catalog["S13_symmetry_certificate"]["anchor_orbit_size"]
        == DOMAIN_SIZE
        and evidence["result_status"]
        == "COMPLETE_EXACT_FULL_MODEL_INFEASIBILITY"
        and solver.get("status") == "INFEASIBLE"
        and int(solver.get("status_code", -1)) == 3
        and classification.get("full_catalog_exhaustive") is True
        and classification.get("counterexample_found") is False
        and classification.get("incomplete") is False
        and evidence.get("witness") is None
        and evidence_model.get("model_textproto_sha256")
        == EXPECTED_MODEL_SHA256
        and evidence_model.get("constraint_count") == 1710
        and evidence_model.get("partition", {}).get("shard_count") == 1
        and evidence_catalog.get("candidate_catalog_sha256")
        == EXPECTED_CATALOG_SHA256
        and evidence_catalog.get("anchored_candidate_catalog_sha256")
        == EXPECTED_ANCHORED_CATALOG_SHA256
    )
    _require(proved, "the exact support-330 classification evidence changed")
    return {
        "scope": "Boolean degree-at-most-two functions of support 330 on J(13,7)",
        "candidate_families": {
            "omitted_pair": {
                "polynomial": "(1-x_i)(1-x_j)",
                "count": 78,
                "anchored_count": 15,
            },
            "all_equal_triple": {
                "polynomial": (
                    "1-x_i-x_j-x_k+x_i*x_j+x_i*x_k+x_j*x_k"
                ),
                "count": 286,
                "anchored_count": 55,
            },
        },
        "full_candidate_count": 364,
        "support_size": SUPPORT330,
        "third_difference_identity_count": 1638,
        "third_difference_identity_sha256": EXPECTED_IDENTITY_SHA256,
        "identity_nullspace_is_exact_degree_at_most_two_space": True,
        "boolean_variable_count": DOMAIN_SIZE,
        "constraint_count": 1710,
        "anchored_nogood_count": 70,
        "model_textproto_sha256": EXPECTED_MODEL_SHA256,
        "candidate_catalog_sha256": EXPECTED_CATALOG_SHA256,
        "anchored_candidate_catalog_sha256": EXPECTED_ANCHORED_CATALOG_SHA256,
        "S13_generator_images_checked": 4368,
        "anchor_orbit_size": DOMAIN_SIZE,
        "solver": {
            "name": solver["name"],
            "version": solver["version"],
            "status": solver["status"],
            "status_code": int(solver["status_code"]),
            "num_search_workers": int(solver["num_search_workers"]),
            "exact_terminal_status": True,
        },
        "classifier_evidence": str(CLASSIFIER_EVIDENCE.relative_to(ROOT)),
        "gpu_cross_check_is_not_a_proof_premise": (
            evidence.get("gpu_catalog_cross_check", {}).get("proof_premise")
            is False
        ),
        "catalog_exhaustive_at_support_330": True,
        "result_status": "exhaustive finite certificate",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def mass12_phase_zero_dichotomy() -> dict[str, object]:
    """Classify the local forms still possible at opposite mean twelve."""
    floors = residual_even_floor_table(P)
    phase_zero = {
        int(b): int(value) for b, value in floors["phase_zero_floors"].items()
    }
    lift = sharp_integral_quadratic_lift_floor(P)
    complement_literal = parity_floor_certificate(P, 1, 1)
    scaled_mass = 12
    stabilizer_coefficient = lift[
        "H_at_least_two_stabilizer_coefficient"
    ]
    possible_H_ge_two = [
        height
        for height in range(2, 1 + scaled_mass)
        if max(
            2 * (P + 1) - 4 * height,
            stabilizer_coefficient * height,
        )
        <= scaled_mass
    ]
    boolean_support = DOMAIN_SIZE * scaled_mass // (4 * P)
    compatible_b = [b for b, floor in phase_zero.items() if floor <= scaled_mass]
    proved = bool(
        floors["proved"]
        and phase_zero == {0: 0, 2: 14, 4: 20, 6: 26, 8: 24, 10: 26, 12: 12}
        and compatible_b == [0, 12]
        and complement_literal["exact_positive_quadrature_certificate"]
        and complement_literal["scaled_floor"] == scaled_mass
        and all(weight > 0 for weight in complement_literal["quadrature_weights"])
        and lift["proved"]
        and possible_H_ge_two == [4]
        and boolean_support == SUPPORT396
    )
    _require(proved, "the phase-zero mass-twelve dichotomy changed")
    return {
        "p": P,
        "phase": 0,
        "scaled_mean": scaled_mass,
        "phase_zero_even_b_floors": {
            str(b): value for b, value in phase_zero.items()
        },
        "floor_compatible_b": compatible_b,
        "literal_branch": {
            "b": 12,
            "pointwise_form": "A=1-x_j",
            "target": "3+2A=4-z_j",
            "coefficient_offset": 3,
            "positive_quadrature_rigidity": True,
        },
        "lift_branch": {
            "b": 0,
            "form": "A=2C",
            "scaled_mass": "4p*E[C]=12",
            "C_is_nonzero_nonnegative_integral_quadratic": True,
            "height_dichotomy": [1, 4],
            "height_one_is_boolean": True,
            "height_one_support_size": boolean_support,
            "H_ge_two_lower_bounds_meet_only_at_H": possible_H_ge_two,
        },
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p3_degree_six_coupling() -> dict[str, object]:
    """Certify the sextic identity and its opposite-cell sign conversion."""
    pairs = tuple(combinations(range(P), 2))
    degrees = (2, 4, 6)
    gauge_rows: dict[str, dict[str, object]] = {}
    for degree in degrees:
        complete_sum = sum((left - right) ** degree for left, right in pairs) % P
        star_sums = tuple(
            sum(
                (center - other) ** degree
                for other in range(P)
                if other != center
            )
            % P
            for center in range(P)
        )
        gauge_rows[str(degree)] = {
            "complete_graph_sum_mod_13": complete_sum,
            "star_sums_mod_13": sorted(set(star_sums)),
            "all_vanish": complete_sum == 0 and set(star_sums) == {0},
        }

    overlap_count = 0
    identity_values = set()
    for baseline_left, baseline_right in pairs:
        for lift_left, lift_right in pairs:
            moments = {
                degree: (
                    (baseline_left - baseline_right) ** degree
                    + (lift_left - lift_right) ** degree
                )
                % P
                for degree in degrees
            }
            identity_values.add(
                (
                    2 * moments[6]
                    + moments[2] ** 3
                    - 3 * moments[2] * moments[4]
                )
                % P
            )
            overlap_count += 1

    # If N_d=(-h)M_d is the opposite cell's own normalization, substituting
    # M_d=-hN_d into the global F6 and multiplying by -1 changes the mixed
    # term from -3 to +3.  Pin that conversion explicitly for the next model.
    sign_conversion_checks = []
    for sign_h in (-1, 1):
        for n2 in range(P):
            for n4 in range(P):
                for n6 in range(P):
                    m2 = (-sign_h * n2) % P
                    m4 = (-sign_h * n4) % P
                    m6 = (-sign_h * n6) % P
                    global_value = (
                        2 * sign_h * m6
                        + sign_h * m2**3
                        - 3 * m2 * m4
                    ) % P
                    opposite_local_value = (
                        2 * n6 + n2**3 + 3 * n2 * n4
                    ) % P
                    sign_conversion_checks.append(
                        global_value == (-opposite_local_value) % P
                    )

    root_count = HARD_DIRECTION_COUNT
    sextic_degree = 6
    proved = bool(
        len(pairs) == 78
        and overlap_count == 78**2 == 6084
        and identity_values == {0}
        and all(row["all_vanish"] for row in gauge_rows.values())
        and all(sign_conversion_checks)
        and root_count > sextic_degree
    )
    _require(proved, "the P=3 degree-six coupling changed")
    return {
        "field": "F_13",
        "hard_family": "baseline pair plus omitted-pair lift",
        "normalized_moment_form": "N_(2r)=h*M_(2r)=alpha^(2r)+beta^(2r)",
        "hard_normalized_Newton_identity": "2*N6+N2^3-3*N2*N4=0",
        "global_sextic": "F6=2*h*M6+h*M2^3-3*M2*M4",
        "pair_count": len(pairs),
        "overlap_inclusive_pair_choices_checked": overlap_count,
        "identity_residue_set": sorted(identity_values),
        "gauge_cancellation": gauge_rows,
        "homogeneous_degree": sextic_degree,
        "distinct_hard_projective_roots": root_count,
        "root_count_exceeds_degree": root_count > sextic_degree,
        "global_F6_is_identically_zero": True,
        "opposite_sign_normalization": "N_d=(-h)*M_d",
        "opposite_local_constraint": "2*N6+N2^3+3*N2*N4=0",
        "opposite_sign_conversion_assignments_checked": len(
            sign_conversion_checks
        ),
        "opposite_sign_conversion_verified": all(sign_conversion_checks),
        "forced_P3_Q5_mass12_cell_satisfies_global_F6_zero": True,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p5_degree_six_no_identity_audit() -> dict[str, object]:
    """Audit weighted-homogeneous even-moment identities through degree six."""
    pairs = tuple(combinations(range(P), 2))
    triples = tuple(combinations(range(P), 3))
    degree_two_rows: list[tuple[int, ...]] = []
    degree_four_rows: list[tuple[int, ...]] = []
    degree_six_rows: list[tuple[int, ...]] = []
    for baseline_left, baseline_right in pairs:
        for i, j, k in triples:
            moments = {}
            for degree in (2, 4, 6):
                moments[degree] = (
                    (baseline_left - baseline_right) ** degree
                    + (i - j) ** degree
                    + (i - k) ** degree
                    + (j - k) ** degree
                ) % P
            degree_two_rows.append((moments[2],))
            degree_four_rows.append((moments[4], moments[2] ** 2 % P))
            degree_six_rows.append(
                (
                    moments[6],
                    moments[2] * moments[4] % P,
                    moments[2] ** 3 % P,
                )
            )
    ranks = {
        2: _mod_rank(degree_two_rows, P),
        4: _mod_rank(degree_four_rows, P),
        6: _mod_rank(degree_six_rows, P),
    }
    expected_ranks = {2: 1, 4: 2, 6: 3}
    proved = bool(
        len(pairs) == 78
        and len(triples) == 286
        and len(degree_six_rows) == 22308
        and ranks == expected_ranks
    )
    _require(proved, "the P=5 degree-six no-identity audit changed")
    return {
        "field": "F_13",
        "hard_family": "baseline pair plus all-equal-triple lift",
        "pattern_count_checked": len(degree_six_rows),
        "weighted_feature_vectors": {
            "2": ["N2"],
            "4": ["N4", "N2^2"],
            "6": ["N6", "N2*N4", "N2^3"],
        },
        "weighted_feature_ranks_mod_13": {
            str(degree): rank for degree, rank in ranks.items()
        },
        "weighted_feature_dimensions": {"2": 1, "4": 2, "6": 3},
        "feature_rank_mod_13": ranks[6],
        "feature_dimension": 3,
        "no_nonzero_universal_weighted_homogeneous_even_moment_identity_through_degree_6": (
            ranks == expected_ranks
        ),
        "scope": (
            "weighted-homogeneous polynomial identities in the even moments "
            "N2,N4,N6 through degree six"
        ),
        "proved": proved,
    }


@lru_cache(maxsize=1)
def t4_u4_catalog_consequence() -> dict[str, object]:
    """Propagate the exact catalog through the common p13 t4 edge ledger."""
    sieve = t4_all_residue_sieve()
    bridge = mass10_boolean_lift_bridge()
    classification = exact_mass10_boolean_classification()
    mass12 = mass12_phase_zero_dichotomy()
    sextic = p3_degree_six_coupling()
    p5_audit = p5_degree_six_no_identity_audit()

    u4_row = next(
        row for row in sieve["residue_rows"] if int(row["u"]) == 4
    )
    family_offsets = {"omitted_pair": 3, "all_equal_triple": 5}
    family_rows: dict[str, dict[str, object]] = {}
    for family, offset in family_offsets.items():
        parallel_candidates = [
            parallel
            for parallel in range(H_EDGE_COUNT // HARD_DIRECTION_COUNT + 1)
            if (parallel - offset) % Q_MODULUS == 0
        ]
        _require(
            len(parallel_candidates) == 1,
            f"the {family} hard parallel count is no longer unique",
        )
        parallel = parallel_candidates[0]
        h_times_T = 14 * parallel - 61
        opposite_sum = H_EDGE_COUNT - HARD_DIRECTION_COUNT * parallel
        minimum_Q = 8 - parallel
        excess_sum = opposite_sum - HARD_DIRECTION_COUNT * minimum_Q
        minimum_mean = 14 * (parallel + minimum_Q) - 100
        minimum_cell_count = HARD_DIRECTION_COUNT - excess_sum
        literal_compatible = (minimum_Q - 3) % Q_MODULUS == 0
        row = {
            "family": family,
            "coefficient_offset": offset,
            "coefficient_congruence": f"6 divides P-{offset}",
            "parallel_candidates_under_7P_le_61": parallel_candidates,
            "common_hard_parallel_count_P": parallel,
            "hard_signed_total_hT": h_times_T,
            "hard_mean_identity": f"22=14*{parallel}-({h_times_T})-39",
            "opposite_parallel_sum": opposite_sum,
            "minimum_opposite_Q": minimum_Q,
            "opposite_excess_definition": "e=P+Q-8",
            "opposite_excess_sum": excess_sum,
            "opposite_mean_formula": "a=12+14e",
            "minimum_opposite_scaled_mean": minimum_mean,
            "directions_at_minimum_at_least": minimum_cell_count,
            "b12_literal_offset": 3,
            "b12_literal_compatible_at_minimum_Q": literal_compatible,
        }
        if family == "omitted_pair":
            row.update(
                {
                    "forced_minimum_cell": "b=0,A=2C,4p*E[C]=12",
                    "minimum_cell_height_dichotomy": [1, 4],
                    "height_one_support_size": SUPPORT396,
                    "global_sextic_constraint": sextic["global_sextic"],
                    "opposite_local_sextic_constraint": sextic[
                        "opposite_local_constraint"
                    ],
                }
            )
        else:
            row.update(
                {
                    "minimum_cell_dichotomy": (
                        "b=12 literal or b=0 mass-12 lift"
                    ),
                    "degree_six_feature_rank": p5_audit[
                        "feature_rank_mod_13"
                    ],
                    "analogous_degree_six_identity_available": False,
                }
            )
        family_rows[family] = row

    hard_parallel_values = {
        row["common_hard_parallel_count_P"] for row in family_rows.values()
    }
    proved = bool(
        sieve["proved"]
        and sieve["surviving_residues_before_prop_15744"] == [0, 3, 4, 6]
        and u4_row["forced_low_direction_count_at_least"] == 7
        and u4_row["forced_low_mean"] == 22
        and u4_row["surviving_low_cells"]
        == [
            {
                "b": 2,
                "floor": 12,
                "excess": 10,
                "status": "sharp lift equality survives",
                "survives": True,
            }
        ]
        and bridge["proved"]
        and classification["catalog_exhaustive_at_support_330"]
        and mass12["proved"]
        and sextic["proved"]
        and p5_audit["proved"]
        and hard_parallel_values == {3, 5}
        and family_rows["omitted_pair"]["opposite_parallel_sum"] == 40
        and family_rows["all_equal_triple"]["opposite_parallel_sum"] == 26
        and all(
            row["opposite_excess_sum"] == 5
            and row["minimum_opposite_scaled_mean"] == 12
            and row["directions_at_minimum_at_least"] == 2
            for row in family_rows.values()
        )
        and not family_rows["omitted_pair"][
            "b12_literal_compatible_at_minimum_Q"
        ]
        and family_rows["all_equal_triple"][
            "b12_literal_compatible_at_minimum_Q"
        ]
    )
    _require(proved, "the p13 t4 u4 catalog consequence changed")
    return {
        "p": P,
        "t": 4,
        "k": 60,
        "u": 4,
        "hard_direction_count": HARD_DIRECTION_COUNT,
        "H_edge_count": H_EDGE_COUNT,
        "hard_cell": {
            "phase": 1,
            "b": 2,
            "scaled_mean": 22,
            "pointwise_form": "A=(1-x_a-x_b)^2+2B",
            "lift_scaled_mass": 10,
            "lift_support_size": SUPPORT330,
        },
        "family_ledgers": family_rows,
        "families_cannot_mix_because_common_P_has_distinct_mod_6_offsets": True,
        "mass12_phase_zero_dichotomy": mass12,
        "P3_degree_six_coupling": sextic,
        "P5_degree_six_no_identity_audit": p5_audit,
        "p13_t4_u4_closed": False,
        "remaining_p13_t4_residues": [4, 6],
        "result_status": "proved open reduction",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15746() -> dict[str, object]:
    """Package the support-330 theorem and exact open reduction for ``u=4``."""
    bridge = mass10_boolean_lift_bridge()
    classification = exact_mass10_boolean_classification()
    consequence = t4_u4_catalog_consequence()
    proved = bool(
        bridge["proved"]
        and classification["proved"]
        and consequence["proved"]
        and not consequence["p13_t4_u4_closed"]
        and consequence["remaining_p13_t4_residues"] == [4, 6]
    )
    _require(proved, "Proposition 15.746 failed")
    return {
        "prop": "15.746",
        "title": "Sharp mass-ten catalog and the p13 t4 u4 reduction",
        "result_status": (
            "exhaustive finite equality classification and proved open reduction"
        ),
        "statement": (
            "the support-330 Boolean catalog is exhaustive; u=4 reduces to "
            "uniform P=3 or P=5 families and remains open"
        ),
        "pointwise_mass10_bridge": bridge,
        "exact_support330_classification": classification,
        "u4_catalog_consequence": consequence,
        "p13_t4_u4_closed": False,
        "p13_k_eq_60_closed": False,
        "remaining_p13_t4_residues": [4, 6],
        "next_exact_gate": (
            "one P=3,Q=5,b=0 mass-12 cell under the opposite-normalized "
            "constraint 2*N6+N2^3+3*N2*N4=0, split into Boolean support "
            "396 or height four"
        ),
        "broad_mass12_or_support396_census_is_not_the_gate": True,
        "residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic Proposition 15.746 package atomically."""
    target = ROOT / "evidence" / "e1_gmin_m4_prop15746.json"
    write_json_atomic(target, proposition_15746())
    return target


def main() -> None:
    result = proposition_15746()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.746 failed")
    target = write_evidence()
    print("Prop 15.746 support-330 catalog: exhaustive")
    print("p=13,t=4,u=4: open reduction; remaining residues [4, 6]")
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
