#!/usr/bin/env python3
"""Exact signed-degree presolver for the p=31 branch-C top endpoint.

This model does not start from the paired-SDR or localized-Mobius gauge.
It starts from an arbitrary simple graph on ``F_31^2`` satisfying the
one-half central-inversion support condition.  It retains the exact balanced
parallel quotas and the signed vertex boundaries forced by the sharp local
branch-C atoms, but eliminates the 461,280 physical edge variables through
the affine-line incidence inverse

    p g(v) = sum_D r[D, D(v)] - G,

where ``g(v)`` is the raw Paley-signed graph degree and

    r[D,s] = eps[D] * (d[D,s] + 2 n[D,s]).

Here ``d[D,s]`` is the normalized transverse signed degree and ``n[D,s]``
is the number of selected direction-D edges inside the affine line D=s.
The identity follows from ``B.T B = p I + J`` for the complete affine-plane
line incidence matrix.  The model additionally enforces exact positive- and
negative-Paley incidence budgets and the capacities imposed on each
antipodal vertex pair by choosing at most one edge from every nonfixed
central-inversion orbit and exactly one fixed antipodal edge.

The target signed-degree projections are exact:

* a hard row is one negative unit star plus ``e_D`` compact atoms;
* an opposite row is six all-positive triangles plus ``Q_D-9`` compact
  atoms.

For six repeated distinct-label triangles, occurrence counts ``A_s`` are
realizable iff ``0 <= A_s <= 6`` and ``sum A_s = 18`` (the corresponding
31-by-6 bipartite degree sequence satisfies Gale--Ryser automatically).
Compact distinguished-label counts are unrestricted weak compositions.

UNSAT therefore excludes the complete p=31, t=177 balanced compact-atom
one-half family at this signed-degree layer.  SAT proves only that this exact
necessary projection survives: it is not a common graph and does not close
residual (ii).  Use ``--export-model`` to write a deterministic CP-SAT proto
for replay on a mesh CPU node.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_inversion_antisymmetric_radon import (  # noqa: E402
    Point,
    _functional_value,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (  # noqa: E402
    paley_direction_sign,
    paley_edge_sign,
)
from io_atomic import write_json_atomic  # noqa: E402


EXPERIMENT = "residual_branch_c_p31_top_signed_degree"
P = 31
R = 7
T = 177
H_SIZE = 125 + 2 * T
E_TOTAL = T + 1
Q_TOTAL = 10 * R + 6 + T
FIXED_EDGE_COUNT = 1


@dataclass(frozen=True)
class Geometry:
    points: tuple[Point, ...]
    directions: tuple[Point, ...]
    epsilon: tuple[int, ...]
    labels: tuple[tuple[int, ...], ...]
    negative_point: tuple[int, ...]
    antipodal_representatives: tuple[int, ...]
    spatial_direction: tuple[int, ...]


@dataclass
class ModelHandles:
    e: dict[int, Any]
    q: dict[int, Any]
    parallel: dict[tuple[int, int], Any]
    centre: dict[tuple[int, int], Any]
    hard_compact: dict[tuple[int, int], Any]
    opposite_positive: dict[tuple[int, int], Any]
    opposite_compact: dict[tuple[int, int], Any]
    transverse_degree: dict[tuple[int, int], Any]
    signed_degree: dict[int, Any]
    positive_degree: dict[int, Any]
    negative_degree: dict[int, Any]
    fixed_edge: dict[int, Any]
    fixed_direction: dict[int, Any]
    auxiliary_direction: dict[int, Any]


def build_geometry() -> Geometry:
    points = tuple((x, y) for x in range(P) for y in range(P))
    point_index = {point: index for index, point in enumerate(points)}
    directions = projective_functionals(P)
    epsilon = tuple(paley_direction_sign(P, direction) for direction in directions)
    if epsilon.count(1) != 16 or epsilon.count(-1) != 16:
        raise ArithmeticError("the p=31 projective Paley split changed")
    labels = tuple(
        tuple(_functional_value(P, direction, point) for point in points)
        for direction in directions
    )
    negative_point = tuple(
        point_index[(-point[0] % P, -point[1] % P)] for point in points
    )
    representatives = tuple(
        index
        for index, point in enumerate(points)
        if point != (0, 0) and index < negative_point[index]
    )
    if len(representatives) != (P * P - 1) // 2:
        raise ArithmeticError("the antipodal point-pair count changed")

    # The spatial direction of {v,-v} is the unique functional annihilating v.
    spatial = [-1] * len(points)
    for index, point in enumerate(points):
        if point == (0, 0):
            continue
        annihilators = [
            direction_index
            for direction_index in range(P + 1)
            if labels[direction_index][index] == 0
        ]
        if len(annihilators) != 1:
            raise ArithmeticError("a nonzero point lost its projective annihilator")
        spatial[index] = annihilators[0]
    return Geometry(
        points=points,
        directions=directions,
        epsilon=epsilon,
        labels=labels,
        negative_point=tuple(negative_point),
        antipodal_representatives=representatives,
        spatial_direction=tuple(spatial),
    )


def geometry_identity_self_test(geometry: Geometry) -> dict[str, int | bool]:
    """Replay the line inverse and sign conventions on deterministic edges."""
    sample_edges: list[tuple[int, int]] = []
    # A nontrivial deterministic mix, including antipodal and nonfixed edges.
    for first in range(0, len(geometry.points), 73):
        for offset in (1, 37, 211):
            second = (first + offset) % len(geometry.points)
            if first != second:
                sample_edges.append(tuple(sorted((first, second))))
    sample_edges = sorted(set(sample_edges))[:37]
    raw_degree = [0] * len(geometry.points)
    parallel = [[0] * P for _ in range(P + 1)]
    transverse = [[0] * P for _ in range(P + 1)]
    sign_sum = 0
    for first, second in sample_edges:
        edge = (geometry.points[first], geometry.points[second])
        tau = paley_edge_sign(P, edge)
        raw_degree[first] += tau
        raw_degree[second] += tau
        sign_sum += 2 * tau
        parallel_rows = []
        for direction_index, eps in enumerate(geometry.epsilon):
            left = geometry.labels[direction_index][first]
            right = geometry.labels[direction_index][second]
            if left == right:
                parallel[direction_index][left] += 1
                parallel_rows.append(direction_index)
            else:
                transverse[direction_index][left] += eps * tau
                transverse[direction_index][right] += eps * tau
        if len(parallel_rows) != 1:
            raise ArithmeticError("an edge lost its unique parallel direction")
        if geometry.epsilon[parallel_rows[0]] != tau:
            raise ArithmeticError("parallel direction and edge Paley signs disagree")

    for point_index in range(len(geometry.points)):
        reconstructed = -sign_sum
        for direction_index, eps in enumerate(geometry.epsilon):
            label = geometry.labels[direction_index][point_index]
            line_sum = eps * (
                transverse[direction_index][label]
                + 2 * parallel[direction_index][label]
            )
            reconstructed += line_sum
        if reconstructed != P * raw_degree[point_index]:
            raise ArithmeticError("B.T B line inverse failed in exact replay")
    return {
        "proved": True,
        "sample_edge_count": len(sample_edges),
        "raw_signed_degree_sum": sum(raw_degree),
        "twice_edge_sign_sum": sign_sum,
    }


def target_row_sum_certificate() -> dict[str, object]:
    """Symbolically check that every directional line-sum total is G."""
    hard_rows = {}
    for e_value in (11, 12):
        transverse_sum = -2 * (P - 1) - 2 * e_value
        parallel_quota = 3 + e_value
        raw_line_sum = transverse_sum + 2 * parallel_quota
        hard_rows[str(e_value)] = {
            "transverse_degree_sum": transverse_sum,
            "parallel_quota": parallel_quota,
            "raw_line_sum": raw_line_sum,
        }
    opposite_rows = {}
    for q_value in (15, 16):
        transverse_sum = 6 * 6 - 2 * (q_value - 9)
        parallel_quota = q_value
        raw_line_sum = -(transverse_sum + 2 * parallel_quota)
        opposite_rows[str(q_value)] = {
            "transverse_degree_sum": transverse_sum,
            "parallel_quota": parallel_quota,
            "raw_line_sum": raw_line_sum,
        }
    expected = 8 - 2 * P
    proved = all(row["raw_line_sum"] == expected for row in hard_rows.values()) and all(
        row["raw_line_sum"] == expected for row in opposite_rows.values()
    )
    if not proved:
        raise ArithmeticError("hard/opposite directional line totals diverged")
    return {
        "hard_rows": hard_rows,
        "opposite_rows": opposite_rows,
        "common_raw_line_sum": expected,
        "proved": True,
    }


def mobius_top_origin_profile(fixed_edge_sign: str) -> dict[str, int]:
    """Return the parity-forced origin degrees in one fixed-sign shard."""
    if fixed_edge_sign == "hard":
        positive, negative = 13, 3
    elif fixed_edge_sign == "opposite":
        positive, negative = 14, 2
    else:
        raise ValueError("the exact profile needs a hard or opposite fixed edge")
    return {
        "positive_degree": positive,
        "negative_degree": negative,
        "signed_degree": positive - negative,
        "ordinary_degree": positive + negative,
    }


def mobius_top_auxiliary_indicator(
    geometry: Geometry,
    e_values: dict[int, int],
    q_values: dict[int, int],
    fixed_direction: int,
) -> tuple[int, ...]:
    """Replay ``z_D=v_D-fdir_D`` for a concrete balanced top profile."""
    if set(e_values) != {
        direction for direction, eps in enumerate(geometry.epsilon) if eps == 1
    }:
        raise ValueError("hard profile has the wrong direction keys")
    if set(q_values) != {
        direction for direction, eps in enumerate(geometry.epsilon) if eps == -1
    }:
        raise ValueError("opposite profile has the wrong direction keys")
    if fixed_direction not in range(P + 1):
        raise ValueError("fixed direction is out of range")
    v = tuple(
        12 - e_values[direction]
        if direction in e_values
        else 16 - q_values[direction]
        for direction in range(P + 1)
    )
    if set(v) - {0, 1} or sum(v) != 17 or v[fixed_direction] != 1:
        raise ValueError("profile/fixed direction violates the exact top parity word")
    z = tuple(
        value - int(direction == fixed_direction)
        for direction, value in enumerate(v)
    )
    if set(z) - {0, 1} or sum(z) != 16:
        raise ArithmeticError("the auxiliary indicator is not a 16-set")
    return z


def build_model(
    geometry: Geometry,
    fixed_edge_sign: str = "any",
    mobius_top_origin: bool = False,
) -> tuple[Any, ModelHandles, dict[str, object]]:
    from ortools.sat.python import cp_model

    if fixed_edge_sign not in {"any", "hard", "opposite"}:
        raise ValueError("fixed_edge_sign must be any, hard, or opposite")
    model = cp_model.CpModel()
    hard = tuple(index for index, eps in enumerate(geometry.epsilon) if eps == 1)
    opposite = tuple(index for index, eps in enumerate(geometry.epsilon) if eps == -1)

    e = {direction: model.new_int_var(11, 12, f"e_{direction}") for direction in hard}
    q = {direction: model.new_int_var(15, 16, f"Q_{direction}") for direction in opposite}
    model.add(sum(e.values()) == E_TOTAL)
    model.add(sum(q.values()) == Q_TOTAL)

    parallel: dict[tuple[int, int], Any] = {}
    centre: dict[tuple[int, int], Any] = {}
    hard_compact: dict[tuple[int, int], Any] = {}
    opposite_positive: dict[tuple[int, int], Any] = {}
    opposite_compact: dict[tuple[int, int], Any] = {}
    transverse_degree: dict[tuple[int, int], Any] = {}

    for direction in range(P + 1):
        quota = 3 + e[direction] if direction in e else q[direction]
        row_parallel = []
        for label in range(P):
            value = model.new_int_var(0, 16, f"n_{direction}_{label}")
            parallel[direction, label] = value
            row_parallel.append(value)
        model.add(sum(row_parallel) == quota)

        if direction in e:
            row_centres = []
            row_compact = []
            for label in range(P):
                # The localized-half slice needs nonzero centers; the full
                # compact target family also permits the zero center.
                is_centre = model.new_bool_var(f"J_{direction}_{label}")
                centre[direction, label] = is_centre
                if mobius_top_origin and label == 0:
                    model.add(is_centre == 0)
                else:
                    row_centres.append(is_centre)
                compact = model.new_int_var(0, 12, f"C_h_{direction}_{label}")
                hard_compact[direction, label] = compact
                row_compact.append(compact)
                degree = model.new_int_var(-54, -1, f"d_{direction}_{label}")
                transverse_degree[direction, label] = degree
                model.add(degree == -1 - (P - 2) * is_centre - 2 * compact)
            model.add(sum(row_centres) == 1)
            model.add(sum(row_compact) == e[direction])
        else:
            row_positive = []
            row_compact = []
            for label in range(P):
                positive = model.new_int_var(0, 6, f"A_o_{direction}_{label}")
                compact = model.new_int_var(0, 7, f"C_o_{direction}_{label}")
                opposite_positive[direction, label] = positive
                opposite_compact[direction, label] = compact
                row_positive.append(positive)
                row_compact.append(compact)
                degree = model.new_int_var(-14, 12, f"d_{direction}_{label}")
                transverse_degree[direction, label] = degree
                model.add(degree == 2 * (positive - compact))
            model.add(sum(row_positive) == 18)
            model.add(sum(row_compact) == q[direction] - 9)

    # Exactly one fixed antipodal edge.  Its spatial direction determines the
    # one central affine line in which it contributes to n[D,0].
    fixed_edge = {
        representative: model.new_bool_var(f"f_{representative}")
        for representative in geometry.antipodal_representatives
    }
    model.add(sum(fixed_edge.values()) == FIXED_EDGE_COUNT)
    if fixed_edge_sign != "any":
        requested_epsilon = 1 if fixed_edge_sign == "hard" else -1
        model.add(
            sum(
                fixed_edge[representative]
                for representative in geometry.antipodal_representatives
                if geometry.epsilon[geometry.spatial_direction[representative]]
                == requested_epsilon
            )
            == FIXED_EDGE_COUNT
        )
    fixed_direction = {
        direction: model.new_bool_var(f"fdir_{direction}")
        for direction in range(P + 1)
    }
    for direction in range(P + 1):
        members = [
            fixed_edge[representative]
            for representative in geometry.antipodal_representatives
            if geometry.spatial_direction[representative] == direction
        ]
        if len(members) != (P - 1) // 2:
            raise ArithmeticError("a spatial direction lost antipodal edges")
        model.add(fixed_direction[direction] == sum(members))
        model.add(parallel[direction, 0] >= fixed_direction[direction])
        model.add(parallel[direction, 0] <= 225 + fixed_direction[direction])
        for label in range(1, (P + 1) // 2):
            # Lines s and -s are exchanged by central inversion and contain
            # 465 paired physical edges in total.
            model.add(parallel[direction, label] + parallel[direction, -label % P] <= 465)

    # Recover raw signed point degrees by exact affine-line inversion.  The
    # global signed degree is fixed by the direction quotas.
    hard_edge_count = 3 * len(hard) + E_TOTAL
    opposite_edge_count = Q_TOTAL
    global_signed_degree = 2 * (hard_edge_count - opposite_edge_count)
    if hard_edge_count + opposite_edge_count != H_SIZE:
        raise ArithmeticError("top endpoint direction quotas lost |H|")
    if global_signed_degree != 8 - 2 * P:
        raise ArithmeticError("top endpoint global signed degree changed")

    signed_degree: dict[int, Any] = {}
    positive_degree: dict[int, Any] = {}
    negative_degree: dict[int, Any] = {}
    for point_index in range(len(geometry.points)):
        signed = model.new_int_var(-480, 480, f"g_{point_index}")
        positive = model.new_int_var(0, 480, f"gplus_{point_index}")
        negative = model.new_int_var(0, 480, f"gminus_{point_index}")
        signed_degree[point_index] = signed
        positive_degree[point_index] = positive
        negative_degree[point_index] = negative
        model.add(signed == positive - negative)
        line_terms = []
        for direction, eps in enumerate(geometry.epsilon):
            label = geometry.labels[direction][point_index]
            line_terms.append(
                eps
                * (
                    transverse_degree[direction, label]
                    + 2 * parallel[direction, label]
                )
            )
        model.add(P * signed == sum(line_terms) - global_signed_degree)

    # These are the ordinary incidences of positive- and negative-Paley
    # direction edges.  They give the exact degree/parity/l1 envelope, not a
    # relaxation with an arbitrary absolute-value slack.
    model.add(sum(positive_degree.values()) == 2 * hard_edge_count)
    model.add(sum(negative_degree.values()) == 2 * opposite_edge_count)

    # At the origin every incident edge orbit is {0,w}<->{0,-w}; there are
    # 15 such orbits in each projective direction and hence 240 per Paley
    # sign.  This is the origin analogue of the antipodal-pair capacities.
    origin = geometry.points.index((0, 0))
    model.add(positive_degree[origin] <= 16 * (P - 1) // 2)
    model.add(negative_degree[origin] <= 16 * (P - 1) // 2)
    auxiliary_direction: dict[int, Any] = {}
    if mobius_top_origin:
        # At the exact top parity slice the 16 auxiliary directions are all
        # distinct, so the unique cancellation cannot be at an origin orbit:
        # a shared {0,u} would force two equal projective kernels ker(M_i).
        # Directionwise, v_D=P_D mod 2 + 1_hard and the auxiliary indicator
        # is z_D=v_D-fdir_D.  Every surviving t=0 edge is parallel to its
        # distinct auxiliary M_D and lies inside the central line D=0.
        for direction in range(P + 1):
            auxiliary = model.new_bool_var(f"zaux_{direction}")
            auxiliary_direction[direction] = auxiliary
            if direction in e:
                model.add(auxiliary + fixed_direction[direction] == 12 - e[direction])
            else:
                model.add(auxiliary + fixed_direction[direction] == 16 - q[direction])
            model.add(parallel[direction, 0] >= auxiliary)
        model.add(sum(auxiliary_direction.values()) == 16)
        model.add(
            positive_degree[origin]
            == sum(auxiliary_direction[direction] for direction in hard)
        )
        model.add(
            negative_degree[origin]
            == sum(auxiliary_direction[direction] for direction in opposite)
        )

    # Central inversion pairs the nonfixed incident edges of v and -v.  In
    # the sign class containing {v,-v} there are 479 such nonfixed orbits;
    # the other sign class has 480.  Selecting the fixed edge contributes
    # two incidences rather than one.
    for representative in geometry.antipodal_representatives:
        negative_index = geometry.negative_point[representative]
        fixed = fixed_edge[representative]
        spatial = geometry.spatial_direction[representative]
        if geometry.epsilon[spatial] == 1:
            model.add(
                positive_degree[representative] + positive_degree[negative_index]
                <= 479 + 2 * fixed
            )
            model.add(
                negative_degree[representative] + negative_degree[negative_index]
                <= 480
            )
            model.add(positive_degree[representative] >= fixed)
            model.add(positive_degree[negative_index] >= fixed)
        else:
            model.add(
                positive_degree[representative] + positive_degree[negative_index]
                <= 480
            )
            model.add(
                negative_degree[representative] + negative_degree[negative_index]
                <= 479 + 2 * fixed
            )
            model.add(negative_degree[representative] >= fixed)
            model.add(negative_degree[negative_index] >= fixed)

    handles = ModelHandles(
        e=e,
        q=q,
        parallel=parallel,
        centre=centre,
        hard_compact=hard_compact,
        opposite_positive=opposite_positive,
        opposite_compact=opposite_compact,
        transverse_degree=transverse_degree,
        signed_degree=signed_degree,
        positive_degree=positive_degree,
        negative_degree=negative_degree,
        fixed_edge=fixed_edge,
        fixed_direction=fixed_direction,
        auxiliary_direction=auxiliary_direction,
    )
    metadata: dict[str, object] = {
        "p": P,
        "r": R,
        "t": T,
        "graph_vertex_count": P * P,
        "physical_edge_count_if_expanded": (P * P) * (P * P - 1) // 2,
        "nonfixed_edge_orbit_count_if_expanded": 230_400,
        "fixed_antipodal_edge_count_if_expanded": 480,
        "selected_edge_count": H_SIZE,
        "selected_fixed_edge_count": FIXED_EDGE_COUNT,
        "fixed_edge_sign_shard": fixed_edge_sign,
        "mobius_top_origin_slice": mobius_top_origin,
        "hard_direction_indices": list(hard),
        "opposite_direction_indices": list(opposite),
        "hard_e_total": E_TOTAL,
        "hard_center_domain": "F_31^*" if mobius_top_origin else "F_31",
        "hard_parallel_total": hard_edge_count,
        "opposite_parallel_total": opposite_edge_count,
        "raw_signed_degree_total": global_signed_degree,
        "target_row_sum_certificate": target_row_sum_certificate(),
        "target_scope": (
            "balanced compact-atom branch-C p31 t177 top endpoint"
            if not mobius_top_origin
            else (
                "localized-Mobius 16-half p31 t177 top construction with exact "
                "origin incidence forced by parallel parity"
            )
        ),
        "projection": "exact target signed degrees plus necessary simple one-half graph degree envelope",
        "sat_semantics": "necessary projection survives; not a common graph",
        "unsat_semantics": (
            "full balanced compact-atom p31 t177 one-half family excluded"
            if not mobius_top_origin
            else (
                "only the localized-Mobius top-construction subfamily is excluded; "
                "not the general compact family"
            )
        ),
        "shard_semantics": (
            "unsharded full family"
            if fixed_edge_sign == "any"
            else (
                f"exact {fixed_edge_sign} fixed-edge half; combine INFEASIBLE hard and opposite shards "
                "to exclude the full fixed-edge partition"
            )
        ),
        "mobius_top_origin_semantics": (
            "no localized-half origin-incidence restriction"
            if not mobius_top_origin
            else (
                "origin cancellation excluded by distinct auxiliary directions; "
                "fixed-hard gives (gplus,gminus,g)=(13,3,10), fixed-opposite gives "
                "(14,2,12)"
            )
        ),
    }
    validation = model.validate()
    if validation:
        raise ArithmeticError(f"invalid signed-degree model: {validation}")
    proto = model.proto
    metadata["cp_sat_variable_count"] = len(proto.variables)
    metadata["cp_sat_constraint_count"] = len(proto.constraints)
    metadata["cp_sat_model_validation"] = validation
    metadata["cp_sat_textproto_sha256"] = hashlib.sha256(
        str(proto).encode("utf-8")
    ).hexdigest()
    return model, handles, metadata


def _value_matrix(solver: Any, variables: dict[tuple[int, int], Any]) -> list[list[int]]:
    return [
        [int(solver.value(variables[direction, label])) for label in range(P)]
        for direction in range(P + 1)
    ]


def realise_six_triangle_occurrences(counts: Sequence[int]) -> list[list[int]]:
    """Construct six distinct-label triples with the prescribed occurrences."""
    if len(counts) != P or any(value < 0 or value > 6 for value in counts):
        raise ValueError("triangle occurrences must be 31 integers in [0,6]")
    if sum(counts) != 18:
        raise ValueError("six triangles require total occurrence 18")
    remaining = list(counts)
    triangles: list[list[int]] = []
    for slot in range(6):
        labels = sorted(range(P), key=lambda label: (-remaining[label], label))[:3]
        if len(set(labels)) != 3 or any(remaining[label] <= 0 for label in labels):
            raise ArithmeticError("Gale--Ryser triangle realization got stuck")
        triangles.append(labels)
        for label in labels:
            remaining[label] -= 1
        slots_left = 5 - slot
        if any(value < 0 or value > slots_left for value in remaining):
            raise ArithmeticError("triangle realization violated its residual bound")
    if any(remaining):
        raise ArithmeticError("triangle realization left occurrences")
    replay = [0] * P
    for triangle in triangles:
        for label in triangle:
            replay[label] += 1
    if replay != list(counts):
        raise ArithmeticError("triangle realization failed replay")
    return triangles


def replay_witness(
    geometry: Geometry,
    solver: Any,
    handles: ModelHandles,
    *,
    mobius_top_origin: bool = False,
) -> dict[str, object]:
    hard = tuple(index for index, eps in enumerate(geometry.epsilon) if eps == 1)
    opposite = tuple(index for index, eps in enumerate(geometry.epsilon) if eps == -1)
    e_values = {direction: int(solver.value(handles.e[direction])) for direction in hard}
    q_values = {direction: int(solver.value(handles.q[direction])) for direction in opposite}
    parallel = _value_matrix(solver, handles.parallel)
    transverse = _value_matrix(solver, handles.transverse_degree)
    signed = [int(solver.value(handles.signed_degree[index])) for index in range(P * P)]
    positive = [int(solver.value(handles.positive_degree[index])) for index in range(P * P)]
    negative = [int(solver.value(handles.negative_degree[index])) for index in range(P * P)]
    fixed_representatives = [
        representative
        for representative in geometry.antipodal_representatives
        if solver.boolean_value(handles.fixed_edge[representative])
    ]
    if len(fixed_representatives) != 1:
        raise ArithmeticError("replay did not find exactly one fixed edge")

    hard_edge_count = sum(3 + value for value in e_values.values())
    opposite_edge_count = sum(q_values.values())
    global_signed_degree = 2 * (hard_edge_count - opposite_edge_count)
    if sum(e_values.values()) != E_TOTAL or set(e_values.values()) - {11, 12}:
        raise ArithmeticError("hard balanced allocation failed replay")
    if sum(q_values.values()) != Q_TOTAL or set(q_values.values()) - {15, 16}:
        raise ArithmeticError("opposite balanced allocation failed replay")
    if sum(positive) != 2 * hard_edge_count or sum(negative) != 2 * opposite_edge_count:
        raise ArithmeticError("ordinary incidence budgets failed replay")
    if any(signed[index] != positive[index] - negative[index] for index in range(P * P)):
        raise ArithmeticError("signed degree split failed replay")
    for direction in range(P + 1):
        expected = 3 + e_values[direction] if direction in e_values else q_values[direction]
        if sum(parallel[direction]) != expected:
            raise ArithmeticError("parallel line allocation failed replay")
    for point_index in range(P * P):
        reconstructed = -global_signed_degree
        for direction, eps in enumerate(geometry.epsilon):
            label = geometry.labels[direction][point_index]
            reconstructed += eps * (
                transverse[direction][label] + 2 * parallel[direction][label]
            )
        if reconstructed != P * signed[point_index]:
            raise ArithmeticError("line inverse failed witness replay")

    centre_values: dict[str, int] = {}
    hard_compact_values: dict[str, list[int]] = {}
    for direction in hard:
        chosen = [
            label
            for label in range(1, P)
            if solver.boolean_value(handles.centre[direction, label])
        ]
        if len(chosen) != 1:
            raise ArithmeticError("hard centre failed replay")
        centre_values[str(direction)] = chosen[0]
        compact = [
            int(solver.value(handles.hard_compact[direction, label]))
            for label in range(P)
        ]
        if sum(compact) != e_values[direction]:
            raise ArithmeticError("hard compact occurrence total failed replay")
        expected_degree = [
            -1 - (P - 2) * int(label == chosen[0]) - 2 * compact[label]
            for label in range(P)
        ]
        if expected_degree != transverse[direction]:
            raise ArithmeticError("hard signed-degree formula failed replay")
        hard_compact_values[str(direction)] = compact

    opposite_positive_values: dict[str, list[int]] = {}
    opposite_compact_values: dict[str, list[int]] = {}
    opposite_triangles: dict[str, list[list[int]]] = {}
    for direction in opposite:
        positive = [
            int(solver.value(handles.opposite_positive[direction, label]))
            for label in range(P)
        ]
        compact = [
            int(solver.value(handles.opposite_compact[direction, label]))
            for label in range(P)
        ]
        if sum(positive) != 18 or any(value < 0 or value > 6 for value in positive):
            raise ArithmeticError("opposite positive occurrences failed replay")
        if sum(compact) != q_values[direction] - 9:
            raise ArithmeticError("opposite compact occurrence total failed replay")
        expected_degree = [2 * (positive[label] - compact[label]) for label in range(P)]
        if expected_degree != transverse[direction]:
            raise ArithmeticError("opposite signed-degree formula failed replay")
        opposite_positive_values[str(direction)] = positive
        opposite_compact_values[str(direction)] = compact
        opposite_triangles[str(direction)] = realise_six_triangle_occurrences(positive)

    fixed_index = fixed_representatives[0]
    fixed_sign = (
        "hard"
        if geometry.epsilon[geometry.spatial_direction[fixed_index]] == 1
        else "opposite"
    )
    if mobius_top_origin:
        expected_origin = mobius_top_origin_profile(fixed_sign)
        origin = geometry.points.index((0, 0))
        actual_origin = {
            "positive_degree": positive[origin],
            "negative_degree": negative[origin],
            "signed_degree": signed[origin],
            "ordinary_degree": positive[origin] + negative[origin],
        }
        if actual_origin != expected_origin:
            raise ArithmeticError("Mobius top origin profile failed replay")
        auxiliary_directions = [
            direction
            for direction in range(P + 1)
            if solver.boolean_value(handles.auxiliary_direction[direction])
        ]
        if len(auxiliary_directions) != 16:
            raise ArithmeticError("Mobius auxiliary direction count failed replay")
        fixed_spatial_direction = geometry.spatial_direction[fixed_index]
        expected_indicator = mobius_top_auxiliary_indicator(
            geometry, e_values, q_values, fixed_spatial_direction
        )
        for direction, z_value in enumerate(expected_indicator):
            if int(direction in auxiliary_directions) != z_value:
                raise ArithmeticError("Mobius directionwise parity failed replay")
            if parallel[direction][0] < z_value:
                raise ArithmeticError("Mobius origin parallel incidence failed replay")
    else:
        auxiliary_directions = []
    witness_core = {
        "hard_e": {str(key): value for key, value in e_values.items()},
        "opposite_Q": {str(key): value for key, value in q_values.items()},
        "hard_centres": centre_values,
        "hard_compact_distinguished_counts": hard_compact_values,
        "opposite_all_positive_occurrences": opposite_positive_values,
        "opposite_all_positive_triangles": opposite_triangles,
        "opposite_compact_distinguished_counts": opposite_compact_values,
        "fixed_antipodal_edge": [
            list(geometry.points[fixed_index]),
            list(geometry.points[geometry.negative_point[fixed_index]]),
        ],
        "fixed_antipodal_edge_sign": fixed_sign,
        "mobius_auxiliary_direction_indices": auxiliary_directions,
        "parallel_line_counts": parallel,
        "transverse_signed_degrees": transverse,
        "raw_signed_point_degrees": signed,
        "positive_point_degrees": positive,
        "negative_point_degrees": negative,
    }
    encoding = json.dumps(witness_core, sort_keys=True, separators=(",", ":")).encode()
    witness_core["witness_sha256"] = hashlib.sha256(encoding).hexdigest()
    witness_core["replay"] = {
        "proved": True,
        "hard_edge_count": hard_edge_count,
        "opposite_edge_count": opposite_edge_count,
        "selected_edge_count": hard_edge_count + opposite_edge_count,
        "raw_signed_degree_total": sum(signed),
    }
    return witness_core


def solve_model(
    geometry: Geometry,
    model: Any,
    handles: ModelHandles,
    metadata: dict[str, object],
    *,
    seconds: float,
    workers: int,
    random_seed: int,
    log_search_progress: bool,
) -> dict[str, object]:
    from ortools import __version__ as ortools_version
    from ortools.sat.python import cp_model

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = random_seed
    solver.parameters.randomize_search = False
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.log_search_progress = log_search_progress
    started = time.monotonic()
    status = solver.solve(model)
    elapsed = time.monotonic() - started
    status_name = solver.status_name(status)
    has_witness = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    fixed_edge_sign = str(metadata["fixed_edge_sign_shard"])
    mobius_top_origin = bool(metadata["mobius_top_origin_slice"])
    full_unsharded = fixed_edge_sign == "any" and not mobius_top_origin
    response = solver.response_proto
    return {
        "experiment": EXPERIMENT,
        "result_status": status_name,
        "model": metadata,
        "solver": {
            "name": "OR-Tools CP-SAT",
            "version": ortools_version,
            "status": status_name,
            "status_code": int(status),
            "max_time_seconds": seconds,
            "num_search_workers": workers,
            "random_seed": random_seed,
            "wall_time_seconds": solver.wall_time,
            "outer_elapsed_seconds": elapsed,
            "deterministic_time": float(response.deterministic_time),
            "branches": solver.num_branches,
            "conflicts": solver.num_conflicts,
        },
        "classification": {
            "model_infeasible": status == cp_model.INFEASIBLE,
            "full_endpoint_excluded": status == cp_model.INFEASIBLE and full_unsharded,
            "fixed_edge_sign_shard_infeasible": (
                status == cp_model.INFEASIBLE and fixed_edge_sign != "any"
            ),
            "localized_mobius_subfamily_only": mobius_top_origin,
            "necessary_projection_survives": has_witness,
            "incomplete": status == cp_model.UNKNOWN,
            "common_graph_constructed": False,
            "residual_ii_closed": False,
        },
        "witness": (
            replay_witness(
                geometry,
                solver,
                handles,
                mobius_top_origin=mobius_top_origin,
            )
            if has_witness
            else None
        ),
    }


def build_only_payload(
    metadata: dict[str, object], self_test: dict[str, int | bool], export_path: Path | None
) -> dict[str, object]:
    return {
        "experiment": EXPERIMENT,
        "result_status": "MODEL_BUILT_NOT_SOLVED",
        "model": metadata,
        "geometry_identity_self_test": self_test,
        "exported_model": str(export_path) if export_path is not None else None,
        "classification": {
            "model_infeasible": False,
            "full_endpoint_excluded": False,
            "fixed_edge_sign_shard_infeasible": False,
            "localized_mobius_subfamily_only": bool(
                metadata["mobius_top_origin_slice"]
            ),
            "necessary_projection_survives": False,
            "incomplete": True,
            "common_graph_constructed": False,
            "residual_ii_closed": False,
        },
        "witness": None,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=0)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument(
        "--fixed-edge-sign", choices=("any", "hard", "opposite"), default="any"
    )
    parser.add_argument(
        "--mobius-top-origin",
        action="store_true",
        help=(
            "restrict to the exact 16-half top construction's parity-forced "
            "origin incidence (a localized-Mobius subfamily only)"
        ),
    )
    parser.add_argument("--export-model", type=Path)
    parser.add_argument("--output", type=Path, default=Path("/tmp/resii_p31_signed_degree.json"))
    parser.add_argument("--log-search-progress", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.seconds <= 0:
        raise ValueError("--seconds must be positive")
    if args.workers <= 0:
        raise ValueError("--workers must be positive")
    geometry = build_geometry()
    self_test = geometry_identity_self_test(geometry)
    model, handles, metadata = build_model(
        geometry, args.fixed_edge_sign, args.mobius_top_origin
    )
    if args.export_model is not None:
        args.export_model.parent.mkdir(parents=True, exist_ok=True)
        if not model.export_to_file(str(args.export_model)):
            raise OSError(f"failed to export CP-SAT model to {args.export_model}")
    if args.build_only:
        payload = build_only_payload(metadata, self_test, args.export_model)
    else:
        payload = solve_model(
            geometry,
            model,
            handles,
            metadata,
            seconds=args.seconds,
            workers=args.workers,
            random_seed=args.random_seed,
            log_search_progress=args.log_search_progress,
        )
        payload["geometry_identity_self_test"] = self_test
        payload["exported_model"] = str(args.export_model) if args.export_model else None
    write_json_atomic(args.output, payload)
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
