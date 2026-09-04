#!/usr/bin/env python3
"""Incremental cut lower bound for p=31 transverse atom rows.

This module supplies a GPU-friendly necessary cost for the live centre
search.  It is strictly stronger than the retired scalar ``l1`` budget.

For a label subset ``A``, a compact atom ``K(a,b;c)=+ab-ac-bc`` has cut
weight either ``-2`` or ``0``.  An all-positive distinct-label triangle has
cut weight either ``0`` or ``2``.  Hence a row made from ``h`` positive
triangles and ``b`` compact atoms obeys

    -2*b <= x(delta(A)) <= 2*h                              (1)

for every ``A``.  The implementation retains all singleton and two-label
cuts.  A cell mutation ``x[u,v] += delta`` changes only two singleton cuts
and ``2*(p-2)`` two-label cuts; the cut cost therefore updates exactly in
``O(p)`` work.

The maximum interval violation is a rigorous lower bound on coefficient
``l1`` edits to any atom-decomposable row.  The additive hinge sum is the
incremental search cost.  Zero cost is necessary, not sufficient: the full
integer atom transport remains a separate exact check.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from itertools import combinations
from typing import Mapping

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    edge_radon_image,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import (
    paley_direction_sign,
    paley_edge_sign,
)
from e1_gmin_m4_p31_direct_mobius_parallel_design import (
    centered_physical_graph,
)


P = 31
LabelEdge = tuple[int, int]


def _label_edge(first: int, second: int, p: int) -> LabelEdge:
    if not 0 <= first < p or not 0 <= second < p or first == second:
        raise ValueError("a label edge needs two distinct labels in range")
    return tuple(sorted((first, second)))


def _interval_hinge(value: int, lower: int, upper: int) -> int:
    if value < lower:
        return lower - value
    if value > upper:
        return value - upper
    return 0


@dataclass(frozen=True)
class AtomRowSpec:
    """The two atom counts defining one transverse target row."""

    positive_triangles: int
    compact_atoms: int

    def __post_init__(self) -> None:
        for name, value in (
            ("positive_triangles", self.positive_triangles),
            ("compact_atoms", self.compact_atoms),
        ):
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0
            ):
                raise ValueError(f"{name} must be a nonnegative integer")

    @property
    def target_edge_sum(self) -> int:
        return 3 * self.positive_triangles - self.compact_atoms

    @property
    def positive_occurrence_budget(self) -> int:
        return 3 * self.positive_triangles + self.compact_atoms

    @property
    def negative_occurrence_budget(self) -> int:
        return 2 * self.compact_atoms

    @property
    def l1_budget(self) -> int:
        return 3 * (self.positive_triangles + self.compact_atoms)

    @property
    def cut_interval(self) -> tuple[int, int]:
        return -2 * self.compact_atoms, 2 * self.positive_triangles

    @classmethod
    def hard(cls, compact_atoms: int) -> "AtomRowSpec":
        return cls(positive_triangles=0, compact_atoms=compact_atoms)

    @classmethod
    def opposite(cls, compact_atoms: int) -> "AtomRowSpec":
        return cls(positive_triangles=6, compact_atoms=compact_atoms)


class IncrementalAtomRowBound:
    """Mutable exact state for the singleton/two-label cut cost.

    ``coefficients`` is the literal-star-subtracted row for a hard target
    and the full transverse row for an opposite target.
    """

    def __init__(
        self,
        coefficients: Mapping[LabelEdge, int],
        spec: AtomRowSpec,
        p: int = P,
    ) -> None:
        if not isinstance(p, int) or isinstance(p, bool) or p < 3:
            raise ValueError("p must be an integer at least three")
        self.p = p
        self.spec = spec
        self.coefficients = {
            (first, second): 0
            for first, second in combinations(range(p), 2)
        }
        for edge, value in coefficients.items():
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError("row coefficients must be integers")
            key = _label_edge(int(edge[0]), int(edge[1]), p)
            self.coefficients[key] += value

        self.degrees = [0] * p
        self.total = 0
        self.positive_mass = 0
        self.negative_mass = 0
        for (first, second), value in self.coefficients.items():
            self.degrees[first] += value
            self.degrees[second] += value
            self.total += value
            self.positive_mass += max(value, 0)
            self.negative_mass += max(-value, 0)
        self.l1 = self.positive_mass + self.negative_mass
        self.pair_cuts = {
            (first, second): (
                self.degrees[first]
                + self.degrees[second]
                - 2 * self.coefficients[(first, second)]
            )
            for first, second in combinations(range(p), 2)
        }
        lower, upper = spec.cut_interval
        self.cut_hinge_sum = sum(
            _interval_hinge(value, lower, upper)
            for value in self.degrees
        ) + sum(
            _interval_hinge(value, lower, upper)
            for value in self.pair_cuts.values()
        )
        self.odd_degree_count = sum(value % 2 != 0 for value in self.degrees)
        self._refresh_scalar_costs()

    def _budget_defects(self) -> dict[str, int]:
        return {
            "edge_sum": abs(self.total - self.spec.target_edge_sum),
            "positive_mass": max(
                0,
                self.positive_mass - self.spec.positive_occurrence_budget,
            ),
            "negative_mass": max(
                0,
                self.negative_mass - self.spec.negative_occurrence_budget,
            ),
            "l1": max(0, self.l1 - self.spec.l1_budget),
        }

    def _degree_projection(self) -> dict[str, object]:
        """Return the exact signed-degree projection for the two row types."""
        lower_cut, upper_cut = self.spec.cut_interval
        even = self.odd_degree_count == 0
        total_ok = self.total == self.spec.target_edge_sum
        singleton_range_ok = all(
            lower_cut <= value <= upper_cut for value in self.degrees
        )
        if self.spec.positive_triangles == 0:
            feasible = even and total_ok and singleton_range_ok
            return {
                "feasible": feasible,
                "type": "compact-only",
                "condition": (
                    "all degrees are nonpositive even and sum to -2b"
                ),
                "interval_gap": 0 if singleton_range_ok else 1,
            }

        if self.spec.positive_triangles != 6:
            return {
                "feasible": even and total_ok and singleton_range_ok,
                "type": "generic-cut-projection-only",
                "condition": "only the cut interval is asserted",
                "interval_gap": 0 if singleton_range_ok else 1,
            }
        if not even:
            return {
                "feasible": False,
                "type": "six-triangle-plus-compact",
                "condition": "degrees must first be even",
                "interval_gap": self.odd_degree_count // 2,
            }

        # Write degree(v)/2=A_v-C_v.  Six repeated distinct-label
        # triangles have occurrence degrees 0<=A_v<=6, sum A_v=18;
        # compact distinguished counts satisfy C_v>=0, sum C_v=b.
        half_degrees = [value // 2 for value in self.degrees]
        lower_a = [max(0, value) for value in half_degrees]
        upper_a = [
            min(6, value + self.spec.compact_atoms)
            for value in half_degrees
        ]
        point_gaps = [
            max(0, lower - upper)
            for lower, upper in zip(lower_a, upper_a, strict=True)
        ]
        lower_sum = sum(lower_a)
        upper_sum = sum(upper_a)
        sum_gap = max(lower_sum - 18, 18 - upper_sum, 0)
        interval_gap = sum(point_gaps) + sum_gap
        feasible = bool(total_ok and interval_gap == 0)
        return {
            "feasible": feasible,
            "type": "six-triangle-plus-compact",
            "condition": (
                "choose A_v in [max(0,d_v/2),min(6,b+d_v/2)] "
                "with sum A_v=18"
            ),
            "lower_occurrence_sum": lower_sum,
            "upper_occurrence_sum": upper_sum,
            "interval_gap": interval_gap,
        }

    def _refresh_scalar_costs(self) -> None:
        defects = self._budget_defects()
        projection = self._degree_projection()
        projection_gap = int(projection["interval_gap"])
        self.search_cost = (
            sum(defects.values())
            + self.odd_degree_count // 2
            + self.cut_hinge_sum
            + projection_gap
        )

    def apply_cell_delta(self, first: int, second: int, delta: int) -> None:
        """Apply one row-cell mutation and update the exact cost in O(p)."""
        if not isinstance(delta, int) or isinstance(delta, bool):
            raise ValueError("delta must be an integer")
        edge = _label_edge(first, second, self.p)
        if delta == 0:
            return
        first, second = edge
        lower, upper = self.spec.cut_interval

        # Remove the affected singleton and pair-cut hinge contributions.
        for vertex in (first, second):
            self.cut_hinge_sum -= _interval_hinge(
                self.degrees[vertex], lower, upper
            )
        affected_pairs = []
        for other in range(self.p):
            if other in (first, second):
                continue
            affected_pairs.append(_label_edge(first, other, self.p))
            affected_pairs.append(_label_edge(second, other, self.p))
        for pair in affected_pairs:
            self.cut_hinge_sum -= _interval_hinge(
                self.pair_cuts[pair], lower, upper
            )

        old_value = self.coefficients[edge]
        new_value = old_value + delta
        self.total += delta
        self.positive_mass += max(new_value, 0) - max(old_value, 0)
        self.negative_mass += max(-new_value, 0) - max(-old_value, 0)
        self.l1 = self.positive_mass + self.negative_mass
        self.coefficients[edge] = new_value
        self.odd_degree_count -= self.degrees[first] % 2 != 0
        self.odd_degree_count -= self.degrees[second] % 2 != 0
        self.degrees[first] += delta
        self.degrees[second] += delta
        self.odd_degree_count += self.degrees[first] % 2 != 0
        self.odd_degree_count += self.degrees[second] % 2 != 0

        # The cut of {first,second} is unchanged: the mutated cell is
        # internal and both degrees changed by delta.  Every other affected
        # two-label cut changes by exactly delta.
        for pair in affected_pairs:
            self.pair_cuts[pair] += delta
        expected_internal = (
            self.degrees[first]
            + self.degrees[second]
            - 2 * self.coefficients[edge]
        )
        if self.pair_cuts[edge] != expected_internal:
            raise ArithmeticError("the internal two-label cut changed")

        for vertex in (first, second):
            self.cut_hinge_sum += _interval_hinge(
                self.degrees[vertex], lower, upper
            )
        for pair in affected_pairs:
            self.cut_hinge_sum += _interval_hinge(
                self.pair_cuts[pair], lower, upper
            )
        self._refresh_scalar_costs()

    def edit_lower_bound(self) -> int:
        """Return a rigorous coefficient-l1 edit lower bound."""
        lower, upper = self.spec.cut_interval
        cut_max = max(
            (
                _interval_hinge(value, lower, upper)
                for value in (*self.degrees, *self.pair_cuts.values())
            ),
            default=0,
        )
        defects = self._budget_defects()
        return max(
            *defects.values(),
            self.odd_degree_count // 2,
            cut_max,
        )

    def summary(self) -> dict[str, object]:
        lower, upper = self.spec.cut_interval
        singleton_violations = sum(
            _interval_hinge(value, lower, upper) > 0
            for value in self.degrees
        )
        pair_violations = sum(
            _interval_hinge(value, lower, upper) > 0
            for value in self.pair_cuts.values()
        )
        return {
            "p": self.p,
            "positive_triangles": self.spec.positive_triangles,
            "compact_atoms": self.spec.compact_atoms,
            "cut_interval": [lower, upper],
            "cut_bank_size": self.p + len(self.pair_cuts),
            "total": self.total,
            "positive_mass": self.positive_mass,
            "negative_mass": self.negative_mass,
            "l1": self.l1,
            "budget_defects": self._budget_defects(),
            "odd_degree_count": self.odd_degree_count,
            "degree_projection": self._degree_projection(),
            "singleton_cut_violations": singleton_violations,
            "two_label_cut_violations": pair_violations,
            "cut_hinge_sum": self.cut_hinge_sum,
            "incremental_search_cost": self.search_cost,
            "coefficient_l1_edit_lower_bound": self.edit_lower_bound(),
            "cell_delta_affected_singletons": 2,
            "cell_delta_affected_two_label_cuts": 2 * (self.p - 2),
            "cell_delta_update_complexity": "O(p)",
            "zero_cost_is_sufficient_for_atom_decomposition": False,
        }


def strict_over_scalar_l1_witness() -> dict[str, object]:
    """Return a top hard-row vector caught only by structural conditions.

    The six-edge seed has exact compact-only masses and nonpositive even
    degrees but cut({0,1})=2.  Nine copies of a compact atom on labels
    6,7,8 pad it to the live top count e=11 without changing that cut.
    """
    coefficients: dict[LabelEdge, int] = {}
    for edge in ((0, 2), (1, 3)):
        coefficients[edge] = coefficients.get(edge, 0) + 1
    for edge in ((0, 1), (2, 4), (4, 5), (3, 5)):
        coefficients[edge] = coefficients.get(edge, 0) - 1
    coefficients[(6, 7)] = 9
    coefficients[(6, 8)] = -9
    coefficients[(7, 8)] = -9
    state = IncrementalAtomRowBound(coefficients, AtomRowSpec.hard(11))
    summary = state.summary()
    pair_cut = state.pair_cuts[(0, 1)]
    scalar_pass = bool(
        summary["budget_defects"]
        == {"edge_sum": 0, "positive_mass": 0, "negative_mass": 0, "l1": 0}
        and summary["degree_projection"]["feasible"]
    )
    proved = bool(
        scalar_pass
        and pair_cut == 2
        and summary["two_label_cut_violations"] >= 1
        and summary["coefficient_l1_edit_lower_bound"] >= 2
    )
    if not proved:
        raise ArithmeticError("the strict two-label-cut witness changed")
    return {
        "p": P,
        "compact_atoms": 11,
        "scalar_l1_and_mass_budgets_pass": scalar_pass,
        "signed_degree_projection_passes": True,
        "violated_subset": [0, 1],
        "violated_cut_value": pair_cut,
        "required_cut_interval": [-22, 0],
        "coefficient_l1_edit_lower_bound": summary[
            "coefficient_l1_edit_lower_bound"
        ],
        "proved": proved,
    }
def normalized_transverse_rows_from_graph() -> tuple[
    str, tuple[tuple[int, AtomRowSpec, dict[LabelEdge, int]], ...]
]:
    exposed = centered_physical_graph()
    directions = projective_functionals(P)
    signs = tuple(paley_direction_sign(P, row) for row in directions)
    centers = {
        int(record["target_direction_index"]): int(
            record["canonical_target_center"]
        )
        for record in exposed["hard_target_centers"]
    }
    source = {
        edge: paley_edge_sign(P, edge) for edge in exposed["edges"]
    }
    image = edge_radon_image(P, source)
    rows = []
    for direction_index, sign in enumerate(signs):
        quota = sign * image.get(("P", direction_index), 0)
        coefficients = {
            (left, right): sign
            * image.get(("K", direction_index, left, right), 0)
            for left, right in combinations(range(P), 2)
        }
        coefficients = {
            edge: value for edge, value in coefficients.items() if value
        }
        if sign == 1:
            center = centers[direction_index]
            for other in range(P):
                if other != center:
                    edge = _label_edge(center, other, P)
                    coefficients[edge] = coefficients.get(edge, 0) + 1
                    if coefficients[edge] == 0:
                        del coefficients[edge]
            spec = AtomRowSpec.hard(quota - 3)
        else:
            spec = AtomRowSpec.opposite(quota - 9)
        rows.append((direction_index, spec, coefficients))
    return str(exposed["graph_sha256"]), tuple(rows)


def frozen_graph_incremental_bound_certificate() -> dict[str, object]:
    """Evaluate the new cost through the exact public 479-edge graph API."""
    graph_sha256, raw_rows = normalized_transverse_rows_from_graph()
    rows = []
    for direction_index, spec, coefficients in raw_rows:
        state = IncrementalAtomRowBound(coefficients, spec)
        row = state.summary()
        row["direction_index"] = direction_index
        rows.append(row)
    hard_rows = [row for row in rows if row["positive_triangles"] == 0]
    opposite_rows = [row for row in rows if row["positive_triangles"] == 6]
    proved = bool(
        graph_sha256
        == "c0b32bdf228401ba5ffe68be543b9e6fddb31f86594ff953e1d290a6faeeae0d"
        and len(rows) == 32
        and len(hard_rows) == len(opposite_rows) == 16
        and all(row["cut_bank_size"] == 496 for row in rows)
        and all(row["coefficient_l1_edit_lower_bound"] > 0 for row in rows)
        and strict_over_scalar_l1_witness()["proved"]
    )
    if not proved:
        raise ArithmeticError("the frozen graph incremental bound changed")
    return {
        "p": P,
        "graph_sha256": graph_sha256,
        "row_count": len(rows),
        "hard_row_count": len(hard_rows),
        "opposite_row_count": len(opposite_rows),
        "cut_bank_per_row": 496,
        "minimum_edit_lower_bound": min(
            int(row["coefficient_l1_edit_lower_bound"]) for row in rows
        ),
        "maximum_edit_lower_bound": max(
            int(row["coefficient_l1_edit_lower_bound"]) for row in rows
        ),
        "rows_with_two_label_cut_violations": sum(
            int(row["two_label_cut_violations"]) > 0 for row in rows
        ),
        "rows": tuple(rows),
        "strictly_stronger_than_scalar_l1": strict_over_scalar_l1_witness(),
        "gpu_state": (
            "per row keep 465 cells, 31 degrees, 465 pair cuts, masses, "
            "odd-degree count, and additive hinge sum"
        ),
        "single_cell_update": (
            "2 singleton cuts plus 58 two-label cuts at p=31"
        ),
        "full_atom_decomposition_certified": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    return {
        "title": "p31 incremental singleton/two-label atom-cut bound",
        "compact_atom_cut_values": [-2, 0],
        "positive_triangle_cut_values": [0, 2],
        "row_cut_interval": "[-2b,2h] for h positive triangles and b compact atoms",
        "strict_witness": strict_over_scalar_l1_witness(),
        "frozen_graph": frozen_graph_incremental_bound_certificate(),
        "incremental_cost_exactly_updated_in_O_p": True,
        "coefficient_edit_lower_bound_rigorous": True,
        "full_atom_transport_still_required_at_zero_cost": True,
        "residual_ii_closed": False,
        "proved": True,
    }


def main() -> dict[str, object]:
    result = theorem_record()
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
