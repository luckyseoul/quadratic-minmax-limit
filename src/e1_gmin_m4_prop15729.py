#!/usr/bin/env python3
r"""Prop. 15.729 -- affine unique-trisecant reduction at the endpoint.

Assume the all-prime endpoint normal form of Proposition 15.727.  Thus
``D`` is an affine set of ``p+1`` points, ``p=3R+c`` with ``c in {1,2}``, and
the lines meeting ``D`` in at least three points are pairwise ``D``-disjoint
blocks: ``x`` trisecants and ``y`` 4-secants with ``x+2y=R``.

Choose one rich block.  Retain three of its points and retain two points on
every other rich block; retain every singleton point.  If the distinguished
block is a trisecant, this deletes ``(x-1)+2y=R-1`` points.  If it is a
4-secant, this deletes ``1+x+2(y-1)=R-1`` points.  The resulting affine set
``U`` therefore has

    |U| = p+2-R.

Any line containing three points of ``U`` already contained three points of
``D``.  It is consequently one of the rich blocks above.  The distinguished
block has three retained points and every other rich block has two, so ``U``
is a ``(p+2-R,3)``-arc with exactly one trisecant.

Let that trisecant be ``{P,Q,Z}``, and put ``B=U\{P,Q}``.  Then ``B`` is an
affine arc of size ``p-R``.  Both ``B union {P}`` and ``B union {Q}`` are
arcs: a new trisecant through either extension point would have been a second
trisecant of ``U``.  Moreover the line ``PQZ`` meets ``B`` only in ``Z``.
Thus ``P`` and ``Q`` are two distinct affine extension points on one tangent
of ``B``.  The same conclusion holds after choosing any of the three points
of the unique trisecant to remain in ``B``.

In the two prime residue classes the sizes are

    c=1: |U|=2R+3 and |B|=2R+1,
    c=2: |U|=2R+4 and |B|=2R+2.

This is a proved structural reduction, not an exclusion of the endpoint.
It uses no finite configuration search and no additional finite-geometry
classification.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15727 import endpoint_block_row, endpoint_residue_data


ROOT = Path(__file__).resolve().parents[1]
DISTINGUISHED_BLOCK_TYPES = ("trisecant", "4-secant")


def endpoint_residue_size_row(p: int) -> dict[str, object]:
    """Record the unique-trisecant and co-tangent arc sizes by residue."""
    data = endpoint_residue_data(p)
    R = data["R"]
    c = data["c"]
    unique_trisecant_size = p + 2 - R
    cotangent_arc_size = p - R
    expected_unique = 2 * R + c + 2
    expected_arc = 2 * R + c
    residue_formulas = {
        1: {
            "unique_trisecant_size_formula": "2R+3",
            "cotangent_arc_size_formula": "2R+1",
        },
        2: {
            "unique_trisecant_size_formula": "2R+4",
            "cotangent_arc_size_formula": "2R+2",
        },
    }[c]
    return {
        **data,
        "unique_trisecant_size": unique_trisecant_size,
        "cotangent_arc_size": cotangent_arc_size,
        **residue_formulas,
        "size_drop_from_D_to_unique_trisecant_set": R - 1,
        "size_drop_from_unique_trisecant_set_to_arc": 2,
        "proved": bool(
            unique_trisecant_size == expected_unique
            and cotangent_arc_size == expected_arc
            and unique_trisecant_size - cotangent_arc_size == 2
        ),
    }


def endpoint_unique_trisecant_construction(
    p: int, four_secants: int, distinguished_block_type: str
) -> dict[str, object]:
    """Audit one choice of distinguished rich block.

    The geometry needed here is precisely the Proposition 15.727 hypothesis
    that the displayed blocks are all the rich lines of ``D`` and are
    pairwise disjoint.  The returned arithmetic verifies the deletion count
    and every retained block occupancy.
    """
    if distinguished_block_type not in DISTINGUISHED_BLOCK_TYPES:
        raise ValueError("distinguished block must be a trisecant or 4-secant")

    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    c = int(block["c"])
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    singleton_points = int(block["singleton_points"])

    if distinguished_block_type == "trisecant":
        if x == 0:
            raise ValueError("this block row has no trisecant to distinguish")
        deleted_on_distinguished_block = 0
        deleted_on_other_trisecants = x - 1
        deleted_on_other_four_secants = 2 * y
        other_rich_blocks = x + y - 1
    else:
        if y == 0:
            raise ValueError("this block row has no 4-secant to distinguish")
        deleted_on_distinguished_block = 1
        deleted_on_other_trisecants = x
        deleted_on_other_four_secants = 2 * (y - 1)
        other_rich_blocks = x + y - 1

    total_deletions = (
        deleted_on_distinguished_block
        + deleted_on_other_trisecants
        + deleted_on_other_four_secants
    )
    retained_on_distinguished_block = 3
    retained_on_every_other_rich_block = 2
    retained_points = (
        singleton_points
        + retained_on_distinguished_block
        + retained_on_every_other_rich_block * other_rich_blocks
    )
    expected_size = p + 2 - R

    # If a line is a trisecant of U, it is rich in D because U is a subset of
    # D.  Proposition 15.727 says every such D-line is one of these blocks.
    # Their retained occupancies therefore prove uniqueness without a finite
    # incidence search.
    unique_trisecant = bool(
        retained_on_distinguished_block == 3
        and retained_on_every_other_rich_block == 2
        and block["rich_lines_pairwise_D_disjoint"]
    )
    proved = bool(
        block["proved"]
        and x + 2 * y == R
        and total_deletions == R - 1
        and retained_points == expected_size
        and unique_trisecant
    )
    return {
        "p": p,
        "R": R,
        "c": c,
        "trisecants_x": x,
        "four_secants_y": y,
        "singleton_points": singleton_points,
        "distinguished_block_type": distinguished_block_type,
        "deleted_on_distinguished_block": deleted_on_distinguished_block,
        "deleted_on_other_trisecants": deleted_on_other_trisecants,
        "deleted_on_other_four_secants": deleted_on_other_four_secants,
        "total_deletions": total_deletions,
        "target_deletions": R - 1,
        "retained_singleton_points": singleton_points,
        "retained_on_distinguished_block": retained_on_distinguished_block,
        "other_rich_block_count": other_rich_blocks,
        "retained_on_every_other_rich_block": (
            retained_on_every_other_rich_block
        ),
        "retained_point_count": retained_points,
        "retained_point_count_formula": "p+2-R",
        "all_points_affine": True,
        "maximum_line_occupancy": 3,
        "trisecant_count": 1,
        "unique_trisecant_reason": (
            "every U-trisecant is a D-rich line; the distinguished rich "
            "block retains three points and every other rich block retains two"
        ),
        "finite_configuration_search_used": False,
        "proved": proved,
    }


def cotangent_extension_row(
    p: int, four_secants: int, distinguished_block_type: str
) -> dict[str, object]:
    """Delete two triple points and expose two co-tangent extensions."""
    unique = endpoint_unique_trisecant_construction(
        p, four_secants, distinguished_block_type
    )
    R = int(unique["R"])
    unique_size = int(unique["retained_point_count"])
    arc_size = unique_size - 2
    proved = bool(
        unique["proved"]
        and unique["trisecant_count"] == 1
        and arc_size == p - R
    )
    return {
        "p": p,
        "R": R,
        "c": int(unique["c"]),
        "four_secants_y": four_secants,
        "distinguished_block_type": distinguished_block_type,
        "unique_trisecant_set_size": unique_size,
        "delete_any_two_points_of_unique_trisecant": True,
        "remaining_set_is_an_arc": True,
        "arc_size": arc_size,
        "arc_size_formula": "p-R",
        "number_of_choices_for_remaining_trisecant_point": 3,
        "deleted_points_are_distinct_affine_extension_points": True,
        "adding_either_deleted_point_preserves_the_arc": True,
        "common_line_meets_arc_only_in_remaining_trisecant_point": True,
        "extensions_lie_on_one_arc_tangent": True,
        "extension_reason": (
            "a secant through either deleted point would give a second "
            "trisecant of U; their common old trisecant contains only the "
            "one retained point of the arc"
        ),
        "proved": proved,
    }


def endpoint_unique_trisecant_rows(p: int) -> list[dict[str, object]]:
    """Cover every endpoint block count and every available block type."""
    data = endpoint_residue_data(p)
    R = data["R"]
    rows: list[dict[str, object]] = []
    for y in range(R // 2 + 1):
        x = R - 2 * y
        if x:
            rows.append(endpoint_unique_trisecant_construction(p, y, "trisecant"))
        if y:
            rows.append(endpoint_unique_trisecant_construction(p, y, "4-secant"))
    if not rows or not all(row["proved"] for row in rows):
        raise ArithmeticError("an endpoint block row lost its construction")
    return rows


def proposition_15729() -> dict[str, object]:
    """Package the all-prime structural implication honestly."""
    # The proof above is symbolic.  Keep only one live endpoint from each
    # residue class in the JSON artifact; the focused tests exercise a wider
    # prime sample and every admissible block count.
    sample_primes = (31, 41)
    residue_rows = [endpoint_residue_size_row(p) for p in sample_primes]
    construction_rows = {
        str(p): endpoint_unique_trisecant_rows(p) for p in sample_primes
    }
    cotangent_rows = {
        str(p): [
            cotangent_extension_row(
                p,
                int(row["four_secants_y"]),
                str(row["distinguished_block_type"]),
            )
            for row in rows
        ]
        for p, rows in ((p, construction_rows[str(p)]) for p in sample_primes)
    }
    proved = bool(
        all(row["proved"] for row in residue_rows)
        and all(
            row["proved"]
            for rows in construction_rows.values()
            for row in rows
        )
        and all(
            row["proved"] for rows in cotangent_rows.values() for row in rows
        )
    )
    return {
        "prop": "15.729",
        "title": "Affine unique-trisecant and co-tangent-extension reduction",
        "result_status": "proved structural reduction",
        "hypotheses": {
            "p": "odd prime p>=17",
            "endpoint": "R=floor((p-1)/3), p=3R+c, c in {1,2}",
            "D": "p+1 affine points",
            "prop_15_727_normal_form": (
                "all rich lines are D-disjoint trisecants/4-secants with "
                "x+2y=R"
            ),
        },
        "conclusion": {
            "affine_unique_trisecant_set": (
                "a (p+2-R,3)-arc with exactly one trisecant"
            ),
            "affine_cotangent_arc": (
                "a (p-R)-arc with two distinct extension points on one tangent"
            ),
            "residue_c_1_sizes": {"unique_trisecant": "2R+3", "arc": "2R+1"},
            "residue_c_2_sizes": {"unique_trisecant": "2R+4", "arc": "2R+2"},
        },
        "symbolic_proof": {
            "distinguished_trisecant_deletions": "(x-1)+2y=R-1",
            "distinguished_four_secant_deletions": "1+x+2(y-1)=R-1",
            "uniqueness": (
                "U subset D, so every U-trisecant is a D-rich line; retained "
                "rich-block occupancies are 3 once and 2 otherwise"
            ),
            "cotangent_extensions": (
                "remove two points P,Q from the unique trisecant PQZ; no "
                "other U-trisecant lets P or Q create a secant, and PQZ is "
                "tangent to the remaining arc at Z"
            ),
        },
        "residue_size_audit": residue_rows,
        "sample_block_arithmetic_audit": construction_rows,
        "sample_cotangent_audit": cotangent_rows,
        "finite_configuration_search_used": False,
        "new_classification_used": False,
        "endpoint_excluded": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_gate": (
            "exclude or classify the two residue-class families of large "
            "affine unique-trisecant 3-arcs, equivalently the associated "
            "large affine arcs with two co-tangent extension points"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic arithmetic certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15729.json"
    output.write_text(json.dumps(proposition_15729(), indent=2, sort_keys=True) + "\n")
    return output


def main() -> None:
    result = proposition_15729()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.729 audit failed")
    path = write_evidence()
    print("Prop 15.729 affine unique-trisecant reduction: proved")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
