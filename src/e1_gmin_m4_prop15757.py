#!/usr/bin/env python3
r"""Prop. 15.757 -- exact binary edge-Radon image and compact-triangle barrier.

Let ``V=F_p^2`` for an odd prime ``p``.  For every affine direction ``L``
push an edge set on ``V`` to the graph on the ``p`` fibres of ``L``: retain
the parity of every off-diagonal fibre-pair count and, in one extra
coordinate, the parity of the total number of edges parallel to ``L``.
Call the direct sum of these ``p+1`` maps ``R_2``.

The exact image dimension is

    rank_F2(R_2) = p^2-1 + (p+1) binom(p-1,2).             (1)

Moreover, its image is exactly the data satisfying the following familiar
conditions.

* The row parities of the off-diagonal fibre graphs are the affine Radon
  transforms of one even point word (the graph boundary).
* In every direction, parallel parity plus off-diagonal edge parity is one
  common bit, and the sum of the parallel parities equals that bit.

There are no further binary cross-direction equations.  To prove (1), split
an off-diagonal ``K_p`` edge space into its row-boundary quotient and its
Eulerian cycle space.  Proposition 15.692 identifies the direct sum of all
row-boundary quotients with the even point space, of dimension ``p^2-1``.
Over an algebraic closure of ``F_2``, a basis of the cycle space in direction
``L`` is

    zeta^(a L(s)+b L(t)) + zeta^(b L(s)+a L(t)),

indexed by unordered distinct ``a,b in F_p^*``.  Its size is
``binom(p-1,2)``.  Different projective directions have disjoint nonzero
Fourier support in ``(V^*)^2``, so all ``p+1`` cycle spaces are independent.
The displayed total-parity identities account for the remaining ``p+1``
codimensions.

This theorem matters for the large residual-(ii) aggregate survivor.  The
sharp scaled-mass ``p+1`` atom need not be an oriented mixed pair with a
large canonical coefficient norm.  For distinct ``a,b,c``,

    B=x_a x_b+x_c-x_a x_c-x_b x_c

is Boolean on the full cube and

    4B=1+z_a z_b-z_a z_c-z_b z_c,     4p E_J[B]=p+1.     (2)

Thus its coefficient graph is a compact signed triangle and is Eulerian
modulo two.  A hard row ``x_j+2 sum B`` has precisely the required odd-fibre
parity: every fibre except ``j`` is odd.  In an opposite row, the omitted-
pair atom, all-equal-triple atoms, and the compact atoms are all Eulerian.
Their parallel/off-diagonal parities also give the common odd edge total.
There are ``m=(p+1)/2`` hard and ``m`` opposite directions; ``I=0`` is the
parallel count at infinity, not an omitted affine direction.  Boundary
compatibility is nevertheless automatic.  Over an algebraic closure of
``F_2``, choose a primitive ``p``-th root ``omega`` and prescribe

    fhat(0)=0,
    fhat(lambda*n_L)=omega^(lambda*j_L)  (L hard, lambda != 0),
    fhat(lambda*n_L)=0                   (L opposite).

Its inverse Fourier transform is ``F_2``-valued by Frobenius and has hard
Radon rows equal to one on every fibre except ``j_L``, and zero opposite
rows.  Finally the exact parallel-count sums over all ``p+1`` directions
equal ``|H|``, so the total-parity conditions also hold.

Consequently every such compact-triangle local survivor has a binary edge
lift.  This is a proved method barrier, not an integer or simple-graph lift:
nonnegativity, exact signed multiplicities, and ``0/1`` edge capacity remain
uncontrolled.  In particular it does not close residual (ii).
"""
from __future__ import annotations

import json
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15751 import atomic_write_json


ROOT = Path(__file__).resolve().parents[1]


def _check_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")


def edge_radon_dimensions(p: int) -> dict[str, int | bool | str]:
    """Return the exact domain, target, rank, and compatibility codimension."""
    _check_prime(p)
    point_count = p * p
    domain = comb(point_count, 2)
    outputs_per_direction = comb(p, 2) + 1
    target = (p + 1) * outputs_per_direction
    boundary_rank = point_count - 1
    cycle_rank_per_direction = comb(p - 1, 2)
    rank = boundary_rank + (p + 1) * cycle_rank_per_direction
    codimension = target - rank
    proved = codimension == p + 1
    if not proved:
        raise ArithmeticError("the edge-Radon dimension identity changed")
    return {
        "p": p,
        "affine_points": point_count,
        "simple_edge_variables": domain,
        "directions": p + 1,
        "outputs_per_direction": outputs_per_direction,
        "target_dimension": target,
        "boundary_component_dimension": boundary_rank,
        "cycle_component_per_direction": cycle_rank_per_direction,
        "image_rank": rank,
        "compatibility_codimension": codimension,
        "compatibility_equations": (
            "p equal-total equations plus one sum-of-parallel-parities equation"
        ),
        "proved": proved,
    }


def _rank_binary_rows(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def exact_edge_radon_rank(p: int) -> dict[str, int | bool]:
    """Build ``R_2`` coefficientwise and row-reduce it exactly over ``F_2``.

    This is intended as a small-prime fail-when-wrong replay.  The general
    proof is the boundary/cycle Fourier decomposition in the module
    docstring, not extrapolation from these ranks.
    """
    _check_prime(p)
    points = tuple(product(range(p), repeat=2))
    edges = tuple(combinations(range(len(points)), 2))
    directions = tuple((1, slope) for slope in range(p)) + ((0, 1),)
    rows: list[int] = []
    for left, right in directions:
        labels = tuple((left * x + right * y) % p for x, y in points)
        for s, t in combinations(range(p), 2):
            bits = 0
            for index, (u, v) in enumerate(edges):
                if tuple(sorted((labels[u], labels[v]))) == (s, t):
                    bits |= 1 << index
            rows.append(bits)
        parallel = 0
        for index, (u, v) in enumerate(edges):
            if labels[u] == labels[v]:
                parallel |= 1 << index
        rows.append(parallel)

    actual = _rank_binary_rows(rows)
    expected = int(edge_radon_dimensions(p)["image_rank"])
    proved = actual == expected
    if not proved:
        raise ArithmeticError("the exact binary edge-Radon rank changed")
    return {
        "p": p,
        "matrix_rows": len(rows),
        "matrix_columns": len(edges),
        "exact_rank_over_F2": actual,
        "formula_rank": expected,
        "proved": proved,
    }


def compact_triangle_atom(p: int) -> dict[str, object]:
    """Audit the Boolean ``p+1`` atom and its three signed coefficients."""
    _check_prime(p)
    m = (p + 1) // 2
    support = comb(p - 3, m - 1) + comb(p - 3, m - 2)
    density = Fraction(support, comb(p, m))
    scaled_mass = 4 * p * density
    values = {}
    for xa, xb, xc in product((0, 1), repeat=3):
        value = xa * xb + xc - xa * xc - xb * xc
        values[f"{xa}{xb}{xc}"] = value
    proved = bool(
        set(values.values()) == {0, 1}
        and sum(values.values()) == 2
        and scaled_mass == p + 1
    )
    if not proved:
        raise ArithmeticError("the compact signed-triangle atom changed")
    return {
        "p": p,
        "formula": "B=x_a*x_b+x_c-x_a*x_c-x_b*x_c",
        "z_formula": "4B=1+z_a*z_b-z_a*z_c-z_b*z_c",
        "cube_values": values,
        "middle_slice_support": support,
        "middle_slice_density": str(density),
        "scaled_mass_4pE": int(scaled_mass),
        "coefficient_l1": 3,
        "coefficient_graph_mod_two": "triangle, hence Eulerian",
        "offset": 1,
        "proved": proved,
    }


def sharp_atom_masses(p: int) -> dict[str, object]:
    """Check the ``p-3`` omitted/all-equal atoms used in opposite rows."""
    _check_prime(p)
    m = (p + 1) // 2
    omitted_density = Fraction(comb(p - 2, m), comb(p, m))
    equal_density = Fraction(
        comb(p - 3, m) + comb(p - 3, m - 3), comb(p, m)
    )
    expected = Fraction(p - 3, 4 * p)
    proved = omitted_density == equal_density == expected
    if not proved:
        raise ArithmeticError("the sharp p-3 atom masses changed")
    return {
        "p": p,
        "omitted_pair_density": str(omitted_density),
        "all_equal_triple_density": str(equal_density),
        "scaled_mass_4pE": p - 3,
        "both_coefficient_graphs_mod_two_are_Eulerian": True,
        "proved": proved,
    }


def compact_survivor_parity_ledger(p: int, t: int) -> dict[str, object]:
    """Verify every parity identity of the ``p=4r+1`` aggregate survivor."""
    _check_prime(p)
    if p % 4 != 1:
        raise ValueError("the compact survivor ledger uses p congruent to 1 mod 4")
    if not isinstance(t, int) or isinstance(t, bool) or t < 0:
        raise ValueError("t must be a nonnegative integer")
    r = (p - 1) // 4
    m = 2 * r + 1
    lower_t = 2 * r * r - 5 * r
    upper_t = 4 * r * r - 6 * r - 3
    if p < 13 or not lower_t <= t <= upper_t:
        raise ValueError(
            f"need the p=4r+1 local-survivor interval {lower_t}<=t<={upper_t}"
        )
    h_edges = 4 * p + 2 * t + 1
    # One balanced choice suffices to audit that the aggregate budgets admit
    # nonnegative integer row parameters.  The parity conclusions below use
    # only their sums and therefore hold for every such choice.
    e_values = [t // m + (index < t % m) for index in range(m)]
    q_total = 6 * r + t
    Q_values = [
        r + (q_total - r * m) // m
        + (index < (q_total - r * m) % m)
        for index in range(m)
    ]
    hard_edge_parities = [e % 2 for e in e_values]
    hard_parallel_parities = [(5 + e) % 2 for e in e_values]
    opposite_edge_parities = [(Q + 1) % 2 for Q in Q_values]
    opposite_parallel_parities = [Q % 2 for Q in Q_values]
    finite_parallel_total = sum(5 + e for e in e_values) + sum(Q_values)
    opposite_mass_checks = [
        (r - 1) * (p - 3) + (Q - r) * (p + 1)
        == (p + 1) * Q - 2 * p + 4
        for Q in Q_values
    ]
    opposite_offset_checks = [
        3 + (-1 + (r - 2) + (Q - r)) == Q for Q in Q_values
    ]
    proved = bool(
        sum(e_values) == t
        and sum(Q_values) == q_total
        and len(e_values) == len(Q_values) == m
        and min(Q_values) >= r
        and finite_parallel_total == h_edges
        and all((edge + parallel) % 2 == 1 for edge, parallel in zip(hard_edge_parities, hard_parallel_parities))
        and all((edge + parallel) % 2 == 1 for edge, parallel in zip(opposite_edge_parities, opposite_parallel_parities))
        and all(opposite_mass_checks)
        and all(opposite_offset_checks)
        and compact_triangle_atom(p)["proved"]
        and sharp_atom_masses(p)["proved"]
    )
    if not proved:
        raise ArithmeticError("the compact survivor parity ledger changed")
    return {
        "p": p,
        "r": r,
        "t": t,
        "hard_direction_count": m,
        "opposite_direction_count": m,
        "all_affine_direction_count": 2 * m,
        "isolated_chart_infinity_parallel_count_I": 0,
        "H_edge_count": h_edges,
        "hard_excess_values": e_values,
        "hard_parallel_counts": [5 + e for e in e_values],
        "opposite_parallel_counts": Q_values,
        "all_direction_parallel_counts_sum_to_H_edges": True,
        "hard_row_parity": "all fibres except the literal fibre j are odd",
        "opposite_row_parity": "Eulerian, so every fibre is even",
        "every_recorded_direction_has_odd_total_edge_parity": True,
        "boundary_fourier_completion": (
            "fhat(0)=0; on each hard dual line fhat(lambda*n_L)="
            "omega^(lambda*j_L); on each opposite dual line fhat=0; "
            "Frobenius makes the inverse transform F2-valued"
        ),
        "opposite_mass_identity": "(r-1)(p-3)+(Q-r)(p+1)=(p+1)Q-2p+4",
        "opposite_offset_identity": "3+[-1+(r-2)+(Q-r)]=Q",
        "binary_edge_lift_exists_by_exact_image_theorem": True,
        "integer_nonnegative_simple_edge_lift_proved": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    ranks = {str(p): exact_edge_radon_rank(p) for p in (3, 5, 7)}
    dimensions = {str(p): edge_radon_dimensions(p) for p in (3, 5, 7, 11, 13)}
    atoms = {str(p): compact_triangle_atom(p) for p in (5, 13, 17, 29)}
    survivors = {
        "p13_t5": compact_survivor_parity_ledger(13, 5),
        "p17_t20": compact_survivor_parity_ledger(17, 20),
        "p29_t70": compact_survivor_parity_ledger(29, 70),
        "p37_t117": compact_survivor_parity_ledger(37, 117),
    }
    proved = bool(
        all(row["proved"] for row in ranks.values())
        and all(row["proved"] for row in dimensions.values())
        and all(row["proved"] for row in atoms.values())
        and all(row["proved"] for row in survivors.values())
    )
    return {
        "prop": "15.757",
        "title": "Exact binary edge-Radon image and compact-triangle barrier",
        "status": "PROVED THEOREM AND PROVED METHOD BARRIER",
        "proved": {
            "edge_Radon_image_dimension_all_odd_primes": proved,
            "boundary_and_total_parity_conditions_are_complete": proved,
            "compact_triangle_scaled_mass_p_plus_1": proved,
            "aggregate_survivor_has_binary_edge_lift": proved,
            "integer_nonnegative_simple_edge_lift": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "rank_formula": "p^2-1+(p+1)*binom(p-1,2)",
        "exact_rank_replays": ranks,
        "dimension_ledgers": dimensions,
        "compact_atom_checks": atoms,
        "compact_survivor_checks": survivors,
        "remaining_obstruction": (
            "Find an integer, nonnegative, simple-edge invariant beyond the "
            "complete F2 boundary/cycle image."
        ),
        "duplicate_work_guards": [
            "Do not use the repeated oriented-pair coefficient l1 as a universal p+1-atom bound.",
            "Do not search for another F2 cross-direction equation after boundary and total parity.",
            "Do not promote a binary edge lift to a nonnegative 0/1 graph realization.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15757.json"
    atomic_write_json(destination, out)
    print("Prop. 15.757 binary edge-Radon image: proved")
    print("  integer/simple-edge lift: OPEN")
    print(f"  wrote {destination}")
    return out


if __name__ == "__main__":
    main()
