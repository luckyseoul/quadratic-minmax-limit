#!/usr/bin/env python3
r"""Prop. 15.759 -- the complete low-degree p-torsion edge-Radon hierarchy.

Fix an odd prime ``p`` and use all ``p+1`` affine projections

    L_lambda(x,y)=x+lambda*y,                 lambda in F_p,
    L_infinity(x,y)=y.

For an integral edge set on ``F_p^2``, let ``K_lambda(s,t)`` be the number
of edges whose distinct endpoint labels are ``s,t`` and let ``P_lambda``
be the total number of edges parallel to ``L_lambda``.  Besides the
``p`` equal-total equations among the ``p+1`` directional totals and the
one parallel-sum equation, reduction modulo ``p`` has the following
relations.  For

    2 <= d <= p-2,     0 <= k < floor(d/2),
    Q_(d,k)(s,t)=(s-t)^2*(s*t)^k*(s+t)^(d-2-2k),

and ``0 <= j <= p-d-2``, one has

    sum_lambda lambda^j sum_(s<t) K_lambda(s,t) Q_(d,k)(s,t)=0 mod p. (1)

Indeed, the contribution of one edge to (1) is the sum over ``lambda``
of a polynomial of degree at most ``d+j<=p-2``.  Every such power sum is
zero in ``F_p``.  Projective completion supplies one more equation,

    sum_lambda lambda^(p-d-1) C_(d,k)(lambda)
       +C_(d,k)(infinity)=0,                               (2)

where ``C_(d,k)(L)=sum K_L(s,t)Q_(d,k)(s,t)``.  The
``Q_(d,k)`` form a basis of the homogeneous symmetric degree-``d``
polynomials divisible by ``(s-t)^2``.

At degree ``d=p-1``, equation (2) still applies.  This gives ``m``
relations, where ``p=2m+1``.  Their one-dimensional intersection with the
ordinary row space is spanned by ``(s-t)^(p-1)``: it is one on every
off-diagonal pair and zero on the diagonal, so its projective sum is the
sum of the off-diagonal totals, namely ``p`` times the common edge total.
The displayed ``Q_(p-1,k)`` basis does not literally contain that vector.
Writing ``D=(s-t)^2=(s+t)^2-4st`` gives

    D^m=sum_(k=0)^(m-1) binom(m-1,k)(-4)^k Q_(p-1,k).

Every coefficient is nonzero modulo ``p``.  Hence omitting
``Q_(p-1,m-1)``, as the implementation does, selects ``m-1`` independent
rows complementary to the ordinary line.

These equations are independent and exhaust the extra characteristic-p
left kernel.  Here is a direct block-rank proof which also controls finite-
field polynomial aliasing.  A symmetric function on an unordered pair, with
one common value on the diagonal, has the unique midpoint/difference basis

    1;  a^i delta^(2b),   0<=i<=p-1, 1<=b<=m,       (3)

where ``a=(s+t)/2`` and ``delta=(s-t)/2``.  Indeed ``delta^2`` takes
``m+1`` values and the terms with ``b>=1`` are
exactly the functions vanishing at ``delta=0``.  Pulling (3) back along a
projective functional ``L`` gives the mutually independent bidegree blocks

    L(a)^i L(delta)^(2b).

There is no hidden function alias between these blocks.  Scaling ``a`` and
``delta`` separately distinguishes their degrees modulo ``p-1``; the value
at ``a=0`` separates degrees zero and ``p-1``.  There is one further
same-character pair which must be checked: the separately displayed constant
``1`` and the block ``(i,b)=(0,m)``, whose ``delta`` degree is ``p-1``.
If a homogeneous form ``F`` of degree ``p-1`` were constant ``c`` on every
nonzero ``delta``, then ``F(1,T)-c`` would have degree at most ``p-1`` and
all ``p`` field elements as roots.  Thus ``F(X,Y)=cX^(p-1)``; evaluation at
``(0,1)`` forces ``c=0``.  This proves that the two blocks are disjoint.
Evaluation of each remaining binary form of degree at most ``p-1`` on
``F_p^2`` (or on its nonzero points for ``delta``) is injective.

For fixed ``i,b``, put ``n=i+2b``.  The ``p+1`` such vectors, as ``L`` runs
over ``(1,lambda)`` and ``(0,1)``, have rank

    min(p+1,n+1).                                      (4)

To see this without assuming characteristic-zero Veronese independence,
expand the finite vectors as ``sum_(k=0)^n lambda^k w_k``.  All binomial
coefficients occurring in the two factors are nonzero because
``i,2b<=p-1``, and the ``w_k`` have disjoint monomial supports.  If
``n<p``, the finite Vandermonde matrix has rank ``n+1``.  If ``n>=p``, its
only aliases in the present range ``n<=2p-2`` are
``lambda^(k+p-1)=lambda^k`` for ``1<=k<=p-1``; it has rank ``p``.  In
particular the finite columns at ``n-(p-1)`` and ``n`` coincide, whereas the
infinity vector is exactly ``w_n``.  It therefore separates that aliased
pair and gives rank ``p+1``.  Thus (4) holds for every block, not merely in
small-prime replays.

Consequently the full image rank is

    1+sum_(b=1)^m sum_(i=0)^(p-1) min(p+1,i+2b+1). (5)

Subtracting (5) from the target gives the ``p+1`` ordinary boundary/total
equations plus

    S(p)=sum_(d=2)^(p-2) floor(d/2)*(p-d) +(m-1)
        =(m-1)(4m^2+7m+6)/6,                p=2m+1.       (6)

The displayed moment equations have exactly this cardinality.  Scaling
the fibre labels separates different degrees modulo ``p-1`` and the
``Q_(d,k)`` are independent within one degree, so the equations are
independent.  The only intersection with the ordinary equations is the
degree-``p-1`` vector ``(s-t)^(p-1)``; it was explicitly removed above.
The rank count (5) then proves exhaustion.

Equivalently, for the target of dimension
``(p+1)*(binom(p,2)+1)``, the characteristic-p rank is

    target-(p+1)-S(p).                                    (7)

The residual normalized rows differ from unsigned edge-Radon rows only by
one sign on every edge column and one sign on every direction block.  These
are unimodular changes.  Consequently (1) remains necessary after inserting
the direction sign in its outer sum.  The same direction signs must be
inserted in every finite term and in the infinity term of (2); the basis (3)
and all ranks are unchanged.

This is a genuine common-edge obstruction beyond Proposition 15.757's
complete binary image.  It does not by itself exclude the compact local
survivor rays of Proposition 15.758: their atom labels have not been forced
through all equations (1).  Residual (ii) therefore remains open.
"""
from __future__ import annotations

import json
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


def p_torsion_codimension(p: int) -> dict[str, object]:
    """Return the full-projective summands and both exact rank formulas."""
    _check_prime(p)
    m = (p - 1) // 2
    by_degree = {
        str(d): (d // 2) * (p - d)
        for d in range(2, p - 1)
    }
    if m > 1:
        by_degree[str(p - 1)] = m - 1
    direct = sum(by_degree.values())
    closed = (m - 1) * (4 * m * m + 7 * m + 6) // 6
    target = (p + 1) * (comb(p, 2) + 1)
    rank = target - (p + 1) - direct
    block_rank = 1 + sum(
        min(p + 1, i + 2 * b + 1)
        for b in range(1, m + 1)
        for i in range(p)
    )
    proved = direct == closed and rank == block_rank
    if not proved:
        raise ArithmeticError("the p-primary codimension identity changed")
    return {
        "p": p,
        "m": m,
        "projective_direction_count": p + 1,
        "target_dimension": target,
        "equal_total_equation_count": p,
        "parallel_sum_equation_count": 1,
        "ordinary_boundary_total_codimension": p + 1,
        "relation_count_by_degree": by_degree,
        "extra_p_primary_codimension_direct_sum": direct,
        "extra_p_primary_codimension_closed_form": closed,
        "rank_in_characteristic_p": rank,
        "independent_bidegree_block_rank": block_rank,
        "block_rank_formula": (
            "1+sum_(b=1)^m sum_(i=0)^(p-1) min(p+1,i+2b+1)"
        ),
        "proved": proved,
    }


def top_degree_basis_complement(p: int) -> dict[str, object]:
    """Record the all-prime change of basis at degree ``p-1``.

    The ordinary off-diagonal-total row is ``D^m``, not one literal member
    of the displayed ``Q`` basis.  Its coefficient on every basis member is
    nonzero, so deleting the last member gives a canonical complement.
    """
    _check_prime(p)
    m = (p - 1) // 2
    coefficients = tuple(
        (comb(m - 1, k) * pow(-4, k, p)) % p
        for k in range(m)
    )
    retained = tuple(range(max(0, m - 1)))
    omitted = m - 1
    proved = bool(
        len(coefficients) == m
        and all(coefficient != 0 for coefficient in coefficients)
        and len(retained) == m - 1
        and omitted not in retained
        and coefficients[omitted] != 0
    )
    if not proved:
        raise ArithmeticError("the top-degree ordinary complement changed")
    return {
        "p": p,
        "m": m,
        "ordinary_vector": "D^m=(s-t)^(p-1)",
        "basis": "Q_(p-1,k), 0<=k<=m-1",
        "expansion": "D^m=sum_k binom(m-1,k)(-4)^k Q_(p-1,k)",
        "expansion_coefficients_mod_p": list(coefficients),
        "retained_new_relation_indices": list(retained),
        "omitted_basis_index": omitted,
        "omitted_coefficient_nonzero": True,
        "constant_top_alias_resolution": (
            "a degree-(p-1) homogeneous form constant on all nonzero delta "
            "is zero: use p finite-chart roots, then the infinity point"
        ),
        "proved": proved,
    }


def moment_relation_rows(p: int) -> tuple[tuple[int, int, int], ...]:
    """Return ``(d,k,j)`` for every new projective relation (1)--(2)."""
    _check_prime(p)
    rows = tuple(
        (d, k, j)
        for d in range(2, p - 1)
        for k in range(d // 2)
        for j in range(p - d)
    ) + tuple((p - 1, k, 0) for k in range(max(0, (p - 3) // 2)))
    expected = int(p_torsion_codimension(p)["extra_p_primary_codimension_closed_form"])
    if len(rows) != expected:
        raise ArithmeticError("the moment-relation indexing changed")
    return rows


def moment_polynomial(p: int, s: int, t: int, d: int, k: int) -> int:
    """Evaluate ``Q_(d,k)`` in ``F_p``."""
    _check_prime(p)
    if not 2 <= d <= p - 1 or not 0 <= k < d // 2:
        raise ValueError("need 2<=d<=p-1 and 0<=k<floor(d/2)")
    return (
        pow((s - t) % p, 2, p)
        * pow((s * t) % p, k, p)
        * pow((s + t) % p, d - 2 - 2 * k, p)
    ) % p


def _full_edge_radon_counts(
    p: int, selected_edges: tuple[tuple[int, int], ...]
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    """Push an explicit edge subset to all ``p+1`` directions."""
    points = tuple(product(range(p), repeat=2))
    pairs = tuple(combinations(range(p), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    offdiag = [[0] * len(pairs) for _ in range(p + 1)]
    parallel = [0] * (p + 1)
    for u, v in selected_edges:
        x1, y1 = points[u]
        x2, y2 = points[v]
        for lam in range(p):
            s = (x1 + lam * y1) % p
            t = (x2 + lam * y2) % p
            if s == t:
                parallel[lam] += 1
            else:
                offdiag[lam][pair_index[tuple(sorted((s, t)))]] += 1
        s = y1
        t = y2
        if s == t:
            parallel[p] += 1
        else:
            offdiag[p][pair_index[tuple(sorted((s, t)))]] += 1
    return tuple(tuple(row) for row in offdiag), tuple(parallel)


def deterministic_edge_sample(p: int) -> tuple[tuple[int, int], ...]:
    """A nontrivial deterministic graph used only to replay all relations."""
    _check_prime(p)
    points = tuple(product(range(p), repeat=2))
    return tuple(
        (u, v)
        for u, v in combinations(range(len(points)), 2)
        if (
            3 * u
            + 5 * v
            + points[u][0] * points[v][1]
            + 2 * points[u][1] * points[v][0]
        )
        % 7
        in (0, 1, 3)
    )


def verify_moment_relations(p: int) -> dict[str, object]:
    """Check the complete hierarchy on one deterministic common graph."""
    _check_prime(p)
    pairs = tuple(combinations(range(p), 2))
    edges = deterministic_edge_sample(p)
    offdiag, _parallel = _full_edge_radon_counts(p, edges)
    residues = []
    for d, k, j in moment_relation_rows(p):
        residue = 0
        for lam in range(p):
            contraction = sum(
                offdiag[lam][index] * moment_polynomial(p, s, t, d, k)
                for index, (s, t) in enumerate(pairs)
            )
            residue += pow(lam, j, p) * contraction
        if j == p - d - 1:
            residue += sum(
                offdiag[p][index] * moment_polynomial(p, s, t, d, k)
                for index, (s, t) in enumerate(pairs)
            )
        residues.append(residue % p)
    proved = all(value == 0 for value in residues)
    if not proved:
        raise ArithmeticError("a p-primary moment relation failed")
    return {
        "p": p,
        "selected_edge_count": len(edges),
        "relations_checked": len(residues),
        "all_residues_zero_mod_p": proved,
        "proved": proved,
    }


def _rank_mod_prime(matrix: list[list[int]], modulus: int) -> int:
    """Exact Gaussian rank over the prime field ``F_modulus``."""
    if not matrix:
        return 0
    rows = [[value % modulus for value in row] for row in matrix]
    row_count = len(rows)
    column_count = len(rows[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, row_count) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, modulus)
        rows[rank] = [(inverse * value) % modulus for value in rows[rank]]
        pivot_row = rows[rank]
        for index in range(row_count):
            if index == rank or rows[index][column] == 0:
                continue
            factor = rows[index][column]
            rows[index] = [
                (value - factor * pivot_value) % modulus
                for value, pivot_value in zip(rows[index], pivot_row)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def exact_full_edge_radon_rank(p: int, modulus: int) -> dict[str, object]:
    """Build the full ``p+1``-direction map and replay its small-prime rank."""
    _check_prime(p)
    if (
        not isinstance(modulus, int)
        or isinstance(modulus, bool)
        or modulus < 2
        or not is_prime(modulus)
    ):
        raise ValueError("modulus must be prime")
    points = tuple(product(range(p), repeat=2))
    edges = tuple(combinations(range(len(points)), 2))
    pairs = tuple(combinations(range(p), 2))
    rows: list[list[int]] = []
    directions = tuple((1, lam) for lam in range(p)) + ((0, 1),)
    for left, right in directions:
        labels = tuple((left * x + right * y) % p for x, y in points)
        for s, t in pairs:
            rows.append(
                [
                    int(tuple(sorted((labels[u], labels[v]))) == (s, t))
                    for u, v in edges
                ]
            )
        rows.append([int(labels[u] == labels[v]) for u, v in edges])
    actual = _rank_mod_prime(rows, modulus)
    ledger = p_torsion_codimension(p)
    target = int(ledger["target_dimension"])
    expected = (
        int(ledger["rank_in_characteristic_p"])
        if modulus == p
        else target - (p + 1)
    )
    proved = actual == expected
    if not proved:
        raise ArithmeticError("the exact full edge-Radon rank changed")
    return {
        "p": p,
        "modulus": modulus,
        "matrix_rows": len(rows),
        "matrix_columns": len(edges),
        "exact_rank": actual,
        "formula_rank": expected,
        "proved": proved,
    }


def signed_residual_transport() -> dict[str, object]:
    """Record the unimodular sign transport to residual normalization."""
    return {
        "raw_edge_variable": "n_e in Z",
        "column_sign": "tau_e, the Paley sign of edge e",
        "direction_row_sign": "epsilon_L, the Paley sign of the parallel class",
        "normalized_parallel_entry": "epsilon_L*tau_e*n_e=n_e",
        "normalized_transverse_entry": "epsilon_L*tau_e*n_e",
        "finite_moment_outer_coefficient": "epsilon_L*lambda^j",
        "projective_infinity_coefficient": "epsilon_infinity",
        "projective_infinity_term_is_signed": True,
        "moment_relation_repair": (
            "multiply every finite outer coefficient and the infinity term by epsilon_L"
        ),
        "row_and_column_sign_changes_are_unimodular": True,
        "p_primary_codimension_unchanged": True,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    ledgers = {str(p): p_torsion_codimension(p) for p in (3, 5, 7, 11, 13)}
    samples = {str(p): verify_moment_relations(p) for p in (3, 5, 7, 11)}
    ranks = {
        f"p{p}_mod{modulus}": exact_full_edge_radon_rank(p, modulus)
        for p in (3, 5, 7)
        for modulus in (2, p)
    }
    signed = signed_residual_transport()
    top_degree = top_degree_basis_complement(13)
    proved = bool(
        all(row["proved"] for row in ledgers.values())
        and all(row["proved"] for row in samples.values())
        and all(row["proved"] for row in ranks.values())
        and signed["proved"]
        and top_degree["proved"]
    )
    return {
        "prop": "15.759",
        "title": "Complete low-degree p-torsion edge-Radon hierarchy",
        "status": "PROVED THEOREM AND PROVED OPEN REDUCTION",
        "hierarchy": (
            "for 2<=d<=p-2: floor(d/2) contractions Q divisible by "
            "(s-t)^2, each with p-d full-projective evaluation checks; "
            "degree p-1 adds m-1 checks"
        ),
        "codimension_formula": "(m-1)(4m^2+7m+6)/6, p=2m+1",
        "proved": {
            "all_endpoint_moment_relations": proved,
            "relations_independent_and_exhaust_p_primary_linear_cokernel": proved,
            "closed_codimension_formula": proved,
            "signed_residual_normalization_transport": proved,
            "compact_aggregate_survivor_excluded": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "codimension_ledgers": ledgers,
        "sample_relation_checks": samples,
        "exact_rank_replays": ranks,
        "signed_transport": signed,
        "all_prime_proof_audit": {
            "equal_total_equation_count": "p among p+1 directional totals",
            "parallel_sum_equation_count": 1,
            "constant_top_scaling_alias": (
                "disjoint by the p finite-chart roots and the infinity point"
            ),
            "top_degree_ordinary_complement": (
                "omit k=m-1 after expanding D^m with all coefficients nonzero"
            ),
            "signed_projective_completion": (
                "epsilon_L multiplies every finite term and epsilon_infinity "
                "multiplies the infinity term"
            ),
            "symbolic_top_degree_replay": top_degree,
            "proved": True,
        },
        "remaining_obstruction": (
            "Test the compact p+1/p-3 atom decompositions against the whole "
            "endpoint-moment hierarchy, then impose nonnegative simple-edge lifting."
        ),
        "duplicate_work_guards": [
            "Do not stop at M2 or M4 when the full endpoint hierarchy is available.",
            "Do not mistake F2 image completeness for saturation of the integer lattice.",
            "Do not claim the p-torsion hierarchy excludes the compact aggregate family without an atom-level proof.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15759.json"
    atomic_write_json(destination, out)
    print("Prop. 15.759 p-primary edge-Radon hierarchy: proved")
    print("  compact aggregate exclusion: OPEN")
    print(f"  wrote {destination}")
    return out


if __name__ == "__main__":
    main()
