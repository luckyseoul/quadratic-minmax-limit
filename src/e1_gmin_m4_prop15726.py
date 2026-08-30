#!/usr/bin/env python3
r"""Prop. 15.726 -- a linear tangent-envelope exclusion for outside slack.

Let ``D`` be the ``p+1`` affine points in an outside chart of the first
general boundary, and let ``R`` be the normalized pair slack from
Proposition 15.722.  Delete an inclusion-minimal set ``T`` so that
``A=D\T`` is an arc, and write ``t=|T|``.  The prior deletion lemma gives
``1<=t<=R``.

For ``z in T``, let ``s_z`` be the number of secants of ``A`` through
``z``.  Minimality gives ``s_z>=1``.  On a line containing two points of
``A`` and ``u`` points of ``T``, the slack contribution satisfies

    h(2+u) >= u.

Indeed, for ``u=2a`` and ``u=2a+1`` the gaps are respectively
``a(a-1)`` and ``a^2``.  Consequently

    I := sum_(z in T) s_z <= R,

so some ``z`` has ``s_z<=floor(R/t)``.

The arc has size ``p+1-t=p+2-tau`` with ``tau=t+1``.  In odd order,
Segre's tangent envelope in Ball--Lavrauw gives a nonzero homogeneous dual
polynomial ``Phi`` of degree ``2*tau`` satisfying

    Phi(X cross P) = f_P(X)^2                         (P in A),

where ``f_P`` is the product of the ``tau`` tangent forms at ``P``.  This
is Theorem 13 in the authors' current manuscript and Theorem 11 in
arXiv v4.  Its size hypothesis is ``|A|>=2*tau+2``.

If ``1<=R<=floor((p-4)/3)``, then for every ``1<=t<=R`` the proof has the
stronger strict inequality ``|A|>2*tau+2``.  Moreover

    3*t + 2*R/t <= 3*R+2 <= p-2 < p-1,

because

    3*R+2-(3*t+2*R/t)=(R-t)(3-2/t)>=0.

Thus the chosen point lies on more than ``2*tau`` tangents to ``A``.
Their dual points are more than ``deg(Phi)`` zeros on the dual line
``z*``, forcing ``z*`` to divide ``Phi``.  But minimality supplies an
``A``-secant ``zPQ``.  At its dual point the tangent-envelope identity
gives ``Phi((zP)*)=f_P(z)^2!=0``, since ``zP`` is a secant rather than a
tangent at ``P``.  This contradicts ``z* | Phi``.

Therefore every positive outside slack

    1 <= R <= floor((p-4)/3)

is impossible for every prime ``p>=17``.  In particular the former seven
``R=4`` cells are closed, with cutoffs ``4,5,6,8,9,11,12`` at
``p=17,19,23,29,31,37,41``.  This is a geometric theorem, not a finite
profile census.  It does not close the remaining ``p+1`` shell,
non-Walsh residual (ii), Type I, or the quadratic-minmax limit.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15722 import occupancy_slack_term


ROOT = Path(__file__).resolve().parents[1]
SMALL_PRIMES = (17, 19, 23, 29, 31, 37, 41)


def _check_prime_parameter(p: int) -> None:
    if p < 17 or not is_prime(p):
        raise ValueError("need an odd prime parameter p>=17")


def linear_tangent_envelope_cutoff(p: int) -> int:
    """Largest integer ``R`` covered by ``3R<=p-4``."""
    _check_prime_parameter(p)
    return (p - 4) // 3


def linewise_slack_incidence_lemma(u: int) -> dict[str, object]:
    """Prove ``h(2+u)>=u`` for an arbitrary nonnegative integer ``u``.

    Here the line contains two retained arc points and ``u`` deleted
    points.  Hence ``u`` is exactly its contribution to
    ``sum_z s_z``, while ``h(2+u)`` is its contribution to ``R``.
    """
    if u < 0:
        raise ValueError("the number of deleted points must be nonnegative")
    slack = occupancy_slack_term(2 + u)
    if u % 2 == 0:
        a = u // 2
        symbolic_gap = a * (a - 1)
        gap_formula = "u=2a: h(2+u)-u=a(a-1)>=0"
    else:
        a = (u - 1) // 2
        symbolic_gap = a * a
        gap_formula = "u=2a+1: h(2+u)-u=a^2>=0"
    exact_gap = slack - u
    return {
        "u": u,
        "line_occupancy": 2 + u,
        "slack_contribution": slack,
        "secant_incidence_contribution": u,
        "gap": exact_gap,
        "gap_formula": gap_formula,
        "formula_matches": exact_gap == symbolic_gap,
        "proved": exact_gap == symbolic_gap and exact_gap >= 0,
    }


def tangent_envelope_dependency() -> dict[str, object]:
    """Record the exact external theorem and the version-number split."""
    return {
        "external_dependency": True,
        "authors": "Simeon Ball and Michel Lavrauw",
        "title": "Planar arcs",
        "journal": "Journal of Combinatorial Theory, Series A 160 (2018), 261-287",
        "doi": "10.1016/j.jcta.2018.06.015",
        "arxiv": "1705.10940v4",
        "current_manuscript_url": (
            "https://web.mat.upc.edu/simeon.michael.ball/planararcs.pdf"
        ),
        "current_manuscript_theorem": 13,
        "arxiv_v4_theorem": 11,
        "odd_order_statement": (
            "if an arc A has size q+2-tau and |A|>=2*tau+2, a nonzero "
            "degree-2*tau dual polynomial Phi satisfies "
            "Phi(X cross P)=f_P(X)^2 for every P in A"
        ),
        "size_hypothesis_is_weak": True,
        "proof_meets_stronger_strict_size_bound": True,
        "tangent_roots_on_each_point_pencil_are_doubled": True,
        "proved": True,
    }


def universal_symbolic_ledger() -> dict[str, object]:
    """Return the parameter-free proof chain behind Proposition 15.726."""
    return {
        "result_status": "proved theorem",
        "minimal_arc_repair": (
            "choose inclusion-minimal T with A=D\\T an arc; 1<=t=|T|<=R"
        ),
        "minimality": "every z in T lies on at least one A-secant",
        "linewise_incidence_bound": {
            "even_deleted_count": "u=2a gives h(2+u)-u=a(a-1)>=0",
            "odd_deleted_count": "u=2a+1 gives h(2+u)-u=a^2>=0",
            "conclusion": "I=sum_(z in T) s_z<=R",
        },
        "average_point": "some z in T has 1<=s_z<=floor(R/t)",
        "arc_parameters": "|A|=p+1-t=p+2-tau, tau=t+1",
        "cutoff": "1<=R<=floor((p-4)/3), equivalently 3R<=p-4",
        "stronger_strict_size_bound": (
            "p>=3t+4 gives |A|=p+1-t>=2t+5>2t+4=2tau+2"
        ),
        "algebraic_identity": (
            "3R+2-(3t+2R/t)=(R-t)(3-2/t)>=0"
        ),
        "strict_tangent_chain": (
            "3t+2s_z<=3t+2R/t<=3R+2<=p-2<p-1, "
            "so |A|-2s_z>2tau"
        ),
        "dual_root_count": (
            "more than deg(Phi)=2tau tangent zeros on z* force z*|Phi"
        ),
        "secant_evaluation": (
            "for an A-secant zPQ, Phi((zP)*)=f_P(z)^2!=0"
        ),
        "contradiction": "the secant dual point lies on z* but Phi is nonzero there",
        "external_theorem": tangent_envelope_dependency(),
        "finite_profile_search_used": False,
        "proved": True,
    }


def deletion_size_incidence_contradiction(
    p: int, R: int, t: int
) -> dict[str, object]:
    """Audit the proof for one possible minimal deletion size ``t``."""
    cutoff = linear_tangent_envelope_cutoff(p)
    if not 1 <= R <= cutoff:
        raise ValueError(f"need 1<=R<={cutoff} at p={p}")
    if not 1 <= t <= R:
        raise ValueError("need 1<=t<=R")

    arc_size = p + 1 - t
    tau = t + 1
    envelope_degree = 2 * tau
    theorem_size_threshold = 2 * tau + 2

    # I<=R and |T|=t give a point with this integer upper bound.
    secant_index_upper = R // t
    tangent_count_lower = arc_size - 2 * secant_index_upper

    # This is the denominator-cleared form of
    # 3*t+2*R/t <= 3*R+2.
    average_inequality_gap_times_t = (R - t) * (3 * t - 2)
    cutoff_gap = p - 4 - 3 * R
    tangent_degree_gap = tangent_count_lower - envelope_degree

    stronger_strict_size_bound = arc_size > theorem_size_threshold
    average_inequality = average_inequality_gap_times_t >= 0
    cutoff_inequality = cutoff_gap >= 0
    more_tangents_than_degree = tangent_degree_gap > 0
    proved = bool(
        stronger_strict_size_bound
        and average_inequality
        and cutoff_inequality
        and secant_index_upper >= 1
        and more_tangents_than_degree
    )
    return {
        "p": p,
        "R": R,
        "minimal_deletion_size_t": t,
        "arc_size": arc_size,
        "tangent_deficiency_tau": tau,
        "envelope_degree": envelope_degree,
        "theorem_size_threshold": theorem_size_threshold,
        "stronger_strict_size_bound_met": stronger_strict_size_bound,
        "incidence_sum_upper": R,
        "chosen_point_secant_index_lower": 1,
        "chosen_point_secant_index_upper": secant_index_upper,
        "tangent_count_lower": tangent_count_lower,
        "tangent_count_minus_envelope_degree": tangent_degree_gap,
        "average_inequality_gap_times_t": average_inequality_gap_times_t,
        "cutoff_gap_p_minus_4_minus_3R": cutoff_gap,
        "dual_line_component_forced": more_tangents_than_degree,
        "minimality_supplies_secant": True,
        "secant_evaluation_nonzero": True,
        "contradiction": proved,
        "proved": proved,
    }


def all_deletion_sizes_excluded(p: int, R: int) -> dict[str, object]:
    """Check every possible ``1<=t<=R`` under the symbolic cutoff."""
    cutoff = linear_tangent_envelope_cutoff(p)
    if not 1 <= R <= cutoff:
        raise ValueError(f"need 1<=R<={cutoff} at p={p}")
    rows = [deletion_size_incidence_contradiction(p, R, t) for t in range(1, R + 1)]
    return {
        "p": p,
        "R": R,
        "possible_minimal_deletion_sizes": list(range(1, R + 1)),
        "rows": rows,
        "every_deletion_size_contradictory": all(row["proved"] for row in rows),
        "proved": all(row["proved"] for row in rows),
    }


def small_prime_cutoff_table() -> dict[int, int]:
    """Exact cutoffs at the seven primes where ``R=4`` was formerly live."""
    table = {p: linear_tangent_envelope_cutoff(p) for p in SMALL_PRIMES}
    expected = {17: 4, 19: 5, 23: 6, 29: 8, 31: 9, 37: 11, 41: 12}
    if table != expected:
        raise ArithmeticError("small-prime linear cutoff table changed")
    return table


def outside_linear_slack_exclusion(p: int) -> dict[str, object]:
    """Prove the complete positive interval through ``floor((p-4)/3)``."""
    cutoff = linear_tangent_envelope_cutoff(p)
    rows = [all_deletion_sizes_excluded(p, R) for R in range(1, cutoff + 1)]
    case_rows = [case for row in rows for case in row["rows"]]
    proved = all(row["proved"] for row in rows)
    return {
        "p": p,
        "excluded_positive_R_interval": [1, cutoff],
        "excluded_positive_R_values": list(range(1, cutoff + 1)),
        "first_possible_positive_R_at_least": cutoff + 1,
        "audited_integer_case_count": len(case_rows),
        "all_integer_cases_verified": proved,
        "minimum_strict_size_margin": min(
            row["arc_size"] - row["theorem_size_threshold"]
            for row in case_rows
        ),
        "minimum_tangent_degree_gap": min(
            row["tangent_count_minus_envelope_degree"] for row in case_rows
        ),
        "symbolic_ledger": universal_symbolic_ledger(),
        "former_R_four_cell_closed": cutoff >= 4,
        "finite_profile_search_used": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "result_status": "proved theorem",
        "proved": proved,
    }


def proposition_15726() -> dict[str, object]:
    """Package the general theorem and the exact former small-prime cells."""
    table = small_prime_cutoff_table()
    prime_rows = {str(p): outside_linear_slack_exclusion(p) for p in SMALL_PRIMES}
    proved = all(row["proved"] for row in prime_rows.values())
    return {
        "prop": "15.726",
        "statement": (
            "for every prime p>=17, outside pair slack "
            "1<=R<=floor((p-4)/3) is impossible at |D|=p+1"
        ),
        "small_prime_cutoffs": {str(p): cutoff for p, cutoff in table.items()},
        "small_prime_rows": prime_rows,
        "universal_symbolic_ledger": universal_symbolic_ledger(),
        "R_four_closed_for_all_primes_p_ge_17": True,
        "remaining_scope": (
            "larger outside slack, the rest of non-Walsh residual (ii), "
            "multi-level Type I, and L remain open"
        ),
        "top_level_gates_changed": False,
        "result_status": "proved theorem",
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the proposition ledger when the coordinated update is ready."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15726.json"
    output.write_text(json.dumps(proposition_15726(), indent=2, sort_keys=True) + "\n")
    return output


def main() -> None:
    result = proposition_15726()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.726 tangent-envelope audit failed")
    path = write_evidence()
    print("Prop 15.726 linear tangent-envelope slack exclusion: proved")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
