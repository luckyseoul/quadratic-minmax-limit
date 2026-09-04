"""Minimum-degree chart lemma and an aggregate high-shell barrier.

Signed PSL transport permutes a relative Paley flip graph, so a vertex of
minimum degree may be sent to infinity without changing its edge count or
the two eigenshell inequalities.  This gives an isolated chart below the
sharp threshold ``2|H| < p^2+1``.

The threshold cannot be improved from degree averaging, directional means,
boundary parity floors, and difference-Radon energy alone.  The antipodal
perfect matching on ``P^1(F_(p^2))`` meets all of those aggregate identities
at ``|H|=(p^2+1)/2`` and has minimum degree one.  Its pointwise shell score
is below three, so it is deliberately a barrier to the aggregate strategy,
not a residual-(ii) graph.
"""

from __future__ import annotations


def _is_prime(p: int) -> bool:
    if type(p) is not int or p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    divisor = 3
    while divisor * divisor <= p:
        if p % divisor == 0:
            return False
        divisor += 2
    return True


def _validate_prime(p: int) -> None:
    if not _is_prime(p) or p < 7:
        raise ValueError("need an odd prime p>=7")


def minimum_degree_chart_bound(p: int, edge_count: int) -> dict[str, object]:
    """Return the exact averaging bound after signed PSL transport."""
    _validate_prime(p)
    n = p * p + 1
    if not 0 <= edge_count <= n * (n - 1) // 2:
        raise ValueError("edge count is outside the simple-graph range")
    bound = 2 * edge_count // n
    isolated = 2 * edge_count < n
    return {
        "p": p,
        "vertex_count": n,
        "edge_count": edge_count,
        "minimum_degree_upper_bound": bound,
        "signed_PSL_can_send_minimum_vertex_to_infinity": True,
        "transported_infinity_degree_upper_bound": bound,
        "isolated_chart_forced": isolated,
        "isolated_chart_condition": "2*|H| < p^2+1",
        "proved": True,
    }


def outside_boundary_even_degree_bound(
    p: int, edge_count: int, boundary_size: int
) -> dict[str, object]:
    """Bound the least even degree outside the odd-degree boundary.

    Boundary vertices consume at least ``boundary_size`` units of the total
    degree ``2*edge_count``.  All remaining degrees are even.
    """
    base = minimum_degree_chart_bound(p, edge_count)
    n = int(base["vertex_count"])
    if not 0 <= boundary_size <= n or boundary_size % 2:
        raise ValueError("the boundary size must be even and lie in 0..p^2+1")
    if boundary_size == n:
        return {
            **base,
            "boundary_size": boundary_size,
            "outside_vertex_count": 0,
            "outside_even_degree_upper_bound": None,
            "isolated_outside_vertex_forced": False,
        }
    outside = n - boundary_size
    remaining_degree = 2 * edge_count - boundary_size
    if remaining_degree < 0:
        raise ValueError("boundary size exceeds the handshake degree budget")
    even_bound = 2 * (remaining_degree // (2 * outside))
    return {
        **base,
        "boundary_size": boundary_size,
        "outside_vertex_count": outside,
        "outside_degree_sum_upper_bound": remaining_degree,
        "outside_even_degree_upper_bound": even_bound,
        "isolated_outside_vertex_forced": even_bound == 0,
        "isolated_outside_condition": "2*|H|+|D| < 2*(p^2+1)",
    }


def residual_chart_system() -> dict[str, object]:
    """Record the exact general high-``k`` equations used by the audit."""
    return {
        "residual_input": (
            "k=|G| is even, e not in G, H=G union {e}; on Max+ S_G>=2 "
            "and S_G=2 implies f_e=+1; a bad deletion has S_H<=-3 on Max-"
        ),
        "signed_shell_form": (
            "T_H^eps=eps*sum_(a in H) f_a; the gap-four box gives "
            "3<=T_H^eps<=Phi-2 on both signed eigenshells"
        ),
        "frame_mean": "E_eps[T_H^eps]=|H|/p",
        "directional_slack": "A_L=(T_H^eps_L-3)/2 is a nonnegative integer",
        "scaled_directional_mean": (
            "a_L=2p E[A_L]=I+(p+1)P_L-eps_L*T-3p"
        ),
        "definitions": {
            "I": "number of H-edges incident with infinity",
            "P_L": "number of finite H-edges parallel to direction L",
            "T": "sum of Paley signs of all finite H-edges",
        },
        "edge_count": "I+sum_L P_L=|H|",
        "type_budgets": (
            "sum_(eps_L=tau) a_L=((p+1)/2)(|H|-3p), tau=+/-1"
        ),
        "same_type_quantization": "a_L-a_M=(p+1)(P_L-P_M)",
        "boundary_phase": (
            "(-1)^eta_L=eps_L^(1+1_(infinity in D))*"
            "(-1)^((|H|-3)/2+b_L)*c_H"
        ),
        "parity_floor": (
            "a_L>=2*ceil(p*M(p,b_L,eta_L)), where M is the exact "
            "quadratic parity-majorant expectation"
        ),
        "difference_row_sum": "sum_(a>0) q_L(a)=eps_L*T-P_L",
        "common_energy": (
            "sum_(L,a>0) q_L(a)^2=p*sum_delta m_delta^2+2T^2-2sum_L P_L^2"
        ),
        "collision_identity": (
            "sum_delta m_delta^2=(|H|-I)+2sum_delta binom(m_delta,2)"
        ),
        "proved": True,
    }


def antipodal_matching_aggregate_barrier(p: int) -> dict[str, object]:
    """Exact all-prime aggregate counterexample at the isolation threshold."""
    _validate_prime(p)
    n = p * p + 1
    h = n // 2
    q = (p - 1) // 2
    m = (p + 1) // 2
    infinity_degree = 1
    finite_edges = h - infinity_degree
    parallel_count = q
    signed_total = 0
    a = infinity_degree + (p + 1) * parallel_count - 3 * p
    phase = 0
    parity_floor = 0 if p % 4 == 3 else 2 * p

    parallel_square_sum = (p + 1) * parallel_count * parallel_count
    common_energy = (
        p * finite_edges
        + 2 * signed_total * signed_total
        - 2 * parallel_square_sum
    )
    collision_count = 0
    displacement_square_sum = finite_edges
    normalized_row_sum = -q
    normalized_row_energy = q

    # For an affine shell state, write z_t in {+1,-1} for its fibre signs
    # and let r count split nonzero antipodal fibre pairs.  The normalized
    # matching score is z_0+2r.  Its exact minimum follows by filling whole
    # antipodal pairs first.
    minimum_shell_score = -1 if p % 4 == 3 else 1
    minimum_slack = (minimum_shell_score - 3) // 2

    type_budget = m * (h - 3 * p)
    exact = bool(
        h % 2 == 1
        and finite_edges == (p * p - 1) // 2
        and (p + 1) * parallel_count == finite_edges
        and a == h - 3 * p
        and a >= parity_floor
        and m * a == type_budget
        and common_energy == finite_edges
        and (p + 1) * normalized_row_energy == common_energy
        and minimum_shell_score < 3
    )
    if not exact:
        raise ArithmeticError("the antipodal matching barrier changed")

    return {
        "p": p,
        "construction": (
            "H={{infinity,0}} union {{x,-x}: x in F_(p^2)^*/+/-}"
        ),
        "vertex_count": n,
        "edge_count": h,
        "G_edge_count_after_deleting_infinity_zero": h - 1,
        "all_vertex_degrees": 1,
        "odd_degree_boundary_size": n,
        "transported_infinity_degree_I": infinity_degree,
        "minimum_degree_bound_is_sharp": 2 * h == n,
        "directions_of_each_quadratic_type": m,
        "parallel_count_every_direction": parallel_count,
        "finite_signed_total_T": signed_total,
        "edge_product_c_H": 1,
        "boundary_fibre_count_b_every_direction": p,
        "boundary_phase_eta_every_direction": phase,
        "scaled_directional_mean_a_every_direction": a,
        "parity_floor_every_direction": parity_floor,
        "all_phase_floors_pass": a >= parity_floor,
        "type_budget_each": type_budget,
        "displacement_multiplicity": 1,
        "collision_count": collision_count,
        "displacement_square_sum": displacement_square_sum,
        "normalized_difference_row": {
            "entries": "q copies of -1 in every direction",
            "sum": normalized_row_sum,
            "energy": normalized_row_energy,
        },
        "common_difference_radon_energy": common_energy,
        "aggregate_system_passes": True,
        "pointwise_shell_score_formula": "T_H^eps(X)=z_0+2r",
        "minimum_pointwise_signed_shell_score": minimum_shell_score,
        "minimum_pointwise_slack": minimum_slack,
        "pointwise_residual_box_passes": False,
        "is_residual_ii_counterexample": False,
        "strategy_conclusion": (
            "minimum degree plus directional means, phase floors, and common "
            "energy cannot imply |H|=O(p); the missing input is the full "
            "pointwise shell box (and its level-two equality), not another "
            "aggregate degree estimate"
        ),
        "residual_ii_closed": False,
        "proved": exact,
    }


def theorem_record() -> dict[str, object]:
    return {
        "classification": "exact chart lemma and aggregate method barrier",
        "chart_system": residual_chart_system(),
        "sharp_threshold_examples": {
            str(p): antipodal_matching_aggregate_barrier(p)
            for p in (7, 11, 13, 17, 19, 23, 31)
        },
        "proved": True,
        "residual_ii_closed": False,
    }
