#!/usr/bin/env python3
r"""Prop. 15.750 -- isolated-chart parity halving closes Type I.

Let ``G`` be a multi-level Type-I bad case for a prime ``p>=5`` and a
distinguished edge ``e``.  In the feature notation of Proposition 15.275,

    |G|=3p-2,  e not in G,
    S_G(y)=3-2f_e(y)                         on Max+,
    S_G(y)<=-1 and S_G(y)<=-3f_e(y)          on Max-.

Put ``H=G union {e}`` and let ``W=H+e=G+2e`` be a multiset.  Then
``|H|=3p-1``, ``|W|=3p``, ``S_W=3`` on Max+, and ``S_H<=-2`` on Max-.
For ``p>=11``, ``H`` has an isolated vertex because

    p^2+1-2(3p-1) = p^2-6p+3 > 0.

Proposition 15.721 transports that vertex to infinity.  The relative-flip
identity preserves the eigenshell equations and inequalities; applying the
same edge permutation to ``W`` also preserves its multiplicities.  Hence
all edges under discussion are finite in the transported affine chart.

Write ``q=(p-1)/2`` and ``m=(p+1)/2``.  In a square direction, let ``P`` be
the total ``W``-multiplicity parallel to the direction and ``K_st`` the
signed off-fibre block sums.  The affine Max+ cylinders give

    P + sum_(s<t) K_st z_s z_t = 3            (sum z_s=1).

The Johnson-slice swap lemma says that a pair quadratic constant on this
middle slice has all pair coefficients equal, say to ``kappa``.  Since
``sum_(s<t) z_s z_t=-q``, and because each block contains at least the
absolute value of its signed sum,

    P=3+q*kappa,
    3p-P >= binom(p,2)|kappa| = p*q*|kappa|.

For ``p>=11`` these relations force ``kappa=0`` and ``P=3``.  Summing over
the ``m`` square directions shows that ``W`` has positive and negative
Paley multiplicities ``3m`` and ``3q`` respectively.

Now fix a nonsquare direction, let ``P_d`` count its parallel ``H``-edges,
and define the nonnegative integer-valued quadratic

    T_d(z)=(-S_H(z)-2)/2.

If ``c=C_e`` and ``tau=sum C_ab H_ab=3-c``, middle-slice averaging gives

    a_d := 2p E[T_d] = (p+1)P_d-2p+tau.        (*)

Nonnegativity forces every ``P_d>=2``.  Their sum is ``3q`` when ``c=1``
and ``3q-1`` when ``c=-1``.  Both are below ``3m``, so one direction has
``P_d=2``.  Equation (*) then gives ``a_d=4`` or ``6``.

The product of the edge features of ``H`` is an affine parity on the fibre
slice.  Since

    (-1)^((|H|-S_H)/2) = product_(g in H) f_g,

``T_d mod 2`` is that affine parity, up to a constant.  For a nonconstant
parity on ``J(p,m)``, the central Krawtchouk estimate

    |K_m(b)| <= binom(p,m)/p,  1<=b<=p-1,

puts at least ``(p-1)/(2p)`` of the slice in either parity class.  Thus an
odd nonnegative ``T_d`` would cost ``2p E[T_d]>=p-1`` (and a constantly odd
one would cost ``2p``).  Since ``a_d<=6<p-1``, ``T_d`` is everywhere even.
Therefore ``B_d=T_d/2`` is a nonzero nonnegative integer-valued quadratic.
Proposition 15.688 gives

    4p E[B_d] >= p-3,

whereas ``4p E[B_d]=a_d in {4,6}<p-3``.  This contradiction closes every
prime ``p>=11``.  The tracked integer Farkas identities in
``e1_type_i_small_prime_exact`` close ``p=5,7`` for the canonical ordered
edge ``(infinity,0)`` without SciPy, an optimizer, or an eigenshell cache.
This loses no generality: ``PSL(2,p^2)`` is 2-transitive.  After sending the
first endpoint to infinity, a determinant-one translation fixing infinity
sends the second endpoint to zero; the signed lift preserves the bad-box
system.  Consequently the multi-level Type-I bad case is empty for every
prime ``p>=5`` and every distinguished edge.

This proposition does not prove residual (ii), E1, or the quadratic
min--max limit.  In particular, it does not revive the obsolete ``3A+B``
or ``Aut_e`` routes in Proposition 15.275.
"""
from __future__ import annotations

import json
from functools import lru_cache
from math import comb
from pathlib import Path

from e1_type_i_small_prime_exact import (
    type_I_badcase_small_primes_exact,
    verify_type_i_badcase_farkas,
)
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def is_prime(p: int) -> bool:
    """Return whether ``p`` is prime, for theorem-domain validation."""
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


def isolated_transport_certificate(p: int) -> dict[str, object]:
    """Replay the isolated-vertex count and the signed transport dependency."""
    if not is_prime(p) or p < 7:
        raise ValueError("need a prime p>=7")
    from e1_gmin_m4_prop15721 import (
        mobius_boundary_normalization,
        signed_relative_flip_transport,
    )

    edge_count_h = 3 * p - 1
    vertex_count = p * p + 1
    isolated_lower_bound = vertex_count - 2 * edge_count_h
    transport = signed_relative_flip_transport()
    normalization = mobius_boundary_normalization(p)
    proved = bool(
        isolated_lower_bound > 0
        and transport["proved"]
        and transport["flip_set_size_preserved"]
        and transport["Max_plus_and_Max_minus_shells_preserved"]
        and normalization["proved"]
    )
    _require(proved, "isolated signed transport dependency changed")
    return {
        "p": p,
        "vertices": vertex_count,
        "H_edges": edge_count_h,
        "isolated_vertex_lower_bound": isolated_lower_bound,
        "relative_flip_transport": {
            "proved": bool(transport["proved"]),
            "flip_set_size_preserved": bool(
                transport["flip_set_size_preserved"]
            ),
            "shells_preserved": bool(
                transport["Max_plus_and_Max_minus_shells_preserved"]
            ),
        },
        "selected_isolated_vertex_maps_to_infinity": bool(
            normalization["selected_boundary_point_maps_to_infinity"]
        ),
        "multiset_W_transport": (
            "The signed-relative identity extends linearly from Boolean "
            "masks to all nonnegative integer edge multiplicities; apply "
            "the same transported edge permutation to H+e."
        ),
        "relative_flip_transport_extends_to_integer_edge_multisets": True,
        "proved": proved,
    }


def distinguished_edge_normalization_certificate(p: int) -> dict[str, object]:
    """Normalize any ordered edge to ``(infinity,0)`` inside signed PSL.

    Translate a finite first endpoint to zero and apply ``z -> -1/z``;
    this determinant-one map sends it to infinity.  The other endpoint is
    then finite.  A translation fixing infinity sends that image to zero.
    If an endpoint is already infinity, omit or reorder the first two steps.
    Both matrices have determinant one, proving 2-transitivity directly.
    """
    if not is_prime(p) or p < 5:
        raise ValueError("need a prime p>=5")
    from e1_gmin_m4_prop15721 import signed_relative_flip_transport

    transport = signed_relative_flip_transport()
    proved = bool(
        transport["proved"]
        and transport["Max_plus_and_Max_minus_shells_preserved"]
    )
    _require(proved, "distinguished-edge normalization changed")
    return {
        "p": p,
        "field": f"F_({p}^2)",
        "first_map": "z -> -1/(z-a), with determinant 1",
        "second_map": "z -> z-u, fixing infinity, with determinant 1",
        "ordered_edge_target": ["infinity", 0],
        "PSL_2_is_two_transitive": True,
        "canonical_matrix_edge_indices": [0, 1],
        "signed_lift_preserves_bad_box_system": bool(transport["proved"]),
        "proved": proved,
    }


def johnson_pair_constancy_certificate(p: int) -> dict[str, object]:
    r"""State the uniform Johnson-slice swap lemma used in square directions.

    If ``sum K_st z_s z_t`` is constant on ``sum z=1``, swap coordinates
    ``a,b`` while fixing a middle set among the other ``p-2`` coordinates.
    The difference is a linear form with coefficients ``K_ak-K_bk``.  Its
    vanishing on that slice makes all those coefficients equal; their dot
    product with a sign vector of sum one makes the common value zero.
    Varying ``a,b`` makes every pair coefficient equal.
    """
    if not is_prime(p) or p < 7:
        raise ValueError("need a prime p>=7")
    remaining_coordinates = p - 2
    remaining_positive_signs = (p - 1) // 2
    _require(
        1 < remaining_positive_signs < remaining_coordinates,
        "swap slice is degenerate",
    )
    return {
        "p": p,
        "slice": "z in {+-1}^p with sum(z)=1",
        "swap_remaining_coordinates": remaining_coordinates,
        "swap_remaining_positive_signs": remaining_positive_signs,
        "linear_swap_forms_vanish": True,
        "all_pair_coefficients_equal": True,
        "proof": (
            "swap a,b; then swap two remaining coordinates; finally use "
            "the remaining sign sum 1 to kill the common difference"
        ),
        "proof_kind": "uniform exact Johnson-slice lemma",
        "proved": True,
    }


def square_direction_kappa_candidates(
    p: int,
    *,
    remove_doubled_edge_endpoint: bool = True,
) -> list[int]:
    """Return all integral ``kappa`` surviving the exact mass inequalities."""
    if not is_prime(p) or p < 7:
        raise ValueError("need a prime p>=7")
    q = (p - 1) // 2
    # These are the exact integral bounds imposed by 0 <= P=3+q*kappa <= 3p.
    kappa_min = -(3 // q)
    kappa_max = (3 * p - 3) // q
    survivors: list[int] = []
    for kappa in range(kappa_min, kappa_max + 1):
        parallel = 3 + q * kappa
        if 3 * p - parallel < p * q * abs(kappa):
            continue
        if (
            remove_doubled_edge_endpoint
            and p == 7
            and kappa == -1
            and parallel == 0
        ):
            # Equality in the l1 estimate makes every off-fibre block have
            # multiplicity one.  But P=0 puts the doubled edge e in one of
            # those blocks, whose multiplicity is at least two.
            continue
        survivors.append(kappa)
    return survivors


def square_direction_rigidity_certificate(p: int) -> dict[str, object]:
    """Package the Johnson and l1 arguments for one square direction."""
    if not is_prime(p) or p < 7:
        raise ValueError("need a prime p>=7")
    q = (p - 1) // 2
    raw = square_direction_kappa_candidates(
        p, remove_doubled_edge_endpoint=False
    )
    survivors = square_direction_kappa_candidates(p)
    johnson = johnson_pair_constancy_certificate(p)
    proved = bool(johnson["proved"] and survivors == [0])
    _require(proved, "square-direction rigidity changed")
    return {
        "p": p,
        "q": q,
        "johnson_pair_constancy": johnson,
        "identities": [
            "P=3+q*kappa",
            "3p-P>=p*q*abs(kappa)",
        ],
        "raw_kappa_candidates": raw,
        "doubled_edge_endpoint_removed": raw != survivors,
        "kappa_candidates": survivors,
        "parallel_W_multiplicity": 3,
        "all_off_fibre_signed_block_sums": 0,
        "proved": proved,
    }


def square_signed_mass_certificate(p: int) -> dict[str, object]:
    """Sum the square-direction rigidity identities over finite slopes."""
    if not is_prime(p) or p < 7:
        raise ValueError("need a prime p>=7")
    q = (p - 1) // 2
    m = q + 1
    rigidity = square_direction_rigidity_certificate(p)
    positive_w = 3 * m
    negative_w = 3 * p - positive_w
    _require(negative_w == 3 * q, "signed W mass identity changed")
    return {
        "p": p,
        "square_directions": m,
        "nonsquare_directions": m,
        "parallel_W_units_per_square_direction": 3,
        "positive_W_multiplicity": positive_w,
        "negative_W_multiplicity": negative_w,
        "signed_W_total": positive_w - negative_w,
        "every_finite_edge_has_one_projective_direction": True,
        "square_direction_rigidity": rigidity,
        "proved": bool(rigidity["proved"] and positive_w - negative_w == 3),
    }


def nonsquare_parallel_two_certificate(p: int, c_e: int) -> dict[str, object]:
    """Force a nonsquare direction with ``P_d=2`` and mean 4 or 6."""
    if not is_prime(p) or p < 11 or c_e not in (-1, 1):
        raise ValueError("need a prime p>=11 and c_e in {-1,+1}")
    q = (p - 1) // 2
    m = q + 1
    signed_mass = square_signed_mass_certificate(p)
    tau_h = 3 - c_e
    negative_h = 3 * q - (1 if c_e == -1 else 0)
    mean_if_parallel_one = (p + 1) - 2 * p + tau_h
    minimum_parallel = 2
    all_at_least_three_mass = 3 * m
    shortfall = all_at_least_three_mass - negative_h
    scaled_mean = 2 * (p + 1) - 2 * p + tau_h
    proved = bool(
        signed_mass["proved"]
        and mean_if_parallel_one < 0
        and shortfall > 0
        and scaled_mean == (4 if c_e == 1 else 6)
    )
    _require(proved, "nonsquare parallel-two argument changed")
    return {
        "p": p,
        "C_e": c_e,
        "tau_H": tau_h,
        "nonsquare_directions": m,
        "sum_parallel_H_over_nonsquare_directions": negative_h,
        "scaled_mean_formula": "2p E[T_d]=(p+1)P_d-2p+tau_H",
        "scaled_mean_at_P_le_1_upper": mean_if_parallel_one,
        "nonnegativity_forces_P_at_least": minimum_parallel,
        "mass_if_all_P_at_least_3": all_at_least_three_mass,
        "mass_shortfall_from_all_P_at_least_3": shortfall,
        "some_direction_has_P": 2,
        "a_equals_2p_E_T": scaled_mean,
        "proved": proved,
    }


def central_krawtchouk(p: int, b: int) -> int:
    """Return the exact central value ``K_m(b)`` for ``p=2m-1``."""
    if type(p) is not int or p < 3 or p % 2 == 0:
        raise ValueError("need odd p>=3")
    if type(b) is not int or not 0 <= b <= p:
        raise ValueError("need 0<=b<=p")
    m = (p + 1) // 2
    lower = max(0, m - (p - b))
    upper = min(b, m)
    return sum(
        (-1) ** j * comb(b, j) * comb(p - b, m - j)
        for j in range(lower, upper + 1)
    )


def parity_bias_exact_replay(p: int) -> dict[str, object]:
    """Replay the Krawtchouk recurrence and bias bound for one odd ``p``."""
    if not is_prime(p) or p < 5:
        raise ValueError("need an odd prime p>=5")
    m = (p + 1) // 2
    denominator = comb(p, m)
    values = tuple(central_krawtchouk(p, b) for b in range(p + 1))
    base = values[1] == values[2] == -denominator // p
    recurrence = all(
        (p - b) * values[b + 1]
        == -values[b] - b * values[b - 1]
        for b in range(1, p)
    )
    bound = all(
        p * abs(values[b]) <= denominator for b in range(1, p)
    )
    proved = bool(base and recurrence and bound)
    _require(proved, "central Krawtchouk replay changed")
    return {
        "p": p,
        "m": m,
        "denominator": denominator,
        "K_1": values[1],
        "K_2": values[2],
        "recurrence": "(p-b)K_(b+1)=-K_b-bK_(b-1)",
        "recurrence_exact": recurrence,
        "max_p_times_absolute_K": max(
            p * abs(values[b]) for b in range(1, p)
        ),
        "bias_bound_denominator": denominator,
        "each_nonconstant_parity_value_probability_at_least": (
            f"{p - 1}/{2 * p}"
        ),
        "proved": proved,
    }


def parity_bias_theorem() -> dict[str, object]:
    """Record the uniform induction proving the central parity bias bound."""
    return {
        "domain": "odd p>=5, m=(p+1)/2, 1<=b<=p-1",
        "base": "K_1=K_2=-binom(p,m)/p",
        "recurrence": "(p-b)K_(b+1)=-K_b-bK_(b-1)",
        "induction_range": "2<=b<=(p-1)/2",
        "induction_ratio": "(b+1)/(p-b)<=1",
        "other_half": "absolute-value complement symmetry",
        "conclusion": "abs(K_m(b))<=binom(p,m)/p",
        "probability_consequence": (
            "each value of a nonconstant affine parity on J(p,m) has "
            "probability at least (p-1)/(2p)"
        ),
        "proof_kind": "uniform exact recurrence induction",
        "proved": True,
    }


def general_prime_row(p: int, c_e: int) -> dict[str, object]:
    """Replay every exact scalar inequality in the ``p>=11`` proof."""
    if not is_prime(p) or p < 11 or c_e not in (-1, 1):
        raise ValueError("need a prime p>=11 and c_e in {-1,+1}")
    from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor

    transport = isolated_transport_certificate(p)
    nonsquare = nonsquare_parallel_two_certificate(p, c_e)
    parity = parity_bias_exact_replay(p)
    floor = sharp_integral_quadratic_lift_floor(p)
    scaled_mean = int(nonsquare["a_equals_2p_E_T"])
    parity_gap = p - 1 - scaled_mean
    lift_gap = int(floor["sharp_scaled_floor"]) - scaled_mean
    proved = bool(
        transport["proved"]
        and nonsquare["proved"]
        and parity["proved"]
        and floor["proved"]
        and parity_gap > 0
        and lift_gap > 0
    )
    _require(proved, "general Type-I contradiction changed")
    return {
        "p": p,
        "C_e": c_e,
        "isolated_transport": transport,
        "nonsquare_parallel_two": nonsquare,
        "parity_bias_replay": parity,
        "a_equals_2p_E_T": scaled_mean,
        "nonconstant_or_odd_parity_floor": p - 1,
        "parity_gap": parity_gap,
        "T_is_constantly_even": parity_gap > 0,
        "B_equals_T_over_2_is_nonzero": scaled_mean > 0,
        "sharp_lift_floor_15_688": int(floor["sharp_scaled_floor"]),
        "lift_floor_gap": lift_gap,
        "contradiction": (
            "4p E[B]=2p E[T]=a is below the Proposition 15.688 floor"
        ),
        "closed": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def uniform_general_prime_certificate() -> dict[str, object]:
    """Prove the inequalities uniformly for every prime ``p>=11``.

    Sample rows are included only as deterministic replays.  Universality
    comes from the displayed monotone symbolic gaps, not from enumeration.
    """
    from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor

    transport = isolated_transport_certificate(11)
    johnson = johnson_pair_constancy_certificate(11)
    parity = parity_bias_theorem()
    lift = sharp_integral_quadratic_lift_floor(11)

    # Each expression is minimized at p=11 over odd p>=11.  The positive
    # kappa gap is the difference between the l1 lower bound and available
    # off-fibre mass at kappa=1; larger kappa only strengthens it.
    minimum_gaps = {
        "isolated_vertices_p2_minus_6p_plus_3": 11 * 11 - 6 * 11 + 3,
        "positive_kappa_l1_gap_(p-1)(p-5)/2": (11 - 1) * (11 - 5) // 2,
        "negative_kappa_parallel_gap_q_minus_3": (11 - 1) // 2 - 3,
        "P_le_1_negative_mean_gap_p_minus_5": 11 - 5,
        "nonsquare_all_three_mass_shortfall": 3,
        "parity_gap_p_minus_1_minus_6": 11 - 1 - 6,
        "lift_gap_p_minus_3_minus_6": 11 - 3 - 6,
    }
    symbolic_monotonicity = {
        "p2_minus_6p_plus_3_increment": "2p-5>0",
        "positive_kappa_gap_increment": "p-2>0",
        "negative_kappa": "P=3-q*abs(kappa)<0 because q>=5",
        "P_le_1_mean": "at most 5-p<0",
        "nonsquare_shortfall": "3 for C_e=1 and 4 for C_e=-1",
        "parity_gap": "p-7>=4",
        "lift_gap": "p-9>=2",
    }
    samples = {
        str(p): {
            str(c): general_prime_row(p, c)
            for c in (-1, 1)
        }
        for p in (11, 13, 17)
    }
    proved = bool(
        transport["proved"]
        and johnson["proved"]
        and parity["proved"]
        and lift["proved"]
        and all(value > 0 for value in minimum_gaps.values())
        and all(
            row["proved"]
            for prime_rows in samples.values()
            for row in prime_rows.values()
        )
    )
    _require(proved, "uniform p>=11 Type-I proof changed")
    return {
        "domain": "every prime p>=11",
        "proof_kind": "uniform theorem; samples are replays only",
        "bad_case_setup": {
            "G_edges": "3p-2",
            "H": "G union {e}",
            "H_edges": "3p-1",
            "W": "G+2e as an edge multiset",
            "W_edges_with_multiplicity": "3p",
            "Max_plus_W_sum": 3,
            "Max_minus_H_upper": -2,
        },
        "signed_transport_dependency": "Proposition 15.721",
        "quadratic_lift_dependency": "Proposition 15.688",
        "johnson_pair_constancy": johnson,
        "parity_bias": parity,
        "minimum_symbolic_gaps_at_p_11": minimum_gaps,
        "symbolic_monotonicity": symbolic_monotonicity,
        "sample_exact_replays": samples,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def type_I_multilevel_bad_case_closed_all_primes() -> bool:
    """Authoritative all-prime Type-I closure predicate."""
    small = type_I_badcase_small_primes_exact()
    all_edges = all(
        distinguished_edge_normalization_certificate(p)["proved"]
        for p in (5, 7)
    )
    general = bool(uniform_general_prime_certificate()["proved"])
    _require(small and all_edges and general, "all-prime Type-I closure failed")
    return True


def type_i_multilevel_bad_case_closed_all_primes() -> bool:
    """PEP-8 alias for the authoritative all-prime predicate."""
    return type_I_multilevel_bad_case_closed_all_primes()


@lru_cache(maxsize=1)
def proposition_15750() -> dict[str, object]:
    """Return the full exact Proposition 15.750 theorem record."""
    small = {
        str(p): {
            "canonical_edge_exact_Farkas": verify_type_i_badcase_farkas(p),
            "all_distinguished_edges": distinguished_edge_normalization_certificate(p),
            "proved": True,
        }
        for p in (5, 7)
    }
    general = uniform_general_prime_certificate()
    closed = type_I_multilevel_bad_case_closed_all_primes()
    proved = bool(
        closed
        and general["proved"]
        and all(
            row["proved"]
            and row["canonical_edge_exact_Farkas"]["proved"]
            and row["all_distinguished_edges"]["proved"]
            for row in small.values()
        )
    )
    _require(proved, "Proposition 15.750 failed")
    return {
        "prop": "15.750",
        "title": "Isolated-chart parity halving closes multi-level Type I",
        "result_status": "proved all-prime theorem with exact base certificates",
        "statement": (
            "the multi-level Type-I bad-case box is empty for every prime p>=5"
        ),
        "small_prime_exact_Farkas_certificates": small,
        "uniform_p_at_least_11": general,
        "type_I_multilevel_bad_case_ND_closed": closed,
        "obsolete_3A_plus_B_route_used": False,
        "finite_graph_census_used": False,
        "scipy_or_optimizer_theorem_dependency": False,
        "eigenshell_cache_theorem_dependency": False,
        "residual_ii_closed": False,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "L_status": "OPEN",
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    """Write deterministic Proposition 15.750 evidence atomically."""
    target = path or ROOT / "evidence" / "e1_gmin_m4_prop15750.json"
    write_json_atomic(target, proposition_15750())
    return target


def main() -> None:
    result = proposition_15750()
    path = write_evidence()
    print(
        json.dumps(
            {
                "prop": result["prop"],
                "result_status": result["result_status"],
                "type_I_multilevel_bad_case_ND_closed": result[
                    "type_I_multilevel_bad_case_ND_closed"
                ],
                "residual_ii_closed": result["residual_ii_closed"],
                "output": str(path),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
