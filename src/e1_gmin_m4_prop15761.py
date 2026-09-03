#!/usr/bin/env python3
r"""Prop. 15.761 -- exact full edge-Radon spectrum and least-norm barrier.

Let ``R`` be the unsigned integral edge-Radon map of Proposition 15.760.
In each of the ``d=p+1`` directions it has one parallel row ``P_L`` and
``C=binom(p,2)=pm`` off-diagonal fibre-pair rows ``K_L(s,t)``, where
``p=2m+1``.  This proposition computes the complete real spectrum of
``RR^t`` on the ordinary image lattice and applies the Moore--Penrose lower
bound to the compact rays of Proposition 15.758.

The row intersections are elementary and exact.  Rows in one direction are
disjoint.  For distinct directions ``L,M``, affine coordinates ``(L,M)``
give

    |P_L cap P_M|=0,
    |P_L cap K_M(s,t)|=p,
    |K_L(a,b) cap K_M(s,t)|=2.                       (1)

The last intersection consists of the two diagonal matchings of the four
coordinate corners.  A parallel row has squared norm ``p^2 m`` and an
off-diagonal row squared norm ``p^2``.

Decompose a real target in the ordinary image ``A_R`` as follows.  Its
common directional total is ``T``.  In direction ``L`` write

    K_L=k_L 1+w_L,       sum w_L=0,
    k_L=(T-P_L)/C.

Then split ``P_L=a+q_L``, where ``a=T/d`` and ``sum q_L=0``.  The three
summands are mutually orthogonal.  Equation (1) gives the eigenvalues

    p^2                 on the w_L spaces,
    p(C+1)              on (P_L,K_L)=(q_L,-q_L/C),
    p^2(m+p)            on (P_L,K_L)=(a,a/m).         (2)

Their multiplicities are respectively

    d(C-1),             p,                            1,

whose sum is ``dC=rank A_R``.  Hence the exact least squared norm of any
real edge vector mapping to the target is

    ||R^+y||^2
      =p^-2 sum_L ||w_L||^2
       +(pC)^-1 sum_L q_L^2
       +2T^2/[p^2(p^2-1)].                           (3)

For the residual normalization, let ``eta_L`` be the direction-row sign,
let ``P_L>=0`` be the actual parallel count, and let ``W_L`` be the
normalized off-diagonal coefficient row.  Undoing the row signs replaces
``(P_L,K_L)`` by ``(eta_L P_L,eta_L W_L)``.  Source and target sign matrices
are orthogonal, so (3) becomes

    Q(W,P,T)=p^-2 sum_L [||W_L||^2-(eta_L T-P_L)^2/C]
       +(pC)^-1 sum_L(eta_L P_L-T/(p+1))^2
       +2T^2/[p^2(p^2-1)].                           (4)

If ``H`` is a simple graph with normalized source coordinates
``z_e=tau_e 1_(e in H)``, then ``||z||^2=|H|``.  Orthogonal projection onto
``(ker R)^perp`` therefore gives the necessary full-midpoint inequality

    Q(W,P,T) <= |H|.                                  (5)

This is strictly finer data than Proposition 15.758's scalar
difference-Radon Parseval expression: (4) retains every individual
fibre-pair coefficient before the midpoint is summed out.

Nevertheless (5) cannot exclude either compact ray.  The conclusion is
uniform in all atom labels.  Each compact or all-equal triangle has
coefficient l1 norm three, a star has norm ``p-1``, and the omitted-pair
row has norm ``2p-3``.  For the balanced allocations used in 15.758, the
whole allowed ``t`` intervals give

    p=4r+1:  e_L<=2r-3, Q_L<=2r-1,
              ||W_hard||_1<=10r-9,
              ||W_opp ||_1<=14r-10;                  (6)

    p=4r+3:  e_L<=2r-2, Q_L<=2r+2,
              ||W_hard||_1<=10r-4,
              ||W_opp ||_1<= 6r-3.                   (7)

Overlaps only decrease these triangle-inequality upper bounds.  Dropping
the negative terms in (4), using ``||W||_2<=||W||_1``, and bounding each of
the last two terms in (4) by one gives

    Q < (2r+1)[(10r-9)^2+(14r-10)^2]/(4r+1)^2+2,
                                                               (8)
    Q < (2r+2)[(10r-4)^2+(6r-3)^2]/(4r+3)^2+2.                 (9)

The minimum edge counts on the two rays are

    h_B=4r^2+6r+5,       h_C=4r^2+8r+9.

After multiplying the gaps ``h_B-(8)`` and ``h_C-(9)`` by their positive
square denominators and putting ``u=r-7``, their numerators are

    64u^4+1328u^3+9796u^2+29864u+30706,
    64u^4+1744u^3+18108u^2+85374u+154867.            (10)

Every coefficient is positive.  Thus (5) holds with strict room for every
``r>=7`` and every target made from the recorded compact atom templates,
even after their labels are adjusted to meet moment equations.

This is a proved method barrier, not a graph construction.  A least-norm
real preimage need not be integral, sign-compatible, nonnegative, or
``0/1``.  Proposition 15.760 reduces the remaining question to the affine
integer fibre intersecting the signed Boolean box; (3)--(10) show that its
Euclidean projection cannot separate the two infinite compact rays.
Residual (ii) remains open.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb


Polynomial = tuple[Fraction, ...]


def _trim(poly: Polynomial) -> Polynomial:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def _add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    return _trim(tuple(
        (left[index] if index < len(left) else Fraction(0))
        + (right[index] if index < len(right) else Fraction(0))
        for index in range(size)
    ))


def _sub(left: Polynomial, right: Polynomial) -> Polynomial:
    return _add(left, tuple(-value for value in right))


def _mul(left: Polynomial, right: Polynomial) -> Polynomial:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return _trim(tuple(out))


def _scale(poly: Polynomial, scalar: int | Fraction) -> Polynomial:
    return _trim(tuple(Fraction(scalar) * value for value in poly))


def _pow(poly: Polynomial, exponent: int) -> Polynomial:
    out: Polynomial = (Fraction(1),)
    for _ in range(exponent):
        out = _mul(out, poly)
    return out


def _shift(poly: Polynomial, amount: int) -> Polynomial:
    """Return ``poly(u+amount)`` exactly."""
    out: Polynomial = (Fraction(0),)
    for degree, coefficient in enumerate(poly):
        term = [Fraction(0)] * (degree + 1)
        for power in range(degree + 1):
            term[power] = coefficient * comb(degree, power) * amount ** (degree - power)
        out = _add(out, tuple(term))
    return _trim(out)


def _serialise(poly: Polynomial) -> list[int | str]:
    return [
        int(value) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
        for value in poly
    ]


def full_radon_spectrum_certificate() -> dict[str, object]:
    """Record the complete real spectrum (1)--(3)."""
    one: Polynomial = (Fraction(1),)
    m_poly: Polynomial = (Fraction(0), Fraction(1))
    p_poly = _add(one, _scale(m_poly, 2))
    d_poly = _add(p_poly, one)
    c_poly = _mul(p_poly, m_poly)
    lambda_directional = _mul(p_poly, _add(c_poly, one))
    lambda_uniform = _mul(_pow(p_poly, 2), _add(m_poly, p_poly))
    checks = {
        "multiplicities_sum_to_ordinary_rank": (
            _add(
                _mul(d_poly, _sub(c_poly, one)),
                _add(p_poly, one),
            )
            == _mul(d_poly, c_poly)
        ),
        "directional_parallel_coordinate_eigenvalue": (
            _add(_mul(_pow(p_poly, 2), m_poly), p_poly)
            == lambda_directional
        ),
        "directional_pair_coordinate_eigenvalue_after_clearing_C": (
            _add(_pow(p_poly, 2), _mul(c_poly, _sub(p_poly, _scale(one, 2))))
            == lambda_directional
        ),
        "uniform_parallel_coordinate_eigenvalue": (
            _add(
                _mul(_pow(p_poly, 2), m_poly),
                _pow(p_poly, 3),
            )
            == lambda_uniform
        ),
        "uniform_pair_coordinate_eigenvalue_after_clearing_m": (
            _mul(_pow(p_poly, 2), _add(one, _scale(m_poly, 3)))
            == lambda_uniform
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError("the symbolic full edge-Radon spectrum changed")
    return {
        "one_direction_row_norms_squared": {
            "parallel": "p^2*m",
            "off_diagonal_pair": "p^2",
        },
        "distinct_direction_row_intersections": {
            "parallel_parallel": 0,
            "parallel_pair": "p",
            "pair_pair": 2,
        },
        "orthogonal_blocks": {
            "within_direction_pair_zero_sum": {
                "eigenvalue": "p^2",
                "multiplicity": "(p+1)(C-1)",
            },
            "directional_aggregate_zero_sum": {
                "coordinates": "(P_L,K_L)=(q_L,-q_L/C), sum q_L=0",
                "eigenvalue": "p(C+1)",
                "multiplicity": "p",
            },
            "uniform": {
                "coordinates": "(P_L,K_L)=(a,a/m)",
                "eigenvalue": "p^2(m+p)",
                "multiplicity": 1,
            },
        },
        "multiplicity_sum": "(p+1)(C-1)+p+1=(p+1)C=rank(A_R)",
        "least_norm_formula": (
            "p^-2 sum_L||w_L||^2+(pC)^-1 sum_L q_L^2+"
            "2T^2/[p^2(p^2-1)]"
        ),
        "symbolic_eigenvalue_checks": checks,
        "complete_full_midpoint_spectrum": True,
        "proved": True,
    }


def signed_least_norm_certificate() -> dict[str, object]:
    """State formula (4) and its signed-Boolean consequence (5)."""
    return {
        "C": "binom(p,2)=p(p-1)/2",
        "signed_formula": (
            "Q=p^-2 sum_L[||W_L||^2-(eta_L*T-P_L)^2/C]+"
            "(pC)^-1 sum_L(eta_L*P_L-T/(p+1))^2+"
            "2T^2/[p^2(p^2-1)]"
        ),
        "source_signing": "z_e=tau_e*1_(e in H)",
        "simple_source_squared_norm": "||z||^2=|H|",
        "necessary_inequality": "Q<=|H|",
        "reason": "orthogonal projection onto (ker R)^perp",
        "strictly_finer_than_difference_aggregate_parseval": True,
        "proved": True,
    }


def compact_ray_norm_bounds() -> dict[str, object]:
    """Record the interval-uniform atom counts and l1 bounds (6)--(7)."""
    return {
        "p_1_mod_4": {
            "parameter": "p=4r+1, r>=7",
            "balanced_count_bounds": {"e_L": "<=2r-3", "Q_L": "<=2r-1"},
            "hard_l1": "<=p-1+3e_L<=10r-9",
            "opposite_l1": "<=2p-3+3(Q_L-2)<=14r-10",
            "least_norm_upper": (
                "(2r+1)[(10r-9)^2+(14r-10)^2]/(4r+1)^2+2"
            ),
            "minimum_H_edge_count": "4r^2+6r+5",
        },
        "p_3_mod_4": {
            "parameter": "p=4r+3, r>=7",
            "balanced_count_bounds": {"e_L": "<=2r-2", "Q_L": "<=2r+2"},
            "hard_l1": "<=p-1+3e_L<=10r-4",
            "opposite_l1": "<=3(Q_L-3)<=6r-3",
            "least_norm_upper": (
                "(2r+2)[(10r-4)^2+(6r-3)^2]/(4r+3)^2+2"
            ),
            "minimum_H_edge_count": "4r^2+8r+9",
        },
        "uniform_in_atom_labels": True,
        "negative_terms_in_exact_formula_discarded": True,
        "each_remaining_aggregate_and_uniform_term_is_less_than_one": True,
        "proved": True,
    }


def symbolic_gap_certificate() -> dict[str, object]:
    """Prove all inequalities after the formal substitution ``r=u+7``.

    This manipulates polynomial coefficients only; it instantiates no prime
    and performs no finite parameter scan.
    """
    one: Polynomial = (Fraction(1),)
    r: Polynomial = (Fraction(0), Fraction(1))

    p_b = _add(_scale(r, 4), one)
    d_b = _add(_scale(r, 4), _scale(one, 2))
    m_b = _add(_scale(r, 2), one)
    hard_b = _sub(_scale(r, 10), _scale(one, 9))
    opposite_b = _sub(_scale(r, 14), _scale(one, 10))
    h_b = _add(_add(_scale(_pow(r, 2), 4), _scale(r, 6)), _scale(one, 5))
    main_gap_b = _sub(
        _mul(_sub(h_b, _scale(one, 2)), _pow(p_b, 2)),
        _mul(m_b, _add(_pow(hard_b, 2), _pow(opposite_b, 2))),
    )
    aggregate_gap_b = _sub(
        _mul(_pow(p_b, 2), _sub(p_b, one)),
        _scale(_mul(d_b, _pow(_add(_scale(r, 2), _scale(one, 4)), 2)), 2),
    )
    uniform_gap_b = _sub(
        _mul(_pow(p_b, 2), _sub(_pow(p_b, 2), one)),
        _scale(_pow(_add(_scale(r, 4), _scale(one, 5)), 2), 2),
    )

    p_c = _add(_scale(r, 4), _scale(one, 3))
    d_c = _add(_scale(r, 4), _scale(one, 4))
    m_c = _add(_scale(r, 2), _scale(one, 2))
    hard_c = _sub(_scale(r, 10), _scale(one, 4))
    opposite_c = _sub(_scale(r, 6), _scale(one, 3))
    h_c = _add(_add(_scale(_pow(r, 2), 4), _scale(r, 8)), _scale(one, 9))
    main_gap_c = _sub(
        _mul(_sub(h_c, _scale(one, 2)), _pow(p_c, 2)),
        _mul(m_c, _add(_pow(hard_c, 2), _pow(opposite_c, 2))),
    )
    aggregate_gap_c = _sub(
        _mul(_pow(p_c, 2), _sub(p_c, one)),
        _scale(_mul(d_c, _pow(_add(_scale(r, 2), _scale(one, 3)), 2)), 2),
    )
    uniform_gap_c = _sub(
        _mul(_pow(p_c, 2), _sub(_pow(p_c, 2), one)),
        _scale(_pow(_sub(_scale(r, 4), one), 2), 2),
    )

    shifted = {
        "p_1_main_gap": _shift(main_gap_b, 7),
        "p_1_aggregate_term_denominator_gap": _shift(aggregate_gap_b, 7),
        "p_1_uniform_term_denominator_gap": _shift(uniform_gap_b, 7),
        "p_3_main_gap": _shift(main_gap_c, 7),
        "p_3_aggregate_term_denominator_gap": _shift(aggregate_gap_c, 7),
        "p_3_uniform_term_denominator_gap": _shift(uniform_gap_c, 7),
    }
    proved = all(all(coefficient > 0 for coefficient in poly) for poly in shifted.values())
    if not proved:
        raise ArithmeticError("the all-r least-norm gap lost coefficientwise positivity")
    return {
        "substitution": "u=r-7>=0",
        "coefficient_order": "constant_to_highest_degree_in_u",
        "shifted_numerators": {
            name: _serialise(poly) for name, poly in shifted.items()
        },
        "every_coefficient_strictly_positive": True,
        "no_prime_or_parameter_scan": True,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    spectrum = full_radon_spectrum_certificate()
    signed = signed_least_norm_certificate()
    bounds = compact_ray_norm_bounds()
    gaps = symbolic_gap_certificate()
    proved = bool(
        spectrum["proved"] and signed["proved"] and bounds["proved"] and gaps["proved"]
    )
    return {
        "prop": "15.761",
        "title": "Exact full edge-Radon spectrum and least-norm method barrier",
        "status": "PROVED STRICTLY STRONGER FULL-MIDPOINT METHOD BARRIER",
        "full_spectrum": spectrum,
        "signed_least_norm": signed,
        "compact_ray_bounds": bounds,
        "symbolic_gap": gaps,
        "proved": {
            "complete_real_spectrum_on_ordinary_image": proved,
            "full_midpoint_least_norm_necessary_inequality": proved,
            "strictly_stronger_than_prior_difference_parseval": proved,
            "both_prop15758_rays_pass_with_strict_room_for_r_ge_7": proved,
            "one_common_simple_graph_constructed": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "remaining_obstruction": (
            "The signed Boolean box, equivalently squarefree nonnegative source "
            "coordinates inside the integral affine fibre of Proposition 15.760."
        ),
        "duplicate_work_guards": [
            "Do not retry a Moore--Penrose or full-target Euclidean norm bound on the compact rays.",
            "Do not identify the least-norm real preimage with an integral or 0/1 graph.",
            "Do not replace the still-live higher moment and affine-box constraints by scalar energy.",
        ],
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    out = theorem_record()
    print("Prop. 15.761 full edge-Radon spectrum: proved")
    print("  least-norm exclusion of compact rays: PROVED IMPOSSIBLE")
    return out


if __name__ == "__main__":
    main()
