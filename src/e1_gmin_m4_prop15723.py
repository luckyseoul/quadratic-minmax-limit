#!/usr/bin/env python3
"""Prop. 15.723 -- paired-cube obstruction to a middle floor-plus-two lift.

Let ``p=2m-1>=17`` be odd and let ``A`` be a nonnegative,
integer-valued quadratic on ``J(p,m)`` with parity

    A(X) = |X intersect B| + eta                         (mod 2).

In the infinity-plus-``p`` boundary shell, ``b=|B|`` is odd.  Suppose the
exact scaled mean lies only two units above the middle parity floor,

    2p E[A] = 2p+2.                                      (1)

For the genuinely middle range

    5 <= min(b,p-b) <= m-1,

(1) is impossible except in the two real equality cells

    (p,b,eta)=(17,5,1), (17,11,0).                       (2)

Both exceptions are attained by ``A=(t-3)^2`` on the smaller parity side,
so no theorem using only nonnegativity, integrality, degree two, parity and
the mean can remove them.

The proof uses paired cubes through a point of the Johnson slice.  If a
nonnegative integral quadratic ``g`` on a Boolean cube has parity
``1+sum_(i in R) z_i``, Fourier orthogonality gives ``E[g]>=1`` for
``|R|>=3``.  For ``|R|>=5``, equality would force ``g-1`` to vanish on one
parity half and hence, by the degree-two Fourier cutoff, on the whole cube.
That contradicts parity on the other half.  Cube means are half-integral,
so

    E[g] >= 3/2.                                         (3)

The paired-cube operator on degree at most two satisfies

    T A(X) = (A(X)+p E[A])/(p+1).

Under (1), this is ``1+A(X)/(p+1)``.  Applying (3) on a contact layer and
then an exact positive three-node quadrature gives a strict lower mean in
every middle cell except (2).  The small phase-zero sides ``k=5,6`` use the
exact fraction of cubes with at least five active parity coordinates.

This proposition repairs the blanket ``excess != 2`` shortcut: it proves
that shortcut in the general middle range while preserving the two genuine
``p=17`` exceptions.  It does not by itself close those exceptions, the
infinity-plus-p shell, residual (ii), Type I, or the limit.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _check_parameters(p: int, b: int | None = None, phase: int | None = None) -> None:
    if p < 17 or p % 2 == 0:
        raise ValueError("need an odd parameter p>=17")
    if b is not None and not (1 <= b <= p and b % 2 == 1):
        raise ValueError("the infinity-plus-p shell requires odd 1<=b<=p")
    if phase is not None and phase not in (0, 1):
        raise ValueError("phase must be zero or one")


def cube_parity_mass_ledger(active_coordinates: int) -> dict[str, object]:
    """Return the sharp cube-mean statement used by the paired operator."""
    if active_coordinates < 0:
        raise ValueError("active-coordinate count must be nonnegative")
    if active_coordinates >= 5:
        lower = Fraction(3, 2)
    elif active_coordinates >= 3:
        lower = Fraction(1)
    else:
        lower = Fraction(0)
    sharp_example = None
    if active_coordinates in (5, 6):
        # E[(S-3)^2] for S~Binomial(r,1/2).
        r = active_coordinates
        sharp_example = Fraction(r, 4) + Fraction(r - 6, 2) ** 2
        if sharp_example != Fraction(3, 2):
            raise ArithmeticError("sharp paired-cube example changed")
    return {
        "active_coordinates": active_coordinates,
        "mean_lies_in_half_integers": True,
        "fourier_orthogonality_floor": str(Fraction(1)) if active_coordinates >= 3 else None,
        "strict_equality_obstruction_for_degree_two": active_coordinates >= 5,
        "mean_lower_bound": str(lower),
        "sharp_example": (
            None
            if sharp_example is None
            else {"polynomial": "(sum_R z_i-3)^2", "mean": str(sharp_example)}
        ),
        "proved": True,
    }


def paired_cube_operator_value(p: int, value_at_x: Fraction, mean: Fraction) -> Fraction:
    """Degree-two paired-cube operator ``(A(X)+p EA)/(p+1)``."""
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd and at least three")
    return (value_at_x + p * mean) / (p + 1)


def paired_cube_operator_audit(p: int) -> dict[str, object]:
    """Brute-force the operator identity on a basis of quadratic monomials."""
    if p not in (5, 7):
        raise ValueError("the finite operator audit is intentionally p=5 or 7")
    m = (p + 1) // 2
    points = tuple(itertools.combinations(range(p), m))
    anchor = frozenset(points[0])
    complement = tuple(i for i in range(p) if i not in anchor)
    cube_points: list[frozenset[int]] = []
    for unmatched in anchor:
        paired_left = tuple(sorted(anchor - {unmatched}))
        for permutation in itertools.permutations(complement):
            pairs = tuple(zip(paired_left, permutation))
            for mask in range(1 << (m - 1)):
                row = {unmatched}
                for index, (left, right) in enumerate(pairs):
                    row.add(right if (mask >> index) & 1 else left)
                cube_points.append(frozenset(row))
    # The operator averages uniformly over unmatched points and bijections;
    # repeated cube vertices retain their natural multiplicity.
    monomials = [(i,) for i in range(p)] + list(itertools.combinations(range(p), 2))
    maximum_error = Fraction(0)
    for monomial in monomials:
        def evaluate(row: frozenset[int]) -> int:
            return int(all(index in row for index in monomial))

        lhs = Fraction(sum(evaluate(row) for row in cube_points), len(cube_points))
        at_anchor = Fraction(evaluate(anchor))
        global_mean = Fraction(sum(evaluate(frozenset(row)) for row in points), len(points))
        rhs = paired_cube_operator_value(p, at_anchor, global_mean)
        maximum_error = max(maximum_error, abs(lhs - rhs))
    return {
        "p": p,
        "m": m,
        "basis_monomials_checked": len(monomials),
        "cube_samples_with_multiplicity": len(cube_points),
        "maximum_error": str(maximum_error),
        "identity": "T A(X)=(A(X)+p*E[A])/(p+1)",
        "proved_for_degree_at_most_two": maximum_error == 0,
    }


def paired_cube_operator_symbolic_ledger() -> dict[str, object]:
    """Give the all-parameter monomial proof of the paired-cube identity.

    The finite audits above are useful regression tests, but they are not the
    proof for arbitrary ``p``.  Constants, coordinates, and squarefree pairs
    span the quadratic functions on a slice.  Conditioning on whether a
    point is the unmatched member gives the displayed probabilities, which
    agree identically with ``(f(X)+p E[f])/(p+1)`` for each monomial type.
    """
    return {
        "basis": "1, x_i, and x_i*x_j for i<j",
        "write_p_as": "p=2m-1 and p+1=2m",
        "cube_probabilities": {
            "constant": "1",
            "coordinate_in_X": "(m+1)/(2m)",
            "coordinate_outside_X": "1/2",
            "pair_both_in_X": "(m+2)/(4m)",
            "pair_not_both_in_X": "1/4",
        },
        "right_hand_side_checks": {
            "coordinate_in_X": "(1+p*(m/p))/(p+1)=(m+1)/(2m)",
            "coordinate_outside_X": "p*(m/p)/(p+1)=1/2",
            "pair_both_in_X": (
                "(1+p*(m(m-1)/(p(p-1))))/(p+1)=(m+2)/(4m)"
            ),
            "pair_not_both_in_X": (
                "p*(m(m-1)/(p(p-1)))/(p+1)=1/4"
            ),
        },
        "linearity_closes_all_degree_at_most_two_functions": True,
        "proved": True,
    }


def _hypergeometric_moments(p: int, k: int) -> tuple[Fraction, Fraction]:
    """First two moments of ``t=|X intersect B|`` on ``J(p,(p+1)/2)``."""
    return (
        Fraction(k * (p + 1), 2 * p),
        Fraction(k * (k + 1) * (p + 1), 4 * p),
    )


def _three_node_moment_weights(
    p: int, k: int, nodes: tuple[int, int, int]
) -> tuple[Fraction, Fraction, Fraction]:
    """Unique weights on three nodes matching moments through degree two."""
    if len(set(nodes)) != 3:
        raise ValueError("quadrature nodes must be distinct")
    mean, second = _hypergeometric_moments(p, k)
    weights = []
    for index, node in enumerate(nodes):
        other = [nodes[j] for j in range(3) if j != index]
        left, right = other
        weights.append(
            Fraction(
                second - (left + right) * mean + left * right,
                (node - left) * (node - right),
            )
        )
    return tuple(weights)  # type: ignore[return-value]


def _quadrature_moment_match(
    p: int,
    k: int,
    nodes: tuple[int, int, int],
    weights: tuple[Fraction, Fraction, Fraction],
) -> bool:
    mean, second = _hypergeometric_moments(p, k)
    targets = (Fraction(1), mean, second)
    return all(
        sum(weight * node**degree for node, weight in zip(nodes, weights))
        == targets[degree]
        for degree in range(3)
    )


def phase_one_contact_quadrature(p: int, k: int) -> dict[str, object]:
    """Exact positive contact quadrature with the upgraded ``t=0`` node."""
    _check_parameters(p)
    m = (p + 1) // 2
    if not 5 <= k <= m - 1:
        raise ValueError("need 5<=k<=m-1")
    near = 2 * (k // 4)
    nodes = (0, near, near + 2)
    weights = _three_node_moment_weights(p, k, nodes)
    zero_weight = weights[0]
    moment_match = _quadrature_moment_match(p, k, nodes, weights)
    positive = all(weight >= 0 for weight in weights)
    contact = all(node % 2 == 0 for node in nodes)
    formula_match = zero_weight == phase_one_zero_weight(p, k)
    gap = (m - 1) * zero_weight - Fraction(1, p)
    exact = bool(moment_match and positive and contact and formula_match)
    if not exact:
        raise ArithmeticError("phase-one three-node quadrature failed")
    return {
        "p": p,
        "k": k,
        "phase": 1,
        "nodes": nodes,
        "weights": weights,
        "moments_through_degree_two_match": moment_match,
        "positive_weights": positive,
        "all_nodes_are_parity_contacts": contact,
        "upgraded_node": 0,
        "active_parity_coordinates_at_upgraded_node": k,
        "paired_cube_forces_A_at_upgraded_node_at_least": m,
        "upgraded_node_weight": zero_weight,
        "weight_formula_matches": formula_match,
        "lower_mean": Fraction(1) + (m - 1) * zero_weight,
        "target_mean": Fraction(p + 1, p),
        "gap": gap,
        "proved": exact,
    }


def _phase_zero_large_nodes(k: int) -> tuple[int, int, int]:
    residue = k % 4
    s = k // 4
    near = 2 * s + 1 if residue == 3 else 2 * s - 1
    far = k if k % 2 else k - 1
    return near, near + 2, far


def phase_zero_large_contact_quadrature(p: int, k: int) -> dict[str, object]:
    """Exact positive quadrature upgrading the far contact, for ``k>=7``."""
    _check_parameters(p)
    m = (p + 1) // 2
    if not 7 <= k <= m - 1:
        raise ValueError("need 7<=k<=m-1")
    nodes = _phase_zero_large_nodes(k)
    weights = _three_node_moment_weights(p, k, nodes)
    moment_match = _quadrature_moment_match(p, k, nodes, weights)
    positive = all(weight >= 0 for weight in weights)
    contact = all(node % 2 == 1 for node in nodes)
    far_weight = weights[2]
    gap = (m - 1) * far_weight - Fraction(1, p)
    denominator = 8 * p * (nodes[2] - nodes[0]) * (nodes[2] - nodes[1])
    formula_match = gap * denominator == phase_zero_large_numerator(p, k)
    exact = bool(moment_match and positive and contact and formula_match)
    if not exact:
        raise ArithmeticError("phase-zero far-contact quadrature failed")
    return {
        "p": p,
        "k": k,
        "phase": 0,
        "nodes": nodes,
        "weights": weights,
        "moments_through_degree_two_match": moment_match,
        "positive_weights": positive,
        "all_nodes_are_parity_contacts": contact,
        "upgraded_node": nodes[2],
        # At t=k (k odd), leaving a B-point unmatched gives k-1 active
        # pairs.  At t=k-1 (k even), the worst case also pairs the unique
        # outside B-point to an inside B-point, leaving k-3 active pairs.
        "minimum_active_parity_coordinates_there": k - 1 if k % 2 else k - 3,
        "paired_cube_forces_A_at_upgraded_node_at_least": m,
        "upgraded_node_weight": far_weight,
        "lower_mean": Fraction(1) + (m - 1) * far_weight,
        "target_mean": Fraction(p + 1, p),
        "gap": gap,
        "positive_gap_numerator": phase_zero_large_numerator(p, k),
        "gap_denominator": denominator,
        "numerator_formula_matches": formula_match,
        "proved": exact,
    }


def phase_zero_small_contact_quadrature(p: int, k: int) -> dict[str, object]:
    """Endpoint quadrature and exact good-cube fraction for ``k=5,6``."""
    _check_parameters(p)
    if k not in (5, 6):
        raise ValueError("the small phase-zero branch is k=5 or k=6")
    m = (p + 1) // 2
    if k > m - 1:
        raise ValueError("need k<=m-1")
    nodes = (1, 3, 5)
    weights = _three_node_moment_weights(p, k, nodes)
    moment_match = _quadrature_moment_match(p, k, nodes, weights)
    positive = all(weight >= 0 for weight in weights)
    endpoint_weight = weights[0] + weights[2]
    expected_endpoint_weight = Fraction(3 * (p - 5), 8 * p)
    good_cube_fraction = Fraction(m - 5, m)
    # Every endpoint cube has at least three active coordinates.  The exact
    # fraction below has at least five, so its mean improves from 1 to 3/2.
    raw_endpoint_floor = 2 * m * good_cube_fraction / 2
    endpoint_floor = int(raw_endpoint_floor)
    if endpoint_floor % 2 == 0:
        endpoint_floor += 1
    endpoint_improvement = endpoint_floor - 1
    gap = endpoint_weight * endpoint_improvement - Fraction(1, p)
    numerator = phase_zero_small_gap_numerator(p)
    formula_match = gap * 16 * p == numerator
    exact = bool(
        moment_match
        and positive
        and endpoint_weight == expected_endpoint_weight
        and formula_match
    )
    if not exact:
        raise ArithmeticError("phase-zero endpoint quadrature failed")
    return {
        "p": p,
        "k": k,
        "phase": 0,
        "nodes": nodes,
        "weights": weights,
        "moments_through_degree_two_match": moment_match,
        "positive_weights": positive,
        "upgraded_nodes": (1, 5),
        "all_endpoint_cubes_have_at_least_three_active_coordinates": True,
        "fraction_with_at_least_five_active_coordinates": good_cube_fraction,
        "forced_odd_endpoint_floor": endpoint_floor,
        "combined_endpoint_weight": endpoint_weight,
        "combined_endpoint_weight_formula": expected_endpoint_weight,
        "lower_mean": Fraction(1) + endpoint_weight * endpoint_improvement,
        "target_mean": Fraction(p + 1, p),
        "gap": gap,
        "positive_gap_numerator": numerator,
        "gap_denominator": 16 * p,
        "numerator_formula_matches": formula_match,
        "proved": exact,
    }


def smaller_parity_side(p: int, b: int, phase: int) -> dict[str, int]:
    """Complement the parity set when that makes it smaller on ``J(p,m)``."""
    _check_parameters(p, b, phase)
    m = (p + 1) // 2
    if b <= m - 1:
        return {"k": b, "phase": phase, "complemented": 0}
    return {"k": p - b, "phase": (phase + m) & 1, "complemented": 1}


def phase_one_zero_weight(p: int, k: int) -> Fraction:
    """Weight at the zero contact node in the positive phase-one quadrature."""
    residue = k % 4
    if residue == 0:
        return Fraction(p - k - 3, p * (k + 4))
    if residue == 1:
        return Fraction(-k * k + k * p - k - 3 * p, p * (k - 1) * (k + 3))
    if residue == 2:
        return Fraction(-k * k + k * p + k - 4 * p, p * (k - 2) * (k + 2))
    return Fraction(p - k, p * (k + 1))


def phase_one_gap(p: int, k: int) -> Fraction:
    """Excess of the contact quadrature over the target mean ``1+1/p``."""
    m = (p + 1) // 2
    return (m - 1) * phase_one_zero_weight(p, k) - Fraction(1, p)


def phase_one_gap_numerator(p: int, k: int) -> int:
    """Numerator of :func:`phase_one_gap` over its positive denominator."""
    residue = k % 4
    if residue == 0:
        return (p - 1) * (p - k - 3) - 2 * (k + 4)
    if residue == 1:
        return (p - 1) * (-k * k + k * p - k - 3 * p) - 2 * (k - 1) * (k + 3)
    if residue == 2:
        return (p - 1) * (-k * k + k * p + k - 4 * p) - 2 * (k - 2) * (k + 2)
    return (p - 1) * (p - k) - 2 * (k + 1)


def phase_one_gap_denominator(p: int, k: int) -> int:
    """Positive denominator paired with :func:`phase_one_gap_numerator`."""
    residue = k % 4
    factor = (
        k + 4
        if residue == 0
        else (k - 1) * (k + 3)
        if residue == 1
        else (k - 2) * (k + 2)
        if residue == 2
        else k + 1
    )
    return 2 * p * factor


def phase_zero_large_numerator(p: int, k: int) -> int:
    """Positive numerator for the phase-zero contact-layer gap, ``k>=7``."""
    residue = k % 4
    return (
        -k * k * p - k * k + k * p * p + 7 * k - 4 * p * p + 4 * p
        if residue == 0
        else -k * k * p - k * k + k * p * p + 2 * k * p - 7 * k - 3 * p * p + 3 * p + 6
        if residue == 1
        else -k * k * p - k * k + k * p * p + 4 * k * p - 5 * k + 8
        if residue == 2
        else -k * k * p - k * k + k * p * p - 2 * k * p + 5 * k - 3 * p * p + 3 * p + 6
    )


def phase_zero_small_gap_numerator(p: int) -> int:
    """Numerator of the ``k=5,6`` endpoint-pairing lower gap."""
    if p % 4 == 1:
        return 3 * (p - 5) * (p - 9) - 16
    return 3 * (p - 5) * (p - 11) - 16


def universal_gap_positivity_certificate() -> dict[str, object]:
    """Elementary all-parameter sign proof for every quadrature gap.

    Put ``d=p-(2k+1)>=0``.  The four phase-one numerators factor as

    ``(d+k-4)(d+2k+2)``,
    ``(d+2k+2)(d(k-3)+k^2-8k+3)``,
    ``(d+2k+2)(d(k-4)+k^2-8k+4)``, and
    ``(d+k-1)(d+2k+2)``.

    In residues one and two, the small values ``k=5,6`` give respectively
    ``d=6,4`` at equality; ``p>=17`` and parity make every later ``d`` at
    least two larger.  For ``k>=9,10`` the remaining quadratic factors are
    already positive and increasing.

    The phase-zero large numerators, after the same substitution, are the
    four polynomials recorded below.  All coefficients relevant on the
    indicated ranges are positive.  The only superficially negative
    constant is the residue-three ``k=7`` value ``-8``; there ``p>=17``
    forces ``d>=2``, giving numerator ``128``.  This is an infinite-range
    sign argument, not a finite prime scan.
    """
    phase_one_factors = {
        "k_mod_4=0": "(d+k-4)(d+2k+2), k>=8",
        "k_mod_4=1": (
            "(d+2k+2)(d(k-3)+k^2-8k+3); "
            "k=5 has equality only at d=6, k>=9 is strict"
        ),
        "k_mod_4=2": (
            "(d+2k+2)(d(k-4)+k^2-8k+4); "
            "k=6 has equality only at d=4, k>=10 is strict"
        ),
        "k_mod_4=3": "(d+k-1)(d+2k+2), k>=7",
    }
    phase_zero_polynomials = {
        "k_mod_4=0": (
            "(k-4)d^2+(3k^2-14k-4)d+2k^2(k-7), k>=8"
        ),
        "k_mod_4=1": (
            "(k-3)d^2+(3k^2-8k-3)d+2k^3-6k^2-10k+6, k>=9"
        ),
        "k_mod_4=2": (
            "k d^2+(3k^2+6k)d+2k^3+10k^2+8, k>=10"
        ),
        "k_mod_4=3": (
            "(k-3)d^2+(3k^2-12k-3)d+2k^3-14k^2-2k+6; "
            "k=7,d>=2 and k>=11"
        ),
    }
    # These are the endpoint/minimum checks used in the monotonicity claims
    # above.  Keeping them executable prevents a transcription from silently
    # changing the sign proof.
    phase_one_thresholds = (
        9 * 9 - 8 * 9 + 3 > 0,
        10 * 10 - 8 * 10 + 4 > 0,
        phase_one_gap_numerator(17, 5) == 0,
        phase_one_gap_numerator(17, 6) == 0,
        phase_one_gap_numerator(19, 5) > 0,
        phase_one_gap_numerator(19, 6) > 0,
    )
    phase_zero_thresholds = (
        3 * 8 * 8 - 14 * 8 - 4 > 0,
        2 * 8 * 8 * (8 - 7) > 0,
        2 * 9**3 - 6 * 9**2 - 10 * 9 + 6 > 0,
        phase_zero_large_numerator(17, 7) == 128,
        2 * 11**3 - 14 * 11**2 - 2 * 11 + 6 > 0,
    )
    small_thresholds = (
        phase_zero_small_gap_numerator(17) > 0,
        phase_zero_small_gap_numerator(19) > 0,
    )
    quadrature_weight_sign_proof = {
        "phase_one": (
            "write k=4s+r; the three weight numerators are positive for "
            "s>=2 after p>=2k+1, while (k,p)=(5,>=17),(6,>=17) "
            "are the two direct endpoint checks; the residue-three middle "
            "weight is zero and the other two are positive"
        ),
        "phase_zero_large": (
            "with the listed nodes, substitution p=2k+1+d gives positive "
            "weight numerators for k>=7; the minima are k=8,9,10 and "
            "(k,d)=(7,2) in the four residue classes"
        ),
        "phase_zero_small": (
            "nodes 1,3,5 have positive weights and combined endpoint "
            "weight 3(p-5)/(8p)"
        ),
        "proved": True,
    }
    proved = bool(
        all(phase_one_thresholds)
        and all(phase_zero_thresholds)
        and all(small_thresholds)
        and quadrature_weight_sign_proof["proved"]
    )
    return {
        "parameterization": "d=p-(2k+1)>=0; d is even",
        "phase_one_gap_factorizations": phase_one_factors,
        "phase_zero_large_gap_polynomials": phase_zero_polynomials,
        "phase_zero_small_monotonicity": {
            "p_mod_4=1": "3(p-5)(p-9)-16, increasing from p=17",
            "p_mod_4=3": "3(p-5)(p-11)-16, increasing from p=19",
        },
        "quadrature_weight_sign_proof": quadrature_weight_sign_proof,
        "only_zero_gaps": [(17, 5, 1), (17, 6, 1)],
        "proved": proved,
    }


def middle_floor_plus_two_cell(p: int, b: int, phase: int) -> dict[str, object]:
    """Classify one odd-profile cell under ``2p E[A]=2p+2``."""
    _check_parameters(p, b, phase)
    m = (p + 1) // 2
    reduced = smaller_parity_side(p, b, phase)
    k = reduced["k"]
    reduced_phase = reduced["phase"]
    applicable = 5 <= k <= m - 1
    if not applicable:
        return {
            "p": p,
            "b": b,
            "phase": phase,
            "smaller_side": reduced,
            "applicable_middle_cell": False,
            "excluded": False,
            "exceptional_equality": False,
            "proved": True,
        }

    method: str
    gap: Fraction | None = None
    numerator: int | None = None
    quadrature: dict[str, object]
    if reduced_phase == 1:
        method = "phase-one zero-contact positive quadrature"
        quadrature = phase_one_contact_quadrature(p, k)
        gap = Fraction(quadrature["gap"])
        excluded = gap > 0
        equality = gap == 0
    elif k >= 7:
        method = "phase-zero far-contact positive quadrature"
        quadrature = phase_zero_large_contact_quadrature(p, k)
        numerator = int(quadrature["positive_gap_numerator"])
        excluded = numerator > 0
        equality = numerator == 0
    else:
        method = "phase-zero k=5,6 mixed active-cube count"
        quadrature = phase_zero_small_contact_quadrature(p, k)
        numerator = int(quadrature["positive_gap_numerator"])
        excluded = numerator > 0
        equality = numerator == 0

    expected_exception = (p, b, phase) in ((17, 5, 1), (17, 11, 0))
    if equality != expected_exception or excluded == expected_exception:
        raise ArithmeticError("paired-cube middle-cell classification changed")
    if not (excluded or expected_exception):
        raise ArithmeticError("an unclassified middle floor-plus-two cell survived")
    return {
        "p": p,
        "b": b,
        "phase": phase,
        "smaller_side": reduced,
        "applicable_middle_cell": True,
        "scaled_mean_assumption": 2 * p + 2,
        "target_unscaled_mean": str(Fraction(p + 1, p)),
        "method": method,
        "quadrature_certificate": quadrature,
        "gap": None if gap is None else str(gap),
        "gap_numerator": numerator,
        "excluded": excluded,
        "exceptional_equality": expected_exception,
        "proved": excluded or expected_exception,
    }


def _check_backward_floor_parameters(p: int, b: int, phase: int) -> None:
    if p < 17 or p % 2 == 0:
        raise ValueError("need an odd parameter p>=17")
    if not 0 <= b <= p:
        raise ValueError("need 0<=b<=p")
    if phase not in (0, 1):
        raise ValueError("phase must be zero or one")


def backward_floor_plus_two_cell(p: int, b: int, phase: int) -> dict[str, object]:
    """Classify an older floor ledger's ``excess == 2`` cell.

    The older all-finite ledgers use even ``b``.  Proposition 15.723 is
    stated in the odd representative appropriate to the infinity-plus-p
    shell, so an even cell must first be complemented.  On ``J(p,m)`` this
    sends ``(b, phase)`` to ``(p-b, phase+m mod 2)``.  Endpoint cells with
    reduced size at most two retain the earlier baseline obstruction;
    reduced sizes three and four are not covered by either that obstruction
    or Proposition 15.723 and therefore must remain in a relaxed ledger.
    """
    _check_backward_floor_parameters(p, b, phase)

    m = (p + 1) // 2
    complemented = b % 2 == 0
    odd_b = p - b if complemented else b
    odd_phase = phase ^ (m & 1) if complemented else phase
    reduced_size = min(b, p - b)

    middle_row: dict[str, object] | None = None
    if reduced_size <= 2:
        forbidden = True
        classification = "endpoint baseline obstruction"
    elif reduced_size <= 4:
        forbidden = False
        classification = "OPEN reduced-size-three-or-four cell"
    else:
        middle_row = middle_floor_plus_two_cell(p, odd_b, odd_phase)
        if not middle_row["applicable_middle_cell"]:
            raise ArithmeticError("normalized middle cell left Proposition 15.723")
        forbidden = bool(middle_row["excluded"])
        classification = (
            "Proposition 15.723 middle-cell obstruction"
            if forbidden
            else "Proposition 15.723 equality exception"
        )

    return {
        "p": p,
        "input_b": b,
        "input_phase": phase,
        "reduced_size": reduced_size,
        "normalized_odd_b": odd_b,
        "normalized_odd_phase": odd_phase,
        "complemented": complemented,
        "classification": classification,
        "floor_plus_two_forbidden": forbidden,
        "floor_plus_two_exclusion_proved": forbidden,
        "admissible_in_relaxed_ledger": not forbidden,
        "exceptional_equality": bool(
            middle_row is not None and middle_row["exceptional_equality"]
        ),
        "middle_cell": middle_row,
        "proved": True,
    }


def floor_excess_admissible(p: int, b: int, phase: int, excess: int) -> bool:
    """Return whether a floor excess survives the proved obstructions.

    Negative excess is below the symbolic floor.  Every nonnegative excess
    other than two is admissible at this stage.  The two-unit case is routed
    through :func:`backward_floor_plus_two_cell`, including complement and
    phase normalization.
    """
    _check_backward_floor_parameters(p, b, phase)
    if excess < 0:
        return False
    if excess != 2:
        return True
    return bool(backward_floor_plus_two_cell(p, b, phase)["admissible_in_relaxed_ledger"])


def exceptional_quadratic_witness(p: int, b: int, phase: int) -> dict[str, object]:
    """Verify the two exact ``p=17`` equality quadratics."""
    row = middle_floor_plus_two_cell(p, b, phase)
    if not row["exceptional_equality"]:
        raise ValueError("this cell is not an equality exception")
    reduced = row["smaller_side"]
    k = int(reduced["k"])
    m = (p + 1) // 2
    first = Fraction(k * m, p)
    second = Fraction(k * (k + 1) * (p + 1), 4 * p)
    mean = second - 6 * first + 9
    scaled = 2 * p * mean
    proved = mean == Fraction(18, 17) and scaled == 2 * p + 2
    if not proved:
        raise ArithmeticError("exceptional quadratic mean changed")
    return {
        "p": p,
        "b": b,
        "phase": phase,
        "smaller_side_size": k,
        "smaller_side_phase": int(reduced["phase"]),
        "quadratic": "A(X)=(|X intersect C|-3)^2",
        "E_t": str(first),
        "E_t_squared": str(second),
        "E_A": str(mean),
        "scaled_mean": int(scaled),
        "required_parity": True,
        "proved": proved,
    }


def theorem_middle_floor_plus_two() -> dict[str, object]:
    """Package the all-parameter proof and finite regression audits."""
    universal = universal_gap_positivity_certificate()
    operator_symbolic = paired_cube_operator_symbolic_ledger()
    sample_parameters = tuple(range(17, 204, 2))
    sample_rows = []
    exceptions = []
    for p in sample_parameters:
        for b in range(1, p + 1, 2):
            for phase in (0, 1):
                row = middle_floor_plus_two_cell(p, b, phase)
                if row["applicable_middle_cell"]:
                    sample_rows.append(row)
                if row["exceptional_equality"]:
                    exceptions.append((p, b, phase))
    expected = [(17, 5, 1), (17, 11, 0)]
    backward_exceptions = [
        backward_floor_plus_two_cell(17, 12, 0),
        backward_floor_plus_two_cell(17, 6, 1),
    ]
    backward_classifier_proved = bool(
        all(bool(row["exceptional_equality"]) for row in backward_exceptions)
        and [
            (int(row["normalized_odd_b"]), int(row["normalized_odd_phase"]))
            for row in backward_exceptions
        ]
        == [(5, 1), (11, 0)]
        and all(
            backward_floor_plus_two_cell(19, b, 0)["floor_plus_two_forbidden"]
            for b in (0, 1, 2)
        )
        and all(
            not backward_floor_plus_two_cell(23, b, 0)["floor_plus_two_forbidden"]
            for b in (3, 4, 19, 20)
        )
    )
    proved = bool(
        universal["proved"]
        and universal["only_zero_gaps"] == [(17, 5, 1), (17, 6, 1)]
        and operator_symbolic["proved"]
        and exceptions == expected
        and backward_classifier_proved
        and all(bool(row["proved"]) for row in sample_rows)
        and all(
            paired_cube_operator_audit(p)["proved_for_degree_at_most_two"]
            for p in (5, 7)
        )
    )
    return {
        "prop": "15.723",
        "title": "Paired-cube obstruction to middle floor-plus-two lifts",
        "proved": proved,
        "scope": (
            "all odd p>=17, odd b from the infinity-plus-p boundary shell, "
            "5<=min(b,p-b)<=m-1"
        ),
        "cube_lemma": {
            "active_at_least_3": "mean at least 1",
            "active_at_least_5": "mean at least 3/2",
            "sharp_at_5_and_6": True,
        },
        "paired_cube_operator": "T A(X)=(A(X)+p*E[A])/(p+1)",
        "paired_cube_operator_symbolic_proof": operator_symbolic,
        "universal_gap_and_quadrature_certificate": universal,
        "scaled_mean": "2p*E[A]=2p+2",
        "classification": {
            "general_middle_cells": "EXCLUDED",
            "exact_exceptions": [list(row) for row in expected],
            "exceptions_are_real_quadratics": True,
        },
        "backward_ledger_classifier": {
            "proved": backward_classifier_proved,
            "endpoint_reduced_sizes": {
                "sizes": [0, 1, 2],
                "status": "EXCLUDED_BY_BASELINE_OBSTRUCTION",
            },
            "uncovered_reduced_sizes": {
                "sizes": [3, 4],
                "status": "ADMISSIBLE_IN_RELAXED_LEDGER",
            },
            "middle_reduced_sizes": "ROUTE_THROUGH_PROPOSITION_15_723",
            "even_p17_exception_cells": [[12, 0], [6, 1]],
            "normalized_odd_exception_cells": [[5, 1], [11, 0]],
        },
        "exception_witnesses": [
            exceptional_quadratic_witness(*row) for row in expected
        ],
        "operator_audits": {
            str(p): paired_cube_operator_audit(p) for p in (5, 7)
        },
        "formula_audit_range": "regression only: all odd parameters 17<=p<=203",
        "formula_audit_cells": len(sample_rows),
        "infinity_plus_p_shell_closed": False,
        "residual_ii": False,
        "type_I": False,
        "limit_exists": False,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main() -> dict[str, object]:
    theorem = theorem_middle_floor_plus_two()
    if theorem["proved"] is not True:
        raise ArithmeticError("Proposition 15.723 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15723.json"
    destination.write_text(json.dumps(_jsonable(theorem), indent=2) + "\n")
    print("Prop 15.723 paired-cube floor-plus-two theorem: proved")
    print("  exact exceptions: p=17 (b,phase)=(5,1),(11,0)")
    print(f"  wrote {destination}")
    return theorem


if __name__ == "__main__":
    main()
