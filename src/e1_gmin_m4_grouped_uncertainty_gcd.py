#!/usr/bin/env python3
"""Homogeneous-form reduction for the grouped uncertainty obstruction.

This module records an exact algebraic reformulation; it does not prove the
still-open even-support inequality.  If ``S`` is a set of antipodal classes
and ``A`` is a projective direction, the nonzero affine blocks in direction
``A`` are the fibres of the squared projection ``ell_i(A)^2``.  The odd
parts of a family of split homogeneous forms detect exactly when all those
fibres have even size.
"""

from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15721 import is_prime


def _check_odd_prime(p: int) -> None:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")


def _check_support_size(p: int, s: int) -> None:
    delta_size = (p * p - 1) // 2
    if (
        not isinstance(s, int)
        or isinstance(s, bool)
        or s < 1
        or s > delta_size
    ):
        raise ValueError("support size must lie in [1,(p^2-1)/2]")


def projection_partition_factor_orders(
    s: int,
    radial_count: int,
    nonzero_block_sizes: tuple[int, ...],
) -> dict[str, object]:
    """Return the local orders of the detecting forms at one direction.

    There are ``radial_count`` classes with zero projection.  The remaining
    classes are partitioned by their nonzero squared projection values, with
    sizes ``nonzero_block_sizes``.  For

        F_i = ell_i * product_(j != i) (ell_i^2-ell_j^2),

    the vanishing order is ``2*radial_count-1`` on a radial class and
    ``block_size-1`` on a class in a nonzero block.  Consequently the
    direction is silent exactly when every nonzero block size is even, and
    this is equivalent to the direction occurring to odd order in every
    ``F_i``.
    """
    if not isinstance(s, int) or isinstance(s, bool) or s < 1:
        raise ValueError("s must be a positive integer")
    if (
        not isinstance(radial_count, int)
        or isinstance(radial_count, bool)
        or radial_count < 0
        or radial_count > s
    ):
        raise ValueError("invalid radial count")
    if any(
        not isinstance(m, int) or isinstance(m, bool) or m < 1
        for m in nonzero_block_sizes
    ):
        raise ValueError("nonzero block sizes must be positive integers")
    if radial_count + sum(nonzero_block_sizes) != s:
        raise ValueError("the radial and nonzero blocks must partition s")

    radial_order = 2 * radial_count - 1 if radial_count else None
    orders = (
        ([radial_order] * radial_count if radial_order is not None else [])
        + [m - 1 for m in nonzero_block_sizes for _ in range(m)]
    )
    silent = all(m % 2 == 0 for m in nonzero_block_sizes)
    common_odd_factor = all(order % 2 == 1 for order in orders)
    if silent != common_odd_factor:
        raise ArithmeticError("the local odd-factor criterion changed")
    return {
        "s": s,
        "radial_count": radial_count,
        "nonzero_block_sizes": nonzero_block_sizes,
        "radial_factor_order": radial_order,
        "factor_orders_by_class": tuple(orders),
        "silent_direction": silent,
        "common_odd_factor": common_odd_factor,
        "proved": True,
    }

def homogeneous_gcd_reduction(p: int, s: int) -> dict[str, object]:
    """Return the exact split-form gcd characterization for support ``s``.

    Choose representatives ``v_i`` of the antipodal support classes and put
    ``ell_i(u)=u(v_i)``.  Let ``odd(F)`` retain the projective linear factors
    of ``F`` having odd multiplicity.  Then the silent directions are
    exactly the roots of

        gcd_i odd(ell_i product_(j != i)(ell_i^2-ell_j^2)).

    In particular the number ``z`` of silent directions is the degree of
    that squarefree gcd.
    """
    _check_odd_prime(p)
    _check_support_size(p, s)
    form_degree = 2 * s - 1
    product_degree = s * form_degree
    paired_factor_degree = s + 4 * comb(s, 2)
    if product_degree != paired_factor_degree:
        raise ArithmeticError("the product-form degree identity changed")
    return {
        "p": p,
        "s": s,
        "linear_forms": "ell_i(u)=u(v_i), with v_i defined modulo sign",
        "detecting_form": (
            "F_i=ell_i*product_(j!=i)(ell_i^2-ell_j^2)"
        ),
        "detecting_form_degree": form_degree,
        "odd_part": (
            "odd(F)=product of projective linear factors of odd multiplicity"
        ),
        "exact_silent_gcd_formula": "z=deg(gcd_i odd(F_i))",
        "local_orders": {
            "radial_class_in_radial_group_r": "2*r-1",
            "class_in_nonzero_squared_projection_block_m": "m-1",
        },
        "product_identity_modulo_nonzero_scalar": (
            "product_i F_i=(product_i ell_i)*"
            "product_(i<j)(ell_i^2-ell_j^2)^2"
        ),
        "product_degree": product_degree,
        "odd_support_bound_proved": s % 2 == 1,
        "even_support_bound_proved": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def odd_support_grouped_uncertainty_theorem(p: int, s: int) -> dict[str, object]:
    """Prove ``z<=s`` when the antipodal support size ``s`` is odd."""
    _check_odd_prime(p)
    _check_support_size(p, s)
    if s % 2 != 1:
        raise ValueError("this theorem is only the odd-support branch")
    reduction = homogeneous_gcd_reduction(p, s)
    if not reduction["odd_support_bound_proved"]:
        raise ArithmeticError("the odd-support gcd consequence changed")
    return {
        "p": p,
        "s": s,
        "conclusion": "z<=s",
        "proof": (
            "If D=gcd_i odd(F_i), every factor of D occurs oddly in every "
            "F_i. Since s is odd it occurs oddly in product_i F_i. The "
            "product identity leaves only odd(product_i ell_i), of degree "
            "at most s. Hence z=deg(D)<=s."
        ),
        "equivalent_radial_parity_proof": (
            "a silent direction has radial count congruent to s modulo 2, "
            "so for odd s it contains a support class; radial directions "
            "of distinct classes account for at most s directions"
        ),
        "proved": True,
    }


def even_support_counterexample_constraints(
    p: int, s: int, z: int
) -> dict[str, object]:
    """Return necessary homogeneous-form constraints on an even counterexample.

    This function deliberately reports necessary conditions only.  In the
    first open case ``s=8, z=9``, division by the common silent factor would
    leave eight homogeneous sextics in a seven-dimensional space, hence a
    forced linear dependence.  Ruling out that dependence (including
    repeated radial factors) is not done here.
    """
    _check_odd_prime(p)
    _check_support_size(p, s)
    if s % 2 != 0:
        raise ValueError("this record is only for even support")
    if (
        not isinstance(z, int)
        or isinstance(z, bool)
        or z <= s
        or z > p + 1
    ):
        raise ValueError("need a putative counterexample with s<z<=p+1")
    form_degree = 2 * s - 1
    quotient_degree = form_degree - z
    quotient_space_dimension = max(0, quotient_degree + 1)
    forced_dependence = quotient_space_dimension < s
    if not forced_dependence:
        raise ArithmeticError("the counterexample dimension cut changed")
    return {
        "p": p,
        "s": s,
        "z": z,
        "necessary_common_squarefree_factor_degree": z,
        "detecting_form_degree": form_degree,
        "quotient_form_degree_at_most": quotient_degree,
        "quotient_form_space_dimension_at_most": quotient_space_dimension,
        "number_of_quotient_forms": s,
        "quotient_forms_forced_linearly_dependent": forced_dependence,
        "reason_even_case_does_not_follow_from_product_identity": (
            "the common factor occurs s times and is therefore a square"
        ),
        "counterexample_excluded": False,
        "grouped_uncertainty_even_branch": "OPEN",
        "residual_ii_closed": False,
        "proved": True,
    }
