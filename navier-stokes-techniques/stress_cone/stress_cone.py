"""Admissible stress cone: positive-cone representation of a target stress.

Implements the technique from OpenAI's "Finite time blowup for Navier-Stokes"
(2026), Section 3.2 / Proposition 7.5 and Figure 4 ("Positive representation
of the target stress").

The construction needs two oscillatory wave families whose averaged quadratic
momentum flux realizes a prescribed annular stress T. Each wave family i
contributes a fixed "covariance direction" v_i in R^2 (its azimuthal- and
axial-momentum-flux ratio per unit squared amplitude); the required stress
must lie in the *interior* of the cone spanned by v_1, v_2 so that both
squared amplitudes c_1, c_2 come out strictly positive:

    T = c_1 v_1 + c_2 v_2,   c_1, c_2 > 0.

If T sits on the cone boundary or outside it, no nonnegative combination of
the two wave families can realize it, and the construction (as posed) fails
for that pair of directions -- exactly the failure mode the paper's profile
construction (Theorem 4.6(ii), Lemma 4.5) is built to avoid.

This module implements the general primitive: given any finite set of
covariance directions in R^d and a target vector, decide cone membership and
return a nonnegative (here: interior, i.e. strictly positive when possible)
representation. The 2-family / R^2 case is solved in closed form, matching
the paper's construction exactly; the general case falls back to a linear
program (scipy.optimize.linprog) so the same tool applies if you want to
explore more than two wave families or higher-dimensional flux targets.

This is a reusable *algorithmic* implementation of the technique. It does not
reproduce the paper's actual covariance-direction formulas (those depend on
the full inner/annulus/exterior profile construction in Sections 4-7); you
supply v_1, v_2 (or v_1..v_k) from whatever wave-family parameterization you
are working with.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


class NotInConeError(ValueError):
    """Raised when the target stress is not a positive combination of the
    supplied covariance directions -- i.e. it is not in the admissible cone."""


@dataclass(frozen=True)
class ConeRepresentation:
    """T = sum_i coefficients[i] * directions[i], with coefficients >= 0."""

    coefficients: np.ndarray
    directions: np.ndarray
    target: np.ndarray
    interior: bool  # True iff every coefficient is strictly positive

    def residual(self) -> np.ndarray:
        return self.target - self.coefficients @ self.directions


def two_family_representation(
    v1: np.ndarray, v2: np.ndarray, target: np.ndarray, *, tol: float = 1e-12
) -> ConeRepresentation:
    """Closed-form solve of T = c1*v1 + c2*v2 in R^2 (the paper's exact case).

    This is Figure 4: v1, v2 in R^2 are the two wave families' covariance
    directions, T is the annular stress the background profile must supply.
    Solves the 2x2 linear system directly; raises NotInConeError if the
    unique solution has a non-positive coefficient (T is outside the open
    cone) or if v1, v2 are (near-)collinear so no cone is spanned at all.
    """
    v1 = np.asarray(v1, dtype=float).reshape(2)
    v2 = np.asarray(v2, dtype=float).reshape(2)
    target = np.asarray(target, dtype=float).reshape(2)

    M = np.column_stack([v1, v2])  # 2x2
    det = np.linalg.det(M)
    if abs(det) < tol:
        raise NotInConeError(
            f"covariance directions are (near-)collinear (det={det:.3e}); "
            "no open cone is spanned, so no target stress off that line is realizable"
        )
    c = np.linalg.solve(M, target)
    interior = bool(np.all(c > tol))
    if not interior:
        raise NotInConeError(
            f"target stress {target} needs coefficients {c}; "
            "at least one is non-positive, so it is not in the interior of the cone"
        )
    return ConeRepresentation(coefficients=c, directions=M.T, target=target, interior=interior)


def general_cone_representation(
    directions: np.ndarray, target: np.ndarray, *, tol: float = 1e-9
) -> ConeRepresentation:
    """General k-direction, d-dimensional nonnegative representation via LP.

    Minimizes the L1 norm of the coefficients subject to
    sum_i c_i * directions[i] == target, c_i >= 0. Feasibility of the LP is
    exactly cone membership; a strictly positive optimal solution (all
    coefficients > tol) certifies interior membership, matching the paper's
    requirement that every wave family contributes with strictly positive
    weight.
    """
    from scipy.optimize import linprog

    directions = np.asarray(directions, dtype=float)
    target = np.asarray(target, dtype=float).reshape(-1)
    k, d = directions.shape
    if target.shape[0] != d:
        raise ValueError(f"target dimension {target.shape[0]} != direction dimension {d}")

    # minimize sum(c) subject to directions.T @ c == target, c >= 0
    res = linprog(
        c=np.ones(k),
        A_eq=directions.T,
        b_eq=target,
        bounds=[(0, None)] * k,
        method="highs",
    )
    if not res.success:
        raise NotInConeError(
            f"target {target} is not a nonnegative combination of the {k} "
            f"supplied directions (LP status: {res.message})"
        )
    coeffs = res.x
    interior = bool(np.all(coeffs > tol))
    return ConeRepresentation(
        coefficients=coeffs, directions=directions, target=target, interior=interior
    )


def _self_check() -> None:
    """Reproduces the qualitative picture in Figure 4 and checks edge cases."""
    # Interior case: two independent directions, target strictly inside the cone.
    v1 = np.array([1.0, 0.3])
    v2 = np.array([0.2, 1.0])
    T = np.array([0.6, 0.5])
    rep = two_family_representation(v1, v2, T)
    assert rep.interior
    assert np.allclose(rep.coefficients @ rep.directions, T)
    assert np.linalg.norm(rep.residual()) < 1e-10

    # Boundary/exterior case: target outside the cone must raise.
    T_bad = np.array([-1.0, 0.5])
    try:
        two_family_representation(v1, v2, T_bad)
    except NotInConeError:
        pass
    else:
        raise AssertionError("expected NotInConeError for a target outside the cone")

    # Collinear directions: no cone at all.
    try:
        two_family_representation(v1, 2.0 * v1, T)
    except NotInConeError:
        pass
    else:
        raise AssertionError("expected NotInConeError for collinear directions")

    # General LP path should agree with the closed form on the same interior case.
    rep_lp = general_cone_representation(np.stack([v1, v2]), T)
    assert rep_lp.interior
    assert np.allclose(rep_lp.coefficients, rep.coefficients, atol=1e-6)

    # Three-direction case in R^2: strictly more slack, should still find an
    # interior representation using a subset with positive weight on all three
    # only if genuinely needed; LP will legitimately zero out a redundant one.
    v3 = np.array([0.5, 0.5])
    rep3 = general_cone_representation(np.stack([v1, v2, v3]), T)
    assert np.allclose(rep3.coefficients @ rep3.directions, T, atol=1e-8)

    print("stress_cone self-check: all assertions passed")
    print(f"  interior 2-family representation: c = {rep.coefficients}")
    print(f"  LP representation matches:        c = {rep_lp.coefficients}")
    print(f"  3-direction LP representation:    c = {rep3.coefficients}")


if __name__ == "__main__":
    _self_check()
