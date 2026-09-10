"""The exact heat-exterior swirl profile from Theorem 4.6(v), eq. (4.29).

Outside the concentrating core and its annulus, the paper's flow is purely
azimuthal and solves the radial heat-type equation for swirl *exactly* --
no approximation, no correction cycle needed there. This is the one part of
the leading-order construction that is a genuine closed form rather than the
output of an existence proof, so it is the part we can faithfully implement
and independently verify, rather than merely illustrate.

Equation (4.29):

    H(Z) = 1/Gamma(1+h) * integral_0^inf e^{-v} v^h (1+Z v)^{-h} dv,
    E(X, eta) = c_inf * X^{-A} * H(2d/X),        d = 1 - h*eta^2... (uses full d)

and, translated out of similarity variables (paper, bottom of p.35):

    K(r, t) = c_inf * (r^2/2)^{-A} * H(4*tau / r^2),   tau = 1 - t,

satisfies the *exact linear PDE*

    d/dt K = (d^2/dr^2 + (1/r) d/dr - 1/r^2) K.

This is the generalized axisymmetric swirl-diffusion equation (the r^{-2}
term is the curvature/centrifugal correction; for h -> the classical value
this family reduces to the familiar Lamb-Oseen vortex). We implement H via
numerical quadrature and independently check the PDE claim by finite
differences -- i.e. we verify one exact ingredient of the paper's
construction, rather than merely restating it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as Gamma
from scipy.special import roots_genlaguerre

# Generalized Gauss-Laguerre nodes/weights for weight v^h * e^{-v}: this
# absorbs the v^h factor (the part of the integrand that is non-smooth at
# v=0 for non-integer h, which made plain Gauss-Laguerre converge only
# algebraically) directly into the quadrature weight, leaving only the
# smooth remaining factor (1+Zv)^{-h} to approximate. That gives fast
# (near machine-precision) convergence with a small, fixed node count,
# computed once and reused as a vectorized dot product instead of one
# scipy.integrate.quad call per (Z, h) pair (that scalar path is kept below
# as `H_scalar_reference` for the finite-difference PDE check, where only a
# handful of evaluations are needed).
_GENLAGUERRE_H = 0.005  # must match the h used by HeatExteriorSwirl below
_GENLAGUERRE_NODES, _GENLAGUERRE_WEIGHTS = roots_genlaguerre(40, _GENLAGUERRE_H)


def H(Z: np.ndarray | float, h: float) -> np.ndarray | float:
    """H(Z) = 1/Gamma(1+h) * int_0^inf e^{-v} v^h (1+Zv)^{-h} dv,  (eq. 4.29).

    Vectorized over Z via fixed generalized-Gauss-Laguerre quadrature with
    weight v^h e^{-v} (matching h), so this accepts and returns arrays
    without a Python-level loop. Falls back to recomputing the quadrature
    rule if called with a different h than the cached one.
    """
    Z = np.asarray(Z, dtype=float)
    if np.any(Z < 0):
        raise ValueError("Z must be >= 0 (Z = 2d/X or 4*tau/r^2, both nonnegative)")
    if h == _GENLAGUERRE_H:
        v, w = _GENLAGUERRE_NODES, _GENLAGUERRE_WEIGHTS
    else:
        v, w = roots_genlaguerre(40, h)
    # Broadcast Z (any shape) against the quadrature nodes on a new last axis.
    Zb = Z[..., np.newaxis]
    smooth_factor = (1.0 + Zb * v) ** (-h)  # v^h * e^{-v} already in w
    value = smooth_factor @ w
    result = value / Gamma(1.0 + h)
    return result if Z.shape else float(result)


def H_scalar_reference(Z: float, h: float) -> float:
    """Reference scalar implementation via scipy.integrate.quad, used only
    to cross-check the vectorized Gauss-Laguerre path in the self-check."""
    if Z < 0:
        raise ValueError("Z must be >= 0")

    def integrand(v: float) -> float:
        return np.exp(-v) * v**h * (1.0 + Z * v) ** (-h)

    value, _err = quad(integrand, 0.0, np.inf, limit=200)
    return value / Gamma(1.0 + h)


@dataclass(frozen=True)
class HeatExteriorSwirl:
    """K(r, t) = c_inf * (r^2/2)^{-A} * H(4*tau/r^2, h),  tau = 1 - t, A = 1/2 + h."""

    h: float
    c_inf: float = 1.0

    @property
    def A(self) -> float:
        return 0.5 + self.h

    def __call__(self, r: np.ndarray | float, t: np.ndarray | float) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        t = np.asarray(t, dtype=float)
        tau = 1.0 - t
        Z = 4.0 * tau / r**2
        Hvals = H(Z, self.h)
        return self.c_inf * (r**2 / 2.0) ** (-self.A) * Hvals


def _swirl_pde_residual(
    swirl: HeatExteriorSwirl, r0: float, t0: float, dr: float = 1e-3, dt: float = 1e-4
) -> tuple[float, float]:
    """Central-difference check of d/dt K == (d^2/dr^2 + (1/r) d/dr - 1/r^2) K
    at a single interior point (r0, t0). Returns (lhs, rhs)."""
    K = swirl

    dKdt = (K(r0, t0 + dt) - K(r0, t0 - dt)) / (2 * dt)

    Kp = K(r0 + dr, t0)
    Km = K(r0 - dr, t0)
    K0 = K(r0, t0)
    d2Kdr2 = (Kp - 2 * K0 + Km) / dr**2
    dKdr = (Kp - Km) / (2 * dr)
    rhs = d2Kdr2 + dKdr / r0 - K0 / r0**2

    return float(dKdt), float(rhs)


def _self_check() -> None:
    h = 0.005  # matches the paper's requirement 0 < h < 1/100
    swirl = HeatExteriorSwirl(h=h, c_inf=1.0)

    # Sanity: H(0) should be 1 (the v-integral reduces to Gamma(1+h)/Gamma(1+h)).
    assert abs(float(H(0.0, h)) - 1.0) < 1e-6, f"H(0) = {H(0.0, h)}, expected 1"

    # H should be positive and decreasing in Z (more diffusion time -> smaller
    # residual swirl at fixed X), matching "H is positive" in Theorem 4.6(v).
    Zs = [0.0, 0.5, 1.0, 2.0, 5.0, 20.0]
    Hs = [float(H(Z, h)) for Z in Zs]
    assert all(v > 0 for v in Hs), Hs
    assert all(Hs[i] >= Hs[i + 1] - 1e-9 for i in range(len(Hs) - 1)), Hs

    # Cross-check the vectorized Gauss-Laguerre path against the scalar
    # scipy.integrate.quad reference implementation.
    for Z in Zs:
        ref = H_scalar_reference(Z, h)
        fast = float(H(Z, h))
        assert abs(ref - fast) < 5e-6, f"Z={Z}: quad={ref}, gauss-laguerre={fast}"

    # Vectorized call over an array must match elementwise scalar calls.
    Z_arr = np.array(Zs)
    H_vec = H(Z_arr, h)
    assert H_vec.shape == Z_arr.shape
    assert np.allclose(H_vec, Hs, atol=1e-9)

    # PDE check at several interior (r, t) points, away from r=0 and t=1.
    max_rel_err = 0.0
    for r0, t0 in [(1.0, 0.0), (1.5, 0.3), (0.7, -0.5), (2.0, 0.6)]:
        lhs, rhs = _swirl_pde_residual(swirl, r0, t0)
        rel_err = abs(lhs - rhs) / max(abs(rhs), 1e-8)
        max_rel_err = max(max_rel_err, rel_err)
    assert max_rel_err < 5e-2, f"PDE residual too large: {max_rel_err:.4f}"

    print("heat_exterior self-check: all assertions passed")
    print(f"  h = {h}, H(0) = {H(0.0, h):.6f} (expect 1.0)")
    print(f"  H(Z) is positive and decreasing on {Zs}: {[round(v,4) for v in Hs]}")
    print(f"  max relative PDE residual over test points: {max_rel_err:.4e}")


if __name__ == "__main__":
    _self_check()
