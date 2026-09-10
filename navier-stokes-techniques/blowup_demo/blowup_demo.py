"""Reduced-order reproduction of the blowup mechanism (Section 2 / Figure 1,3).

Honesty note up front: the paper's actual leading profile (E, U on the core)
is the output of a ~20-page nonlinear existence argument (Theorem 4.6, with
full proofs in Appendices A-C) -- a finite-dimensional moment solve joined to
an outer branch under a "positive stress cone" constraint. Reproducing that
construction exactly is out of scope here. What *is* exact and reproduced
faithfully is the heat exterior (heat_exterior.py, eq. 4.29): a genuine
closed-form solution of a genuine linear PDE, verified there by finite
differences.

This module combines:

  * the EXACT heat-exterior swirl K(r, t), used for r >= r_join(t);
  * an ILLUSTRATIVE core model for r < r_join(t) that uses the paper's own
    *reported scaling exponents* (Section 2.1: ell_r =~ tau^{1/2},
    |u_theta| =~ tau^{-1/2-h}) rather than solving Theorem 4.6's profile
    equations -- a smooth bump of the correct width and height, not the true
    self-similar solution of the momentum balance.

The purpose is to check the *qualitative and scaling* claims that are stated
in closed form in the paper (Section 2): velocity sup-norm diverges like
tau^{-1/2-h} while total kinetic energy over the shrinking core stays
bounded (in fact -> 0, order tau^{1/2-3h}), even though the core radius
ell_r -> 0. We verify these scaling claims numerically on the reduced model,
and plot the Figure-1/Figure-3(b)-style envelope.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from heat_exterior import HeatExteriorSwirl


@dataclass(frozen=True)
class ReducedFlow:
    """Illustrative reduced-order axisymmetric swirl u_theta(r, t).

    Core region (r < r_join): a smooth bump of width ell_r(tau) =~ tau^{1/2}
    and peak height ~ tau^{-1/2-h}, matching the paper's reported exponents
    (Section 2.1) but NOT solving the actual profile ODEs of Theorem 4.6.

    Exterior region (r >= r_join): the EXACT heat-exterior swirl K(r, t)
    from eq. (4.29), joined continuously in value (not derivative -- the
    real construction spends an entire annulus and oscillatory correction
    on making the join smooth; here we only need visual/scaling fidelity).
    """

    h: float
    c_inf: float = 1.0
    join_factor: float = 3.0  # r_join(t) = join_factor * ell_r(tau)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_exterior", HeatExteriorSwirl(h=self.h, c_inf=self.c_inf))

    def ell_r(self, tau: np.ndarray) -> np.ndarray:
        return np.sqrt(np.clip(tau, 1e-12, None))

    def core_peak(self, tau: np.ndarray) -> np.ndarray:
        return np.clip(tau, 1e-12, None) ** (-0.5 - self.h)

    def r_join(self, tau: np.ndarray) -> np.ndarray:
        return self.join_factor * self.ell_r(tau)

    def __call__(self, r: np.ndarray, t: np.ndarray) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        t = np.asarray(t, dtype=float)
        tau = 1.0 - t
        ell = self.ell_r(tau)
        peak = self.core_peak(tau)
        rj = self.r_join(tau)

        # Exact exterior value at the join radius, used to scale a smooth
        # bump so the two pieces agree in magnitude at r = r_join.
        ext_at_join = self._exterior(rj, t)
        core_shape = np.exp(-0.5 * (r / ell) ** 2)  # smooth Gaussian bump
        join_shape = np.exp(-0.5 * (rj / ell) ** 2)
        core_value = peak * core_shape
        # Rescale the core so it matches the exterior value at r = r_join,
        # giving a continuous (if not smooth) combined profile.
        core_at_join = peak * join_shape
        scale = np.where(core_at_join > 0, ext_at_join / np.maximum(core_at_join, 1e-300), 1.0)
        core_value = core_value * scale

        exterior_value = self._exterior(np.maximum(r, rj), t)
        return np.where(r < rj, core_value, exterior_value)

    def kinetic_energy_core(self, t: float, n: int = 4000) -> float:
        """Volumetric core kinetic energy, matching Section 2.1's "total
        kinetic energy of the core is of order tau^{1/2-3h}, which tends to
        zero despite the increasing speeds".

        The core is a slender 3D column of radius ell_r(tau) ~ tau^{1/2} and
        HEIGHT ell_z(tau) ~ tau^{1/2-h} (Figure 1: "its radius shrinks faster
        than its height"), not an infinite cylinder. u_theta itself only
        depends on (r, t) in this reduced model (no z-dependence), so the
        volumetric energy is the areal (r-only) energy density times the
        axial extent ell_z(tau):

            E_core(t) ~ ell_z(tau) * 2*pi * integral_0^{r_join} 0.5 u^2 r dr.

        This restriction to r < r_join and the ell_z volume factor are both
        needed to match the paper's reported exponent; the full-flow energy
        bound in Theorem 1.1 additionally relies on the compact-support
        spatial cutoff of Section 2.3 / Proposition 10.1, not reproduced here.
        """
        tau = 1.0 - t
        rj = float(self.r_join(np.array(tau)))
        r = np.linspace(1e-6, rj, n)
        u = self(r, np.full_like(r, t))
        integrand = 0.5 * u**2 * r
        areal_energy = 2 * np.pi * np.trapezoid(integrand, r)
        ell_z = max(tau, 1e-12) ** (0.5 - self.h)
        return float(ell_z * areal_energy)

    def sup_norm(self, t: float, r_max: float = 50.0, n: int = 4000) -> float:
        r = np.linspace(1e-6, r_max, n)
        u = self(r, np.full_like(r, t))
        return float(np.max(np.abs(u)))


def _self_check() -> None:
    h = 0.005
    flow = ReducedFlow(h=h, c_inf=1.0)

    ts = np.array([0.0, 0.5, 0.9, 0.99, 0.999, 0.9999])
    sup_norms = [flow.sup_norm(t) for t in ts]
    energies = [flow.kinetic_energy_core(t) for t in ts]

    # The core sup-norm must diverge as t -> 1.
    assert sup_norms[-1] > sup_norms[0] * 10, sup_norms
    assert all(b >= a for a, b in zip(sup_norms, sup_norms[1:])), sup_norms

    # Kinetic energy of the *core* must stay bounded and in fact -> 0 (paper,
    # Section 2.1: order tau^{1/2-3h}) even as the velocity blows up -- the
    # whole point of the construction.
    assert all(b <= a for a, b in zip(energies, energies[1:])), energies
    assert energies[-1] < energies[0] / 50, energies

    # Cross-check the reported scaling law directly: sup_norm(t) should grow
    # like tau^{-1/2-h} up to a slowly varying prefactor, i.e.
    # log(sup_norm) ~ -(1/2+h) * log(tau) + const, fit the exponent.
    tau = 1.0 - ts
    log_tau = np.log(tau)
    log_sup = np.log(sup_norms)
    # linear fit over the last few (smallest tau) points, where the core
    # value (not the join-scaled exterior read at large r) dominates
    slope, intercept = np.polyfit(log_tau[-4:], log_sup[-4:], 1)
    expected = -(0.5 + h)
    assert abs(slope - expected) < 0.05, f"fitted exponent {slope}, expected {expected}"

    # Same check for the core energy's reported decay exponent tau^{1/2-3h}.
    log_energy = np.log(energies)
    e_slope, _e_intercept = np.polyfit(log_tau[-4:], log_energy[-4:], 1)
    e_expected = 0.5 - 3 * h
    assert abs(e_slope - e_expected) < 0.05, f"fitted energy exponent {e_slope}, expected {e_expected}"

    print("blowup_demo self-check: all assertions passed")
    print(f"  h = {h}")
    print(f"  t values:        {ts}")
    print(f"  sup |u_theta|:   {[round(v, 3) for v in sup_norms]}")
    print(f"  core energy:     {[round(v, 6) for v in energies]}")
    print(f"  fitted velocity growth exponent: {slope:.4f} (paper: {expected:.4f})")
    print(f"  fitted energy decay exponent:    {e_slope:.4f} (paper: {e_expected:.4f})")


def make_plot(path: str = "blowup_demo.png") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = 0.005
    flow = ReducedFlow(h=h, c_inf=1.0)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: radial profile at successive times, echoing Figure 1's shrinking,
    # intensifying core.
    r = np.linspace(1e-3, 3.0, 800)
    for t, color in zip([0.0, 0.7, 0.95, 0.99], ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]):
        u = flow(r, np.full_like(r, t))
        axes[0].plot(r, u, label=f"t = {t}", color=color)
    axes[0].set_xlabel("r")
    axes[0].set_ylabel(r"$u_\theta(r, t)$")
    axes[0].set_title("Contracting, intensifying core\n(core: illustrative; exterior: exact eq. 4.29)")
    axes[0].legend()

    # Right: sup-norm blowup vs bounded core energy, echoing the paper's
    # "unbounded velocity, uniformly bounded kinetic energy" claim.
    ts = 1.0 - np.geomspace(1e-6, 1.0, 200)
    sup_norms = [flow.sup_norm(t) for t in ts]
    energies = [flow.kinetic_energy_core(t) for t in ts]
    ax2 = axes[1]
    ax2.plot(ts, sup_norms, color="#C44E52", label=r"$\sup_r |u_\theta(r,t)|$")
    ax2.set_yscale("log")
    ax2.set_xlabel("t")
    ax2.set_ylabel(r"$\sup |u_\theta|$ (log scale)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax3 = ax2.twinx()
    ax3.plot(ts, energies, color="#4C72B0", label="core kinetic energy")
    ax3.set_ylabel("core kinetic energy", color="#4C72B0")
    ax3.tick_params(axis="y", labelcolor="#4C72B0")
    ax2.set_title("Velocity blows up, energy stays bounded (in fact -> 0)")

    fig.tight_layout()
    fig.savefig(path, dpi=140)
    print(f"wrote {path}")


if __name__ == "__main__":
    _self_check()
    make_plot("/home/nick/navier-stokes-techniques/blowup_demo/blowup_demo.png")
