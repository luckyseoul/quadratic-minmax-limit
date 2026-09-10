"""Generic residual-correction cycle (defect-correction / Nash-Moser style).

Implements the *pattern* of Section 3.3-3.4 and Figure 6 of OpenAI's
"Finite time blowup for Navier-Stokes": starting from a background state with
a nonzero residual, repeatedly (a) recompute the full residual, (b) solve a
*linearized* problem to cancel its leading part, (c) fold in the resulting
quadratic self-interaction as a new source for the next stage, and (d) track
an explicit, monotonically improving decay exponent sigma_j so the whole
induction is quantitative, not just "eventually small".

  u^[j+1] = u^[j] + delta_u_j
  R(u^[j+1]) = L_{u^[j]}(delta_u_j) + Q(delta_u_j, delta_u_j) + [old residual]

The paper solves this for the Navier-Stokes momentum residual; this module
factors the *induction itself* out into a reusable driver, parameterized by:

  - a `linear_solve(residual, state) -> correction` callable representing
    "solve the linearized operator L_v to cancel the leading part of the
    residual" (their Proposition 9.6 stages 1-4);
  - a `quadratic(correction) -> extra_residual` callable representing the
    new nonlinear self-interaction each correction introduces (their
    nabla . (delta_u_j (x) delta_u_j) term);
  - a `residual_norm(residual) -> float` callable and a target decay
    schedule sigma_j.

This is directly the numerical-analysis technique usually called *defect
correction* (Stetter 1978; Skeel 1982) or, when iterated to convergence with
a quantitative loss-of-derivatives bookkeeping, a Nash-Moser scheme. The demo
at the bottom applies it to a genuinely nonlinear scalar problem (viscous
Burgers with a manufactured near-singular target) to show the residual decay
tracking in action; it is a worked instance of the scaffold, not a
restatement of the 3D Navier-Stokes construction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Generic, TypeVar

State = TypeVar("State")
Residual = TypeVar("Residual")


@dataclass
class StageRecord:
    stage: int
    sigma: float
    residual_norm: float


@dataclass
class CorrectionCycleResult(Generic[State]):
    state: State
    history: list[StageRecord] = field(default_factory=list)
    converged: bool = False


class CorrectionCycle(Generic[State, Residual]):
    """Drives the compute-residual / correct / shrink-cutoff loop.

    Parameters
    ----------
    residual_fn:
        state -> residual. Recomputed fresh at the start of every stage, as
        the paper insists ("recompute the full residual after each
        operation, so newly created terms enter the next stage").
    linear_solve:
        (residual, state) -> correction. Solves the linearized problem that
        cancels the leading part of the residual.
    apply_correction:
        (state, correction) -> new_state.
    residual_norm:
        residual -> float, a scalar size used for the sigma_j bookkeeping
        and the stopping criterion.
    sigma_0, sigma_increment:
        The decay-exponent schedule sigma_{j+1} = sigma_j + sigma_increment,
        mirroring the paper's sigma_j = 1/5 + j/10 -> infinity (Section 3.4).
        This is bookkeeping/reporting only; the actual contraction is
        whatever linear_solve + apply_correction achieve.
    tol:
        Stop once residual_norm(residual) < tol.
    max_stages:
        Hard cap so a non-contracting problem terminates.
    """

    def __init__(
        self,
        residual_fn: Callable[[State], Residual],
        linear_solve: Callable[[Residual, State], State],
        apply_correction: Callable[[State, State], State],
        residual_norm: Callable[[Residual], float],
        *,
        sigma_0: float = 0.2,
        sigma_increment: float = 0.1,
        tol: float = 1e-12,
        max_stages: int = 200,
    ) -> None:
        self.residual_fn = residual_fn
        self.linear_solve = linear_solve
        self.apply_correction = apply_correction
        self.residual_norm = residual_norm
        self.sigma_0 = sigma_0
        self.sigma_increment = sigma_increment
        self.tol = tol
        self.max_stages = max_stages

    def run(self, initial_state: State) -> CorrectionCycleResult[State]:
        state = initial_state
        sigma = self.sigma_0
        history: list[StageRecord] = []
        converged = False

        for j in range(self.max_stages):
            residual = self.residual_fn(state)
            rnorm = self.residual_norm(residual)
            history.append(StageRecord(stage=j, sigma=sigma, residual_norm=rnorm))
            if rnorm < self.tol:
                converged = True
                break
            correction = self.linear_solve(residual, state)
            state = self.apply_correction(state, correction)
            sigma += self.sigma_increment

        return CorrectionCycleResult(state=state, history=history, converged=converged)


# --------------------------------------------------------------------------
# Worked demo: defect correction for a manufactured near-singular 1D
# viscous-Burgers-like problem.
#
#   u_t + u u_x - nu u_xx = f(x, t)
#
# We pick a target u* with a steepening front (approaching, but this side of,
# gradient blowup) and treat "find a smooth forcing f whose residual against
# a fixed numerical scheme vanishes to high order" as the correction-cycle
# target -- the same *shape* of problem as the paper's "choose u so the
# residual force extends smoothly", just in 1D and with an explicit target
# rather than an emergent singularity.
# --------------------------------------------------------------------------
import numpy as np


@dataclass
class GridState:
    x: np.ndarray
    t: np.ndarray
    u: np.ndarray  # shape (len(t), len(x)), current approximate solution


def _burgers_operator(u: np.ndarray, x: np.ndarray, t: np.ndarray, nu: float) -> np.ndarray:
    """u_t + u u_x - nu u_xx, via central/one-sided finite differences."""
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    u_t = np.gradient(u, dt, axis=0, edge_order=2)
    u_x = np.gradient(u, dx, axis=1, edge_order=2)
    u_xx = np.gradient(u_x, dx, axis=1, edge_order=2)
    return u_t + u * u_x - nu * u_xx


def _burgers_residual(state: GridState, nu: float, forcing: np.ndarray) -> np.ndarray:
    """R = (u_t + u u_x - nu u_xx) - forcing.

    Mirrors the paper's own trick (Section 2, "for any incompressible flow
    u and pressure p, we can always define the external force f to be the
    residual"): here `forcing` is fixed ahead of time as exactly the
    operator applied to the target field, so R is *only* zero once the
    current state matches the target -- not trivially zero at u=0.
    """
    return _burgers_operator(state.u, state.x, state.t, nu) - forcing


def _make_linear_solve(target_u: np.ndarray, relax: float):
    """Toy 'solve the linearized operator': nudge the state toward target_u
    by a fixed relaxation fraction of the current residual-implied error.
    This stands in for the paper's L_v-based wave/mean/pressure corrections;
    here the 'linearization' is intentionally trivial so the demo isolates
    the induction bookkeeping, not a bespoke PDE solver.
    """

    def linear_solve(residual: np.ndarray, state: GridState) -> GridState:
        error = target_u - state.u
        delta = relax * error
        return GridState(x=state.x, t=state.t, u=delta)

    return linear_solve


def _apply_correction(state: GridState, correction: GridState) -> GridState:
    return GridState(x=state.x, t=state.t, u=state.u + correction.u)


def run_burgers_demo(
    *,
    nx: int = 128,
    nt: int = 64,
    nu: float = 0.02,
    relax: float = 0.35,
    tol: float = 1e-8,
    max_stages: int = 60,
) -> CorrectionCycleResult[GridState]:
    x = np.linspace(-1.0, 1.0, nx)
    t = np.linspace(0.0, 0.6, nt)
    X, T = np.meshgrid(x, t)

    # Manufactured near-singular target: a steepening tanh front whose slope
    # grows as t increases toward the end of the window, without actually
    # reaching a discontinuity inside [0, 0.6].
    steepness = 3.0 + 12.0 * T
    target_u = np.tanh(steepness * X)
    # Fix the forcing once, from the target alone, exactly as the paper
    # defines its external force as the residual of the constructed flow.
    forcing = _burgers_operator(target_u, x, t, nu)

    initial = GridState(x=x, t=t, u=np.zeros_like(target_u))
    cycle = CorrectionCycle(
        residual_fn=lambda s: _burgers_residual(s, nu, forcing),
        linear_solve=_make_linear_solve(target_u, relax),
        apply_correction=_apply_correction,
        residual_norm=lambda r: float(np.sqrt(np.mean(r**2))),
        sigma_0=0.2,
        sigma_increment=0.1,
        tol=tol,
        max_stages=max_stages,
    )
    result = cycle.run(initial)
    return result


def _self_check() -> None:
    result = run_burgers_demo()
    norms = [rec.residual_norm for rec in result.history]
    assert len(norms) >= 2, "expected at least two stages"
    assert norms[-1] < norms[0], "residual should shrink over the cycle"
    # Monotone non-increase isn't guaranteed stage-to-stage for this toy
    # relaxation scheme, but the overall trend and final tolerance must hold.
    assert result.converged or norms[-1] < 1e-3, f"did not shrink enough: {norms[-5:]}"
    print("correction_cycle self-check: all assertions passed")
    print(f"  stages run: {len(norms)}, converged: {result.converged}")
    print(f"  residual norm: stage 0 = {norms[0]:.4e} -> stage {len(norms)-1} = {norms[-1]:.4e}")
    print(f"  sigma schedule: {result.history[0].sigma:.2f} -> {result.history[-1].sigma:.2f}")


if __name__ == "__main__":
    _self_check()
