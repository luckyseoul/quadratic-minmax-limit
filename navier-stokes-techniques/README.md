# Techniques from OpenAI's Navier–Stokes blowup proof

Implementations of techniques borrowed from OpenAI's "Finite time blowup for
Navier–Stokes" (2026) and its companion Lean 4 formalization.

Source material (downloaded locally, not vendored here):

* paper: `/home/nick/openai-navier-stokes/navier-stokes.pdf` (166 pp.)
  — <https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf>
* Lean: `/home/nick/NavierStokesAndEuler` — <https://github.com/openai/NavierStokesAndEuler>

## What the paper does, in one paragraph

For every viscosity `nu > 0` it builds a smooth, compactly supported force and
a Navier–Stokes solution starting from **rest** whose velocity becomes
unbounded in finite time while kinetic energy stays bounded — establishing
Clay's alternative (C), and (D) on the torus. The construction is a
self-similar, contracting vortex **core** (radial and axial scales shrinking
at *different* rates, `ell_r ~ tau^{1/2}` vs `ell_z ~ tau^{1/2-h}`), joined
across an **annulus** to an exactly-solvable **heat exterior**. The join leaves
a momentum residual; rather than forcing it away, two families of oscillatory
**pulses**, amplified by the background shear, supply it through their
*nonlinear* averaged momentum flux, as a strictly positive combination of two
covariance directions. A Nash-Moser-style **correction cycle** then drives the
leftover residual to be flat at the singular time.

## Contents

`quadratic_repair.py` applies the cone and correction-cycle primitives to
actual complete signings: fractional cone coefficients rank integral edge
repairs, and exact full Boolean norms decide every accepted correction.
The [mesh experiment](../evidence/fixed_repair_path_20260911/README.md)
records two escaped single-edge stalls among ten prescribed inputs, with
matched proposal controls and independent integer replay. This is finite
evidence of local usefulness, with no all-orders convergence claim.

| Directory | What it is | Fidelity |
|---|---|---|
| `stress_cone/` | Positive-cone representation `T = c1 v1 + c2 v2`, `c1,c2 > 0` (§3.2, Prop 7.5, Fig 4) | Faithful algorithm, generic inputs |
| `correction_cycle/` | Reusable residual-correction / defect-correction driver with `sigma_j` bookkeeping (§3.4, Fig 6) | Faithful *pattern*, toy 1D demo |
| `blowup_demo/` | `heat_exterior.py`: the **exact** exterior swirl profile (eq. 4.29) + PDE verification; `blowup_demo.py`: reduced-order core | Exterior exact; core illustrative |
| `lean_verification/` | Building and independently re-checking the Lean proofs (Comparator + nanoda dual-kernel) | Reproduction of their own pipeline |

Run any module directly to execute its self-check:

```sh
python3 stress_cone/stress_cone.py
python3 correction_cycle/correction_cycle.py
cd blowup_demo && python3 heat_exterior.py && python3 blowup_demo.py
```

## Honest scoping

Two things here are genuinely faithful to the paper: the **heat exterior**
(eq. 4.29 is a real closed form, and `heat_exterior.py` independently verifies
by finite differences that it solves its claimed swirl-diffusion PDE
`d_t K = (d_rr + r^-1 d_r - r^-2) K`), and the **Lean verification**, which is
their actual artifact re-checked with an independent kernel.

Everything else implements a *technique*, not the paper's specific objects:

* The **stress cone** solver takes `v1, v2, T` as inputs. The paper's actual
  covariance directions come out of the full profile construction
  (§4–§7); deriving them is not reproduced.
* The **correction cycle** captures the induction's shape (recompute residual
  → linear solve → fold in the quadratic self-interaction → shrink cutoff,
  tracking a monotone `sigma_j`). The demo applies it to 1D viscous Burgers
  with a manufactured forcing, not to the 3D momentum residual.
* The **core** in `blowup_demo.py` is a smooth bump carrying the paper's
  *reported scaling exponents*, not a solution of the profile equations.
  Constructing the real profiles is Theorem 4.6 — roughly 20 pages plus
  Appendices A–C, including a finite-dimensional moment solve joined to an
  outer branch under the cone constraint. That is not redone here.

The reduced model does reproduce both reported exponents numerically:
fitted velocity growth `tau^{-0.5050}` against the paper's `-(1/2+h)`, and
fitted core-energy decay `tau^{0.4850}` against the paper's `1/2-3h`, at
`h = 0.005`.

## A numerical-methods note

`heat_exterior.py` evaluates `H(Z) = 1/Gamma(1+h) * int_0^inf e^-v v^h (1+Zv)^-h dv`.
Plain Gauss–Laguerre converges only *algebraically* here, because `v^h` is not
smooth at `v=0` for non-integer `h` (64 nodes → ~5e-5 error; 160 nodes still
~2e-5). Switching to **generalized** Gauss–Laguerre with weight `v^h e^{-v}`
absorbs the singular factor into the quadrature weight, leaving only the
smooth `(1+Zv)^{-h}`: 40 nodes then give ~1e-14 at small `Z`. That change also
took the demo from >5 minutes to ~5 seconds.
