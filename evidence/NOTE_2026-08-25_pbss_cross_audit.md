# Cross-audit: Perry--Beurling Spectral Sieve versus the two live leftovers

Date: 2026-08-25.  Read-only audit of
`/home/nick/perry-beurling-spectral-sieve`.  The PBSS worktree was already
dirty and was not modified.

## Verdict

There is no drop-in PBSS theorem or implementation that closes either the
non-Walsh affine residual or R1.  One abstraction is worth transplanting to
the R1 attack: explicitly separate a known positive channel from an
adversarial cancellation operator, then demand a quantitative coercivity
bound before treating the positive channel as evidence of the whole sign.

## The useful R1 translation

Proposition 15.631 gives the exact positive first-dual-shell contribution

\[
 Q_{0,p}(W;t)=\frac{\|W\|_F^2}{8(d+2)}e^{-\pi/(8t)}.
\]

Write the higher dual shells as a quadratic operator on admissible `W`:

\[
 R_p(W;t)=\sum_{\substack{u\in L^*\\
                 \|u\|^2\ge(p-1)/p}}
 (-1)^{2p\|u\|^2}H_W(u/2)e^{-\pi\|u\|^2/(4t)}.
\]

The PBSS `HSTAR_ANTICANCELLATION` pattern says to test the actual hinge

\[
 \|R_p(t)\|_{\mathrm{op}}
 <\frac{1}{8(d+2)}e^{-\pi/(8t)},                           \tag{C}
\]

not merely to observe that the first term is positive.  Its forced
projection-cancellation control is the right sanity test: if the available
higher-shell model can cancel `Q_0`, first-shell positivity has bought no
global sign.

For R1, (C) at one Gaussian parameter is still insufficient.  The target is
the first *primal* odd-coset coefficient, while small `t` isolates the first
dual shell and large `t` isolates the first primal shell.  The actionable
variant is therefore a finite signed combination of Gaussian parameters:

1. choose coefficients `c_j,t_j` so the primal radial kernel
   `sum c_j exp(-pi t_j r)` isolates `r=n` and suppresses `r>n`;
2. use Poisson summation to obtain the reciprocal dual kernel;
3. minimize the higher-dual-shell operator norm subject to a fixed positive
   first-shell coefficient;
4. run an explicit forced-cancellation negative control;
5. promote only a uniform analytic inequality, or record a general
   countermechanism if the LP remains cancellable.

This is a multi-Gaussian theta-window / Cohn--Elkies-style linear program,
not code already present in PBSS.  PBSS contributes the experiment design,
not the mathematics of the Paley lattice.

## Non-Walsh residual translation

PBSS's Legendre projection code is continuous, floating point, and tied to a
logarithmic interval.  The corresponding exact basis for Proposition 15.632
would be the degree-0/1/2 Hahn basis under the hypergeometric law on a Johnson
slice.  The current exact three-variable parity-majorant LP already uses all
pointwise degree-two constraints and is stronger than simply computing that
projection energy.  Copying PBSS's NumPy/CuPy projector would therefore add
machinery without strengthening the bound.

## Components not to import as proofs

- The repository's own project record says the original spectral-sieve
  operator collapsed to a weighted `L2` norm.  Its classification code does
  not address either live inequality here.
- `HSTAR_ANTICANCELLATION` is explicitly a model statement with an injected
  off-critical bump.  The arithmetic anti-cancellation step remains open; it
  is an analogy, not a theorem transferable to R1.
- The weighted-mode Lemma M6 is not valid under its stated `L-infinity`-only
  hypothesis.  The displayed step

  \[
  |\langle wq,\phi\rangle|
  \le\|w\|_\infty|\langle q,\phi\rangle|
  \]

  is false for signed integrals.  A valid integration-by-parts estimate needs
  regularity/variation control of `w*phi`, with the derivative term included.
  No M6 constant should be imported.
- The generic GPU projector performs floating-point dense products.  It does
  not preserve the exact parity, integrality, or finite-field structure needed
  by Proposition 15.632.
- The Jensen/moment note is a useful warning that low-order bulk diagnostics
  can miss a remote binding feature, but its parameter-count slogan is not a
  substitute for a rigorous theta-shell tail estimate.

## Concrete next use

If the current exact jobs do not close their branches, prototype the
multi-Gaussian theta-window locally in this repository.  Start with exact or
interval-certified low dual shells and operator norms at the smallest live
prime, then ask whether the optimized margin scales in a form provable from
the all-prime dual norm gap.  Do not add PBSS as a dependency.
