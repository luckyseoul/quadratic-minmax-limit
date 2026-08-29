# Exact Max+ dilation energy and the character/PSD no-go

**Date:** 2026-08-29
**Status:** exact normalization proved; proposed character-only and PSD-only
proof routes are insufficient; R1 remains open

Put `q=p^2`, let `H=(F_q^*)^2`, `K=H/{+1,-1}`, and write
`q-1=8R`. If `lambda_e,lambda_m` are the exceptional and principal scalars
of the full-Max+ frame operator and

```text
lambda_bar = 8(q-1)/(q-5),
mu_e = lambda_e-lambda_bar,
mu_m = lambda_m-lambda_bar,
```

then the trace constraint is

```text
mu_e/2 + sum_m mu_m = 0.
```

The centered square-torus character is

```text
gamma(l) = mu_e(-1)^l
         + sum_m mu_m(exp(pi*i*m*l/R)+exp(-pi*i*m*l/R)).
```

Thus `gamma(1)=gamma(-1)=0`, and cyclic orthogonality gives exactly

```text
(1/(q-1)) sum_{a in H} |gamma(a)|^2
  = V/n
  = 24||delta||^2/n.
```

For a Max+ vector, let `D` be its negative finite coordinates,
`N(u)=|D intersect (D-u)|`, and

```text
L(a) = E_y sum_{u != 0} N_y(u)N_y(au).
```

Affine-map averaging gives

```text
Gamma(a) = (16/q)L(a) - (q-1)(q+5),
L_bar = q(q-1)(q+5)/16.
```

If `S_K` is the left side of equation (5) in the cold strategy note, then

```text
S_K = (q-1)V/(2n) = 12(q-1)||delta||^2/n.
```

Therefore

```text
S_K <= q-1  iff  V <= 2n  iff  ||delta||^2 <= n/12.
```

Equation (5) is exactly strong R1 in dilation coordinates, not a separate
lemma that representation theory can prove downstream.

There are two rigorous route obstructions.

1. Positivity, equivariance, trace, constituent multiplicities, and fixed
   vector norm permit abstract invariant ensembles whose spectrum violates
   the bound by a factor of order one or larger. What they lack is the
   Boolean rank-one identity `B_y+2P_+=yy^T`.
2. On `K`, with `m=|K|` and `c=lambda_bar`, additive Fourier positivity
   rewrites equation (5) as

   ```text
   sum_{chi != 1}(lambda_chi-c)^2 <= 4.
   ```

   A quartic-residue PSD autocorrelation has the same nonnegativity,
   diagonal, mass, and trace constraints but violates this inequality by a
   factor `Theta(q^2)`. At `p=7` an integral instance has energy
   `1105920/11` against target `48`.

These are not counterexamples to the actual uniform full-Max+ ensemble.
They prove that torus character orthogonality, Bochner positivity, and
autocorrelation constraints alone cannot close R1. A successful proof must
use phase-sensitive Boolean identities and exact cancellation between all
Max+ orbit types.
