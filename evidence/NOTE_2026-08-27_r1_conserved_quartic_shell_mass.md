# Conserved positive quartic mass on R1 dual shells

Date: 2026-08-27. This is Proposition 15.665. It supplies the nonlinear
shell coupling absent from Proposition 15.641's linear modular constraints.
It does not, by itself, prove R1.

## Setup

Let

\[
 n=p^2+1,\qquad d=n/2,\qquad P=(I+C/p)/2,
\]

and regard `P` as the orthogonal projector onto the `+p` conference
eigenspace. Put

\[
 \mathcal Z=\{W=W^T:PWP=W,\ \operatorname{diag}W=0\},
 \qquad z=\dim\mathcal Z={n(n-6)\over8}.
\]

For `x` in the `+p` eigenspace, let

\[
 b_x=\Pi_{\mathcal Z}(xx^T).
\]

Then `x^T W x=<b_x,W>` for every `W` in `Z`. For a complete dual shell
`X_s` of size `N_s` and common squared radius `r_s`, define

\[
 R_s=\sum_{x\in X_s} b_x\otimes b_x.
\]

This operator is positive semidefinite before any theta phase is inserted.

## Harmonic operator equals raw mass minus one scalar

For the degree-four harmonic polynomial used in Proposition 15.631,

\[
 H_W(x)=(x^TWx)^2-{4\over d+4}\|x\|^2x^TW^2x
 +{2\|W\|_F^2\over(d+2)(d+4)}\|x\|^4,
\]

the shell is a tight frame in the irreducible `+p` space:

\[
 \sum_{x\in X_s}xx^T={N_s r_s\over d}P.
\]

Consequently

\[
 \sum_{x\in X_s}\|x\|^2x^TW^2x={N_s r_s^2\over d}\|W\|_F^2.
\]

Combining the last two radial terms gives exactly

\[
 A_s=R_s-\rho_s I_{\mathcal Z},\qquad
 \boxed{\rho_s={2N_s r_s^2\over d(d+2)}}.                 \tag{1}
\]

Thus the previously signed harmonic eigenvalue `a_(s,c)` becomes the
nonnegative raw eigenvalue

\[
 q_{s,c}=a_{s,c}+\rho_s.                                  \tag{2}
\]

## Closed scalar trace series

The diagonal map on `Sym^2(P R^n)` has Gram matrix

\[
 K=P\circ P={(p^2-1)I+J\over4p^2},\qquad
 K^{-1}={4p^2\over p^2-1}I-{2\over p^2-1}J.               \tag{3}
\]

Writing `a_i=x_i^2` and `r=||x||^2`, orthogonal projection onto the kernel
of that diagonal map yields

\[
 \boxed{\|b_x\|_F^2
 =r^2-{4p^2\over p^2-1}\sum_i x_i^4+{2\over p^2-1}r^2.}  \tag{4}
\]

If `{W_alpha}` is an orthonormal basis of `Z`, equivariance and a trace
calculation give

\[
 \sum_\alpha W_\alpha^2={z\over d}P.                      \tag{5}
\]

Therefore the single polynomial whose shell coefficient is `tr(A_s)` is

\[
 H_{\rm tr}(x)=\|b_x\|_F^2
 -{4z\over d(d+4)}r^2+{2z\over(d+2)(d+4)}r^2.             \tag{6}
\]

No basis of the `z`-dimensional space is needed to construct this weighted
theta series. Its shell coefficient `h_s` determines the complete raw mass:

\[
 \boxed{\tau_s=\operatorname{tr}R_s=h_s+z\rho_s.}         \tag{7}
\]

There is also a computationally useful symmetry reduction. For the
coordinate functional `ell_i(x)=x_i` on the `P`-space,
`||ell_i||^2=P_ii=1/2`, so

\[
 Z_i(x)=x_i^4-{3\over d+4}r x_i^2
        +{3\over4(d+2)(d+4)}r^2
\]

is harmonic. Coordinate transitivity and `sum_i x_i^2=r` simplify (6) to

\[
 \Theta_{H_{\rm tr}}
 =-{4p^2(p^2+1)\over p^2-1}\Theta_{Z_i}.           \tag{7a}
\]

Thus the exact trace series can be computed from one compact zonal quartic
instead of expanding all `p^2+1` coordinate fourth powers.

## Conserved channel mass

Every complete shell is PSL-invariant. The real PSL decomposition of
`Z` used in Propositions 15.597 and 15.640 is multiplicity-free, so `R_s`
is scalar on each constituent `c` of dimension `m_c`. Positivity and (7)
give

\[
 \boxed{q_{s,c}\ge0,\qquad
        \sum_c m_cq_{s,c}=\tau_s,\qquad
        q_{s,c}\le{\tau_s\over m_c}.}                     \tag{8}
\]

The smallest constituent dimension is `d=(p^2+1)/2`; every ordinary
principal constituent has dimension `n=p^2+1`. Thus (8) improves the
generic rank-one shell upper bound by a factor of order `p^2` and couples
all channels through one conserved scalar coefficient.

For reference, the constituent partition has one Weil constituent and
`(p^2-9)/8` principal constituents. Their square-circle grouping is:

- if `p=1 (mod 4)`: kernel `(p-1)(p-3)/8` principal; low `(p-1)/4`
  principal; high one Weil plus `(p-5)/4` principal;
- if `p=3 (mod 4)`: kernel `(p-1)(p-3)/8` principal; low one Weil plus
  `(p-3)/4` principal; high `(p-3)/4` principal.

Their weighted dimensions sum to `z` in both cases.

## Exact p=11 checks

At `p=11`, the channel dimensions are `1220`, `305`, and `244`. For the
four proved nonempty shells, equations (1), (2), and (7) give:

| scaled norm | `N_s` | `tr(A_s)` | `rho_s` | `tau_s` |
|---:|---:|---:|---:|---:|
| 11 | 244 | `-3538/63` | `2/63` | `0` |
| 20 | 16,104 | `-5368/21` | `1600/231` | `923784/77` |
| 24 | 14,762 | `-15921/28` | `64/7` | `436943/28` |
| 27 | 442,860 | `-527406/7` | `2430/7` | `538752` |

The channel-weighted raw eigenvalues reproduce each `tau_s` exactly and
are all nonnegative. The complete rational audit is emitted by
`src/e1_gmin_m4_prop15665.py`.

## Computational consequence and remaining scope

For `p=11`, every constituent's harmonic theta series lies in the exact
66-dimensional Kohnen space already cached for Proposition 15.641. The
cusp gaps and all coefficients through exponent 28 leave a 32-dimensional
affine space per square-circle channel. Applying (8) at every later
coefficient produces a rational LP for the half-cusp target. This is a
strict strengthening of the falsified coefficient-determination route:
the old nullspace witness need not satisfy raw-shell positivity.

The LP is still an outer relaxation of actual channel theta series, and a
`p=11` certificate is not an all-prime theorem. Closing R1 requires either
a uniform version of the resulting dual certificate or a multi-scale theta
inequality that transports (8) to the odd-coset target. No R1, global-QVAR,
Type-I, or final-limit flag is flipped by Proposition 15.665 alone.

The mechanism is compatible with the weighted harmonic-theta and lattice
shadow framework already discussed in `HISTORY_AND_REFERENCES.md`, especially
Rains--Sloane. Their strongly modular conclusions still do not apply
verbatim to this growing-level `4p` family; equations (3)--(8) are the
Paley-specific input.
