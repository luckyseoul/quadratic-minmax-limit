# Actual cross sign matrices: singular-moment rounding constraints

2026-09-05. Analytic draft. This supplies moment constraints for an
actual cross matrix and its shell-width resolvents. It is not an
evaluation of the final conditional Gaussian maximum.

## 1. Exact finite-order inequality

Let B be an n by n sign matrix, let sigma_1,...,sigma_n be its singular
values, and put kappa=2/pi. Write

\[
 \beta(B)=\max_{x,y\in\{-1,1\}^n}|x^TBy|.
\]

Let f be any real function on the finite singular spectrum satisfying
`0<=f(sigma)<=sigma`. Functional calculus is applied to
`|B|=(B^T B)^(1/2)`. Then

\[
 \boxed{\displaystyle
 \beta(B)\ge {\kappa\over n}\sum_j\sigma_j^2 f(\sigma_j)
       -{1-\kappa\over n^2}\sum_j\sigma_j^2 f(\sigma_j)^2.}
                                                               \tag{1}
\]

To prove this, let u_i be row i of B divided by sqrt(n), and let v_j
be column j of f(|B|) divided by sqrt(n). Each u_i is a unit vector.
Moreover `f(|B|)^2 <= |B|^2=B^T B` in the PSD order, so each v_j
has norm at most one. Add mutually orthogonal auxiliary coordinates
to the v_j, orthogonal to every u_i, to make every norm exactly one.
The cross Gram matrix remains

\[
                   Z={1\over n}Bf(|B|),\qquad |Z_{ij}|\le1.
                                                               \tag{2}
\]

Using one standard Gaussian vector on the resulting common Euclidean
space, set `x_i=sign(<u_i,g>)`, `y_j=sign(<v_j,g>)`. The Gaussian
arcsine identity gives

\[
 \beta(B)\ge\mathbb E[x^TBy]
                    =\kappa\sum_{ij}B_{ij}\arcsin Z_{ij}.
                                                               \tag{3}
\]

The positive Taylor coefficients of arcsin imply, for every z in
[-1,1],

\[
 |\arcsin z-z|\le(\pi/2-1)|z|^2.                             \tag{4}
\]

Indeed the powers in the remainder are at least three and hence at
most |z|^2 on this interval; the sum of their coefficients is
pi/2-1. The endpoint cases follow by continuity. Thus (3) is at least

\[
 \kappa\langle B,Z\rangle-(1-\kappa)\|Z\|_F^2.
\]

Cyclic trace identities give

\[
 \langle B,Z\rangle={1\over n}\sum_j\sigma_j^2f(\sigma_j),
 \qquad
 \|Z\|_F^2={1\over n^2}\sum_j\sigma_j^2f(\sigma_j)^2,
\]

proving (1). No SDP optimality, operator cap, randomness of B, or
conditional optimality is assumed.

## 2. Cubic and clipped consequences

Taking f(sigma)=sigma gives the exact inequality

\[
 \boxed{\displaystyle
 \beta(B)\ge {\kappa\over n}\sum_j\sigma_j^3
       -{1-\kappa\over n^2}\sum_j\sigma_j^4.}                \tag{5}
\]

For any fixed K>0, the choice f(sigma)=min(sigma,K sqrt(n)) gives

\[
 \boxed{\displaystyle
 \beta(B)\ge {\kappa\over n}
       \sum_j\sigma_j^2\min(\sigma_j,K\sqrt n)
                    -(1-\kappa)K^2n.}                      \tag{6}
\]

The error follows from `sum sigma_j^2=n^2`. Formula (6) applies to
every B, including matrices with operator norm much larger than
sqrt(n). It does not replace B or alter its entries.

## 3. The actual normalized resolvent constraint

Let d>=||B||op, write d=q sqrt(n), and define the empirical moments

\[
  y_j=\sigma_j^2/d^2\in[0,1],\qquad
  \langle h(y)\rangle={1\over n}\sum_j h(y_j),\qquad
  b={\beta(B)\over n^{3/2}}.
\]

The sign-entry Frobenius identity gives `q>=1` and

\[
             \langle y\rangle={1\over q^2}.
                                                               \tag{7}
\]

Equation (5) gives, with no asymptotic approximation,

\[
 \boxed{\displaystyle
 \langle y^{3/2}\rangle
 \le {b\over\kappa q^3}
       +{1-\kappa\over\kappa}\,{q\over\sqrt n}
                                             \langle y^2\rangle.}
                                                               \tag{8}
\]

For q bounded by a fixed constant, its error is O(n^(-1/2)).
These are constraints on the actual singular spectrum, not an
arbitrary probability measure introduced without source conditions.

If an actual shell has cross energy c=beta(B), and
`u=c/(nd)=b/q`, then (7)--(8) give

\[
 \langle y^{3/2}\rangle
       \le {u\over\kappa}\langle y\rangle+O(n^{-1/2}).       \tag{9}
\]

The same statement holds after replacing a negatively attaining
cross energy by its absolute value. If |c|<beta(B), the correct
numerator in (8) is still beta(B), not |c|. Conditional optimality
only gives beta(B)<=F_A^*; it does not authorize this substitution
on every shell.

In particular, if u<kappa in (9), the whole Frobenius-weighted mass
cannot concentrate at y=1. This is the concrete extra constraint
available when evaluating the Boolean-sensitive cross resolvent.

## 4. Exact partial-isometry check and remaining scope

If every nonzero singular value equals d, equation (5) reduces to

\[
              \beta(B)\ge\kappa nd-(1-\kappa)d^2.             \tag{10}
\]

Here the natural cross Gram is `(d/n)B`. It attains vector objective
nd, and the operator bound gives the matching SDP upper nd.
Consequently, when d=O(sqrt(n)), this exact partial-isometry case
has asymptotic SDP/integer ratio at most pi/2. This observation does
not extend that constant to arbitrary sign matrices or to every
actual conditional optimizer.

Equations (7)--(9) may be combined with the previously reviewed full
Boolean ellipsoid remainder. They do not by themselves evaluate its
supremum over all actual spectra and attainable joint shells, and
they do not establish the original doubling target or convergence.

No numerical experiment or finite-order spectral search is used.
