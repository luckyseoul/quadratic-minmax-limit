# A cap-independent paired-field bound for the original signing minimum

2026-09-17. Analytic theorem; convergence remains open.

Write `kappa=2/pi` and

    B = [81 kappa/2 + 9 sqrt(kappa(2/3+kappa/2))]/101.

For the original complete symmetric zero-diagonal signing minimum,

    liminf_n m_n/n^(3/2) >= B > 13/40.                     (1)

The bound is approximately 0.32584. It is not a proof of convergence,
an identification of the limit, or a finite-order cutoff. The new step
is a paired local-field estimate whose limiting constant is independent
of a fixed source operator cap. This permits the existing same-order
regularization theorem to remove the cap without losing the improvement.

## 1. Reused inputs and the new implication

Put `M=A/sqrt(n)`, `q=(n-1)/n`, `Q_M(x)=x^T M x/2`, and
`alpha=Phi(A)/n^(3/2)`. Completeness gives

    diag M=0,  |M_ij|=1/sqrt(n) (i!=j),  (M^2)_ii=q.

We first assume `||M||op<=L` for a fixed `L>=1`. This assumption is
removed in Section 6, not imposed on arbitrary minimizing signings.

The following existing results are reused:

- CORE, Section 4: the two Gaussian phases based on `I+/-M` and the
  Gaussian sign identity. The construction and its baseline are not new.
- [Uniform scalar Gaussianization, Section 6](NOTE_2026-09-06_ORIGINAL_SOURCE_NEAR_FLAT_STRICT_GAIN.md):
  if `R` is a correlation matrix with fixed operator cap `C`, and
  `F=sum_j d_j sign(G_j)` with `max|d_j|<=n^(-1/2)` and
  `sum_j d_j^2<=1`, then uniformly
  `E|F|=sqrt(kappa)*sqrt(E F^2)+o_C(1)`.
- The exact original-Boolean-norm mean-update identity already in the
  September 15 explicit-lower note (Section 5, preserved in commit
  `e8e2af56ee028adb098f0c63443f1054672b0cdd`). Section 5 below repeats its
  short derivation to keep the normalization self-contained; it is not
  claimed as a new result.
- [Same-order spectral regularization, Section 5](NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md):
  for `Gamma=4 pi/log(1+sqrt(2))` and any fixed `K>0`,
  `0<=alpha_n^(K+8)-alpha_n<=2 sqrt(2 Gamma/K)`.

The new implication is the uniform paired-field lower (9). Unlike the
earlier spectral-sign defect estimate, its limiting constant does not
depend on `L` or closeness to a flat spectrum. It uses no family scan,
operator-flatness premise, or limiting spectral law.

## 2. Keep the even covariance and linearize only the small odd part

For `s in {+1,-1}`, define

    R_s=(I+sM)^2/(1+q),   R_0=(I+M^2)/(1+q),
    G_s ~ N(0,R_s),       X_s=sign(G_s),
    C_s=E X_s X_s^T=kappa arcsin[R_s].

Square brackets on arcsine denote entrywise application. Each `R_s`
is a correlation matrix and has operator norm at most `(1+L)^2`.
The two phases can be sampled separately. Define the exact even part

    C_0=kappa arcsin[R_0].

For `n>=64` there are symmetric zero-diagonal error matrices `E_s`
such that

    C_s=C_0+s [2 kappa/(1+q)] M+E_s,
    ||E_s||_F <= L+8.                                    (2)

Here is an explicit estimate. Off the diagonal, write

    a_ij=(M^2)_ij/(1+q),    b_ij=s 2M_ij/(1+q).

By Cauchy--Schwarz on rows of M, `|a_ij|<=q/(1+q)<=1/2`.
Also `|b_ij|<=2/sqrt(n)<=1/4`. For `f=arcsin`,

    |f'(a)-1|<=a^2 (|a|<=1/2),
    |f''(z)|<=3   (|z|<=3/4).

Taylor's formula therefore gives

    |arcsin(a+b)-arcsin(a)-b| <= |b| a^2+2b^2.

Now `sum_(i!=j) a_ij^2<=||M^2||_F^2<=L^2 n` and
`sum a_ij^4<=sum a_ij^2/4`. Hence the two terms have Frobenius
norm at most `L` and `8`, respectively. Multiplication by `kappa<1`
proves (2). The diagonal is exact, so no endpoint Taylor expansion
at a diagonal correlation of one is used.

Set

    b=kappa/(1+q),   c=1-kappa arcsin(q/(1+q)).

The following is an exact PSD inequality:

    C_0 >= c I+b M^2.                                    (3)

Indeed, `C_0=c I+kappa arcsin[M^2/(1+q)]` entrywise, including
the diagonal. The odd arcsine series has nonnegative coefficients;
every Schur power of the PSD matrix `M^2/(1+q)` is PSD. Removing
its linear term leaves a PSD sum. Entries have absolute value at
most one half, so the series converges. Since `q/(1+q)<=1/2`,

    c>=1-kappa arcsin(1/2)=2/3.                          (4)

The entire even covariance is retained. No variance-only assertion
about absolute local fields is made at this step.

## 3. Pair the local variances before taking a lower bound

Let `F_(s,i)=(sM X_s)_i` and `sigma_(s,i)^2=E F_(s,i)^2`.
By (2)--(3),

    sigma_(s,i)^2 >= c q+b w_i+s 2b r_i+delta_(s,i),
    w_i=(M^4)_ii,  r_i=(M^3)_ii,
    delta_(s,i)=(M E_s M)_ii,
    (1/n)sum_i |delta_(s,i)| <= L^2(L+8)/sqrt(n).          (5)

The zero diagonal, not a spectral-law assumption, gives the useful
sharpened moment inequality

    w_i>=q^2,        r_i^2<=q(w_i-q^2).                   (6)

For the second statement, write
`r_i=sum_(j!=i) M_ij (M^2)_ji` and use Cauchy--Schwarz. The omitted
diagonal of `M^2` is exactly q, while `sum_(j!=i) M_ij^2=q`.

Put `y_i=sqrt(w_i-q^2)` and `a=cq+bq^2`. Then
`|r_i|<=sqrt(q)y_i` and

    a-bq=q[c-b(1-q)]>0,

because `q>=1/2`, `c>=2/3`, and `b(1-q)<=kappa/3<1/3`.
Thus both `a+b y_i^2 +/- 2b r_i` are nonnegative. Moreover,

    [sqrt(a+b y_i^2+2b r_i)
       +sqrt(a+b y_i^2-2b r_i)]/2 >= sqrt(a).             (7)

For clarity, this is a scalar inequality for every `a>=bq`, `b,q>=0`
and `r^2<=qy^2`. The sum of roots decreases as `|r|` increases,
so it suffices to take `|r|=sqrt(q)y`. The two resulting roots are
the Euclidean norms

    sqrt(b(y+/-sqrt(q))^2+a-bq).

Their average is an even convex function of y and is therefore minimized
at y=0. Equivalently, in the case requiring a second squaring, the
remaining polynomial is `4b y^2(a-bq)>=0`.

Use `sqrt(v)>=sqrt(u)-sqrt(|v-u|)` for nonnegative u,v, or the
corresponding one-sided version when only `v>=u-|delta|` is known.
Cauchy--Schwarz applied to (5) and then (7) proves

    (1/(2n))sum_(s,i) sigma_(s,i)
      >=sqrt(cq+bq^2)-L sqrt(L+8)n^(-1/4).                (8)

This pairing is the new estimate: it does not discard the shared even
term or replace the two cubic moments by two unrelated worst cases.

## 4. Transfer to the actual local fields and their phase energy

Apply the reused scalar Gaussianization lemma to each actual row of
`sM` and to its own `R_s`. Its coefficient hypotheses hold exactly,
and its frame cap depends only on the fixed L. No Gaussian limit for
the entire vector of local fields is required. With

    f=(1/(2n))sum_(s,i) E|F_(s,i)|,
    e=(1/(2n))sum_s E Q_(sM)(X_s),

equation (8) gives, uniformly over the fixed-L class,

    f >= sqrt(kappa(cq+bq^2))-o_L(1)
      = sqrt(kappa(2/3+kappa/2))-o_L(1).                  (9)

Here q tends to one, `c` tends to `2/3`, and b tends to `kappa/2`.
The error may depend on L, but the limiting lower bound does not.

For the same actual phases the standard arcsine difference gives exactly

    e >= kappa q/(1+q).                                  (10)

In detail, `R_+-R_-=4M/(1+q)`, and `(arcsin)' >=1` gives
`M_ij(C_+-C_-)_ij>=4 kappa M_ij^2/(1+q)`. Sum ordered pairs and
divide by `4n`. Since `tr M^2=nq`, this is (10).

## 5. Reuse the original-norm mean update

For either phase, put `Y_s=sign(sM X_s)` and
`z_s=(1-p)X_s+pY_s`, `0<=p<=1`. Multilinearity and zero diagonal
give `|Q_(sM)(z_s)|<=n alpha`. Expanding the quadratic gives

    Q_(sM)(z_s)=(1-p)^2 Q_(sM)(X_s)
        +p(1-p)||(sM)X_s||_1+p^2 Q_(sM)(Y_s).

Since `Q_(sM)(Y_s)>=-n alpha`, average over both phases to obtain
the previously established exact inequality

    (1+p^2)alpha >= (1-p)^2 e+p(1-p)f.                   (11)

At `p=1/10`, (9)--(11) prove for every fixed `L>=1`

    alpha >= B-o_L(1),
    B=[81 kappa/2+9 sqrt(kappa(2/3+kappa/2))]/101.         (12)

There is no penalty proportional to L in B.

## 6. Remove the operator cap at the original minimum

For every fixed `K>0`, same-order regularization gives complete signings
with operator norm at most `(K+8)sqrt(n)` and normalized norm at most
`alpha_n+2 sqrt(2 Gamma/K)`. Apply (12) at the fixed cap `L=K+8`:

    liminf_n alpha_n >= B-2 sqrt(2 Gamma/K).

Now let K tend to infinity. This proves the first inequality in (1).
The order of limits is essential: no uniform Gaussianization error for
a growing cap is assumed, and the old source energy is not replaced
by the energy of a small perturbation.

For an exact conservative certificate of the second inequality, reuse
`7/11<kappa<16/25`. Then

    kappa(2/3+kappa/2)>455/726>(791/1000)^2,
    B>[81(7/22)+9(791/1000)]/101>13/40.                   (13)

The exact differences and the scalar polynomial underlying (7) are
checked on the mesh in `paired_polynomial_field_20260917/`. That
computation does not formalize Gaussianization, regularization, or the
matrix argument above. Root author-review is recorded there; no
independent human or proof-assistant verification is claimed.

## 7. Duplication check and scope

The exact-object check included the ten local branch heads, all four
worktree locations, the source Gaussian and adaptive-update notes, and
the September 6 campaign's unpublished analytic frontier. The old
midpoint arcsine linearization concerns tensor-indexed cross-block laws,
not (2)--(9). The earlier polynomial Gaussian baseline and mean-update
identity are explicitly reused, and the already-completed September 15
constant check is not new evidence for this result.

The contribution is (2)--(9) followed by the cap-independent limit in
(12). It upgrades the unconditional lower bound, but supplies no
cross-order upper comparison and does not settle convergence. Neither
the optional two-ray criterion nor the optimized-pressure transport
inequality is claimed here.
