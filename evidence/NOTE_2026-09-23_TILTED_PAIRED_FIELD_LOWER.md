# Tilted paired Gaussian phases improve the unconditional lower bound

2026-09-23. This strengthens the September 17 paired-field lower bound.
It introduces one fixed scalar tilt in the already-used Gaussian phase
construction, then optimizes the same mean-update inequality. No new family
hypothesis, cross-order assumption, or finite census is used. The original
convergence problem remains OPEN.

Let

    kappa = 2/pi,
    M = A/sqrt(n),
    q = (n-1)/n,
    alpha = Phi(A)/n^(3/2),

and first assume the same fixed normalized operator cap ||M||op<=L used
in the September 17 proof. For a fixed 0<t<=1 and s in {+1,-1}, define

    R_s(t) = (I+s t M)^2 / (1+t^2 q),
    G_s ~ N(0,R_s(t)),
    X_s = sign(G_s).

Every R_s(t) is a correlation matrix.

## 1. Even covariance and the tilted odd signal

Put

    d = 1+t^2 q,
    R_0(t) = (I+t^2 M^2)/d,
    C_s = E[X_s X_s^T] = kappa arcsin[R_s(t)],
    C_0 = kappa arcsin[R_0(t)].

Entrywise arcsine is meant throughout. Since 0<t<=1, every off-diagonal

    a_ij = t^2 (M^2)_ij/d

has |a_ij|<=1/2, while

    b_ij = s 2t M_ij/d

is O(n^(-1/2)). The same Taylor estimate used in the September 17 proof
therefore gives symmetric zero-diagonal E_s with

    C_s = C_0 + s [2 kappa t/d] M + E_s,
    ||E_s||F = O_(L,t)(1).                                  (1)

No new CLT is used here; this is the same bounded-cap arcsine expansion
with the scalar t retained.

Set

    z = t^2 q/d,
    c = 1-kappa arcsin(z),
    b = kappa t^2/d.

Because t^2 M^2/d is PSD and the odd arcsine series has nonnegative
coefficients,

    C_0 >= c I + b M^2.                                     (2)

## 2. Paired local-field variance

Let

    F_(s,i) = (s M X_s)_i,
    sigma_(s,i)^2 = E F_(s,i)^2,
    w_i = (M^4)_ii,
    r_i = (M^3)_ii.

From (1)-(2),

    sigma_(s,i)^2
      >= c q + b w_i + s [2 kappa t/d] r_i + delta_(s,i),   (3)

with

    (1/n) sum_i |delta_(s,i)| = o_(L,t)(1).

Completeness gives the same moment inequality as before,

    w_i >= q^2,
    r_i^2 <= q(w_i-q^2).                                   (4)

Write y_i=sqrt(w_i-q^2) and

    a_q(t)=c q+b q^2.

The average of the two square roots in (3), before delta-errors, is at
least sqrt(a_q(t)). To see this, put h=2 kappa t/d. It suffices that

    h^2 q <= 4 a_q(t) b.                                   (5)

Now h^2 q/(4b)=kappa q/d, while

    a_q(t)/q = c+bq
             >= 2/3
             > kappa
             >= kappa/d,

because z<=1/2 gives c>=2/3. Hence (5) holds. The same scalar argument as
in the September 17 note then yields

    (1/(2n)) sum_(s,i) sigma_(s,i)
       >= sqrt(a_q(t))-o_(L,t)(1).                          (6)

The reused uniform scalar Gaussianization lemma gives

    f_t := (1/(2n)) sum_(s,i) E|F_(s,i)|
       >= sqrt(kappa a(t))-o_(L,t)(1),                      (7)

where, after n->infinity,

    z_t = t^2/(1+t^2),
    a(t)=1-kappa arcsin(z_t)+kappa z_t.                     (8)

## 3. Tilted phase energy

For every off-diagonal entry the two arcsine arguments are a_ij+/-b_ij.
Since arcsin' >=1 and b_ij has the sign of M_ij,

    M_ij (C_+-C_-)_ij >= [4 kappa t/d] M_ij^2.

Therefore

    e_t := (1/(2n)) sum_s E Q_(sM)(X_s)
       >= kappa t q/d,                                     (9)

and hence asymptotically

    e(t)=kappa t/(1+t^2).                                  (10)

At t=1, (7) and (10) reduce exactly to the September 17 constants.

## 4. Mean-update optimization

The already-proved Boolean mean-update inequality gives, for every fixed
p in [0,1],

    (1+p^2) alpha
      >= (1-p)^2 e_t + p(1-p) f_t.                         (11)

Thus, after n->infinity at fixed L and then removing the cap by the same
same-order regularization theorem,

    liminf alpha_n >= B(t)                                 (12)

for every fixed 0<t<=1, where

    e=e(t),
    f=sqrt(kappa a(t)),

and optimization over p gives the exact closed form

    B(t)
      = e-f/2 + (1/2)*sqrt(f^2+(f-2e)^2),                  (13)

with

    p_*(t) = (f-2e)/(f+sqrt(f^2+(f-2e)^2))                 (14)

whenever f>2e. Near t=1 this lies in (0,1), so it is admissible.

## 5. Strict improvement over t=1

The earlier optimized bound is B(1)=0.3258474004377943....

The phase-energy term has e'(1)=0. On the other hand, if
z=t^2/(1+t^2), then

    a'(1)
      = (kappa/2) [1-1/sqrt(1-(1/2)^2)] < 0,

so f'(1)<0. At the t=1 optimizer p_*(1)>0, the envelope derivative is

    B'(1)
      = p_*(1)(1-p_*(1)) f'(1)/(1+p_*(1)^2) < 0.           (15)

Hence decreasing t slightly below one strictly raises the unconditional
lower bound. This is an analytic strict improvement, not a numerical
search artifact.

For the explicit rational choice

    t = 993/1000,

equations (8), (10), and (13) give

    z = 0.4964877503022332...,
    e = 0.3183020328103728...,
    f = 0.7920062144464151...,
    p_* = 0.0971801303973218...,

and

    B(993/1000)
      = 0.3258530333538241....                              (16)

Thus

    liminf_(n->infinity) alpha_n
      >= 0.3258530333538241...,                             (17)

an increase of about 5.63291603e-6 over the already optimized t=1 bound,
and about 1.227789825e-5 over the original p=1/10 September 17 constant.

The displayed decimal is only a numerical evaluation of the exact
expression (13) at t=993/1000; the theorem is the exact expression itself.

## 6. Scope

This is a stronger unconditional asymptotic lower bound for the original
minimum. It does not prove convergence, give a finite-order cutoff, identify
the limit, or impose a special minimizing family. The only new ingredient is
retaining and optimizing a fixed tilt t in the same paired Gaussian source
used by the September 17 proof. The cap-removal order of limits is unchanged.
