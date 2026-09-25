# Restriction-bootstrapped cubic interior gap

2026-09-24. Positive same-order strengthening. This couples the optimized
two-sided product floor to principal restriction BEFORE applying the
fractional-spin interior argument. The original convergence problem remains
OPEN.

Let A be a complete symmetric zero-diagonal signing of order n, put

    F=Phi(A)=alpha n^(3/2),
    M=A/sqrt(n),

and assume the universal lower bound alpha>=1/pi-o(1).

For k>=2 let

    Delta_k^*
      =(1/4) max_(1<=b<=k-1) (k-b)^2 mu_b^2,

where mu_b=E|eps_1+...+eps_b|. As proved in
NOTE_2026-09-24_OPTIMIZED_PRODUCT_SPLIT,

    Delta_k^*=[c_*+o(1)]k^3,
    c_*=2/(27 pi).                                         (1)

## 1. Bootstrap the one-sided floor through the complement

Fix T subset [n], |T|=k, and S=T^c. Write

    H_T=Phi(A[T]).

The one-block strict restriction theorem, applied with T retained and S
deleted, gives

    H_T <= F-Delta_(n-k)^*/F.                               (2)

The optimized product theorem inside T gives

    P_T N_T >= Delta_k^*.

Since max(P_T,N_T)=H_T,

    boxed:
    min(P_T,N_T)
      >= Delta_k^*/(F-Delta_(n-k)^*/F).                     (3)

This is strictly stronger than the previous denominator F whenever both
blocks have at least two vertices.

If k/n -> u in [0,1], (1)--(3) give

    min(P_T,N_T)
      >= n^(3/2) [ f_alpha(u)-o(1) ],                       (4)

where

    f_alpha(u)
      = c_* alpha u^3 /
        [alpha^2-c_*(1-u)^3].                               (5)

The denominator is positive throughout the relevant range because

    alpha^2/c_* >= 27/(2 pi)+o(1) > 4.

## 2. Convexity of the bootstrapped profile

The function f_alpha is increasing and convex on [0,1].

Indeed, put

    A0=alpha^2/c_* > 2,
    g(u)=u^3/[A0-(1-u)^3],

so f_alpha=alpha g. Direct differentiation gives

    g'(u)
      =3u^2[A0-(1-u)^2]/[A0-(1-u)^3]^2 >=0,                (6)

and

    g''(u)
      =6u N(A0,u)/[A0-(1-u)^3]^3,                           (7)

where with v=1-u,

    N(A0,u)
      =A0^2+A0(2v^3-5v^2+v)+v^4(2-v).

For 0<=v<=1,

    2v^3-5v^2+v >= -2,
    v^4(2-v)>=0,

hence

    N(A0,u)>=A0(A0-2)>0.                                    (8)

Thus Jensen can be used on the ACTUAL random changed-coordinate fraction;
no concentration theorem is needed.

## 3. Strengthened fractional-update inequality

Fix one orientation sigma and a Boolean source X. Let Y be its full best
response, let K be the number of changed coordinates, and for 0<=p<=1/2 set

    z=(1-p)X+pY.

The K changed coordinates have softness 2p. The symmetric perturbation
identity and (3) give, uniformly after normalization,

    Q_(sigma M)(z)
      <= n alpha-4p^2 n f_alpha(K/n)+o(n).                  (9)

Average over the two phases and source laws. Let

    e=(1/(2n)) sum_sigma E Q_(sigma M)(X_sigma),
    f=(1/(2n)) sum_sigma E ||M X_sigma||_1,
    mu=(1/(2n)) sum_sigma E K_sigma.

Using convexity and monotonicity,

    (1/2)sum_sigma E f_alpha(K_sigma/n)
      >= f_alpha(mu).                                       (10)

The exact quadratic interpolation therefore yields

    boxed:
    (1+p^2)alpha
      >=(1-p)^2 e+p(1-p)f
        +4p^2 f_alpha(mu)-o(1).                             (11)

The old optimized-product interior correction was

    4c_* p^2 mu^3/alpha.

Equation (11) replaces it by

    4c_* p^2 alpha mu^3/
      [alpha^2-c_*(1-mu)^3],                                (12)

which is strictly larger for every 0<mu<1.

## 4. Explicit improved source certificate

For the tilted paired Gaussian source, retain the already proved functions

    kappa=2/pi,
    z(t)=t^2/(1+t^2),
    a(t)=1-kappa asin(z(t))+kappa z(t),
    e(t)=kappa t/(1+t^2),
    f(t)=sqrt(kappa a(t)),
    beta(t)=a(t)-kappa,
    C(t)=2sqrt(kappa)t/(1+t^2),

and

    H(t)=min{1/4, atan(sqrt(beta(t))/C(t))/pi},

valid for 99/100<=t<=1.

For fixed p in (0,1/2], put

    R(t,p)=(1-p)^2 e(t)+p(1-p)f(t).

Since mu>=H(t), (11) gives the scalar inequality

    boxed:
    (1+p^2)alpha
      >= R(t,p)
       + 4c_* p^2 alpha H(t)^3/
         [alpha^2-c_*(1-H(t))^3].                           (13)

Equivalently,

    [(1+p^2)alpha-R(t,p)]
    [alpha^2-c_*(1-H(t))^3]
      >=4c_*p^2 alpha H(t)^3.                               (14)

Let B_boot(t,p) be the largest positive root above
sqrt(c_*(1-H(t))^3) of equality in (14). Then

    boxed:
    liminf alpha_n
      >= sup_(.99<=t<=1, 0<p<=1/2) B_boot(t,p).             (15)

At the explicit rational parameters

    t=4967/5000=0.9934,
    p=61/625=0.0976,

the root is

    B_boot
      =0.3258785815293334647646857... .                     (16)

Thus (15) supplies the explicit unconditional certificate

    boxed:
    liminf alpha_n >= 0.3258785815.                         (17)

The decimal is only a locator for the exact algebraic-transcendental
certificate (14) at the displayed rational parameters.

## 5. Scope

The advance is not another stable-state geometry statement. It consumes the
new cross-order restriction theorem inside the same-order interior argument:
the complement's mandatory product payment reduces the norm available to the
changed principal block, which raises its opposite one-sided extremum. The
convexity calculation then survives averaging without any concentration
assumption.

No finite signing census, Gaussian reapplication after an update, operator
cap, or optimizer classification is used.
