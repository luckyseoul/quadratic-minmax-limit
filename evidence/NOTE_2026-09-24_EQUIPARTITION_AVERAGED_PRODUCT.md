# Equipartition-averaged two-sided product bound

2026-09-24. Positive strengthening of the finite two-sided Boolean product
bound. The original convergence problem remains OPEN.

For an order-k complete symmetric zero-diagonal signing B, write

    P=max_x Q_B(x),   N=-min_x Q_B(x),

and let

    mu_m = E|eps_1+...+eps_m|.

The earlier response/seed argument bounded the internal response block by the
FULL opposite one-sided extremum. Averaging over a random balanced
equipartition cuts that penalty in half.

## 1. Even order

Let k=2m. Choose a uniformly random equipartition V=S disjoint union T with
|S|=|T|=m. Conditional on S,T, sample X_T uniformly and define the response

    Y_i=sign(sum_(j in T) B_ij X_j),   i in S.

For 0<=s<=1, the point (sY_S,X_T) lies in the continuous cube, so

    P >= E Q_B(sY_S,X_T)
      >= s L - s^2 E N(B[S]),                              (1)

where

    L=m mu_m.

For every equipartition,

    N(B[S])+N(B[T]) <= N(B)                                (2)

by the disjoint-phase inequality. Since S and T have identical marginal
laws,

    E N(B[S]) <= N/2.                                      (3)

Thus, for EVERY 0<=s<=1,

    P >= sL-(s^2/2)N.                                      (4)

Applying the same argument to -B gives

    N >= sL-(s^2/2)P.                                      (5)

Normalize p=P/L and n=N/L. Put

    f(u)=max_(0<=s<=1) [s-(u/2)s^2]
        = {1-u/2,  0<=u<=1;
           1/(2u), u>=1}.                                  (6)

Equations (4)--(5) say

    p>=f(n),   n>=f(p).                                    (7)

If p>=1 or n>=1, (6) immediately gives pn>=1/2.
If p,n<1, then

    p+n/2>=1,   n+p/2>=1.

The feasible region in the unit square has minimum product at

    p=n=2/3,

so in every case

    boxed:
    P N >= (4/9)L^2
        = (4/9)m^2 mu_m^2.                                 (8)

The same inequalities at s=1 also give the phase-sum bound

    boxed:
    P+N >= (4/3)m mu_m.                                    (9)

## 2. Odd order

Let k=2m+1. Restrict B to any 2m principal vertices. One-sided principal
restriction gives

    P(B)>=P(B[U]),   N(B)>=N(B[U]).

Therefore (8) on B[U] yields

    boxed:
    P(B)N(B) >= (4/9)m^2 mu_m^2.                           (10)

Thus for every k>=2, with m=floor(k/2), define

    Gamma_k=(4/9)m^2 mu_m^2.

Then

    boxed:
    P(B)N(B) >= Gamma_k.                                   (11)

## 3. Asymptotic constant

Since

    mu_m^2=(2/pi)m+O(1),

and m=k/2+O(1),

    boxed:
    Gamma_k=[1/(9pi)+o(1)]k^3.                             (12)

The immediately preceding optimized asymmetric-split theorem gave only

    [2/(27pi)+o(1)]k^3.

Hence the universal cubic product coefficient improves by the exact factor

    boxed:
    [1/(9pi)]/[2/(27pi)] = 3/2.                            (13)

Relative to the original balanced non-averaged coefficient 1/(16pi), the
improvement factor is 16/9.

## 4. Cross-order consequences

Every theorem that used only a universal lower bound on P(B)N(B) may now use
Gamma_k instead.

In particular, if A is a complete host of norm F and T is any deleted
principal block of size t, then

    min(P(A[T]),N(A[T])) >= Gamma_t/Phi(A[T]) >= Gamma_t/F,

so

    boxed:
    Phi(A[S]) <= F-Gamma_t/F.                              (14)

The stronger denominator min(F,binom(t,2)) can be retained:

    boxed:
    Phi(A[S])
      <= F-Gamma_t/min(F,binom(t,2)).                      (15)

Likewise the coupled two-block quadratic restriction theorem remains valid
with Delta_k^* replaced everywhere by Gamma_k.

For exact minima,

    boxed:
    m_N-m_(N-t)
      >= Gamma_t/min(m_N,binom(t,2)).                      (16)

At proportional deletion t=lambda N with m_N=alpha_N N^(3/2), the
host-scale branch becomes

    m_N-m_(N-t)
      >= [lambda^3/(9pi alpha_N)+o(1)]N^(3/2),             (17)

whenever the host norm is the active denominator.

This is a full N^(3/2)-scale improvement for proportional restriction, not a
linear-in-t correction.

## 5. Duplication audit

Repository searches for the constant 1/(9pi), an equipartition-averaged
response penalty, and the inequality E N(B[S])<=N(B)/2 in this product
argument returned no prior instance. The ingredients (disjoint-phase
restriction and response-spin construction) were known separately; their
averaged combination and the resulting 4/9 product factor are new here.
