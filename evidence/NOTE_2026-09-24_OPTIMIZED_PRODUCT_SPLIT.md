# Optimized asymmetric product split

2026-09-24. Positive constant improvement to the two-sided Boolean product
bound and every theorem that consumes it. The original convergence problem
remains OPEN.

For an order-k complete signing B, k>=2, write

    P=max_x Q_B(x),   N=-min_x Q_B(x),

and let

    mu_b = E|eps_1+...+eps_b|.

The September-23 product proof used the balanced split
a=ceil(k/2), b=floor(k/2). Balance is not optimal.

For ANY integers a,b>=1 with a+b=k, the identical response-spin argument
gives

    P N >= a^2 mu_b^2 / 4.                                  (1)

Therefore define the optimized finite constant

    Delta_k^*
      := (1/4) max_(1<=b<=k-1) (k-b)^2 mu_b^2.              (2)

Then

    boxed: P N >= Delta_k^*.                                (3)

The old constant Delta_k is one admissible term in (2), so

    Delta_k^* >= Delta_k                                    (4)

at every order, with a strict asymptotic improvement.

## Proof

Fix a response set S of size a and an independent uniform seed set T of
size b. For X_T uniform put

    Y_i=sign(sum_(j in T) B_ij X_j),  i in S.

Then

    L=a mu_b.

For 0<=s<=1, multilinearity on the continuous cube gives

    P >= sL-s^2 N.

Applying the same construction to -B gives

    N >= sL-s^2 P.

At s=1, P+N>=L, hence H=max(P,N)>=L/2. Taking s=L/(2H)<=1 in the inequality
whose quadratic coefficient is -H yields

    PN>=L^2/4=a^2 mu_b^2/4.

Since the split was arbitrary, maximize over b.

## Asymptotic optimizer

The exact Rademacher mean satisfies

    mu_b^2=(2/pi)b+O(1).

Write b=lambda k. Then

    (1/4)(k-b)^2 mu_b^2
      = [lambda(1-lambda)^2/(2pi)+o(1)] k^3.

The cubic profile lambda(1-lambda)^2 is maximized at

    lambda=1/3.

Thus the optimal split is asymptotically

    response : seed = 2 : 1,

and

    boxed:
    Delta_k^*
      = [2/(27 pi)+o(1)] k^3.                               (5)

The old balanced coefficient was

    1/(16 pi).

Hence the exact product coefficient improves by the factor

    boxed:
    [2/(27pi)]/[1/(16pi)] = 32/27 = 1.185185... .           (6)

## Cross-order consequence

Every occurrence of the old Delta_t in the September-24 strict principal
restriction theorem can be replaced by Delta_t^*. In particular, for a
complete host A of norm F and any principal deletion of t>=2 vertices,

    boxed:
    Phi(A[S])
      <= F - Delta_t^*/min(F,binom(t,2)).                   (7)

and for exact minima,

    boxed:
    m_N-m_(N-t)
      >= Delta_t^*/min(m_N,binom(t,2)).                     (8)

The host-scale macroscopic coefficient therefore improves from

    1/(16pi)

to

    2/(27pi).

The independent half-per-deleted-vertex theorem remains stronger in the
small-deletion regime; (8) is the improved proportional-scale branch.

## Interior-gap consequence

The same substitution improves every finite interior-gap statement that was
derived only from the product lower bound. At proportional subset size k,
the asymptotic cubic coefficient changes from

    1/(16pi)

to

    2/(27pi),

without any new distributional, spectral, or optimizer assumption.

This note is a constant optimization of an already-proved mechanism, but it
strictly strengthens its quantitative output. Repository search found no
previous optimization of the response/seed split.
