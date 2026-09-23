# One-vertex extension from near-extremal partition mass

2026-09-23. Author derivation. This is a positive reduction toward the
original convergence problem. It does not assert the partition-mass
hypothesis below.

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), and define the Boolean deficit

    D_A(x)=F-|Q_A(x)| >= 0

and its low-temperature deficit partition sum

    S_A(lambda)=sum_(x in {+/-1}^n) exp(-lambda D_A(x)),   lambda>0.

## 1. Exact extension objective

If a in {+/-1}^n is the incident row of one new vertex, then maximizing
over the new spin y gives

    Phi(A extended by a)
      = max_x ( |Q_A(x)| + |a.x| )
      = F + max_x ( |a.x|-D_A(x) ).                       (1)

Thus the one-vertex problem is exactly a competition between row
correlation and deficit from the old Boolean norm.

## 2. Partition extension lemma

For every lambda>0 there exists a sign row a such that

    Phi(A extended by a)-F
      <= [ n log cosh(lambda) + log(2 S_A(lambda)) ]/lambda. (2)

Proof. For fixed a, write R_a=max_x(|a.x|-D_A(x)). Then

    exp(lambda R_a)
      <= sum_x exp(-lambda D_A(x)) exp(lambda |a.x|)
      <= sum_x exp(-lambda D_A(x))
             [exp(lambda a.x)+exp(-lambda a.x)].

Average over a with independent uniform signs. For every fixed x,

    E_a exp(+/- lambda a.x)=cosh(lambda)^n.

Hence

    E_a exp(lambda R_a)
      <= 2 cosh(lambda)^n S_A(lambda).

At least one a is no larger than this average; taking logarithms proves
(2). No independence or structural assumption on the maximizers of A is
used.

Using log cosh(lambda)<=lambda^2/2 gives the simpler bound

    Phi(A extended by a)-F
      <= n lambda/2 + log(2 S_A(lambda))/lambda.          (3)

## 3. A concrete monotonicity certificate

Write alpha=F/n^(3/2) and choose

    lambda_n = 3 alpha/(2 sqrt(n)).

If

    log(2 S_A(lambda_n)) <= 9 alpha^2/8,                  (4)

then (3) gives a row a with

    Phi(A extended by a)-F <= (3/2) alpha sqrt(n).        (5)

Since

    (n+1)^(3/2)-n^(3/2) > (3/2) sqrt(n),

(5) is strictly below the neutral normalized increment. Therefore, if A
is an exact order-n minimizer satisfying (4),

    alpha_(n+1) < alpha_n.                                (6)

More generally, any excess epsilon_n in the right side of (4) translates
directly through (3) into a quantitative one-vertex excess. This turns the
previous requirement for a single row controlling every near-maximizer
into a scalar bound on the weighted near-extremal mass.

## 4. Why this advances the bottleneck

The finite interior-gap theorem controls softness of fractional
near-extremizers but did not select one row that simultaneously controls
all Boolean near-maximizers. Equations (2)-(4) avoid selecting or covering
those maximizers individually. Their entire contribution is compressed
into S_A(lambda).

At the relevant lambda=Theta(n^(-1/2)), states with a macroscopic
Theta(n^(3/2)) deficit receive exp(-Theta(n)) weight. Thus (4) asks for a
low-temperature entropy estimate near the true Boolean edge, not a bound
on all 2^n states and not a classification of optimizers.

A sufficient all-orders route is now precise: prove (4), or a
Dini-summable relaxation of it, for a choice of exact minimizer at each
large order. Exact (4) would make alpha_n eventually strictly decreasing
and hence convergent. A summable relaxation plugs directly into the
one-vertex neutral-increment criterion already recorded in HANDOFF.

No route is closed here. The next positive target is an upper bound on
S_A(3 alpha/(2 sqrt(n))) for exact minimizers, using pressure/cavity or
near-extremal entropy information already present in the repository.
