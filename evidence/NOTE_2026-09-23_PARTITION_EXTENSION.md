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

## 3. Correction: the first scalar certificate is impossible

The initially proposed specialization
`log(2 S_A(lambda)) <= 9 alpha^2/8` cannot hold and is withdrawn.
Spin reversal gives at least two Boolean maximizers, so `S_A(lambda)>=2`.
Meanwhile `alpha<=1/2+o(1)` makes the proposed right side at most about
9/32. Thus that specialization was vacuous. The partition-extension lemma
(2) itself remains valid.

A sharper formulation avoids paying the exponential-moment factor at the
maximizers. Identify x and -x, since both D_A and |a.x| are unchanged, and
let H be one representative from each pair. For r>=0, a uniformly random
sign row obeys the exact union estimate

    P_a{ max_x (|a.x|-D_A(x)) >= r }
      <= sum_(x in H) P{|S_n| >= D_A(x)+r},                (4)

where S_n is a sum of n independent uniform signs. Consequently, whenever

    sum_(x in H) P{|S_n| >= D_A(x)+r} < 1,                 (5)

there exists a row with

    Phi(A extended by a)-F < r.                            (6)

This is stronger than applying Hoeffding termwise and, crucially, retains
the exact central-binomial tail at the zero-deficit states.

Proof: for each fixed x, a.x has exactly the S_n distribution. The event
that the maximum exceeds r is the union of the displayed tail events.
Quotienting by spin reversal removes duplicate events. No independence
between different x is asserted or needed.

## 4. Neutral-increment consequence

Put r_n=F*((1+1/n)^(3/2)-1), the exact neutral increment. If the left side
of (5), evaluated at r=r_n, is <1 for an exact minimizer A, then

    alpha_(n+1) < alpha_n.                                (7)

More generally, replacing r_n by r_n+rho_n gives the existing
one-vertex summable-error criterion when sum rho_n n^(-3/2)<infinity.

## 5. Positive content and next implication

Equation (5) is a finite, directly checkable sufficient condition expressed
in the ACTUAL deficit spectrum of a minimizer. Unlike the withdrawn
partition-mass specialization, it is not automatically impossible at the
zero-deficit shell: one maximizer pair contributes the exact binomial tail
P{|S_n|>=r_n}, which is strictly below one for every positive r_n.


The finite interior-gap theorem controls softness of fractional
near-extremizers but did not select one row that simultaneously controls
all Boolean near-maximizers. Equations (2)-(4) avoid selecting or covering
those maximizers individually. Their entire contribution is compressed
into S_A(lambda).

At the relevant lambda=Theta(n^(-1/2)), states with a macroscopic
Theta(n^(3/2)) deficit receive exp(-Theta(n)) weight. Thus (4) asks for a
low-temperature entropy estimate near the true Boolean edge, not a bound
on all 2^n states and not a classification of optimizers.

The remaining mathematical issue is to control the combined exact binomial tails of the near-extremal shells. The lemma above is unconditional and strictly sharpens the earlier scalar reduction by preserving the discrete correlation law instead of an exponential-moment relaxation.
