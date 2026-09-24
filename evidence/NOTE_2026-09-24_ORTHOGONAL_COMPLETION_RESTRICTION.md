# Orthogonal-completion restriction identity

2026-09-24. Deterministic cross-order structure for complete signings.
The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N, let
S be any principal set of size n, let T=S^c have size t, and put F=Phi(A).
Write B=A[T,S] for the t-by-n cross block.

For EVERY Boolean x on S,

    Q_(A[S])(x)^2 + binom(t,2) + ||B x||_2^2 <= F^2.       (1)

This is stronger than ordinary restriction because it retains the entire
orthogonal completion variance.

Proof. Complete x by an independent uniform Boolean vector Y on T. Then

    Q_A(x,Y)=q + L(Y)+Q_(A[T])(Y),

where q=Q_(A[S])(x), L(Y)=(Bx).Y is degree one in Y, and Q_(A[T])(Y)
is degree two. Under the uniform cube these three Walsh levels are
orthogonal. Hence

    E_Y Q_A(x,Y)^2
      = q^2 + ||Bx||_2^2 + binom(t,2).                     (2)

Since |Q_A(x,Y)|<=F pointwise, (1) follows.

Consequences:

1. For an extremizer x_* of the principal restriction,

       Phi(A[S])^2
       <= F^2-binom(t,2)-||B x_*||_2^2.                    (3)

   Thus a restriction can remain close to the host norm only if its
   extremizing Boolean vector is an approximate Boolean kernel vector of
   the omitted cross block.

2. Dropping the cross term gives the universal square-gap

       Phi(A[S]) <= sqrt(F^2-binom(t,2)).                  (4)

3. For one-vertex deletion T={v},

       Phi(A-v)^2 + |a_v.x_*|^2 <= F^2.                    (5)

   Therefore exact norm preservation Phi(A-v)=F is possible only when
   the deleted row is exactly orthogonal to a Boolean extremizer of A-v.

4. Plateau rigidity. If m_(n+1)=m_n and A is ANY exact order-(n+1)
   minimizer, then every principal n-vertex restriction has norm m_n:
   admissibility gives m_n<=Phi(A-v)<=Phi(A)=m_n. Hence for every vertex v
   there is a Boolean extremizer x^(v) of A-v satisfying

       a_v . x^(v) = 0.                                    (6)

   In particular n must be even, recovering the parity obstruction, but
   (6) is stronger: every row of every plateau minimizer has an exact
   Boolean zero against an extremizer of its complementary principal minor.

5. Combining (1) with the September-24 two-sided product restriction gives
   two independent mandatory payments under deletion: internal omitted
   two-sided energy and cross-block completion variance. Neither payment
   was present in ordinary monotonicity.

Equation (1) is exact and finite; no Gaussian approximation, operator cap,
or optimizer assumption is used.
