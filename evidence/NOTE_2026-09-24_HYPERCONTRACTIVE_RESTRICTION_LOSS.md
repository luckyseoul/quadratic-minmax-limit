# Hypercontractive completion gives a stronger linear restriction loss

2026-09-24. Deterministic cross-order theorem. The original convergence
problem remains OPEN.

This strengthens the host-free t/32 principal-restriction loss proved
earlier today.

Let A be a complete symmetric zero-diagonal signing of order N, let
S be any principal set of size n=N-t with t>=2, let T=S^c, and put

    F=Phi(A).

Write B=A[T,S] for the t-by-n cross block.

## 1. A one-sided lemma for degree-two cube polynomials

LEMMA. If Z is a nonzero mean-zero real polynomial of degree at most 2 on a
uniform Boolean cube, then

    max Z >= ||Z||_2/(2 e^2),
    -min Z >= ||Z||_2/(2 e^2).                              (1)

Proof. For every q>2, Bonami--Beckner hypercontractivity gives

    ||Z||_q <= (q-1)||Z||_2.                                (2)

Interpolating L1 and Lq to L2, with
theta=(q-2)/(2(q-1)), gives

    ||Z||_1 >= ||Z||_2 (q-1)^(-q/(q-2)).                   (3)

Let q decrease to 2. Since
(q-1)^(-q/(q-2)) -> e^(-2),

    ||Z||_1 >= e^(-2)||Z||_2.                              (4)

Because EZ=0,

    E Z_+ = E Z_- = ||Z||_1/2.

Each expectation is at most the corresponding one-sided maximum, proving
(1). QED.

## 2. Apply the lemma to every principal completion

Choose x on S and sigma in {+1,-1} so that

    sigma Q_(A[S])(x)=Phi(A[S])=:G.

Complete x by a uniform Boolean vector Y on T. Then

    sigma Q_A(x,Y)=G+Z(Y),                                  (5)

where

    Z(Y)=sigma[(Bx).Y+Q_(A[T])(Y)].

The linear and quadratic Walsh levels are orthogonal, hence

    EZ=0,
    ||Z||_2^2 = ||Bx||_2^2 + binom(t,2).                    (6)

Since sigma Q_A(x,Y)<=F for EVERY Y,

    G + max_Y Z(Y) <= F.                                    (7)

Using (1) and (6) gives the sharpened restriction theorem

    boxed:
    Phi(A[S])
      <= F - sqrt( binom(t,2)+||Bx||_2^2 )/(2 e^2).         (8)

In particular, dropping the cross-block term,

    boxed:
    Phi(A[S])
      <= F - sqrt(binomial(t,2))/(2 e^2).                   (9)

This is finite, deterministic, and holds for EVERY principal S. No optimizer
assumption, Gaussian approximation, operator cap, or census is used.

## 3. New cross-order inequality

Apply (9) to an exact order-N minimizer. Every retained S is an admissible
order-n signing, so for t=N-n>=2,

    boxed:
    m_N-m_n >= sqrt(binomial(t,2))/(2 e^2).                 (10)

Since the left side is integral,

    boxed:
    m_N-m_n >= ceil( sqrt(binomial(t,2))/(2 e^2) ).         (11)

As t->infinity,

    m_N-m_(N-t)
      >= [1/(2 sqrt(2) e^2)+o(1)] t
      = (0.04784...+o(1)) t.                                (12)

The previous universal host-free theorem gave only t/32=0.03125 t.
Thus the guaranteed all-scale linear restriction slope improves by more
than fifty percent.

The improvement is already strict at t=2:

    1/(2e^2) > 1/16 = 2/32,

and the ratio only increases toward the asymptotic constant above.

## 4. Cross-block rigidity near equality

Equation (8) retains information discarded by (9). If a principal
restriction loses only delta=F-Phi(A[S]), then its extremizer x must satisfy

    ||Bx||_2^2
      <= (2e^2 delta)^2 - binom(t,2).                       (13)

Thus near equality in the cross-order bound forces an approximate Boolean
kernel vector for the omitted cross block. In particular no restriction can
come within sqrt(binomial(t,2))/(2e^2) of the host norm at all.

This complements the orthogonal-completion square identity: the second
moment there gave a quadratic gap, while hypercontractivity converts the
same completion variance into a direct one-sided LINEAR norm loss.
