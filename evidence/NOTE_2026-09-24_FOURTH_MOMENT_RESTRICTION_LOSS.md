# Sharp fourth-moment one-sided completion loss

2026-09-24. Deterministic cross-order theorem. This strictly supersedes the
same-day hypercontractive constant 1/(2e^2). The original convergence problem
remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order N. Let S be any
principal set of size n=N-t, with T=S^c and t>=2. Put F=Phi(A), and let
B=A[T,S].

Define

    c4 = (sqrt(18)-sqrt(14))/2
       = (3-sqrt(7))/sqrt(2)
       = 0.250491650... .

Then for an extremizer x of the principal restriction,

    boxed:
    Phi(A[S])
      <= F-c4 sqrt( binom(t,2)+||Bx||_2^2 ).                (1)

In particular,

    boxed:
    Phi(A[S]) <= F-c4 sqrt(binom(t,2)).                     (2)

Hence for the exact minima,

    boxed:
    m_N-m_(N-t) >= c4 sqrt(binom(t,2)),                     (3)

and therefore, by integrality,

    boxed:
    m_N-m_(N-t) >= ceil(c4 sqrt(binom(t,2))).               (4)

As t->infinity,

    boxed:
    m_N-m_(N-t)
      >= ((3-sqrt(7))/2+o(1)) t
      = (0.177124344...+o(1)) t.                            (5)

The previous host-free slope was 1/32=0.03125; (5) is more than 5.6 times
larger.

## 1. Fourth moment lemma for degree at most two

LEMMA 1. If Z is a nonzero mean-zero multilinear polynomial of degree at most
2 on a uniform Boolean cube, then

    E Z^4 <= 15 (E Z^2)^2.                                  (6)

Proof. Write

    Z(eps)=sum_i b_i eps_i + sum_(i<j) c_ij eps_i eps_j.

Replace the Rademacher coordinates one at a time by independent standard
Gaussians. Conditional on all other coordinates the polynomial is affine,
a X+b, and

    E_eps (a eps+b)^4 = a^4+6a^2b^2+b^4,
    E_g   (a g+b)^4   =3a^4+6a^2b^2+b^4.

Thus the fourth moment can only increase under the replacement. The second
moment is unchanged because the degree-one and degree-two monomials remain
orthogonal with the same L2 norms.

For the resulting Gaussian polynomial, diagonalize the symmetric quadratic
coefficient matrix. Since its trace is zero, the centered polynomial becomes
a sum of independent variables

    Z_i = a_i g_i + d_i(g_i^2-1).

For one component,

    v_i=E Z_i^2=a_i^2+2d_i^2,

and direct Gaussian moments give

    E Z_i^4
      =3(a_i^4+20a_i^2d_i^2+20d_i^4)
      <=15(a_i^2+2d_i^2)^2
      =15v_i^2.

For independent centered components,

    E(sum_i Z_i)^4
      =sum_i E Z_i^4 + 6 sum_(i<j) v_i v_j
      <=15(sum_i v_i)^2.

This proves (6).

## 2. Convert the fourth moment into a one-sided extremum

LEMMA 2. If X is mean zero, E X^2=1, E X^4<=15, then

    max X >= c4,
    -min X >= c4.                                           (7)

Proof. Let M=max X. If M>=1 there is nothing to prove. Assume M<1. For every
real a,

    (M-X)(X+a)^2 >=0.

Taking expectations gives

    E X^3 <= M(1+a^2)-2a.

Choose a=1/M:

    E X^3 <= M-1/M <0.                                      (8)

Pearson's moment inequality follows from
E[(X^2-(E X^3)X-1)^2]>=0:

    E X^4 >= 1+(E X^3)^2.                                  (9)

Therefore

    15 >= 1+(1/M-M)^2.

Solving the positive quadratic inequality gives

    M >= (sqrt(18)-sqrt(14))/2=c4.

Apply the same argument to -X for the lower extremum. QED.

Combining Lemmas 1 and 2, every nonzero mean-zero degree-at-most-two Boolean
polynomial Z satisfies

    max Z >= c4 ||Z||_2,
    -min Z >= c4 ||Z||_2.                                  (10)

## 3. Principal completion

Choose x on S and sigma in {+1,-1} with

    sigma Q_(A[S])(x)=Phi(A[S])=:G.

Complete x by a uniform Boolean vector Y on T. Then

    sigma Q_A(x,Y)=G+Z(Y),

where

    Z(Y)=sigma[(Bx).Y+Q_(A[T])(Y)].

The linear and quadratic Walsh levels are orthogonal, so

    E Z=0,
    ||Z||_2^2=||Bx||_2^2+binom(t,2).                       (11)

Since sigma Q_A(x,Y)<=F for every Y,

    G+max Z<=F.

Equation (10) now gives (1), and the remaining statements follow.

## 4. Rigidity near equality

If a restriction loses only delta=F-Phi(A[S]), then (1) gives

    ||Bx||_2^2
      <= (delta/c4)^2-binom(t,2).                           (12)

Thus near equality forces an approximate Boolean kernel vector for the whole
omitted cross block. Unlike the earlier square-gap identity, this produces a
direct linear-in-t norm loss.

This theorem uses only finite Boolean/Gaussian moment identities and Pearson's
inequality. No Gaussian approximation, operator cap, optimizer classification,
or finite census is used.
