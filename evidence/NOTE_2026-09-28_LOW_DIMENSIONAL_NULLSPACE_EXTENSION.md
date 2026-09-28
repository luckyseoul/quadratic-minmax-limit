# Low-dimensional nullspace rounding for one-vertex extension

2026-09-28. Direct one-vertex extension theorem under the scope freeze.
This gives a new compression route: the active stable family need not have
small cardinality, sparse signature overlap, or small reciprocal-square mass.
It is enough that the active states lie close to a low-dimensional linear
subspace.

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), fix a desired increment r>0, and let

    H={s : D_s+r<=n}

be the active antipodal stable family. Choose Boolean representatives
x^s in {+-1}^n and write

    t_s=D_s+r.

Let V be ANY d-dimensional linear subspace of R^n, with 0<=d<n, and define

    delta_s = dist_2(x^s,V)
            = ||P_(V^perp) x^s||_2.                        (1)

Put

    m_s = t_s-sqrt(n) delta_s.                              (2)

## Theorem 1: approximate-rank extension certificate

Assume m_s>0 for every active state and

    boxed:
    sum_(s in H) 2 exp(-m_s^2/(2d)) < 1                    (3)

when d>=1. Then there exists an actual incident row
a in {+-1}^n such that

    boxed:
    Phi(A extended by a) < F+r.                             (4)

For d=0 the statement is interpreted separately: if every delta_s=0 then H
is empty, since no nonzero Boolean state lies in the zero subspace.

## Proof

Consider the polytope

    P = V^perp intersect [-1,1]^n.

It contains 0. Choose an extreme point z of P. Let

    J={i:|z_i|<1}

be its fractional coordinates and f=|J|.

Because V^perp is defined by d independent linear equations, one has

    f<=d.                                                    (5)

Indeed, if f>d, there is a nonzero vector h supported on J and orthogonal to
V. For sufficiently small epsilon>0, both z+epsilon h and z-epsilon h remain
inside the cube and in V^perp, contradicting extremality.

Set a_i=z_i on every integral coordinate. Independently round each
fractional z_i, i in J, to a_i in {+-1} with

    E a_i=z_i.

Fix an active state s. Since z lies in V^perp,

    x^s.z
      =(P_(V^perp)x^s).z,

hence

    |x^s.z|
      <=delta_s ||z||_2
      <=sqrt(n) delta_s.                                    (6)

The rounding error is

    R_s=sum_(i in J) (a_i-z_i)x_i^s.                        (7)

Its summands are independent, mean zero, and each lies in an interval of
length two. Hoeffding's lemma gives

    E exp(theta R_s) <= exp(f theta^2/2)
                      <= exp(d theta^2/2),

so

    P(|R_s|>=u) <= 2 exp(-u^2/(2d)).                        (8)

If |R_s|<m_s, then by (2) and (6),

    |a.x^s|
      <=|x^s.z|+|R_s|
      <sqrt(n)delta_s+m_s
      =D_s+r.                                               (9)

The union bound and (3) give positive probability that (9) holds for every
active stable state simultaneously. Inactive states are automatically
harmless because D_s+r>n>=|a.x^s|. The exact stable-state extension identity
therefore proves (4). QED.

## Theorem 2: exact-rank corollary

Let d_H be the real linear rank of the active stable-state matrix whose rows
are the x^s. Taking

    V=span{x^s:s in H}

gives delta_s=0 for every active state. Hence

    boxed:
    sum_(s in H)
      2 exp(-(D_s+r)^2/(2d_H)) < 1                         (10)

is sufficient for one-vertex extension below F+r.

In shell form, if M_D active states have deficit D,

    boxed:
    sum_D 2 M_D exp(-(D+r)^2/(2d_H)) < 1.                  (11)

Thus cardinality is exponentially discounted by deficit and by the ACTUAL
linear dimension of the stable family rather than by n.

## Useful uniform corollary

If M=|H| and every active state has threshold at least T, then

    boxed:
    d_H < T^2/[2 log(2M)]                                  (12)

implies extension below F+r.

At the neutral increment of an exact minimizer,

    r_n=m_n[(1+1/n)^(3/2)-1]=Theta(sqrt(n)),

so neutral extension follows whenever, for example,

    d_H = o(n/log(2|H|)).                                   (13)

More generally, the approximate-rank theorem permits d much smaller than
d_H provided every active state is within

    (D_s+r)/sqrt(n)

of V, with the remaining margin satisfying (3).

## Spectral discriminator

Let X be the matrix of active stable representatives. A natural finite test
is to choose V as the span of the top d right singular vectors of X. Then
delta_s is the rowwise residual norm after the rank-d approximation.

For each d=1,...,n-1 compute

    S_d
      =sum_(s in H)
         2 exp(
           -[D_s+r-sqrt(n)delta_s]_+^2/(2d)
         ),

discarding d if any D_s+r<=sqrt(n)delta_s.

Any d with

    S_d<1

is already a complete one-vertex extension certificate.

This is a much smaller computation than searching 2^n incident rows and is
directly aligned with the frozen target.

## Neutral and summable-error consequence

For an exact order-n minimizer, use r=r_n. If (3) holds for some V, then

    boxed:
    alpha_(n+1)<alpha_n.                                    (14)

More generally use r=r_n+rho_n. If such a subspace certificate exists for
all large n and

    sum_n rho_n/(n+1)^(3/2)<infinity,

then the frozen one-vertex argument proves convergence of alpha_n.

This theorem is finite and deterministic apart from the final elementary
rounding existence argument. It does not assume a spectral cap on A, a
Gaussian approximation, a bound on the number of stable states, or any
optimizer classification.
