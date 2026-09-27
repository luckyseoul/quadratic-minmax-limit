# Balanced-core hypergeometric one-vertex extension

2026-09-27. Direct one-vertex extension theorem under the scope freeze.
This replaces worst-pair residual control by a per-state impurity tail.
The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order n,
F=Phi(A), and let H be the active antipodal stable family for a desired
increment r>0:

    H={s : D_s+r<=n}.

Choose a core C subset H. For each coordinate i, let its projective core
signature be the column (x_i^c)_(c in C), modulo global sign. In each
projective signature class G choose canonical signs lambda_i in {+-1} so
that

    lambda_i x_i^c

is independent of i in G for every c in C.

If |G| is odd, remove one arbitrary coordinate from G. Let U be the set of
all removed coordinates and write

    u=|U|.

Thus every remaining class G' has even size.

For a residual state s in H\C define the canonical residual signs

    y_i^s=lambda_i x_i^s,  i notin U,

and the within-class minority count

    q_s
      = sum_(G')
          min( #{i in G': y_i^s=+1},
               #{i in G': y_i^s=-1} ).                     (1)

Put

    tau_s=D_s+r-u.                                          (2)

## Theorem

If

    u<r                                                      (3)

and

    sum_(s in H\C, q_s>0)
      2 exp( -tau_s^2/(8 q_s) ) < 1,                       (4)

then there exists an actual incident sign row a in {+-1}^n such that

    Phi(A extended by a) < F+r.                             (5)

Thus, for an exact order-n minimizer and the neutral increment

    r_n=m_n[(1+1/n)^(3/2)-1],

condition (3)--(4) implies

    alpha_(n+1)<alpha_n.                                    (6)

More generally the same statement with r=r_n+rho_n and
sum rho_n/(n+1)^(3/2)<infinity gives convergence through the existing
one-vertex argument.

## Proof

Fix arbitrary signs on the leftover coordinates U.

On each even class G', independently choose a uniformly random balanced
sign vector b=(b_i)_(i in G') satisfying

    sum_(i in G') b_i=0.

Set

    a_i=lambda_i b_i

on paired coordinates.

For a core state c,

    a_i x_i^c
      = b_i (lambda_i x_i^c),

and the parenthesized factor is constant on G'. Therefore every even class
contributes exactly zero. Only U remains, so

    |a.x^c|<=u<r<=D_c+r.                                    (7)

Now fix a residual active state s. The contribution from U has absolute
value at most u. On one even class G', after possibly replacing all y_i^s
by their negatives, let M be the minority set and q=|M|. Since sum b_i=0,

    sum_(i in G') a_i x_i^s
      =sum_i b_i y_i^s
      =-2 sum_(i in M) b_i.                                (8)

Under a uniformly balanced b, the q entries (b_i)_(i in M) are a sample
without replacement from a population containing equally many +1 and -1.
Hoeffding's comparison for sampling without replacement gives

    E exp(theta sum_(i in M)b_i)
      <= exp(q theta^2/2).                                  (9)

Hence the class contribution in (8) is subgaussian with variance proxy 4q.
The choices on distinct classes are independent. Therefore the total
non-leftover contribution Z_s is subgaussian with variance proxy 4q_s:

    P(|Z_s|>=t)
      <=2 exp(-t^2/(8q_s)),   q_s>0.                        (10)

If q_s=0 then Z_s=0 identically.

Taking t=tau_s and summing (10) over residual states, condition (4) and the
union bound show that with positive probability

    |Z_s|<tau_s

for every residual active s simultaneously. Adding the leftover contribution,

    |a.x^s|
      <tau_s+u
      =D_s+r.                                               (11)

Inactive stable states satisfy D_s+r>n>=|a.x^s| automatically. The exact
stable-state extension identity now gives (5). QED.

## Why this is a real strengthening of the pairing certificate

The September-25 core-signature theorem first fixed a pairing inside every
signature class and then required EVERY chosen pair to have small weighted
residual disagreement mass.

The present theorem does not commit to a pairing. It uses the full balanced
slice of each even signature class. A residual state pays only its total
within-class minority count q_s, and the payment is exponential:

    exp(-(D_s+r-u)^2/(8q_s)).

Consequently:

- a class on which the residual state is constant costs exactly zero;
- different residual states may use different minority locations;
- no single bad coordinate pair controls the certificate;
- large deficits are exponentially discounted rather than only by an
  inverse-square weight.

This is a direct construction of the extension row and is therefore inside
the frozen scope.

## Finite discriminator

For a recorded exact minimizer, the test is now:

1. enumerate active stable antipodal classes and deficits;
2. choose a small core C;
3. form projective signature classes;
4. remove one coordinate from each odd class;
5. compute q_s for every residual active state;
6. evaluate the left side of (4).

A value below 1 is already a complete one-vertex extension certificate.
No search over pairings is required.
