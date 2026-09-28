# Deterministic low-rank one-vertex extension

2026-09-28. Direct one-vertex extension theorem under the scope freeze.
This strengthens NOTE_2026-09-28_LOW_DIMENSIONAL_NULLSPACE_EXTENSION by
removing the union bound entirely when the active stable family is captured
by a low-dimensional subspace.

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), fix r>0, and let

    H={s : D_s+r<=n}

be the active antipodal stable family. Choose representatives
x^s in {+-1}^n and thresholds

    t_s=D_s+r.

Let V be any d-dimensional subspace of R^n and define

    delta_s=dist_2(x^s,V).

## Theorem

There exists an actual sign row a in {+-1}^n satisfying

    boxed:
    |a.x^s| <= sqrt(n) delta_s + d
    for every active state s.                               (1)

Consequently, if

    boxed:
    sqrt(n) delta_s + d < D_s+r
    for every s in H,                                      (2)

then

    boxed:
    Phi(A extended by a) < F+r.                            (3)

The sharper statement replaces d in (1)--(2) by a fractional-rounding cost
L<=d defined in the proof.

## Proof

Consider the polytope

    P=V^perp intersect [-1,1]^n.

Choose an extreme point z of P. Let

    J={i:|z_i|<1}.

Exactly as in the preceding low-dimensional-nullspace note,

    |J|<=d,                                                 (4)

because otherwise there would be a nonzero perturbation supported on J that
stays inside V^perp.

Round z deterministically to its nearest Boolean vertex:

    a_i=sign(z_i),

with either sign allowed when z_i=0. Define

    L=sum_(i in J) (1-|z_i|).                               (5)

Then

    L<=|J|<=d.                                              (6)

For an active state s,

    |x^s.z|
      =|(P_(V^perp)x^s).z|
      <=sqrt(n) delta_s.                                    (7)

Also

    |x^s.(a-z)|
      <=sum_(i in J)|a_i-z_i|
      =L.                                                   (8)

Therefore

    |a.x^s|
      <=sqrt(n)delta_s+L
      <=sqrt(n)delta_s+d,                                   (9)

which proves (1). If (2) holds, every active stable constraint is strictly
satisfied. Inactive states are automatic because D_s+r>n>=|a.x^s|. The exact
stable-state extension identity gives (3). QED.

## Exact-rank corollary

Let

    d_H=rank_R span{x^s:s in H}.

Taking V equal to the exact active-state span gives delta_s=0 for every
active state. Hence

    boxed:
    d_H < min_(s in H)(D_s+r)                               (10)

is sufficient for extension below F+r.

Since min_s(D_s+r)=r whenever a zero-deficit active state exists, the common
case reduces to

    boxed:
    rank(H) < r.                                            (11)

At the neutral increment

    r_n=m_n[(1+1/n)^(3/2)-1]=Theta(sqrt(n)),

any exact minimizer whose active stable family has real rank below r_n admits
a neutral one-vertex extension and therefore satisfies

    boxed:
    alpha_(n+1)<alpha_n.                                    (12)

## Approximate-rank corollary

More generally, define the rowwise d-dimensional approximation error

    epsilon_d
      = min_(dim V=d) max_(s in H)
          [ sqrt(n) dist_2(x^s,V) - D_s ].                  (13)

Then

    boxed:
    d + epsilon_d < r                                      (14)

is sufficient for extension below F+r.

This is a deterministic rank-width criterion with no dependence on the
number of active stable states.

## Why this is stronger than the probabilistic low-rank certificate

The previous theorem rounded the at most d fractional coordinates randomly
and paid a tail/union-bound condition over the active family.

Nearest-vertex rounding instead pays only the total fractional l1 movement

    L=sum_(i in J)(1-|z_i|)<=d,

simultaneously for every state. Therefore cardinality disappears completely.

The criterion can succeed even when the active stable family is enormous,
provided it is linearly compressible at the threshold scale.

## Finite discriminator

For each d, choose a candidate d-dimensional subspace V, for example from
the leading right singular vectors of the active stable-state matrix, and
compute

    max_s [sqrt(n) dist_2(x^s,V)-D_s].

Any d satisfying (14) is already a complete deterministic one-vertex
certificate. If V is the exact row span, only its rank is needed.

No search over incident rows is required.
