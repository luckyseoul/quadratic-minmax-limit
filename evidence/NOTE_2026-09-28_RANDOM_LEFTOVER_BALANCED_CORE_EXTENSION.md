# Random-leftover balanced-core one-vertex extension

2026-09-28. Direct one-vertex extension theorem under the scope freeze.
This removes the deterministic leftover penalty u from the balanced-core
criterion. The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), fix r>0, and let

    H={s : D_s+r<=n}

be the active antipodal stable family.

Choose any core C subset H. Form the projective core-signature classes of
coordinates exactly as in the September-27 balanced-core note. In every odd
signature class remove one coordinate; let U be the set of removed
coordinates and put

    u=|U|.

On every remaining even class G choose canonical signs lambda_i so that
lambda_i x_i^c is constant over i in G for every c in C.

For every active state s and even class G define

    q_(s,G)
      = min( #{i in G: lambda_i x_i^s=+1},
             #{i in G: lambda_i x_i^s=-1} ),

and put

    q_s=sum_G q_(s,G).                                      (1)

For core states c in C one has q_c=0.

## Theorem

If

    boxed:
    sum_(s in H)
      2 exp( -(D_s+r)^2 / [2(u+4q_s)] ) < 1,               (2)

with the convention that a term is zero when u=q_s=0, then there exists an
actual incident sign row a in {+-1}^n such that

    boxed:
    Phi(A extended by a) < F+r.                             (3)

In particular, no condition u<r is required.

## Proof

On every even class G choose independently a uniformly random balanced sign
vector b^G with sum zero and set

    a_i=lambda_i b_i^G.

On every leftover coordinate i in U choose a_i independently and uniformly
from {+-1}. All class choices and leftover signs are mutually independent.

Fix an active state s.

For an even class G, after possibly flipping all canonical residual signs,
let M be the minority set, |M|=q_(s,G). Because sum_i b_i^G=0,

    sum_(i in G) a_i x_i^s
      = -2 sum_(i in M) b_i^G.                              (4)

Sampling without replacement from a balanced +-1 population gives the mgf
bound

    E exp(theta * class contribution)
      <= exp( 4 q_(s,G) theta^2 / 2 ).                     (5)

Thus the even-class contribution is subgaussian with variance proxy
4q_(s,G).

The leftover contribution is

    sum_(i in U) a_i x_i^s,

a sum of u independent Rademachers, hence subgaussian with variance proxy u.

By independence across all classes and leftovers, the full dot product
a.x^s is subgaussian with variance proxy

    v_s=u+4q_s.                                             (6)

Therefore

    P( |a.x^s| >= D_s+r )
      <= 2 exp( -(D_s+r)^2/(2v_s) ).                        (7)

If v_s=0 then a.x^s=0 deterministically.

Summing (7) over all active stable classes and using (2), the union bound
gives positive probability that

    |a.x^s| < D_s+r

for every active s simultaneously.

Inactive states satisfy D_s+r>n>=|a.x^s| automatically. The exact
stable-state extension identity therefore gives (3). QED.

## Why this is strictly more flexible than the previous balanced-core theorem

The September-27 theorem fixed arbitrary signs on U and paid their worst-case
absolute contribution u. Its residual threshold was

    tau_s=D_s+r-u,

and it required u<r.

Here the leftovers are randomized and cost only subgaussian variance u.
The threshold remains the FULL

    D_s+r.

Thus a large number of odd signature classes is no longer fatal.

For M=|H|, if every active state satisfies q_s<=q and D_s+r>=T, the simple
uniform certificate is

    boxed:
    T^2 > 2(u+4q) log(2M).                                  (8)

For core states q_s=0, so their cost is only

    2 exp(-(D_s+r)^2/(2u)).

At the neutral increment r_n=Theta(sqrt(n)), this allows

    u = o(n/log |C|)

for the core constraints, rather than the old deterministic requirement

    u<r_n=Theta(sqrt(n)).

Hence one may take a core large enough to create far more than sqrt(n)
projective signatures while still retaining a possible neutral extension.
This materially enlarges the signature-refinement regime available to the
one-vertex route.

## Neutral and summable-error consequences

For an exact order-n minimizer take

    r_n=m_n[(1+1/n)^(3/2)-1].

If some core C satisfies (2) with r=r_n, then

    boxed:
    alpha_(n+1)<alpha_n.                                    (9)

More generally use r=r_n+rho_n. If (2) holds for all sufficiently large n
and

    sum_n rho_n/(n+1)^(3/2)<infinity,

then the frozen one-vertex argument proves convergence of alpha_n.

## Finite discriminator

Given an exact minimizer, one now needs only:

1. the active stable classes and deficits;
2. a candidate core C;
3. the number u of odd projective signature classes;
4. the impurity q_s of every active state.

The scalar in (2) is then a complete extension certificate. No pairing
optimization, dependency graph, or 2^n row search is required.
