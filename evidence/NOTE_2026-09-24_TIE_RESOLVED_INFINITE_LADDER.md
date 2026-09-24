# Tie-resolved successor dynamics and an infinite strict lower-bound ladder

2026-09-24. Positive structural strengthening of the actual-successor route.
The original convergence problem remains OPEN. No independent human or
proof-assistant review is claimed.

The September 24 successor-noise theorem deliberately kept the old spin at
a zero local field. That convention is harmless for energy, but it caused
the all-orders changed-coordinate recursion to lose a second factor of p:
mu'>=p^2 mu (with the stronger p mu only when parity forbids zero fields).

Zero-field ties are FREE best-response choices. Choosing the opposite spin
at every zero field removes that loss completely.

The result is more than a better constant: the actual non-Gaussian update
can now be iterated to arbitrary fixed depth with an explicit scalar
recursion, and every finite depth can be extended to a strictly stronger
unconditional liminf certificate. Thus the lower-bound mechanism has no
finite-depth stopping point.

## 1. Tie-flip best response

Let A be a complete symmetric zero-diagonal signing, M=A/sqrt(n), fix an
orientation s, and let X be any Boolean source law. Put F=sMX.

Choose the full best response Y coordinatewise by

    Y_i =  X_i   if X_i F_i > 0,
    Y_i = -X_i   if X_i F_i <= 0.                           (1)

Thus a zero field is deliberately flipped. Let

    S={i:Y_i!=X_i}={i:X_i F_i<=0},   k=|S|,
    R=sum_(i in S)|F_i|.

The zero-field coordinates contribute zero to R. Therefore all identities
using weighted disagreement are unchanged:

    ||MX||_1-2 Q_(sM)(X)=2R,                               (2)

and, because Y_i F_i=|F_i| also at a zero field,

    X^T(sM)Y=||MX||_1.                                     (3)

Independently flip each coordinate in S with probability 0<p<=1/2 and
leave the other coordinates fixed. Call the actual successor X'. At the
successor, use the SAME tie-flip convention (1), producing S' and k'.

## 2. Exact count retention: mu' >= p mu for every n

Fix a source state x and i in S. The successor field F'_i does not depend
on i's own flip coin because M_ii=0. Condition on all other coins and put

    a=x_i F'_i.

If i is not flipped, it belongs to S' exactly when a<=0. If i is flipped,
it belongs to S' exactly when a>=0. Hence

    P(i in S' | all other coins)
      =(1-p) 1_(a<=0)+p 1_(a>=0)
      >=p.                                                  (4)

At a=0 the left side is one, rather than zero. Summing over the ORIGINAL
S gives

    E[k'|x] >= p k.                                         (5)

Averaging the two phases therefore proves the finite cap-free estimate

    boxed:  mu' >= p mu                                    (6)

for every order n, with no parity condition and no error term.

This strictly strengthens the previous universal p^2 mu estimate. The old
estimate remains correct under the old-sign tie convention; (6) uses the
better admissible convention.

## 3. The successor field-noise estimate survives unchanged

For i in S, conditional on all coins except i,

    E[(-X'_i F'_i)_+ | other coins] >= p |F'_i|             (7)

for p<=1/2. This remains true whether the source coordinate entered S
through a negative field or a zero field.

Conditioned on x, for i in S,

    F'_i=a_i+sum_(j in S,j!=i) u_ij (xi_j-p),
    |u_ij|=2/sqrt(n),

with independent Bernoulli(p) coins. The shifted fourth-moment lemma from
NOTE_2026-09-24_SUCCESSOR_FLIP_NOISE therefore applies verbatim. With

    d_p=2sqrt(p(1-p)/3),  E_p=d_p+4/3,

it gives

    g' >= 2p d_p mu^(3/2)-2p E_p/sqrt(n),                  (8)

where g=f-2e. The energy recursion is likewise unchanged,

    e' >= (1-p^2)e+p(1-p)g-p^2 alpha.                      (9)

For the interior count,

    eta' >= (mu')^3-1/n^3 >= p^3 mu^3-1/n^3.              (10)

The old all-orders p^6 propagation in eta is therefore replaced by p^3.

## 4. Arbitrary fixed-depth scalar dynamics

Use the initial tilted paired source at a fixed 99/100<=t<=1. Retain

    e(t), f(t), H(t)>1/5

from the both-phase theorem, where H(t) is the phase-averaged changed-count
floor.

Choose a first update probability p_0 in (0,1/2]. Define

    a_1=(1-p_0)^2 e(t)+p_0(1-p_0)f(t),
    b_1=p_0^2,
    M_1=p_0 H(t),
    G_1=2p_0 d_(p_0) H(t)^(3/2).                           (11)

For j>=1 and any further fixed p_j in (0,1/2], recursively put

    a_(j+1)=(1-p_j^2)a_j+p_j(1-p_j)G_j,
    b_(j+1)=(1-p_j^2)b_j+p_j^2,
    M_(j+1)=p_j M_j,
    G_(j+1)=2p_j d_(p_j) M_j^(3/2).                        (12)

For every FIXED depth j, after n->infinity at the initial fixed source cap,

    e_j >= a_j-b_j alpha-o_L(1),
    mu_j >= M_j-o_L(1),
    g_j >= G_j-o_L(1).                                    (13)

No successor Gaussian approximation occurs at any stage. The only Gaussian
law is the initial source.

The sharp deficit theorem implies the convenient rational closure

    alpha-e_j >= g_j^2/(8alpha+4g_j).                      (14)

Consequently every state (a_j,b_j,G_j) gives the scalar certificate

    [(1+b_j)alpha-a_j][alpha+G_j/2] >= G_j^2/8-o_L(1).     (15)

Let C(a,b,G) be the positive root of

    [(1+b)x-a][x+G/2]=G^2/8.                               (16)

Then

    liminf alpha_n >= C(a_j,b_j,G_j)                       (17)

for every fixed depth j. Same-order regularization removes the initial
operator cap after the n-limit exactly as in the preceding source proofs.

Thus arbitrary finite actual-update depth is reduced to the scalar
recursion (11)--(12).

## 5. Every finite depth can be improved strictly

This is the key structural consequence.

Suppose a valid state certificate has G>0 and let

    x=C(a,b,G),
    X=(1+b)x-a=G^2/[4(2x+G)] > 0.                          (18)

Append one more actual partial update with probability v. Dropping only
nonnegative terms, the new mean certificate is

    a'=(1-v^2)a+v(1-v)G,
    b'=(1-v^2)b+v^2.                                      (19)

At alpha=x its slack is exactly

    (1+b')x-a'
      = X-Gv+(2x+G-X)v^2.                                 (20)

Put K=2x+G-X. In the present chain G<=1/sqrt(3) and x>1/pi, so X<2x,
hence K>G and

    v_*=G/(2K)                                             (21)

lies in (0,1/2). At this value,

    (1+b')x-a'
      =X-G^2/(4K)
      <X-G^2/[4(2x+G)]
      =0.                                                   (22)

Therefore the new MEAN bound already satisfies

    a'/(1+b') > x.                                         (23)

The tie-resolved count theorem gives M'=v_* M>0, so the next successor
also has G'>0. Applying closure (16) only strengthens (23).

Hence the construction can be repeated indefinitely: every finite-depth
certificate has a strictly stronger finite-depth successor certificate.

## 6. Strict improvement beyond the current B_2

The existing two-update constant B_2 is the supremum of the two-update
MEAN certificate T(t,u,v). Its parameter rectangle is compact and the
maximum is attained. A maximizing triple has u>0 and v>0:

- u=0 gives at most e(t)<=1/pi;
- v=0 reduces to the pre-successor mean bound, strictly below the already
  proved one-successor closure B_opt;
- B_2>B_opt.

Therefore at a B_2-maximizing triple, the tie-resolved first count gives

    M_1=u H(t)>0,

and the SECOND actual successor has the strictly positive field floor

    G_2=2v d_v [u H(t)]^(3/2)>0.                            (24)

Apply closure (16) to that second successor. Since its mean bound is B_2
and G_2>0, the closure root is strictly larger than B_2. Thus

    boxed: liminf alpha_n >= B_tie > B_2.                  (25)

Now apply Section 5 repeatedly from that same certified state. It produces
an explicit infinite sequence

    B_2 < B_tie = x_2 < x_3 < x_4 < ... < 1/2,             (26)

where every x_j is an unconditional liminf lower bound coming from a
FINITE number of actual Boolean updates.

The sequence is increasing and bounded, so

    B_infty := sup_j x_j = lim_j x_j

exists, and because liminf alpha_n>=x_j for every fixed j,

    boxed: liminf alpha_n >= B_infty > B_2.                (27)

No interchange of n and an infinite update trajectory is used: each x_j
is first proved at fixed finite depth, and only the scalar certified
constants are then supremized.

## 7. What this does and does not solve

This removes the previous finite-depth wall in the lower-bound dynamics.
The actual non-Gaussian successor mechanism can be propagated to arbitrary
fixed depth, and there is no finite stage at which the certified lower
bound saturates exactly.

It does NOT prove that alpha_n converges, provide a cross-order comparison,
or show that B_infty is near the true liminf. The strict increments become
very small under the elementary recursion. What is new is the theorem that
the repeated-update mechanism itself supports an infinite strict hierarchy
of unconditional lower bounds rather than only one or two isolated updates.
