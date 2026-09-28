# Hybrid Komlos--hypergeometric one-vertex extension

2026-09-28. Direct one-vertex extension theorem under the scope freeze.
This removes the odd-signature-class penalty entirely by balancing the
class means with Komlos and leaving only within-class impurity for the
probabilistic tail. The original convergence problem remains OPEN.

Let A be a complete symmetric zero-diagonal signing of order n, put
F=Phi(A), fix r>0, and let

    H={s : D_s+r<=n}

be the active antipodal stable family.

Choose a core C subset H and form the projective core-signature classes
G exactly as in the balanced-core notes. In each class choose canonical
lambda_i in {+-1} so that lambda_i x_i^c is constant over i in G for
every c in C.

For every active state s and class G define

    y_i^s=lambda_i x_i^s,

    mu_(s,G)=|G|^(-1) sum_(i in G) y_i^s,                  (1)

    q_(s,G)
      =min( #{i in G:y_i^s=+1},
            #{i in G:y_i^s=-1} ),                          (2)

and

    q_s=sum_G q_(s,G).                                     (3)

For core states c in C, every class is pure:

    q_(c,G)=0,   |mu_(c,G)|=1.                             (4)

Let O be the set of odd-cardinality signature classes. Put

    K=3 sqrt(2 pi),

the explicit Guo--Fang--Lu Komlos constant.

Choose margins h_s satisfying

    0<h_s<D_s+r.                                           (5)

## Theorem

Assume that every odd signature class G obeys

    boxed:
    K^2 sum_(s in H) mu_(s,G)^2 / h_s^2 <= 1.              (6)

Also assume

    boxed:
    sum_(s in H:q_s>0)
      2 exp( -(D_s+r-h_s)^2/(8q_s) ) < 1.                 (7)

Then there exists an actual incident sign row a in {+-1}^n such that

    boxed:
    Phi(A extended by a) < F+r.                            (8)

No bound on the NUMBER of odd signature classes is required.

## Proof

For an even class G choose a uniformly random balanced sign vector
b^G in {+-1}^G with

    sum_(i in G) b_i^G=0.                                  (9)

For an odd class G we will first choose a class imbalance sign
eps_G in {+-1}; conditional on eps_G, choose b^G uniformly among sign
vectors satisfying

    sum_(i in G) b_i^G=eps_G.                              (10)

Finally set

    a_i=lambda_i b_i^G,   i in G.                          (11)

### Step 1: choose all odd-class imbalance signs by Komlos

For each odd class G form the vector indexed by active states

    v_G(s)=K mu_(s,G)/h_s.                                 (12)

Condition (6) says ||v_G||_2<=1 for every G. The Guo--Fang--Lu Komlos
theorem therefore gives signs eps_G such that simultaneously for every
active state s,

    |sum_(G in O) eps_G mu_(s,G)| < h_s.                   (13)

Fix these signs from now on.

### Step 2: center the internal class fluctuations

For a class G and state s define

    X_(s,G)=sum_(i in G) b_i^G y_i^s.                      (14)

If G is even, E X_(s,G)=0. If G is odd, symmetry under the uniform slice
with sum b_i=eps_G gives

    E X_(s,G)=eps_G mu_(s,G).                              (15)

Orient y^s inside G so that its minority set M has size
q=q_(s,G). Since sum_i b_i is fixed,

    X_(s,G)-E X_(s,G)
      = -2[ sum_(i in M)b_i^G - E sum_(i in M)b_i^G ].     (16)

The q coordinates in M are sampled without replacement from a fixed
+-1 population. Hoeffding's comparison theorem therefore gives the mgf

    E exp(theta[X_(s,G)-E X_(s,G)])
      <= exp(4q theta^2/2).                                (17)

Thus the centered class fluctuation is subgaussian with variance proxy 4q.

The class choices are independent after the eps_G have been fixed, so the
total centered fluctuation

    Z_s=sum_G (X_(s,G)-E X_(s,G))                           (18)

is subgaussian with variance proxy

    4q_s.                                                   (19)

By (13), the total deterministic mean obeys

    |E sum_G X_(s,G)| < h_s.                               (20)

Hence

    P( |a.x^s| >= D_s+r )
      <=
      2 exp( -(D_s+r-h_s)^2/(8q_s) )                       (21)

when q_s>0.

If q_s=0, every class is pure for s, there is no centered fluctuation at
all, and (13) alone gives

    |a.x^s|<h_s<D_s+r.                                     (22)

Condition (7) and the union bound now give positive probability that every
active state satisfies

    |a.x^s|<D_s+r.                                         (23)

Inactive stable states are automatically harmless because D_s+r>n.
The exact stable-state extension identity proves (8). QED.

## Why this is stronger than the random-leftover theorem

The random-leftover theorem paid every odd signature class as one full unit
of variance:

    v_s=u+4q_s.

Here the odd classes are not discarded and are not independently charged.
Their deterministic class means are balanced ALL AT ONCE by Komlos.
After that step, only the genuine within-class impurity remains:

    boxed: variance proxy = 4q_s.                           (24)

Thus the odd-class count u disappears completely.

The price is the column condition (6), which depends on the weighted mean
profiles mu_(s,G). Pure states have |mu|=1 but zero fluctuation; highly
impure states have smaller |mu| and are cheaper in the Komlos step.

## Linear-size core corollary

For a core state c, |mu_(c,G)|=1 on every odd class. If we momentarily
ignore residual contributions to (6), a common core margin h gives the
necessary core-side budget

    K^2 |C|/h^2 <= 1.                                      (25)

So with h on the neutral sqrt(n) scale, the construction can support a
core of LINEAR size in n, independent of how many odd signature classes
that core creates.

More generally, condition (6) is exactly

    sum_(c in C) 1/h_c^2
      + sum_(s in H\C) mu_(s,G)^2/h_s^2
      <= 1/(18pi)                                          (26)

for every odd class G.

This is the quantitative core-size/impurity tradeoff missing from the
previous notes: large cores are now allowed when their margins are large,
and the residual family is charged according to class purity rather than
cardinality.

## Neutral consequence

For an exact order-n minimizer take

    r_n=m_n[(1+1/n)^(3/2)-1].

If some core C and margins h_s satisfy (5)--(7) with r=r_n, then

    boxed: alpha_(n+1)<alpha_n.                            (27)

More generally, replacing r_n by r_n+rho_n and obtaining such a
certificate for all large n with

    sum_n rho_n/(n+1)^(3/2)<infinity

proves convergence through the frozen one-vertex route.

## Finite discriminator

For a candidate exact minimizer, the data to compute are now:

1. active stable classes and deficits D_s;
2. a candidate core C;
3. signature classes G;
4. mu_(s,G) and q_s;
5. margins h_s.

Conditions (6) and (7) are a complete extension certificate. In
particular, there is no longer any reason to reject a core merely because
it creates many odd signature classes.
