# Hypercontractive row/column bridge without the one-half loss

2026-11-04. Direct progress on the current ProofAtlas/finite-temperature
bridge target. This note proves the row/column double-counting inequality
for independently noised one-point tests from a common center. It does not
by itself identify those tests with the full MOP-DINF-CLOSURE kernel.

Let E=I x J be the edge set of a rectangular sign matrix and let
X=(X_ij)_(i,j in E) be independent uniform Rademachers.

For every row i let f_i:{+-1}^J -> [0,1], and for every column j let
g_j:{+-1}^I -> [0,1]. Let T_rho be the usual Boolean noise operator,
0<=rho<=1, acting on the coordinates of the relevant row or column.

Then

    boxed:
    E_X [ prod_i (T_rho f_i)(X_i,*)
          prod_j (T_rho g_j)(X_*,j) ]
      <=
      prod_i ||f_i||_(1+rho^2)
      prod_j ||g_j||_(1+rho^2).                            (1)

In particular, if f_i=1_(A_i), g_j=1_(B_j), with
p_i=P(A_i), q_j=P(B_j), then

    boxed:
    E_X [ prod_i T_rho 1_(A_i)(X_i,*)
          prod_j T_rho 1_(B_j)(X_*,j) ]
      <=
      ( prod_i p_i prod_j q_j )^(1/(1+rho^2)).             (2)

Equivalently, writing the row and column one-point budgets

    B_row=-sum_i log p_i,
    B_col=-sum_j log q_j,

the joint same-center smoothed test has negative log probability at least

    boxed:
    (B_row+B_col)/(1+rho^2).                               (3)

Thus the generic row/column Finner loss of one half is replaced, after a
rho-noise layer, by the factor 1/(1+rho^2). The unrecovered bridge budget is

    boxed:
    [rho^2/(1+rho^2)] (B_row+B_col).                       (4)

Hence whenever the total one-point budget is O(L_N), choosing rho_N with

    rho_N^2 L_N=o(R_N)

makes the double-counting loss o(R_N). In particular, any polylogarithmic
excess in L_N can be absorbed by letting rho_N tend to zero at a matching
inverse-polylogarithmic rate, provided the surrounding MOP Markov-kernel
identification preserves the same common temperature set.

## Proof

Put

    F_i=T_rho f_i,  G_j=T_rho g_j.

Every edge coordinate X_ij occurs in exactly two factors: the row factor
F_i and the column factor G_j. Finner's inequality for the bipartite
edge-product space, with fractional cover weight 1/2 on every row and every
column factor, gives

    E prod_i F_i prod_j G_j
      <= prod_i ||F_i||_2 prod_j ||G_j||_2.                (5)

Bonami--Beckner hypercontractivity at q=2 and
p=1+rho^2 gives exactly

    ||T_rho h||_2 <= ||h||_(1+rho^2).                      (6)

Applying (6) to every factor in (5) proves (1).

For an indicator 1_A,

    ||1_A||_(1+rho^2)=P(A)^(1/(1+rho^2)),

which proves (2), and taking minus logarithms gives (3)--(4).

## Markov interpretation

For each row and column factor, T_rho f(x) is the success probability of a
rho-correlated BSC observation started from the common center x. Thus the
left side of (1) is the probability that all row and column one-point tests
succeed when their observation noises are conditionally independent given
that common center.

This is precisely the regime in which generic checkerboard/Finner accounting
would charge half of each bridge budget. The extra BSC layer plus
hypercontractivity recovers asymptotically all of both budgets when
rho->0.

## What remains before this closes MOP-DINF

The remaining identification is theorem-shaped and narrow:

1. represent the source-reported row and column one-point tests in
   MOP-DINF-CLOSURE as conditionally independent rho_N-noised observations
   from one common deepest output (or dominate them pointwise by such tests);
2. keep both dilation factors r=2,3 on one common unbounded temperature set;
3. choose rho_N so the loss (4), plus the Markov synchronization cost, is
   sublinear on the logarithmic scale required by the source reduction.

No generic row/column entropy inequality can improve the 1/2 exponent
without extra structure. Equation (1) shows that BSC smoothing is enough
structure to replace 1/2 by 1/(1+rho^2).
