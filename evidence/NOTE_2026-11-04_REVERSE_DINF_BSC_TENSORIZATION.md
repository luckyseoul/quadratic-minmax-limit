# Reverse-D-infinity same-center BSC tensorization

2026-11-04. This is the reverse-likelihood normalization of the same-center
BSC synchronization identity. In this orientation the synchronization term
is FAVORABLE: k bridge budgets are recovered with no loss.

Let U be the uniform law on {+-1}, and let

    K_beta(y|x)=exp(beta*x*y)/(2 cosh beta)

be the binary symmetric channel in exponential coordinates. Define the
reverse one-point likelihood ratio

    R_beta(y|x)
      := U(y)/K_beta(y|x)
       = cosh(beta) exp(-beta*x*y).                         (1)

This is the pointwise quantity entering first-argument / reverse
D_infinity.

For every integer k>=2,

    boxed:
    R_(k beta)(y|x)
      = d_(k,beta) R_beta(y|x)^k,                           (2)

where

    d_(k,beta)=cosh(k beta)/cosh(beta)^k >= 1.              (3)

Consequently

    boxed:
    log R_(k beta)(y|x)
      >= k log R_beta(y|x).                                (4)

## Proof

Directly,

    R_beta(y|x)^k
      =cosh(beta)^k exp(-k beta xy),

while

    R_(k beta)(y|x)
      =cosh(k beta) exp(-k beta xy).

This proves (2). To prove (3), let S_1,...,S_k be independent Rademachers.
Then

    cosh(beta)^k
      = E exp(beta sum_i S_i)
      <= cosh(k beta),

because |sum_i S_i|<=k and cosh is increasing on [0,infinity);
equivalently the standard multiple-angle polynomial has nonnegative
coefficients in sinh(beta)^2. Hence (4).

## Product-coordinate form

For a finite coordinate set E and common center x, common output y,

    R_beta^E(y|x)
      =prod_(e in E) R_beta(y_e|x_e).

Equation (2) tensorizes exactly:

    boxed:
    log R_(k beta)^E(y|x)
      =
      k log R_beta^E(y|x)
      + |E| [log cosh(k beta)-k log cosh beta]
      >= k log R_beta^E(y|x).                              (5)

Thus if k row/column one-point reverse-information tests use the SAME center
and the SAME output coordinate on every shared matrix entry, their k bridge
budgets are dominated pointwise by one deepest-output reverse budget at
temperature k beta.

There is no double-counting loss.

## Row/column matrix version

For a rectangular matrix E=I x J, the product of all row one-point reverse
likelihood ratios and all column one-point reverse likelihood ratios counts
each shared entry exactly twice. Therefore, under exact same-center and
same-output synchronization,

    boxed:
    B_cross(2 beta)
      >= B_rows(beta)+B_cols(beta),                         (6)

where each B is the corresponding pointwise reverse log-likelihood budget.

For a k-fold synchronized construction the identical statement is

    boxed:
    B_cross(k beta)
      >= sum_(r=1)^k B_r(beta).                            (7)

This is precisely the inequality shape requested by the current public
ProofAtlas description of MOP-DINF-CLOSURE: cross reverse max-information
must dominate k bridge budgets.

## Common temperature set

No n-dependent temperature choice is needed. Equation (7) holds for EVERY
beta>=0. Hence one may choose any unbounded set S that is closed under
multiplication by 2 and 3, for example

    S={2^a 3^b beta_0 : a,b>=0}.                            (8)

The k=2 and k=3 inequalities then live on the same unbounded temperature
semigroup, with zero synchronization error.

## Scope / exact remaining identification

The public ProofAtlas page states that at each fixed temperature exact bridge
and cross reverse-max-information identities reduce the two- and
three-dilation bounds to a same-center inequality MOP-DINF-CLOSURE, and that
the current priority is to combine row and column one-point tests without
double-counting shared entries.

Equations (6)--(8) close that algebraic double-counting step PROVIDED the
source's one-point bridge factors are exactly the reverse BSC ratios (1), or
factor pointwise into them, and PROVIDED 'same-center' also identifies the
same output bit on each shared entry.

Those source definitions are not reproduced on the public page. Therefore
this note does not yet register a completion proof in original_mo_status.py.
The remaining task is an identification audit against the exact hidden/source
MOP-DINF definitions, not a new analytic inequality.
