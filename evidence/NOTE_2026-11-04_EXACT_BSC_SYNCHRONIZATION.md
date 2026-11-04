# Exact same-center BSC synchronization for k bridge budgets

2026-11-04. Direct theorem for the current ProofAtlas MOP-DINF bridge target.
It strengthens NOTE_2026-11-04_HYPERCONTRACTIVE_ROW_COLUMN_BRIDGE when the
one-point tests share the same center AND the same output coordinate.

Write the binary symmetric channel in exponential form

    K_beta(y|x)=exp(beta*x*y)/(2 cosh beta),
    x,y in {+-1}, beta>=0.

For every integer k>=2,

    boxed:
    K_beta(y|x)^k = c_(k,beta) K_(k beta)(y|x),             (1)

where

    c_(k,beta)
      =2^(1-k) cosh(k beta)/cosh(beta)^k.                   (2)

The synchronization loss is

    s_k(beta):=-log c_(k,beta)
      =k log(1+e^(-2 beta))
        -log(1+e^(-2k beta))                               (3)

and therefore

    boxed:
    0<=s_k(beta)<=k e^(-2 beta).                           (4)

## Proof

Equation (1) is immediate:

    K_beta(y|x)^k
      =exp(k beta xy)/(2^k cosh(beta)^k)

and

    K_(k beta)(y|x)
      =exp(k beta xy)/(2 cosh(k beta)).

Their ratio is (2).

Using

    cosh beta=(e^beta/2)(1+e^(-2 beta))

in (2), every exponential and power-of-two factor cancels, giving

    c_(k,beta)
      =(1+e^(-2k beta))/(1+e^(-2 beta))^k.

Taking minus logarithms gives (3). Positivity follows from
(1+u)^k>=1+u^k for 0<=u<=1. Finally
log(1+u)<=u gives (4).

## Matrix form

Let E be ANY finite coordinate set and suppose k one-point likelihood tests
use the same center x in {+-1}^E and the same output y in {+-1}^E.
Then

    prod_(r=1)^k prod_(e in E) K_beta(y_e|x_e)
      = c_(k,beta)^|E|
        prod_(e in E) K_(k beta)(y_e|x_e).                 (5)

If

    B_r=-log prod_e K_beta(y_e|x_e)

is the r-th bridge budget and

    B_deep=-log prod_e K_(k beta)(y_e|x_e),

then

    boxed:
    B_deep
      = sum_(r=1)^k B_r - |E| s_k(beta).                   (6)

Thus all k bridge budgets are recovered, not a 1/k fraction, with
synchronization cost at most

    boxed:
    k |E| e^(-2 beta).                                     (7)

## Common-temperature consequence for dilation 2 and 3

For a rectangular bridge with |E|=O(n^2), choose ONE common unbounded
temperature sequence

    beta_n=log n.                                           (8)

Then for both k=2 and k=3,

    k |E| e^(-2 beta_n)=O(1).                              (9)

More generally beta_n=(1/2+epsilon)log n gives
O(n^(1-2 epsilon)), already sublinear.

Therefore the specific double-counting obstruction stated in the current
ProofAtlas task -- recovering two/three same-center bridge budgets up to
sublinear logarithmic cost on one common unbounded temperature set -- is
removed whenever its one-point tests are exact BSC likelihoods of the form
above at a shared output.

## D-infinity orientation

For a point output y, reverse max-information uses a negative log likelihood
(or a likelihood ratio against a reference output law). Multiplying by a
reference factor independent of x merely adds the corresponding reference
normalization to both sides of (6). Thus the identity preserves the
FIRST-ARGUMENT/reverse orientation: the deepest-output surprise equals the
sum of the k same-center bridge surprises minus the explicit synchronization
penalty.

No total-variation approximation is used; (1) is pointwise and
multiplicative even for exponentially small output probabilities.

## Remaining audit before claiming convergence

ProofAtlas states that MOP-DINF-CLOSURE would imply the two dilation bounds
if the row/column one-point tests can be combined without double-counting.
Equation (6) proves that combination for exact same-center BSC kernels and
(9) supplies one common unbounded temperature sequence for k=2,3.

A completion claim still requires checking, against the exact source
definitions not reproduced on the public ProofAtlas page, that:

1. each bridge one-point factor in MOP-DINF-CLOSURE is exactly K_beta
   (or has a pointwise likelihood-ratio representation that cancels to it);
2. all k tests bind to the SAME output coordinate y, not merely the same
   orbit or center;
3. the source's native temperature parameter is beta in (1), or its
   reparameterization preserves beta_n=log n as an admissible unbounded
   common set;
4. no additional error outside the stated double-counting term has exponent
   >=1 after synchronization.

Subject to those source-identification checks, (6)--(9) close the priority
row/column bridge with O(1) synchronization cost for both dilation factors.
