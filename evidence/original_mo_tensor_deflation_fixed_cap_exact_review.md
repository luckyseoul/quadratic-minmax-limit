# Independent exact-worker review: fixed-cap tensor-deflation rate

2026-09-05. Status: PASS, complete analytic review. No mathematical
computation, signing construction run, numerical experiment, or source
edit was performed by this reviewer.

Reviewed the complete 342-line source
`/tmp/original_mo_tensor_deflation_fixed_cap_rate.md`, SHA256
`22febfa722afb3e18878f23f8e140895da90a3eb41fe0179356b08232d44f27a`.
The proof worker is the source author; this receipt is an independent
full-source review by the exact worker. No corrections are requested.

## 1. Spectral repair estimates

Both scalar positive-part estimates are correctly normalized:
the maxima of `(t-1)/t^2` and `(t-1)/t^3` are `1/4` and `4/27`.
Counting ordered eigenvalue products gives the three rank bounds;
using `S_2=n(n-1)` gives the displayed trace bounds. The operator-norm
bound uses the maximal absolute eigenvalue product and is valid for
either tensor sign. Spectral calculus also gives
`R <= (A^2 tensor A^2)/(4K^2 n^2)`, and hence the stated coordinate
variance bound because `diag(A^2)=n-1`.

The quoted norm-only cubic estimate is used only for the ancillary
`O_C(n^2/K^3)` trace comparison. The main construction and rate
obstruction in Sections 3--7 do not depend on that quoted estimate.

## 2. Finite sign template and localized overlap

For `H=F tensor^(2t)`, the order, diagonal, regular eigenvector,
square identity, and positive-eigenspace dimension are all correct.
Replacing the principal S block by all ones preserves symmetry and
every diagonal sign. The identity
`C_k=J_S+(H-P_S H P_S)` gives perturbation norm at most `2 sqrt(k)`.
The Rayleigh quotient at `v=1_S/sqrt(s)` is exactly s, not an
approximation. Projecting an actual top eigenvector away from v
therefore proves the overlap lower bound `sqrt(3s)/2` for q>=4.

The subspace `V=ker(H-sqrt(k)I) intersect {z:z_S=0}` is genuinely
untouched by the block perturbation. Its dimension is at least k/4.
Its projection annihilates `1_S` and also the top eigenvector u,
since the corresponding eigenvalues differ. The inequality
`sum sqrt(Q_jj)>=tr(Q)` needs no uniform diagonal hypothesis.

## 3. Actual sign source and uniform cap

`A=H_m tensor C_k-I` is an actual complete symmetric zero-diagonal
signing. The cap does not illegitimately insert the large template
operator norm: it bounds `tau(C_k)` by `tau(H)+2s^2`.
The grouped vectors `x_b/sqrt(m)` and `H_m y_c/m` are both unit
vectors, establishing the exact amplification inequality
`beta(H_m tensor C_k)<=m^(3/2) tau(C_k)`.
Consequently `Phi(A)<=(1/2+a^2)n^(3/2)+n/2` as claimed.

For every fixed C>1/2 a positive dyadic amplitude can be selected
independently of K with `1/2+a^2<C`. Holding k, and therefore K,
fixed while increasing m absorbs the identity term. Thus the cap
holds for arbitrarily large constructed n at each selected K.

## 4. Spectral subcovariances and Gaussian norm direction

The three displayed projections are actual mutually appropriate
eigenspace projections of A. Their eigenvalues include the subtracted
identity exactly. The lower products in (5.2) are at least `9/4`,
so the positive-part covariance contains the indicated unit
high-by-bulk product projection for each tensor sign.
This uses commuting spectral projections, not operator monotonicity
of a positive-part compression.

For `G_0=P Z Q_epsilon`, the covariance is the indicated ordered
product with the consistent row-vectorization convention. The
Loewner comparison permits an independent centered Gaussian
increment. Conditional Jensen for the convex norm beta gives the
required LOWER bound on the larger covariance.
The sparse vector `x_0=1_m tensor 1_S` is legitimately admissible
in the cube formulation. Its projected norm and the sum of square
root bulk diagonal variances give (6.2), with
`s/k=a^2/q=a^2/(4K)`. The exact lower constant is
`a sqrt(6/pi)/32`, and weakening it to `a/32` is valid.

## 5. Symmetric and diagonal-free version

The symmetric restriction contains both high-by-bulk product spaces.
`(PZQ+QZ^T P)/sqrt(2)` has unit Gaussian coefficients in the
orthonormal symmetric basis, so its covariance normalization is
correct. Since `Q_epsilon x_0=0`, the fixed-vector argument incurs
exactly the stated factor `1/sqrt(2)`.

Diagonal covariance coordinates are unchanged by restricting to the
symmetric space. Their expected absolute sum is at most
`n/(K sqrt(2pi))`. The norm triangle inequality gives the correct
subtraction direction when removing the diagonal. Cube polarization
for a symmetric zero-diagonal matrix gives `beta<=4 Phi`, proving
(7.2). Its subtracted term is lower order for each fixed K.

## 6. Quantifiers and limits of the conclusion

Dividing the constructed lower bound by `n^(3/2)` rules out every
uniform fixed-cap rate with exponent greater than 1/2 along the
unbounded threshold sequence. The order is fixed amplitude, then
K, then arbitrarily large m and n; the cap constant is independent
of K. This also gives the stated asymptotic-in-n obstruction when
its estimate is uniform over the indicated sources and thresholds.

The source correctly makes no claim about C=1/2, exact original
minimizers, slack coupled adaptively to K, a directly analyzed clipped
law, or the full cross-coupled operator. It also does not identify
the repair norm with the difference of two Gaussian-law norms or
invoke a Gaussian covariance for an indefinite unrepaired matrix.
