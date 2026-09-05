# Actual-diagonal stability of the weighted two-trace upper

2026-09-05. Analytic all-shell stability theorem. No mathematical
computation or numerical scalar evaluation was run. This is separate
from the frozen weighted-field theorem and does not alter it.

The comparison retains the ACTUAL positive field covariance and the
ACTUAL normalized source contraction. It does not replace the source
contraction by K divided by a scalar, which need not be a contraction
when a small set of diagonal weights is large.

## 1. Setup and the diagonal-dispersion parameter

Use the complete source, positive diagonal D, weighted matrices, and
Gaussian coefficients of
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
Thus N=2n, D+-K>=0, k,v>=0, w=k+v<=1, and

    S=tr D,       dbar=S/N,
    T=D^(-1/2)KD^(-1/2),                 ||T||op<=1,
    L_D=D^(-1/2)H_BD^(-1/2),             ||L_D||op<=1.

The diagonal in this theorem need not be trace-optimal. Define its
actual arithmetic-harmonic dispersion by

    delta=S tr(D^(-1))/N^2-1 >=0.                        (1.1)

For any actual representative theta=(p_D,q_D,c_D), chosen WITHIN
its final original/weighted refined cell, let M=M_theta
be the positive base field covariance of the cited theorem. In
particular

    M>=0,       diag M=wn,       tr M=wN^2/2.             (1.2)

The original energies p,q,c always remain those of that representative.
For |eta|<1 put

    a=|eta|,        b=1-a>0,       F_eta=(I-eta L_D)^(-1),
    P_eta=D-eta H_B,
    T_eta=tr(MP_eta^(-1)),
    R_eta=tr(DP_eta^(-1)MP_eta^(-1)),
    u=2c/S,        u_D=2c_D/N.                           (1.3)

Both u and u_D lie in [-1,1]. The first uses D+-H_B>=0; the second
uses contraction of L_D. The theorem below applies to all original
shells, not only p=q=0.

## 2. Uniform energy compatibility from this actual dispersion

Write t_i=d_i/dbar, so `N^(-1)sum t_i=1` and

    N^(-1)sum (t_i-1)^2/t_i=delta.

The pointwise inequalities

    (sqrt(t)-1)^2 <= (t-1)^2/t,
    (1/sqrt(t)-1)^2 <= (t-1)^2/t

therefore give

    sum_i (sqrt(d_i)-sqrt(dbar))^2 <= S delta,
    sum_i (sqrt(dbar/d_i)-1)^2 <= N delta.                (2.1)

No maximum-diagonal bound is used. Put E=K-dbar T. Since
`K=D^(1/2)T D^(1/2)` and ||T||op<=1, for every full Boolean z,

    |z^T E z|
      <= ||(D^(1/2)-sqrt(dbar)I)z||
             (||D^(1/2)z||+sqrt(dbar)||z||)
      <=2S sqrt(delta).

Consequently

                            Phi(E)<=S sqrt(delta).      (2.2)

The matrix E is symmetric with zero diagonal. Its principal-block
quadratic forms on a cube have absolute value at most 2Phi(E):
extend the unused coordinates by zero and use multilinearity of a
zero-diagonal quadratic form to maximize on Boolean vertices.
Its cross block has Boolean bilinear norm at most Phi(E), because
flipping the second Boolean block gives quadratic half-values
`s+c_E` and `s-c_E`, whose maximum absolute value is at least |c_E|.
Thus (2.2) proves, uniformly over ACTUAL states,

    |p_D-p/dbar| <=2N sqrt(delta),
    |q_D-q/dbar| <=2N sqrt(delta),
    |c_D-c/dbar| <= N sqrt(delta),
    |u_D-u|      <=2 sqrt(delta).                        (2.3)

These are quantitative compatibility inequalities, not an assertion
that the original and weighted shell triples coincide.

## 3. Trace-norm stability of the Gaussian covariance congruence

Let V_0=sqrt(dbar)D^(-1/2), and define

    E_delta=(wN^2/2)sqrt(delta)(1+sqrt(1+delta)).

Then the nuclear norm obeys

                     ||V_0 M V_0-M||_1<=E_delta.         (3.1)

Here ||.||_1 denotes the sum of singular values, not an entrywise norm.
To prove this, factor through the positive square root M^(1/2).
The constant diagonal (1.2) and (2.1) give the exact or bounded
Frobenius identities

    ||M^(1/2)||_F^2=wN^2/2,
    ||V_0M^(1/2)||_F^2=(wN^2/2)(1+delta),
    ||(V_0-I)M^(1/2)||_F^2<=(wN^2/2)delta.

Expand
`V_0 M V_0-M=(V_0-I)M V_0+M(V_0-I)` and apply the product
inequality `||AB||_1<=||A||_F||B||_F` to each term. This proves
(3.1) without commutation, a bound on ||D||, or a rank restriction.

The natural-D field traces satisfy

    dbar T_eta=tr(V_0 M V_0 F_eta),
    dbar R_eta=tr(V_0 M V_0 F_eta^2).

Define scalar-I REFERENCE traces using the same actual M and L_D:

    bar T_eta=tr(MF_eta),    bar R_eta=tr(MF_eta^2).

Since ||F_eta||op<=1/b, (3.1) implies

    |dbar T_eta-bar T_eta|<=E_delta/b,
    |dbar R_eta-bar R_eta|<=E_delta/b^2.                  (3.2)

For the combined trace one can retain its cancellation. The spectral
function `F_eta-bF_eta^2` is positive and has norm at most 1/(4b).
Indeed, with t=1-eta lambda>=b, its eigenvalue is
`(t-b)/t^2`, whose maximum over t>=b is 1/(4b). Therefore

    |dbar(T_eta-bR_eta)-(bar T_eta-b bar R_eta)|
                                           <=E_delta/(4b),
    0<=bar T_eta-b bar R_eta<=wN^2/(8b).                 (3.3)

The same combined natural-D trace is nonnegative by its exact
ellipsoid representation. Thus all square roots below are real.

## 4. Stability of the complete two-trace expression

Write the VALID natural-D upper expression at fixed eta as

    B_D(theta,c;eta)=sqrt(S) {
       sqrt((a-eta u)(T_eta-bR_eta))
                        +sqrt(kappa)b sqrt(R_eta)},
    kappa=2/pi.

Define the numerical scalar-I reference expression

    B_flat(theta;eta)=sqrt(N) {
       sqrt((a-eta u_D)(bar T_eta-b bar R_eta))
                        +sqrt(kappa)b sqrt(bar R_eta)}. (4.1)

The reference still uses the actual M_theta and actual contraction
L_D. It is not the Gaussian law obtained by replacing those matrices
with unweighted scalar multiples.

For 0<=delta<=1 and every fixed |eta|<1,

    |B_D(theta,c;eta)-B_flat(theta;eta)|
          <=3sqrt(w)N^(3/2)delta^(1/4)/sqrt(1-|eta|).   (4.2)

The w=0 or delta=0 cases follow directly. For the remaining cases,
rewrite B_D with outer factor sqrt(N), replacing its two traces by
dbar times the natural-D traces. Each radius coefficient lies in
[0,2], and (2.3) bounds their difference by 2a sqrt(delta).
Equations (3.3) give

    |(a-eta u)dbar(T_eta-bR_eta)
               -(a-eta u_D)(bar T_eta-b bar R_eta)|
      <=E_delta/(2b)+a wN^2 sqrt(delta)/(4b)
      <=wN^2 sqrt(delta)/b.                             (4.3)

For the final inequality use
`1+sqrt(1+delta)+a<=2+sqrt(2)<4`.
The elementary inequality `|sqrt(x)-sqrt(y)|<=sqrt(|x-y|)`
for nonnegative x,y bounds the first square-root difference by
`sqrt(w)N delta^(1/4)/sqrt(b)` before the outer factor sqrt(N).

For the second term (3.2) gives

    sqrt(kappa)b |sqrt(dbar R_eta)-sqrt(bar R_eta)|
                                      <=sqrt(kappa E_delta).

Because delta<=1, kappa<1, and b<=1, the sum of these two bounds
is at most `3sqrt(w)N delta^(1/4)/sqrt(b)` before the outer factor.
This proves (4.2). The proof covers arbitrary attained p_D,q_D,c_D;
no zero-internal-energy cancellation is used.

## 5. Correct use on real-valued cells and in limits

For a cell of the frozen weighted-field theorem, its representative
theta has the same ORIGINAL c as every cell state. Hence B_D is a
genuine field-width upper at each eta. In contrast, the representative
c_D need not equal the weighted cross value at every state in the
cell. Accordingly B_flat in (4.1) is only a COMPARED NUMERICAL
FUNCTIONAL for that cell, not a claim of an exact scalar-I ellipsoid
constraint on the entire cell. Equation (4.2) is what justifies using
it as a proxy in an upper bound.

Fix any 0<b_0<=1 and restrict to `|eta|<=1-b_0`. Applying (4.2)
uniformly on this compact interval, and then the already proved
cell/noise comparison, gives

    E max_z |a_z+X_z|
      <= max_(actual cells j) [|a_j|
              +inf_(|eta|<=1-b_0) B_flat(theta_j;eta)]
         +3sqrt(w)N^(3/2)delta^(1/4)/sqrt(b_0)
         +6sqrt(k n log(2))
         +n sqrt(2(2k+v)log(2m)),                       (5.1)

where `m<=(2n^2+1)^6`, exactly as in the frozen field theorem.
The original drift is still `(p_j-q_j)/2+s c_j`.
Add its recorded O_(C_D)(n) padding cost if the full Gaussian
covariance is retained instead of the positive base covariance.

For delta_N tending to zero, at every fixed b_0>0 the added metric
comparison cost is o(N^(3/2)). One may subsequently study endpoint
limits of the actual reference traces. Equation (4.2) does not give
a uniform error as |eta| tends to one with N, and no unrestricted
endpoint-infimum interchange is asserted in (5.1).

## 6. Connection to a small canonical gap and remaining scope

A separate residual/commutator theorem for a trace-optimal diagonal
of the literal complete K, source
`original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`,
establishes the bound

    g=S-tr|K|^3/(N-1)>=0,
    delta<=4Sg/[(N-1)N^2].                              (6.1)

That implication uses the actual optimal diagonal and canonical SDP
residual; it is not assumed merely from D+-K>=0 in Sections 1--5.
Under a fixed trace cap S=O(N^(3/2)), the additional actual hypothesis
g=o(N^(3/2)) implies delta=o(1). Thus (5.1) gives a rigorously stable
all-shell numerical trace comparison in that source-compatible range.
The residual theorem supplies (6.1); the present theorem supplies the
general resolvent and metric stability that follows from delta.

No assertion is made that conditional original minimizers satisfy this
small-gap hypothesis. Nor does the note evaluate the supremum in
(5.1), even in the small-gap range. It identifies a valid scalar-I
metric comparison while retaining the actual weighted source matrices;
their trace evaluation and the complementary large-gap case remain
separate requirements. Original MO convergence is not claimed.
