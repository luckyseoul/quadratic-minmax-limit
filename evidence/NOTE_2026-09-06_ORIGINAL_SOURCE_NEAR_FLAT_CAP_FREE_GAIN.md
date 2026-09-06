# Cap-free original-source gain for the actual near-flat spectral law

2026-09-06. Analytic extension for ACTUAL complete symmetric signings.
No mathematical program, checker, solver, integral evaluation, spectral
scan, optimization run, or construction search was executed.

## 1. Statement and scope

Let A be an actual real symmetric zero-diagonal signing of order n,
with A_ij in {+1,-1} for i!=j. Write

    M=A/sqrt(n),       Q_M(x)=x^T M x/2,
    Phi(A)=max_(x in {+1,-1}^n)|x^T A x/2|,
    alpha(A)=Phi(A)/n^(3/2),
    kappa=2/pi,       rho=16/25,       lambda0=5/4.

Suppose only that the FULL empirical eigenvalue law of M converges to

    nu=(9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4)).        (1.1)

Then, without an actual operator-norm cap,

    liminf alpha(A)>=5kappa/8+16/3125
                   >2/5+3/1100.                           (1.2)

The constant is deliberately the same as in the frozen bounded
near-flat theorem. The cutoff C=5/3 below is an AUXILIARY bulk cutoff,
not a bound on ||M||op. Every candidate mean is evaluated in the full
original Q_M and independently rounded on the original n coordinates.
No source coordinate, outlying eigenspace, or actual objective is
discarded. This theorem concerns (1.1) only; it does not remove the
operator premise from the separate all-law theorem or settle global MO.

## 2. Actual second moment controls the spectral tails

Completeness gives exactly

    diag M=0,       (M^2)_ii=1-1/n,       tr M^2=n-1.       (2.1)

The second moment of nu is rho lambda0^2=1. For every fixed C>lambda0,
weak convergence at the two non-atomic cutoffs gives

    (1/n)tr[M^2 1_(|M|<=C)] ->1.

Subtracting from (2.1) proves the required second-moment tail control.
In particular, for the fixed C=5/3, put

    M_b=M 1_(|M|<=C),       M_o=M-M_b.

Then

    ||M_b||op<=C,       tr M_o^2=o(n),
    ||M_o||_1<=sqrt(n)||M_o||_F=o(n).                      (2.2)

Here ||.||_1 denotes the nuclear norm, not the entrywise norm. Weak
convergence alone would not imply (2.2); equality of the limiting
second moment with the exact complete-source second moment is essential.

Keep the FULL projectors, including any outlying eigenvectors:

    P_+=1_(M>lambda0/2),       P_-=1_(M<-lambda0/2),
    P=P_++P_-,       T=P_+-P_-.

Using a fixed bounded spectral interval and (2.2) on its complement
in the corresponding scalar error functions gives

    ||M/lambda0-T||_F=o(sqrt(n)),
    ||M^2/lambda0^2-P||_1=o(n).                            (2.3)

The second estimate is a TRACE-NORM estimate: squaring that spectral
error would require an unavailable fourth-moment tail bound. The
square of the first scalar error and the absolute second scalar error
are bounded by a constant times 1+x^2, so second-moment tails suffice
for both stated
norms. Their bounded-interval limits vanish at the three atoms.

## 3. Full positive-projector frame and actual phase baseline

Since diag M=0, (2.3) gives mean-square convergence of T_ii to zero.
The diagonal trace-norm inequality, (2.1), and (2.3) give

    (1/n)sum_i |P_ii-rho| ->0.

As 0<=P_ii<=1, this also gives mean-square convergence. Thus, with

    R0=(2/rho)P_+,       r_i=(R0)_ii,

we have ||R0||op<=2/rho and h_n=(1/n)sum_i(r_i-1)^2->0. Choose
eta_n=max(n^(-1/8),h_n^(1/4)) and I={i:|r_i-1|<=eta_n}.
Then |I^c|=o(n). Normalize R0[I,I] by its diagonal to define R[I,I],
put R[I^c,I^c]=I, and set the two off-blocks to zero. Consequently

    R>=0,       diag R=1,       ||R||op<=B_n=2/rho+o(1),
    ||R-R0||_F=o(sqrt(n)),       ||R-R0||_1=o(n).           (3.1)

Indeed deleting o(n) rows and columns of a bounded operator costs
o(sqrt(n)) in Frobenius norm; the inserted identity has the same cost.
The good diagonal rescaling differs from identity by O(eta_n).
No coordinates are deleted from A, M, or any Boolean variable.

Let G be centered Gaussian with covariance R and X_i=sign(G_i).
Every G_i is standard normal. Schur multiplication by a correlation
matrix preserves positive order and fixes I; hence for every q>=1,

    0<=R^(circ q)<=B_n I.

The normalized Hermite expansion of sign gives

    C_X=E[XX^T]=kappa R+C_tail,
    C_tail=sum_(odd q>=3)c_q^2 R^(circ q),
    sum_(odd q>=3)c_q^2=1-kappa,
    0<=C_tail<=(1-kappa)B_n I,       C_X<=B_n I.             (3.2)

All series converge in operator norm. Actual complete entries imply,
uniformly for odd q>=3,

    |tr(M R^(circ q))|
       <=n^(-1/2)sum_(i!=j)|R_ij|^q
       <=n^(-1/2)tr R^2<=B_n sqrt(n).                     (3.3)

By (3.1), ||M||_F=sqrt(n-1), and Cauchy--Schwarz,
tr[M(R-R0)]=o(n). First-moment convergence on P_+, which also follows
from (2.2), gives tr(M R0)/n->lambda0. Therefore

    E Q_M(X)/n ->kappa lambda0/2=5kappa/8.                 (3.4)

This is the actual positive original-source phase, not a bilinear
quantity or the energy of a restricted source.

## 4. Higher-chaos variance and first-chaos alignment without a cap

The projector identity P+T=rho R0 and (2.3), (3.1) give

    E_n=M^2/lambda0^2-rho R+M/lambda0,
    ||E_n||_1=o(n).                                       (4.1)

For odd q>=3, pair (4.1) with R^(circ q). Since
tr(R R^(circ q))=sum_(i,j)R_ij^(q+1)>=n, equations (3.3) and
||R^(circ q)||op<=B_n give, UNIFORMLY in q,

    tr(M^2 R^(circ q))>=lambda0^2 rho n-o(n)=n-o(n).

For the ACTUAL local field F=MX define

    F_1=sqrt(kappa)MG,
    v_i=E(F_i-F_(i,1))^2=(M C_tail M)_ii,
    mu_n=(1/n)sum_i v_i.

Summing the preceding uniform estimate with the nonnegative Hermite
weights, and using the exact squared complete-row lengths, yields

    mu_n>=1-kappa-o(1),
    0<=v_i<=(1-kappa)B_n(1-1/n)<=V+o(1),
    V=2(1-kappa)/rho.                                     (4.2)

Put a=lambda0 sqrt(kappa). The first-chaos alignment is

    (1/n)E||F_1-aG||^2
       =(kappa/n)tr[(M-lambda0 I)^2 R] ->0.                (4.3)

For R0 this follows directly from the positive limiting atom and the
second-moment tails. For Delta R=R-R0, split (M-lambda0 I)^2 on
|M|<=C and its complement. The first operator is bounded and
||Delta R||_1=o(n). The second is positive with trace o(n), by
(2.2) and rank 1_(|M|>C)=o(n), while ||Delta R||op=O(1).
Both pairings are o(n). This proves (4.3) without an actual cap or a
fourth moment. Complete row lengths and (3.2) also uniformly bound
E F_i^2 and E F_(i,1)^2, as required by the local-field estimates.

## 5. Reused local Gaussianization and scalar estimates

Sections 6, 7, and 9 of the frozen 612-line near-flat source identified
in Section 8 give the following consequences of (3.1)--(4.3):

    G_n=(1/n)E[||MX||_1-X^T M X]>=(8/25)mu_n-o(1),
    p_n=(1/n)sum_i P(sign(F_i)!=X_i)
                             <=(8/25)sqrt(mu_n)+o(1).      (5.1)

For clarity, that distinguished-coordinate Gaussianization lemma
assumes only bounded ||R||op, max_j |d_j|<=n^(-1/2), and sum_j d_j^2<=1
for F_i=sum_j d_j sign(G_j). We use the ACTUAL row d_j=M_ij, not a
row of M_b; its hypotheses hold exactly by completeness. No cap on M
occurs in that lemma. Alignment (4.3) and (4.2) give, before the
elementary scalar bounds,

    G_n=sqrt(kappa)(1/n)sum_i[sqrt(a^2+v_i)-a]+o(1),
    p_n=(1/pi)(1/n)sum_i arctan(sqrt(v_i)/a)+o(1).

The chord estimate on [0,V] and arctan x<=x prove (5.1), using only
the frozen enclosure 7/11<kappa<16/25. The same uniform rowwise
Gaussianization and scalar proofs apply unchanged; no new CLT,
field-covariance assertion, integral computation, or certificate is
needed. In particular no covariance bound for sign(MX) is asserted.

## 6. Adaptive bulk means and a full-source smoothing transfer

For r>=0 define the nonnegative, globally 2-Lipschitz function

    g_C(r)=r^2/(2C)-(r-2C)_+^2/(2C)
          =max_(0<=p<=1)(2pr-2Cp^2).                      (6.1)

For a realization of X set

    Y=clip((I+M_b/C)X,[-1,1]^n),
    r_(b,i)=(-X_i(M_b X)_i)_+.

Coordinatewise Y_i-X_i=-2X_i p_i where
p_i=min(r_(b,i)/(2C),1). The quadratic identity and ||M_b||op<=C
give deterministically

    Q_(M_b)(Y)-Q_(M_b)(X)>=sum_i g_C(r_(b,i)).              (6.2)

This identity permits a nonzero diagonal of M_b. It does not use a
Boolean-rounding identity for M_b. By (2.2) and C_X<=B_n I,

    E||M_o X||^2<=B_n tr M_o^2=o(n),
    |E Q_(M_o)(X)|<=B_n ||M_o||_1/2=o(n).

Since g_C and positive part are Lipschitz, if r_i=(-X_i F_i)_+ then

    E Q_(M_b)(Y)>=E Q_M(X)+sum_i E g_C(r_i)-o(n).          (6.3)

It remains to pay for M_o at the OUTPUT mean. We do not assume
Cov(Y) is bounded. For fixed epsilon>0 use the odd smoothing

    q_epsilon(t)=clip(t/epsilon,[-1,1]),
    X_epsilon=q_epsilon(G),
    Y_epsilon=clip((I+M_b/C)X_epsilon,[-1,1]^n).

For any fixed B with B_n<=B eventually, the map from a standard
Gaussian Z with G=R^(1/2)Z to Y_epsilon is odd and has Euclidean
Lipschitz constant at most 2sqrt(B)/epsilon. Gaussian Poincare
applied to every linear functional consequently gives

    E Y_epsilon=0,       Cov(Y_epsilon)<=4B epsilon^(-2)I.

One may use the finite-dimensional Gaussian Poincare inequality
Var f(Z)<=E||grad f(Z)||^2, obtained by comparing the positive-degree
Hermite coefficients on its two sides. It extends to Lipschitz maps
by smooth approximation and needs no nonsingularity of R. Thus

    |E Q_(M_o)(Y_epsilon)|
                        <=2B epsilon^(-2)||M_o||_1=o_epsilon(n).

Also each G_i is standard Gaussian and
P(|G_i|<=epsilon)<=sqrt(kappa)epsilon. Therefore

    E||X_epsilon-X||^2<=sqrt(kappa)n epsilon,
    E||Y_epsilon-Y||^2<=4sqrt(kappa)n epsilon.

Both output means have norm at most sqrt(n), so the BOUNDED bulk
energy obeys

    |E Q_(M_b)(Y_epsilon)-E Q_(M_b)(Y)|
                         <=2C kappa^(1/4)n sqrt(epsilon).

Since the FULL original M has zero diagonal, independent Boolean
rounding of any mean y in [-1,1]^n has expected Q_M equal to Q_M(y).
Thus alpha(A)>=E Q_M(Y_epsilon)/n. Combining this with (6.3) gives

    alpha(A)>=E Q_M(X)/n+(1/n)sum_i E g_C(r_i)
                     -2C kappa^(1/4)sqrt(epsilon)-o_epsilon(1).
                                                               (6.4)

The limit order is n->infinity at each FIXED epsilon, followed by
epsilon->0. There is no uncontrolled pairing of M_o with Y or with
sign(MX), and no change of the original source or objective.

## 7. The existing strict gain survives unchanged

For any fixed 0<=p<=1, (6.1) gives the sharper zero-field form

    g_C(r)>=2pr-2Cp^2 1_(r>0).

Here 2r_i=|F_i|-X_i F_i, and {r_i>0} is contained in the
sign-disagreement event, including the convention sign(0)=+1.
Consequently, with p=1/10 and C=5/3, (5.1) yields

    (1/n)sum_i E g_C(r_i)>=pG_n-2Cp^2 p_n
       >=(4/125)mu_n-(4/375)sqrt(mu_n)-o(1).

By (4.2), liminf sqrt(mu_n)>=sqrt(1-kappa)>3/5. The function
3s^2-s is increasing for s>=3/5, so the last liminf is at least

    (4/375)[3(3/5)^2-3/5]=16/3125.

Apply (3.4) and (6.4), with the stated limit order, to prove the
first inequality in (1.2). Finally 5kappa/8>35/88,
16/3125>1/200, and 2/5-35/88=1/440 give its strict comparison,
since 1/200-1/440=3/1100. This proves the theorem.

## 8. Frozen dependencies and contribution boundaries

The unchanged local Gaussianization, empirical scalar reductions,
and coarse kappa enclosure are imported from Sections 6, 7, and 9 of

    original_mo_original_source_near_flat_strict_gain.md, 612 lines
    SHA256 7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298.

The adaptive mean mechanism in Section 6 restates and specializes
the deterministic clipping mechanism from Section 7 of

    original_mo_all_law_adaptive_nuclear_gain.md, 553 lines
    SHA256 0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51.

Those full sources were directly read and remain frozen. Their
underlying lemma and enclosure provenance remains as disclosed there.
This extension does not claim a new proof of their unchanged CLT.

Root supplied the actual-second-moment tail argument and the strategy
of smoothing a bounded-bulk update to control its full-source outlier
energy by Gaussian Poincare. The exact worker supplied the full-projector
frame/baseline repair, checked the trace-norm and smoothing links,
and authored this complete extension. The docs-gate worker contributed
the trace-norm higher-chaos/alignment repair and the direct smoothed
adaptive-mean gain transfer, with contributing independent checks of
the frame links. These are contributions, not independent whole-source
reviews. A full frozen-source reviewer must disclose any overlap.

Only this /tmp source and its author receipt were authored in this
step. No canonical file, frozen predecessor, covariance, cross block,
active cell, or paired-source definition was changed. Publication,
documentation gates, commits, and backups remain root's workflow.
The full near-flat law (1.1) is still a genuine conditional premise;
all other laws and the global original MO target remain OPEN.
