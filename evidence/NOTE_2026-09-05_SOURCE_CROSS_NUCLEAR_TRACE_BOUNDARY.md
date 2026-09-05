# Actual source/cross nuclear coupling and a precise trace-relaxation boundary

2026-09-05. Analytic bounded note. No mathematical computation, checker,
optimization, signing search, or numerical metric evaluation was run.

Part I proves a finite necessary inequality for ACTUAL complete paired
signings. Part II exhibits a FORMAL spectral/block relaxation satisfying
the stated necessary inequalities but failing the desired drift-plus-
ellipsoid certificate for every shifted Gaussian sign threshold and
every signed metric, including their endpoint limits. The formal parameters are NOT
claimed to arise from complete signings or an actual active optimizer.
This is not a counterexample to the desired original inequality.

## 1. Actual matrices and a cap-free finite coupling inequality

Let n>=2, N=2n, and let the actual complete signing be

    K=[[A,B],[B^T,-A]],       D=diag(D_L,D_R)>0,
    D-K>=0,                  D+K>=0.

Here A is a complete symmetric zero-diagonal signing and B is a complete
cross sign matrix. Put

    S=tr D,      dbar=S/N,      delta=S tr(D^(-1))/N^2-1,
    T=D^(-1/2)KD^(-1/2)=[[A_L,W],[W^T,-A_R]],
    A_L=D_L^(-1/2)AD_L^(-1/2),
    A_R=D_R^(-1/2)AD_R^(-1/2),
    W=D_L^(-1/2)BD_R^(-1/2),
    alpha=Phi(A)/n^(3/2),      kappa=2/pi.

The actual T is a contraction. Let nu be the empirical probability
measure of the n squared singular values of W, including zeros, and
let m=integral y dnu(y). Thus nu is supported in [0,1]. Then

    integral sqrt(1-y) dnu(y)
      >= [kappa sqrt(m)/(2alpha)] (1-1/n)/(1+delta)
          -sqrt[(2N/dbar^2)(2delta+delta^2)].                       (1.1)

No optimality of D, small canonical gap, or fixed trace cap is needed
for (1.1). In particular, as delta tends to zero,

    integral sqrt(1-y) dnu(y)
                   >=kappa sqrt(m)/(2alpha)-O(sqrt(delta)+1/n),    (1.2)

with an absolute error constant for 0<=delta<=1. The measure and cross
matrix in both statements are the actual original ones.

The original nuclear prerequisite is the complete 262-line source
`original_mo_original_phase_spectral_moment.md`, SHA256
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
Its finite bound is

                    Phi(A)>=kappa n^2(n-1)/(2tr|A|).              (1.3)

## 2. Direct nuclear transfer: no auxiliary trimming is needed

Write t_i=d_i/dbar and r_i=1/t_i. Exactly,

    sum_i t_i=N,       sum_i r_i=N(1+delta),
    sum_i sqrt(r_i)>=N.

The last inequality is Jensen for the convex function t^(-1/2).
Literal off-diagonal sign squares in K, followed only by adding
nonnegative diagonal terms, give

    ||K/dbar-T||_F^2
       =dbar^(-2) sum_(i!=j)(1-sqrt(r_i r_j))^2
       <=dbar^(-2)[N^2+(sum r_i)^2-2(sum sqrt(r_i))^2]
       <=(N^2/dbar^2)(2delta+delta^2).                            (2.1)

The trace-norm triangle inequality and ||E||_1<=sqrt(n)||E||_F
therefore imply, separately for either internal block,

    tr|A_L|/n, tr|A_R|/n
      >=tr|A|/(n dbar)-sqrt[(2N/dbar^2)(2delta+delta^2)].           (2.2)

Set m_0=n/dbar^2. From (1.3), the first term in (2.2) is at least

                         kappa sqrt(m_0)(1-1/n)/(2alpha).         (2.3)

The exact original cross sign squares give

    m=[sum_L d_i^(-1)][sum_R d_i^(-1)]/n=m_0 ell h,
    ell=(1/n)sum_L r_i,       h=(1/n)sum_R r_i,
    ell+h=2(1+delta).

Arithmetic-geometric mean yields ell h<=(1+delta)^2 and hence

                         sqrt(m_0)>=sqrt(m)/(1+delta).            (2.4)

Finally, T^2<=I implies the two ACTUAL block inequalities

    A_L^2+WW^T<=I,             A_R^2+W^TW<=I.

The square root is operator monotone on positive semidefinite matrices,
so no commutation is required to conclude

    |A_L|<=(I-WW^T)^(1/2),     |A_R|<=(I-W^TW)^(1/2).             (2.5)

For completeness, square-root monotonicity follows from
sqrt(X)=(1/pi)integral_0^infinity t^(-1/2)X(tI+X)^(-1)dt:
inverse Loewner reversal makes each integrand increasing in X.
Positive-definite regularization handles singular endpoints.
Taking traces in (2.5), then using (2.2)--(2.4), proves (1.1).

The contraction row-square inequalities give
S>=(N-1)tr(D^(-1)), so Cauchy--Schwarz gives dbar^2>=N-1.
Also (1.3) and tr|A|<=n sqrt(n-1) imply
alpha>=kappa sqrt(1-1/n)/2. These facts, together with m<=1,
bound all coefficients in (1.1) and prove the cap-free error in (1.2).

This proof is a direct comparison of nuclear moments. It does not
replace K by K/dbar in the Gaussian construction, does not assert that
K/dbar is contractive, and does not alter any final original/weighted
cell or its representative.

## 3. The particular retained relaxation to be tested

The additional actual near-scalar normalization previously proved is
in the complete 280-line source
`original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
In its delta-to-zero fixed-cap limit, with mu the full weighted
spectral measure and mu_j=integral |lambda|^j dmu, it requires

    mu supported in [-1,1],       mu_2=r,
    u>=max{kappa(1+r)mu_3/(2r), kappa r/mu_1}.                    (3.1)

Here u=c_D/n may be substituted for the original norm ratio only on
the separately assumed active original face p=q_A=0, c=Phi(K).
With f=Phi(K)/n^(3/2), actual scalar/near-scalar compatibility also
gives, at leading order,

                          m=r/2,       u=f sqrt(m).               (3.2)

The relaxation considered below retains (3.1)--(3.2), the algebraic
block contraction and its internal second moments, the source nuclear
and cubic necessary bounds explicitly checked in Section 5, the cross
cubic necessary bound, and (1.2) at zero error. It does NOT retain
entrywise complete-signing realizability or the existence of actual
Boolean states witnessing p=q_A=0 and c=Phi(K). Those are missing
conditions, not conclusions inferred from the listed trace data.

## 4. Explicit formal spectral/block data

Choose the formal numerical parameters

    alpha=2/5,       f=4/3,       u=kappa,
    m=9kappa^2/16,   r=2m,       s=kappa u=kappa^2,
    a=sqrt[m/(1-m)].

Then 0<m<1/4, a<1, and u=f sqrt(m), as required in (3.2).
Define the formal cross and full spectral measures by

    nu=(1-m)delta_0+m delta_1,
    mu=(m/2)(delta_1+delta_(-1))
                      +((1-m)/2)(delta_a+delta_(-a)).             (4.1)

The internal absolute-value law is

                         chi=m delta_0+(1-m)delta_a.               (4.2)

These laws are mutually consistent with an algebraic commuting block
contraction: let P be an orthogonal projection with normalized trace m,
let W=P, and let H be zero on ran(P) and have eigenvalues +a and -a
on equal halves of its complement. Set

                          T_form=[[H,W],[W,-H]].

Then HW=WH=0, W^2=P, and H^2=a^2(I-P). Thus
T_form^2=diag(H^2+P,H^2+P)<=I, its full law is (4.1), and its
internal squared moments and cross squared moment all equal m.
For irrational m this is interpreted as a weighted trace model or as
a limit of rational-rank models, not as a finite empirical equality.

Even zero internal diagonals and constant diagonal second moments can
be imposed in rational approximations: take a real normalized Sylvester
Hadamard basis, a projection onto an even number of its columns, and
equal numbers of +a and -a eigenvalues on the complement. Conjugation
in that basis gives diag H=0, diag H^2=m, and diag W^2=m.
This observation realizes ONLY these real-matrix algebraic data.
The entries are not asserted to have the magnitudes of complete
signings. In particular alpha and f have NOT been identified with
the Boolean norms of these real matrices, and u has NOT been realized
as an active original cross value. No actual-signing example is claimed.

## 5. Exactly which necessary inequalities the formal data pass

Directly from (4.1),

    mu_1=m+sqrt[m(1-m)],
    mu_2=2m=r,
    mu_3=m(1+a).

The first full normalization term is

    kappa(1+r)mu_3/(2r)=kappa(1+2m)(1+a)/4<=kappa=u,

because m<=1/2 and a<=1. The nuclear term is likewise at most u:

    kappa r/mu_1
       =2kappa sqrt(m)/(sqrt(m)+sqrt(1-m))<=kappa.                 (5.1)

The scalar canonical-gap parameter of these formal moments is

                          gamma=1-mu_3/r=(1-a)/2>0.

This is consistent with a positive-gap near-scalar branch. It is not
claimed to be the gap of an actual complete signing.

The source/cross necessary bound is satisfied with substantial slack:

    integral sqrt(1-y)dnu=1-m,
    kappa sqrt(m)/(2alpha)=15kappa^2/16,
    1-m-15kappa^2/16=1-(3/2)kappa^2>1/3.                         (5.2)

For the internal law (4.2), the normalized source nuclear bound is

                         alpha>=kappa/[2sqrt(1-m)].               (5.3)

It holds because m<1/4 and kappa<2/3 imply the right side is less
than 2/(3sqrt(3))<2/5. The internal nonzero absolute eigenvalue is a,
with balanced signs. Its common-zero-odd-diagonal cubic phase bound
requires the SAME inequality (5.3); the weaker row-variance phase
bounds are therefore satisfied as well. This checks the relevant
internal nuclear/cubic necessary inequalities, not an actual upper
bound Phi(A)<=alpha n^(3/2) for a realized signing.

The cross cubic necessary bound also passes: its normalized form is
u>=kappa[integral y^(3/2)dnu]/m, whose right side equals kappa=u.
Thus it is not omitted in declaring the listed relaxation consistent.
Its prerequisite is `original_mo_cross_singular_moment_rounding.md`,
SHA256
`6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756`.
For a genuine active pure-cross state, c=Phi(K) and beta(B)<=Phi(K)
would also imply c=beta(B), as required by that cross normalization.

No claim is made that these moment checks supply actual optimal
frames, their tensor-feature correlations, or any realization theorem
for finite templates. Additional frame or entrywise constraints are
not silently folded into the relaxation. The stopping statement is
limited to the explicitly retained constraints above.

## 6. The actual-measure ellipsoid functional being tested

The full derivation of the numerical reference functional is in
`original_mo_small_gap_pure_cross_upper.md`, final 312-line source,
SHA256
`035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
For the centered sign law w=1, k=kappa, a positive metric 0<=t<1
and a supplied cross measure, the normalized reference upper is

    U_nu(t)=t sqrt[(1-u) integral A_s,t(y)dnu(y)]
                +(1-t)sqrt(kappa) sqrt[integral B_s,t(y)dnu(y)],

    A_s,t(y)=[1+(t^2-2t+s(1-2t))y+s t^2y^2]/(1-t^2y)^2,
    B_s,t(y)=[1+(t^2-2st)y]/(1-t^2y)^2.                          (6.1)

The original target is U_nu<=sqrt(2)alpha, NOT U_nu<=f/2.
Here alpha=2/5, so the target equals 2sqrt(2)/5.

For the formal endpoint law (4.1), formula (6.1) is exactly

    U_nu(t)=t sqrt(1-kappa)
                  sqrt[1-m+m(1+s)/(1+t)^2]
             +(1-t)sqrt(kappa)
                  sqrt[1-m+m(1+t^2-2st)/(1-t^2)^2].              (6.2)

We do not replace the two expectations by separate Jensen/chord
bounds. Equation (6.2) evaluates the SAME formal law in both terms.
After absorbing 1-t into the second square root, (6.2) has a finite
continuous limit at t=1. The proof below covers that limit too.

## 7. Uniform analytic lower bound for every positive metric

Let A_*=sqrt[(1-kappa)(1+s)]. View the two terms in (6.2) as
Euclidean norms of two nonnegative two-coordinate vectors. The triangle
inequality gives U_nu(t)>=sqrt[a(t)^2+b(t)^2], where

    a(t)=sqrt(1-m)[t sqrt(1-kappa)+(1-t)sqrt(kappa)],
    b(t)=sqrt(m)[t A_*+sqrt(kappa)sqrt(1+t^2-2st)]/(1+t).

Since sqrt(kappa)>=sqrt(1-kappa), and since
sqrt(1+t^2-2st)>=1-st on [0,1],

    a(t)^2>1/4,
    b(t)>=sqrt(m) min{sqrt(kappa),
                                  [A_*+(1-s)sqrt(kappa)]/2}.     (7.1)

The first strict bound uses m<1/4 and kappa<2/3. The minimum in
the second bound comes from the two endpoints of the affine ratio
[sqrt(kappa)+t(A_*-s sqrt(kappa))]/(1+t); this is an exact scalar
identity, not a numerical optimization.

Only the following coarse interval is needed:

                              63/100<kappa<16/25.                 (7.2)

It follows directly from the ALREADY VERIFIED pi enclosure
31415926/10000000<pi<31415927/10000000, whose baseline checker
SHA256 is
`d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`
and baseline result SHA256 is
`fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
The two integer comparisons for (7.2) are
1979203401<2000000000 and 500000000<502654816.
That baseline was reused, not rerun.

All remaining comparisons can be written explicitly as rationals:

    A_*^2>(9/25)(1+(63/100)^2)=125721/250000>49/100,
    sqrt(kappa)>79/100,       1-s>369/625,
    [A_*+(1-s)sqrt(kappa)]/2>72901/125000>29/50,
    m>35721/160000>11/50.

Thus the minimum in (7.1) exceeds 29/50. Uniformly for 0<=t<=1,

    U_nu(t)^2>1/4+(11/50)(29/50)^2
                  =40501/125000>40000/125000=8/25.               (7.3)

Consequently every positive metric and its endpoint limit satisfy

                              U_nu(t)>2sqrt(2)/5.                 (7.4)

The positive squared margin in (7.3) is exactly 501/125000. This
is an analytic rational comparison; there is no new checker result.

## 8. Negative metrics cannot evade this particular example

For a negative metric -t, t>=0, even spectral pairing changes s to
-s and replaces the first prefactor 1-u by 1+u in (6.1).
For the present endpoint law, the squared first term cannot decrease:
at y=0 the relevant coefficient changes from 1-kappa to 1+kappa;
at y=1 its numerator difference is

    (1+kappa)(1-kappa^2)-(1-kappa)(1+kappa^2)
                                   =2kappa(1-kappa)>0.

The denominator there is the same positive (1+t)^2. Also

    B_(-s),t(y)-B_s,t(y)=4st y/(1-t^2y)^2>=0.

Thus the negative-metric functional is at least U_nu(t). Formula
(7.4) therefore holds for every signed interior metric and both
continuous endpoint limits. No metric scan or optimization is used.

## 9. Every shifted Gaussian sign threshold has the same obstruction

The preceding centered-law calculation extends analytically to every
shifted threshold h in the SAME formal trace/block model. The threshold
coefficients below are those of Section 4.1 of
`original_mo_diagonal_majorizer_cross_covariance.md`, SHA256
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`;
the original drift is from Section 4.6 of
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
Write

    z=|s_h|=|2 Phi_Gauss(h)-1|,   w=1-z^2,
    k=kappa exp(-h^2),           a_h=k/w,       s=a_h kappa.

For finite h one has w>0. If X,Y are independent standard Gaussians,
the square [-|h|,|h|]^2 is contained in the disk of radius sqrt(2)|h|.
Consequently

    z^2=Pr(|X|<=|h|, |Y|<=|h|)
       <=Pr(X^2+Y^2<=2h^2)=1-exp(-h^2),

so w>=exp(-h^2), a_h<=kappa, and 0<=s<=kappa^2.

The noise covariance is now w times the reference covariance having
spectral coefficient s. Hence its normalized ellipsoid functional is
sqrt(w) U_(nu,s)(eta), where (6.1)--(6.2) define U_(nu,s) with u=kappa
still fixed. The second-term sqrt(kappa) is the Boolean completion-square
constant, not the threshold coefficient; it is unchanged.

For positive metrics, (7.1) remains valid with this value of s. Set

    C(s)=[sqrt((1-kappa)(1+s))+(1-s)sqrt(kappa)]/2.

Then

    C'(s)=sqrt(1-kappa)/(4sqrt(1+s))-sqrt(kappa)/2<0,

using kappa>1/2. Thus C(s)>=C(kappa^2)>29/50 by the already proved
centered-law rational bound. Together with sqrt(kappa)>79/100,
a(t)^2>1/4 and m>11/50, this proves

    U_(nu,s)(t)>L_0:=sqrt(40501/125000)>2sqrt(2)/5               (9.1)

uniformly for 0<=t<=1 and 0<=s<=kappa^2. We use monotonicity of
C as a whole: the separate bound A_*>7/10 from Section 7 is NOT
asserted for all s. For negative metrics the first unit-atom numerator
difference becomes

    (1+kappa)(1-s)-(1-kappa)(1+s)=2(kappa-s)>0,

and the second-term difference remains 4st y/(1-t^2y)^2>=0.
Thus (9.1) also holds for all negative metrics and their endpoint limits.

The original drift is (p-q_A)/2+s_h c. On the formal active pure-cross
face p=q_A=0, c=f n^(3/2), its absolute value divided by 2n^(3/2)
is exactly z f/2. This is the drift of the original unweighted energies,
not a substituted weighted drift. Therefore the entire normalized
formal drift-plus-ellipsoid certificate is at least

    z f/2+sqrt(1-z^2)L_0
       >=z f/2+(1-z)L_0
       >=min{f/2,L_0}>sqrt(2)alpha,                              (9.2)

since 0<=z<=1, f/2=2/3, and alpha=2/5. At infinite-threshold limits
z=1 the noise vanishes and the drift alone is f/2, so no endpoint is
excluded. The statement concerns this reference certificate on the
explicit formal data, not a realized signing or an actual optimizer.

## 10. Exact stopping point and non-claims

The finite inequality (1.1) is a genuine new necessary coupling for
actual source and cross matrices. Its proof keeps all original
matrices, norms, and cells. It is not retracted by the formal example.

The example shows only that this coupling, the retained full/source
nuclear and cubic inequalities, and the specified block trace data do
not certify the desired smaller-alpha ellipsoid upper. The original
target is kept throughout; it has not been replaced by f/2 or by a
universal sqrt(2) normalization. Section 9 extends the centered-law
calculation to every shifted Gaussian sign threshold and signed metric
within this SAME explicitly retained formal relaxation.

The missing implication must impose additional information from actual
complete-signing entries, Boolean active states, conditional optimality,
or their frame geometry. In particular the formal example does not
supply zero-internal-energy states whose cross value saturates Phi(K).
It is NOT a counterexample to the original inequality, an impossibility
theorem for all spectral methods, or a closure result. This specific
trace-only route stops at the explicitly identified relaxation boundary.

The exact worker authored the direct Frobenius/nuclear transfer and
the formal relaxation calculation. The docs-gate reviewer independently
checked the transfer and supplied the directional AM-GM refinement
(2.4) before this note was written. That contribution must be disclosed
in the subsequent review record. The root authored the all-threshold
extension in Section 9; the exact worker independently checked it and
integrated it here. No mathematical computation was run.
