# A strengthened formal trace profile still exceeds the original certificate target

2026-09-05. Bounded analytic diagnostic. No mathematical computation,
numerical metric evaluation, scan, checker, or parameter search was run.
This is a formal analytic diagnostic, not an actual counterexample.

The profile below is FORMAL. It passes the explicitly listed source,
full-spectrum, and cross-entry necessary inequalities, including the
new complete-cross gain. Nevertheless its SAME original drift-plus-
ellipsoid reference certificate exceeds the desired target for every
shifted Gaussian sign threshold and every signed metric. No complete
signing, actual Boolean norm, or active optimizer is claimed to realize it.

## 1. Frozen prerequisites and the fixed profile

The earlier actual-coupling/formal-boundary source is
`original_mo_source_cross_nuclear_trace_boundary.md`, 444 lines, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.
It supplies the general reference functional, original drift conventions,
and already verified coarse enclosure

                         63/100<kappa=2/pi<16/25.                (1.1)

The general measure formula used here is (2.2) of
`original_mo_small_gap_pure_cross_upper.md`, final SHA256
`035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
It is tested here as a numerical reference on formal data, not invoked
as a theorem asserting that the formal profile is an actual source.

The near-scalar full-spectrum source is
`original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
The new actual complete-cross gain, transferred to full weighted moments,
is the complete 364-line source
`original_mo_near_scalar_cross_spectral_gain.md`, SHA256
`ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9`.

Fix the proposed formal parameters, without searching over parameters:

    alpha=2/5,       f=4/3,       u=4/5,
    m=9/25,         r=18/25,     a=3/4,
    nu=(16/25)delta_0+(9/25)delta_1,
    mu=(9/50)(delta_1+delta_(-1))
                        +(8/25)(delta_(3/4)+delta_(-3/4)).       (1.2)

Thus u=f sqrt(m), r=2m, and a=sqrt(m/(1-m)). The same algebraic
block construction as in the earlier diagnostic gives these moments:
W=P, H=0 on ran(P), H has balanced eigenvalues +a and -a on its
complement, and

    T_form=[[H,P],[P,-H]],       HP=PH=0,
    H^2=a^2(I-P),               T_form^2<=I.

The normalized trace of P is m. The internal absolute-value law is
(9/25)delta_0+(16/25)delta_(3/4). These are consistent trace/block
data, not assertions about complete-signing entries. As before, a
weighted trace model or limiting rational-rank models suffice for this
FORMAL test; no finite-entry or Boolean realization is inferred.

## 2. The strengthened retained inequalities all pass

The full absolute moments are exactly

    mu_1=21/25,       mu_2=18/25=r,       mu_3=63/100.

Consequently the two retained full-normalization lower bounds on u are

    kappa(1+r)mu_3/(2r)=301kappa/400<4/5,
    kappa r/mu_1=6kappa/7<4/5.                                (2.1)

Both strict comparisons follow already from kappa<16/25. The scalar
canonical moment-gap label is 1-mu_3/r=1/8; this is not the asserted
gap of a realized complete signing.

The internal source nuclear and common-zero-odd-diagonal cubic
conditions coincide for this internal law. Their required lower bound is

                alpha>=kappa/[2sqrt(1-m)]=5kappa/8,

which holds strictly by kappa<16/25. The retained source/cross coupling
also passes:

    integral sqrt(1-y)dnu=16/25,
    kappa sqrt(m)/(2alpha)=3kappa/4<12/25.                      (2.2)

The old cross cubic condition is u>=kappa because
integral y^(3/2)dnu=m; it passes. The NEW cross-entry condition at
zero dispersion and this endpoint spectrum is

    u>=kappa+(sqrt(kappa)-kappa)m.

It too passes strictly, since (1.1) gives sqrt(kappa)<4/5 and hence

    kappa+(sqrt(kappa)-kappa)m
      =(16/25)kappa+(9/25)sqrt(kappa)
      <436/625<4/5=u.                                        (2.3)

This checks precisely the stated retained inequalities. It does not
assert that the original source has Boolean norm alpha n^(3/2),
that the paired norm is f n^(3/2), or that an actual state attains u.
Those omitted realizability and active-state conditions remain omitted.

## 3. Exact formal reference for every positive metric

Write the shifted-sign coefficients as

    z=|2 Phi_Gauss(h)-1|,       w=1-z^2,
    k=kappa exp(-h^2),         s=(k/w)u

for finite real h. The square [-|h|,|h|]^2 lies in the disk of
radius sqrt(2)|h|. For two independent standard Gaussians this gives
z^2<=1-exp(-h^2), hence w>=exp(-h^2). Therefore

                      0<=s<=kappa u<64/125.                    (3.1)

After factoring out sqrt(w), the exact positive-metric reference on
the SAME endpoint law nu is, for 0<=t<1,

    U_s(t)=t sqrt(1-u)
                    sqrt[1-m+m(1+s)/(1+t)^2]
       +(1-t)sqrt(kappa)
                    sqrt[1-m+m(1+t^2-2st)/(1-t^2)^2].          (3.2)

Both terms evaluate the same formal law. The square-root completion
constant sqrt(kappa) is independent of the threshold coefficient k.
Canceling 1-t inside the second square root gives its continuous
endpoint expression

    T_2=sqrt(kappa)
          sqrt[(1-m)(1-t)^2+m(1+t^2-2st)/(1+t)^2].              (3.3)

Thus (3.2), interpreted through (3.3), has a finite t=1 limit.

## 4. Uniform rational lower bound, without any metric scan

Denote the first term in (3.2) by T_1. Since s>=0 and 1+t<=2,

    T_1>=t sqrt[(1/5)(16/25+9/100)]
         =t sqrt(73/500)>=(19/50)t,                           (4.1)

where 73/500=365/2500>361/2500=(19/50)^2.

For the ratio in (3.3), t/(1+t)^2<=1/4 gives

    (1+t^2-2st)/(1+t)^2
      =1-2(1+s)t/(1+t)^2
      >=(1-s)/2>=61/250.                                    (4.2)

Using kappa>63/100 in (3.3) therefore yields, with x=1-t,

    T_2>=sqrt(A x^2+B),
    A=252/625,       B=34587/625000.                          (4.3)

The following single rational comparison supplies a global supporting
line for the square root:

    B[1-(19/50)^2/A]
      =355203/10000000
      >353440/10000000=(47/250)^2.                            (4.4)

Indeed the vector ((19/50)/sqrt(A),(47/250)/sqrt(B)) has norm
strictly less than one. Its scalar product with (sqrt(A)x,sqrt(B))
and Cauchy--Schwarz imply, for every x in [0,1],

                      sqrt(A x^2+B)>(19/50)x+47/250.

Combining this with (4.1)--(4.3) proves the uniform bound

    U_s(t)>19/50+47/250=71/125>2sqrt(2)/5                    (4.5)

for every 0<=t<=1 and every s in (3.1). The final target comparison
is the exact rational square inequality

    (71/125)^2=5041/15625>5000/15625=(2sqrt(2)/5)^2.            (4.6)

The squared margin is 41/15625. These are direct analytic rational
inequalities, not numerical evaluations or newly executed certificates.

## 5. Negative metrics, original drift, and threshold endpoints

For metric -t, the first prefactor changes from 1-u to 1+u and
the spectral coefficient changes from s to -s. Its zero-atom first
coefficient increases. Its unit-atom numerator difference is

    (1+u)(1-s)-(1-u)(1+s)=2(u-s)>0,

because s<=kappa u<u. The second term cannot decrease, since its
unit-atom numerator changes from 1+t^2-2st to 1+t^2+2st.
Thus every negative-metric reference is at least the positive one.
Continuous endpoint limits are covered by the cancellation in (3.3).

The ORIGINAL formal pure-cross drift is |s_h|c with p=q_A=0 and
c=f n^(3/2). After division by 2n^(3/2), it is exactly zf/2.
Consequently the entire formal drift-plus-ellipsoid certificate obeys

    zf/2+sqrt(1-z^2) U_s(eta)
      >=z(2/3)+(1-z)(71/125)
      >=71/125>sqrt(2)alpha,                                 (5.1)

where sqrt(1-z^2)>=1-z and f/2=2/3>71/125 were used.
This applies to both metric signs. At infinite thresholds the noise
vanishes and the original drift alone is 2/3, still above target.

The original target is sqrt(2)alpha=2sqrt(2)/5 throughout.
It has not been replaced by f/2; f/2 is only the actual drift scale.

## 6. Exact diagnostic conclusion and non-claims

The proposed changed profile passes the listed strengthened trace
constraints but its same reference certificate remains above the
original target for every signed metric and shifted threshold.
Thus the new actual-entry gain eliminates the earlier u=kappa
profile without, by itself, closing this strengthened trace relaxation.

This is not a lower bound on the actual Gaussian width and is not an
actual-signing or conditional-optimizer counterexample. It does not
invalidate the complete-cross gain or its near-scalar transfer.
Additional source-entry, Boolean-active-state, conditional-optimality,
or frame information could exclude this profile; none is supplied by
the listed moment checks alone. No original-problem closure or general
impossibility claim is made.

The root supplied this particular parameter profile and requested the
bounded check. The exact worker checked the retained inequalities and
derived the uniform supporting-line bound (4.1)--(4.6). The previously
verified kappa enclosure was reused, not rerun. During derivation, no
canonical repository file was changed and no publication was performed.
