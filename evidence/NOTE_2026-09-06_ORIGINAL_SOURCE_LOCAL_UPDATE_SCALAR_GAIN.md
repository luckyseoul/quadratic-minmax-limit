# Scalar local-update gain with an explicit operator cap

2026-09-06. Analytic scalar implication. No mathematical computation,
scan, checker, numerical integral, optimization, or parameter search.

This note proves ONLY the scalar implication requested by root. It
does not prove the proposed Gaussian local-field representation, its
variance constraints for an actual source, or the claimed original-
energy improvement from independent Boolean updates. Those premises
are being checked separately by root and the proof worker.

This revision retains the operator cap C in the update penalty. It
must not be silently replaced by the nonzero limiting spectral value
1/sqrt(rho)=5/4. The live bounded-operator premise is C<=5/3+o(1).
The leading result uses an admissible FIXED update probability, not
an unconstrained quadratic optimizer whose probability might exceed one.

## 1. Scalar setting and the fixed-update conclusions

Let kappa=2/pi, rho=16/25, and define

    a=sqrt(kappa/rho),
    V=2(1-kappa)/rho.

Let v have any probability distribution supported on [0,V], with
mean mu=E v>0. Empirical averages are a special case. Set

    G=sqrt(kappa) E[sqrt(a^2+v)-a],
    p=(1/pi) E arctan(sqrt(v)/a).

For an explicit C>0 define the proposed scalar update lower

    J_C(epsilon)=epsilon G-2C epsilon^2 p,   0<=epsilon<=1.

C here is only the operator-cap coefficient in this penalty, not a
diagonal trace-cap parameter. No homogeneity of v is assumed.

Under the STRONGER proposed mean premise `mu>=1-kappa`:

1. If C<=5/3, the fixed probability epsilon=1/10 satisfies

       J_C(1/10)>16/3125>1/200,
       5kappa/8+J_C(1/10)>2/5+3/1100.                     (1.1)

2. Even if only C<=3 is available, epsilon=1/20 satisfies

       J_C(1/20)>9/3125>1/400,
       5kappa/8+J_C(1/20)>2/5+1/4400.                     (1.2)

The rational thresholds 1/200 and 1/400, and their resulting objective
margins, remain eventually valid if the respective cap is enlarged
by o(1) and the mean lower bound is weakened by o(1), with the stated
a, V and scalar formulas retained. Additional o(1) errors in a
separately justified objective lower are absorbed by the strict margins.

Section 5 records the unconstrained quadratic value as an OPTIONAL
scalar bound. Section 6 records the earlier half-mean calculation at
C=5/4 only as an exact-scale specialization, not the live robust premise.

## 2. Reused kappa interval

Only

                       7/11<kappa<16/25                       (2.1)

is needed. This follows from the already verified pi enclosure
`31415926/10000000<pi<31415927/10000000` recorded in
`original_mo_source_cross_nuclear_trace_boundary.md`, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.
Its upper endpoint is below 22/7 because

                    7*31415927=219911489<220000000.

Hence kappa>7/11. The same source already gives kappa<16/25.
No pi certificate is rerun. The comparisons below are direct analytic
rational inequalities, not newly executed certificate outputs.

## 3. A distribution-free chord for G

The function `sqrt(a^2+v)-a` is concave on [0,V] and vanishes
at zero. It is therefore bounded below by its endpoint chord:

    sqrt(kappa)[sqrt(a^2+v)-a]>=c v,
    c=sqrt(kappa)/(sqrt(a^2+V)+a)
      =(4/5)sqrt(kappa)/(sqrt(kappa)+sqrt(2-kappa)).

Consequently G>=c mu. Moreover

                             c>8/25,                           (3.1)

because this is equivalent to
`3sqrt(kappa)>2sqrt(2-kappa)`, or `13kappa>8`.
The latter follows already from kappa>7/11>8/13.
All quantities squared in this comparison are positive.

## 4. Fixed probabilities give the required robust gain

For every nonnegative x, arctan(x)<=x. Cauchy--Schwarz gives
E sqrt(v)<=sqrt(mu). Hence

                         p<=sqrt(mu)/(pi a).

The exact relation pi=2/kappa and rho=16/25 imply

    pi a=5/(2sqrt(kappa))>25/8.

Thus, using (3.1),

                 G>(8/25)mu,       p<(8/25)sqrt(mu).            (4.1)

Suppose C<=5/3 and take epsilon=1/10, which is in [0,1]. Then

    J_C(1/10)>(4/125)mu-(4/375)sqrt(mu)
              =(4/375)(3r^2-r),       r=sqrt(mu).              (4.2)

The strong mean premise and kappa<16/25 give r>3/5. The polynomial
3r^2-r is increasing on r>=3/5, since 6r-1>0 there. Therefore

    J_C(1/10)>(4/375)[3(3/5)^2-3/5]
              =16/3125>1/200,
    16*200=3200>3125.                                       (4.3)

The baseline is strictly above 35/88 by (2.1), and

    2/5-35/88=1/440,
    1/200-1/440=3/1100.

This proves (1.1). The choice epsilon=1/10 is fixed; no condition
on the location of a quadratic maximizer is required.

For the more permissive cap C<=3 choose epsilon=1/20. In the same way,

    J_C(1/20)>(2/125)mu-(3/625)sqrt(mu)
              =(10r^2-3r)/625.

This polynomial is increasing for r>=3/5. At r=3/5 its value is
9/3125>1/400, since 9*400=3600>3125. Together with
`1/400-1/440=1/4400`, this proves (1.2).

All probability bounds are distribution-free. In particular neither
the variances nor their distribution are assumed homogeneous.

## 5. Optional unconstrained quadratic value

Because mu>0, one has p>0. Maximizing the scalar quadratic over
all nonnegative real epsilon, without the constraint epsilon<=1, gives

    epsilon_*=G/(4Cp),       I_C=G^2/(8Cp).

The same chord and probability bounds imply

    I_C>=pi a c^2 mu^(3/2)/(8C)>mu^(3/2)/(25C).              (5.1)

For the strong mean premise and C<=5/3, this exceeds
81/15625>1/200, with the latter comparison given by 16200>15625.
For C<=3 it exceeds 9/3125>1/400.

These are scalar bounds on the UNCONSTRAINED optimum. They are not
asserted to be implementable update gains unless epsilon_*<=1 is
proved separately. The fixed probabilities of Section 4 avoid that
issue and are the relevant conclusions for the proposed actual updates.

## 6. Earlier half-mean specialization: C=5/4 only

The original scalar request used C=1/sqrt(rho)=5/4 and the weaker
mean mu>=(1-kappa)/2. At this EXACT-SCALE specialization, (5.1) gives

    I_(5/4)>(4/125)mu^(3/2)
            >(4/125)(9/50)^(3/2)=27sqrt(2)/15625>1/420.

The last comparison follows from sqrt(2)>7/5 and
`189/78125>1/420`, since 189*420=79380>78125.
Thus `5kappa/8+I_(5/4)>2/5+1/9240`, because
`1/420-1/440=1/9240`.

This remains a correct optional scalar calculation. It does not
replace the actual cap 5/3 by the limiting nonzero spectral value
5/4, and it does not establish admissibility of epsilon_*.
No stronger-operator version of this half-mean conclusion is claimed.

## 7. Asymptotic stability and exact contribution scope

For a sequence with mu>=1-kappa-o(1), its eventual square-root
mean is at least 3/5 up to a vanishing error; in fact the strict
inequality 1-kappa>9/25 leaves fixed room at that endpoint.
If C<=5/3+o(1), the fixed-epsilon lower (4.2) has only an
additional o(1) penalty: p is uniformly bounded by 1/2 and epsilon
is fixed. Its liminf is at least 16/3125>1/200. Thus the gain
is eventually above 1/200. The analogous argument for C<=3+o(1)
gives liminf at least 9/3125>1/400. This proves the stated stability
without convergence or homogeneity of the individual variances.

The robust fixed-update implication uses the STRONG mean premise.
It does not itself assert a theorem about the original symmetric
Boolean quadratic norm. In particular the local Gaussian
approximation, actual variance inequalities, covariance/dependence
effects of updates, and the objective-improvement formula must still
be established in the separate actual-source argument.

Root supplied the scalar setting, corrected the relevant operator cap,
and relayed the proof worker's stronger proposed mean constraint.
The docs-gate worker derived the distribution-free chord/probability
bounds and the earlier exact-scale margins. Root proposed the live
C=5/3 extension and fixed-epsilon=1/10 corollary; the docs-gate
worker checked its coefficients and stability and supplied the
C<=3, epsilon=1/20 fallback. This is a contribution-disclosed
authored scalar lemma, not an independent full-source review of the
underlying actual-source update construction.
No canonical repository file was edited or mathematical job executed.
