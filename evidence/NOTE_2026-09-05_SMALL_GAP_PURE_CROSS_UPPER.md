# Evaluated actual-measure upper on the small-gap pure-cross face

2026-09-05. Analytic fixed-metric evaluation. The narrow final rational
enclosure is supplied by one fixed certificate, not by optimization or a
spectral-law ansatz. No mathematical computation was run by the author.

This note evaluates the ACTUAL singular-measure upper at the diagnostic
normalization f=sqrt(2). It does not treat that normalization as sufficient
for original convergence when Phi(A)/n^(3/2) may be below one half. It also
does not replace the maximum over all original/weighted cells by the
pure-cross face considered here.

## 1. Actual source, pure-cross face, and exact reference field

Use an actual complete paired signing

    K=[[A,B],[B^T,-A]],       N=2n,

an actual positive optimal same-diagonal majorizer D, and

    S=tr D,   dbar=S/N,   gamma=(S-tr|K|^3/(N-1))/S,
    delta=S tr(D^(-1))/N^2-1,
    eta_source=S/[N sqrt(N-1)].

Assume a fixed original norm cap, so eta_source is bounded, and consider
sequences with gamma tending to zero. The reviewed canonical-gap theorem
gives delta<=4 eta_source^2 gamma, hence delta tends to zero. Its source is
`original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.

Retain the actual matrices

    W=diag(D_L)^(-1/2) B diag(D_R)^(-1/2),
    L_D=[[0,W],[W^T,0]],                 ||L_D||op<=1.

Choose a representative WITHIN a final refined cell whose ORIGINAL
internal energies satisfy p=q_A=0. Its weighted cross value is c_D, and
put u=c_D/n. The positive pure-cross reference covariance is

    M_0=wn I_(2n)-k c_D L_D
       =wn(I-a u L_D),       a=k/w,                    (1.1)

when w>0 and 0<=k<=w. For w=0 the field vanishes. The original cell is
compared to this positive field by the reviewed gap theorem and the
weighted-field theorem. The all-shell metric-stability theorem permits
its scalar-I NUMERICAL trace upper at any fixed interior metric parameter,
with an additional o(n^(3/2)) error. Its final source is
`original_mo_diagonal_majorizer_metric_stability.md`, SHA256
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.

Its covariance-congruence proof applies unchanged to M_0: the only field
properties used there are positivity and the constant diagonal wn, both
of which M_0 has. The exact natural-D ellipsoid bound likewise holds for
this positive field on the same original cell.

The comparison does not assert that K/dbar is contractive, and does not
assert that every state in a weighted bin has cross value c_D. The exact
original cross constraint enters the natural-D upper first; its uniform
metric comparison justifies the numerical reference expression below.

For the standard centered sign law the parameters are

                  w=1,       k=kappa=2/pi,       a=kappa. (1.2)

The independent Gaussian padding, bin error, and cell-selection error of
the reviewed construction remain separate and are o(n^(3/2)). The present
note evaluates a specified face, not the supremum over other cells.

## 2. Exact measure functional and the actual first-moment relation

Let nu be the empirical probability measure of the n squared singular
values of the ACTUAL W, including zeros. Its support lies in [0,1]. Put

    m=integral y dnu(y),        s=a u.

For a fixed ellipsoid parameter 0<t<1 define

    A_s,t(y)=[1+(t^2-2t+s(1-2t))y+s t^2 y^2]
                                                 /(1-t^2 y)^2,
    B_s,t(y)=[1+(t^2-2s t)y]/(1-t^2 y)^2.             (2.1)

Evenness of the spectrum of L_D and the commutation in (1.1) give the
exact normalized numerical reference upper

    W_flat(t)/(2 n^(3/2))
      =sqrt(w) {t sqrt((1-u) integral A_s,t dnu)
                    +(1-t)sqrt(kappa) sqrt(integral B_s,t dnu)}.
                                                               (2.2)

For example, if F=(I-t L_D)^(-1), its combined trace is
`tr[M_0(F-(1-t)F^2)]=wn t tr[(I-sL_D)(I-L_D)F^2]`.
The two spectral signs at each singular value give exactly A_s,t;
the trace of M_0 F^2 similarly gives B_s,t. Both functions are
nonnegative because their unsymmetrized factors are nonnegative.

The mean m is not chosen independently. Since B_ij^2=1 exactly,

    m=[sum_(i<=n) d_i^(-1)]
          [sum_(j>n) d_j^(-1)]/n.                       (2.3)

Let r_i=dbar/d_i. The reviewed dispersion bounds imply

    sum_i (sqrt(r_i)-1)^2<=N delta,
    sum_i r_i=N(1+delta).

Cauchy--Schwarz gives

    sum_i |r_i-1|
       <=N sqrt(delta)(1+sqrt(1+delta)).                (2.4)

Consequently each of the two sums of inverse diagonals in (2.3) is
`(n/dbar)(1+o(1))`, uniformly in the chosen cell. Therefore

    m=n/dbar^2+o(1)=1/(2 eta_source^2)+o(1).             (2.5)

If its original c=f_n n^(3/2), with f_n tending to a positive f,
uniform original/weighted energy compatibility gives

    u=c/(n dbar)+o(1),         m=u^2/f_n^2+o(1).        (2.6)

The actual-law statement (2.6), not a Dirac measure assumption, is the
only spectral constraint used in the following fixed-metric evaluation.

## 3. One fixed metric: concavity for A and convexity for B

Set t=3/5, r=t^2=9/25, and write A_s=A_s,t, B_s=B_s,t. For every
0<=s<=2/3, A_s is concave and B_s is convex on [0,1].

Here is a direct rational proof. Write

    A_s(y)=1-y(alpha+beta y)/(1-r y)^2,
    alpha=3/25+s/5,       beta=(9/25)(9/25-s).

The second derivative of 1-A_s is

    2[2alpha r+beta+r(alpha r+2beta)y]/(1-r y)^4.

The numerator inside square brackets is affine in y. At y=0 it is
`(27/125)(1-s)>0`; at y=1 it is

                    5076/15625-(7020/15625)s>0          (3.1)

for s<=2/3. This proves concavity. Also A_s(0)=1 and
`A_s(1)=(1+s)/(1+t)^2`, so throughout the interval

                              A_s(y)>=25/64.            (3.2)

For B_s put d=r-2s t. Its second derivative is

    2r[5r-4s t+r(r-2s t)y]/(1-r y)^4.

Again the numerator inside square brackets is affine in y. At y=0
it is `9/5-(12/5)s>0`; at y=1 it is

                              1206/625-(1770/625)s>0.  (3.3)

Thus Jensen for A and the endpoint chord for B prove the ACTUAL-measure
bound

    integral A_s dnu<=A_s(m),
    integral B_s dnu<=1+m[B_s(1)-1]
                    =1+m(297-375s)/128.                (3.4)

The two separate bounds do not assume that one measure attains both.
Their simultaneous use is valid for every actual probability measure.

Combining (2.2) and (3.4) yields the evaluated algebraic expression

    W_flat(3/5)/(2n^(3/2))
      <=sqrt(w) {(3/5)sqrt((1-u)A_(a u)(m))
               +(2/5)sqrt(kappa)
                        sqrt(1+m(297-375a u)/128)},     (3.5)

whenever 0<=a u<=2/3. In a small-gap active face one can substitute
`m=u^2/f_n^2+o(1)` at this fixed metric, using continuity on the compact
support. Formula (3.5) retains general f and w. It is not automatically
at most f/2, or at most sqrt(2)Phi(A)/n^(3/2), when f<sqrt(2).

## 4. Uniform diagnostic evaluation at f=sqrt(2)

Take w=1, a=kappa, m=u^2/2, and kappa<=u<=1. Define the right side of
(3.5) by

    U(u)=(3/5)sqrt((1-u)A_(kappa u)(u^2/2))
        +(2/5)sqrt(kappa)
                    sqrt(1+u^2(297-375kappa u)/256).    (4.1)

This is a continuous function on [kappa,1], and it is strictly
decreasing. The following elementary argument avoids optimization.

Let C(u)=A_(kappa u)(u^2/2). Then

    1-C(u)=G(u)/(1-9u^2/50)^2,
    G(u)=3u^2/50+kappa u^3/10+81u^4/2500-9kappa u^5/100.

The numerator determining the sign of (1-C)' after multiplication
by the positive denominator is

    u[3/25+(3kappa/10)u+(189/1250)u^2
                        -(54kappa/125)u^3
                        +(81kappa/5000)u^5].           (4.2)

For 0<u<=1, the two possibly competing kappa terms are at least
`-(33kappa/250)u`. Since 3/5<kappa<2/3, the bracket is strictly
positive: already `3/25-33(2/3)/250>0`. Thus C'(u)<0.
Together with C(u)>=25/64, the derivative of the first term of (4.1)
is at most -3/16 for u<1.

Put V(u)=1+u^2(297-375kappa u)/256. It is at least one, because
297-375kappa u>0. Also

    V'(u)=u(594-1125kappa u)/256<=189/256                (4.3)

on u>=kappa: use kappa u>=kappa^2>9/25 and u<=1 whenever the
bracket is positive; a negative bracket satisfies the bound directly.
The derivative of the second term of (4.1) is consequently at most
189/1280, using sqrt(kappa)<1 and sqrt(V)>=1. Since

                              3/16>189/1280,

the sum is strictly decreasing, as claimed. Continuity covers u=1.
Thus every actual measure in this diagnostic range satisfies

                              U(u)<=U(kappa).           (4.4)

## 5. Predetermined exact rational enclosure

The fixed certificate is
`original_mo_small_gap_pure_cross_fraction_certificate.py`.
Its frozen 81-line source has SHA256
`10d76c46fbdf75d8b856d06bae07a3d6304c78ce2d5b17de225567435f63fdf2`.
It contains no numerical optimization, signing construction, SDP,
search, simulation, or parameter scan. It reuses the already-proved
interval

                         3.1415926<pi<3.1415927,

from baseline checker SHA256
`d3af3d3bac9ba4d73a7589ba9ed4ff6261fde3263c64d04de36da7f36a1c65d3`
and baseline result SHA256
`fbc10c4760d963f9364dca586cca3d8df5692ab786cd155634a651fac3a62d9d`.
It does not rerun that prior Machin calculation.

Let k_- and k_+ be twice the inverse endpoints of this pi interval,
in increasing order. Put s_-=k_-^2, s_+=k_+^2 and m_-=s_-/2,
m_+=s_+/2. Positive rational upper bounds for the two means at u=kappa
are

    A_up=[1-(21/25+s_-/5)m_-+(9/25)s_+m_+^2]
                                              /(1-(9/25)m_+)^2,
    B_up=1+m_+(297-375s_-)/128.

The certificate checks by exact rational squaring that

    (9/25)(1-k_-)A_up < (35317/100000)^2,
    (4/25)k_+ B_up    < (35391/100000)^2.                (5.1)

It also checks the elementary strict comparison

    (35317+35391)/100000=70708/100000,
                    (70708/100000)^2<1/2.              (5.2)

Therefore the certified conclusion is

    U(u)<70708/100000<1/sqrt(2),
    W_flat(3/5)<(70708/50000)n^(3/2)<sqrt(2)n^(3/2).    (5.3)

Root ran these eleven new rational comparisons exactly once on the
authorized offload host; all passed and the process exited zero.
The result `/tmp/original-mo-pure-cross-rational.w3EsHK/result.json`
has SHA256
`0ea064435322e698b8e33a4d9bce8ab29156e3cfe013c9885f1f35e205156e41`.
It records the exact squared margin
`1/2-(17677/25000)^2=23671/625000000`. The author has read this
result and verified its hash but has run no mathematical computation
on the local coordinator. The reused baseline was not rerun.

## 6. Actual small-gap application and exact limits

The separate two-phase full-K rounding theorem supplies, under the
fixed norm cap and gamma tending to zero,

                             Phi(K)>=kappa S/2-o(n^(3/2)).

Its complete 274-line source is
`original_mo_full_sdp_gap_original_phase_bound.md`, SHA256
`1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`.

It is needed here to infer u>=kappa-o(1) when the ORIGINAL pure-cross
face is active: p=q_A=0 and c=Phi(K), with the sign chosen so c>=0.
This implication is not inferred merely from D+-K>=0 or from a norm cap.

If in addition Phi(K)/n^(3/2) tends to sqrt(2), equations (2.5)--(2.6)
give m=u^2/2+o(1). Fixed t=3/5 makes the functional uniformly continuous
on the relevant compact domain, so the small discrepancies in u and m,
the positive-field comparison error, and the actual-diagonal metric
comparison error are all o(n^(3/2)). The resulting actual pure-cross
cell bound is

    limsup [E max_(z in that cell) X_z]/(2n^(3/2))
                              <=70708/100000<1/sqrt(2). (6.1)

This is a strict, evaluated upper on the actual singular-measure
functional in the diagnostic small-gap range. It does not select a
formal Dirac measure or silently vary the actual cross signing.

For other f, the general expression (3.5), with m=u^2/f^2+o(1), remains
the stated result. The desired original conditional comparison is
F<=2sqrt(2)Phi(A); it is not replaced by F<=sqrt(2)n^(3/2). Neither
the other original internal-energy cells nor that smaller-alpha target
is settled here. The complementary large-gap case is separate work.
Original MO convergence is not claimed.
