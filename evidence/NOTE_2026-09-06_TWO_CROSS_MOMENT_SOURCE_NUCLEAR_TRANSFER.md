# Two actual cross moments control the original source nuclear cap

2026-09-06. Conditional analytic theorem for ACTUAL complete paired
signings. No mathematical computation, scan, solver, or search was run.
No full cross or internal spectral law is assumed.

## 1. Hypotheses and conclusions

Let n tend to infinity. Let K=[[A,B],[B^T,-A]] be a complete symmetric
zero-diagonal signing of order 2n, and let D=diag(D_L,D_R)>0 be diagonal
with D+-K>=0. Define the ACTUAL weighted blocks and diagonal dispersion

    dbar=tr(D)/(2n),       delta=tr(D)tr(D^(-1))/(2n)^2-1,
    H=D_L^(-1/2) A D_L^(-1/2),
    W=D_L^(-1/2) B D_R^(-1/2),       Y=WW^T,
    m_D=tr Y/n,       Delta_D=tr[Y(I-Y)]/n.

Assume only delta->0, m_D->m>0, and Delta_D->Delta. Then

    0<m<=1/2,       0<=Delta<=m(1-m).

Put C_m=1/sqrt(m) and define the continuous nuclear-cap function

    L(m,Delta)=sqrt[(1-m)(1-Delta/m)]+sqrt(Delta),
                                            if 0<=Delta<m^2,
    L(m,Delta)=1,                            if Delta>=m^2. (1.1)

There is one common original index set J, of order q with q/n->1,
for which the SAME actual source signing A_J satisfies

    dbar/sqrt(q)->C_m,
    limsup ||A_J||op/sqrt(q)<=C_m,
    limsup tr|A_J|/q^(3/2)<=L(m,Delta),
    Phi(A_J)/q^(3/2)<=Phi(A)/n^(3/2)+o(1),                 (1.2)

where Phi(A)=max_(x in {+-1}^n)|x^T A x|/2. No diagonal optimality,
separate trace cap, active-state premise, or internal shape is needed.
The full actual weighted internal moments also satisfy

    tr H^2/n->m,
    limsup tr|H|/n<=sqrt(m)L(m,Delta).                     (1.3)

Using the separately established all-law source theorem stated in
Section 5, these facts imply the functional original-source bound

    liminf Phi(A)/n^(3/2)>=F_(C_m)(L(m,Delta)).             (1.4)

The exact region where these CERTIFIED caps give C_m<=5/3 and L<=4/5 is

    R={9/25<=m<=1/2, 0<=Delta<=Delta_crit(m)},
    Delta_crit(m)=m[4sqrt(m)-3sqrt(1-m)]^2/25.              (1.5)

Throughout R, the original-source lower is at least 2/5+7/55000.
This is a sufficient conditional exclusion, not a necessary condition
for actual source exclusion and not a proof of the global MO limit.

## 2. Finite noncommuting Schatten bound

More generally suppose H is symmetric and 0<=H^2<=I-Y, with 0<=Y<=I.
Write M=tr H^2/n, v=tr Y/n, d=tr[Y(I-Y)]/n, and z=tr(H^2Y)/n.
Trace positivity gives 0<=z<=min(M,d). The Schatten Cauchy--Schwarz
inequality applied to H=H(I-Y)+HY gives, without any commutation,

    tr|H|/n<=sqrt[(M-z)(1-v)]+sqrt(z v).                  (2.1)

Indeed factor the summands as (H sqrt(I-Y))sqrt(I-Y) and
(H sqrt(Y))sqrt(Y). Their squared Frobenius factors are respectively
n(M-z), n(1-v), nz, and nv. The trace-norm triangle inequality proves
(2.1). The right side is concave in z and maximized at z=vM before
the constraint z<=d. Consequently, with z_*=min(d,vM),

    tr|H|/n<=S(M,v,d)
       :=sqrt[(M-z_*)(1-v)]+sqrt(z_* v).                  (2.2)

This formula includes the boundary cases by continuity. For d>=vM
it is sqrt(M); otherwise it is sqrt[(M-d)(1-v)]+sqrt(d v).
The matrix inequality H^2<=I-Y yields z<=d because the trace of
[(I-Y)-H^2]Y is nonnegative. No simultaneous diagonalization is used.

For the actual paired blocks, T=D^(-1/2)KD^(-1/2) is a contraction.
Its upper block row therefore gives H^2+Y<=I, so (2.2) applies to
the FULL actual H and Y, not substituted spectral surrogates.

## 3. Near-scalar scale and one common source compression

The following scale and common-label argument restates the endpoint
transfer's mechanism, but now uses only convergence of m_D.
Set t_i=d_i/dbar, m_0=n/dbar^2, and

    e_L=(1/n)sum_(i=1)^n 1/t_i,
    e_R=(1/n)sum_(i=1)^n 1/t_(n+i).

Literal complete cross sign squares give m_D=m_0 e_L e_R, while
e_L+e_R=2(1+delta). The arithmetic means of t in the two halves sum
to 2. Cauchy--Schwarz in each half therefore gives

    1/e_L+1/e_R<=2,
    1+delta<=e_L e_R<=(1+delta)^2,
    m_D/(1+delta)^2<=m_0<=m_D/(1+delta).                   (3.1)

Thus m_0->m and dbar/sqrt(n)->C_m without a separate trace cap.
The exact dispersion identity is

                 (1/(2n))sum_i (t_i-1)^2/t_i=delta.

For delta>0 take epsilon=delta^(1/3), eventually at most 1/8, and
retain the common labels J with both |t_i-1|<=epsilon and
|t_(n+i)-1|<=epsilon. Write a=q/n and b=1-a. Since outside this
interval (t-1)^2/t>=epsilon^2/(1+epsilon),

    b<=2delta(1+epsilon)/epsilon^2->0.                    (3.2)

For delta=0 take epsilon=0 and all labels. This convention also
handles mixed zero/positive-dispersion sequences.
Let H_J=H[J,J] and Q=diag(sqrt(t_i):i in J). Principal feasibility
and exact diagonal congruence give

    ||H||op, ||H_J||op<=1,       A_J/dbar=Q H_J Q,
    ||A_J||op<=(1+epsilon)dbar,
    ||A_J/dbar-H_J||op<=3epsilon.                          (3.3)

Here ||Q-I||<=epsilon and ||Q||<=sqrt(1+epsilon) prove the last bound.
The scale and operator assertions in (1.2) follow.
Put M_H=tr H^2/n and M_J=tr H_J^2/q. Interlacing on [-1,1] bounds
the loss of each positive and negative squared-eigenvalue sum by n-q.
Thus 0<=M_H-aM_J<=2b and |M_H-M_J|<=2b. The square function is
3-Lipschitz on the common interval [-(1+epsilon),1+epsilon]. Weyl's
inequality in (3.3), and the exact completeness identity
tr A_J^2=q(q-1), now give

    |M_H-(q-1)/dbar^2|<=2b+9epsilon.                      (3.4)

It follows that M_H->m. Since M_H+m_D<=1, necessarily m<=1/2.
Also 0<=Delta_D<=m_D(1-m_D), by 0<=Y<=I and the scalar second-moment
inequality. This proves the asserted range of Delta. Continuity of
(2.2) yields

    limsup tr|H|/n<=S(m,m,Delta)=sqrt(m)L(m,Delta).

Principal compression decreases the UNNORMALIZED trace norm. Weyl
comparison with (3.3) then gives

    tr|A_J/dbar|/q<=a^(-1)tr|H|/n+3epsilon.

Multiply by dbar/sqrt(q)->C_m to obtain the nuclear claim in (1.2).
For every fixed Boolean signing on J, independent unbiased extension
to the removed labels leaves the expected original energy unchanged.
Hence Phi(A_J)<=Phi(A). Finally feasibility gives
Phi(A)<=tr(D_L)/2<=n dbar, so Phi(A)/n^(3/2) remains bounded and

    Phi(A_J)/q^(3/2)<=a^(-3/2)Phi(A)/n^(3/2)
                          =Phi(A)/n^(3/2)+o(1).

The compression is used only for this ORIGINAL norm lower. It does
not replace the paired covariance, full W, active field, or optimizer.

## 4. Exact cap-threshold curve and a simple interior rectangle

For fixed 0<m<=1/2, L(m,Delta) increases from sqrt(1-m) to 1 as
Delta runs from 0 to m^2, and then stays 1. Put

    beta=arcsin(sqrt(m)),       theta=arcsin(sqrt(Delta/m))

on 0<=Delta<=m^2. Then 0<=theta<=beta and L=cos(beta-theta).
Consequently L<=4/5 is possible exactly when m>=9/25. In that range
it holds exactly when

    sqrt(Delta/m)<=[4sqrt(m)-3sqrt(1-m)]/5,

which proves (1.5); C_m<=5/3 gives the same lower threshold on m.
This is exact only for the two displayed certified caps, not a claim
that all actual sources outside R fail an exclusion. In particular
Delta_crit(9/25)=0 and Delta_crit(1/2)=1/100.

A convenient sufficient rectangle inside R is

                 2/5<=m<=1/2,       0<=Delta<=1/1600.     (4.1)

On Delta<m^2, the first squared term in L equals
1-m-Delta/m+Delta; its m derivative is -1+Delta/m^2<0.
Thus L decreases in m and increases in Delta in this rectangle. At
its worst corner,

    L(2/5,1/1600)=sqrt(1917/3200)+1/40<4/5,

because (31/40)^2=1922/3200. Also C_(2/5)<5/3.

## 5. Monotone composition with the separately proved all-law gain

For kappa=2/pi and a,w>=0 define

    Psi(a,w)=[(a+w)arctan(sqrt(w/a))-sqrt(a w)]/pi, a>0,
    Psi(0,w)=w/2,       Psi(a,0)=0,
    F_C(u)=kappa/(2u)+(1-exp(-1))/(2C)
                          Psi(kappa C^2,(1-kappa)u^2).

The frozen 553-line all-law theorem proves, uniformly for actual
complete source signings with ||A||op/sqrt(n)<=C, C>=1,

    Phi(A)/n^(3/2)>=F_C(tr|A|/n^(3/2))-e_C(n),
    e_C(n)->0.                                            (5.1)

It also proves the consequence 2/5+7/55000 under the two limsup
caps C<=5/3 and nuclear moment <=4/5. This is an imported, separately
proved theorem, not an additional unproved premise about our sources.

The new composition needs monotonicity in the retained nuclear cap.
For a>0, partial_w Psi(a,w)=arctan(sqrt(w/a))/pi<=1/2. Therefore for
C>=1 and 0<u<=1,

    F_C'(u)<=-kappa/(2u^2)
                    +(1-exp(-1))(1-kappa)u/(2C)<0,        (5.2)

using kappa>1/2. On any subsequence relevant to the original liminf,
extract a further subsequence where u_J=tr|A_J|/q^(3/2) converges.
Completeness and the operator cap bound u_J away from zero:
u_J>=[(q-1)/q]/(||A_J||op/sqrt(q)); also u_J<=sqrt((q-1)/q).
Its limit u thus lies in (0,L(m,Delta)]. Apply (5.1) with any fixed
C=C_m+eta, eta>0, then let eta decrease to zero. Continuity and
(5.2) yield the bound F_(C_m)(u)>=F_(C_m)(L(m,Delta)). The original
norm comparison in (1.2) proves (1.4). For (m,Delta) in R, the
separately proved two-cap consequence gives 2/5+7/55000 directly.

The constant region also has a no-moment-convergence formulation.
If delta->0 and EVERY accumulation point of (m_D,Delta_D) belongs
to the compact region R, then the same original liminf is at least
2/5+7/55000. Indeed every subsequence has a further moment-convergent
subsequence in R, where the preceding proof applies. A simple
sufficient condition is liminf m_D>=2/5 and
limsup Delta_D<=1/1600: (3.1)--(3.4) force all its moment limits to
have m<=1/2, so (4.1) applies. No empirical spectral law is required.

## 6. Frozen prerequisites, credit, and scope

The common-good-label scale transfer is restated from the 230-line
`original_mo_near_scalar_cross_endpoint_source_nuclear_transfer.md`,
SHA256 `6a486df0fd46aa76259e3f02e3734eb2529162500f98f89af58e90562e6a2187`.
The all-law theorem imported in Section 5 is
`original_mo_all_law_adaptive_nuclear_gain.md`, 553 lines,
SHA256 `0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51`.
Both were read fully and remain frozen, unchanged by this work.

Root contributed the finite two-moment Schatten inequality. The exact
worker independently checked it, derived the threshold inversion,
monotone composition, and source transfer, and authored this note.
Root and docs supplied contributing or pre-writing checks of new
links. Those checks do not replace a full frozen-source review, and
authorship is not claimed as independent whole-proof verification.

Only this /tmp source and its author receipt are written by its author.
No canonical source, publication, or backup was changed or performed.
The proof supplies actual-source consequences of near-scalar diagonal
feasibility and two cross moments. It does not show those hypotheses
for arbitrary original minimizers or settle complementary active cells.
The original MO convergence question remains OPEN.
