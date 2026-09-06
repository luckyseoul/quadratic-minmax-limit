# All-law original-source nuclear gain from adaptive Boolean updates

2026-09-06. Analytic theorem for ACTUAL complete symmetric zero-diagonal
signings. No mathematical program, checker, spectral scan, numerical
integral, optimization run, signing construction, or search was executed.

This theorem uses the actual coordinate-normalized |A|+/-A phases and
coordinate-dependent Boolean updates. It imposes NO limiting spectral
law, spectral flatness, spectral symmetry, diagonal homogeneity, or
near-scalar SDP premise. Its separate normalized operator cap is explicit.
All energies are the ORIGINAL same-source quadratic energies.

## 1. Uniform all-law theorem and consequential excluded region

Let A be a complete symmetric zero-diagonal signing of order n>=2, and set

    Q_A(x)=x^T A x/2,
    Phi(A)=max_(x in {+1,-1}^n)|Q_A(x)|,
    alpha(A)=Phi(A)/n^(3/2),
    kappa=2/pi,
    M=A/sqrt(n),       q=(n-1)/n,
    ell=tr|M|/n=tr|A|/n^(3/2).

Fix C>=1 and suppose the ACTUAL operator norm satisfies ||M||op<=C.
For a,w>=0 define the continuous function

    Psi(a,w)=[(a+w)arctan(sqrt(w/a))-sqrt(a w)]/pi,
                      when a>0,
    Psi(0,w)=w/2,       Psi(a,0)=0.                         (1.1)

There is a function e_C(n) tending to zero such that, uniformly over
ALL these actual A,

    alpha(A)>=kappa/(2ell)+D_C(ell)-e_C(n),
    D_C(ell)=(1-exp(-1))/(2C)
                         Psi(kappa C^2,(1-kappa)ell^2).     (1.2)

In particular the explicit elementary lower envelope is

    alpha(A)>=kappa/(2ell)+J_C(ell)-e_C(n),
    J_C(ell)=5[kappa(1-kappa)]^(3/2)ell^3
                    /[48(kappa C^2+(1-kappa)ell^2)].        (1.3)

The error function may be enlarged between formulas. No finite-n rate
for the Gaussianization error is claimed. For fixed C and ell>0,
both displayed gain functions are strictly positive.

Consequently, for ANY sequence of actual complete source signings with

    limsup ||A||op/sqrt(n)<=5/3,
    limsup tr|A|/n^(3/2)<=4/5,                              (1.4)

one has

    liminf Phi(A)/n^(3/2)
       >=35/88+3/1250
        =2/5+7/55000.                                      (1.5)

Thus the entire actual bounded-operator/nuclear-moment region (1.4)
is incompatible with alpha(A) tending to 2/5. This is not restricted
to one three-atom profile; the empirical spectral laws need not even
converge. Neither an arbitrary norm cap nor original near-optimality is
asserted to imply the operator cap in (1.4). The global MO limit and
the remaining all-profile/all-active-cell implication remain OPEN.

## 2. Existing-object audit and what is new

The previously proved nuclear inequality uses the same |A|+/-A phases
and gives alpha(A)>=kappa q/(2ell), without an operator cap. It is
recorded in original_mo_nuclear_spectral_budget.md and in Section 4
of original_mo_original_phase_spectral_moment.md.

The separate original_mo_original_source_near_flat_strict_gain.md
proves a stronger numerical gap at one specified near-flat spectral law.
Its distinguished-coordinate Gaussianization proof is reused below,
but its projector approximation and first-chaos eigenvalue alignment
are NOT assumed here. The present changes are:

1. Actual, potentially unequal |M| diagonals still give a bounded frame
   under the explicit operator cap, and a NEW uniform higher-chaos
   mean lower follows from |M|=H^(1/2)R_s H^(1/2)-sM.
2. The update probability depends on each ACTUAL local field. A radial
   Gaussian argument controls its clipping cost for every covariance.
3. A convex perspective averages the resulting negative-field second
   moments without first-chaos alignment or homogeneous local variances.

No old variance-only claim, scalar moment scan, or finite profile hunt
is being rerun or relabelled. The new inequality depends only on the
actual source's retained nuclear moment and its separately given cap.

## 3. Coordinate-normalized actual spectral phases

Put L=|M|, h_i=L_ii, H=diag(h_i), and ell=(1/n)sum_i h_i.
Completeness gives M_ii=0 and (M^2)_ii=q. Spectral calculus and
Cauchy--Schwarz in each coordinate's spectral measure give

    M^2<=C L,
    q/C<=h_i<=sqrt(q),
    q/C<=ell<=sqrt(q).                                     (3.1)

In particular every h_i is positive, even though the h_i need not be
close to their mean. For s in {+1,-1} define

    R_s=H^(-1/2)(L+sM)H^(-1/2),
    K_n=2C^2/q.

These are genuine Gaussian correlation matrices, since L+sM>=0 and
diag(L+sM)=h. Their operator norms obey

    diag R_s=1,       0<=R_s<=K_n I,       K_n<=4C^2.        (3.2)

Indeed ||L+sM||op<=2C and ||H^(-1/2)||op^2<=C/q.
This controls the actual frame WITHOUT replacing H by ell I.

Let G_s be a centered Gaussian vector with covariance R_s, and put
X_s=sign(G_s). The two phases may be sampled separately; no joint
coupling between them is used. Gaussian singularity is allowed.
Use sign(0)=+1 throughout.

## 4. Exact original nuclear baseline of these same phases

Define the two oriented original-phase energies

    E_s=E Q_(sA)(X_s)/n^(3/2).

The common diagonal h makes the pair-correlation difference exact:

    (R_+)_ij-(R_-)_ij=2M_ij/sqrt(h_i h_j).

The Gaussian sign identity is E sign(U)sign(V)=kappa arcsin r.
Arcsine has derivative at least one on (-1,1), with its continuous
endpoint extension. Since M_ij=A_ij/sqrt(n) off the diagonal,

    (E_++E_-)/2
       >=(kappa/n^2)sum_(i<j)1/sqrt(h_i h_j)
       >=kappa q/(2ell).                                   (4.1)

For the last inequality, write k=n(n-1)/2 and use

    sum_(i<j)1/sqrt(h_i h_j)
       >=2sum_(i<j)1/(h_i+h_j)
       >=2k^2/[(n-1)sum_i h_i]
        =n(n-1)/(2ell).

This is an exact lower on the average of the two ACTUAL phase energies,
not merely an already maximized nuclear-norm inequality. That distinction
allows us to improve the same phases and add the new gain to (4.1).

## 5. All-law higher-chaos noise lower, with actual diagonals retained

Fix either phase and abbreviate R=R_s. For every integer k>=1,

    R^{circ k}>=0,       diag(R^{circ k})=1,
    ||R^{circ k}||op<=||R||op<=K_n.                          (5.1)

This follows by iterating the positive Schur multiplier associated with
a correlation matrix; it preserves PSD order and fixes scalar diagonals.
For every odd k>=3, actual complete entries give

    |tr(M R^{circ k})|
       <=(1/sqrt(n))sum_(i!=j)|R_ij|^k
       <=tr(R^2)/sqrt(n)<=K_n sqrt(n).                      (5.2)

The EXACT matrix identity L=H^(1/2)R H^(1/2)-sM yields

    tr(L R^{circ k})
      =sum_(i,j)sqrt(h_i h_j)R_ij^(k+1)
                                      -s tr(M R^{circ k})
      >=n ell-K_n sqrt(n),                                 (5.3)

because k+1 is even and its diagonal terms sum to n ell. In particular
the noncommutation of H with M causes no lost or unbounded error here.

For every real t>=0, scalar spectral calculus gives

    M^2=L^2>=t L-(t^2/4)I.

Pair with the PSD matrix R^{circ k}, use (5.3), and optimize this
scalar quadratic by taking t=2(ell-K_n/sqrt(n))_+. This gives

    (1/n)tr(M^2 R^{circ k})
       >=(ell-K_n/sqrt(n))_+^2,
                    uniformly for every odd k>=3.           (5.4)

Let p_k denote the normalized Gaussian Hermite polynomials, and write
sign(z)=sum_(k odd>=1)c_k p_k(z), with

    c_1=sqrt(kappa),       sum_k c_k^2=1.

The sign covariance and its higher-chaos tail are

    C_X=E[X_s X_s^T]=kappa R+C_tail,
    C_tail=sum_(k odd>=3)c_k^2 R^{circ k},
    0<=C_tail<=(1-kappa)K_n I.                              (5.5)

Operator convergence of this series follows from (5.1) and the scalar
coefficient tails. Define local fields, first-chaos fields, and variances

    F_i=(sM X_s)_i,
    F_(i,1)=sqrt(kappa)(sM G_s)_i,
    v_i=E(F_i-F_(i,1))^2=(M C_tail M)_ii,
    b_i=E F_(i,1)^2,
    c_i=E[(G_s)_i F_i]=E[(G_s)_i F_(i,1)],
    sigma_i^2=E F_i^2=b_i+v_i,
    w_i=sigma_i^2-c_i^2>=v_i.                              (5.6)

The equalities and inequality use orthogonality of different Gaussian
chaoses and Cauchy--Schwarz in the first chaos. Summing (5.4) with the
nonnegative tail coefficients proves the finite bound

    (1/n)sum_i w_i>=(1/n)sum_i v_i
      >=(1-kappa)(ell-K_n/sqrt(n))_+^2
       =(1-kappa)ell^2-O_C(n^(-1/2)).                       (5.7)

No error depending on the chaos order was inserted into an infinite
series. This lower holds for BOTH phases separately, regardless of
their own original-energy baselines.

We also have the crucial finite bounds

    sigma_i^2=(M C_X M)_ii<=K_n q=2C^2,
    (1/n)sum_i c_i^2<=(1/n)sum_i b_i
       =(kappa/n)tr(M^2 R)<=kappa C^2.                      (5.8)

The first uses the exact row-square sum q and C_X<=K_n I. The second
uses M^2<=C^2 I and tr R=n. It does not assume c_i is positive or
close to a common value. This replaces the near-flat alignment premise.

## 6. Uniform Gaussianization for the clipped adaptive objective

For a correlation matrix R with bounded operator norm, a Gaussian G
of covariance R, coefficients d_j with max|d_j|<=1/sqrt(n) and
sum_j d_j^2<=1, and any chosen i, put F=sum_j d_j sign(G_j).
The complete proof in Section 6 of the frozen near-flat source theorem
gives the following JOINT marginal assertion: along any subsequence
where the covariance matrix of (G_i,F) converges, this pair converges
in distribution to the centered Gaussian pair with that covariance.
The approximation is uniform over the admissible inputs by compactness.

Here is its proof mechanism, including the relevant coefficient and
distinguished-coordinate checks. Represent G_j by unit Gram vectors u_j
with Gram matrix R. The unscaled chaos kernels are
g_k=sum_j d_j u_j^{tensor k}. Positive Schur powers have a common
operator bound B. For a contraction of r indices, with
1<=r<min(p,k), its squared norm is

    tr(T R^{circ(k-r)} T R^{circ(p-r)})<=B^4/n,
    T=diag(d)R^{circ r}diag(d).

For a full contraction of the smaller, unequal order, writing
z=R^{circ p}d gives the bound

    z^T diag(d)R^{circ(k-p)}diag(d)z<=B^3/n.

The additional distinguished-coordinate contraction satisfies

    ||u_i tensor_1 g_k||^2
       <=(B/n)sum_j R_ij^2<=B^2/n,       k>=3.

For a fixed finite Hermite truncation F_Q, every nonconstant contraction
of the linear combination aG_i+bF_Q is covered by these estimates.
The equal-order full contractions are exactly constant covariances.
Gaussian polynomial integration by parts, with U obtained by dividing
each chaos component by its order, gives

    Gamma=grad(aG_i+bF_Q) dot grad U,
    E Gamma=Var(aG_i+bF_Q),       Var Gamma=O_(Q,B,a,b)(1/n).

The resulting characteristic-function equation has Gaussian error at
argument z bounded by z^2 sqrt(Var Gamma)/2, without inverse variance.
The omitted sign tail has L2 norm at most sqrt(B tau_Q), where
tau_Q=sum_(k>Q)c_k^2->0, and is orthogonal to both G_i and all retained
chaoses. Taking n first and then Q proves the claim for every fixed
linear combination. Compact covariance subsequences and Cramer--Wold
give the stated two-variable joint limit, including singular limits.
The full identities and normalization are proved in the cited source;
no new external limit theorem or growing-dimensional CLT is invoked.

For C>0 and r>=0 set

    g_C(r)=r^2/(2C)-(r-2C)_+^2/(2C)
          =r^2/(2C)                  if 0<=r<=2C,
           2r-2C                    if r>=2C.              (6.1)

This is continuous, nonnegative, and at most 2r. Apply the joint limit
to the function g_C((-sign(u)v)_+). Its only possible discontinuities
lie on u=0, a set of zero probability under every limiting Gaussian
pair because its first marginal is standard. Its linear growth is
uniformly integrable under the bounded second moments of F. Therefore,
uniformly over the above data,

    E g_C((-sign(G_i)F)_+)
      =E g_C((-sign(Z)(cZ+sqrt(w)Z'))_+)+o_(B,C)(1),        (6.2)

where Z,Z' are independent standard Gaussians,
c=E[G_i F], and w=E F^2-c^2>=0. Degenerate second marginals cause
no difficulty: the function vanishes at v=0. Uniformity follows by
choosing a violating subsequence, extracting its compact covariance
limit, and applying the preceding convergence and uniform integrability.
No uniform integrability of F^2 is inferred from its bounded expectation;
only the linearly growing clipped function is passed through this limit.

The actual row coefficients d_j=sA_ij/sqrt(n), including their one
zero diagonal coefficient, satisfy both required bounds. Thus (6.2)
applies uniformly to all n rows of BOTH phases with B=4C^2.

## 7. Exact coordinate-adaptive update of the original Boolean source

Fix either phase s and condition on X=X_s. Let

    F=sM X,       Y=sign(F),       Delta=Y-X,
    r_i=(-X_i F_i)_+,
    epsilon_i=min(r_i/(2C),1).

Conditionally on X, choose INDEPENDENT Bernoulli(epsilon_i) variables
xi_i and set X'_i=X_i+xi_i Delta_i. Every X' is Boolean on the same
original source. With z_i=epsilon_i Delta_i, zero diagonal gives exactly

    E_xi Q_(sM)(X')
       =Q_(sM)(X)+(sM X)^T z+z^T(sM)z/2.

There is no diagonal correction because (sM)_ii=0. No independence of
the original signs or of the adaptive probabilities is asserted.

On a mismatch, F_i Delta_i=2r_i and Delta_i^2=4. On a match,
Delta_i=0. A zero field has r_i=epsilon_i=0 and causes no exception.
Since ||sM||op<=C,

    E_xi Q_(sM)(X')-Q_(sM)(X)
       >=sum_i[2epsilon_i r_i-2C epsilon_i^2]
        =sum_i g_C(r_i).                                   (7.1)

The equality uses the chosen clipped optimizer epsilon_i, which is
always in [0,1]. It is an admissible actual update, not an unconstrained
scalar optimizer. Since Q_(sA)(X')<=Phi(A) for every realization,

    alpha(A)>=E_s+(1/n)sum_i E g_C(r_i).                    (7.2)

Thus the task is to lower-bound the new gain for each phase, retaining
the actual operator penalty C and the local joint Gaussian covariance.

## 8. Uniform Gaussian clipping control by radial integration

Let (Z,F) be ANY centered Gaussian pair with Var Z=1 and
sigma^2=Var F<=2C^2. Put r=(-sign(Z)F)_+. Represent the pair as
linear forms of two independent standard Gaussians, permitting a
rank-one or zero second form. In polar coordinates their radius T
is independent of their uniform angle and has density t exp(-t^2/2).
The mismatch set depends only on the angle. At each mismatching angle,
r=bT with 0<=b<=sigma; at other angles r=0.

For a>=0 direct one-dimensional integration gives

    E(T-a)_+^2 / E T^2
       =exp(-a^2/2)-a integral_a^infinity exp(-u^2/2)du
       <=exp(-a^2/2),       E T^2=2.                       (8.1)

For example, expand (t-a)^2 t and integrate its three terms; integration
by parts gives the displayed expression. No numerical integral is used.
At each angle with b>0, take a=2C/b. Then

    E[(bT-2C)_+^2 | angle]
       <=exp(-2C^2/b^2) E[(bT)^2 | angle]
       <=exp(-1) E[(bT)^2 | angle].                        (8.2)

The last inequality uses b^2<=sigma^2<=2C^2. Zero b and nonmismatch
angles contribute zero. Averaging over angles and using (6.1) proves

    E g_C(r)>=(1-exp(-1))/(2C) E r^2.                      (8.3)

This works for negative, positive, or zero covariance, and for degenerate
Gaussian pairs. It does not divide by the mismatch probability and does
not discard a potentially large unquantified Gaussian tail.

## 9. A convex perspective eliminates the unknown local frame geometry

Write F=cZ+sqrt(w)Z' with independent standard Gaussians Z,Z' and w>=0.
If c<0, replacing c by |c| decreases r pointwise in the equivalent
representation

    r distributed as (-c|Z|+sqrt(w)Z')_+.

For c>=0, a planar Gaussian integral gives

    E r^2=Psi(c^2,w)
       =[(c^2+w)arctan(sqrt(w)/c)-c sqrt(w)]/pi.            (9.1)

Explicitly, set theta=arctan(sqrt(w)/c) and sigma^2=c^2+w.
The two angular mismatch sectors each have length theta. With the
radial second moment two, their integral is

    (2sigma^2/pi)integral_0^theta sin(u)^2 du
       =(sigma^2/pi)[theta-sin(2theta)/2],

which equals (9.1) because sin(2theta)/2=c sqrt(w)/sigma^2.
At c=0 it equals w/2, and at w=0 it equals zero. Thus for every sign
of c,

    E r^2>=Psi(c^2,w).                                     (9.2)

We require the JOINT convexity of Psi in its two nonnegative arguments,
not convexity in w alone. For a,w>0 write

    Psi(a,w)=w f(a/w),
    f(x)=[(1+x)arctan(x^(-1/2))-sqrt(x)]/pi.

Direct differentiation gives

    f'(x)=[arctan(x^(-1/2))-x^(-1/2)]/pi,
    f''(x)=1/[2pi x^(3/2)(1+x)]>0.                         (9.3)

The perspective of this convex function is jointly convex: its Hessian
is f''(a/w)/w times the outer product of (1,-a/w) with itself.
The continuous extension in (1.1) preserves convexity on the closed
nonnegative quadrant. Moreover Psi decreases in a, since f'<=0, and
increases in w; for a>0 its w derivative is

    partial_w Psi(a,w)=arctan(sqrt(w/a))/pi>=0.              (9.4)

For either phase, (5.7)--(5.8) and Jensen therefore imply

    (1/n)sum_i Psi(c_i^2,w_i)
       >=Psi((1/n)sum_i c_i^2,(1/n)sum_i w_i)
       >=Psi(kappa C^2,
                    (1-kappa)(ell-K_n/sqrt(n))_+^2).        (9.5)

There is no assumption on the distribution, signs, homogeneity, or
convergence of the local parameters c_i and w_i. This is the missing
averaging step that a first-chaos-alignment argument would not supply
for arbitrary spectral laws.

Combining (6.2), (8.3), and (9.5) gives a common lower on each phase's
adaptive gain. Insert it in (7.2) for s=+1 and s=-1, and average the
two valid lower bounds on alpha(A). Equation (4.1) yields the finite-
parameter form

    alpha(A)>=kappa q/(2ell)
       +(1-exp(-1))/(2C)
          Psi(kappa C^2,(1-kappa)(ell-K_n/sqrt(n))_+^2)
       -o_C(1).                                           (9.6)

Uniformity follows from the rowwise uniform Gaussianization errors.
The ranges (3.1), K_n<=4C^2, and 0<=partial_w Psi<=1/2 show that
replacing q by one and the second Psi argument by (1-kappa)ell^2
costs only O_C(n^(-1/2)). This proves the uniform theorem (1.2).

## 10. Elementary gain envelope and the complete nuclear region

For c>0, integrate (9.4) from zero to w to obtain the useful identity

    Psi(c^2,w)=(1/pi)integral_0^w arctan(sqrt(t)/c)dt.

For x>=0, arctan x>=x/(1+x^2), as follows by differentiation from
equality at zero. Bounding the denominator t+c^2 above by w+c^2 gives

    Psi(c^2,w)>=2c w^(3/2)/[3pi(c^2+w)].                   (10.1)

The integral step matters: directly replacing the arctangent in the
closed form (9.1) by that lower bound would lose the useful positive term.

The elementary exponential series gives exp(1)>1+1+1/2+1/6=8/3,
so 1-exp(-1)>5/8. Taking c=C sqrt(kappa),
w=(1-kappa)ell^2 in (10.1), and using pi=2/kappa, proves
D_C(ell)>=J_C(ell) as displayed in (1.3).

Only the previously verified interval 7/11<kappa<16/25 is needed for
the following rational region. No pi certificate is rerun.
At C0=5/3 and 3/4<=ell<=4/5,

    [kappa(1-kappa)]^(3/2)>(12/25)^3,
    kappa C0^2+(1-kappa)ell^2
       <=16/9+144/625
        =11296/5625<81/40.                                (10.2)

The first uses that kappa(1-kappa) decreases for kappa>1/2.
For the second, the expression increases in ell and, at ell=4/5,
in kappa; evaluate at kappa=16/25. The final comparison is
451840<455625. Consequently

    J_(5/3)(ell)
       >=5(12/25)^3(3/4)^3/[48(81/40)]
        =3/1250,
    kappa/(2ell)>=5kappa/8>35/88.                          (10.3)

For 0<ell<=3/4 the ORIGINAL nuclear baseline alone gives

    kappa/(2ell)>=2kappa/3>14/33
       >35/88+3/1250.                                    (10.4)

Indeed 14/33-35/88=7/264>3/1250. Thus the smaller-nuclear-moment
part of the region is covered with ample room; a gain estimate is not
needed there. Finally

    35/88+3/1250=2/5+7/55000.                              (10.5)

To justify (1.5) with only the two limsup hypotheses in (1.4), no moment
convergence is imposed. On any subsequence relevant to the liminf,
extract a further subsequence where ell converges; (3.1) and the cap
place its limit in [3/5,4/5]. Apply the uniform theorem with any fixed
C=5/3+epsilon for all sufficiently large n, then let epsilon decrease
to zero. Continuity of the displayed functions returns C0=5/3.
Equations (10.3)--(10.4) then give (1.5) at every possible limit ell.
This also covers sequences alternating across ell=3/4 or approaching
4/5 from above by a vanishing error.

## 11. Provenance, independent-check boundaries, and remaining work

Frozen prerequisites directly read for this argument include:

    original_mo_original_source_near_flat_strict_gain.md, 612 lines
    SHA256 7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298

    original_mo_nuclear_spectral_budget.md
    SHA256 ee8ad5ff3dbf9aa9e251c4190e98ee1671c9a2140c759ba6f768f8c9c03ef13d

    original_mo_complete_cross_flat_spectral_gain.md, 411 lines
    SHA256 b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20

    original_mo_original_phase_spectral_moment.md, 262 lines
    SHA256 7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600

    original_mo_source_cross_nuclear_trace_boundary.md, 444 lines
    SHA256 106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556

The older scalar209 fixed-probability note is NOT a logical prerequisite
for the adaptive clipping or convexity arguments; those are derived here.
Its reusable coarse kappa interval comes from the listed older pi source.

Root posed the all-law actual-source extension and requested an entire
consequential bad-region exclusion. The proof worker and exact worker
independently derived the bounded normalized |M|+/-M frame and the
all-law higher-chaos mean. The proof worker proposed the field-adaptive
clipped update, derived the radial clipping bound, convex averaging,
and the explicit region, and wrote this full note. The exact worker
independently checked those proposals before writing and supplied the
finite uniform form (5.4), along with explicit checks of the adaptive
quadratic identity and scalar constants. These are contributing checks,
not a claim of whole-new-source independence by either contributor.
The older Gaussianization and nuclear prerequisites retain their own
disclosed provenance. A full new-source reviewer must state such roles.

Only this /tmp source and its review artifacts are new. The reviewed
612-line theorem and all earlier sources remain frozen and unchanged.
No canonical repository file was edited by this source's author, and
no mathematical computation was run. Publication, documentation gates,
commits, and backups are root's separate workflow.

The theorem removes the three-atom and diagonal-homogeneity premises,
not the ACTUAL bounded-operator premise. It neither proves that arbitrary
original minimizers satisfy that cap nor bounds the complementary
nuclear-moment region sufficiently for the paired all-cell comparison.
No source, covariance, cross block, or active state is silently replaced.
Original MO convergence and its possible limiting value remain OPEN.
