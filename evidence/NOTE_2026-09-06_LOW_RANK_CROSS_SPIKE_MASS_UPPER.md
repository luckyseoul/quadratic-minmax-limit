# Low-rank cross spikes: retain their state mass inside the Boolean upper

2026-09-06. General actual-cell upper and a weak-middle-law corollary.
No mathematical program, scan, solver, signing construction, or numerical
optimization was run. The original source, covariance, and cross energy
are retained. A low-rank Gaussian projection is paid for explicitly;
the Boolean state's mass in that subspace is NOT assumed small.

## 1. Actual objects and the reusable allocation upper

Let N=2n, K=[[A,B],[B^T,-A]] be a complete paired signing, and let
D>0 be diagonal with D+-K>=0. Set

    S=tr D,       dbar=S/N,       delta=S tr(D^(-1))/N^2-1,
    H_B=[[0,B],[B^T,0]],       L=D^(-1/2)H_BD^(-1/2).

Then ||L||op<=1 by block-sign conjugation and averaging. Let C be any
nonempty actual Boolean cell with a common ORIGINAL cross energy
c=x^TBy, and put u=2c/S in [-1,1]. Its positive Gaussian linear field
g has the weighted-shell covariance M, or more generally any covariance

    0<=M<=2wn I,       M_ii=wn,       0<=w<=1.              (1.1)

The w=0 case is immediate. Fix a BULK cutoff 0<r0<1 and define

    Pi=1_(|L|>r0),       Q=I-Pi,       d=rank(Pi).

This is not an operator-norm bound on the full actual L. Fix eta>=0
with eta r0<1, and put

    e=1-eta r0,
    F=Q[(I-eta L)|_(ran Q)]^(-1)Q,
    Mhat=D^(-1/2)M D^(-1/2),
    T=tr(Mhat F),       R=tr(Mhat F^2).                   (1.2)

F is zero on the spike space; no inverse is taken there. In particular
these definitions remain valid if I-eta L is singular on that space.
One has 0<=eF<=I and T-eR>=0. Define

    I_u=[max(0,(u-r0)/(1-r0)),1].

For sequences with delta->0 and d=o(n), the general actual-cell upper is

    E max_(z in C) g^Tz
      <=sqrt(S) max_(t in I_u) {
          sqrt(eta[r0-u+(1-r0)t](T-eR))
                      +sqrt(kappa)e sqrt(R)sqrt(1-t)}
          +o(n^(3/2)),       kappa=2/pi.                  (1.3)

The error is uniform over the fields and cells satisfying (1.1), at
fixed r0,eta. No convergence of T,R,u is required. The retained t is
the actual weighted spike-mass fraction. Formula (1.3) does not replace
u by a bulk cross energy or assert that t vanishes. Negative metric
orientations follow by replacing L,u by -L,-u in the metric calculation.

## 2. Project only the Gaussian field; retain the full state allocation

For z=(x,y) in the Boolean cube write

    zeta=D^(1/2)z,       ||zeta||^2=S,
    t_z=||Pi zeta||^2/S,
    v_z=zeta^T Pi L Pi zeta/S.

Since Pi commutes with L, the exact allocation satisfies

    |v_z|<=t_z,       |u-v_z|<=r0(1-t_z).                (2.1)

In particular t_z belongs to I_u. Set

    g0=D^(1/2)Q D^(-1/2)g.

The removed Gaussian field has the finite uniform bound

    E sup_(z Boolean) |(g-g0)^Tz| <= Lambda_sp,
    Lambda_sp=2n sqrt(w[d+N sqrt(delta(1+delta))]).         (2.2)

Indeed the left side is at most sqrt(S tr(Pi Mhat)). For
q_i=dbar/d_i, the exact identities give

    sum_i q_i=N(1+delta),       sum_i 1/q_i=N,
    sum_i |q_i-1|<=N sqrt(delta(1+delta)).

The last bound is weighted Cauchy--Schwarz using
sum_i (q_i-1)^2/q_i=N delta. Thus
tr(Pi D^(-1))<=dbar^(-1)[d+N sqrt(delta(1+delta))].
Together with M<=2wnI and S=N dbar, this proves (2.2).
No state-mass estimate is used in this Gaussian projection loss.

Use the positive metric and diagonal remainder

    P=D^(1/2)[e Pi+Q(I-eta L)Q]D^(1/2),       E=eD.

They satisfy P>=E>0. The ORIGINAL cell has exact radius

    z^TPz=S[1-eta r0 t_z-eta(u-v_z)],
    z^T(P-E)z=S eta[r0(1-t_z)-u+v_z]
                    <=S eta[r0-u+(1-r0)t_z].            (2.3)

Only the inequality v_z<=t_z was used in the last step; no state with
v_z=t_z has been asserted to exist. This relaxation is useful because
only t_z, not a second real allocation coordinate, must be binned.

Let h=EP^(-1)g0. A pointwise positive-semidefinite Cauchy--Schwarz
split retains the same constrained states:

    g0^Tz <=h^Tz+sqrt(z^T(P-E)z)
              sqrt(g0^TP^(-1)(P-E)P^(-1)g0).             (2.4)

The expectation of the last squared Gaussian factor is T-eR.
Also, by cyclic trace identities,

    sum_i Var(h_i)/d_i=e^2R.                              (2.5)

These facts do not require a constant diagonal after projection,
commutation of Mhat with L, or any replacement of the paired covariance.

## 3. A translated Gaussian inequality and constrained Boolean width

For every sigma>=0 and real b,

                    E|sigma G+b|<=sqrt(kappa sigma^2+b^2),
                    G standard Gaussian.                (3.1)

Here is a proof. Let H(b)=E|G+b|. Rotating two independent standard
Gaussians to independent U,V gives exactly

    H(b)^2-b^2=E[(V^2-(U+sqrt(2)b)^2)_+].

For each fixed V the integrand, as a function of U+sqrt(2)b, is even
and decreasing in absolute value. Gaussian translation decreases its
integral: by layer cake, each superlevel set is a centered interval,
whose Gaussian measure is maximal at zero translation. Consequently
H(b)^2-b^2<=H(0)^2=kappa. Scaling proves (3.1), including sigma=0.

Fix a unit vector a in the spike space, and restrict to Boolean states
with a^Tzeta>=sqrt(tS). For every lambda>=0, Fenchel's elementary
linear-constraint bound and (3.1) give

    E max h^Tz
      <=sum_i E|h_i+lambda sqrt(S d_i)a_i|-lambda S sqrt(t)
      <=sqrt(S)sqrt(kappa sum_i Var(h_i)/d_i+lambda^2 S)
                                               -lambda S sqrt(t).

The second line is Cauchy--Schwarz with weights d_i; correlations
between coordinates do not matter. Minimizing over lambda and using
(2.5) gives the exact fixed-direction estimate

             E max h^Tz<=sqrt(kappa S)e sqrt(R)sqrt(1-t). (3.2)

For t=1 the same result follows by letting lambda tend to infinity.
This is the spike-mass loss INSIDE the Boolean remainder in (2.4),
not merely an unrelated upper to be minimized with the whole width.

## 4. Uniform weighted conditioning and the low-rank net

The direction a in (3.2) may depend on z. We now justify the union
over directions without imposing pointwise conditioning on D.
Put A0=eF, K0=D^(1/2)/sqrt(dbar), J0=K0^(-1). Then

    h=K0 A0 J0 g,       hbar=A0g,       ||A0||op<=1.

Their uniform Boolean-width coupling satisfies

    E||h-hbar||_1<=Lambda_delta,
    Lambda_delta=N sqrt(wn)sqrt(delta)(1+sqrt(1+delta)).   (4.1)

To prove it, split
h-hbar=(K0-I)A0J0g+A0(J0-I)g. The dispersion bounds are

    ||K0-I||_F^2<=N delta,       ||J0-I||_F^2<=N delta.

The constant diagonal M_ii=wn gives
E||J0g||^2=wnN(1+delta) and E||(J0-I)g||^2<=wnN delta.
Apply Cauchy--Schwarz to the first summand's coordinatewise product
and use ||x||_1<=sqrt(N)||x||_2 on the second. This proves (4.1).
Thus hbar is only a uniformly controlled auxiliary field, not a free
change of covariance. Its covariance has operator norm at most 2wn.

For 0<xi<1 take a Euclidean xi-net of the unit sphere in ran(Pi),
of cardinality at most (1+2/xi)^d. If t_z>=t, some net direction
satisfies a^Tzeta>=(1-xi)sqrt(tS). The net concerns the WEIGHTED
state zeta, whose norm is exactly sqrt(S).

For any fixed subset of the Boolean cube, the supremum of hbar^Tz
has Gaussian Lipschitz constant at most sqrt(2wnN)=2n sqrt(w).
Gaussian concentration selects among the net's constrained subsets
at cost at most 2n sqrt(2w d log(1+2/xi)). Transfer to and from hbar
using two copies of (4.1), and apply (3.2) to h itself on each subset.
The resulting finite estimate is

    E max_(z in C, t_z>=t) h^Tz
      <=sqrt(kappa S)e sqrt(R)sqrt(1-(1-xi)^2 t)
          +2Lambda_delta+2n sqrt(2w d log(1+2/xi)).        (4.2)

Empty subsets are ignored. When d=0 only t=0 is nonempty and the
ordinary weighted Boolean bound gives the same assertion directly.
No Gaussian independence between subsets or between projected fields
is assumed. The conditioning cost is O(n^(3/2)sqrt(delta)), while
the net cost is o(n^(3/2)) for fixed xi and d=o(n).

## 5. Finite mass bins, selection costs, and the order of limits

For a nonempty mass bin t_-<=t_z<=t_+, equations (2.2)--(2.4) and
(4.2) give the finite bound

    E max_bin g^Tz
      <=sqrt(S) {
          sqrt(eta[r0-u+(1-r0)t_+](T-eR))
                +sqrt(kappa)e sqrt(R)sqrt(1-(1-xi)^2t_-)}
          +Lambda_sp+2Lambda_delta
                          +2n sqrt(2w d log(1+2/xi)).     (5.1)

The first coefficient is nonnegative in every nonempty bin, by
(2.3). Partition I_u into mass bins of mesh at most h0. Selecting
among their maxima for the ACTUAL g costs at most
2n sqrt(2w log(ceil(1/h0)+1)), since M<=2wnI. This selection does
not require an operator bound on the projected or transformed field.

All normalized trace factors are uniformly bounded at fixed eta,r0:

    S tr(Mhat)=wnN^2(1+delta),       ||F||op<=1/e.

Changing mass endpoints by h0 changes the normalized envelope by
O_(eta,r0)(sqrt(h0)); replacing (1-xi)^2t by t costs O(sqrt(xi)).
First take n to infinity at fixed xi,h0, using delta->0 and d=o(n).
Then let xi,h0 decrease to zero. Equations (2.2), (4.1), and (5.1)
prove precisely (1.3), including sequences with varying u,T,R.
No rate at a moving spectral cutoff or moving resolvent endpoint is used.

For the weighted-shell Gaussian cross-process upper, refine its actual
original/weighted cells by these mass bins before selecting their actual
representatives. The old increment comparison applies unchanged, and
at fixed h0 the number of cells remains polynomial in n. Its original
cell-selection, bin, and padding errors remain as recorded in the
381-line theorem. In particular the deterministic ORIGINAL internal
and cross drift is unchanged. A finite mass mesh can be sent to zero
only after the n-limit, as above.

## 6. Consequence: weak middle law, with arbitrary low-rank outliers

Consider actual zero-original-source cells p=q_A=0 with S=O(n^(3/2)),
delta->0, and representative weighted cross value u_D=c_D/n->4/5.
Assume ONLY that the empirical law of all squared singular values of
W=D_L^(-1/2)B D_R^(-1/2) converges weakly to

                    (delta_0+delta_(4/5))/2.              (6.1)

There is no assumption that the actual operator radius tends to
2/sqrt(5). It may remain as large as one because of spectral outliers.
The centered sign parameters here are w=1 and k=kappa.

The frozen 279-line actual-radius theorem supplies the uniform
o(n^(3/2)) comparison to the genuine positive pure-cross field

    M0=n(I-kappa u_D L),       diag M0=n,

and |u-u_D|<=2sqrt(delta). Its constant-diagonal congruence estimate
also gives, for J0=sqrt(dbar)D^(-1/2),

    ||J0 M0 J0-M0||_1
              <=nN sqrt(delta)(1+sqrt(1+delta)).           (6.2)

It applies to the truncated F in (1.2): both F^2 and F-eF^2 have
bounded operator norms at fixed e>0. Thus the exact weighted traces
in (1.3) have the same normalized limits as their scalar-diagonal
reference traces for the SAME actual M0,L. This step does not replace
the actual original cell radius by a weighted shell identity.

Put rho=2/sqrt(5). Fix r0 in (rho,1), take beta=9/10,
eta=beta/r0, and e=1/10. By (6.1), rank(1_(|L|>r0))=o(n), so
(1.3) applies. The normalized reference traces are those of the
279-line theorem with L/r0, including all retained bulk eigenvalues.
They use the original dimension as denominator, not a renormalization
of the retained spectrum. The omitted spectral mass tends to zero.
First let n tend to infinity with this cutoff fixed. Then let
r0 decrease to rho. The inverse gap e=1/10 stays fixed throughout;
the limiting trace functions are continuous and uniformly bounded.

For clarity, their exact endpoint formulas are as follows. Set

    s=kappa rho^3,
    IA=(461+100s)/722,       IB=(18461-18000s)/722,
    F0=(9/10)sqrt((1-rho)IA),
    G0=(1/10)sqrt(kappa IB).

After normalization by 2n^(3/2), the general allocation envelope
therefore tends, uniformly for 0<=t<=1, to

                 U(t)=F0 sqrt(1+t/rho)+G0 sqrt(1-t).      (6.3)

The term 1+t/rho comes from the retained original radius:
1-u/rho+(1/rho-1)t=(1-rho)(1+t/rho). No outlier state mass was
discarded, and the relaxation v_z<=t_z was only an upper inequality.

The frozen 279-line rational calculation proves

    rho>25/28,       s<1/2,       F0<1/4,
    G0<31/100,       F0+G0<14/25.                         (6.4)

Only one new elementary derivative check is needed. Since
IB>9461/722>13 and kappa>7/11,

    F0/rho<7/25,
    G0^2>91/1100>49/625,       hence G0>7/25.              (6.5)

The final squared comparison is 56875>53900. Consequently
U'(0)=(F0/rho-G0)/2<0. The function U is concave on [0,1], and
strictly concave in its interior, so U(t)<=U(0)<14/25 for EVERY
allocation. This proves, including the prior cross-process errors,

    limsup E max_cell X_z/(2n^(3/2))
                       <=14/25<2sqrt(2)/5.              (6.6)

The last inequality follows from 196/625<200/625. Thus the middle-law
upper is below sqrt(2)alpha at alpha=2/5 using the WEAK empirical law
alone, with no actual-radius convergence premise. This is a conditional
actual-cell theorem; neither existence of the diagnostic source nor a
bound over all other internal-energy cells is asserted.

## 7. Prerequisites, collaboration, and scope

The 381-line weighted-shell source was fully read:
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
The complete 279-line prerequisite was fully read and remains unchanged:
`original_mo_actual_cross_radius_shell_upper.md`, SHA256
`44fa3e7361e2142b20dce58d2dde727458db786529690f15e752390b8081725f`.
Its cited compatibility and Boolean-shell prerequisites retain their
disclosed earlier provenance. The 312-line pure-cross source was also
fully read; no old certificate or mathematical computation was rerun.
A targeted existing-object search found no prior spike-mass remainder.

Root supplied the translated-Gaussian inequality and the decisive idea
of retaining spike mass inside the Boolean remainder. The exact worker
derived the finite allocation metric, Gaussian projection loss, weighted
Fenchel version, conditioning coupling, and the full source. The proof
worker independently checked these links and supplied an independent
middle-law derivative enclosure. The exact worker's simpler derivative
check is written in (6.5). These are contributing roles, not independent
whole-new-source reviews. A frozen-source reviewer must disclose them.

Only this /tmp source and its author receipt are written by its author.
No canonical file, computation, publication, or backup was changed or
performed. The actual-r theorem remains a valid separate result. This
extension removes its operator-radius convergence premise in the stated
weak-middle-law corollary, not the near-scalar or zero-original-source
premises. The original MO convergence question remains OPEN.
