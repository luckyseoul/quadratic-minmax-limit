# Whole-profile cross-bulk upper at the original source target

2026-09-06. Analytic actual-cell theorem. The full spectral profile is
arbitrary; only a bulk edge and its first moment enter the upper.
Low-rank operator outliers are allowed and their Boolean state mass
is retained through the spike-remainder theorem. No mathematical
program, parameter scan, solver, census, or numerical optimization ran.

## 1. Actual hypotheses and the whole-profile bound

Consider sequences of actual paired complete signings

    K=[[A,B],[B^T,-A]],       N=2n,

with feasible positive diagonal D satisfying D+-K>=0. Throughout assume

    S=tr D=O(n^(3/2)),       dbar=S/N,
    delta=S tr(D^(-1))/N^2-1 ->0.                         (1)

The fixed trace cap is part of the theorem; it is retained for the
existing Gaussian cross-process padding and transfer errors. Set

    W=D_L^(-1/2) B D_R^(-1/2),       ||W||op<=1,
    m_n=tr(W^TW)/n.

Choose actual final original/weighted cells whose original internal
energies satisfy p=q_A=0. Their common original cross energy is c,
and their representatives are chosen WITHIN those final cells. Write

    u_n=2c/S,       u_D,n=c_D/n,       c_D=x^TWy.

The actual dispersion comparison gives |u_n-u_D,n|<=2sqrt(delta).
Suppose u_n->u in [0,1] and m_n->m. Fix a bulk edge 0<r<1 and assume

    for every fixed R in (r,1),
    #{j: singular_value_j(W)>R}/n ->0.                   (2)

This is an asymptotic bulk condition, NOT ||W||op<=r+o(1).
A weak limiting squared-singular law supported in [0,r^2] implies
(2), but no full empirical-law limit is required for this theorem.

Let X_z be the centered-sign base Gaussian cross process of the
weighted-shell theorem, with parameters w=1,k=kappa=2/pi. All of
its original cell, padding, and comparison errors remain present.
Then

    limsup E max_cell X_z/(2n^(3/2)) <= Gamma_r(m,u),      (3)

where, with d_r=1-(1-kappa)r,

    Gamma_r(m,u)^2
      = kappa(1-u)/d_r *
        [1+kappa^2 m(1-u d_r)
                /((1-r)(1-r+2kappa r))].                (4)

No optimality of A, B, or D is used to prove (3). Its applications
below compare the result with the ORIGINAL Phi(A), not a substituted
source normalization or a prescribed internal spectral law.

## 2. Eliminate every spike allocation by one quadratic inequality

Fix first a cutoff R in (r,1), set

    d_R=1-(1-kappa)R,       beta=kappa R/d_R,
    eta=beta/R,            e=1-beta=(1-R)/d_R>0.

Condition (2) makes the associated cross-spike projection rank o(n).
The frozen 347-line spike-remainder theorem applies to the actual
field, with all state allocations retained. The pure-cross comparison
and constant-diagonal congruence from the frozen 279-line theorem
justify its scalar numerical trace upper, with o(n^(3/2)) errors at
this FIXED cutoff and inverse gap e. The covariance still uses the
actual W and u_D,n, and the original radius still uses u_n.

For the reference trace calculation, let nu_Q be the actual retained
squared-singular measure, with weight 1/n per retained singular value.
It is a subprobability measure, NOT renormalized after removing spikes.
For y in [0,1] define

    A_s,beta(y)=[1+(beta^2-2beta+s(1-2beta))y
                            +s beta^2 y^2]/(1-beta^2 y)^2,
    B_s,beta(y)=[1+(beta^2-2s beta)y]/(1-beta^2 y)^2,
    s=kappa u R.

Put I_A=integral A_s,beta(y/R^2)dnu_Q(y) and define I_B similarly.
Write A0=beta sqrt(I_A), B0=(1-beta)sqrt(kappa I_B).
After normalization by 2n^(3/2), the spike theorem's envelope is

    A0 sqrt(1-u/R+(1/R-1)t)+B0 sqrt(1-t),                (5)

over its feasible interval max(0,(u-R)/(1-R))<=t<=1. Terms in (5)
are nonnegative there. Cauchy--Schwarz with weights R and 1-R gives

    (5)^2 <= (1-u)[A0^2/R+B0^2/(1-R)].                  (6)

Indeed R[1-u/R+(1/R-1)t]+(1-R)(1-t)=1-u exactly.
Thus (6) holds simultaneously for EVERY allocation, including u>R.
It does not assert that the relaxed assignment of spike cross energy
equal to t is attainable. No search over t is needed.

## 3. A matched metric turns the two traces into one convex resolvent

The chosen beta has the exact properties

    beta^2/R=kappa^2 R/d_R^2,
    kappa(1-beta)^2/(1-R)=kappa(1-R)/d_R^2.

Before pairing the two spectral signs, A_s,beta is represented by
(1-sv)(1-v)/(1-beta v)^2 and B_s,beta by
(1-sv)/(1-beta v)^2. The identity

    kappa R(1-v)+(1-R)=d_R(1-beta v)

therefore collapses the right side of (6) to

    kappa(1-u) integral f_R,u(y)dnu_Q(y),
    f_R,u(y)=(d_R-kappa^2 u y)/(d_R^2-kappa^2 y).        (7)

Pairing v and -v is legitimate because the cross operator has an
even spectrum. This is a trace identity for the numerical reference,
not a replacement of the actual field or cross law. At finite n its
two u parameters are u_n and u_D,n; their difference tends to zero,
uniformly in these bounded reference integrands at fixed R. Hence
the displayed single-u limiting expression is justified.

The function f_R,u is positive, increasing, and convex on [0,R^2]:

    f'_R,u(y)=kappa^2 d_R(1-u d_R)/(d_R^2-kappa^2 y)^2,
    f''_R,u(y)=2kappa^4 d_R(1-u d_R)
                                      /(d_R^2-kappa^2 y)^3.

The denominator is positive because d_R-kappa R=1-R>0.
Its endpoint chord, the bulk mass at most one, and the bulk first
moment at most m_n consequently give

    integral f_R,u dnu_Q
      <=1/d_R+kappa^2 m_n(1-u d_R)
                    /[d_R(d_R^2-kappa^2 R^2)].           (8)

Here d_R^2-kappa^2R^2=(1-R)(1-R+2kappa R). The omitted spectral
mass is o(1); no hidden normalization changes the first moment.
All finite trace comparisons and mass/net binning errors are uniform
over the actual cells at fixed R. First let n tend to infinity,
then let R decrease to r. Its inverse gap remains bounded away from
zero near this fixed r<1. Equations (6)-(8) prove (3)-(4).

In particular, (4) needs only the asymptotic first moment and the
bulk-edge condition (2). Neither individual atom masses nor existence
of a limiting spectral law has been assumed.

## 4. Complete-source moment cap and a strict ORIGINAL-target region

Actual completeness and near-scalar dispersion force limsup m_n<=1/2.
Here is the short finite proof. Diagonal entries of
(D^(-1/2)KD^(-1/2))^2 are at most one, giving

    (N-1)tr(D^(-1))<=S,       dbar^2>=N-1.

Literal squared cross entries and arithmetic--geometric mean give

    m_n=(sum_L d_i^(-1))(sum_R d_i^(-1))/n
        <=n(1+delta)^2/dbar^2
        <=n(1+delta)^2/(2n-1).                           (9)

Now assume (2) with r=9/10 and liminf u_n>=7/8. No extra first-moment
hypothesis is needed: use (9), subsequences, and m<=1/2 in (4).
The expression there decreases with u and increases with m. The
elementary bounds 7/11<kappa<16/25 imply, at r=9/10,

    d_r>37/55,
    (1-r)(1-r+2kappa r)>137/1100.

Therefore its prefactor and correction satisfy

    kappa(1-u)/d_r < (16/25)(1/8)/(37/55)=22/185,
    1-u d_r < 1-(7/8)(37/55)=181/440,
    kappa^2 m(1-u d_r)/((1-r)(1-r+2kappa r))
        < (256/625)(1/2)(181/440)/(137/1100)
        =11584/17125.

Consequently the ENTIRE profile class has

    Gamma_(9/10)(m,u)^2
      <(22/185)(1+11584/17125)
       =631598/3168125<1/5,                              (10)

where the final comparison is 3157990<3168125.
The verified original-source lower from CORE.md is

    Phi(A)>=n sqrt(n-1)/pi.

The actual paired target, normalized in precisely the convention of
(3), is (2sqrt(2)Phi(A))/(2n^(3/2))=sqrt(2)Phi(A)/n^(3/2).
Its limiting square is at least 2/pi^2. Since pi<22/7,

                         2/pi^2>49/242>1/5.             (11)

Thus (3), (10), and (11) give a strict upper below the ORIGINAL
2sqrt(2)Phi(A) scale, uniformly for this pure-cross class. For example,

    limsup [E max_cell X_z-2sqrt(2)Phi(A)]/n^(3/2)
                   <=2/sqrt(5)-2sqrt(2)/pi<0.            (12)

No assumption that Phi(A)/n^(3/2) tends to 2/5, or that A is an
original minimizer, is needed for (12). Fixed trace cap (1), zero
ORIGINAL internal energies, dispersion, and bulk condition (2) remain.

## 5. A second wider diagnostic-ratio region

At the same bulk edge r=9/10, suppose instead liminf u_n>=4/5 and
limsup m_n<=2/5. Using d_r>2/3 and the denominator in (4) greater
than 3/25 gives

    kappa(1-u)/d_r<24/125,
    1-u d_r<7/15,
    kappa^2 m(1-u d_r)/((1-r)(1-r+2kappa r))<3584/5625.

Hence

    Gamma_(9/10)(m,u)^2<221016/703125<8/25,               (13)

with squared margin 3984/703125. This entire profile class is below
sqrt(2)alpha at alpha=2/5. Unlike (12), comparison (13) alone uses
that separate source-scale target; it is not asserted to imply the
original target at every smaller alpha. This distinguishes a useful
diagnostic region from the genuine uniform-source region in Section 4.

## 6. Dependencies, provenance, and remaining scope

The complete frozen 347-line prerequisite and its 63-line author
receipt were read before this proof:

    original_mo_low_rank_cross_spike_mass_upper.md
    SHA256 30347140ecf9fb2458444fb152490c601fe81d8a1733e90f31be692126ecdf1c.

The author also read its complete 279-line and 381-line prerequisites,
and the complete 237-line CORE.md (including the proof of the actual
Gaussian source lower used in (11)-(12)), SHA256
f70cb5120d6554c75bc9b10311d993087099796cce79eee4effe5342c6dfde26.
The 279-line and 381-line hashes are recorded in the 347-line source.

Root supplied the spike-mass remainder idea in that prerequisite;
the exact worker authored it, and this author contributed checks.
Here the proof author derived the allocation elimination, matched
metric/resolvent identity, all-profile first-moment bound, and the
two rational regions. The exact worker independently checked the new
resolvent/chord algebra and Section 5 rational constants before this
source was written; that preliminary check is not a complete frozen
source review. No duplicate arithmetic program or certificate ran.

The condition concerns an actual BULK edge, permitting operator
outliers of vanishing relative rank. It does not derive that condition,
the high-u condition, or delta->0 from conditional optimality. Other
pure-cross ratios, higher bulk edges, and general internal-energy cells
remain to be controlled in the full weighted-shell maximum. Therefore
the conditional paired-norm theorem and original all-orders convergence
remain OPEN. Only this separate /tmp proof was written here; no
canonical source, reviewed prerequisite, or backup was changed.
