# Full-SDP gap and two actual original quadratic phases

2026-09-05. Analytic original-norm theorem. No mathematical computation,
numerical phase evaluation, signing census, or solver was run.

This proves the small-gap source-scale consequence
`Phi(K)>=(kappa/2-o(1))tau(K)` under a fixed original norm cap.
The proof uses two actual Gaussian sign phases for the SAME complete
source K, with a shared coordinate-dependent diagonal normalization.
It does not obtain this constant by transferring a rectangular bound
through `beta(K)<=4Phi(K)`. The final section proves that the displayed
gap-only estimate itself is insufficient for a specified large-gap range.

## 1. Actual source, canonical gap, and a finite evaluated inequality

Let K be a complete symmetric zero-diagonal signing of order N>=3.
For the conditional application it is the literal paired matrix
`K=[[A,B],[B^T,-A]]`, with N=2n. Throughout,

    Phi(K)=max_(z in {+-1}^N) |z^T K z|/2,
    q=N-1,           L=||K||op,          S_j=tr|K|^j,
    kappa=2/pi,      rho=1-kappa.

Choose an actual trace-optimal same-diagonal SDP majorizer

    D-K>=0,       D+K>=0,       S=tr D=tau(K).

Set

    g=S-S_3/q>=0,       gamma=g/S in [0,1),
    eta=S/[N sqrt(q)]>=1,
    delta=S tr(D^(-1))/N^2-1,
    c_*=1-1/sqrt(2),
    b_*=min(1/2, eta(1+N/q+c_*)sqrt(gamma)).             (1.1)

Then the following exact finite-dimensional original-norm inequality
holds:

    Phi(K)/S
      >=(kappa/2)[1-gamma-b_*-sqrt(2gamma b_*)]
                       -rho(L/q)(1-gamma).              (1.2)

No small-gap hypothesis is required for (1.2). Its right-hand side
is permitted to be negative. Under a fixed original norm cap
`Phi(K)<=C N^(3/2)`, its useful small-gap consequence is

    Phi(K) >= (kappa/2)S
                  -O_C(N^(3/2)sqrt(gamma)+N^(5/4)).      (1.3)

In particular gamma tending to zero gives the claimed leading
coefficient kappa/2. Small gamma is an extra actual condition; the
theorem does not infer it from conditional or original minimality.

The canonical residual and dispersion inputs are proved in
`original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
The earlier original-phase spectral-moment theorem supplies the
Gaussian sign identity and elementary arcsine remainder mechanism;
the coordinatewise normalization below is a new construction that
avoids the earlier maximum-diagonal phase denominator.

## 2. Two actual phases with a shared coordinate normalization

Write

    H=K|K|,       h_i=H_ii,       v_i=q+|h_i|,
    P=diag(v_i^(-1/2)).

Since `K^2+H=2K_+^2` and `K^2-H=2K_-^2` are PSD and have
diagonals q+h_i and q-h_i, respectively, one has |h_i|<=q.
Both matrices

    C_+=P[K^2+H+diag(|h_i|-h_i)]P,
    C_-=P[K^2-H+diag(|h_i|+h_i)]P                       (2.1)

are positive semidefinite and have diagonal one. Thus they are actual
Gaussian correlation matrices, including when singular. Their
coordinate normalization is the SAME v_i in the two phases, but it
need not be constant in i. Let z_+ and z_- be their Gaussian signs.
Each realization is an actual signing state of the original K.

For any correlation matrix C,

    E Q_K(sign G)=kappa sum_(i<j) K_ij arcsin(C_ij),
    Q_K(z)=z^T K z/2.

The elementary remainder bound

    |arcsin(t)-t| <= (pi/2-1)t^2,       -1<=t<=1

shows that the absolute error beyond the linear trace term is at
most `(rho/2)||C_off||_F^2`. Since v_i>=q, the two off-diagonal
Frobenius norms in (2.1) satisfy

    ||(C_+)_off||_F^2 <=4 tr(K_+^4)/q^2,
    ||(C_-)_off||_F^2 <=4 tr(K_-^4)/q^2.

Padding diagonals do not enter either estimate or the trace pairing
with K. Therefore, with

    J=<K,PHP>=sum_(i!=j) K_ij H_ij/sqrt(v_i v_j),

the difference of the two expected ORIGINAL quadratic energies obeys

    E Q_K(z_+)-E Q_K(z_-)
                       >=kappa J-2rho S_4/q^2.

Each of the two expectations lies in [-Phi(K),Phi(K)], so

                    Phi(K)>=(kappa/2)J-rho S_4/q^2.     (2.2)

This is the normalization responsible for the factor kappa/2.
The two phases are not cross states of a rectangular objective,
and no source signing is replaced.

## 3. The actual gap controls the diagonal-normalization loss

The canonical residual

    R_1=DK-H

satisfies

    sum_(i,j) |(R_1)_ij|^2/d_i <=4qg.                   (3.1)

This is one nonnegative summand of the complete weighted residual
bound in the cited compatibility theorem. Because K has zero
diagonal, `(R_1)_ii=-h_i`. Consequently

    sum_i h_i^2/d_i<=4qg,
    sum_i |h_i|<=2sqrt(qSg).                            (3.2)

The second inequality is Cauchy--Schwarz with the actual positive
weights d_i. It does not assert a uniform coordinate bound better
than |h_i|<=q.

Put

    f_i=sqrt(q/v_i)=(1+|h_i|/q)^(-1/2),
    a_i=1-f_i,       m_ij=1-f_i f_j,
    A_0=sum_i a_i,
    B_0=(1/q)sum_(i!=j) d_i m_ij.                       (3.3)

Then

    1/sqrt(2)<=f_i<=1,
    0<=a_i<=c_*,        0<=m_ij<=1/2,
    A_0<=sum_i |h_i|/(2q)<=sqrt(Sg/q),
    0<=B_0<=S/2.                                       (3.4)

The bound for a_i follows from
`1-(1+x)^(-1/2)<=x/2` for x>=0. The last bound uses the q
off-diagonal entries in each row, and only those entries.

To control the weighted a_i sum, put dbar=S/N. The dispersion
identity gives

    sum_i |d_i-dbar|<=S sqrt(delta).

Indeed for t_i=d_i/dbar use Cauchy--Schwarz on
`sum |t_i-1|` with factors sqrt(t_i) and |t_i-1|/sqrt(t_i).
The total positive deviation of d_i from dbar is half this absolute
sum. Thus

    sum_i d_i a_i
       <=dbar A_0+(c_*/2)S sqrt(delta).

Using `m_ij<=a_i+a_j` and dropping a nonpositive diagonal correction
now gives

    B_0<=sum_i d_i a_i+(S/q)A_0
       <=(S/N+S/q)sqrt(Sg/q)+(c_*/2)S sqrt(delta).        (3.5)

The actual canonical-gap dispersion theorem says
`delta<=4eta^2 gamma`. Dividing (3.5) by S and retaining (3.4)
therefore proves

                              B_0/S<=b_*.              (3.6)

All diagonal-distribution costs have been retained in this estimate;
no d_max or maximum phase variance is substituted.

## 4. Masked residual and the original source-scale bound

For i!=j the literal sign identity K_ij^2=1 gives

    K_ij H_ij/q=d_i/q-K_ij(R_1)_ij/q.

Since the unnormalized canonical objective is
`<K,H/q>=S-g`, equations (3.3) yield the exact relation

    J=S-g-B_0+(1/q)sum_(i!=j)K_ij(R_1)_ij m_ij.          (4.1)

The diagonal is omitted throughout this identity; the sign-square
identity is not falsely applied to K_ii=0. From (3.1), weighted
Cauchy--Schwarz, and `m_ij^2<=m_ij/2`,

    |(1/q)sum_(i!=j)K_ij(R_1)_ij m_ij|
      <=(1/q)sqrt(4qg)sqrt(sum_(i!=j)d_i m_ij^2)
      <=sqrt(2gB_0).                                    (4.2)

Thus

    J>=S-g-B_0-sqrt(2gB_0)
      >=S[1-gamma-b_*-sqrt(2gamma b_*)].                 (4.3)

The last step uses that `b+sqrt(2gamma b)` is increasing for b>=0.
Finally S_4<=L S_3=LqS(1-gamma). Substitution into (2.2) proves
the finite bound (1.2).

For the norm-only asymptotic interpretation use the reviewed facts

    L^2<=8Phi(K),
    S=tau(K)<=4G Phi(K),     G=pi/[2log(1+sqrt(2))].

Under Phi(K)<=C N^(3/2), eta is bounded and L/q=O_C(N^(-1/4)).
The quantity b_* is at most a bounded multiple of sqrt(gamma), and
`sqrt(2gamma b_*)=O_C(gamma^(3/4))`. For 0<=gamma<1, both gamma
and gamma^(3/4) are at most sqrt(gamma). Equation (1.3) follows.
The bilinear SDP estimate in this paragraph only bounds the source
scale S under a norm cap; it does not supply the factor kappa/2 in
the actual original-phase lower bound.

## 5. What this evaluates for the paired original problem

For a sequence of actual paired sources with a fixed norm cap and
gamma tending to zero, (1.3) gives

    Phi(K)>=(kappa/2-o_C(1))S.                           (5.1)

This is a statement about the original norm of the literal K, not
about its rectangular bilinear norm or an artificial moment model.
For an active positive original cross value c it implies the correctly
normalized source-scale constraint

    S<=(2/kappa+o_C(1))Phi(K),
    2c/S >= (kappa-o_C(1)) c/Phi(K).                    (5.2)

The latter is interpreted along a bounded source-scale sequence and
an active cross value; it does not assert c=Phi(K) on an arbitrary
shell. At zero original internal values, active c=Phi(K) is a
separate, explicit saturation condition. In that case (5.2) and
the uniform compatibility theorem give

                      2c_D/N >= kappa-o_C(1).           (5.3)

For nonzero internal values the ratio c/Phi(K) must be retained.
Equation (5.3) does not evaluate the Gaussian field width by itself;
it supplies the actual source-scale input for that separate analysis.

## 6. A precise limitation on the complementary large-gap branch

The explicit right-hand side of (1.2), at fixed eta and L/q, is not
a positive gain with gamma. Its main bracket decreases as gamma
increases. In fact it is already vacuous for every gamma>=1/4:
eta>=1 and N/q>1 imply

    eta(1+N/q+c_*)sqrt(gamma)>1/2,
    b_*=1/2,
    1-gamma-b_*-sqrt(2gamma b_*)
                  =1/2-gamma-sqrt(gamma)<=-1/4.        (6.1)

The subtracted finite fourth-moment term cannot improve this bound.
Thus this PARTICULAR fully evaluated uniform phase estimate cannot
exclude, or settle the original upper in, that large-gap range.
This is a proved insufficiency of the displayed estimate, not an
impossibility theorem for all optimized two-phase constructions.

The useful new conclusion is the actual small-gap source normalization
with coefficient kappa/2 and controlled finite error. A strategy that
splits on the canonical gap still needs a separate evaluated argument
for the complementary positive-gap sources. Nothing here proves that
every conditional optimizer has a small gap or establishes original
MO convergence.
