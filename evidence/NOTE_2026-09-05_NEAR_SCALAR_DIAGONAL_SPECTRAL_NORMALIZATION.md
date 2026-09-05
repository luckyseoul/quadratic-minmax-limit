# Near-scalar diagonals: actual full-spectrum normalization without small gap

2026-09-05. Analytic original-norm theorem. No mathematical computation,
signing construction, search, or numerical spectral evaluation was run.

Small arithmetic-harmonic diagonal dispersion is sufficient for the
normalization below. The full canonical-primal gap need not be small.
The full weighted spectrum always belongs to the original actual source.
An auxiliary principal submatrix is used only to LOWER-bound its original
quadratic norm; it never replaces the source covariance or cross block.

## 1. Actual source, diagonal, and spectral quantities

Let K be an actual complete symmetric zero-diagonal signing of order
N>=4. Let D be a positive diagonal satisfying D-K>=0 and D+K>=0.
Write

    S=tr D,       dbar=S/N,       q=N-1,
    delta=S tr(D^(-1))/N^2-1,
    T=D^(-1/2) K D^(-1/2),       ||T||op<=1,
    r=q/dbar^2,                  kappa=2/pi,   rho=1-kappa.

Let mu be the empirical probability measure of the N eigenvalues of
the ACTUAL T, and put

    mu_j=integral |lambda|^j dmu(lambda),       j=1,2,3.

In particular mu has support in [-1,1]. Feasibility and the canonical
complete-signing lower bound imply S>=N sqrt(q), hence 0<r<=1.
Optimality of D is not needed below if its source-scale cap is supplied
separately. Define the actual original-norm ratio

    U=2Phi(K)/S,       Phi(K)=max_(z in {+-1}^N)|z^TKz|/2.

The asymptotic conclusion, uniform under S<=C N^(3/2), is

    U >= max{ kappa(1+r)mu_3/(2r), kappa r/mu_1 }
                                  -O_C(delta^(1/3)+N^(-1/2)),    (1.1)

as delta tends to zero and N tends to infinity. Also

                         mu_2=r+O(delta^(1/3)).                   (1.2)

The first term of (1.1) is an actual original two-phase cubic-moment
bound. The second is the actual nuclear-moment bound. Both remain
available with positive canonical gap; neither was obtained by discarding
that gap inside the earlier small-gap masked-residual estimate.

The full original-phase prerequisite is
`original_mo_original_phase_spectral_moment.md`, SHA256
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
Its row-moment and nuclear facts are restated and justified below.

## 2. Good coordinates and exact finite comparison parameters

Put t_i=d_i/dbar. Then the exact dispersion identity is

             (1/N)sum_i (t_i-1)^2/t_i=delta.

Choose 0<epsilon<=1/2 such that

                   delta(1+epsilon)/epsilon^2<=1/4.

Let I contain exactly those coordinates with |t_i-1|<=epsilon, and set

    M=|I|,       b=(N-M)/N,       a=M/N=1-b,
    q_I=M-1,     r_I=q_I/dbar^2,
    L_*=(1+epsilon)dbar,
    C_j=3j(3/2)^(j-1),       j=1,2,3.

Since (t-1)^2/t is at least epsilon^2/(1+epsilon) outside the good
interval, one has

    0<=b<=delta(1+epsilon)/epsilon^2<=1/4,
    M>=3,                    (2/3)r<=r_I<=r.                      (2.1)

The last lower bound uses N>=4 and M>=3N/4. Denote by K_I and T_I
the corresponding ACTUAL principal submatrices. Their roles are
different: K_I is a complete signing for a lower-bound argument, whereas
T_I is used to compare its spectral moments with the full actual T.

The inequalities D_I+-K_I>=0 and D_I<=(1+epsilon)dbar I_M give

                              ||K_I||op<=L_*.                      (2.2)

For every fixed signing on I, extend it by independent unbiased signs
on the remaining coordinates. The conditional expected full quadratic
energy equals its K_I energy: all crossing and removed internal terms
have zero mean. Thus, exactly,

                              Phi(K)>=Phi(K_I).                    (2.3)

No assertion that K_I is itself an optimizer is required.

## 3. Explicit rank and congruence errors for the full actual spectrum

For j=1,2,3 define

                    E_j=(1/N)tr|K_I/dbar|^j.

Principal-eigenvalue interlacing gives

    0<=mu_j-(1/N)tr|T_I|^j<=2b.                                  (3.1)

Indeed deleting N-M coordinates loses at most N-M from the sum of
the positive j-th powers and at most N-M from the sum of the negative
absolute j-th powers, because the full spectrum is in [-1,1]. This
proves (3.1) without an operator norm bound on K/dbar.

On I let Q=diag(sqrt(t_i)). Then, exactly,

             K_I/dbar=Q T_I Q,
             ||Q-I||op<=epsilon,       ||Q||op<=sqrt(1+epsilon).

Consequently

    ||Q T_I Q-T_I||op
       <=epsilon(1+sqrt(1+epsilon))<=3epsilon.                    (3.2)

Both symmetric matrices in (3.2) have spectrum in
[-(1+epsilon),1+epsilon]. Weyl's ordered-eigenvalue comparison and the
scalar Lipschitz bound for |x|^j therefore imply

    |E_j-(1/N)tr|T_I|^j|<=C_j epsilon,
    E_j>=mu_j-2b-C_j epsilon,
    E_j<=mu_j+C_j epsilon.                                      (3.3)

In particular C_1=3, C_2=9, and C_3=81/4. Literal sign squares in
the ACTUAL K_I give the additional identity

                           E_2=a r_I.                            (3.4)

The difference r-a r_I equals b(N+M-1)/dbar^2 and is at most
(5/2)b for N>=4, using dbar^2>=N-1. Combining this with (3.3)
gives the convenient finite bound

                           |mu_2-r|<=5b+9epsilon.                (3.5)

All moments retained on the right sides belong to the FULL original
weighted matrix T. The removed coordinates are not discarded from mu.

## 4. Actual Gaussian phases on the auxiliary complete signing

For any complete signing J of order M with ||J||op<=L_*, let
q_I=M-1 and h_i=(J|J|)_ii. Its local spectral measure has mean zero,
second moment q_I, and support in [-L_*,L_*]. These three scalar facts
give

                    |h_i|<=q_I(L_*^2-q_I)/(L_*^2+q_I).            (4.1)

For clarity, the one-sided argument bounds P_2=E X_+^2 and
N_2=E X_-^2, with P_2+N_2=q_I. If t=E X_+=E X_-, then
P_2<=L_*t, while the two nonzero-side probabilities are at least
t/L_* and t^2/N_2. Their sum is at most one, giving
q_I P_2<=L_*^2 N_2. This proves the upper bound for h_i=P_2-N_2;
apply it to -X for the other side. Necessarily L_*^2>=q_I.

Set the COMMON constant phase variance

              v_*=2q_I L_*^2/(L_*^2+q_I)>0.

By (4.1), v_*>=q_I+|h_i|. Therefore

    C_+=[J^2+J|J|+diag(v_*-q_I-h_i)]/v_*,
    C_-=[J^2-J|J|+diag(v_*-q_I+h_i)]/v_*

are genuine PSD unit-diagonal Gaussian correlations. The common
constant variance makes their linear original-energy difference
kappa tr|J|^3/v_*. The two arcsine remainder magnitudes sum to at
most 2rho tr|J|^4/v_*^2. Since both expected original energies lie
in [-Phi(J),Phi(J)],

              2Phi(J)>=kappa tr|J|^3/v_*
                                      -2rho tr|J|^4/v_*^2.       (4.2)

This uses actual original quadratic states, not a rectangular-norm
transfer. The diagonal padding contributes to neither trace pairing
with J nor the off-diagonal remainder bounds.

Apply (4.2) to J=K_I, use (2.3), (3.3), and
tr|K_I|^4<=M L_*^4. With [x]_+=max(x,0), the exact finite bound is

    U >= kappa[(1+epsilon)^2+r_I]/[2r_I(1+epsilon)^2]
                            *[mu_3-2b-(81/4)epsilon]_+
          -rho a[(1+epsilon)^2+r_I]^2/[2dbar r_I^2].              (4.3)

The fourth-moment term in (4.3) has been explicitly retained. It is
O_C(N^(-1/2)) under the stated source-scale cap and (2.1).

## 5. Actual nuclear phases and their exact finite consequence

The original nuclear bound for the same complete J is

                  2Phi(J)>=kappa M^2 q_I/tr|J|.                  (5.1)

Here is the common-diagonal mechanism behind (5.1). The PSD matrices
|J|+-J have the same positive diagonal ell_i=|J|_ii. Normalize both
by diag(ell_i^(-1/2)). On each original edge their correlation
difference, multiplied by J_ij, is 2/sqrt(ell_i ell_j).
The slope of arcsine is at least one, so the difference of their
expected original energies is at least
2kappa sum_(i<j)(ell_i ell_j)^(-1/2). Both expectations lie in
[-Phi(J),Phi(J)]. Finally

    sum_(i<j)(ell_i ell_j)^(-1/2)
       >=sum_(i<j)2/(ell_i+ell_j)
       >=M^2(M-1)/(2sum_i ell_i),

by scalar arithmetic-geometric mean and Cauchy--Schwarz. This gives
(5.1). Positivity of each ell_i follows because a zero diagonal of
the PSD |J| would give a zero row and hence (J^2)_ii=0, contrary to
q_I>0. Thus no singular normalization is hidden here.

Applying (5.1) to K_I, followed by (2.3) and the UPPER E_1 estimate
in (3.3), gives the second exact finite bound

                         U>=kappa a^2 r_I/(mu_1+3epsilon).        (5.2)

Its denominator is positive. Equations (4.3) and (5.2) hold
simultaneously, so one may take their maximum.

## 6. Uniform near-scalar limit and paired active-cell use

For 0<delta<=1/512 choose epsilon=delta^(1/3). Then
epsilon<=1/8 and b<=(9/8)delta^(1/3)<1/4. Under S<=C N^(3/2),
r is bounded below by a positive constant depending only on C;
the same is true for r_I. Formula (3.5) proves (1.2). Since
mu_1>=mu_2, its denominator is uniformly separated from zero for
sufficiently small delta. The coefficients in (4.3) and (5.2) can
therefore be compared uniformly with those at a=1, r_I=r, epsilon=0.
Using mu_j<=1, their combined error is
O_C(delta^(1/3)+N^(-1/2)), proving (1.1). If delta=0, D=dbar I,
and the same assertion follows directly with no trimming or dispersion
error. Alternatively one may first take delta to zero at fixed epsilon
and then take epsilon to zero; no uncontrolled endpoint is used.

In the EXACT scalar case, the full canonical gap
gamma=(S-tr|K|^3/q)/S satisfies mu_3=r(1-gamma). Hence the first
term of (1.1) becomes

                        U>=kappa(1+r)(1-gamma)/2-O_C(N^(-1/2)).

This remains positive in the positive-gap range where the earlier
masked-residual estimate was vacuous. The nuclear term still retains
the actual mu_1. The displayed identity between gamma and mu_3 is
asserted here only for exactly scalar D; it is not silently transferred
to a near-scalar diagonal with potentially exceptional coordinates.

For the literal paired source K=[[A,B],[B^T,-A]] with N=2n, retain
the original actual weighted cross matrix

                 W_D=D_L^(-1/2) B D_R^(-1/2).

On an actual final original/weighted refined cell with original
p=q_A=0 and active positive c=Phi(K), let c_D be its actual within-cell
representative's weighted cross value and u_D=c_D/n. The reviewed
uniform energy compatibility gives |u_D-2c/S|<=2sqrt(delta).
Its source is
`original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`;
the compatibility estimate itself uses delta, not small gamma.
Thus (1.1) also holds with U replaced by this actual u_D, with the
same stated error order. On a general shell c/Phi(K) must be retained;
the active saturation assumption is not silently imposed there.

The exact natural-D field and the measure of W_D are unchanged. The
separate reviewed metric-stability theorem has source
`original_mo_diagonal_majorizer_metric_stability.md`, SHA256
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.
It can therefore use this new delta-only normalization in its actual
weighted-cell upper, without interpreting the auxiliary K_I as a new
source or as a covariance replacement.

The cap Phi(A)/n^(3/2)<=1/2+o(1) and conditional optimality do not
themselves become a proof that delta tends to zero. The theorem is a
conditional normalization for the scalar/near-scalar diagonal branch,
including positive canonical gap. The full actual moment constraints
must still be coupled to the actual cross measure and all original
cells to evaluate the desired upper. No large-gap closure, all-cell
width bound, or original MO convergence is claimed.
