# Near-scalar diagonals retain the actual complete-cross gain

2026-09-05. Analytic transfer theorem. No mathematical computation,
solver, signing search, spectral evaluation, or numerical run was used.

The complete-cross gain below is expressed in the FULL ACTUAL weighted
cross spectrum. A balanced complete submatrix is only an auxiliary
lower bound on the original cross norm. It never replaces the original
source, cross covariance, active cell, or spectral measure in the result.

## 1. Actual setting and conclusion

Let n>=2, N=2n, and let

    K=[[A,B],[B^T,-A]]

be an actual complete zero-diagonal paired signing. Let D>0 be any
positive diagonal such that D-K>=0 and D+K>=0. Write

    S=tr D,       dbar=S/N,       delta=S tr(D^(-1))/N^2-1,
    T=D^(-1/2) K D^(-1/2),       ||T||op<=1,
    W=D_L^(-1/2) B D_R^(-1/2),
    beta(B)=max_(x,y in {+1,-1}^n) x^T B y.

Let nu be the empirical probability measure of the n squared singular
values of this ACTUAL W, including zeros, and put

    m=integral y dnu(y),       v_2=integral y^2 dnu(y),
    m_0=n/dbar^2,             kappa=2/pi,
    g_kappa=sqrt(kappa)-kappa>0.

The symbol g_kappa is only this fixed positive constant, not an SDP gap.
Assume the separate source-scale cap S<=C N^(3/2), for fixed C>0.
There is a nonnegative R_C(n) tending to zero such that, uniformly for
0<=delta<=1/512,

    beta(B)/(n dbar)
      >=kappa v_2/m+g_kappa m
           -[25kappa C^2+6g_kappa]delta^(1/3)-R_C(n).       (1.1)

The full actual moments satisfy 0<v_2<=m<=1 and

                    1/(2C^2)<=m_0<=m.                         (1.2)

No trace optimality of D, small canonical-primal gap, or actual
unweighted bound ||B||op=O(sqrt(n)) is assumed. In particular this
removes that extra unweighted operator hypothesis from the earlier
near-flat gain, within the stated near-scalar branch.

The Gaussian absolute-moment prerequisite is the complete 411-line
source `original_mo_complete_cross_flat_spectral_gain.md`, SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
Its actual complete-cross theorem is restated in Section 5. The
balanced trimming argument below is given fully, including its first
and fourth singular-moment comparisons.

## 2. Exact normalization and directional inverse means

The contraction row-square inequalities give

    (N-1)tr(D^(-1))<=S,       dbar^2>=N-1.

The second assertion follows from the first and Cauchy--Schwarz.
The source-scale cap gives dbar^2<=C^2 N. Hence

             1/(2C^2)<=m_0<=n/(2n-1)<=1.                       (2.1)

Put t_i=d_i/dbar. Exactly,

    (1/N)sum_i t_i=1,
    (1/N)sum_i 1/t_i=1+delta,
    (1/N)sum_i (t_i-1)^2/t_i=delta.                            (2.2)

Define directional inverse means

    ell=(1/n)sum_L 1/t_i,       h=(1/n)sum_R 1/t_i.

Literal squared cross entries give

    m=m_0 ell h,       ell+h=2(1+delta).                        (2.3)

Thus ell h<=(1+delta)^2 by arithmetic-geometric mean. For the other
direction let t_L and t_R be the respective averages of t_i in the
two halves. Cauchy--Schwarz gives ell>=1/t_L and h>=1/t_R;
also t_L+t_R=2, so t_L t_R<=1. Consequently

                      m_0<=m<=m_0(1+delta)^2.                   (2.4)

This is a comparison of actual complete-signing moments; no scalar
replacement for W is made. Since W is a cross block of the contraction
T, it is itself contractive, proving v_2<=m<=1. Positivity follows
because the complete cross block is nonzero.

## 3. Balanced good coordinates and the auxiliary complete block

For a finite comparison choose 0<eta<=1/2 with

                    b_0=delta(1+eta)/eta^2<=1/4.

Call a coordinate good if |t_i-1|<=eta. Outside that interval,
(t_i-1)^2/t_i>=eta^2/(1+eta). By (2.2), at most b_0 N
coordinates are bad in the two halves combined.

Let q be the smaller of the two good-coordinate counts. Choose q
good coordinates I in the left half and q good coordinates J in the
right half, removing extra good coordinates from the larger half.
Write

    a=q/n,       theta=1-a,       d'=(1+eta)dbar,
    B_J=B[I,J],       W_J=W[I,J],       Y=B_J/d',
    m'=q/d'^2=a m_0/(1+eta)^2.

The particular choices within the good sets may be deterministic and
arbitrary. They satisfy

    0<=theta<=2b_0,       1/2<=a<=1,
    q>=n/2,              m'>=1/(9C^2).                         (3.1)

The last bound uses (2.1) and (1+eta)^2<=9/4.

The auxiliary B_J is a genuine q by q COMPLETE sign matrix. Moreover

    B_J=D_I^(1/2) W_J D_J^(1/2),
    ||W_J||op<=1,       ||B_J||op<=d'.                         (3.2)

The norm inequality follows because all retained diagonal entries of
D are at most (1+eta)dbar. It bounds only the auxiliary B_J; the
original B may still have exceptional large singular values.

For any fixed left and right signing on I,J, extend it by independent
unbiased signs on the removed coordinates. The expected full cross
energy is its B_J energy, since every omitted term has zero mean.
Taking the full maximum therefore gives the exact original-norm bound

                              beta(B)>=beta(B_J).               (3.3)

No optimality or active-cell property of B_J is needed.

## 4. First and fourth moments return to the full actual W

For k=1,2 let v_k=integral y^k dnu(y); thus v_1=m. Let nu_Y
be the q-point empirical squared-singular-value measure of Y, and put

                         v_k'=integral y^k dnu_Y(y).

First compare W_J with W. Their symmetric dilations are

    L=[[0,W],[W^T,0]],       L_J=[[0,W_J],[W_J^T,0]].

The second is a principal submatrix of the first, after deleting
2(n-q) coordinates. Both are contractions. Their positive eigenvalues
are their cross singular values. Principal-eigenvalue interlacing and
monotonicity of (max(x,0))^(2k) show that

    0<=v_k-(1/n)tr[(W_J W_J^T)^k]<=2theta.                    (4.1)

More explicitly, deleting r coordinates from a symmetric contraction
loses between zero and r from the sum of positive eigenvalues to the
power 2k. Here r=2(n-q). Each dilation has paired eigenvalues, so its
positive-power sum is exactly the indicated cross singular-power sum.
Zero eigenvalues cause no difficulty in this comparison.

On the retained coordinates define

    P_L=diag(sqrt(t_i/(1+eta))),
    P_R=diag(sqrt(t_j/(1+eta))).

Every diagonal entry of either P lies in [1-eta,1]. Indeed
sqrt((1-eta)/(1+eta))>=1-eta. Thus

    Y=P_L W_J P_R,       ||Y||op<=1,
    ||Y-W_J||op<=||P_L-I||op+||P_R-I||op<=2eta.                (4.2)

Weyl's singular-value bound and the 2k-Lipschitz property of
x^(2k) on [0,1] now give

    |a v_k'-(1/n)tr[(W_J W_J^T)^k]|<=4k a eta,
    |a v_k'-v_k|<=2theta+4k a eta.                            (4.3)

In particular, literal squared signs in B_J give the exact first-
moment identity v_1'=m', and hence

    |a m'-m|<=2theta+4a eta,
    a v_2'>=[v_2-2theta-8a eta]_+.                            (4.4)

Here [x]_+=max(x,0). The moment v_2 is the FULL actual fourth
singular moment of W, not the fourth moment of a discarded or
independently supplied cross matrix. The directional identity (2.4)
is additionally available and is stronger than needed for the first
moment transfer in (4.4).

## 5. Apply the actual complete-cross theorem only as a norm lower

The prerequisite 411-line theorem states: for any complete q by q
sign matrix J, scale d'>=||J||op, and m'=q/d'^2>=mu>0,

    beta(J)/(q d')
       >=kappa+(sqrt(kappa)-kappa)m'-kappa epsilon'-e_mu(q),
    epsilon'=1-tr[(J^T J)^2]/(q^2 d'^2),

where e_mu(q) tends to zero uniformly over J,d'. We choose its
error function bounded and nonnegative, as permitted by its uniform
absolute-moment proof. With J=B_J, its epsilon' satisfies

                       1-epsilon'=v_2'/m'.                     (5.1)

Equations (3.1)--(3.3) therefore imply

    beta(B)/(n dbar)
       >=a(1+eta)[kappa v_2'/m'+g_kappa m'-e_mu(q)],
    mu=1/(9C^2).                                              (5.2)

Both factors a=q/n and 1+eta=d'/dbar are retained explicitly.
Using (4.4) gives the following finite transferred inequality:

    beta(B)/(n dbar)
       >=[kappa(1+eta)/m'][v_2-2theta-8a eta]_+
             +g_kappa a(1+eta)m'-a(1+eta)e_mu(q).              (5.3)

Every quantity on the right is either a controlled scalar comparison
parameter or a moment of the full original W.

## 6. One-sided coefficient comparison and explicit dispersion error

Since m'=a m_0/(1+eta)^2<=m_0<=m, the first term in (5.3)
is at least

    (kappa/m)[v_2-2theta-8a eta]_+
          >=kappa v_2/m-kappa(2theta+8a eta)/m.                (6.1)

The positive part makes this one-sided comparison valid even when
the bracket before truncation is negative; no subtraction is multiplied
by an uncontrolled enlarged coefficient.

The second term in (5.3) is exactly

                     g_kappa a^2 m_0/(1+eta).

Using (2.4), m_0<=1, and a=1-theta gives

    g_kappa m-g_kappa a^2 m_0/(1+eta)
      <=g_kappa m_0[(1+delta)^2-a^2/(1+eta)]
      <=g_kappa(2delta+delta^2+2theta+eta).                    (6.2)

To see the last bound, use 1-a^2<=2theta and
1-1/(1+eta)<=eta. Combining (5.3)--(6.2), we have

    beta(B)/(n dbar)
       >=kappa v_2/m+g_kappa m
          -kappa(2theta+8a eta)/m
          -g_kappa(2delta+delta^2+2theta+eta)
          -a(1+eta)e_mu(q).                                  (6.3)

For 0<delta<=1/512 choose eta=delta^(1/3). Then eta<=1/8,
b_0=eta(1+eta)<=9eta/8<1/4, and theta<=9eta/4. Thus

    2theta+8a eta<=25eta/2,
    2delta+delta^2+2theta+eta<=6eta,
    1/m<=2C^2,       a(1+eta)<=3/2.                           (6.4)

Define the uniform tail envelope

    R_C(n)=(3/2)sup_(integer j>=ceil(n/2)) e_(1/(9C^2))(j).

It is finite and tends to zero. It controls the last term of (6.3)
because q>=n/2. Substitution of (6.4) proves (1.1).

If delta=0, identity (2.2) forces D=dbar I. Take q=n, a=1,
eta=theta=0 and d'=dbar directly; there is no trimming. The earlier
complete-cross theorem applies to the original B, with nu_Y=nu and
m'=m=m_0. This proves (1.1) also at the exact scalar endpoint.
No rate in n beyond R_C(n)->0 is asserted or needed.

## 7. The actual weighted active-state ratio

For any actual Boolean pair x,y write

    p=x^T A x,       q_A=y^T A y,       c=x^T B y,
    c_D=x^T W y,     u_D=c_D/n.

The finite uniform compatibility bound is

                   |u_D-c/(n dbar)|<=2sqrt(delta).               (7.1)

For clarity it holds for every feasible D here, without its trace
being optimal. Let Q=diag(sqrt(t_i)). By (2.2),
sum_i(sqrt(t_i)-1)^2<=N delta. For any full Boolean z,
||Qz||=||z||=sqrt(N). Since K/dbar=Q T Q and T is contractive,

    |z^T(K/dbar-T)z|<=2N sqrt(delta),
    Phi(K/dbar-T)<=N sqrt(delta).

For a symmetric block matrix, compare the quadratic energies at
(x,y) and (x,-y): the cross pairing is bounded by its Phi norm.
The cross block here is B/dbar-W. Dividing the resulting bound
|c/dbar-c_D|<=N sqrt(delta) by n proves (7.1).

Now assume the actual original state has p=q_A=0 and positive
c=Phi(K). For every paired signing beta(B)<=Phi(K), by comparing
the full objective at (x,y) and (x,-y); the stated state also gives
beta(B)>=c. Hence c=beta(B)=Phi(K). From (1.1) and (7.1),

    u_D>=kappa v_2/m+g_kappa m
       -[25kappa C^2+6g_kappa]delta^(1/3)
       -2sqrt(delta)-R_C(n).                                 (7.2)

This applies equally to an actual within-cell representative of a final
original/weighted refined cell satisfying those original active-state
conditions. On general states the original ratio c/beta(B) cannot be
dropped, and (7.2) is not asserted there.

## 8. Consequence and exact remaining scope

Along any sequence of actual matrices with the fixed source-scale cap,
n tending to infinity, delta tending to zero, and

                         (m-v_2)/m tending to zero,

the original cross norm has

    beta(B)/(n dbar)>=kappa+g_kappa m-o_C(1).

Because m>=1/(2C^2), the gain above kappa is uniformly positive
at leading order. In an actual active pure-cross state the same is
true for u_D. In particular a full weighted cross law tending to
(1-m)delta_0+m delta_1 with a positive limiting m cannot coexist
with u_D=kappa+o(1) under these near-scalar, fixed-cap assumptions.

This excludes the previously displayed kappa-floor formal trace data
from this actual branch using complete-signing information. It does
not assert that the entire formal relaxation is realizable, and it
does not change or retract the earlier relaxation-only statement.

The final measure is always the original full nu. Original and
weighted energies are compared explicitly; the auxiliary B_J is not
inserted into the Gaussian field, ellipsoid functional, or active-cell
definition. No global bound on ||B||op or maximum d_i is claimed.
The theorem does not prove delta->0 for every conditional optimizer,
does not show that all active cells are pure-cross, and does not
evaluate the original ellipsoid upper at every remaining moment law.
The original all-cell inequality and convergence remain open.

## 9. Prior arguments and contribution record

The good-coordinate method was already justified for full symmetric
moments in the complete 280-line source
`original_mo_near_scalar_diagonal_spectral_normalization.md`, SHA256
`c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
The full 303-line compatibility source is
`original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
Both were read in full for this transfer. Their needed good-coordinate
and uniform-energy arguments have been rederived explicitly above;
the only substantive imported inequality is the frozen 411-line
complete-cross theorem cited in Section 1.

The root proposed applying that theorem to balanced good coordinates
and transferring its first and fourth moments to the actual W.
The optimized-profile proof worker derived the finite factors and
one-sided error bounds and authored this note. The exact worker
independently checked the proposed formulas before final-source review;
no optional stronger refinement supplied during that check is used
by this artifact. A separate full-source receipt is required for any
claim of completed review. No mathematical computation was run.
