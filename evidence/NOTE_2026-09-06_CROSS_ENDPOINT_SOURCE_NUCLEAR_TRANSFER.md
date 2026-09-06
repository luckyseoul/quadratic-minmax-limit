# Actual cross-endpoint law transfers to an original source nuclear cap

2026-09-06. Supporting analytic theorem for ACTUAL complete paired signings.
No mathematical computation, solver, construction, scan, or search was run.

The FULL actual weighted cross law is the input. No internal spectral
shape is assumed. The proof first bounds the nuclear moment of the FULL
actual weighted internal matrix, without assuming commutation. Only then
is one common large principal submatrix of the SAME original A used for
an original quadratic-norm comparison. It does not replace the paired
covariance, cross block, active field, or any optimizer.

## 1. Hypotheses and conclusion

Let n tend to infinity and let

    K=[[A,B],[B^T,-A]]

be a complete symmetric zero-diagonal signing of order 2n. Thus A is a
complete symmetric zero-diagonal signing and every entry of B is a sign.
Let D=diag(D_L,D_R)>0 be diagonal with D+-K>=0, and write

    dbar=tr(D)/(2n),
    delta=tr(D)tr(D^(-1))/(2n)^2-1,
    H=D_L^(-1/2) A D_L^(-1/2),
    W=D_L^(-1/2) B D_R^(-1/2),
    Phi(A)=max_(x in {+-1}^n)|x^T A x|/2.

Assume delta->0 and the empirical law of ALL n squared singular
values of the actual W converges weakly to

                         (1-m)delta_0+m delta_1,             (1.1)

where 0<m<=1 is fixed. Then in fact m<=1/2. There is one common
original index set J, of order q=|J| with q/n->1, such that

    dbar/sqrt(q) -> 1/sqrt(m),
    limsup ||A_J||op/sqrt(q) <= 1/sqrt(m),
    limsup tr|A_J|/q^(3/2) <= sqrt(1-m),
    Phi(A_J)/q^(3/2) <= Phi(A)/n^(3/2)+o(1).                (1.2)

Consequently m>=9/25 implies the actual source caps

    limsup ||A_J||op/sqrt(q)<=5/3,
    limsup tr|A_J|/q^(3/2)<=4/5.                            (1.3)

No separate trace cap, diagonal optimality, small canonical-primal gap,
pure-cross active state, or internal flat-law hypothesis is required.
The proof also establishes the full weighted internal moment facts

    (1/n)tr H^2 -> m,
    limsup (1/n)tr|H| <= sqrt[m(1-m)].                      (1.4)

This transfer theorem does not itself prove an additional original
source-norm lower bound. Its role is to supply (1.3) to a separately
established theorem that uses those two actual source caps.

## 2. Exact inverse-mean identity determines the source scale

Feasibility makes T=D^(-1/2) K D^(-1/2) a contraction. In particular
the actual W is a contraction, so its squared singular values lie in
[0,1]. Let

    m_D=(1/n)tr(WW^T),       m_0=n/dbar^2,
    t_i=d_i/dbar,
    ell=(1/n)sum_(i=1)^n 1/t_i,
    h=(1/n)sum_(i=1)^n 1/t_(n+i).

Literal complete cross sign squares give EXACTLY

    m_D=m_0 ell h,       ell+h=2(1+delta).                  (2.1)

To control the product, put
L=(1/n)sum_(i=1)^n t_i and R=(1/n)sum_(i=1)^n t_(n+i).
Then L+R=2. Cauchy--Schwarz gives L ell>=1 and R h>=1.
Therefore

    1/ell+1/h<=2,
    1+delta=(ell+h)/2<=ell h<=(1+delta)^2.                  (2.2)

Equivalently,

    m_D/(1+delta)^2<=m_0<=m_D/(1+delta).                    (2.3)

Because the actual cross spectra have common compact support, (1.1)
implies m_D->m. Thus m_0->m, proving dbar/sqrt(n)->1/sqrt(m).
The assumed positive endpoint mass supplies the scale; no separate
bound on tr(D) is needed.

## 3. Common good labels and the full internal second moment

The exact diagonal dispersion identity is

                 (1/(2n))sum_i (t_i-1)^2/t_i=delta.

For delta>0 choose epsilon=delta^(1/3), eventually at most 1/8, and set

    J={i in {1,...,n}: |t_i-1|<=epsilon,
                            |t_(n+i)-1|<=epsilon},
    q=|J|,       a=q/n,       b=1-a.

Outside the good interval, (t-1)^2/t>=epsilon^2/(1+epsilon). Each
excluded original label has at least one bad coordinate, so

             b<=2delta(1+epsilon)/epsilon^2 ->0.            (3.1)

If delta=0 take epsilon=0 and every original label. This convention
also handles sequences containing both zero and positive dispersion.
The same literal A_J is retained in both diagonal halves.

Let H_J=H[J,J] be the actual principal compression and
Q=diag(sqrt(t_i):i in J). Principal feasibility gives

    ||H||op<=1,       ||H_J||op<=1,
    A_J/dbar=Q H_J Q,
    ||A_J||op<=(1+epsilon)dbar,
    ||A_J/dbar-H_J||op<=3epsilon.                           (3.2)

The last bound follows from ||Q-I||<=epsilon and
||Q||<=sqrt(1+epsilon). No norm bound on the full untrimmed A/dbar
is inferred. Equation (2.3) and a->1 now give the scale and operator
claims in (1.2).

For completeness, write M_H=(1/n)tr H^2 and M_J=(1/q)tr H_J^2.
Positive and negative squared eigenvalue sums each lose between zero
and n-q under principal compression, by interlacing on [-1,1]. Thus

    0<=M_H-a M_J<=2b,       |M_H-M_J|<=2b.

The square function is Lipschitz with constant at most 3 on the common
spectral interval [-(1+epsilon),1+epsilon]. Weyl comparison in (3.2)
and the EXACT identity tr A_J^2=q(q-1) therefore give

              |M_H-(q-1)/dbar^2|<=2b+9epsilon.              (3.3)

Since q/n->1 and n/dbar^2->m, this proves M_H->m. The actual block
row inequality H^2+WW^T<=I then yields M_H+m_D<=1. Passing to the
limit proves the asserted feasibility restriction 2m<=1.

The good-coordinate and congruence mechanism was established in
`original_mo_near_scalar_diagonal_spectral_normalization.md` (280 lines),
SHA256 `c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2`.
The common-label version and the normalized 2b moment estimate also
appear in `original_mo_near_scalar_internal_flat_law_transfer.md`
(141 lines), SHA256
`f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`.
They are restated above; the new input here is the cross endpoint law,
not the earlier source's prescribed full internal spectral shape.

## 4. Endpoint cross mass bounds the full internal nuclear moment

Fix 0<tau<1 and define the spectral projection of the ACTUAL cross
matrix

    P_tau=1_(WW^T>=1-tau),       r_tau=rank(P_tau),
    p_tau=r_tau/n.

For each fixed tau, (1.1) gives p_tau->m; the threshold 1-tau
lies strictly between the two limiting atoms. The actual block row
inequality gives, without any commutation assumption on H,

    tr(H^2 P_tau)<=tr[(I-WW^T)P_tau]<=tau r_tau.            (4.1)

Use only the trace-norm triangle inequality and the rank--Frobenius
inequality on the matrix decomposition

                       H=H(I-P_tau)+H P_tau.

The first summand has rank at most n-r_tau and squared Frobenius
norm at most tr H^2. The second has rank at most r_tau and squared
Frobenius norm tr(H^2 P_tau). Hence the finite estimate is

    (1/n)tr|H|
       <=sqrt[(1-p_tau)M_H]+p_tau sqrt(tau).                (4.2)

This is a matrix spectral-moment estimate only. No original Boolean
objective is changed into a rectangular norm, and no subspace or
compression is substituted into the paired field.

First let n tend to infinity in (4.2), using M_H->m and p_tau->m.
Then let tau decrease to zero. This proves

                 limsup (1/n)tr|H|<=sqrt[m(1-m)],

the second assertion in (1.4). The two limits are taken in the stated
order; no unproved convergence rate at a moving spectral cutoff is used.

## 5. Transfer the nuclear cap and retain the original norm

Principal compression decreases the UNNORMALIZED trace norm of a
symmetric matrix, by positive/negative eigenvalue interlacing. Thus
tr|H_J|<=tr|H|. The absolute value function is 1-Lipschitz, so
Weyl comparison in (3.2) also gives

    (1/q)tr|A_J/dbar|
       <=(1/q)tr|H_J|+3epsilon
       <=a^(-1)(1/n)tr|H|+3epsilon.

Multiplying by dbar/sqrt(q)->1/sqrt(m), and using Section 4, proves

                  limsup tr|A_J|/q^(3/2)<=sqrt(1-m).

For each fixed Boolean signing on J, independent unbiased extension
to the other original labels leaves its expected original quadratic
energy equal to the A_J energy. Every crossing or removed internal
edge has zero expectation. Therefore

    Phi(A_J)<=Phi(A),
    Phi(A_J)/q^(3/2)<=a^(-3/2)Phi(A)/n^(3/2).                (5.1)

The normalized original norm on the right is bounded: feasibility
gives Phi(A)<=tr(D_L)/2<=n dbar, and dbar/sqrt(n)->1/sqrt(m).
Since a->1, (5.1) proves the final claim of (1.2).

For m>=9/25 the elementary inequalities
1/sqrt(m)<=5/3 and sqrt(1-m)<=4/5 give (1.3). They require no
additional full internal spectral law. This is precisely the extra
premise removed relative to the earlier 141-line transfer.

This note does not assert that arbitrary active sequences have small
dispersion or an endpoint cross law. It supplies a conditional actual
source transfer, not an all-profile exclusion or an original MO closure.
Any further original-source lower bound applied to (1.3) must be
established and reviewed separately.

The exact worker derived and authored the cross-to-nuclear transfer.
Root requested the bounded full proof after receiving the argument.
The older good-coordinate machinery is explicitly credited above.
During derivation no canonical repository file was changed and no
publication was performed.
