# Full-SDP canonical gap controls weighted source compatibility

2026-09-05. Analytic conditional theorem, separate from the preceding
four-note weighted-covariance milestone. No mathematical computation,
signing census, numerical SDP, search, or optimization was run.

The gap below belongs to the vector SDP of the literal COMPLETE matrix K,
not to the Boolean norm and not just to its cross block B. The theorem
proves uniform weighted/unweighted energy compatibility when that gap is
small, without a maximum-diagonal bound. It does not prove that all
original or conditional minimizers have a small gap.

## 1. Setting and the quantitative conclusions

Let K be a complete symmetric zero-diagonal signing of order N>=3 and put

    q=N-1,                 S_3=tr |K|^3.

Let D=diag(d_1,...,d_N)>0 be ANY trace-optimal same-diagonal majorizer
for the real bipartite vector SDP of K:

    D-K>=0,       D+K>=0,       S=tr D=tau(K).

Existence, attainment, and the same-diagonal normalization are proved in
`original_mo_diagonal_majorizer_cross_covariance.md`, final SHA256
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
Write

    g=S-S_3/q,
    delta=S tr(D^(-1))/N^2-1,
    dbar=S/N,             T=D^(-1/2) K D^(-1/2).

Then g>=0, delta>=0, ||T||op<=1, and

    delta <= 4q ||D^(-1)||op g/N^2
      <= 4Sg/(qN^2).                                      (1.1)

In addition

    sum_i (sqrt(d_i)-sqrt(dbar))^2 <= S delta,
    Phi(K-dbar T) <= S sqrt(delta).                        (1.2)

Here Phi(E)=max_{z in {+-1}^N}|z^T E z|/2, also for real symmetric
zero-diagonal E. Define the normalized quantities

    eta=S/[N sqrt(q)] >= 1,          gamma=g/S in [0,1).

Then (1.1) reads

    delta <= 4 eta^2 gamma.                                (1.3)

Under a fixed S=O(N^(3/2)) cap, eta is bounded. Thus a vanishing
RELATIVE canonical gap gamma implies the uniform comparison (1.2),
with no assumption that D is scalar, uniformly close to scalar in
operator norm, or bounded above by a multiple of sqrt(N) coordinatewise.

For the paired complete matrix

    K=[[A,B],[B^T,-A]],            N=2n,

set U=diag(d_1^(-1/2),...,d_n^(-1/2)),
V=diag(d_(n+1)^(-1/2),...,d_(2n)^(-1/2)), and

    A_L=UAU,        A_R=VAV,        W_D=UBV.

For EVERY actual Boolean pair x,y, write

    p=x^T A x,       q_A=y^T A y,       c=x^T B y,
    p_D=x^T A_L x,   q_D=y^T A_R y,     c_D=x^T W_D y.

The uniform, finite bounds are

    |p_D-p/dbar| <= 2N sqrt(delta),
    |q_D-q_A/dbar| <= 2N sqrt(delta),
    |c_D-c/dbar| <= N sqrt(delta).                          (1.4)

In particular the weighted source energies are controlled on the actual
ORIGINAL zero-source slice p=q_A=0. Section 6 turns that fact into a
positive-Gaussian pure-cross field comparison, with normalized error
O_eta(gamma^(1/4)). Neither that field comparison nor (1.4) assumes
that the small-gap premise is already known for all minimizers.

## 2. The actual canonical primal and its weighted residual

Because K is a complete signing, every diagonal entry of K^2 is q.
The same holds for |K|^2. Consequently the rows of both K/sqrt(q)
and |K|/sqrt(q) are unit vectors. Their bipartite SDP objective is
S_3/q, proving g>=0.

For an explicit residual calculation put

    Z=(1/sqrt(q)) [[K],[|K|]],
    Dcal=diag(D,D),
    Qcal=[[D,-K],[-K,D]].

Here Z has 2N rows and N columns, and its Gram matrix has diagonal one.
The two majorizations of K imply

    0<=Qcal<=2Dcal,
    tr(Z^T Qcal Z)=2S-2S_3/q=2g.

Whitening by Dcal and squaring the resulting PSD contraction gives

    Qcal Dcal^(-1) Qcal <= 2Qcal.

Define the N by N residual matrices

    R_1=DK-K|K|,
    R_2=D|K|-K^2.

Substitution of Qcal Z=q^(-1/2)(R_1,R_2)^T proves

    ||D^(-1/2)R_1||_F^2+||D^(-1/2)R_2||_F^2 <= 4qg.       (2.1)

No unweighted residual estimate, maximum d_i, or operator bound on K
has been inserted into (2.1).

## 3. The inverse-weighted commutator identifies diagonal spread

Choose an orthogonal symmetric polar factor O for K, assigning either
sign on any zero eigenspace. Then K=O|K|, O commutes with K and |K|,
and O K^2=K|K|. Hence

    O R_2^T=KD-K|K|,
    [D,K]=DK-KD=R_1-O R_2^T.                              (3.1)

Let b=||D^(-1)||op. The two weighted Frobenius estimates needed for
(3.1) are

    ||D^(-1/2)R_1D^(-1/2)||_F^2
      <= b ||D^(-1/2)R_1||_F^2,

    ||D^(-1/2)O R_2^T D^(-1/2)||_F^2
      <= b ||R_2^T D^(-1/2)||_F^2
       = b ||D^(-1/2)R_2||_F^2.

The second estimate bounds the LEFT factor D^(-1/2) in operator norm
and uses orthogonality of O. It does not commute O through D.
The squared triangle inequality and (2.1) give

    ||D^(-1/2)[D,K]D^(-1/2)||_F^2 <= 8q b g.               (3.2)

On the other hand, K_ij^2=1 for i!=j, so the same squared norm is
exactly

    sum_{i,j} (d_i-d_j)^2/(d_i d_j)
      = 2[S tr(D^(-1))-N^2].                              (3.3)

The diagonal summands vanish, so adding them has not changed the
sum. Equations (3.2)-(3.3) prove the first part of (1.1).

The row-square bound for the contraction T gives, as in the covariance
note,

    d_i>=sum_{j!=i}1/d_j>=q^2/S,
    b<=S/q^2,
    tr D^(-1)<=S/q.                                      (3.4)

This proves the second part of (1.1). Cauchy--Schwarz also gives
S tr D^(-1)>=N^2; together with (3.4), it implies S>=N sqrt(q).
This verifies delta>=0 and eta>=1, and then (1.3) follows exactly.

The small-gap conclusion therefore controls an inverse-weighted spread,
not merely an unweighted variance multiplied by d_max+||K||op.

## 4. Uniform rescaling on the actual cube

Put t_i=d_i/dbar and let E denote uniform averaging over i. Then

    E t=1,
    E(1/t)=1+delta,
    E[(t-1)^2/t]=delta.

For every t>0,

    (sqrt(t)-1)^2 <= (t-1)^2/t.

Thus V:=sum_i(sqrt(d_i)-sqrt(dbar))^2<=S delta. For every real cube
vector z, let a=D^(1/2)z and b=sqrt(dbar)z. Their norms are at most
sqrt(S), and ||a-b||<=sqrt(V). Since K=D^(1/2)TD^(1/2),

    |z^T(K-dbar T)z|
      = |a^T T a-b^T T b|
      <= (||a||+||b||)||a-b||
      <= 2S sqrt(delta).                                (4.1)

Both K and T have zero diagonal. Taking the Boolean quadratic maximum
in (4.1) proves (1.2). The estimate is uniform over every cube state,
not averaged only over an SDP frame or over a Gaussian rounding.

## 5. Actual weighted source and cross energies

For the paired matrix put E=K-dbar T. Its diagonal blocks are

    E_L=A-dbar A_L,          E_R=-A+dbar A_R,

and its cross block is E_B=B-dbar W_D. All diagonal entries are zero.
For each fixed Boolean x, averaging the other block's coordinates as
independent fair signs gives

    |x^T E_L x|/2 <= Phi(E).

The same argument holds for E_R. For the cross block, comparing the
quadratic energies at (x,y) and (-x,y) gives

    |x^T E_B y| <= Phi(E).

These facts use zero-diagonal multilinearity and a block sign flip;
they do not require the two diagonal blocks of E to be negatives of
one another. Substituting Phi(E)<=S sqrt(delta) and dividing by
dbar=S/N proves all three inequalities in (1.4).

For mixed source pairings, ordinary cube polarization further gives

    |x^T(A_L-A/dbar)x'| <= 4N sqrt(delta),
    |y^T(A_R-A/dbar)y'| <= 4N sqrt(delta).                 (5.1)

The cross bound in (1.4) already applies to arbitrary mixed x,y'.
In terms of the relative canonical gap, the unprimed source bounds
in (1.4) are at most 4eta N sqrt(gamma), and the cross bound is at
most 2eta N sqrt(gamma).

## 6. Evaluated pure-cross field stability on an original zero-source cell

We use the actual weighted linear-field theorem in
`original_mo_diagonal_majorizer_weighted_shell_upper.md`, initial frozen
381-line version, SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
Its covariance normalization is the base model after the independently
controlled diagonal padding is removed. Its symbol d_0 denotes S;
the present dbar=S/N is the AVERAGE diagonal and is different.

Let 0<=k<=w=k+v<=1. Choose a nonempty actual cell with original energies
p=q_A=0. If weighted binning is used, choose its representative from
WITHIN THIS FINAL REFINED CELL, not from a larger weighted bin whose
other states might have different original energies. Such a choice
is permitted by the field note's comparison proof and binning error.

Let theta=(p_D,q_D,c_D) be that representative's actual weighted triple.
The field theorem gives the genuine PSD covariance

    M_theta=wnI_N
       +k[[q_D A_L,-c_D W_D],[-c_D W_D^T,p_D A_R]].        (6.1)

Define its pure-cross counterpart

    M_0=[[wnI_n,-k c_D W_D],[-k c_D W_D^T,wnI_n]].        (6.2)

It is also PSD: ||W_D||op<=1 and |c_D|<=n||W_D||op<=n, so
the off-diagonal block norm is at most kn<=wn. No scalar-D
assumption or Gaussian with an indefinite covariance is used.

Since the representative has ORIGINAL p=q_A=0, (1.4) gives

    |p_D|,|q_D| <= 2N sqrt(delta),
    ||M_theta-M_0||op <= 2kN sqrt(delta).                (6.3)

Let g_theta and g_0 be Gaussians with these covariances. For any
nonempty subset of the actual cell and any fixed identical offsets
a_z on its states, the Gaussian finite-maximum comparison gives

    |E max_z(a_z+g_theta^T z)-E max_z(a_z+g_0^T z)|
      <= 2sqrt(k log(2)) N^(3/2) delta^(1/4).            (6.4)

Indeed there are at most 2^N states, each has squared norm N, and
the covariance operator difference is bounded by (6.3). For an
absolute-value maximum the safe constant 2sqrt(2k log(2)) applies,
by including both augmenting signs and using N+1<=2N.

Combining (6.4) with (1.3), the non-augmented bound is at most

    2sqrt(2k eta log(2)) N^(3/2) gamma^(1/4).             (6.5)

Thus a bounded-eta, vanishing-relative-gap regime makes the internal
weighted-source contribution negligible on an actual zero-original-
source cell. The pure-cross covariance M_0 still uses the ACTUAL
weighted W_D and c_D. It is not replaced by an unweighted scalar
covariance, and (6.4) is not a final evaluation of its width.

The separate weighted-cell approximation and cell-selection errors
from the field theorem remain present if this corollary is applied
to binned, rather than exact, weighted cells. The representative
choice above ensures that the original zero-source premise is not
silently transferred from a cell to an unrelated representative.

## 7. Exact scope and remaining implication

The theorem applies to every actual trace-optimal D, including a
nonunique optimum and a K with singular spectrum. The polar factor
is completed orthogonally on the zero eigenspace; no inverse of K
or |K| was used.

The canonical gap g is a separate, nonnegative SDP slack. A fixed
original or conditional norm cap bounds eta but does not, in this
proof, force gamma=g/S to vanish. The result therefore provides a
quantified ACTUAL small-gap range, not uniform scalar compatibility
for every optimizer and not a proof of original convergence.

In particular small gamma does not imply a conference-scale operator
cap on K, justify a Gaussian with I-K tensor K/dbar^2 when that
matrix is indefinite, or permit deletion of the remaining weighted
cross metric. All Gaussian covariances invoked in Section 6 are
individually proved positive semidefinite.
