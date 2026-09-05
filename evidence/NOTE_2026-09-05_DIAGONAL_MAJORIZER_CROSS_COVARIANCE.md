# Diagonal-majorizer cross covariance with the original source unchanged

2026-09-05. Analytic covariance and Gaussian-reduction theorem. No
mathematical computation, signing census, simulation, or optimization was
run for this proof. An SDP optimum is used as an existence object, not as
a numerical calculation.

The construction keeps both internal blocks A and -A, the cross drift B,
and the original conditional objective unchanged. It uses a diagonal
majorizer of the literal complete signing K, rather than dividing by its
largest tensor eigenvalue. It proves an admissible local-scale covariance
and a same-source Gaussian reduction. An evaluated upper comparison with
the original optimum is not proved here.

## 1. Literal complete source and existence of a diagonal majorizer

Let n>=2, N=2n, and let

    K = [[A,B],[B^T,-A]],

where A is a symmetric zero-diagonal order-n signing and B is an n by n
sign matrix. Thus K is an actual complete symmetric zero-diagonal signing
of order N. Write L=||K||op and

    Phi(K) = max_{z in {+-1}^N} |z^T K z|/2.

For a real square array T write beta(T)=max_{x,y}|x^T T y| and let tau(T)
be its real vector SDP with unit vectors on both sides. The already proved
finite tensor-rounding inequality is

    tau(T) <= g beta(T),       g=pi/[2 log(1+sqrt(2))].       (1.1)

This is the same finite real-array inequality used in the reviewed
cross-only and same-order regularization notes. Cube polarization for a
symmetric zero-diagonal matrix gives beta(K)<=4Phi(K).

There is an attained positive diagonal D=diag(d_1,...,d_N) such that

    D-K >= 0,          D+K >= 0,
    S:=tr D=tau(K) <= 4g Phi(K).                            (1.2)

Here and below matrix inequalities are Loewner inequalities. To verify the
same diagonal in (1.2), the bipartite SDP dual for tau(K) has a feasible
block matrix [[D_r,-K],[-K,D_c]], with objective
(tr D_r+tr D_c)/2. Since K is symmetric, swapping its two blocks preserves
feasibility. Average the two feasible matrices and set D=(D_r+D_c)/2.
The resulting block matrix [[D,-K],[-K,D]] is orthogonally equivalent to
diag(D-K,D+K), and tr D is exactly the dual objective. Strong duality and
attainment follow from strict primal and dual feasibility and compact
nonnegative diagonal trace sublevels. No diagonal entry of an attained
feasible D can be zero, because a zero diagonal of a PSD block matrix
would force the corresponding nonzero row of K to vanish.

More generally everything below applies to ANY positive diagonal D with
D+-K>=0 and

    S <= C_D N^(3/2),                                      (1.3)

where C_D is fixed. A norm cap on Phi(K) supplies such a D by (1.2).
Conversely (1.3) already gives Phi(K)<=S/2.

D and K are selected before drawing any Gaussian or rounded signs.

## 2. Exact local lower bounds and the weighted contraction

Put

    T=D^(-1/2) K D^(-1/2).

Then T is symmetric, has zero diagonal, and ||T||op<=1. In particular
the i-th diagonal entry of T^2 is at most one, so

    d_i >= sum_{j != i} 1/d_j
      >= (N-1)^2/(S-d_i) >= (N-1)^2/S.                    (2.1)

Summing the first inequality also gives

    tr D^(-1) <= S/(N-1).                                 (2.2)

For every distinct i,j one has d_i d_j>1. Indeed a value at most one
would make the j-th term in the i-th row sum for T^2 at least one,
and there are N-2 other strictly positive terms. Here N>=4.

On cross coordinates e=(i,j), representing the edge {i,n+j}, define

    q_ij=1/sqrt(d_i d_(n+j)),       Q=diag((q_ij)_{ij}),
    epsilon=max_ij q_ij^2 < 1.

The trace-cap bound gives the quantitative estimate

    epsilon <= S^2/(N-1)^4 = O_{C_D}(1/N).                 (2.3)

Thus local cross-edge scales are of order at least sqrt(N), regardless
of whether L/sqrt(N) is bounded. No upper bound on individual d_i at
the sqrt(N) scale is asserted.

## 3. Symmetric-edge compression and exact unit variance

On the N^2 ordered matrix coordinates, I-T tensor T is PSD and at most
2I. Compress it to the orthonormal symmetric cross-edge matrices

    F_ij=(E_(i,n+j)+E_(n+j,i))/sqrt(2).

The compressed T tensor T entry is exactly

    T_ik T_(n+j,n+l) + T_(i,n+l) T_(n+j,k)
      = q_ij q_kl [-A_ik A_jl+B_il B_kj].                  (3.1)

Define the self-adjoint exchange operator

    (S_B X)_ij=(B X^T B)_ij,
    H=A tensor A-S_B+I.

The tensor has entry A_ik A_jl on coordinates ij,kl, and H has diagonal
zero. The compressed covariance is therefore

    C_base = I+Q(A tensor A-S_B)Q,
    0 <= C_base <= 2I,
    (C_base)_(ij,ij)=1-q_ij^2.                             (3.2)

Restore the missing variances by an independent diagonal Gaussian:

    R_D = C_base+Q^2 = I+Q H Q.                            (3.3)

This is an exact correlation matrix:

    diag R_D=1,          0 <= R_D <= (2+epsilon)I < 3I.     (3.4)

For distinct cross coordinates,

    (R_D)_(ij,kl)=q_ij q_kl(A_ik A_jl-B_il B_kj),
    |(R_D)_(ij,kl)| <= 2epsilon=O_{C_D}(1/N).              (3.5)

In particular (3.3) is PSD before any Gaussian is invoked. It is not a
positive-part repair of an indefinite covariance.

For comparison, literal variance normalization of C_base is also valid.
With s_ij=sqrt(d_i d_(n+j)-1) its exact correlation matrix is

    R_norm=I+diag(1/s_ij) H diag(1/s_ij).

If epsilon<=1/2, conjugating C_base by
J=diag((1-q_ij^2)^(-1/2)) shows

    ||R_norm-R_D||op <= 6epsilon.                          (3.6)

Indeed ||J||<=sqrt(2), ||J-I||<=epsilon, ||C_base||<=2, and
||Q^2||<=epsilon. The separable construction (3.3) is sufficient here;
no change of source is involved in choosing it.

## 4. The shifted-sign law and the Gaussian models

Let G_D have covariance R_D on its actual n by n matrix entries, and
let W have independent standard Gaussian entries. For a deterministic
real h chosen before the disorder put

    B_h,ij=sign(G_D,ij+h B_ij),
    s_h=2 Phi_Gauss(h)-1,
    k_h=4 phi_Gauss(h)^2,
    v_h=1-s_h^2-k_h >= 0,
    Z_h=s_h B+sqrt(k_h)G_D+sqrt(v_h)W.                    (4.1)

The Gaussian matrices G_D and W are independent. Every B_h is an
actual complete cross signing, and E B_h=s_h B. For arbitrary fixed
real internal energies I(x,y) and |theta|<=1 define

    M_I(U)=max_{x,y in {+-1}^n}|I(x,y)+theta x^T U y|.

Under (1.3), with a constant depending only on C_D,

    |E M_I(B_h)-E M_I(Z_h)| <= C n^(16/11),               (4.2)

uniformly in h, I, theta, A, B, and the admissible D.

There is a second, slightly simpler Gaussian model with covariance
C_base instead of R_D. If G_base has covariance C_base, independently
of W, set

    Z_h^base=s_h B+sqrt(k_h)G_base+sqrt(v_h)W.

Then the same bound (4.2) holds with Z_h^base in place of Z_h. Only
the Gaussian model drops the variance padding; the actual sign law in
(4.1) always uses the unit-variance R_D. Sections 5-6 prove both claims.

## 5. Complete weighted Hermite remainder

Throughout the separated series decomposition in this section assume
epsilon<=1/2. Section 6 separately treats the remaining bounded set of
orders under a fixed C_D.

Let c_j(h) be the normalized Hermite coefficients of
sign(t+h)-s_h in standard Gaussian L2. Then c_1^2=k_h and
sum_{j>=2}c_j^2=v_h. Replacing h by -h multiplies c_j by
(-1)^(j+1). Thus threshold orientation by B_ij affects exactly the
even coefficients. Let b=vec B and D_b=diag b.

We use two real matrices indexed by all unordered edges of the actual
complete signing K. Let E_N be its line-graph adjacency matrix, and let

    (Q_N)_{ij,kl}=K_ik K_jk K_il K_jl.

The latter is zero for intersecting edges. Their algebraic operator
bounds, valid without any conference-scale norm assumption, are

    ||E_N||=2(N-2),        ||Q_N||<=L(N-1)/2.               (5.1)

For the second bound, square the full ordered-pair four-cycle matrix.
Its squared entries are (sum_r K_ir K_jr K_pr K_qr)^2>=0. With
V_(ij),r=K_ir K_jr one has V^T V=(K^2) circ (K^2)<=L^2(N-1)I
by the positive Schur map, and each row of V has squared norm at most
N-1. Every squared-matrix row sum is at most L^2(N-1)^2. Its operator
norm is consequently at most L(N-1). Compression to the orthonormal
symmetric off-diagonal basis gives twice Q_N, proving (5.1).

Let Pi select the cross edges. For every even integer p>=2, the
entrywise power of the raw zero-diagonal cross numerator H is exactly

    H^{circ p} = 2^(p-1) [11^T+Pi Q_N Pi^T]
      +(1-2^(p-1))Pi E_N Pi^T-2^(p-1)I.                  (5.2)

Indeed on disjoint edges the numerator has magnitude two exactly when
Q_N=1, and otherwise zero. On adjacent cross edges it has magnitude one;
on the diagonal it is zero. This checks all entries of (5.2).

For each p define Q_p=diag(q_e^p). Let C_h be the exact centered
covariance of B_h and C_0=k_h R_D+v_h I. The full even-Hermite
off-diagonal correction decomposes as P_h+E_even, where

    P_h=sum_{p>=2, p even} c_p^2 2^(p-1)
            (b circ q^p)(b circ q^p)^T >= 0,              (5.3)

and

    E_even=sum_{p>=2, p even} c_p^2 D_b Q_p
       Pi[2^(p-1)Q_N+(1-2^(p-1))E_N-2^(p-1)I]Pi^T
       Q_p D_b.                                         (5.4)

The diagonal cancellation between (5.3) and (5.4) is exact: all
diagonal Hermite mass was already placed in v_h I. The positive
covariance P_h may have more than one rank-one component; it is not
claimed to be a single scalar multiple of bb^T.

Since sum_p c_p^2<=1,

    sum_{p>=2, p even} c_p^2 2^(p-1) epsilon^p
      <= 2epsilon^2.

Therefore (5.1) and (5.4) give

    ||E_even|| <= epsilon^2[L(N-1)+4(N-2)+2].             (5.5)

The odd remainder starts at degree three. Its entrywise absolute value
is at most |(R_D)_(e,f)|^3, hence its row sum and operator norm are
at most 8n^2 epsilon^3. Combining the two remainders yields

    C_h=(C_0+P_h)+E_h,
    ||E_h|| <= 8n^2 epsilon^3
      +epsilon^2[L(N-1)+4(N-2)+2]
      =O_{C_D}((L+1)/N).                                 (5.6)

Both C_h and C_0+P_h are PSD. The error E_h need not be PSD.

The retained positive covariance has a small ACTUAL Boolean norm cost.
If Y_P is a centered Gaussian with covariance P_h, its coordinate
variance is at most 2q_e^4. The series (5.3) converges absolutely on
these finite coordinates and defines such a Gaussian covariance. Thus

    E beta(Y_P) <= E sum_e |(Y_P)_e|
      <= (2/sqrt(pi)) sum_ij q_ij^2
      = (2/sqrt(pi))(sum_{i<=n}1/d_i)(sum_{j>n}1/d_j)
      <= S^2/[2sqrt(pi)(N-1)^2] = O_{C_D}(N).             (5.7)

The last step uses (2.2) and the product bound for two nonnegative
partial sums. This is a norm estimate on all weighted rank-one
components together, not an inference from low rank alone.

## 6. Gaussian comparison, universality, and removal of padding

The Gaussian finite-maximum comparison with arbitrary deterministic
offsets says that a covariance operator error delta in dimension M
costs at most sqrt(2delta M log J) for J coefficient states of norm
at most sqrt(M). It follows by adding independent delta-I noise and
Gaussian convex order in both directions. Here M=n^2 and J<=2^(2n+1),
including the sign augmenting the absolute value. Equations (5.6)-(5.7)
therefore compare the matched Gaussian of B_h to Z_h with error

    O_{C_D}(n sqrt(L+1)+n).                               (6.1)

The mean-preserving shifted-sign universality theorem applies directly
to (4.1): (3.4) supplies the uniform latent covariance bound three,
the diagonal is one, and the observable is the augmented family
sigma theta(x_i y_j) with the arbitrary fixed internal energy absorbed
into its deterministic prior. Its growing-temperature estimate has
absolute error O(n^(16/11)), uniformly in the deterministic threshold h.
This is exactly the theorem used in the reviewed whole-edge and direct
cross reductions; its hypotheses have been verified anew here.

Finally L^2<=8Phi(K)<=4S, by the reviewed interpolation/polarization
bound and (1.3). Hence n sqrt(L+1)=O_{C_D}(n^(11/8)), which is smaller
than n^(16/11). This proves (4.2) for large N. The finitely many smaller
orders under a fixed C_D are absorbed by increasing the constant:
the signing norm is at most n^2 and the latent Gaussian covariance
is always bounded by three. No covariance formula is singular there.

For removal of the padding, write in distribution
G_D=G_base+W_q, independently, where (W_q)_ij=q_ij Z_ij and the Z_ij
are independent standard Gaussians. The variance of every fixed Boolean
pairing is

    V_q=sum_ij q_ij^2
      <= S^2/[4(N-1)^2]=O_{C_D}(N).

The standard Gaussian finite-maximum bound gives

    E beta(W_q) <= sqrt(4n log(2) V_q)=O_{C_D}(n).          (6.2)

Since sqrt(k_h)<=1, the actual objective is Lipschitz under this
addition with cost at most beta(W_q). Equation (6.2) proves the
claimed version of (4.2) for Z_h^base. C_base is PSD throughout;
there is no Gaussian associated to an indefinite formal covariance.

Equation (3.6) likewise costs O_{C_D}(n) in the expected Gaussian
maximum if one prefers literal variance normalization. It is not
needed for either stated model.

## 7. Exact weighted process covariance and the remaining compatibility

Set

    U=diag(d_1^(-1/2),...,d_n^(-1/2)),
    V=diag(d_(n+1)^(-1/2),...,d_(2n)^(-1/2)),
    A_L=U A U,          A_R=V A V,          W_D=U B V.

Thus T has blocks [[A_L,W_D],[W_D^T,-A_R]], with ||T||<=1. For two
Boolean pairs (x,y),(x',y'), the centered Gaussian cross process from
Z_h^base has the EXACT covariance

    (k_h+v_h)(x^T x')(y^T y')
      +k_h[(x^T A_L x')(y^T A_R y')
             -(x^T W_D y')(x'^T W_D y)].                  (7.1)

For Z_h, add the padding term

    k_h(x^T U^2 x')(y^T V^2 y').                           (7.2)

In particular (7.1) contains two differently weighted copies of A,
not the original unweighted source energy. The deterministic mean is
still s_h B and the internal energy is still Q_A(x)-Q_A(y). No change
to those quantities was made in obtaining the local-scale covariance.

The trace bound gives no uniform upper bound d_i=O(sqrt(N)) at each
coordinate, and no coercive lower comparison of the weighted forms in
(7.1) with the original Q_A is asserted. The earlier scalar attenuation
by a possibly large mu is absent; its removal is not itself an evaluated
Gaussian upper. Weighted source/drift compatibility remains to be proved.

## 8. Actual conditional-optimizer consequence and scope

For a prescribed original source A let

    F_A^*=min_B Phi([[A,B],[B^T,-A]]).

Choose any exact conditional minimizer B_*, form the literal K with
that B_*, and choose D before the Gaussian. The elementary independent
cross-sign estimate gives

    F_A^* <= 2Phi(A)+2sqrt(log(2)) n^(3/2).

Thus every fixed original-source norm cap supplies (1.3) with a fixed
constant via (1.2). Every rounded B_h in (4.1) is an admissible competitor
with those same internal blocks. Applying (4.2) with
I(x,y)=Q_A(x)-Q_A(y) gives the valid floor

    F_A^* <= E Phi([[A,Z_h^base],[(Z_h^base)^T,-A]])
                 +O_C(n^(16/11)).                        (8.1)

The identical statement holds with Z_h. This is conditional ORIGINAL
norm optimality, not pressure optimality or optimality of a modified
source. The signing A has never been trimmed, deflated, or replaced.

Equation (8.1) does not reverse to give an order upper bound. A useful
upper analysis must retain the weighted forms of (7.1), the unweighted
internal energy and drift, and their actual coupling. No conclusion
about the original MO convergence problem follows from this note alone.
