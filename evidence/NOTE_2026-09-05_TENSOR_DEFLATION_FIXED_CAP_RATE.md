# Fixed-cap obstruction to a fast tensor-deflation repair rate

2026-09-05. Analytic source-only theorem. No mathematical computation,
matrix search, numerical optimization, or signing census was run.

This note keeps the original signing A unchanged. It studies the Gaussian
cost of the positive-part repair of an n-scale tensor covariance. It does
not replace A by a regularized signing, does not assume that A is an exact
minimizer, and does not assert a result for the full coupled operator
A tensor A - S_B + I. Its obstruction applies to uniform estimates based
only on a fixed source-norm cap C>1/2.

## 1. The repair and the conclusions

Let A be a complete symmetric zero-diagonal signing of order n. Write

    Phi(A) = max_{x in {+-1}^n} |x^T A x|/2,
    beta(T) = max_{x,y in {+-1}^n} |x^T T y|,
    S_j = tr |A|^j,          L = ||A||op.

For a threshold K>0 and a choice sigma in {+1,-1}, define the positive
semidefinite covariance on ordered matrix coordinates

    R_{K,sigma} = [sigma (A tensor A)/(K n) - I]_{+}.          (1.1)

Here the positive part is spectral functional calculus. Thus
I-sigma(A tensor A)/(K n)+R_{K,sigma} is positive semidefinite. Let G_R be
the centered real n by n Gaussian matrix whose vectorized covariance is R.

The elementary exact bounds are

    rank R <= min(n^2, (n-1)^2/K^2, S_3^2/(K^3 n^3)),
    tr R <= min((n-1)^2/(4K^2), 4S_3^2/(27K^3 n^3)),         (1.2)
    ||R||op <= (L^2/(Kn)-1)_{+}.

The reviewed norm-only cubic bound gives S_3=O_C(n^(5/2)) under
Phi(A)<=C n^(3/2), so (1.2) yields rank R,tr R=O_C(n^2/K^3).
This trace bound does not by itself give the needed Boolean Gaussian norm.

The main result is a rate obstruction using actual sign sources:

For every fixed C>1/2 there are a constant c_C>0, an unbounded sequence of
thresholds K, and, for each such K, arbitrarily large orders n and actual
signings A with

    Phi(A) <= C n^(3/2),
    E beta(G_{R_{K,+}}) >= c_C n^(3/2)/sqrt(K),
    E beta(G_{R_{K,-}}) >= c_C n^(3/2)/sqrt(K).               (1.3)

Consequently no uniform estimate

    E beta(G_{R_{K,sigma}}) <= M_C n^(3/2)/K^a

can hold under that fixed cap for any a>1/2. The same obstruction applies
to an asymptotic-in-n estimate uniform in A and K, by first fixing K and
letting n grow in the constructed family. In particular an exponent a>1,
which would make this additive repair cost smaller than a 1/K signal,
cannot follow from the sole fixed cap.

Section 7 gives the corresponding same-order, symmetric zero-diagonal
Gaussian quadratic-norm statement, with an absolute constant-factor loss.
The result is not an impossibility theorem for a directly analyzed clipped
law, an adaptively coupled choice of source slack and K, or an exact-source
minimizer. It only rules out the displayed uniform repair estimate.

## 2. Exact rank, trace, and diagonal bounds

The tensor eigenvalues are lambda_i lambda_j. For t>=0 one has

    (t-1)_{+} <= t^2/4,
    (t-1)_{+} <= 4t^3/27.                                  (2.1)

The respective maxima of (t-1)/t^2 and (t-1)/t^3 for t>=1 occur at
t=2 and t=3/2. Counting tensor products with absolute value greater
than Kn, and using S_2=n(n-1), proves (1.2).

The same scalar inequality, applied spectrally, gives the Loewner bound

    0 <= R_{K,sigma} <= (A^2 tensor A^2)/(4K^2 n^2).          (2.2)

Since every diagonal entry of A^2 is n-1, every ordered-coordinate
variance obeys

    (R_{K,sigma})_{ij,ij} <= (n-1)^2/(4K^2 n^2) <= 1/(4K^2).
                                                                    (2.3)

These statements are valid for both tensor signs and every actual A.
They do not use, or imply, a small operator norm for A/sqrt(n).

For comparison, the immediate Frobenius estimate is only
E beta(G_R)<=n sqrt(tr R). Even with the cubic-moment bound this is
O_C(n^2/K^(3/2)), not a uniform n^(3/2)-scale estimate.

## 3. A finite template with one large mode and many untouched bulk modes

Put F=J_4-2I_4. This symmetric sign matrix satisfies

    F^2=4I_4,     F 1_4=2 1_4,     diag F=-1.

Fix a dyadic amplitude a=2^(-h), with h a nonnegative integer. For an
integer t sufficiently large, set

    k=16^t,          s=a k^(3/4)=2^(3t-h),
    q=s/sqrt(k)=2^(t-h),       K=q/4.                        (3.1)

Choose t large enough that q>=4 and s<=k/4. Then s is an integer and
K tends to infinity with t. Define

    H=F^{tensor 2t}.

It has order k, constant diagonal one, H^2=kI, and H1=sqrt(k)1.
Let S be any subset of [k] of size s, and overwrite the S by S principal
block of H by J_s, leaving every other entry unchanged. Call the resulting
matrix C_k. It is symmetric, all its entries are signs, and its diagonal
is still one.

Write P_S for the coordinate projection onto S and v=1_S/sqrt(s). Then

    C_k = J_S + E,
    E = H-P_S H P_S,             ||E||op <= 2 sqrt(k).        (3.2)

Its largest eigenvalue lambda satisfies lambda>=s, because
v^T C_k v=s. Choose a unit eigenvector u for lambda, with v^T u>=0.
Projecting C_k u=lambda u orthogonally to v gives

    ||(I-vv^T)u|| <= ||E||op/lambda <= 2/q.

Hence

    1_S^T u = sqrt(s) v^T u
      >= sqrt(s) sqrt(1-4/q^2) >= sqrt(3s)/2.               (3.3)

This is the needed localized cube overlap. It does not assume that the
entire eigenvector is flat or supported on S.

There are also many exact, unperturbed bulk eigenvectors. Let

    V = ker(H-sqrt(k)I) intersect {z:z_S=0},

and let Q be its orthogonal projection. Every z in V is also a
sqrt(k)-eigenvector of C_k, since the perturbation is supported on S by S.
The positive eigenspace of H has dimension (k+sqrt(k))/2, by its trace.
Imposing z_S=0 has codimension at most s. Therefore

    rank Q >= (k+sqrt(k))/2-s >= k/4,
    Q1_S=0,               Qu=0.                             (3.4)

The last identity follows from symmetry and lambda>=s>sqrt(k). Each
diagonal entry of an orthogonal projection lies in [0,1], so

    sum_{j=1}^k sqrt(Q_jj) >= tr Q >= k/4.                  (3.5)

No uniform-diagonal assertion about Q is made.

## 4. Amplification gives actual signings with a fixed original-norm cap

Independently let m=16^j with j>=1, put H_m=F^{tensor 2j}, and set

    n=mk,                  A=H_m tensor C_k - I_n.           (4.1)

The tensor has constant diagonal one and every entry is a sign. Thus A
is an actual complete symmetric zero-diagonal signing of order n.

We verify its norm cap without using ||C_k||op as an upper bound. For a
real array T let tau(T) be its real vector SDP with unit row vectors.
Since ||H||op=sqrt(k), tau(H)<=k^(3/2). The perturbation C_k-H is
supported on s^2 entries and every entry has magnitude at most two, so

    tau(C_k) <= tau(H)+sum_ij |(C_k-H)_ij|
      <= k^(3/2)+2s^2 = (1+2a^2) k^(3/2).                 (4.2)

For Boolean vectors x,y of order mk, group their m coordinates according
to the template coordinate, writing x_b,y_c in {+-1}^m. The vectors
x_b/sqrt(m) and H_m y_c/m all have norm one, and

    x_b^T H_m y_c = m^(3/2) <x_b/sqrt(m),H_m y_c/m>.

Consequently

    beta(H_m tensor C_k) <= m^(3/2) tau(C_k).

Applying this bound to x=y and restoring the subtracted identity gives

    Phi(A) <= (1/2+a^2)n^(3/2)+n/2.                        (4.3)

Given any fixed cap C>1/2, choose the fixed dyadic a so that
1/2+a^2<C. For all sufficiently large m, (4.3) proves Phi(A)<=C n^(3/2).
The amplitude and cap are independent of K; m may grow arbitrarily for
each fixed k, hence for each fixed K.

## 5. Both tensor signs contain an expensive high-by-bulk subcovariance

Let

    P_m^+ = (I+H_m/sqrt(m))/2,
    P_m^- = (I-H_m/sqrt(m))/2.

These are orthogonal projections. The vector 1_m lies in the plus
space, and their constant diagonals are (1+1/sqrt(m))/2 and
(1-1/sqrt(m))/2. Because m>=16, both diagonals are at least 1/4.

Define projections on R^n by

    P = P_m^+ tensor uu^T,
    Q_+ = P_m^+ tensor Q,
    Q_- = P_m^- tensor Q.

They are spectral projections of A, with respective eigenvalues

    alpha = sqrt(m)lambda-1,
    b_+ = sqrt(mk)-1,
    b_- = -sqrt(mk)-1.                                     (5.1)

The P space is orthogonal to either Q space. Both sqrt(m)lambda and
sqrt(mk) are at least four. Therefore

    alpha >= (3/4)sqrt(m)lambda,
    b_+ >= (3/4)sqrt(mk),
    |b_-| >= (3/4)sqrt(mk).

Using lambda>=s, q=s/sqrt(k), K=q/4, and n=mk gives

    alpha b_+/(Kn) >= (9/16) lambda/(K sqrt(k)) >= 9/4,
    alpha |b_-|/(Kn) >= 9/4.                               (5.2)

Thus the eigenvalue of the corresponding positive-part covariance on
each indicated product space is at least 5/4, in particular at least one:

    R_{K,+} >= P tensor Q_+,
    R_{K,-} >= P tensor Q_-.                               (5.3)

These are genuine Loewner inequalities between commuting spectral
projections and the spectral positive part. No general compression
inequality for the positive-part function is being assumed.

## 6. A fixed cube vector proves the Boolean Gaussian lower bound

Fix either Q_epsilon=Q_+ or Q_-. For an n by n matrix Z of independent
standard Gaussians, the matrix

    G_0 = P Z Q_epsilon

has vectorized covariance P tensor Q_epsilon, with the row/column
vectorization convention used consistently. By (5.3), G_R can be
realized as G_0 plus an independent centered Gaussian. Convexity of beta
and conditional Jensen give

    E beta(G_R) >= E beta(G_0).                             (6.1)

The cube vector x_0=1_m tensor 1_S has coordinates zero or one, and

    ||P x_0|| = sqrt(m) (1_S^T u) >= sqrt(3ms)/2.

The maximum defining beta is unchanged if either Boolean cube is
replaced by [-1,1]^n, because a bilinear objective attains a maximum at
vertices. We may therefore hold x_0 fixed. The vector G_0^T x_0 is
Gaussian with covariance ||P x_0||^2 Q_epsilon. It follows exactly that

    E beta(G_0) >= E||G_0^T x_0||_1
      = sqrt(2/pi) ||P x_0|| sum_i sqrt((Q_epsilon)_ii).

By the lower diagonal bound for P_m^epsilon and (3.5),

    sum_i sqrt((Q_epsilon)_ii) >= (m/2)(k/4)=mk/8.

Consequently

    E beta(G_R) >= sqrt(6/pi) (mk sqrt(ms))/16
      = sqrt(6/pi) n^(3/2) sqrt(s/k)/16.

The exact parameter identities give s/k=a^2/q=a^2/(4K). Thus

    E beta(G_{R_{K,sigma}})
      >= [a sqrt(6/pi)/32] n^(3/2)/sqrt(K)
      >= (a/32) n^(3/2)/sqrt(K),       sigma in {+1,-1}.    (6.2)

The last inequality uses only pi<6. This proves (1.3), with a positive
constant depending solely on the selected fixed cap.

For any putative exponent a_exp>1/2 and constant M_C, dividing (6.2)
by n^(3/2) would force (a/32)K^(a_exp-1/2)<=M_C for an unbounded
sequence of K, an impossibility. The order of limits is explicit: choose
the fixed amplitude, then K, then arbitrarily large m and n.

## 7. Same-order symmetric Gaussian quadratic norms

The tensor operator and its spectral positive part preserve the real
symmetric matrix space. Give that space the Frobenius inner product,
and let G_R^sym have covariance equal to this restriction of R. Both
product spaces P tensor Q_epsilon and Q_epsilon tensor P occur with the
same tensor eigenvalue. Their symmetric part has the genuine Gaussian
representative

    G_0^sym = (P Z Q_epsilon + Q_epsilon Z^T P)/sqrt(2).

It is a subcovariance of G_R^sym. In addition Q_epsilon x_0=0, by (3.4).
Hence G_0^sym x_0=Q_epsilon Z^T P x_0/sqrt(2), and the calculation of
Section 6 proves

    E beta(G_R^sym) >= (a/(32 sqrt(2))) n^(3/2)/sqrt(K).    (7.1)

This covariance convention is on orthonormal symmetric-matrix
coefficients; an off-diagonal matrix entry is its coefficient divided
by sqrt(2). The normalization of G_0^sym above uses precisely that
convention.

Let G_R^0 be G_R^sym with its matrix diagonal removed. Bound (2.3)
applies also to the diagonal coordinates of the symmetric restriction,
so

    E sum_i |(G_R^sym)_ii| <= n/(K sqrt(2 pi)).

For every zero-diagonal symmetric matrix T, cube polarization gives
beta(T)<=4Phi(T). The diagonal removal and (7.1) therefore yield

    E Phi(G_R^0)
      >= [a/(128 sqrt(2))] n^(3/2)/sqrt(K)
          - n/(4K sqrt(2 pi)).                             (7.2)

For each fixed K the last term is lower order as n grows. Thus the
failure of a uniform exponent greater than 1/2 is also present in a
genuine same-order, zero-diagonal quadratic Gaussian norm, not just
in an unrestricted Frobenius bound or a sphere relaxation.

## 8. Scope of the obstruction

The construction proves a lower bound for the Gaussian repair itself.
It does not assert that the norm of a full repaired covariance law is
the sum of an unrepaired norm and the repair norm. Before repair, the
n-scale covariance can be indefinite, so no Gaussian for that
indefinite object is invoked.

The source cap is any fixed C>1/2. The theorem says nothing at C=1/2,
does not establish that the constructed A are exact original
minimizers, and does not prevent a useful argument in which a source
slack tends to zero in a controlled relation to K. The constant in
(6.2) tends to zero with the chosen dyadic amplitude.

Finally, no cross block B has been specified. A conclusion for the
coupled positive-part repair of A tensor A - S_B + I requires a
separate analysis of that actual cross perturbation. The present
source-only theorem must not be transferred to it without that proof.
