# A finite interior gap for complete Boolean quadratic forms

2026-09-23. Author-reviewed supporting theorem; original convergence OPEN.
No spectral cap, optimizer assumption, limiting law, or Gaussian closure
is used. This does not prove an order-extension or pressure comparison.

For a complete symmetric zero-diagonal signing A of order n>=2, write

    Q_A(z)=sum_(i<j) A_ij z_i z_j,
    F=Phi(A)=max_(x in {-1,1}^n) |Q_A(x)|.

The new finite estimate is: if z belongs to [-1,1]^n and at least k>=2
coordinates satisfy |z_i|<=1-epsilon, with 0<epsilon<=1, then

    F-|Q_A(z)| >= epsilon^2 k^3/(96 F).                 (1)

In particular, with s(z)=sum_i(1-|z_i|),

    s(z) <= 1+3[96 F(F-|Q_A(z)|)]^(1/3),                (2)
    F-|Q_A(z)| >= [(s(z)-1)_+]^3/(2592 F).              (3)

Thus, for ANY sequence with Phi(A_n)<=C n^(3/2), an o(n^(3/2))
norm deficit implies s(z_n)=o(n). Equivalently, rounding each coordinate
of z_n to a nearest sign changes it by o(n) in l1. This is a statement
about a near-extremal fractional spin vector, not a bound on the number,
Hamming separation, or covering radius of the Boolean maximizers.

## 1. Reused material and exact additional scope

Multilinearity gives the same one-sided extrema on the continuous cube
as on its vertices. Also, either one-sided extremum of a principal
submatrix is at most the corresponding full-matrix extremum: average
the omitted spins independently and uniformly. These are the existing
restriction facts in CORE, Section 2.

The September 6 archive already contains a cubic one-sided bound and
its macroscopic pairing consequence. Specifically,
`original_mo_broad_campaign_20260906/independent_and_hadamard.tar.gz`,
member `original-mo-campaign-root.Z0Um0rAh/UNPUBLISHED_WORKING_FRONTIER.md`,
records a spectral-bootstrap estimate on macroscopic induced sets under
a fixed normalized norm bound. That estimate and the pairing idea are
NOT new here; the member labels itself unreviewed scratch.

Sections 2--3 instead give a self-contained finite one-sided product
bound at every order and an interior-gap estimate at every subset size.
Section 4 applies the gap to the actual mean update, giving an additional
nonnegative term with no post-update distributional assumption. We do not
claim that a cubic pairing bound now supplies the linear posterior bound
missing in the archived Gaussian comparison.

The positive/negative trade-off was suggested by Bollobas and Scott,
[Discrepancy in graphs and hypergraphs, Theorem 1](https://people.maths.ox.ac.uk/~scott/Papers/disc.pdf).
Their induced-subgraph discrepancy theorem is not an input to the proof
below. The direct Boolean argument retains both one-sided extrema, so
no density-centering error or asymptotic spectral cutoff is needed.

## 2. A direct finite product bound

For an order-k complete signing B, k>=2, set

    P=max_x Q_B(x),   N=-min_x Q_B(x).

Both are positive because the uniform mean is zero and the variance is
binom(k,2)>0. Partition the vertices into S,T of sizes
a=ceil(k/2), b=floor(k/2). Sample X_T uniformly and put

    Y_i=sign(sum_(j in T) B_ij X_j),  i in S,
    L=a E|epsilon_1+...+epsilon_b|.

Ties in sign may be resolved either way. For any 0<=t<=1, the vector
(tY_S,X_T) is in the continuous cube. Its expected energy is

    t^2 E Q_(B[S])(Y_S) + t L + E Q_(B[T])(X_T)
       >= -t^2 N+t L.

The last internal expectation is zero; the principal lower bound is
-N. Therefore P>=t L-t^2 N. Applying the same argument to -B gives
N>=t L-t^2 P. At t=1 these imply P+N>=L. Let H=max(P,N), so H>=L/2.
Use t=L/(2H)<=1 in the inequality whose quadratic coefficient is -H:

    P N >= L^2/4.                                      (4)

For R=sum_(j=1)^b epsilon_j, independence gives

    E R^2=b,   E R^4=b+6 binom(b,2)=3b^2-2b<=3b^2.

Holder, applied to |R|^2=|R|^(2/3)|R|^(4/3), implies
E|R| >= (E R^2)^(3/2)/(E R^4)^(1/2) >= sqrt(b/3).
Thus L^2/4>=a^2 b/12. The elementary bound a^2 b>=k^3/8
is exact for even k; for k=2m+1, m>=1, its numerator difference is
8(m+1)^2 m-(2m+1)^3=4m^2+2m-1>0. Hence

    P N >= k^3/96.                                     (5)

In particular, on any k-vertex principal submatrix of A,

    min(max Q_(A[S]), max(-Q_(A[S]))) >= k^3/(96 F).      (6)

The stronger exact bound (4) is retained, but no constant optimization
or new signing census is needed for (1)--(3).

## 3. Symmetric perturbations inside the cube

Choose k coordinates with |z_i|<=1-epsilon. For any vector v supported
on these coordinates with entries +/-1 there, both z+epsilon v and
z-epsilon v lie in the cube. Their quadratic average is exactly

    [Q_A(z+epsilon v)+Q_A(z-epsilon v)]/2
       =Q_A(z)+epsilon^2 Q_A(v).                        (7)

Bounding this average above by F and maximizing Q_A(v) gives
F-Q_A(z)>=epsilon^2 max Q_(A[S]). Apply the same argument to -A
for F+Q_A(z). Equation (6) proves (1).

Let Delta=F-|Q_A(z)| and k(u)=#{i:1-|z_i|>=u}. If k(u)>=2,
(1) gives k(u)<=[96 F Delta]^(1/3) u^(-2/3); otherwise k(u)<=1.
Use the valid common bound

    k(u)<=1+[96 F Delta]^(1/3) u^(-2/3),  0<u<=1.

Integrating and using integral_0^1 u^(-2/3) du=3 proves (2).
Cubing its nonnegative part proves (3). The additive 1 is retained:
the argument cannot penalize an isolated fractional coordinate with
zero local field. No spectral estimate is used at any stage.

## 4. Extra term in the original mean-update inequality

Set M=A/sqrt(n), alpha=F/n^(3/2). For either phase sigma=+/-1,
let X_sigma have ANY law on the Boolean cube, and define

    Y_sigma=sign(sigma M X_sigma),
    k_sigma=#{i:Y_(sigma,i)!=X_(sigma,i)},
    z_sigma=(1-p)X_sigma+pY_sigma,   0<=p<=1/2.

Take sign(0) equal to the old coordinate, so zero-field ties are not
counted as updates. Set g(k)=k^3 for k>=2 and g(0)=g(1)=0. The k_sigma
changed coordinates of z_sigma have softness 2p. Equation (1), including
the p=0 case by its zero right side, gives the pointwise upper

    Q_(sigma M)(z_sigma)
       <= n alpha-p^2 g(k_sigma)/(24 alpha n^2).         (8)

The existing exact expansion is

    Q_(sigma M)(z_sigma)
      =(1-p)^2 Q_(sigma M)(X_sigma)
        +p(1-p)||M X_sigma||_1+p^2 Q_(sigma M)(Y_sigma).

Use Q_(sigma M)(Y_sigma)>=-n alpha. Average over the two phases
with equal weight (the laws need not be related), and put

    e=(1/(2n)) sum_sigma E Q_(sigma M)(X_sigma),
    f=(1/(2n)) sum_sigma E||M X_sigma||_1,
    eta=(1/(2n^3)) sum_sigma E g(k_sigma).

Then the strengthened finite inequality is

    (1+p^2)alpha >= (1-p)^2 e+p(1-p)f+p^2 eta/(24 alpha). (9)

The first two terms are reused from the September 17 paired-field note,
Section 5. The last term is the additional interior-gap correction.
It is valid for actual non-Gaussian post-update laws too, but does NOT
give values of their e, f, or eta, or authorize reapplying a Gaussian
source lemma to them. No new numerical lower-bound constant is claimed.

## 5. Scope and verification

The implication (9) is a same-order constraint on correction steps, not
the order-extension estimate needed for convergence. In particular,
the one-vertex formula still needs a single row that controls all
near-extremal Boolean states. Equations (1)--(3) do not select that row.
Nor do they lower-bound the signed pressure-defect integral.

The exact-object check covered the ten local branch heads, the four
worktree locations, the active source/update notes, and the September 6
analytic scratch member above. Its older cubic estimate is explicitly
retained as prior work, not presented as a new pairing route.

`boolean_interior_gap_20260923/check_identities.py` checks local algebra,
normalization, and boundary cases only. It does not formalize Holder,
the probabilistic selection argument, or this all-orders proof. Review
is by the author; no independent human or proof-assistant review is
claimed. The separate receipt records the mesh run and input hashes.
