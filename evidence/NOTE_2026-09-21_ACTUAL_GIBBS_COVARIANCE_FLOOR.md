# A covariance floor from actual Ising conditional scores

2026-09-21. Author-reviewed analytic lemma and actual-optimizer corollaries.
This is a LOWER bound on unsigned covariance, not an upper fluctuation
bound or a comparison between orders. The original MO limit remains OPEN.

## 1. Conditional-score lemma, including external fields

Let J be any real symmetric zero-diagonal matrix, f any real field, and
let mu(x) be proportional to exp(x^T J x/2 + f^T x) on the sign cube.
Write

    h_i = f_i + sum_j J_ij X_j,
    d_i(X) = sech(h_i)^2,   D = diag(E_mu d_i),
    C = Cov_mu(X),
    K(J) = 1 + ||J||op + 2 max_i sum_j J_ij^2.

Then, in positive-semidefinite order,

    C >= D^2 / K(J).                                         (1)

In particular, if E_mu |h_i| <= H for every i, then

    C >= exp(-4H) I / K(J).                                  (2)

There is no small-coupling or high-temperature assumption in (1).
The expectations, including those in D, are under the actual law with
the specified field; f is not discarded or replaced by its mean.

### Proof

Set R_i = X_i - tanh(h_i). Conditional expectation given all spins
except i gives E R_i=0 and

    Cov(X,R) = D,   E R_i^2 = E d_i.                         (3)

For i != j put t=J_ij and b_j=f_j+sum_(k != i,j) J_jk X_k.
As a function of X_i,

    tanh(b_j+t X_i) = u_j + v_j X_i,
    v_j = [tanh(b_j+t)-tanh(b_j-t)]/2.

Both u_j and v_j are measurable without X_i and X_j. Conditional
centering of R_i consequently gives the EXACT identity

    E R_i R_j = - E[v_j d_i].                               (4)

Indeed E R_i X_j=E R_i u_j=0, and conditioning out X_i in
E[v_j X_i R_i] replaces X_i R_i by d_i.

The function d(z)=sech(z)^2 is 2-Lipschitz. For either sign of t,

    v_j = (t/2) integral_(-1)^1 d(b_j+t s) ds,
    |v_j - t d(b_j+t X_i)| <= 2 t^2.                       (5)

Here integral_(-1)^1 |s-X_i| ds=2 for X_i in {-1,1}; thus
(5) is uniform in every field and spin configuration. Since d_i<=1,
(4) implies

    Cov(R) = D - J o T + E,
    T = E[d(X)d(X)^T],
    E_ii=0,   |E_ij| <= 2 J_ij^2.                          (6)

The circle denotes entrywise product. T is positive semidefinite and
T_ii<=1. Applying its positive Schur map to
-||J||op I <= J <= ||J||op I proves

    ||J o T||op <= ||J||op.

E is symmetric, so its maximum absolute row sum bounds its operator
norm. Together with 0<=D<=I, this proves Cov(R)<=K(J)I.
The covariance block matrix of (X,R) is positive semidefinite. More
explicitly, for any real vector z, choose a=Dz and apply Cauchy--Schwarz
to z^T(X-E X) and a^T R. Equations (3) and Cov(R)<=K(J)I give

    z^T C z >= ||Dz||_2^2 / K(J),

with the zero-vector case immediate. This is (1).
Finally sech(h)^2>=exp(-2|h|); Jensen gives
E d_i>=exp(-2 E|h_i|)>=exp(-2H), proving (2).

## 2. Both phases of a half-product pressure minimum

Let A be a complete order-N signing, beta=c/sqrt(N)>0, and suppose A
is edge-local minimal for

    a_A(beta) = [log Z_A(beta)+log Z_A(-beta)]/2,
    Z_A(beta) = E_uniform exp(beta Q_A).

Let U,V be its actual zero-field covariances at +beta and -beta.
Put

    K_A = 1 + beta ||A||op + 2 beta^2(N-1).

Then BOTH phases satisfy

    U,V >= eta_A I,
    eta_A = exp[-4(2c^2+1)] / K_A.                         (7)

Proof: write r_e^sigma=sigma A_ij E_sigma X_i X_j. The exact
edge-flip partition ratio in phase sigma is
cosh(2beta)-r_e^sigma sinh(2beta). Edge-local optimality of a_A,
followed by concavity of logarithm, implies

    (r_e^+ + r_e^-)/2 <= tanh(beta).                       (8)

For a fixed row put g_i=beta sum_(j != i) A_ij X_j.
In either phase, integrating X_i gives

    E_sigma[g_i tanh(g_i)]
       = beta sum_(j != i) r_ij^sigma >= 0.

The average over the two phases is at most
(N-1)beta tanh(beta)<=c^2 by (8). Each nonnegative phase quantity
is therefore at most 2c^2. The scalar bound
|g|<=g tanh(g)+1 proves E_sigma |g_i|<=2c^2+1: for u>=0,
u(1-tanh(u))=2u/(exp(2u)+1)<=1.
Apply (2) with J=sigma beta A and f=0 to obtain (7).

## 3. A symmetric-pressure minimum also has an unsigned floor

Instead suppose A is edge-local minimal for

    F_A(beta)=log E_uniform cosh(beta Q_A).

Use its ACTUAL augmented Gibbs law on (sigma,X), proportional to
exp(sigma beta Q_A). Its phase probabilities need not be equal.
Let C_sym=E[XX^T], which is unsigned; E X=0 by spin reversal. Then

    C_sym >= exp[-4(c^2+1)] I / K_A.                       (9)

To check the mixture step, define R_i=X_i-tanh(sigma g_i).
Equations (3)--(6) remain valid after conditioning on sigma, except
that the Schur matrix is the signed matrix

    T_signed=E[sigma d(X)d(X)^T].

It is the difference of the two phase positive Schur maps, with
nonnegative weights summing to one. Each is a contraction in the
operator estimate used above, so ||beta A o T_signed||op
is still at most beta||A||op. The remainder bound is unchanged.
Thus (1) holds for C_sym and D=diag(E d_i) under this mixture.

The exact edge-local inequality from the existing variational-control
note gives r_e=A_ij E[sigma X_i X_j]<=tanh(beta). Integrating X_i
again gives E[g_i tanh(g_i)]<=c^2 for each row, hence
E|g_i|<=c^2+1. Jensen now proves (9). This argument does NOT assert
a positive-semidefinite lower bound on the signed matrix E[sigma XX^T].

## 4. Actual single-row cavity fields: a finite-field consequence

Assume the stronger hypothesis that A is a GLOBAL half-product minimum.
Fix i, let T omit i, and consider the two actual cavity laws

    nu_(sigma,t)(x_T) proportional to
      exp(sigma beta Q_(A_T)(x_T) + t sigma g_i(x_T)),
    -1<=t<=1.

Their field direction is the actual removed row, not an arbitrary
subsequently chosen vector. Set

    H_c = 2 exp(c^2)(2c^2+1+c)+c,
    eta_cavity = exp(-4 H_c) / K_A.

Uniformly in sigma and t,

    Cov_(nu_(sigma,t))(X_T) >= eta_cavity I.                (10)

Here is the needed change-of-measure proof. Whole-row random replacement
and global minimality give

    0 <= p_i^sigma := log E_(nu_(sigma,0)) cosh(g_i),
    (p_i^+ + p_i^-)/2 <= (N-1)log cosh(beta) <= c^2/2.

Thus p_i^sigma<=c^2 for each phase. The full phase marginal on T is
nu_(sigma,0) tilted by exp(-p_i^sigma)cosh(g_i). Spin reversal gives
E_(nu_(sigma,0)) exp(t sigma g_i)=E cosh(t g_i)>=1, while
exp(t sigma g_i)<=2cosh(g_i). Consequently

    d nu_(sigma,t) / d(mu_sigma)_T <= 2 exp(c^2).          (11)

For j in T, write g_j^T=beta sum_(k in T, k != j) A_jk X_k.
Under the full phase law, E|g_j^T|<=2c^2+1+beta by Section 2.
The local field in nu_(sigma,t) is sigma g_j^T+t sigma beta A_ji.
Equation (11) bounds its expected absolute value by H_c, using beta<=c.
The principal interaction matrix sigma beta A_T has operator norm at
most beta||A||op and squared row sums at most beta^2(N-1).
Apply (2) to obtain (10).

The point of (10) is a simultaneous matrix bound for EVERY test vector
under these actual nonzero-field laws. A scalar lower response in the
removed-row direction was already supplied by the signed-field theorem
in NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md, with the field amplitude
set to beta. That scalar consequence is not a new result here.

## 5. Scope, reuse, and remaining implication

Reused results: the exact edge-flip identity and whole-row replacement
comparison from NOTE_2026-09-05_GLOBAL_OPTIMIZER_VARIATIONAL_CONTROL.md,
and the exact half-product deletion comparison from
NOTE_2026-09-05_EXACT_HALFPRODUCT_SUBCRITICAL_SPECTRAL.md.
The new step is the conditional-score covariance bound (1), followed
by the actual-law first-moment and domination arguments above.

At fixed c and a fixed normalized operator cap, (7), (9), and (10)
give dimension-independent positive floors. Without such a cap their
printed denominators must be retained. Same-order regularization alone
does NOT transfer edge-local or global pressure optimality to its output.
The earlier exact-half-product theorem supplies only
||A_N||op=o(N^(3/4)), not a fixed normalized cap.

These are lower susceptibility bounds. They do not bound fourth moments
from above, establish uniform integrability, control fields from a growing
set of pinned rows, or lower-bound the signed internal/cross pressure-defect
integral. That last order-comparison implication remains unresolved.

The two-spin conditional identities have a separate exact symbolic check
in actual_gibbs_covariance_floor_20260921/. It is an algebra regression,
not a signing census or formal verification of the all-orders proof.
