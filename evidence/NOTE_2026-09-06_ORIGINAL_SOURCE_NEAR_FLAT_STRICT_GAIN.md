# Strict original-source gain in the bounded near-flat spectral regime

2026-09-06. Analytic theorem for ACTUAL complete symmetric signings.
No mathematical computation, solver, checker, spectral scan, numerical
integral, optimization run, or construction search was executed.

The argument improves an actual positive Gaussian phase by a fixed,
independent-coordinate Boolean update. Every objective below is the
ORIGINAL same-source quadratic energy. No bilinear norm is substituted,
and no polarization or change-of-order factor is suppressed.

## 1. Statement and scope

For a real complete symmetric zero-diagonal signing A of order n, put

    Q_A(x)=x^T A x/2,
    Phi(A)=max_(x in {+1,-1}^n)|Q_A(x)|,
    kappa=2/pi,
    rho=16/25,
    lambda=(5/4)sqrt(n)=sqrt(n/rho).

Let n tend to infinity through any sequence of such actual A. Assume

    limsup ||A||op/sqrt(n)<=5/3,                               (1.1)
    empirical eigenvalue law of A/sqrt(n)
        -> (9/25)delta_0+(8/25)(delta_(5/4)+delta_(-5/4)).      (1.2)

Then

    liminf Phi(A)/n^(3/2)
        >=5kappa/8+16/3125
        >2/5+3/1100.                                         (1.3)

The conclusion concerns a near-flat limiting law, not exact finite-n
flatness. In particular (1.2) does not impose a large exact kernel.
The operator cap in (1.1) is retained in the update penalty; weak
empirical convergence does NOT allow replacing it by 5/4.

The theorem is conditional on this specific actual source spectral
regime. It does not prove that all candidate extremizers enter the
regime, exclude every formal trace profile, close every active cell,
or solve the global original MO target.

## 2. Positive spectral projector and a genuine correlation matrix

Let

    P_+=1_(A>lambda/2),
    P_-=1_(A<-lambda/2),
    P=P_++P_-,
    T=P_+-P_-.

The fixed cutoffs lie away from the three limiting spectral atoms.
The spectra of A/lambda are uniformly bounded by (1.1). Therefore
weak convergence in (1.2), applied on this common compact interval,
gives

    rank(P_+)/n ->rho/2,       rank(P_-)/n ->rho/2,
    ||A/lambda-T||_F=o(sqrt(n)),
    ||A^2/lambda^2-P||_F=o(sqrt(n)).                          (2.1)

For example, each squared Frobenius error divided by n is the
empirical integral of the corresponding bounded, piecewise continuous
spectral error squared. Its only discontinuities are the two cutoffs,
neither of which carries limiting mass; it vanishes at 0 and +/-1.

Completeness and zero diagonal give EXACTLY

    A_ii=0,       (A^2)_ii=n-1,       lambda^2 rho=n.

Taking diagonal parts in (2.1) consequently gives

    (1/n)sum_i T_ii^2 ->0,
    (1/n)sum_i (P_ii-rho)^2 ->0.

Define the positive semidefinite matrix

    R0=(2/rho)P_+,       r_i=(R0)_ii.

Since 2P_+=P+T, the preceding diagonal estimates imply

    ||R0||op=2/rho when P_+ is nonzero,
    (1/n)sum_i(r_i-1)^2 ->0.                                (2.2)

The raw diagonal r_i need not be constant or bounded below. We therefore
normalize only good coordinates and fill bad coordinates independently.
Choose eta=eta_n tending to zero slowly enough that

    |I^c|=o(n),       I={i:|r_i-1|<=eta}.

For instance, if h_n=(1/n)sum_i(r_i-1)^2, the choice
eta=max(n^(-1/8),h_n^(1/4)) has this property. Monotone decrease of
eta is unnecessary; only eta->0 and h_n/eta^2->0 are used.
For all sufficiently large n, eta<1/2. On I define

    R[I,I]=diag(r_i^(-1/2):i in I) R0[I,I]
                                      diag(r_i^(-1/2):i in I).

Set R[I^c,I^c]=I and R[I,I^c]=R[I^c,I]=0. Then R is an ACTUAL
correlation matrix, possibly singular, and

    R>=0,       diag R=1,
    ||R||op<=max(2/[rho(1-eta)],1)=2/rho+o(1),
    ||R-R0||_F=o(sqrt(n)).                                  (2.3)

Here is an explicit justification of the last estimate. Removing the
rows and columns in I^c from R0 costs O(sqrt(|I^c|)) in Frobenius norm,
because R0 has bounded operator norm and those removed matrices have
rank at most 2|I^c|. Inserting I on I^c has the same size. The good
diagonal normalization differs from the identity by O(eta) in operator
norm, and ||R0||_F=O(sqrt(n)), so it costs O(eta sqrt(n)). All three
costs are o(sqrt(n)).

Let G be a centered Gaussian vector with covariance R, and put

    X_i=sign(G_i).

Each G_i is a standard Gaussian, so its zero event has probability
zero. The convention sign(0)=+1 will also be used for local fields.
All bad coordinates remain genuine Boolean variables; none is deleted
from A or from the original energy.

## 3. Uniform Schur powers and the original phase baseline

Fix a constant C_R such that ||R||op<=C_R for all sufficiently large
n. For every integer q>=1,

    R^{circ q}>=0,       diag(R^{circ q})=1,
    ||R^{circ q}||op<=||R||op<=C_R.                           (3.1)

Indeed Schur multiplication by a correlation matrix preserves positive
semidefinite order and sends cI to cI. Iteration proves (3.1), including
when some entries of R are negative. In particular

    ||R^{circ q}||_F<=C_R sqrt(n),       tr R^2<=C_R n.       (3.2)

Write the normalized Gaussian Hermite expansion

    sign(z)=sum_(q odd>=1)c_q h_q(z),
    c_1=sqrt(kappa),       sum_(q odd>=1)c_q^2=1.

Orthogonality and the correlated Hermite identity give

    C_X=E[XX^T]=kappa R+C_tail,
    C_tail=sum_(q odd>=3)c_q^2 R^{circ q},
    0<=C_tail<=(1-kappa)||R||op I.                            (3.3)

The series converges in operator norm, since its tails are bounded by
C_R times the corresponding scalar coefficient tails.

The actual complete entries and zero diagonal imply, uniformly for
all odd q>=3,

    |tr(A R^{circ q})|
        <=sum_(i!=j)|R_ij|^q
        <=sum_(i,j)R_ij^2
        =tr R^2<=C_R n.                                     (3.4)

Consequently the original positive phase has expectation

    E Q_A(X)=(kappa/2)tr(AR)+O(n).                           (3.5)

By (2.3), ||A||_F=sqrt(n(n-1)), and Cauchy--Schwarz,

    tr[A(R-R0)]=o(n^(3/2)).

Also (2.1) implies

    tr(AR0)=(2/rho)tr(AP_+)=n lambda+o(n^(3/2)).

Thus

    E Q_A(X)/n^(3/2) ->5kappa/8.                             (3.6)

This is an expected value of Q_A itself. We will improve this actual
positive energy, not a cross norm or a surrogate norm.

## 4. The higher-chaos variance has the stronger mean lower bound

The projector relations P+T=2P_+=rho R0 and (2.1)--(2.3) imply

    E_n=A^2/lambda^2-rho R+A/lambda,
    ||E_n||_F=o(sqrt(n)).                                   (4.1)

For every odd q>=3, take its trace pairing with R^{circ q}:

    tr(A^2 R^{circ q})/lambda^2
      =rho tr(R R^{circ q})
         -(1/lambda)tr(A R^{circ q})+tr(E_n R^{circ q}).

The first trace is

    tr(R R^{circ q})=sum_(i,j)R_ij^(q+1)>=n,                 (4.2)

because q+1 is even and all n diagonal terms are one. Equations
(3.2), (3.4), and (4.1) bound the other two terms uniformly in q.
It follows that

    tr(A^2 R^{circ q})/lambda^2>=rho n-o(n)
                              uniformly for odd q>=3.       (4.3)

Define each normalized actual local field and its first-chaos part by

    F_i=(AX)_i/sqrt(n),
    F_(i,1)=sqrt(kappa)(AG)_i/sqrt(n),
    v_i=E(F_i-F_(i,1))^2
        =(1/n)(A C_tail A)_ii>=0,
    mu_n=(1/n)sum_i v_i.

The coefficient tail in (3.3) has total mass 1-kappa. Since
lambda^2 rho=n, summing (4.3) with those nonnegative coefficients
gives the STRONGER lower bound

    mu_n=(1/n^2)tr(A^2 C_tail)>=1-kappa-o(1).                (4.4)

No interchange of a nonuniform error and an infinite series occurs:
the error in (4.3) was uniform for every odd q>=3. On the other hand,
each complete row of A has squared length n-1, so (3.3) gives

    0<=v_i<=(1-kappa)||R||op (n-1)/n<=V+o(1)
                    uniformly in i,
    V=2(1-kappa)/rho.                                       (4.5)

The mean is not divided by two. Identity (4.1), which retains the
A/lambda term and uses the actual complete entries in (3.4), is what
supplies the full 1-kappa lower bound. The v_i need not be homogeneous.

## 5. First-chaos alignment in average squared mean

Set

    a=sqrt(kappa/rho)=(5/4)sqrt(kappa).

We claim

    (1/n)sum_i E|F_(i,1)-a G_i|^2 ->0.                       (5.1)

The left side is exactly

    (kappa/n^2)tr[(A-lambda I)^2 R].

For R0 the trace equals

    (2/rho)||(A-lambda I)P_+||_F^2=o(n^2),

by the bounded spectrum and the flat positive limiting atom in (2.1).
For the difference, the ACTUAL operator cap in (1.1) implies

    ||(A-lambda I)^2||_F=O(n^(3/2)),
    |tr[(A-lambda I)^2(R-R0)]|=o(n^2)

using (2.3). This proves (5.1) without inverting R and without
inferring variance control merely from a Frobenius estimate for AR.

Write

    b_i=E F_(i,1)^2,
    c_i=E[G_i F_i]=E[G_i F_(i,1)],
    sigma_i^2=E F_i^2=b_i+v_i.                              (5.2)

The equality for c_i uses orthogonality of the higher chaoses to the
first chaos. The row lengths and the bounded operator norm of R give
uniform bounds on b_i, sigma_i^2, and |c_i|. From (5.1) and
Cauchy--Schwarz it follows that

    (1/n)sum_i |c_i-a| ->0,
    (1/n)sum_i |b_i-a^2| ->0.                               (5.3)

For the second statement, factor b_i-a^2 as the expectation of
(F_(i,1)-aG_i)(F_(i,1)+aG_i) and average; the second factor has a
uniformly bounded second moment. In particular all but o(n) rows have
c_i bounded away from zero, and b_i-a^2 and c_i-a tend to zero in
empirical probability.

## 6. Uniform local Gaussianization, including the distinguished sign

We use the following extension of the delocalized Gaussian-sign lemma
in `original_mo_complete_cross_flat_spectral_gain.md`, Sections 2--5.
Its additional point is a JOINT limit with one chosen input coordinate,
not a growing-dimensional joint limit of all local fields.

### Lemma

Fix C>=1. Let R be any correlation matrix of order n with ||R||op<=C,
let G be centered Gaussian with covariance R, and let

    F=sum_j d_j sign(G_j),
    max_j |d_j|<=1/sqrt(n),       sum_j d_j^2<=1.

Choose any index i. Set sigma^2=E F^2 and c=E[G_i F]. Uniformly over
these data,

    E|F|=sqrt(kappa) sigma+o_C(1).                           (6.1)

For every fixed c_*>0, uniformly over the data also satisfying c>=c_*,

    P(sign(F)!=sign(G_i))
       =(1/pi) arccos(c/sigma)+o_(C,c_*)(1).                 (6.2)

Here sign(0)=+1. The covariance inequality gives sigma>=c_*, so the
ratio is well defined and in [-1,1]. Degenerate Gaussian limits with
correlation one are permitted. All errors tend to zero as n increases.

### Proof of the extension

For every q>=1 the Schur power has operator norm at most C, by (3.1).
Let u_j be unit vectors with Gram matrix R. In the Gaussian polynomial
normalization I_q(u^{tensor q})=H_q(<Z,u>), define

    g_q=sum_j d_j u_j^{tensor q},
    f_q=(c_q/sqrt(q!))g_q,
    F_Q=sum_(q odd<=Q) I_q(f_q).

The proof of the cited scalar lemma uses only the two displayed bounds
on d, not that every coefficient is nonzero or has maximal magnitude.
In detail, for 1<=r<min(p,q), set

    M=diag(d)R^{circ r}diag(d),
    N_p=R^{circ(p-r)},       N_q=R^{circ(q-r)}.

Then

    ||g_p tensor_r g_q||^2=tr(M N_q M N_p)<=C^4/n.            (6.3)

If r=p<q, write z=R^{circ p}d; then

    ||g_p tensor_p g_q||^2
       =z^T diag(d)R^{circ(q-p)}diag(d)z<=C^3/n.             (6.4)

The case r=q<p follows by symmetry. Full equal-order contractions
are scalar covariance terms, not fluctuation terms. These estimates
follow from ||diag(d)||op<=1/sqrt(n), ||d||<=1, and (3.1).

For a fixed linear combination s G_i+t F_Q, the first-order kernel is
s u_i+t sqrt(kappa)g_1; its higher-order kernels are t f_q. The NEW
mixed contraction, for q>=3, is bounded by

    ||u_i tensor_1 g_q||^2
      =(d circ R_i)^T R^{circ(q-1)}(d circ R_i)
      <=(C/n)sum_j R_ij^2<=C^2/n.                           (6.5)

In the last step R^2<=C R gives sum_j R_ij^2<=C. Together with
(6.3)--(6.4), this controls every nonconstant contraction for the
fixed linear combination, including its first/higher-chaos interaction.

For completeness, write L_Q=sG_i+tF_Q as a finite sum sum_q I_q(k_q),
and U_Q=sum_q I_q(k_q)/q. Gaussian integration by parts and the
Ornstein--Uhlenbeck eigenvalue identity give

    Gamma_Q=grad L_Q dot grad U_Q,
    E[L_Q h(L_Q)]=E[h'(L_Q)Gamma_Q],
    E Gamma_Q=E L_Q^2=omega_Q^2.

The polynomial product identity

    I_p(f)I_q(g)=sum_(r=0)^min(p,q)
       r! binom(p,r) binom(q,r)
          I_(p+q-2r)(sym(f tensor_r g))

expands Gamma_Q. Its full equal-order contractions are exactly the
constant omega_Q^2. Every other term is bounded by (6.3)--(6.5).
Gaussian polynomial isometry and the finite number of terms yield

    Var Gamma_Q<=K_(Q,C,s,t)/n.                              (6.6)

These identities follow directly from the Hermite generating function
and finite-dimensional Gaussian integration by parts, as detailed in
the cited scalar lemma. They also apply to singular R through its
Gram representation; no nonsingularity or inverse variance is needed.

If phi_Q(z)=E exp(izL_Q), the exact differential identity is

    phi_Q'(z)+omega_Q^2 z phi_Q(z)
       =-z E[(Gamma_Q-omega_Q^2)exp(izL_Q)].

Solving with phi_Q(0)=1 and using (6.6) gives, for every real z,

    |phi_Q(z)-exp(-omega_Q^2 z^2/2)|
                 <=(z^2/2)sqrt(K_(Q,C,s,t)/n).               (6.7)

Let tau_Q=sum_(q>Q)c_q^2. Orthogonality and the Schur-power cap give

    E(F-F_Q)^2<=C tau_Q,
    E F^2<=C.                                              (6.8)

The omitted tail is also orthogonal to G_i and all retained chaoses.
Thus the characteristic function of sG_i+tF differs from that of
L_Q by at most |zt|sqrt(C tau_Q), and their variances differ by at
most t^2 C tau_Q. Taking n to infinity and then Q to infinity in
(6.7) proves a Gaussian characteristic-function approximation for
EVERY fixed linear combination sG_i+tF, with its own actual variance.

The pair covariances lie in a compact set: Var G_i=1, Var F<=C,
and |c|<=sqrt(C). On every covariance-convergent subsequence the
Cramer--Wold criterion now gives convergence of (G_i,F) to the centered
Gaussian pair with that limiting covariance. The uniform second-moment
bound makes |F| uniformly integrable and proves (6.1). Uniformity
follows by contradiction: any violating sequence has such a compact
covariance-convergent subsequence.

If c>=c_*, both limiting marginals have positive variance. The boundary
of the sign-disagreement set lies in the two coordinate axes, each of
zero limiting Gaussian probability, even when the pair is degenerate.
Convergence of probabilities therefore applies, and the Gaussian-pair
disagreement probability is arccos(c/sigma)/pi. The same compactness
contradiction proves (6.2) uniformly. The Gaussian formula follows by
representing the pair as two linear forms in two independent standard
Gaussians and measuring the angular sectors where their signs differ;
the degenerate endpoint follows by continuity. This completes the lemma.

## 7. Local gain and sign-disagreement frequency

Apply the lemma to each actual row of A with d_j=A_ij/sqrt(n), including
its one zero diagonal coefficient. These coefficients meet both bounds
in Section 6. Define

    G_n=E[||AX||_1-X^T A X]/n^(3/2),
    p_n=(1/n)E|{i:sign((AX)_i)!=X_i}|.

By (6.1), (5.2)--(5.3), and uniform continuity of square root on a
fixed compact interval,

    (1/n)sum_i E|F_i|
       =sqrt(kappa)(1/n)sum_i sqrt(a^2+v_i)+o(1).

Equation (3.6) gives E[X^TAX]/n^(3/2)=5kappa/4+o(1), and
sqrt(kappa)a=5kappa/4. Hence

    G_n=sqrt(kappa)(1/n)sum_i[sqrt(a^2+v_i)-a]+o(1).         (7.1)

All but o(n) rows have c_i>=a/2 by (5.3). On these rows (6.2) applies
uniformly. The other rows contribute at most o(1) to p_n. Equations
(5.3) and (4.5), followed by uniform continuity of the Gaussian angle
on the compact region with sigma bounded away from zero, give

    p_n=(1/pi)(1/n)sum_i arctan(sqrt(v_i)/a)+o(1).           (7.2)

Indeed for b=a^2 and c=a,
arccos(c/sqrt(b+v))=arctan(sqrt(v)/a), including at v=0. The
empirical mean errors in (5.3) tend to zero, so one can first restrict
to rows with both errors at most any fixed tolerance, then let that
tolerance go to zero. This handles the possibly non-Lipschitz endpoint
without assuming a uniform lower bound on v_i.

No joint limit of n fields was used. Nor was an operator bound on the
covariance of sign(AX) asserted or required.

## 8. An actual independent-coordinate Boolean update

For each realization of X, let

    Y=sign(AX),       Delta=Y-X.

Conditionally on X, choose independent Bernoulli(epsilon) variables
xi_i and set X'_i=X_i+xi_i Delta_i. For every 0<=epsilon<=1, X'
is an actual Boolean vector on the SAME n coordinates. Since A_ii=0,
expansion of its original quadratic form gives exactly

    E_xi Q_A(X')
       =Q_A(X)+epsilon(AX)^T Delta
                    +(epsilon^2/2)Delta^T A Delta
       =Q_A(X)+epsilon[||AX||_1-X^T A X]
                    +(epsilon^2/2)Delta^T A Delta.           (8.1)

Independence is used only for distinct-coordinate terms; the diagonal
terms vanish exactly. At a zero local field the chosen sign convention
still gives (AX)_i Y_i=|(AX)_i|.

Each disagreement has Delta_i^2=4 and other coordinates have Delta_i=0.
Thus

    Delta^T A Delta>=-||A||op ||Delta||^2
       =-4||A||op |{i:Y_i!=X_i}|.

Every realization satisfies Q_A(X')<=Phi(A). Averaging (8.1), with
C_n=||A||op/sqrt(n), yields the actual-source lower bound

    Phi(A)/n^(3/2)
       >=E Q_A(X)/n^(3/2)+epsilon G_n-2epsilon^2 C_n p_n.    (8.2)

The penalty keeps C_n<=5/3+o(1) from (1.1). The nonzero limiting
atom 5/4 cannot replace that actual cap: a vanishing fraction of
outlying eigenvalues could survive (1.2).

## 9. A fixed admissible probability gives a strict rational margin

The following scalar estimates are the fixed-probability argument in
`original_mo_original_source_local_update_scalar_gain.md`, reproduced
here with the vanishing errors justified by (4.4)--(4.5).
Only the already established interval

    7/11<kappa<16/25                                         (9.1)

is used. It follows from the frozen pi enclosure recorded in
`original_mo_source_cross_nuclear_trace_boundary.md`; no certificate
was rerun. In particular its upper endpoint 31415927/10000000 is
less than 22/7 because 7*31415927=219911489<220000000.

On [0,V], concavity and the zero value at v=0 give the chord bound

    sqrt(kappa)[sqrt(a^2+v)-a]>=c0 v,
    c0=sqrt(kappa)/(sqrt(a^2+V)+a)
       =sqrt(kappa rho)/(sqrt(kappa)+sqrt(2-kappa))>8/25.     (9.2)

The strict comparison is equivalent to 13kappa>8, which follows
from (9.1). The uniform enlargement v_i<=V+o(1) changes this bound
by only o(1) after empirical averaging. Thus (7.1) gives

    G_n>=(8/25)mu_n-o(1).                                   (9.3)

The elementary inequalities arctan x<=x and
(1/n)sum_i sqrt(v_i)<=sqrt(mu_n) give from (7.2)

    p_n<=sqrt(mu_n)/(pi a)+o(1)
        <=(8/25)sqrt(mu_n)+o(1),
    pi a=5/(2sqrt(kappa))>25/8.                             (9.4)

Choose the FIXED probability epsilon=1/10. It is admissible regardless
of the variance profile. Equations (8.2), (9.3)--(9.4), and
C_n<=5/3+o(1) give

    epsilon G_n-2epsilon^2 C_n p_n
       >=(4/125)mu_n-(4/375)sqrt(mu_n)-o(1)
        =(4/375)(3r_n^2-r_n)-o(1),       r_n=sqrt(mu_n).

By (4.4) and kappa<16/25, liminf r_n>=sqrt(1-kappa)>3/5. The
polynomial 3r^2-r increases for r>=3/5. Consequently

    liminf[epsilon G_n-2epsilon^2 C_n p_n]
       >=(4/375)[3(3/5)^2-3/5]=16/3125>1/200.                (9.5)

Together with (3.6), this proves the first bound in (1.3). For its
strict comparison, kappa>7/11 gives 5kappa/8>35/88, and

    16/3125>1/200,
    2/5-35/88=1/440,
    1/200-1/440=3/1100.

Thus 5kappa/8+16/3125>2/5+3/1100, proving (1.3). Neither a numerical
integration nor an unconstrained quadratic optimizer is used. The
heterogeneous empirical distribution of v_i need not converge.

## 10. Original near-scalar paired-source corollary

The separate supporting lemma
`original_mo_near_scalar_internal_flat_law_transfer.md` states the
following. For an ACTUAL paired complete signing

    K=[[A,B],[B^T,-A]],
    D=diag(D_L,D_R)>0 diagonal,       D+/-K>=0,
    dbar=tr(D)/(2n),
    delta=tr(D)tr(D^(-1))/(2n)^2-1,
    H_L=D_L^(-1/2) A D_L^(-1/2),

if delta->0 and the empirical law of the FULL actual H_L tends to

    (9/25)delta_0+(8/25)(delta_(3/4)+delta_(-3/4)),

then one common original principal source A_J, of order q with q/n->1,
satisfies (1.1)--(1.2) with q in place of n and

    Phi(A_J)/q^(3/2)<=Phi(A)/n^(3/2)+o(1).

Applying (1.3) to that actual source proves the same strict lower for
liminf Phi(A)/n^(3/2). In particular this excludes the actual
near-scalar internal-law regime with Phi(A)/n^(3/2)->2/5 underlying
the specific strengthened FORMAL profile. This is an original-source
contradiction; the paired covariance, cross block, and active field
are not replaced or reoptimized. No separate claim of optimality of
D is needed, and the transfer lemma does not need a separate trace cap.

This corollary uses the independently written transfer lemma as a
separate premise. It does NOT assert that an arbitrary active sequence
has delta->0 or this internal law. Other actual profiles, the remaining
active-state implication, and the global original MO target stay OPEN.

## 11. Dependencies, contribution scope, and publication boundary

The frozen predecessor Gaussian polynomial argument is
`original_mo_complete_cross_flat_spectral_gain.md` (411 lines), SHA256
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
Section 6 extends that argument by allowing zero/smaller coefficients
and proving the distinguished-coordinate contractions and joint limit.

The authored scalar support is
`original_mo_original_source_local_update_scalar_gain.md` (209 lines),
SHA256 `7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`.
The separate authored original-source transfer is
`original_mo_near_scalar_internal_flat_law_transfer.md` (141 lines),
SHA256 `f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8`.
The original-phase normalization context is
`original_mo_original_phase_spectral_moment.md` (262 lines), SHA256
`7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600`.
The reused pi-enclosure source is
`original_mo_source_cross_nuclear_trace_boundary.md`, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

Root supplied the positive-projector/local-update strategy, the
relevant actual operator cap, and the fixed-probability specialization,
and independently emphasized the trace-of-square first-chaos estimate.
The proof worker derived the robust correlation normalization, the
stronger higher-chaos mean via (4.1), the distinguished-coordinate
Gaussianization extension, and this full actual-source proof. The
docs-gate worker authored the distribution-free scalar support with
the operator cap and fixed probabilities retained. The exact-saturation
worker authored the separate original-source transfer lemma. These
roles are contributions, not claims of independence for those same
derivations. A full-source reviewer must disclose any such contribution.

At the time of authorship this is a /tmp proof artifact only. No
canonical repository file was edited, publication was performed, or
mathematical job was executed by this artifact's author.
