# Actual complete-cross gain beyond the cubic spectral floor

2026-09-05. Analytic proof. No mathematical computation, solver,
signing search, numerical integration, or spectral scan was run.

The theorem concerns ACTUAL complete cross sign matrices. Its gain uses
their delocalized entries, not only their singular-value measure. The
logical core is a self-contained Gaussian polynomial argument; the
standard multiple-chaos limit theorems were independently checked as
backup, not assumed in the proof below.

## 1. Statement, normalization, and exact scope

Write kappa=2/pi. For an actual n by n matrix B with every entry in
{+1,-1}, put

    beta(B)=max_(x,y in {+1,-1}^n) x^T B y
           =max_x ||B^T x||_1.

Choose a real d>=||B||op and define

    m=n/d^2,
    epsilon=1-tr[(B^T B)^2]/(n^2 d^2).

Then 0<m<=1 and 0<=epsilon<=1-m. For every fixed m_0>0 there is a
function e_(m_0)(n) tending to zero such that, uniformly over all these
actual B and d with m>=m_0,

    beta(B)/(nd)
      >=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-e_(m_0)(n).  (1.1)

In particular, if all nonzero singular values of B equal d, then
epsilon=0 and

    beta(B)/(nd)
      >=kappa+(sqrt(kappa)-kappa)m-o_(m_0)(1).              (1.2)

Thus exact flatness with m bounded away from zero forces a strictly
positive leading-order gain over kappa. Formula (1.1) also retains a
quantified leading-order penalty for failure of flatness. No finite-n
rate for e_(m_0) is claimed.

These are actual unweighted cross-norm statements. They do not assume
or produce an internal matrix A, a pure-cross active optimizer, or any
conditional-optimizer property. The hypothesis d>=||B||op together
with m>=m_0 is an actual O(sqrt(n)) operator bound. It is NOT silently
deduced from a near-scalar diagonal majorizer or its trace alone.

## 2. A uniform absolute-moment lemma for delocalized Gaussian signs

Fix C>=1. Let R be any n by n real correlation matrix with

                       0<=R<=C I,       diag R=1.

Let G be a centered Gaussian with covariance R, possibly singular,
and let a_i in {+1/sqrt(n),-1/sqrt(n)}. Set

    F=sum_i a_i sign(G_i),       sigma^2=E F^2.

Then, uniformly over R and all such coefficient vectors a,

                |E|F|-sqrt(kappa) sigma|<=e_C(n),
                e_C(n) tends to zero as n tends to infinity.   (2.1)

This is a scalar marginal assertion, uniform over the input data.
It does NOT assert joint convergence of a growing vector of n column
sums. Such a growing-dimensional limit is unnecessary here.

For every integer q>=1 the entrywise power R^{circ q} is a
correlation matrix and satisfies

                       0<=R^{circ q}<=C I.                       (2.2)

Indeed, taking a Schur product with a correlation matrix preserves
positive semidefinite order and sends C I to C I. Induction starting
from R proves the assertion, including for negative entries of R.

Let h_q be the normalized probabilists' Hermite polynomials. In the
one-dimensional Gaussian L2 space write

    sign(z)=sum_(q odd>=1) c_q h_q(z),
    sum_q c_q^2=1,       c_1=sqrt(kappa).

The first coefficient is E|Z| for a standard Gaussian Z. The even
coefficients and the constant coefficient vanish by parity.
For every integer Q>=1 put

    F_Q=sum_(q odd<=Q) c_q sum_i a_i h_q(G_i),
    tau_Q=sum_(q>Q) c_q^2.

Orthogonality of distinct chaoses, the correlated Hermite identity,
and (2.2) give

    E F^2<=C,       E F_Q^2<=C,
    E(F-F_Q)^2
        =sum_(q>Q) c_q^2 a^T R^{circ q} a<=C tau_Q.        (2.3)

In particular tau_Q tends to zero, independently of n and R.
Sections 3--5 prove the finite-Q Gaussian approximation needed to
turn this L2 tail estimate into (2.1).

## 3. Tensor kernels and all mixed contraction estimates

Choose unit vectors v_1,...,v_n with Gram matrix R and a standard
Gaussian vector Z in their finite-dimensional span, so G_i=<Z,v_i>.
For a symmetric tensor f of order q, use the Gaussian polynomial
normalization

    I_q(v^{tensor q})=H_q(<Z,v>)       when ||v||=1,
    E[I_q(f)I_q(g)]=q! <f,g>.

Here H_q=sqrt(q!)h_q. Define the unscaled kernels

    g_q=sum_i a_i v_i^{tensor q},
    f_q=(c_q/sqrt(q!))g_q,

so F_Q=sum_(q odd<=Q) I_q(f_q). Only finitely many orders occur
for any fixed Q. The coefficients in passing between g_q and f_q
depend only on Q, not on the dimension or the correlation matrix.

For 1<=r<=min(p,q), contraction of r tensor indices gives

    g_p tensor_r g_q
      =sum_(i,j) a_i a_j R_ij^r
                    v_i^{tensor(p-r)} tensor v_j^{tensor(q-r)}.

Suppose first r<min(p,q). Write

    M=diag(a) R^{circ r} diag(a),
    N_p=R^{circ(p-r)},       N_q=R^{circ(q-r)}.

Expanding the squared tensor norm exactly yields

    ||g_p tensor_r g_q||^2=tr(M N_q M N_p)
                =||N_q^(1/2) M N_p^(1/2)||_F^2.

All the entrywise powers here have positive exponents, so (2.2)
and ||diag(a)||op=1/sqrt(n) imply

    ||M||_F^2<=||R^{circ r}||_F^2/n^2<=C^2/n,
    ||g_p tensor_r g_q||^2<=C^4/n.                       (3.1)

Next suppose r=p<q. Collapse the fully contracted p indices and
write b=R^{circ p}a. Then

    g_p tensor_p g_q=sum_j a_j b_j v_j^{tensor(q-p)},
    ||g_p tensor_p g_q||^2
      =b^T diag(a) R^{circ(q-p)} diag(a)b
      <=(C/n)||b||^2<=C^3/n.                            (3.2)

The case r=q<p follows by reversing the tensor factors. This also
covers the interaction of the first chaos with every higher chaos;
it is not dismissed merely because the first chaos is Gaussian.

The only remaining case is r=p=q. That contraction is the scalar
<g_p,g_p>. In the argument below it contributes exactly the constant
variance term and is removed before bounding fluctuations.
Symmetrizing any of these tensors can only decrease its norm.

## 4. Self-contained finite-chaos Gaussian approximation

Two Gaussian polynomial identities suffice. First, the product rule is

    I_p(f)I_q(g)=sum_(r=0)^min(p,q)
       r! binom(p,r) binom(q,r)
           I_(p+q-2r)(sym(f tensor_r g)).                       (4.1)

For completeness, multiply the two generating functions
exp(t<Z,v>-t^2||v||^2/2) and exp(u<Z,w>-u^2||w||^2/2).
Their product is the joint Wick generating function times
exp(tu<v,w>). Comparing coefficients gives (4.1) for tensor powers;
polarization and linearity give it for all symmetric tensors.
The same generating function proves the Hermite covariance identity
used in (2.3), as well as the derivative identity
grad I_q(f)=q I_(q-1)(f with one tensor index left free).

Second, let L=Delta-Z dot grad be the finite-dimensional Gaussian
Ornstein--Uhlenbeck operator. Gaussian integration by parts gives

    E[(-L U)V]=E[grad U dot grad V].                            (4.2)

Set U_Q=sum_(q odd<=Q) I_q(f_q)/q. Since -L acts as multiplication
by q on degree-q Hermite polynomials, -L U_Q=F_Q. Consequently

    Gamma_Q=grad F_Q dot grad U_Q,
    E[F_Q h(F_Q)]=E[h'(F_Q) Gamma_Q],
    E Gamma_Q=E F_Q^2=:sigma_Q^2.                               (4.3)

The identity for h(x)=exp(itx) follows directly by integration by
parts: the polynomial factors and their derivatives are integrable
against the Gaussian density, and the exponential has bounded modulus.

Applying (4.1) to the two gradients expands Gamma_Q as

    sum_(p,q odd<=Q) sum_(r=1)^min(p,q)
      p(r-1)! binom(p-1,r-1) binom(q-1,r-1)
           I_(p+q-2r)(sym(f_p tensor_r f_q)).                   (4.4)

The terms p=q=r are constants p!||f_p||^2; together they are
sigma_Q^2. Every other contraction in (4.4) is covered by (3.1)
or (3.2). There are only finitely many terms for fixed Q, and the
Gaussian polynomial isometry and the triangle inequality therefore give

                     Var(Gamma_Q)<=K_(Q,C)/n,                   (4.5)

where K_(Q,C) is finite and independent of n, R, and a.
No independence of distinct chaos components is assumed.

Let phi_Q(t)=E exp(itF_Q). From (4.3), exactly,

    phi_Q'(t)+sigma_Q^2 t phi_Q(t)
       =-t E[(Gamma_Q-sigma_Q^2)exp(itF_Q)].                    (4.6)

Solving this scalar differential equation with phi_Q(0)=1 and using
Cauchy--Schwarz gives, for every real t,

    |phi_Q(t)-exp(-sigma_Q^2 t^2/2)|
           <=(t^2/2)sqrt(Var Gamma_Q)
           <=(t^2/2)sqrt(K_(Q,C)/n).                            (4.7)

For t>=0 the integrating-factor kernel is
exp[-sigma_Q^2(t^2-s^2)/2]<=1 for 0<=s<=t; negative t follows by
complex conjugation. Thus no inverse variance or dimension-dependent
Gaussian approximation constant is used, including when variances
approach zero.

## 5. From characteristic functions to uniform absolute moments

For fixed Q, the variance sigma_Q^2 stays in [0,C]. Given any
sequence of admissible inputs with n tending to infinity, select a
subsequence on which sigma_Q^2 converges to some v in [0,C].
Equation (4.7) and the continuity theorem for characteristic functions
give convergence in distribution along that subsequence to N(0,v).
The uniform second-moment bound in (2.3) implies uniform integrability
of |F_Q|: E[|F_Q|; |F_Q|>M]<=C/M. Hence

                    E|F_Q|-sqrt(kappa) sigma_Q ->0.             (5.1)

This conclusion is uniform over all R,a at the given n. Otherwise
one could choose a sequence violating uniformity, select its variance-
convergent subsequence, and contradict the preceding argument. Write
its uniform error as e_(Q,C)(n), which tends to zero for each fixed Q.

Since the omitted chaoses are orthogonal to those retained,

    0<=sigma^2-sigma_Q^2=E(F-F_Q)^2<=C tau_Q,
    |sigma-sigma_Q|<=sqrt(C tau_Q).

The elementary inequality ||x|-|y||<=|x-y| and (2.3) now yield

    |E|F|-sqrt(kappa) sigma|
      <=e_(Q,C)(n)+(1+sqrt(kappa))sqrt(C tau_Q).                 (5.2)

First take n to infinity with Q fixed, and then Q to infinity.
This proves the uniform lemma (2.1). In particular the proof does
not attempt to bound E|F| below using variance alone: the Gaussian
approximation is the necessary bridge between these quantities.

## 6. Apply the lemma to the actual cross columns

Return to B,d,m in Section 1, and put

    R=BB^T/n,       E=BB^T/d^2=mR,
    G Gaussian with covariance R,       X=sign(G),
    h_j=(B^T X)_j.

The complete sign entries imply diag R=1. Also
0<=R<=I/m<=I/m_0. The covariance matrix of X is

    C_X=sum_(q odd>=1) c_q^2 R^{circ q},
    0<=C_X<=I/m.                                              (6.1)

For the jth column b_j of B define

    ell_j=b_j^T R b_j/d^2,
    t_j=Var(h_j)/d^2=b_j^T C_X b_j/d^2.

Every column has exactly n entries of squared magnitude one.
It follows from ||B||op<=d that

    0<=ell_j<=1,       (1/n)sum_j ell_j=1-epsilon,
    kappa ell_j<=t_j<=1.                                     (6.2)

For the upper bound on ell_j use
b_j^T R b_j=||B^T b_j||^2/n<=d^2||b_j||^2/n=d^2.
The lower bound on t_j retains the q=1 summand in (6.1);
all the other summands are positive semidefinite. The upper bound
uses C_X<=I/m and ||b_j||^2=n.

Taking the average of these variances retains more information than
the q=1 bound alone:

    (1/n)sum_j t_j=tr(E C_X)/n
       =sum_(q odd>=1) c_q^2 (m/n)sum_(i,j) R_ij^(q+1)
       >=kappa(1-epsilon)+(1-kappa)m.                          (6.3)

The q=1 term is exactly kappa(1-epsilon). For every odd q>=3,
q+1 is even, so all displayed summands are nonnegative and the n
diagonal summands already contribute m. This is the retained higher-
Hermite variance gain, derived from the actual complete sign entries.
The infinite sum is justified either by its nonnegative scalar terms
or by the L2 covariance limit in (6.1).

For each column j, h_j/sqrt(n) has exactly the coefficient pattern
required in Section 2. Its uniform absolute-moment conclusion gives

    |E|h_j|-sqrt(kappa) d sqrt(t_j)|<=sqrt(n)e_(1/m_0)(n),

with an error independent of j. Since beta(B)>=E||B^T X||_1,

    beta(B)/(nd)
       >=(1/n)sum_j sqrt(kappa)sqrt(t_j)
                         -sqrt(m)e_(1/m_0)(n).                 (6.4)

Only the separate uniform marginal estimates are summed here.
Neither independence of the h_j nor a joint column-sum CLT is used.

## 7. The explicit flatness penalty and the strict gain

Put b=sqrt(kappa)/(1+sqrt(kappa)). For every t in [0,1],

    sqrt(kappa)sqrt(t)
       >=kappa+b(t-kappa)-(1-b)(kappa-t)_+.                    (7.1)

If t<=kappa the right side equals t, which is at most
sqrt(kappa)sqrt(t). If t>=kappa, (7.1) is the chord bound for the
concave square root between t=kappa and t=1.

By (6.2),

    (kappa-t_j)_+<=kappa(1-ell_j),
    (1/n)sum_j (kappa-t_j)_+<=kappa epsilon.                   (7.2)

Averaging (7.1) and applying (6.3)--(7.2) therefore gives

    (1/n)sum_j sqrt(kappa)sqrt(t_j)
       >=kappa+b[(1-kappa)m-kappa epsilon]
                                      -(1-b)kappa epsilon
       =kappa+(sqrt(kappa)-kappa)m-kappa epsilon.              (7.3)

Insert this in (6.4), use sqrt(m)<=1, and rename the uniform error
to obtain (1.1). For exact flat nonzero singular values,
B B^T B=d^2 B and epsilon=0, proving (1.2).

Finally, the elementary spectral inequalities
tr[(B^T B)^2]<=d^2 tr(B^T B)=d^2 n^2 and
tr[(B^T B)^2]>=tr(B^T B)^2/n=n^3 verify the claimed range
0<=epsilon<=1-m. They also show that epsilon=0 is exactly flatness
at the permitted endpoint d, apart from zero singular values.

## 8. What this says about a pure-cross active normalization

If an ACTUAL paired signing has p=q_A=0 and c=Phi(K), then
c=beta(B): the cross-only norm is at most Phi(K), whereas that
actual state attains c through its cross term. If its scalar scale
is d with d>=||B||op, the corresponding ratio is u=c/(nd).
Under the additional hypotheses of Section 1, (1.1) then requires

        u>=kappa+(sqrt(kappa)-kappa)m-kappa epsilon-o(1).        (8.1)

For exact flatness, or epsilon tending to zero with m bounded below,
this excludes u=kappa+o(1). In particular the formal trace model
with cross law (1-m)delta_0+m delta_1 and u=kappa cannot be realized
in such an actual scalar, operator-bounded pure-cross setting.

This is a new actual-entry restriction beyond its retained cubic
cross moment inequality. It does not yet transfer (8.1) to every
near-scalar weighted cross measure or to an arbitrary diagonal
majorizer with potentially exceptional large diagonal entries.
Nor does it evaluate the original ellipsoid upper at every possible
active cell or prove the original inequality. Those are separate
remaining implications, not consequences hidden in the scalar CLT.

## 9. Primary-source cross-check and contribution record

The following primary sources were consulted to verify the standard
chaos-limit route independently of the self-contained argument above:

- Nualart and Peccati, *Central limit theorems for sequences of multiple
  stochastic integrals*, Theorem 1, page 3:
  https://arxiv.org/pdf/math/0503598 . Its normalized fixed-chaos
  criterion equates vanishing nontrivial self-contractions with a
  Gaussian limit and convergence of fourth moments.
- Noreddine and Nourdin, *On the Gaussian approximation of vector-valued
  multiple integrals*, Theorem 1.1, pages 1--2:
  https://arxiv.org/pdf/1009.1310 . Fixed chaos orders may include one;
  convergent covariances and vanishing marginal fourth cumulants imply
  joint Gaussian convergence. Singular limiting covariances are allowed.

Both theorem statements were read from the actual PDFs. Downloaded
copies and their SHA256 hashes are:

    original_mo_flat_cross_clt_primary_nualart_peccati.pdf
    b8c80cb9b03b016da921316239d0acaed6594fc9303a6176c860fe25db173ea9

    original_mo_flat_cross_clt_primary_noreddine_nourdin.pdf
    bb07986ac769c63ae02596b470903ba3c4f89c5d83d494c4d3eff2c98d6cb259

Those sources are optional cross-checks; Sections 3--5 establish the
needed scalar limit through Gaussian polynomial identities and the
characteristic-function equation without importing their theorems.

The optimized-profile proof worker originated the actual flat-cross
variance/chord gain and delocalized-chaos CLT bridge and authored this
note. The exact worker independently checked that proposal and supplied
the mixed-contraction/characteristic-equation route and the robust
epsilon extension. The root independently checked these contributions
and requested their integration. Consequently root/exact checks of
this combined artifact are contributing-author checks, not substitutes
for a genuinely independent full-source review. No mathematical
computation was run by the proof worker in deriving this result.
