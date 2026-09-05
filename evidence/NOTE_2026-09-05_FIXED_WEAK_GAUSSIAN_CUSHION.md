# Fixed-strength weak Gaussian rounding: a conditional-noise obstruction

2026-09-05. All-orders, method-scoped theorem; original convergence OPEN.

This strengthens the earlier vanishing-strength exclusion. It concerns
actual log pressure and excludes a nonempty interval of FIXED positive
rounding strengths at suitable fixed temperatures. It does not use an
o(n)-relative-entropy hypothesis, which is false at fixed strength on
norm-capped hosts. The singular full-strength law rho=1 remains open.
No computation or numerical spin-glass constant is used.

## 1. Objects and statement

Use the proved quenched note
`NOTE_2026-09-05_IID_QUENCHED_CROSS_OBSTRUCTION.md`, including its
Gaussian pure-cross lower bound and its explicit constant

    K0 = 4/(3 sqrt(pi)) > 1/sqrt(2).

Fix c>0, 0<t<=1, and 0<=rho<1 independently of n. Set

    beta=c/sqrt(n), eta=beta sqrt(1-t/2),
    gamma=beta sqrt(t/2), d=n^2,
    L(B)=log E_(x,y) exp(gamma x^T B y).

Here x,y are independent uniform sign vectors of length n. For every
complete internal host A, the paired pressure is

    F_(A,B)=log E_(x,y) cosh[eta(Q_A(x)-Q_A(y))+gamma x^T B y].

The block-flip identity in the quenched note proves POINTWISE that

    min_A F_(A,B) >= L(B).                                      (1)

Let S be ANY d by d positive-semidefinite matrix with diagonal one,
possibly singular. Generate U~N(0,S) and independent iid standard normal
g in R^d, and reshape the signs of

    sqrt(rho) U + sqrt(1-rho) g

into B. Its Gaussian covariance is (1-rho)I+rho S. Then

    E min_A F_(A,B)
      >= [c sqrt(2t) K0 - 2 log2
                        - c^2 t arcsin(rho)/(2 pi)]n - o_(c,t)(n). (2)

The error is uniform in S and rho. The covariance-generating host,
if any, is held fixed while producing this law; the internal host in
(1) can still be selected after seeing B.

Define

    delta(t)=sqrt(2t)K0-1,
    Deltarho(c,t)=c delta(t)-2 log2-c^2 t arcsin(rho)/(2 pi).       (3)

Since 2R_n(beta)<=cn+o_c(n), a positive Deltarho excludes the desired
mean comparison with the optimized paired endpoint. For the ACTUAL
centered Gaussian construction, it also gives exponentially small
probability of a successful outcome, uniformly over generating hosts.
The probability claim is proved in Section 4.

For any fixed t with delta(t)>0, there is a suitable fixed c whenever

    0 < arcsin(rho) < pi delta(t)^2/(4t log2).                   (4)

Indeed the quadratic in c in (3) has a positive maximum exactly under
(4), attained at c=pi delta(t)/(t arcsin(rho)). This is a nonempty
fixed-positive-strength interval at t=1. At rho=0 the original iid
criterion c delta(t)>2log2 applies. No decimal approximation is needed.

## 2. Conditional independent replacement and removal of the mean

Conditional on U, the d signs B_e are independent with means

    m_e = m(U_e),
    m(u)=2 Phi_std(sqrt(rho/(1-rho)) u)-1,
    w_e=1-m_e^2,

where Phi_std is the standard normal distribution function. Replace
them, one coordinate at a time, by independent Gaussian variables with
the SAME conditional means and variances, m_e+sqrt(w_e)z_e.

For any real matrix argument, the third coordinate derivative of L
has absolute value at most 2 gamma^3: the corresponding spin observable
x_i y_j takes values +-1 and its third centered moment is
-2a(1-a^2), where a is its Gibbs mean. The centered third absolute
moment of a sign with mean m is 1-m^4<=1. The Gaussian counterpart
has third absolute moment at most E|z|^3=2sqrt(2/pi). Taylor's theorem
and matching the first two moments therefore give, uniformly in U,

    |E[L(B)|U] - E_z L(M+R)|
      <= (gamma^3/3)(1+2sqrt(2/pi)) d = O_(c,t)(sqrt(n)),        (5)

where M=(m_e) and R=(sqrt(w_e)z_e).

The pure-cross pressure L is convex and globally even, L(-B)=L(B).
Since R and -R have the same law, the function

    M -> E_z L(M+R)

is also convex and even, with its minimum at M=0. Consequently

    E_z L(M+R) >= E_z L(R).                                    (6)

This is an actual Gaussian replacement followed by convexity. It does
not treat the original centered Bernoulli residual as symmetric, and
does not move a logarithm outside the expectation over cross blocks.

For independent centered Gaussian entries with coordinate variances v_e,
Gaussian integration by parts gives

    0 <= partial_(v_e) E_z L((sqrt(v_e)z_e)_e)
       = (gamma^2/2) E_z[1-<x_i y_j>^2] <= gamma^2/2.

The formula extends to zero variances by continuity. Increasing all w_e
to one and integrating the upper bound yields

    E_z L(R) >= E_z L(z) - (gamma^2/2) sum_e m_e^2.             (7)

Every U_e is standard normal, regardless of S. Two independent noise
copies conditional on U_e have Gaussian correlation rho, so the arcsine
sign identity gives

    E m(U_e)^2 = a_rho := (2/pi) arcsin(rho).                   (8)

The already reviewed Gaussian pure-cross bound is

    E_z L(z) >= [c sqrt(2t) K0 - 2 log2] n - o_(c,t)(n).

Combining (1), (5)--(8), and gamma^2=c^2 t/(2n) proves (2).
No positive spectral bound on S is used for this mean conclusion.

## 3. Concentration of the conditional variance loss

For 0<rho<1 put

    K_rho=2sqrt(2/pi) sqrt(rho/(1-rho)).

The derivative of f(u)=m(u)^2 has absolute value at most K_rho.
Writing U=S^(1/2)z shows that

    T(z)=sum_e f((S^(1/2)z)_e)

is K_rho sqrt(d ||S||_op)-Lipschitz. Therefore, for epsilon>0,

    P(T >= a_rho d + epsilon d)
       <= exp[-epsilon^2 d/(2 K_rho^2 ||S||_op)].               (9)

For completeness, the Gaussian concentration inequality needed here
follows directly from the heat martingale. For a smooth ell-Lipschitz
function h of a standard Gaussian vector, let W_s be standard Brownian
motion and M_s=E[h(W_1)|F_s]=(P_(1-s)h)(W_s). Heat averaging preserves
the gradient bound ell, so Ito's formula gives dM_s=v_s dW_s with
||v_s||<=ell. The exponential martingale implies
E exp(lambda(M_1-M_0))<=exp(lambda^2 ell^2/2). Exponential Markov
and lambda=r/ell^2 give the tail exp(-r^2/(2ell^2)). Smooth bounded
approximation, or directly the bounded smooth function T here, handles
the endpoints. This argument requires no independent coordinates in U.

Conditional on U, flipping any one of the independent B_e changes L
by at most 2gamma. The bounded-difference MGF, valid also for biased
independent signs, thus gives

    P(L(B)<=E[L(B)|U]-r | U)
       <= exp[-r^2/(2d gamma^2)].                             (10)

No concentration of F or of an adaptively selected host is necessary;
the pointwise host-free lower bound (1) is retained throughout.

## 4. Actual centered laws and rare-outcome consequence

For any complete generating host A0 and the admissible actual centering
alpha from the Gaussian-sign information note, put

    H=A0 tensor A0-alpha(A0 tensor I+I tensor A0),
    mu=-lambda_min(H), S=I+H/mu.

The reviewed exact bounds give

    mu >= (n-1)/2, ||H||_op <= 2||A0||_op^2,
    ||A0||_op <= n-1,

and hence the UNIFORM estimate

    ||S||_op <= 1+4(n-1)=4n-3.                                (11)

Thus this S is admissible in Section 1, singularity allowed, and the
generated covariance is precisely I+rho H/mu. No norm-minimizer,
pressure-minimizer, or extra spectral hypothesis is required for (11).

Fix parameters with Deltarho>0 and choose

    epsilon=Deltarho/(c^2 t).

On the complement of the event in (9), equations (5)--(7) give

    E[L(B)|U] >= (c+3Deltarho/4)n-o_(c,t)(n).

For every deterministic e_n=o(n), sufficiently large n therefore obey

    P(min_A F_(A,B) <= 2R_n(beta)+e_n)
      <= exp[-epsilon^2 n^2/(2K_rho^2(4n-3))]
             + exp[-Deltarho^2 n/(16c^2 t)].                  (12)

Indeed 2R_n(beta)+e_n<=cn+o(n), and a downward conditional deviation
of at least Deltarho n/4 is then necessary. Substitute this in (10).
Both terms in (12) are exponentially small in n at the fixed parameters.
At rho=0 the conditional variance loss vanishes deterministically and
the first term is omitted.

The bound is uniform over generating hosts and their admissible
centerings. Consequently exp(o(n)) proposals whose marginals have
these laws, or mixtures of them formed by choosing a generating host
before drawing its conditional law, have vanishing success probability
by a union bound. Independence among proposals is not needed. The
internal host may be optimized separately after each draw by (1).
This is not a lower entropy assertion for mixtures over hosts.

The theorem does not exclude a specially selected exponentially rare
outcome, a stronger fixed rho outside (3), the actual Gram--Schmidt law,
or rho=1. In particular, it does not prove a cross-order inequality
or settle the original MO limit.
