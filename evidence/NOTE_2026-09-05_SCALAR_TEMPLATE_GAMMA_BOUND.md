# A quantitative completion bound for scalar-optimal finite templates

Status: proved finite-template obstruction. This is a restriction on the
particular upper-certificate functional Gamma below. It is not a lower bound on
the Boolean norm of a constructed large sign matrix, is not an attainability
claim for Gamma, and does not settle the original MO convergence question.

The proof is analytic. Its displayed finite rational inequalities were
verified once on soulkiller using only fractions.Fraction: all 28 checks
passed. No floating-point arithmetic, numerical integration, matrix search,
or optimization was performed for this note. The terminating decimals below
denote exact rationals, with finite-sum enclosures given in Section 6.

## 1. Setting and conclusions

Let C be a real p by p matrix, let q = ||C||op > 1, and suppose its real vector
SDP has the scalar-dual optimal value

    tau(C) = max sum_ij C_ij <u_i,v_j> = p q,

where the maximum is over unit vectors u_i and v_j in a common Euclidean
space. Define the finite soft-spin completion functional

    Gamma(C) = max_{a,b in [-1,1]^p}
      [a^T C b / p + sqrt((1-||a||^2/p)(1-||b||^2/p))].          (1.1)

In particular Gamma(C) >= 1, by choosing a=b=0. All maxima here are attained:
the vector SDP can be written over the compact set of correlation matrices of
size 2p, and the cube in (1.1) is compact.

Let G be a standard real Gaussian, and let f: R -> [-1,1] be a measurable odd
function. Put

    c = E[G f(G)],                v = E[f(G)^2],
    k = E[(G^2-1) f(G)^2] / 2,
    R = Var(f(G)^2) - 2 k^2 >= 0,
    D(q) = q(2 c^2-v) - v.

Then

    Gamma(C) >= 1 + D(q)
      - sqrt(R (Gamma(C)-1) / (2(q-1))).                       (1.2)

When D(q) >= 0 this gives the explicit lower bound

    Gamma(C) >= 1 + ((sqrt(h+4D(q))-sqrt(h))/2)^2,
    h = R/(2(q-1)).                                            (1.3)

The concrete choice f(x)=clip(x,-1,1) proves

    q >= 5/2  ==>  Gamma(C) > sqrt(2).                         (1.4)

Consequently every such finite template satisfying Gamma(C) <= sqrt(2) has
q < 5/2. If a template normalization also writes u=sqrt(2)/q, it necessarily
has u > 2 sqrt(2)/5. No conclusion about a matrix's actual Boolean norm is
drawn merely from a failure of its Gamma upper certificate.

## 2. Equality in the scalar-dual SDP gives matched frames

Choose an optimal vector representation and write U,V for the p by r matrices
whose rows are u_i^T,v_j^T. Their row norms are one, so

    ||U||F^2 = ||V||F^2 = p.

Equality holds throughout

    p q = tr(U^T C V) <= ||U||F ||C V||F <= p q.

Equality in the first Cauchy--Schwarz inequality, together with the positive
value p q, gives C V = q U. Applying the same argument to the transposed
objective gives C^T U = q V. Hence

    M := U^T U = V^T V,
    tr M = p.                                                 (2.1)

Indeed U^T C V = q U^T U and also U^T C V = q V^T V.

Set

    mu = lambda_max(M)/p,
    s = tr(M^2)/p^2.

Because M is positive semidefinite and has trace p,

    0 < s <= mu <= 1.                                         (2.2)

For a unit top eigenvector z of M choose a=Uz and b=Vz. Every coordinate of
these vectors lies in [-1,1], their squared norms both equal p mu, and
a^T C b=q p mu. Substitution into (1.1) proves

    Gamma(C) >= 1+(q-1) mu.

Thus

    s <= mu <= (Gamma(C)-1)/(q-1).                            (2.3)

This elementary linear rounding will control the fluctuation penalty of the
nonlinear rounding. No bound on the rank of U or V is assumed.

## 3. The Gaussian cross term

Let g be a standard Gaussian vector in R^r, and put X=Ug, Y=Vg. Every X_i and
Y_j is a standard one-dimensional Gaussian. The matched-frame identity gives
the pointwise identity

    ||X||^2 = g^T M g = ||Y||^2.                              (3.1)

Take the admissible random cube vectors a_i=f(X_i), b_j=f(Y_j). Write
r(x)=f(x)-c x. Then E[G r(G)]=0 and E[r(G)^2]=v-c^2.

For any jointly standard Gaussian pair (X_i,Y_j), including correlations
plus or minus one, Gaussian conditional expectation gives

    E[X_i r(Y_j)] = E[r(X_i) Y_j] = 0.

Therefore

    E[a^T C b] = c^2 E[X^T C Y] + E[r(X)^T C r(Y)].

The first term is c^2 tr(U^T C V)=c^2 p q. The operator norm and ordinary
Cauchy--Schwarz give

    E[r(X)^T C r(Y)]
      >= -q E[||r(X)|| ||r(Y)||]
      >= -q sqrt(E||r(X)||^2 E||r(Y)||^2)
      = -p q (v-c^2).

Consequently

    E[a^T C b]/p >= q(2c^2-v).                               (3.2)

No assertion that the nonlinear rounding itself preserves optimal SDP
correlations is needed.

## 4. Exact quadratic cancellation and the completion term

Put

    A = ||a||^2/p,             B = ||b||^2/p.

These random variables lie in [0,1], and E A = E B = v. For every A,B in
[0,1],

    sqrt((1-A)(1-B)) >= 1-max(A,B)
      = 1-(A+B)/2-|A-B|/2.

It follows that

    E sqrt((1-A)(1-B)) >= 1-v-(1/2) E|A-B|.                  (4.1)

Define

    H(x) = f(x)^2-v-k(x^2-1).

The function H is even and orthogonal in Gaussian L2 to both 1 and x^2-1.
Its squared L2 norm is exactly R. In the normalized probabilists' Hermite
basis it therefore has the expansion

    H(x) = sum_{j>=2} alpha_j He_{2j}(x)/sqrt((2j)!),
    sum_{j>=2} alpha_j^2 = R.                                (4.2)

For jointly standard Gaussian Z,W of correlation rho, Hermite orthogonality
gives

    E[H(Z)H(W)] = sum_{j>=2} alpha_j^2 rho^(2j).

In particular, even when rho is negative,

    0 <= E[H(Z)H(W)] <= R rho^2.                             (4.3)

For completeness, the identity used here follows by expanding
E[exp(tZ-t^2/2) exp(sW-s^2/2)] = exp(rho t s): the coefficient comparison
gives E[He_j(Z)He_l(W)] = 1_{j=l} j! rho^j. Apply it first to finite Hermite
sums, then pass to the Gaussian L2 limit by Cauchy--Schwarz. The normalized
Hermite polynomials form the usual complete orthonormal basis of Gaussian
L2. The endpoint correlations cause no difficulty in this limit argument.

The constant terms and the entire quadratic terms in A-B cancel, the latter
by (3.1). Hence, pointwise,

    A-B = (1/p) sum_i H(X_i) - (1/p) sum_j H(Y_j).             (4.4)

By (4.3) the covariance of the two averages on the right is nonnegative. Also

    Var((1/p)sum_i H(X_i))
      <= (R/p^2) sum_{i,l} <u_i,u_l>^2
      = R tr(M^2)/p^2 = R s.

The same bound holds for the Y average, with the same M. Thus

    E[(A-B)^2] <= 2 R s,
    (1/2) E|A-B| <= sqrt(R s/2).                             (4.5)

Equations (4.1) and (4.5) yield the quantitative completion estimate

    E sqrt((1-A)(1-B)) >= 1-v-sqrt(R s/2).                   (4.6)

The nonnegative cross covariance in this argument is a termwise consequence
of the even Hermite expansion, not a use of a qualitative Gaussian
correlation principle. No unproved equality of the clipped norms is used.

## 5. Combination and the explicit threshold

Since the maximum in (1.1) dominates the expected objective of every random
admissible pair, (3.2), (4.6), and (2.3) prove (1.2).

To derive (1.3), set x=Gamma(C)-1 >= 0 and h=R/(2(q-1)). Then (1.2) says
x+sqrt(hx) >= D(q). When D(q)>=0, solving this quadratic inequality in
sqrt(x)>=0 gives (1.3).

Now choose f(x)=clip(x,-1,1), and write

    P = P(|G|<=1),             phi = exp(-1/2)/sqrt(2 pi).

Integration by parts in the standard Gaussian density gives

    E[G^2 1_{|G|<=1}] = P-2 phi,
    E[G^4 1_{|G|<=1}] = 3P-8 phi.

Also E[|G|1_{|G|>1}]=2 phi. These identities imply

    c = P,
    v = 1-2 phi,
    E[f(G)^4] = 1+2P-8 phi,
    k = P-2 phi,
    R = 1+2P-8 phi-(1-2 phi)^2-2(P-2 phi)^2.                 (5.1)

The exact enclosures established in Section 6 are

    0.68268 < P < 0.68270,
    0.24197 < phi < 0.24198.                                 (5.2)

They imply 2P^2-v>0, so D(q) is strictly increasing in q. At q=5/2,

    D(5/2) = 5P^2-7/2+7 phi
      > 5(0.68268)^2-7/2+7(0.24197)
      = 0.524049912 > 0.524.                                 (5.3)

For the variance remainder, (5.1)-(5.2) give

    R < 0.42964-(0.51604)^2-2(0.19872)^2
      = 0.0843634416 < 0.0844.                               (5.4)

Suppose q>=5/2 and Gamma(C)<=sqrt(2). Put x=Gamma(C)-1. Since
1.415^2>2, one has x<0.415. Moreover

    R x/(2(q-1)) < (0.0844)(0.415)/3 < (0.109)^2.

Using (5.3) in (1.2) now forces

    x > 0.524-0.109 = 0.415,

which is a contradiction. This proves (1.4).

In fact the same arithmetic gives the uniform bound Gamma(C)>283/200 for
q>=5/2. Under Gamma(C)<=283/200 one has x<=0.415, and the strict bound
R<0.0844 still gives sqrt(R x/(2(q-1)))<0.109. Equation (1.2) again forces
x>0.415. Since (283/200)^2>2, this is a fixed positive margin above sqrt(2).

## 6. Exact elementary enclosures, with no numerical oracle

All endpoints in this section are rationals. Define the finite sums

    T_m(z) = sum_{j=0}^m (-1)^j z^(2j+1)/(2j+1),
    E_m = sum_{j=0}^m (-1)^j/(2^j j!),
    I_m = sum_{j=0}^m (-1)^j/(2^j j!(2j+1)).

Alternating-series bounds for arctan, exp(-1/2), and the integrated power
series of exp(-x^2/2) on [0,1] show that odd-indexed sums are strict lower
bounds and even-indexed sums are strict upper bounds.

The elementary Machin identity pi=16 arctan(1/5)-4 arctan(1/239) therefore
gives

    16 T_5(1/5)-4 T_2(1/239) < pi
      < 16 T_6(1/5)-4 T_1(1/239).

The following are finite rational inequalities, checkable by multiplying
positive denominators:

    3.1415926 < 16 T_5(1/5)-4 T_2(1/239),
    16 T_6(1/5)-4 T_1(1/239) < 3.1415927;

    0.6065306 < E_9 < E_10 < 0.6065307;
    0.8556243 < I_7 < I_8 < 0.8556245.                        (6.1)

In addition, direct rational multiplication gives

    (0.398942)^2 (2*3.1415927) < 1,
    (0.398943)^2 (2*3.1415926) > 1.

It follows that, with Z=1/sqrt(2 pi) and
I=integral_0^1 exp(-x^2/2) dx,

    0.398942 < Z < 0.398943,
    0.6065306 < exp(-1/2) < 0.6065307,
    0.8556243 < I < 0.8556245.

The final products satisfy the following exact rational inequalities:

    0.24197 < (0.398942)(0.6065306),
    (0.398943)(0.6065307) < 0.24198;

    0.68268 < 2(0.398942)(0.8556243),
    2(0.398943)(0.8556245) < 0.68270.

Since phi=Z exp(-1/2) and P=2 Z I, these prove (5.2). Thus the explicit
threshold is based on displayed finite rational inequalities and elementary
alternating remainders, not floating-point evaluation.

## 7. Exact scope

The theorem applies to every finite real template C satisfying both
||C||op=q and tau(C)=p q. It does not require C to have sign entries, does not
assume symmetric singular spectrum, and does not assert that every such C
comes from a realizable large sign-matrix construction.

In a construction where Gamma(C) is only a proved asymptotic upper bound
for an actual Boolean norm, Gamma(C)>sqrt(2) means that this upper
certificate cannot establish the desired sqrt(2) cap. It does not imply the
actual Boolean norm exceeds sqrt(2). Even a uniform gap for Gamma does not
by itself exclude those actual matrices from the original weak-repair
regime. A transfer from Gamma to an actual Boolean lower bound would need a
separate argument and is not part of this note.
