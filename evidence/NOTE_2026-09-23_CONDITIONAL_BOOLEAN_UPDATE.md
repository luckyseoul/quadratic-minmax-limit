# Exact conditional moments after a Boolean correction layer

2026-09-23. Author-reviewed finite identities and polynomial-arithmetic
evaluation, not a convergence theorem or a new universal numerical bound.

The September 6 [adaptive-source proof](NOTE_2026-09-06_ALL_LAW_ADAPTIVE_NUCLEAR_GAIN.md),
Section 7, already uses independent coordinate masks after conditioning on
the actual source. We reuse that update. The additional result here is an
exact conditional evaluation of its next field magnitudes, instability
count, and target quadratic self-action. In particular, the updated law
is not identified with a Gaussian sign phase.

The exact-object check was carried forward from the preceding research
turn: active source/update notes, ten preserved local branch heads, four
worktrees, and the indexed September 6 archive. The independent-mask
identity is prior work. No prior grouped row/pair evaluator was located;
the September 23 interior-gap lemma is a separate, explicitly reused input.

## 1. A fixed conditional source and exact row distributions

Let A be a complete symmetric zero-diagonal signing, x a fixed Boolean
vector, and B_j independent Bernoulli variables with probabilities p_j.
Set

    X'_j = x_j(1-2B_j),   m_j = x_j(1-2p_j),   h_i = (Ax)_i,
    H'_i = (AX')_i = h_i-2 sum_(j!=i) A_ij x_j B_j.       (1)

The low-level identity allows all p_j in [0,1]. For the existing mean
update, use p_j=0 where x_j h_j>=0; the nonzero probabilities may differ.

Partition the sum in (1) according to A_ij x_j=+1 or -1. The two counts
are independent Poisson-binomial variables. For rational p_j=u_j/v_j,
their probability-generating polynomials have integer factors

    ((v_j-u_j)+u_j z)/v_j.                               (2)

Thus convolution, followed by subtraction of the two counts, gives the
entire rational law of H'_i. Zero and unit probabilities are allowed.

Define Y'_i=sign(H'_i) at a nonzero field and Y'_i=X'_i at a zero field.
Since A_ii=0, X'_i is independent of H'_i. Writing P_i^+, P_i^- and P_i^0
for the three field-sign probabilities gives

    E Y'_i = P_i^+ - P_i^- + m_i P_i^0,
    q_i := P(X'_i H'_i<0)
         = [(1+m_i)P_i^-+(1-m_i)P_i^+]/2.                (3)

The row law also determines E|H'_i| and its first two moments exactly.
Independence of the *updated spins*, not of the fields, gives

    a := E Q_A(X') = Q_A(m),
    b := E ||AX'||_1 = sum_i E|H'_i|,
    q := E K' = sum_i q_i,  K'=#{i:X'_i H'_i<0}.          (4)

## 2. Pair fields: four classes retain dependence

Fix distinct i,j and exclude them from the remaining coordinates.
Partition k by (A_ik x_k,A_jk x_k) in {++, +-, -+, --}. Let

    U = sum_(++) B_k - sum_(--) B_k,
    V = sum_(+-) B_k - sum_(-+) B_k.

The two variables U,V and the endpoint masks B_i,B_j are mutually
independent. Their distributions are given by (2), including for unequal
probabilities. Conditional on the endpoint masks,

    H'_i = h_i-2[U+V+A_ij x_j B_j],
    H'_j = h_j-2[U-V+A_ij x_i B_i].                      (5)

The fallback signs at zero are x_i(1-2B_i) and x_j(1-2B_j). Summing the
product of the two target signs in (5) over U,V and the at most four
endpoint choices gives E[Y'_i Y'_j] exactly. Endpoint conditioning is
essential: a field omits its *own* mask, not the other endpoint's mask.
In general E[Y'_i Y'_j] is not E[Y'_i]E[Y'_j]. Consequently,

    c := E Q_A(Y') = sum_(i<j) A_ij E[Y'_i Y'_j]          (6)

is now evaluated, instead of replaced by the unknown lower -Phi(A).

Each univariate convolution uses O(n^2) integer arithmetic operations;
each pair has O(n^2) support pairs. The n(n-1)/2 pair evaluations and n
row evaluations therefore cost O(n^4) such operations, rather than an
enumeration of 2^k masks for k nontrivial probabilities. Rational bit
complexity is not constant: denominators divide products of the supplied
v_j, and the integer arithmetic cost depends on their bit lengths.
The implementation bounds its reusable distribution cache to 512 entries.
No GPU speedup or measured exponential timing ratio is asserted.

## 3. A second mean step with evaluated quadratic self-action

For 0<=t<=1 put Z=(1-t)X'+tY'. Quadratic expansion and
(AX')^T Y'=||AX'||_1 give exactly

    C(t) := E Q_A(Z) = (1-t)^2 a+t(1-t)b+t^2 c.          (7)

Independent Boolean rounding of Z has this expected energy. This is a
valid actual second correction, not a reapplication of the old Gaussian
source lemma. It does not assert that simultaneous full correction t=1
always raises the energy.

Write F=Phi(A)>0. For 0<=t<=1/2, the [finite interior gap](NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md)
applied pointwise to Z gives

    F >= C(t)+t^2 E g(K')/(24F),
    g(0)=g(1)=0,  g(k)=k^3 for integers k>=2.             (8)

Let G be the linear interpolation of g between consecutive integers.
Its slopes are 0,8,19,37,... and increase; hence G is convex. Because
K' is integer-valued, Jensen's inequality gives E g(K')>=G(E K')=G(q).
This step does not suppose independence of the instability indicators.
Thus all coefficients in the following finite norm certificate are
conditionally evaluated rationals when p and t are rational:

    F^2-C(t)F-D(t)>=0,   D(t)=t^2 G(q)/24.               (9)

The nonnegative root is a lower bound on F; the code encloses it by exact
rational bisection. The elementary bounds F>=|a| and F>=|c| remain valid.

For comparison, replacing c by -F in (7) gives the older certificate

    (1+t^2)F^2-[(1-t)^2 a+t(1-t)b]F-D(t)>=0.             (10)

The root of (9) alone need not dominate that of (10) if c is negative.
Retaining F>=|c| fixes the comparison: if R is the larger of |c| and
the root of (9), then the left side of (10) at R is the left side of
(9) plus t^2 R(R+c), hence is nonnegative. Therefore R is at least the
nonnegative root of (10). This is a comparison of available finite
certificates, not an improved universal constant.

## 4. Outer laws, actual limitation, and verification scope

For any source law X, the above identities hold conditional on X=x,
including when p=p(x). Taking the outer expectation preserves (7)--(9).
In (9), either average the conditional G(q(x)) or apply its convexity
again to the outer mean count. This is an exact tower-expectation
statement; it does not evaluate the outer expectation for free.

After further adaptive steps, the unconditional spins generally form a
mixture of product laws. The first two target moments do not specify that
mixture or close subsequent sign-of-field expectations. No compressed
outer-law evolution, uniform asymptotic gain, Dini error, one-vertex
extension, or cross-order comparison follows here. These remain actual
research requirements, not hypotheses supplied by the tests.

`src/boolean_mean_update.py` implements (1)--(10) with integer/Fraction
arithmetic. `tests/test_boolean_mean_update.py` compares it to independently
implemented direct mask enumeration on fixed order-2, order-5 and order-6
fixtures, including nonuniform probabilities, all four pair classes,
deterministic masks, zero fields, switching and root brackets. This is not
a census of signings or a formal verification of the all-orders argument.
The separate receipt records a single retained order-26 input evaluation,
without rescoring its already certified norm or attempting a repair scan.
Review is by the author, not an independent human or proof assistant.
Original MO convergence remains OPEN; the universal lower B is unchanged.
