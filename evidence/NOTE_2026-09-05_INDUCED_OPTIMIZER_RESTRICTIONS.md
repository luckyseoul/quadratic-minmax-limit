# Quantitative induced restrictions of low-norm complete signings

Status: proved quantitative sampling and universality statements, followed by
an obstruction to typical or modest-many uniform restrictions. This does not
prove or disprove convergence of the original MathOverflow sequence. It does
not assume that any conference or Paley signing is globally optimal.

No computation, finite census, or solver is used in this proof.

## 1. Definitions and the bilinear estimate

Let A be a symmetric N by N matrix with zero diagonal and off-diagonal
entries in {−1,+1}. Write

    Q_A(x) = sum_{i<j} A_ij x_i x_j = (1/2) x^T A x,
    M = Phi(A) = max_{x in {−1,+1}^N} |Q_A(x)|,
    beta(A) = max_{x,y in {−1,+1}^N} |x^T A y|.

For every z in [−1,1]^N, choose independent signs X_i with E X_i = z_i.
The zero diagonal gives Q_A(z) = E Q_A(X), so |Q_A(z)| <= M.

For signs x,y, put u=(x+y)/2 and v=(x−y)/2. Symmetry of A gives

    x^T A y = u^T A u − v^T A v = 2(Q_A(u) − Q_A(v)).

Both u and v belong to [−1,1]^N. Consequently

    beta(A) <= 4 M.                                             (1)

A bilinear form achieves its maximum absolute value on a product of cubes
at a pair of vertices. Thus, in particular, for any f,g in [0,1]^N,

    |sum_{i,j} A_ij f_i g_j| <= beta(A) <= 4 M.                  (2)

This argument does not require A to minimize Phi.

## 2. The complete induced-table law, including collisions

Fix integers 2 <= n <= N and put k=binom(n,2). A labelled signing B of order
n is specified by its k signs B_ab, 1 <= a < b <= n. Let U_n denote the
uniform distribution on all 2^k such signings, equivalently independent fair
signs on their edges.

Let (I_1,...,I_n) be uniform among ordered injective tuples from [N]. Define
the induced labelled signing by B^I_ab=A_{I_a I_b}. Let nu denote its law.
This is a uniform n-set with an independent uniform labelling. It is not a
claim about the full table law obtained by sorting the sampled vertex labels.
Its Phi, however, is exactly the Phi of the underlying induced submatrix,
because Phi is invariant under relabelling.

To handle collisions rigorously, first take independent uniform indices
J_1,...,J_n in [N]. Define an auxiliary random signing as follows. For each
position pair a<b with J_a != J_b, use A_{J_a J_b}. For each pair with
J_a=J_b, use a fresh independent fair sign; these auxiliary signs are also
independent of the indices. Let mu be this auxiliary law. For every labelled
pattern B, its probability is exactly

    mu(B) = E_J product_{a<b} (1 + B_ab A_{J_a J_b})/2,           (3)

where A_ii=0. The factors 1/2 at collisions in (3) are precisely the
probabilities supplied by the auxiliary fair signs. Repeated nonloop edges
retain their actual dependence; no independence of those edge values is
being asserted.

Let

    pi_(N,n) = P(some J_a=J_b)
             = 1 − (N)_n / N^n <= k/N,                         (4)

where (N)_n=N(N−1)...(N−n+1). Conditional on no collision, the auxiliary law
is nu. Therefore mu=(1−pi_(N,n))nu + pi_(N,n)kappa for another probability
law kappa, and

    TV(mu,nu) <= pi_(N,n).                                    (5)

Here TV(P,Q)=(1/2)sum_b |P(b)−Q(b)|. The correction in (5) is at law level:
it must not be multiplied by the number 2^k of labelled patterns.

### Proposition 1: quantitative full-law estimate

For every A,N,n as above,

    |mu(B) − 2^(−k)| <= 2 k M / N^2       for every B,          (6)

    |nu(B) − 2^(−k)| <= 2 k M / N^2 + pi_(N,n),               (7)

and

    TV(nu,U_n) <= k 2^k M / N^2 + pi_(N,n).                   (8)

Proof. Order the k position pairs in (3), write their factors as F_1,...,F_k,
and telescope product F_r against product (1/2). Each term has one factor

    F_r − 1/2 = (B_ab/2) A_{J_a J_b}.

Condition on all indices except J_a,J_b. Every other retained edge factor
is either constant, a function of J_a alone, or a function of J_b alone:
the position graph is simple and has only one edge joining a and b. Thus
the remaining product has the form c f(J_a)g(J_b), with c in [0,1] and
f,g taking values in [0,1]. The independent uniform conditional law of
J_a,J_b and (2) bound the absolute conditional expectation of this term by

    beta(A)/(2 N^2) <= 2 M/N^2.

Summing the k terms proves (6). Equation (7) follows from (5). Summing (6)
over the 2^k patterns and dividing by two gives

    TV(mu,U_n) <= k 2^k M/N^2.

Triangle inequality with (5), used only once, proves (8). This controls the
entire joint induced table, not merely its edge marginals or fixed moments.
End of proof.

For later use, set

    D_(N,n)(A) = min(1, k 2^k M/N^2 + pi_(N,n)).                (9)

## 3. All-labelled-pattern universality

### Corollary 2: finite sufficient condition

If

    2 k M/N^2 + pi_(N,n) < 2^(−k),                            (10)

then every labelled complete signing of order n occurs as an induced
labelled submatrix of A. Indeed (7) makes each of its probabilities strictly
positive. The same inequality holds uniformly for every pattern; no union
bound over failures is needed.

In particular, let A_N be any sequence with Phi(A_N) <= C N^(3/2) for all
sufficiently large N, for a fixed finite C. Let n=n(N) satisfy

    n -> infinity,              n^2 = o(log N).                (11)

Then (10) holds for all sufficiently large N. To check its scale, multiply
its left side by 2^k and use (4):

    2^k(2 k M/N^2 + pi_(N,n))
        <= 2 C k 2^k / sqrt(N) + k 2^k/N -> 0.                 (12)

Indeed k log 2 = o(log N). This is stronger than TV convergence alone;
TV=o(1) by itself would not prove that every pattern occurs.

Define the original extremal quantity

    m_n = min_B Phi(B),          alpha_n = m_n/n^(3/2),

where B ranges over all complete signings of order n. Under (10),

    min_{S subset [N], |S|=n} Phi(A[S]) = m_n.                 (13)

The lower bound is the definition of m_n; universality embeds a minimizing
pattern and proves the upper bound. This applies, in particular, to every
sequence of exact globally minimizing A_N. It also applies to every other
sequence satisfying the same O(N^(3/2)) norm bound.

Equation (13) is an existence statement about all restrictions. It does not
assert that a uniform sample finds those minimizing patterns with useful
probability, and it gives no new normalized comparison between m_n and m_N.

## 4. An elementary exponential lower tail for random-signing norm

Let B have law U_n. Expose its edges by fresh columns: at stage i expose
precisely B_ji for j<i, after all edges among [i−1] have been exposed. Set

    x_1 = 1,
    H_i = sum_{j<i} B_ji x_j,
    x_i = sign(H_i),      with sign(0)=1,
    G_i = x_i H_i = |H_i|.                                   (14)

Every x_j with j<i is measurable with respect to the preceding history.
Conditional on that history, the fresh B_ji are independent fair signs, so
the products B_ji x_j are independent fair signs. The conditional law of
H_i is therefore the law of a sum of i−1 independent Rademacher variables,
and it does not depend on the history. In particular, G_i is independent
of the preceding history. Induction proves that G_2,...,G_n are mutually
independent, although the states x_i themselves need not be independent.

The energy of the resulting state is exactly

    Q_B(x) = sum_{i=2}^n G_i.

Write mu_j=E|epsilon_1+...+epsilon_j|, with mu_0=0, and

    b_n = sum_{j=0}^{n−1} mu_j.                               (15)

Thus E Q_B(x)=b_n. This notation b_n is local to this proof and does not
denote the source matrix or any previously defined extremal sequence.

### The bounded-differences MGF, with its normalization

If F is a function of j independent fair signs and changing any one sign
changes F by at most 2, its Doob martingale reveals one input at a time.
Each centered martingale increment, conditional on the previous inputs,
lies in an interval of length at most 2. Hoeffding's elementary lemma gives
conditional MGF at most exp(t^2/2) for that increment. Iteration yields

    E exp(t(F−EF)) <= exp(t^2 j/2),       t in R.               (16)

For completeness, the lemma used here is: if E Z=0 and Z is supported on
an interval of length L, then E exp(tZ) <= exp(t^2 L^2/8). It follows by
bounding the exponential by its endpoint chord, followed by the elementary
two-point bound; equivalently, the second derivative of the log MGF of an
endpoint-valued variable is its tilted variance, at most L^2/4, and two
integrations from zero give the displayed inequality.

Apply (16) to F=|sum_{a=1}^j epsilon_a|, which is 2-Lipschitz in each
input. Independence of the G_i, or conditional iteration of the same bound,
then gives

    E exp(t(sum_i G_i − b_n)) <= exp(t^2 k/2),                 (17)

because sum_{j=0}^{n−1} j=k. Thus, for every s>=0, optimizing the negative
MGF parameter at t=s/k gives

    P(sum_i G_i <= b_n−s) <= exp(−s^2/(2k)).                  (18)

This uses the individual Rademacher inputs to obtain proxy j at stage j.
Using only the coarse range 0<=G_i<=i−1 would instead give a sum of squared
ranges and would lose the exponential-in-n rate needed below.

Since Phi(B)>=Q_B(x), for every real t<=b_n,

    P_(U_n)(Phi(B)<=t) <= exp(−(b_n−t)^2/(2k)).               (19)

For arbitrary real t, the valid unified upper bound is

    P_(U_n)(Phi(B)<=t)
        <= exp(−((b_n−t)_+)^2/(2k)),                         (20)

where r_+=max(r,0).

### The leading constant

The exact elementary Rademacher identities are

    mu_(2r)   = 2r binom(2r,r)/2^(2r),
    mu_(2r+1) = (2r+1) binom(2r,r)/2^(2r).                   (21)

For r>=1, the first identity follows by pairing binomial terms on the two
sides of zero and telescoping the first moment; the second also follows
from E(|S+epsilon| | S)=|S|+1_{S=0}. The case r=0 in the second identity
gives mu_1=1. Stirling's formula in (21) gives

    mu_j = sqrt(2j/pi) + O(j^(−1/2))          for j>=1.

Summing and comparing sum sqrt(j) with its integral yields

    b_n = c_g n^(3/2) + O(sqrt(n)),
    c_g = (2/3) sqrt(2/pi) > 1/2.                            (22)

The strict inequality is equivalent to pi<32/9 (for example pi<22/7
suffices). Numerically c_g=0.5319230405..., but the proof uses its exact
expression. No result about the precise random-signing ground-state limit
is invoked. Equations (19) and (22) imply, for each fixed epsilon>0,

    P_(U_n)(Phi(B)<=(c_g−epsilon)n^(3/2))
        <= exp(−(epsilon^2+o(1))n).                          (23)

## 5. Transfer to induced restrictions, with the TV error retained

For any fixed source A of order N, let S be a uniform n-set. Since Phi is
invariant under relabelling, (8), (20), and the ordered-injection definition
give the following finite statement for every real threshold t:

    P(Phi(A[S])<=t)
      <= exp(−((b_n−t)_+)^2/(2k)) + D_(N,n)(A).              (24)

The right side may of course be replaced by its minimum with 1. The TV
error in (24) is additive; one must control its rate, not just know that it
tends to zero.

Under the bounded-norm and order-separation assumptions (11),

    D_(N,n)(A_N)
      <= C k 2^k/sqrt(N) + k/N
       = exp(−omega(n)).                                    (25)

More explicitly, the logarithm of the first positive bound on the right is
at most −(1/2−o(1))log N, and log N/n -> infinity under (11). This verifies
the superexponential-in-n rate required when exponentially many event
bounds are accumulated.

Now suppose additionally that the source sequence satisfies

    limsup_N Phi(A_N)/N^(3/2) <= 1/2.                         (26)

Exact global optimizers satisfy (26) by the already proved universal
construction upper bound for m_N. That upper bound uses only the existence
of suitable signings at asymptotically dense orders; it does not assert
Paley optimality or identify the unknown limit.

Let eta_N be any real sequence tending to zero and consider the
source-normalized threshold

    t_(N,n) = (n/N)^(3/2) Phi(A_N) + eta_N n^(3/2).           (27)

Put delta_g=c_g−1/2>0. Equations (22), (26), and (27) give

    (b_n−t_(N,n))/n^(3/2) >= delta_g−o(1).

It follows from (24) and (25) that, for every fixed a with
0<a<delta_g^2, all sufficiently large N satisfy

    P(Phi(A_N[S])<=t_(N,n)) <= exp(−a n)
                                      + exp(−omega(n)).     (28)

Equivalently, the upper exponential rate is at most −delta_g^2. The less
precise bound exp(−Omega(n)) is already enough for the consequence below.

For comparison, (23)–(25) also show that a uniform restriction has

    Phi(A_N[S])/n^(3/2) >= c_g−o_P(1),
    liminf E[Phi(A_N[S])/n^(3/2)] >= c_g.                     (29)

For the expectation assertion, fix any epsilon>0 and use the high-probability
lower bound c_g−epsilon, together with nonnegativity of Phi; then let
epsilon decrease to zero. No transfer of an unbounded statistic by a
bare TV bound is being used.

## 6. Many restrictions: arbitrary dependence, uniform marginals

Fix the host A_N. Let S_1,...,S_T be random n-subsets on any common
probability space, each with the uniform marginal distribution on the
n-subsets of [N]. No independence between them is required; they may overlap
or be arbitrarily correlated. Applying (24) to each marginal and taking
a union bound gives the exact finite inequality

    P(exists j<=T: Phi(A_N[S_j])<=t)
      <= T [exp(−((b_n−t)_+)^2/(2k)) + D_(N,n)(A_N)].         (30)

If n satisfies (11), the source satisfies (26), t is (27), and

    log T = o(n),                                            (31)

then (25), (28), and (30) imply that this success probability tends to
zero. Indeed T exp(−a n)->0 and T D_(N,n)(A_N)->0 separately.

The hypothesis is uniform marginals for the fixed host, not merely that
the realized family has small cardinality. Inspecting the host and then
selecting its best restriction can produce a singleton family which is not
uniformly distributed; (30) does not apply to that selection. An adaptive
procedure qualifies only if the required marginal condition really holds
(or a separate suitable conditional bound is proved). If the host itself
is random, the displayed argument applies conditionally whenever each
candidate has the stated uniform marginal conditional on that host.

## 7. Exact scope of the result

For actual global optimizers A_N, with n=n(N) tending to infinity and
n^2=o(log N), all of the following are simultaneously true:

1. Every labelled n-vertex signing occurs as an induced submatrix.
2. The minimum norm over all n-vertex restrictions is exactly m_n.
3. A uniform restriction has typical normalized norm at least
   c_g=(2/3)sqrt(2/pi)>1/2.
4. Any exp(o(n)) candidate restrictions with uniform marginals fail, with
   probability tending to one, to achieve the source-normalized threshold
   (27), despite arbitrary dependence between those candidates.

Universality proves existence of m_n-optimal restrictions. It does NOT prove
existence of restrictions meeting the source-normalized threshold (27):
alpha_n may exceed Phi(A_N)/N^(3/2) + eta_N. In fact, under universality,
such a threshold-(27) restriction exists if and only if

    alpha_n <= Phi(A_N)/N^(3/2) + eta_N.                       (32)

For a globally minimizing source, the right side is alpha_N + eta_N.
This cross-order comparison is not proved here; the threshold-(27) event
may even be empty. The proved existence statement concerns m_n-optimal
restrictions only, and must not be promoted to source-constant preservation.

What fails is typical, mean, or modest-many uniform sampling as a way to
transfer the source's optimal leading constant. Biased selection, adaptive
methods without uniform marginals, all-subset arguments, and comparisons
at other order ratios are not excluded. No conclusion about existence,
nonexistence, or value of lim alpha_n follows from this result.
