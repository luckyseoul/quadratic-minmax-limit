# A deterministic selected-half restriction inequality

2026-09-05. Analytic theorem, no computation. This gives a genuine but
non-sharp cross-order inequality. It neither proves the requested sharp
proportional paving factor nor proves convergence of the original MO
sequence. No hypothesis of global optimality is needed for the restriction
theorem; applying it to an actual global minimizer supplies the order
comparison without identifying its restrictions as smaller minimizers.

## 1. Definitions and disjoint-phase inequality

For any real symmetric zero-diagonal matrix A on a finite vertex set V,
put Q_A(x)=sum_(i<j) A_ij x_i x_j and Phi(A)=max_x |Q_A(x)|.
For S contained in V put

    p(S)=max_x Q_(A[S])(x),
    q(S)=-min_x Q_(A[S])(x),
    Phi(A[S])=max(p(S),q(S)).

Both p(S) and q(S) are nonnegative, because the uniform spin mean is zero.
If S and T are disjoint, then

    p(S)+p(T) <= p(V) <= Phi(A),
    q(S)+q(T) <= q(V) <= Phi(A).                         (1)

To prove the first inequality, choose separately maximizing assignments
on S and T. Multiply each entire block assignment by an independent fair
sign and assign independent fair spins to unused vertices. Internal
energies remain fixed, whereas the expected energy of every other edge
is zero. The expected full energy is p(S)+p(T), which is at most p(V).
Choose minimizing assignments to obtain the second inequality.

## 2. Exact odd-order selected restriction

For every n>=1 and every real A of order 2n+1,

    min_(|S|=n) Phi(A[S]) <= Phi(A)/2.                  (2)

Label the vertices by the integers modulo 2n+1, and consider the sets

    S_i={i,i+2,...,i+2(n-1)},   i=0,...,2n.

Consecutive sets S_i,S_(i+1), including S_(2n),S_0, are disjoint.
Indeed an intersection would give 2(a-b)=1 modulo 2n+1 with
|a-b|<=n-1, which is impossible. Thus these sets form an explicitly
specified odd cycle of disjoint pairs.

Write M=Phi(A). If every S_i had norm greater than M/2, color S_i
positive if p(S_i)>M/2 and negative otherwise. A negatively colored
set has q(S_i)>M/2. Equation (1) prohibits two consecutive sets with
the same color. But a cycle of odd length has no such two-coloring:
alternation along all its edges would return the opposite starting
color. This contradiction proves (2). If both phases exceed M/2,
choosing the positive color still makes the argument valid.

The proof selects from only 2n+1 explicitly listed n-sets after seeing
their norms. It is not a uniform-marginal sampling statement. Its threshold
is half the full norm, not the source-normalized threshold from the
previous tiny-n restriction theorem. There is no conflict with that theorem.

The same conclusion holds for every order N>=2n+1: restrict first to
any 2n+1 vertices, whose norm is at most the full norm, then use (2).

## 3. Even-order selected restriction with an explicit boundary error

Suppose now that A has order 2n and |A_ij|<=a. Then

    min_(|S|=n) Phi(A[S]) <= [Phi(A)+a(n-1)]/2.           (3)

For n=1 the assertion is immediate. For general n set L=a(n-1).
If n-sets S,T differ by replacing exactly one vertex, their common
(n-1)-set R obeys

    p(R) <= p(S),p(T) <= p(R)+L,

and consequently

    |p(S)-p(T)|<=L.                                    (4)

The lower bound follows by averaging over the added spin; the upper
bound follows because its row contributes at most a(n-1). The same
estimate also holds for q, although only (4) is needed below.

Write M=Phi(A) and t=(M+L)/2. Suppose every n-set has norm greater
than t. Color a set positive if p(S)>t, and negative otherwise, so
q(S)>t in the latter case. By (1), a set and its complement must
have opposite colors. (They cannot both be positive or both negative,
since 2t>=M.) In particular, the complement of a negative set T is
positive, so

    p(T) <= M-p(T^c) < M-t.

If a positive set S and a negative set T differed by one vertex, then

    p(S)-p(T) > t-(M-t)=L,

contradicting (4). Thus adjacent n-sets cannot have different colors.
Any n-set can be changed into its complement by successively replacing
one vertex, so the opposite colors at those endpoints yield a
contradiction. This proves (3).

No rounding of a dense coefficient block and no control of a selected
Gibbs covariance is used in this argument.

## 4. Consequences for actual global norm minima

For complete signings define m_N=min_A Phi(A) and alpha_N=m_N/N^(3/2).
Apply (2) and (3) to an ACTUAL global minimizer at the larger order. Each
selected restriction is admissible at order n, and therefore

    m_(2n+1) >= 2 m_n,
    m_(2n)   >= 2 m_n-(n-1).                            (5)

The selected restriction need not itself minimize at order n. Its
nonnegative gap above m_n is retained implicitly in the admissibility
inequality m_n<=Phi(A[S]); discarding that gap only weakens the bound
in the valid direction.

At balanced even orders (5) gives

    alpha_n <= sqrt(2) alpha_(2n)+(n-1)/(2n^(3/2)).       (6)

This improves the factor 2^(3/2) from ordinary restriction monotonicity
to sqrt(2), but the latter remains a fixed factor greater than one.
For a norm-minimizing source of order 2n, the sharp target would instead
be

    min_(|S|=n) Phi(A[S]) <= 2^(-3/2) m_(2n)+o(n^(3/2)). (7)

Equation (3) has leading factor 1/2, not 2^(-3/2). The difference is
on the leading n^(3/2) scale, so the linear boundary error is not the
live issue. Neither (6) nor its iteration proves alpha convergence.

## 5. Scope and duplication audit

The current induced-optimizer note proves a complete tiny-n signing law,
universality, and failure of uniform-marginal sampling at a different,
source-normalized threshold. Those arguments are not used or repeated.
The sharp-influence dimension audit establishes ordinary restriction
monotonicity and explains the wrong averaging direction; (2)--(3) use
a deterministic phase coloring instead of reversing that inequality.

The older arbitrary fixed block partition bound in solution.md gives
k m_n<=2 m_(kn) by adding positive and negative block extrema. At k=2
that is only m_n<=m_(2n). Equations (2)--(5) instead exploit the full
family of induced subsets and its disjointness/exchange structure.

No proof or refutation of the sharp proportional restriction statement
under the exact-global-minimizer hypothesis is supplied. There is no
new assertion about norm optimality of conference constructions and no
asymptotic conclusion based on a fixed-order counterexample.
