# Cut-code, spherical, and sharp-influence audit of the original MO limit

Date: 2026-09-02

Status: **rigorous reformulations and dimension-comparison obstructions; this
does not prove the limit.**  No finite census, Paley multiplier computation, or
residual-(ii) search is used here.

Put

\[
 Q_A(x)=\sum_{1\leq i<j\leq n}a_{ij}x_ix_j,\qquad
 \Phi(A)=\max_{x\in\{\pm1\}^n}|Q_A(x)|,
 \qquad m_n=\min_A\Phi(A),
\]

where every `a_ij` is in `{-1,1}`, and put `N=binom(n,2)`.

## 1. Exact sharp-\(L^1\)-influence reformulation

For a function on the cube use

\[
 \partial_i f(x)={f(x)-f(x^{\oplus i})\over2},\qquad
 I_1(f)=\sum_i\mathbb E|\partial_i f|,
\]

and set

\[
 \mu_r=\mathbb E|\epsilon_1+\cdots+\epsilon_r|,
 \quad
 K_n=\max_{A\in\{\pm1\}^{E(K_n)}}
 I_1\!\left({Q_A\over\Phi(A)}\right).
\]

### Proposition 1 (exact identity)

For every `n>=2`,

\[
 \boxed{\quad K_n={n\mu_{n-1}\over m_n}.\quad}                 \tag{1}
\]

Consequently,

\[
 {m_n\over n^{3/2}}={\mu_{n-1}/\sqrt n\over K_n},
 \qquad {\mu_{n-1}\over\sqrt n}\longrightarrow\sqrt{2/\pi}.
                                                               \tag{2}
\]

**Proof.**  Direct differentiation gives

\[
 \partial_iQ_A(x)=x_i\sum_{j\ne i}a_{ij}x_j.
\]

Under uniform `x`, the summands `a_ij x_j` are independent Rademacher
variables.  Hence `||partial_i Q_A||_1=mu_(n-1)`, independently of `A`, and
`I_1(Q_A)=n mu_(n-1)`.  Maximizing after normalization by `Phi(A)` is therefore
exactly the same operation as minimizing `Phi(A)`, proving (1).  Finally

\[
 \mu_r={r\over2^{r-1}}
 \binom{r-1}{\lfloor(r-1)/2\rfloor},
\]

and Stirling's formula proves (2).  QED.

The word "consequently" can be strengthened to a genuine equivalence, not
merely a one-way implication:

### Corollary 2 (the limit is exactly a sharp-constant limit)

The finite limit of `m_n/n^(3/2)` exists if and only if the finite limit of
`K_n` exists.  If either exists, then

\[
 \lim {m_n\over n^{3/2}}
 ={\sqrt{2/\pi}\over\lim K_n}.                                \tag{3}
\]

**Positivity audit.**  This reciprocal statement cannot silently fail at zero.
The degree-`d` theorem of Filmus--Hatami--Keller--Lifshitz gives
`I_1(f)<=d^2 ||f||_infty`, hence `K_n<=4`.  Conversely, choose all coefficients
of `A` independently.  Hoeffding plus a union bound over the cube shows that
some `A` satisfies

\[
 \Phi(A)\leq\sqrt{2N(n+2)\log2}.
\]

Khintchine's inequality `mu_r>=sqrt(r/2)` then gives, uniformly for `n>=2`,

\[
 K_n\geq\sqrt{{n\over2(n+2)\log2}}>0.
\]

Thus `K_n` stays in a fixed compact subinterval of `(0,infinity)`, and taking
reciprocals in (2) is legitimate in both directions.

Combining (2) with the canonical sandwich already proved in `solution.md`,

\[
 {1\over\pi}\leq\liminf {m_n\over n^{3/2}}
 \leq\limsup {m_n\over n^{3/2}}\leq{1\over2},
\]

places the sharp constants in the narrower interval

\[
 2\sqrt{2/\pi}\leq\liminf K_n
 \leq\limsup K_n\leq\sqrt{2\pi}.                             \tag{3a}
\]

In particular, the repository's candidate limit `1/2` is exactly the candidate
sharp-influence limit `2 sqrt(2/pi)`.

This is a new sharp target, but it is not supplied by the existing influence
literature.  That literature bounds the supremum over all bounded degree-two
functions uniformly in the dimension.  It gives no comparison between the
much narrower, equimodular, complete-support constants `K_n` and `K_(n+1)`.

## 2. Adding one vertex: exact interpolation identity and obstruction

Write a signing of `K_(n+1)` as an old signing `A` and a new-edge vector
`b in {+-1}^n`.  Since

\[
 Q_{A,b}(x,t)=Q_A(x)+t\sum_i b_ix_i
\]

and `max_(t=+-1)|u+tv|=|u|+|v|`, one has the exact identity

\[
 \boxed{\quad
 m_{n+1}=\min_{A,b}\max_x
 \left(|Q_A(x)|+\left|\sum_i b_ix_i\right|\right).
 \quad}                                                       \tag{4}
\]

In particular,

\[
 m_n\leq m_{n+1}\leq m_n+n.                                  \tag{5}
\]

The first inequality also follows by averaging an `(n+1)`-variable form over
the last variable.  If `n` is odd then the linear form in (4) is never zero,
so `m_(n+1)>=m_n+1`.  Equivalently, if
`Delta_n=m_(n+1)-m_n`, then

\[
 0\leq\Delta_n\leq n,
 \qquad \Delta_n\equiv n\pmod2.                              \tag{6}
\]

For the congruence, every value of an `n`-vertex form is a sum of `N` odd
integers and hence `m_n` has the parity of `N`; subtract the two consecutive
parities.

The exact Rademacher-mean formula gives

\[
 r_n:={(n+1)\mu_n\over n\mu_{n-1}}
 =\begin{cases}
  (n+1)/n,&n\text{ even},\\
  (n+1)/(n-1),&n\text{ odd}.
 \end{cases}
\]

Therefore

\[
 \boxed{\quad {K_{n+1}\over K_n}
 ={r_n\over1+\Delta_n/m_n}.\quad}                            \tag{7}
\]

This locates the missing information exactly.  For even `n`, monotonicity of
`K_n` changes direction according as `Delta_n` lies below or above `m_n/n`;
for odd `n`, the threshold is `2m_n/(n-1)`.  These are order-`sqrt(n)`
increment questions, while (5) has an order-`n` error.  Thus vertex padding
does not yield a summable interpolation error.

There is no hidden monotonicity of `K_n` even before asymptotics.  This can be
seen without a census: `m_2=1`, `m_3=3`, and `m_4=4`.  The middle identity
follows because the three edge products have product one.  For `n=4`,
Parseval and parity give `m_4>=4`, while

\[
 Q=x_1(x_2+x_3+x_4)+x_2x_3-x_2x_4-x_3x_4
\]

has norm at most four: if `x_2=x_3=x_4`, the two displayed absolute
contributions are `3` and `1`; otherwise they are at most `1` and `3`.
Since `mu_1=mu_2=1` and `mu_3=3/2`,

\[
 (K_2,K_3,K_4)=(2,1,3/2),
\]

which rules out monotonicity in either direction.

## 3. Two blocks: exact decomposition, but no tensorization

For positive `r,s` define the rectangular switching value

\[
 b_{r,s}=\min_{C\in\{\pm1\}^{r\times s}}
 \max_{x\in\{\pm1\}^r,y\in\{\pm1\}^s}|x^TCy|.
\]

Splitting a complete signing into its two diagonal blocks and its rectangular
cross block proves

\[
 \boxed{\quad
 \max\{m_r,m_s,b_{r,s}\}\leq m_{r+s}
 \leq m_r+m_s+b_{r,s}.
 \quad}                                                       \tag{8}
\]

For the lower bound, conditional expectation over the other block extracts
each diagonal quadratic form.  Also

\[
 x^TCy={Q(x,y)-Q(x,-y)\over2},
\]

so the cross-block norm is at most the full norm.  For the upper bound, choose
three minimizers separately and use the triangle inequality.

The quantity `b_(r,s)` is precisely the discrepancy in the rectangular
Gale--Berlekamp switching game.  Known results give its correct order
`sqrt(rs(r+s))`, but the term is of the same leading order as `m_(r+s)`.
Thus (8) is not a Fekete inequality.  In the equal split, even a sharp bound
of the form `m_(2n)<=2m_n+b_(n,n)` still requires a leading-constant estimate
on the cross term; treating it as an error cannot prove convergence.

This also explains the coding obstruction.  Covering radii add exactly for an
abstract direct sum of binary codes:

\[
 \rho(D_1\oplus D_2)=\rho(D_1)+\rho(D_2).
\]

Consequently their central deficits add.  But a quantity scaling as `n^(3/2)`
must multiply by `2sqrt(2)`, not by `2`, under an equal doubling.  The missing
leading contribution is exactly the complete bipartite cross block.  The
augmented cut code of `K_(2n)` is not the direct sum of two augmented cut
codes, so direct-sum covering-radius theory cannot supply it.

## 4. Exact cut-code and spherical formulations

For a projective Boolean vector `[x] in {+-1}^n/{+-1}`, define

\[
 c_x=(x_ix_j)_{i<j}\in\{\pm1\}^N.
\]

These vectors obey

\[
 \langle c_x,c_y\rangle={(x\cdot y)^2-n\over2},              \tag{9}
\]

and form an isotropic tight frame:

\[
 2^{-(n-1)}\sum_{[x]}c_xc_x^T=I_N.                           \tag{10}
\]

Let `C_n` be the binary cut space of `K_n`, and let
`D_n=C_n+span{1}` be its augmentation by the all-one word.  Translating signs
to binary words gives the exact covering-radius identity

\[
 \boxed{\quad m_n=N-2\rho(D_n).\quad}                         \tag{11}
\]

Indeed, for a sign word `a=(-1)^z`, maximizing the absolute correlation with
`c_x` is `N-2 d(z,D_n)`; minimizing that maximum over `z` maximizes the distance
to `D_n`.

The augmentation is essential.  The classical maximum-frustration theorem
computes the covering radius of the *ordinary* cut code `C_n` (equivalently a
one-sided switching objective).  It does not compute `rho(D_n)`: adjoining
the all-one word identifies every switched signing with its global negative
and enforces the absolute value in the MO objective.  Confusing these two
codes changes the central deficit from order `n` to the unresolved order
`n^(3/2)` quantity.

Equation (10) alone only yields

\[
 \max_x|\langle a,c_x\rangle|\geq\sqrt N,
\]

because the mean squared correlation is `N`.  This is order `n`, a factor
`sqrt(n)` below the order relevant to the MO problem.  The missing strength
comes from restricting the ambient normal `a` itself to a hypercube vertex;
ordinary spherical covering or tight-frame estimates discard that arithmetic
restriction.

## 5. Exact fixed-density cut-discrepancy equivalence

Define

\[
 H_n=\min_{\substack{G\text{ on }[n]\\
                     e(G)=\lfloor N/2\rfloor}}
 \max_{S\subseteq[n]}
 \left|e_G(S,S^c)-\tfrac12|S||S^c|\right|.
\]

Then, for every `n>=2`,

\[
 \boxed{\quad m_n-1\leq4H_n\leq m_n+\sqrt N.\quad}           \tag{12}
\]

For the first inequality, let `a_e=2 1_(e in G)-1`.  Its total
`T=2e(G)-N` is zero or minus one, and for the sign vector associated with `S`,

\[
 Q_A(x_S)=T-4\left(e_G(S,S^c)-\tfrac12|S||S^c|\right).
\]

Thus `m_n<=4H_n+1`.

Conversely, take an optimal `A`.  Walsh orthogonality gives
`E Q_A^2=N`, so some `x` has `|Q_A(x)|<=sqrt(N)`.  Switch by this `x`; the new
total coefficient sum `T` has `|T|<=sqrt(N)` and the same norm.  If `N` is
even, flip `|T|/2` majority edges to make the total zero.  If `N` is odd,
globally negate first so `T<0`, then flip `(|T|-1)/2` negative edges to make
the total minus one.  Each edge flip changes the sup norm by at most two.  The
resulting half-density graph has discrepancy at most `(m_n+|T|)/4`, proving
the other half of (12).

Hence the original limit exists if and only if `H_n/n^(3/2)` has a limit, and
the two limits differ by the factor four.

## 6. Why ordinary graphons erase the question

Regard a signing as a symmetric step kernel.  For subsets `S,T`, polarization
and multilinearity on `[-1,1]^n` bound its rectangular sum by a fixed multiple
of `Phi(A)`.  Explicitly, for `u=1_S`, `v=1_T`,

\[
 u^TAv={1\over2}\{Q_A(u+v)-Q_A(u-v)\},
\]

while `|Q_A((u+v)/2)|<=Phi(A)` and `|Q_A(u-v)|<=Phi(A)`.
Thus the normalized cut norm is `O(Phi(A)/n^2)`.  Since signings with
`Phi(A)=O(n^(3/2))` exist, every optimizing sequence converges in ordinary
dense cut metric to the constant half-density graphon, at rate
`O(n^(-1/2))`.

Therefore any functional continuous in the ordinary Lovasz--Szegedy cut
topology has the same limit on every candidate optimizer and cannot see the
`n^(3/2)` coefficient.  A graph-limit proof would need a second-order
fluctuation object, not an ordinary graphon.

## 7. Literature boundary

The following primary sources were checked against the exact targets above:

1. Y. Filmus, H. Hatami, N. Keller, and N. Lifshitz, *On the sum of the
   L1 influences of bounded functions*, Israel J. Math. 214 (2016), 167--192,
   DOI `10.1007/s11856-016-1355-0`, arXiv:1404.3396.  Its dimension-free
   `I_1(f)<=d^2` theorem supplies compactness of the numerical constants, not
   comparison or convergence of the restricted `K_n`.
2. A. Backurs and M. Bavarian, *On the sum of L1 influences*, CCC 2014,
   132--143, arXiv:1302.4625.  Its graph application and homogeneous bounds
   likewise control scale only.
3. P. Sole and T. Zaslavsky, *A Coding Approach to Signed Graphs*, SIAM J.
   Discrete Math. 7 (1994), 544--553, DOI `10.1137/S0895480189174374`.
   This identifies switching/frustration with cut-code covering radius and
   gives covering bounds, but no second-order limit theorem for the augmented
   cut codes of complete graphs.
4. T. A. Brown and J. H. Spencer, *Minimization of +/-1 matrices under line
   shifts*, Colloq. Math. 23 (1971), 165--171, and D. Pellegrino and
   A. Raposo Jr., *Upper bounds for the constants of Bennett's inequality and
   the Gale--Berlekamp switching game*, arXiv:2111.00445.  These control the
   rectangular cross-block problem in (8), but do not turn it into a lower
   order error or provide the needed complete-graph composition theorem.
5. L. Lovasz and B. Szegedy, *Limits of dense graph sequences*, J. Combin.
   Theory Ser. B 96 (2006), 933--957, DOI `10.1016/j.jctb.2006.05.002`.
   Its topology operates at the `n^2` scale and therefore collapses the
   second-order distinction described in Section 6.
6. A. Volberg, *An estimate of Sidon constant for complex polynomials with
   unimodular coefficients*, arXiv:2205.04936, and A. Defant,
   M. Mastylo, A. Perez, *On the Fourier spectrum of functions on Boolean
   cubes*, Math. Ann. 374 (2019), 653--680, arXiv:1706.03670.  The relevant
   Sidon/Bohnenblust--Hille estimates recover order of growth, not an
   `n -> n+1` or block interpolation for the flat real quadratic class.

**Audit conclusion.**  No theorem in these sources implies convergence of
`K_n`, `H_n/n^(3/2)`, or the covering deficit in (11).  The exact fresh target
is now (7): control the optimal one-vertex increment on the `sqrt(n)` scale,
or, equivalently in a two-half attack, control the cross block in (8) at its
leading constant.  Dimension-free influence bounds, direct-sum covering
radii, tight-frame geometry, and ordinary graphon compactness do not do this.
