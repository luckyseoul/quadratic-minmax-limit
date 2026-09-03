# Random vertex orders: exact variance and the missing minimizer statistic

**Status:** proved all-orders averaging identities and an open reduction.  A
random vertex order is not shown to satisfy the multiplier-two bound.  The
calculation identifies the exact joint statistic that is absent from
`Phi(A)=M` and from the immediate variational consequence of global
minimality.  No finite signing or order census is used.

Let `A` be a symmetric signing and `M=Phi(A)`.  Regard `pi(i)` as the
position of vertex `i` in the order, and define its transitive tournament by
`(T_pi)_ij=+1` exactly when `pi(i)<pi(j)`.  Put

\[
 R^\pi=A\circ T_\pi,
 \qquad
 \Lambda_\pi(A)={1\over2}B(A,R^\pi)
 =\max_U\Phi\left(A^{E^+_\pi(U)}\right).            \tag{1}
\]

The last equality is the directed half-cut identity.  The desired ordered
construction would be

\[
 \min_\pi\Lambda_\pi(A)\le\sqrt2M+O(n).             \tag{2}
\]

## 1. Exact average of the prefix variation

Fix a Boolean state `y` and write

\[
 \ell_i=y_i(Ay)_i=\sum_{j\ne i}A_{ij}y_iy_j.
\]

For a uniformly random order, the local field when vertex `i` is crossed in
the prefix path is

\[
 L_i^\pi=\sum_{j\ne i}A_{ij}y_iy_j
          \operatorname{sgn}(\pi(j)-\pi(i)).         \tag{3}
\]

Among three uniformly ordered vertices, two comparisons with the same
endpoint have product mean `1/3`.  Hence

\[
 \mathbb E_\pi L_i^\pi=0,
 \qquad
 \mathbb E_\pi(L_i^\pi)^2
 ={2(n-1)\over3}+{\ell_i^2\over3}.                  \tag{4}
\]

If `q_k` is the `A`-energy after flipping the first `k` vertices of the
order, then (q_{k-1}-q_k=2L^pi_{v_k}).  Thus the prefix total-variation
identity, with its normalization explicit, is

\[
 \max_x|x^TR^\pi y|=\sum_i|L_i^\pi|
 ={1\over2}\sum_{k=1}^n|q_k-q_{k-1}|.              \tag{5}
\]

Consequently

\[
 \mathbb E_\pi\max_x|x^TR^\pi y|
 \le\sum_i\sqrt{{2(n-1)\over3}+{\ell_i^2\over3}}.  \tag{6}
\]

The degree-two best-response inequality gives

\[
 \sum_i|\ell_i|
 \le2M+2\sqrt{M^2-Q_A(y)^2},                        \tag{7}
\]

so a completely scalar consequence of (6) is

\[
 \mathbb E_\pi\max_x|x^TR^\pi y|
 \le n\sqrt{{2(n-1)\over3}}
 +{2\over\sqrt3}\left(M+\sqrt{M^2-Q_A(y)^2}\right).\tag{8}
\]

This is a `for every y, there exists an order depending on y` estimate.  It
does not permit the quantifier reversal needed in (2).  In fact the first
average over `y` contains no order-selection information at all.  For every
skew signing `R`, each `(Ry)_i` is a sum of `n-1` independent signs under
uniform `y`; therefore

\[
 \boxed{
 \mathbb E_y\max_x|x^TRy|
 =\mathbb E_y\|Ry\|_1=n\mu_{n-1}.
 }                                                     \tag{9}
\]

Thus every order has exactly the same mean prefix variation.  Only its
upper tail and its correlation with the two endpoint energies can matter.

## 2. Exact fixed-event permutation variance

Fix a nontrivial cut `U,V` of sizes `k,l`, a Boolean `y`, and define the
signed cross matrix

\[
 W_{ij}=A_{ij}y_i y_j\quad(i\in U,j\in V),
 \qquad h=kl.
\]

If `F_pi,G_pi` are the cross-edge sums directed respectively from `U` to
`V` and from `V` to `U`, put

\[
 Z_\pi=G_\pi-F_\pi.
\]

Then the corresponding directed half-cut value is

\[
 Q_{A^{E^+_\pi(U)}}(y)=I(U,y)+Z_\pi,                \tag{10}
\]

where `I` is the noncrossing energy.  Let

\[
 r=W\mathbf1,qquad c=W^T\mathbf1,qquad
 D_A(U,y)=\|r\|_2^2+\|c\|_2^2.                     \tag{11}
\]

Comparisons on disjoint edges have covariance zero; comparisons sharing
their `U` endpoint or their `V` endpoint have covariance `1/3`.  Expanding
the square therefore gives the exact formula

\[
 \boxed{
 \mathbb E_\pi Z_\pi=0,qquad
 \operatorname{Var}_\pi Z_\pi={h+D_A(U,y)\over3}.
 }                                                     \tag{12}
\]

There is an equivalent Hoeffding decomposition.  Generate the order using
independent uniform ranks `s_i in [0,1]`.  With `i in U,j in V`, the sign
appearing in `Z_pi` is `sgn(s_i-s_j)`, and

\[
 \operatorname{sgn}(s_i-s_j)
 =(2s_i-1)+(1-2s_j)+\kappa_0(s_i,s_j),              \tag{13}
\]

where `kappa_0` has zero conditional mean in either variable, the first two
terms contribute `D_A/3` to (12), while the degenerate kernel contributes
the irreducible `h/3`.

The `Phi` bound controls `D_A`, but not at the scale needed for an
exponential union bound.  For all Boolean `a,b`, polarization gives

\[
 |a^TWb|\le M.                                       \tag{14}
\]

In particular `||W 1||_1<=M` and `||W^T 1||_1<=M`.
Since individual row and column sums have magnitudes at most `l` and `k`,
respectively,

\[
 \boxed{D_A(U,y)\le lM+kM=nM,}
 \qquad
 \operatorname{Var}_\pi Z_\pi\le{h+nM\over3}.       \tag{15}
\]

At `M=Theta(n^(3/2))`, (15) permits variance of order `n^(5/2)`, not
`n^2`.  A deviation of order `n^(3/2)` then has only a `sqrt(n)` quadratic
variance exponent, whereas (2) must control exponentially many `(U,y)`
events.

The degree term is genuinely present on average.  For every fixed cut,
orthogonality of the Boolean characters gives

\[
 \mathbb E_yD_A(U,y)=2h,
 \qquad
 \mathbb E_{y,\pi}Z_\pi^2=h.                        \tag{16}
\]

Thus a typical event has the same second moment as the independent-edge
orientation model.  Even the fictitious perfectly balanced case `D_A=0`
has variance `h/3`.  If it also had an ideal Gaussian tail and `I=0`, the
central-cut entropy calculation would require

\[
 12\alpha^2>\log4,qquad
 \alpha>{\sqrt{\log2/6}}=0.3398889967\ldots,         \tag{17}
\]

still above the universal lower endpoint `1/pi=0.3183098861...`.  Equation
(17) is a barrier for a purely quadratic-variance Chernoff argument, not a
lower-tail theorem for every ordering statistic.

## 3. Exact random-order criterion and its missing input

Let `OrdTail_W(t)` denote

\[
 \Pr_\pi\{|Z_\pi|>t\}.
\]

The same eightfold orbit quotient used for random independent tournaments
shows that a sufficient condition for an order satisfying (2), without its
`O(n)` allowance, is

\[
 \boxed{
 {1\over8}\sum_{(U,y)}
 \operatorname{OrdTail}_{W(U,y)}
 \left(\sqrt2M-|I(U,y)|\right)<1.
 }                                                     \tag{18}
\]

Unlike the independent-tournament tail, `OrdTail_W` is not determined by
`h`; already its exact variance depends on (11).  The first missing
minimizer statistic is therefore the hereditary switched cross-degree
energy

\[
 \mathcal D(A)=\max_{U,y}
 \left(\|A_{U,V}y_V\|_2^2+
       \|A_{U,V}^Ty_U\|_2^2\right).                 \tag{19}
\]

The scalar norm gives only `mathcal D(A)<=nM`.  A bound
`mathcal D(A)=O(n^2)` would put the linear Hoeffding component on the
necessary variance scale, but would still have to be paired with a
speed-`n` tail bound for the degenerate component in (13) and with the joint
energy enumerator of `I(U,y)`.  Equivalently, the genuinely sufficient
missing object is a uniform joint large-deviation estimate for

\[
 \bigl(I(U,y),D_A(U,y),Z_\pi\bigr),                 \tag{20}
\]

strong enough that the sum in (18), not merely each fixed event, is small.

Spectral flatness would control the first new statistic:

\[
 \mathcal D(A)\le n\|A\|_{op}^2.                   \tag{21}
\]

Thus `||A||op=O(sqrt(n))` is a concrete sufficient input for the variance
scale, though not by itself for (18).  It is not a consequence of the known
cube norm estimate.  A direct randomized-rounding argument gives only

\[
 \|A\|_{op}^3\le2nM,
 \qquad
 \|A\|_{op}=O(n^{5/6})                              \tag{22}
\]

at the relevant scale.  For completeness, if `v` is a unit extremal
eigenvector and `a=||v||_infty`, independent signs with means `v_i/a` give
`2M>=||A||op/a^2`; the eigenvector equation gives
`a<=sqrt(n)/||A||op`, proving (22).

## 4. A norm-scale countermechanism

The gap in (15) is not an artifact of the estimate.  Fix
`r=floor(c sqrt(n))` vertices and set every edge of their rectangular cut to
`+1`.  Fill all remaining edges with independent signs.  For any fixed
Boolean state the random part is a sum of at most `N` independent signs, so
a union bound shows that, for every fixed `C>sqrt(log 2)`, there are fillings
with

\[
 \Phi(A)\le(c+C+o(1))n^{3/2}.                       \tag{23}
\]

Indeed the deterministic rectangle contributes at most
`r(n-r)=(c+o(1))n^(3/2)`, while Hoeffding gives failure probability at most
`2 exp((log 2-C^2+o(1))n)` for the random remainder.

For the rectangular cut and `y=1`, however,

\[
 D_A(U,1)=r(n-r)^2+(n-r)r^2
          =r(n-r)n=\Theta(n^{5/2}),                 \tag{24}
\]

and (12) attains that order.  Moreover the event that all `r` vertices
precede all other vertices has probability

\[
 {1\over{n\choose r}}=\exp(-O(\sqrt n\log n)),       \tag{25}
\]

and on it `|Z_pi|=r(n-r)`.  Hence no speed-`n` ordering tail bound follows
from the correct `n^(3/2)` cube-norm scale alone.

This construction is not a global minimizer and is not a counterexample to
(2).  Its role is exact: it proves that any successful use of (18) must
extract new structure from **global** minimality, rather than use only its
numerical value `M=m_n`.

## 5. What global minimality currently says

If `A` is globally minimizing, every directed half-cut neighbor is another
signing of order `n`, so the immediate variational statement is

\[
 \Phi(A^{E^+_\pi(U)})\ge m_n=M.                     \tag{26}
\]

This is a lower bound; (2) needs a uniform upper bound.  Averaging (26) over
orders preserves the wrong direction, and none of (4), (9), (12), or (16)
changes when `A` is minimizing.  No bound on (19) or (20) follows from the
definition alone.

Therefore the ordered construction remains an honest open reduction.  The
specific next lemma is not another choice of order or a finite test: prove
that global minimizers have a speed-`n` joint ordering-tail estimate of the
form (18), with (19) as its first necessary structural statistic.  Without
such a minimizer theorem, averaging over vertex orders does not yield
`Lambda_pi<=sqrt(2)M+O(n)`.
