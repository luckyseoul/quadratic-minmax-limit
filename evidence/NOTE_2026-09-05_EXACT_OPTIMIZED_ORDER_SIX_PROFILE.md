# Exact optimized balanced profile at order six

2026-09-05. **Analytic fixed-order theorem; original convergence OPEN.**

This note determines an actual global minimum over all complete order-six
signings for every point of the balanced two-block profile. Its proof does
not require a signing census, floating-point optimization, or a grid in
temperature. The finite integer catalog and GPU profile evaluation are
independent regressions, not premises of the theorem.

The conclusion is specific to order six. It does not compare arbitrarily
large orders, establish a fixed-temperature small-oh error in the order,
or prove convergence of the original normalized optimum.

## 1. Definitions and statement

Let \(\mathcal S_6\) denote the symmetric \(6\times6\) matrices with zero
diagonal and off-diagonal entries in \(\{-1,1\}\). Partition the six
vertices into two labelled sets of size three. For independent uniform
spin vectors \(x,y\in\{-1,1\}^3\), write

\[
 A=\begin{pmatrix}D_L&B\\B^T&D_R\end{pmatrix},\qquad
 I_A(x,y)=Q_{D_L}(x)+Q_{D_R}(y),\qquad C_A(x,y)=x^TBy,
\]

where \(Q_D(x)=\sum_{i<j}D_{ij}x_ix_j\). Thus \(I_A\) is the sum of
six signed internal edge characters, while \(C_A\) is the sum of nine
signed cross-edge characters. Define

\[
 Z_A(u,v)=\mathbb E_{x,y}\cosh(uI_A+vC_A),\qquad
 Z_*(u,v)=\min_{A\in\mathcal S_6}Z_A(u,v).
\]

All expectations in this note use all \(64\) equally weighted spin
states. Since logarithm is increasing, \(\log Z_*\) is also the minimum
of the individual log pressures.

**Theorem.** For every \(u\ge v\ge0\), putting
\(X=\cosh(2u)\) and \(Y=\cosh(2v)\),

\[
 \boxed{Z_*(u,v)=\frac{\cosh v}{4}
       \bigl(3X^2+3Y^2+2Y-4\bigr).}                       \tag{1}
\]

If \(v>0\), a signing minimizes precisely when \(A^2=5I_6\): the
minimizers are exactly the order-six symmetric conference signings.
There are twelve minimizers after normalizing every entry of the first
row off the diagonal to \(+1\) by vertex switching.

If \(v=0<u\), the minimizers are precisely the signings whose two
internal triangle edge-products have opposite signs; the cross block
is arbitrary. If \(u=v=0\), every signing minimizes.

The balanced path of the preceding variational note is

\[
 u=c\sqrt{\frac{2-t}{6}},\qquad
 v=c\sqrt{\frac t6},\qquad c\ge0,\quad0\le t\le1,
 \qquad f_6(c,t)=\log Z_*(u,v).                             \tag{2}
\]

For \(c>0\), this parameterization covers exactly \(u\ge v\ge0\),
apart from the origin: its inverse is
\(c^2=3(u^2+v^2)\), \(t=2v^2/(u^2+v^2)\).

## 2. Exact polynomial representation

Replacing \(y\) by \(-y\) preserves \(I_A\) and reverses \(C_A\).
Consequently

\[
 Z_A(u,v)=\mathbb E[\cosh(uI_A)\cosh(vC_A)].                \tag{3}
\]

This block-reversal identity is already established in the variational
note; here it is used to produce a particular fixed-order polynomial.
The possible absolute values satisfy

\[
 |I_A|\in\{0,2,4,6\},\qquad |C_A|\in\{1,3,5,7,9\}.
\]

Let \(h_A(i,j)\) be the number of the \(64\) states with
\((|I_A|,|C_A|)=(i,j)\). Define polynomials by

\[
 T_0(X)=1,\quad T_1(X)=X,\quad
 T_{a+1}(X)=2XT_a(X)-T_{a-1}(X),
\]
\[
 R_0(Y)=1,\quad R_1(Y)=2Y-1,\quad
 R_{b+1}(Y)=2YR_b(Y)-R_{b-1}(Y).
\]

The elementary hyperbolic cosine addition identity gives
\(T_a(\cosh2u)=\cosh(2au)\) and
\(\cosh v\,R_b(\cosh2v)=\cosh((2b+1)v)\). It follows that

\[
 64Z_A(u,v)=\cosh v\,P_A(X,Y),\qquad
 P_A(X,Y)=\sum_{i,j}h_A(i,j)T_{i/2}(X)R_{(j-1)/2}(Y).       \tag{4}
\]

Put

\[
 p=X-Y,\qquad q=Y-1,
 \qquad X=1+p+q,\quad Y=1+q.
\]

The domain \(u\ge v\ge0\) is exactly \(p,q\ge0\). Every coefficient
of \(P_A(1+p+q,1+q)\) is nonnegative. To see this without any general
Chebyshev coefficient assertion, the complete bases needed here are

\[
\begin{array}{c|l}
a&T_a(1+s)\\\hline
0&1\\
1&1+s\\
2&1+4s+2s^2\\
3&1+9s+12s^2+4s^3
\end{array}
\qquad
\begin{array}{c|l}
b&R_b(1+q)\\\hline
0&1\\
1&1+2q\\
2&1+6q+4q^2\\
3&1+12q+20q^2+8q^3\\
4&1+20q+60q^2+56q^3+16q^4.
\end{array}                                                \tag{5}
\]

These expansions follow directly from the displayed recurrences. Their
coefficients are nonnegative, substitution \(s=p+q\) preserves that
property, and all the histogram weights in (4) are nonnegative.

## 3. A candidate attaining the formula

Let \(J\) be the \(3\times3\) all-ones matrix and let \(P\) be any
\(3\times3\) permutation matrix. Choose

\[
 D_L=J-I_3,\qquad D_R=-(J-I_3),\qquad B=J-2P.               \tag{6}
\]

This is a complete signing. Since

\[
 (J-I_3)^2=I_3+J,\qquad
 (J-2P)(J-2P)^T=4I_3-J,\qquad JP=PJ=J,
\]

its full matrix satisfies \(A^2=5I_6\).

Its joint histogram can be obtained directly. Modulo overall sign, the
four equally likely three-spin representatives are

\[
 w_0=(1,1,1),\quad w_1=(1,-1,-1),\quad
 w_2=(-1,1,-1),\quad w_3=(-1,-1,1).
\]

The positive internal triangle has energy three at \(w_0\) and minus
one at the other representatives; the negative triangle has the opposite
energies. Each ordered pair of representatives accounts for four of the
full spin states. Both row and column sums of \(B\) equal one.

For the six pairs with exactly one representative equal to \(w_0\),
\(|I_A|=4\) and \(|C_A|=1\). For the pair \((w_0,w_0)\),
\(I_A=0\) and \(C_A=3\). For \(i,j>0\),

\[
 w_i^TBw_j=1-2w_i^TPw_j
 =\begin{cases}-5,&w_i=Pw_j,\\3,&w_i\ne Pw_j.\end{cases}
\]

The permutation \(P\) permutes \(w_1,w_2,w_3\), so exactly three
of these nine pairs give minus five and the other six give three. Thus
the only nonzero absolute histogram entries are

\[
 h(4,1)=24,\qquad h(0,3)=28,\qquad h(0,5)=12.              \tag{7}
\]

By (4), the candidate polynomial is

\[
\begin{aligned}
 P_*(X,Y)
 &=24T_2(X)+28R_1(Y)+12R_2(Y)\\
 &=16(3X^2+3Y^2+2Y-4),\\
 P_*(1+p+q,1+q)
 &=64+96p+224q+48p^2+96pq+96q^2.                           \tag{8}
\end{aligned}
\]

In particular it is exactly quadratic in \(p,q\). This fact is
essential to the following comparison.

## 4. Universal moments and coefficient comparison

The edge characters \(x_ix_j\) are pairwise orthogonal under uniform
independent spins. Hence, for every complete order-six signing,

\[
 \mathbb EI_A^2=6,\qquad \mathbb EC_A^2=9,
 \qquad \mathbb E(I_A+C_A)^2=15.                           \tag{9}
\]

For any signed simple graph with \(K\) edges, expansion of the fourth
moment of its edge sum gives

\[
 \mathbb E Q^4=3K^2-2K+
 24\sum_{\text{unoriented simple four-cycles }\gamma}
                  \prod_{e\in\gamma}A_e.                 \tag{10}
\]

Indeed, a product of four edge characters has nonzero expectation only
when every vertex has even total degree, with edge repetitions counted.
The terms consisting of one edge repeated four times or two edges each
repeated twice contribute \(K+6\binom K2=3K^2-2K\). Four distinct
edges contribute precisely when they form a simple four-cycle; each such
cycle has \(4!=24\) orders in the expansion. No other case contributes.

The internal graph is the disjoint union of two triangles and has no
four-cycle. Taking \(K=6\) in (10) yields

\[
 \mathbb EI_A^4=96.                                       \tag{11}
\]

For the complete order-six graph, let \(S_4\) be the signed four-cycle
sum in (10). With \(Q=I_A+C_A\),

\[
 \mathbb E Q^4=645+24S_4,
 \qquad \operatorname{tr}A^4=270+8S_4,
 \qquad \mathbb E Q^4=3\operatorname{tr}A^4-165.             \tag{12}
\]

For clarity, the trace counts closed walks of length four. Walks on two
vertices contribute \(6\cdot5=30\); walks on three vertices contribute
\(2\cdot6\cdot5\cdot4=240\). Their edge products are all one.
Each simple four-cycle contributes its signed product eight times, once
per starting vertex and direction. This proves the middle identity.

Write \([p^aq^b]\) for polynomial coefficient extraction after the
substitution \((X,Y)=(1+p+q,1+q)\). For a fixed state with
\(i=|I_A|\), \(j=|C_A|\), (5) gives

\[
 T_{i/2}(1+s)
 =1+\frac{i^2}{4}s+\frac{i^4-4i^2}{96}s^2+
                    \text{terms of degree at least three},
\]
\[
 R_{(j-1)/2}(1+q)
 =1+\frac{j^2-1}{4}q+
       \frac{(j^2-1)(j^2-9)}{96}q^2+
                    \text{terms of degree at least three}.             \tag{13}
\]

Multiplying, substituting \(s=p+q\), and using the factor \(64\) in
(4), equations (9)--(11) give

\[
\begin{aligned}
 [1]P_A&=64,\\
 [p]P_A&=16\mathbb EI_A^2=96,\\
 [q]P_A&=16\mathbb E(I_A^2+C_A^2-1)=224,\\
 [p^2]P_A&=\frac23\mathbb E(I_A^4-4I_A^2)=48,\\
 [pq]P_A&=96+4\mathbb E[I_A^2(C_A^2-1)].                  \tag{14}
\end{aligned}
\]

The last excess is nonnegative pointwise because \(C_A\) is odd,
so \(C_A^2\ge1\).

Block reversal also gives \(\mathbb E I_A^3C_A=
\mathbb E I_AC_A^3=0\). Therefore the remaining quadratic coefficient
is

\[
\begin{aligned}
 [q^2]P_A
 &=\frac23\mathbb E\bigl[I_A^4+6I_A^2C_A^2+C_A^4
                    -10(I_A^2+C_A^2)+9\bigr]\\
 &=\frac23(\mathbb E Q^4-141)\\
 &=96+2\|A^2-5I_6\|_F^2.                                 \tag{15}
\end{aligned}
\]

The last equality uses (12) and

\[
 \|A^2-5I_6\|_F^2
 =\operatorname{tr}A^4-10\operatorname{tr}A^2+150
 =\operatorname{tr}A^4-150,
\]

since \(A\) is symmetric and \(\operatorname{tr}A^2=30\).

Equations (8), (14), and (15) prove that every coefficient of
\(P_A-P_*\) of total degree at most two is nonnegative. Every higher
coefficient is nonnegative by (5), because \(P_*\) has no higher
terms at all. Thus

\[
 P_A(1+p+q,1+q)-P_*(1+p+q,1+q)
 \quad\text{has only nonnegative coefficients}.            \tag{16}
\]

This proves the claimed dominance throughout \(p,q\ge0\). The
explicit candidate (6) attains equality, proving (1).

## 5. Equality and the two boundaries

Suppose first that \(v>0\), so \(q>0\). If a signing attains (1),
(16) and (15) imply \(\|A^2-5I_6\|_F^2=0\). This necessity also
holds when \(u=v>0\), since the strictly positive \(q^2\) term alone
detects every nonconference signing.

Conversely assume \(A^2=5I_6\). Vertex switching replaces \(A\) by
\(DAD\) for a diagonal sign matrix \(D\); it corresponds to a
bijective change of spin variables and preserves \(Z_A(u,v)\).
Global negation also preserves \(Z_A\), by evenness of cosine hyperbolic.

Let \(\tau_L,\tau_R\in\{-1,1\}\) be the products of the three
edge signs in each internal triangle. By switching separately inside
the two blocks, their matrices can be made

\[
 D_L=\tau_L(J-I_3),\qquad D_R=\tau_R(J-I_3).
\]

The off-diagonal block of \(A^2=5I_6\) is the Sylvester equation

\[
 D_LB+BD_R=0.                                              \tag{17}
\]

If \(\tau_L=\tau_R=\tau\), each internal matrix has spectrum
\(\{2\tau,-\tau,-\tau\}\). No eigenvalue of the first is the
negative of an eigenvalue of the second: all possible sums are
\(4\tau,\tau,-2\tau\), which are nonzero. Diagonalizing the two
real symmetric matrices separately in (17) therefore forces every
entry of the transformed \(B\) to be zero. This contradicts its sign
entries. Thus the triangle products are opposite.

After optional global negation, take \(D_L=J-I_3\) and
\(D_R=-(J-I_3)\). Equation (17) becomes \(JB=BJ\), so every row
and every column sum of \(B\) has a common value \(r\). The first
diagonal block of \(A^2=5I_6\) gives

\[
 BB^T=4I_3-J.
\]

Applying this to the all-ones vector gives \(r^2=1\). If necessary,
switching all three vertices of the second block replaces \(B\) by
\(-B\), without changing its internal matrix, and makes \(r=1\).
Each row and column then contains exactly one negative entry. Hence
\(B=J-2P\) for a permutation matrix \(P\), exactly the candidate
form (6). This proves sufficiency and the classification for \(v>0\).

The count of twelve switching-normalized minimizers also has a short
analytic proof. Normalize the first row off the diagonal to be positive.
For each of the other five vertices, \((A^2)_{0i}=0\) says that among
its four incident edges within those five vertices, two are positive
and two are negative. The positive edges consequently form a simple
two-regular graph on five vertices, necessarily a five-cycle. There
are \((5-1)!/2=12\) labelled undirected five-cycles. They all occur:
the explicit conference candidate yields one after normalization, and
permuting the remaining five vertices yields every labelled five-cycle.
Thus exactly twelve normalized signings satisfy \(A^2=5I_6\).

At \(v=0\), the cross block is absent. For opposite triangle products,
the absolute internal histogram is \(h_I(0)=40,h_I(4)=24\), giving

\[
 P_{\rm opp}(X,1)=16+48X^2.
\]

For equal triangle products, the corresponding counts are
\(h_I(2)=60,h_I(6)=4\), giving

\[
 P_{\rm same}(X,1)=16X^3+48X,
 \qquad P_{\rm same}(X,1)-P_{\rm opp}(X,1)=16(X-1)^3.       \tag{18}
\]

For \(u>0\), this difference is strictly positive. This proves the
boundary classification. At \(u=v=0\), every partition function is one.

## 6. Exact endpoint comparison and its unique crossing

Equation (1) gives the two optimized endpoint partition functions

\[
 Z_0(c):=e^{f_6(c,0)}
       =\frac{3\cosh(4c/\sqrt3)+5}{8},
\]
\[
 Z_1(c):=e^{f_6(c,1)}
       =\frac{3\cosh(5c/\sqrt6)+5\cosh(3c/\sqrt6)}{8}.      \tag{19}
\]

There is exactly one \(c_*>0\) at which these are equal. More precisely,
\(f_6(c,1)>f_6(c,0)\) for \(0<c<c_*\), while the inequality is
reversed for \(c>c_*\).

To prove this, put \(b=c/\sqrt6\) and write

\[
 G(b):=8(Z_1(c)-Z_0(c))
 =3\cosh5b+5\cosh3b-3\cosh(4\sqrt2\,b)-5.                \tag{20}
\]

Its constant term is zero and its \(b^2\) coefficient is twelve.
For every integer \(k\ge2\), its \(b^{2k}\) coefficient is strictly
negative, since

\[
 3\,25^k+5\,9^k<3\,32^k.
\]

At \(k=2\), the two sides are \(2280<3072\); after division by
\(32^k\), both summands on the left decrease with \(k\). Thus
\(G(b)/b^2\) is strictly decreasing for \(b>0\): its power series
has constant term twelve and every higher nonzero coefficient is
negative, with termwise differentiation valid on compact intervals.
It starts positive and tends to minus infinity, because
\(4\sqrt2>5\). Continuity proves the unique positive zero and all the
stated signs. Taking logarithms preserves those signs.

This is a fixed-order endpoint transition. It is not a classification
of all stationary points of the interior profile.

## 7. An unbounded-in-temperature interior excursion

The same exact minimum can be written with positive weights as

\[
 Z_*(u,v)=\frac1{16}\bigl[
 3\cosh(4u+v)+3\cosh(4u-v)+3\cosh5v+7\cosh3v\bigr].       \tag{21}
\]

This follows by expanding the squared cosines hyperbolic in (1) and
then applying the product-to-sum identity. When \(u\ge v\ge0\), the
largest argument in (21) is \(4u+v\). Along (2), it is

\[
 M(c,t)=\frac c{\sqrt6}\bigl(4\sqrt{2-t}+\sqrt t\bigr).
\]

Its maximum on \([0,1]\) occurs at \(t=2/17\), as follows by
setting its derivative to zero in the interior; the derivative changes
from positive to negative there. The value is

\[
 M(c,2/17)=c\sqrt{17/3},\qquad M(c,0)=4c/\sqrt3.           \tag{22}
\]

Only this maximal-energy function is being differentiated in (22), not
the finite-temperature log pressure. At the fixed point \(t=2/17\),
we have \(u=4c/\sqrt{51}\), \(v=c/\sqrt{51}\). The four arguments
in (21) are then \(17v,15v,5v,3v\), so as \(c\to\infty\),

\[
 f_6(c,2/17)=c\sqrt{17/3}+\log(3/32)+o(1).
\]

At \(t=0\), (19) similarly gives
\(f_6(c,0)=4c/\sqrt3+\log(3/16)+o(1)\). Hence

\[
 \boxed{f_6(c,2/17)-f_6(c,0)
 =\frac{\sqrt{17}-4}{\sqrt3}\,c-\log2+o(1).}              \tag{23}
\]

Since \(\sqrt{17}>4\), no bound uniform in \(c\) can keep the
whole exact optimized order-six profile within a bounded additive amount
above its \(t=0\) endpoint.

An endpoint derivative must not be inferred by differentiating the
largest argument alone. For example, with \(b=c/\sqrt6>0\), direct
differentiation of (21) along (2) gives

\[
 \partial_t f_6(c,1-)
 =\frac{3b(\sinh5b+\sinh3b)}
        {6\cosh5b+10\cosh3b}>0.                           \tag{24}
\]

At this endpoint, the two distinct terms with arguments \(4u+v\)
and \(5v\) tie at the largest value and have different derivatives;
both contribute. Thus (23) neither asserts monotonicity of the profile
nor identifies its finite-temperature maximum with \(t=2/17\).

## 8. Exact scope and the obstruction to immediate extension

The proof of (16) succeeds because the candidate polynomial is exactly
quadratic in \(p,q\). Matching the low moments and comparing the two
remaining quadratic coefficients exhaust every coefficient that the
candidate asks a competing signing to dominate. All higher coefficients
of the competing signing are nonnegative, and the candidate has zero
coefficients there.

For larger orders, a candidate generally has nonzero higher-degree
coefficients. Nonnegativity of each competitor's coefficients does not
compare them with those positive candidate coefficients. No such
higher-coefficient comparison has been proved here. In particular, this
argument does not say that fourth moments determine critical pressure
in general, that conference signings minimize at other orders, or that
extending the numerical grid or signing census supplies a missing
all-orders implication.

Equation (23) excludes a temperature-uniform bounded-excursion claim
already for an actual global minimum at this one order. It does not
refute an estimate with an arbitrary dependence on fixed \(c\), an
\(o(N)\) error as \(N\to\infty\), an integrated order comparison,
or convergence by another method. The original question about
\(m_n/n^{3/2}\), and the value of a possible limit, remain open.
