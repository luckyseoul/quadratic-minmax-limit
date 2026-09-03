# Finite-anchor signature-cell tournament shielding

**Status:** proved all-orders construction and an exact conditional cover
criterion.  Any fixed family of Boolean anchor states can be simultaneously
shielded from the skew bilinear form at `O(n)` cost.  This rules out a finite
near-maximizer list as the true obstruction, but it does not control pairs
outside that list and does not settle the original MO problem.

## 1. Signature-cell construction

Let

\[
 \mathcal X=\{x^{(1)},\ldots,x^{(k)}\}\subseteq\{\pm1\}^n
\]

be deduplicated modulo global sign.  Put
`z^(a)=x^(1) circ x^(a)` and partition `[n]` into the `L` nonempty cells on
which

\[
 (z_i^{(2)},\ldots,z_i^{(k)})
\]

is constant.  Thus

\[
 L\le\min(n,2^{k-1}).                                      \tag{1}
\]

There is a skew signing `R` such that, simultaneously for every anchor,

\[
 \boxed{\|Rx^{(a)}\|_1\le Ln,\qquad
 |(x^{(a)})^TRy|\le Ln\quad\hbox{for every Boolean }y.}     \tag{2}
\]

### Proof

Gauge the first anchor to the all-one vector: let
`D=diag(x^(1))`.  It is enough to construct a tournament matrix `T` in the
gauged coordinates and put `R=DTD`, because

\[
 \|Rx^{(a)}\|_1=\|Tz^{(a)}\|_1.
\]

Order the `L` signature cells and locally index the vertices in each cell.
Inside an odd cell use a regular tournament, whose row sums are all zero.
Inside an even cell, take a regular tournament on one additional vertex and
delete that vertex; every remaining row sum is `+1` or `-1`.

Between two cells `V_c,V_d`, with local indices `r,s`, set

\[
 T_{i_rj_s}=(-1)^{r+s},                                    \tag{3}
\]

and define the reverse block by skew symmetry.  Every row and column sum of
this alternating rectangular block is zero or `+/-1`.

Each `z^(a)` is constant, say `chi_c^(a)`, on `V_c`.  Therefore, for
`i in V_b`,

\[
 (Tz^{(a)})_i
 =\sum_{c=1}^L\chi_c^{(a)}
   \sum_{\substack{j\in V_c\\j\ne i}}T_{ij},
\]

and every inner sum has magnitude at most one.  Thus every coordinate has
magnitude at most `L`, and summing proves the first part of (2).  Finally,

\[
 |(x^{(a)})^TRy|
 =|(R^Tx^{(a)})\cdot y|\le\|Rx^{(a)}\|_1,
\]

which proves the second part.

## 2. Exact parity boundary

For one anchor the optimum is

\[
 \min_R\|Rx\|_1=\begin{cases}0,&n\text{ odd},\\n,&n\text{ even}.
 \end{cases}                                               \tag{4}
\]

After gauging `x` to one, every tournament row sum has parity `n-1`;
regular and nearly regular tournaments attain the stated floors.

Two projectively distinct Boolean anchors can never both lie exactly in the
kernel.  Gauge the first to one and split the second signature into
nonempty cells `P,Q` of sizes `p,q`.  For a vertex `i`, let `a_i` be its
row sum into its own cell and `b_i` its row sum into the other.  The two
kernel equations force `a_i=b_i=0`.  On `P` this requires `p` odd and `q`
even; on `Q` it requires `q` odd and `p` even, a contradiction.

More precisely,

\[
 \|T\mathbf1\|_1+\|T\tau\|_1
 =2\sum_i\max(|a_i|,|b_i|).                                \tag{5}
\]

Parity makes (5) at least `2n` for even `n`.  For odd `n` it gives only
twice the size of the even signature cell, which can be `O(1)`; there is no
uniform odd-order `Omega(n)` two-anchor floor.

## 3. Conditional high-difference cover

Let

\[
 \Delta_A(x,y)=|Q_A(x)-Q_A(y)|,qquad
 \Gamma(R)=\max_{x,y}|x^TRy|.
\]

The outgoing-half identity is

\[
 4D_{\to}(A,S)=\max_{x,y}
 \bigl(\Delta_A(x,y)+|x^TRy|\bigr),qquad S=A\circ R.       \tag{6}
\]

For every pair incident to an anchor in `mathcal X`, (2) and
`Delta_A<=2M` give

\[
 \Delta_A(x,y)+|x^TRy|\le2M+Ln.                            \tag{7}
\]

Define the high-difference graph

\[
 E_\Gamma=\left\{\{x,y\}:
 \Delta_A(x,y)>2\sqrt2M-\Gamma(R)\right\}.                 \tag{8}
\]

If `mathcal X` is a vertex cover of `E_Gamma` and

\[
 Ln\le2(\sqrt2-1)M,                                        \tag{9}
\]

then the target follows exactly.  Covered pairs satisfy (7) and (9), while
uncovered pairs satisfy

\[
 \Delta_A(x,y)+|x^TRy|
 \le2\sqrt2M-\Gamma(R)+\Gamma(R)=2\sqrt2M.
\]

Thus

\[
 \boxed{D_{\to}(A,S)\le M/\sqrt2.}                         \tag{10}
\]

Proposition 6.5e applied to the complete signing `A` gives the refined
universal finite-order bound

\[
 M\ge {n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}}.
\]

Thus (9) is automatic for every `A` whenever

\[
 L\le {2(\sqrt2-1)(n-1)\over\pi}
       \arcsin{1\over\sqrt{n-1}}
 =\left({2(\sqrt2-1)\over\pi}+o(1)\right)\sqrt n
 =(0.2636965\ldots+o(1))\sqrt n.                           \tag{10a}
\]

In the worst case `L<=2^(k-1)`, so arbitrary anchor families through about
`(1/2)log_2 n-O(1)` states still fit inside the critical margin on every
anchor-incident pair.

This is a conditional theorem, not a closure: the checkerboard tournament
can have `Gamma(R)=Theta(n^2)`, making `E_Gamma` far too dense.  The
construction controls anchor-incident pairs but supplies no global bound on
the complement.

For asymptotic use, the normalized anchor cost is `L/sqrt(n)`.  A
Dini-admissible error requires

\[
 \sum_{j\ge0}{L(2^jn_0)\over\sqrt{2^jn_0}}<\infty.         \tag{11}
\]

Fixed `k` satisfies this automatically, as does
`L(n)=O(n^(1/2-delta))`; the mere condition `L=o(sqrt(n))` is not enough.

## 4. Consequence for the live obstruction

A finite list of extremizers, or more generally any family with few
coordinate-signature cells, can be removed at lower-order cost.  Therefore a
failure of the directed target cannot be explained by one maximizer orbit or
any fixed collection of anchors.  It must come from a diffuse energy layer
with growing signature complexity, or from uncontrolled non-anchor pairs.
The remaining task is to control `Gamma(R)` jointly with the high-difference
graph, not to enumerate more finite anchors.
