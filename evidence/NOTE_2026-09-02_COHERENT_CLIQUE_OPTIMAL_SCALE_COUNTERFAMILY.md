# An optimal-order counterfamily to coherent clique-flip control

**Status:** proved theorem / asymptotic counterexample.  There is an infinite
family of complete signings `A` with

\[
 \Phi(A)=\Theta(n^{3/2})
\]

for which

\[
 L_{\rm cl}(A)>\sqrt2\,\Phi(A)+c n^{3/2}
\]

with an absolute constant `c>0`.  Thus the coherent fourth-phase estimate
cannot follow from the order-optimal hypothesis
`Phi(A)=O(n^(3/2))`.  This does **not** settle the same estimate for global
minimizers `Phi(A)=m_n`: the matrices below have a deliberately large fixed
leading constant.

This note advances the order-four and quadratic-scale obstructions in
`NOTE_2026-09-02_COMPLEXIFICATION_OPPOSITE_DIAGONAL_AUDIT.md`.  It uses the
Hadamard saddle gadget proved there, but dilutes it into many growing blocks
and fills every missing edge.  The resulting matrices are complete and have
the correct `n^(3/2)` order.  No finite-order census is involved.

## 1. Statement

For a complete signing `A`, recall

\[
 \Phi(A)=\max_{x\in\{\pm1\}^n}|Q_A(x)|,
 \qquad
 L_{\rm cl}(A)=\max_{T\subseteq[n]}\Phi(A^{K_T}),
\]

where `A^(K_T)` reverses every edge with both endpoints in `T`.

Let `b>=2`, and suppose that a Hadamard matrix of order `r>=4` exists.  Put

\[
 n=4rb,
 \qquad
 N_\times={b\choose2}(4r)^2=8r^2b(b-1).
\]

Then there is a complete signing `A=A_(b,r)` of order `n` such that

\[
 \boxed{
 \Phi(A)\le
 b(2r^2+2)
 +\sqrt{16r^2b(b-1)\bigl(4rb\log2+1\bigr)}
 }                                                        \tag{1}
\]

and

\[
 \boxed{
 L_{\rm cl}(A)\ge b(4r^2-2r).
 }                                                        \tag{2}
\]

In particular, fix `kappa=256`, take `b=2^j`, and take
`r=kappa b`.  Sylvester Hadamard matrices exist at every such order.  Since
`n=4 kappa b^2`, (1)--(2) imply

\[
 \liminf_{j\to\infty}{L_{\rm cl}(A)\over n^{3/2}}
 \ge {\sqrt\kappa\over2},                                 \tag{3}
\]

whereas

\[
 \limsup_{j\to\infty}{\Phi(A)\over n^{3/2}}
 \le {\sqrt\kappa\over4}+\sqrt{\log2}.                    \tag{4}
\]

Consequently

\[
 \liminf_{j\to\infty}
 {L_{\rm cl}(A)-\sqrt2\,\Phi(A)\over n^{3/2}}
 \ge
 {\sqrt\kappa(2-\sqrt2)\over4}-\sqrt{2\log2}
 =1.1657\ldots>0.                                         \tag{5}
\]

The ratio certified by these bounds tends to at least

\[
 {2\sqrt\kappa\over\sqrt\kappa+4\sqrt{\log2}}
 =1.6554\ldots>\sqrt2.                                    \tag{6}
\]

Proposition 5.2's universal bound
`Phi(A)>=n sqrt(n-1)/pi`, together with (1), gives
`Phi(A)=Theta(n^(3/2))`.  Thus this is an order-optimal family, although not
a family known (or claimed) to have the optimal leading constant.

## 2. The internal saddle blocks

We use the exact gadget already proved in the complexification audit.  If
`C` is an order-`r` Hadamard matrix, set `k=2r` and

\[
 P=J_k-I_k,
 \qquad
 H=C\otimes\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
 \qquad
 G_r=\begin{pmatrix}P&H\\H^T&-P\end{pmatrix}.              \tag{7}
\]

This is a complete signing of order `4r`, and the exact calculation in the
earlier audit gives

\[
 \Phi(G_r)=2r^2+2.                                         \tag{8}
\]

Partition its vertices into the first and second `2r` coordinates.  Put the
second half into a clique `T_r` and evaluate the clique-flipped signing at
the all-one state.  The first diagonal block contributes
`2r^2-r`; after the flip, the second diagonal block contributes the same
amount.  The cross block contributes zero because every row and column sum
of the displayed order-two kernel is zero.  Hence

\[
 Q_{G_r^{K_{T_r}}}({\bf1})=4r^2-2r.                        \tag{9}
\]

The new point is to retain the sum of (9) over a growing number of blocks
while completing all inter-block edges at only `O(n^(3/2))` cost.

## 3. A complete low-norm inter-block filler

Take `b` disjoint copies of `G_r`.  On every edge joining two different
copies, independently choose a Rademacher sign, and call the resulting
zero-within-block weighted signing `W`.  There are exactly `N_\times` such
random edges.

For each fixed Boolean state `x`, `Q_W(x)` is a sum of `N_\times` independent
Rademacher variables.  Hoeffding's inequality gives

\[
 \Pr\{|Q_W(x)|\ge t\}
 \le2\exp\left(-{t^2\over2N_\times}\right).               \tag{10}
\]

The energies for `x` and `-x` agree, so there are only `2^(n-1)` distinct
states.  With

\[
 t=\sqrt{2N_\times(n\log2+1)},                             \tag{11}
\]

the union bound in (10) is at most `e^(-1)<1`.  Therefore a deterministic
choice of the inter-block signs exists with

\[
 \Phi(W)<t.                                                 \tag{12}
\]

Fill the `b` diagonal blocks with `G_r` and all inter-block positions with
this `W`; call the resulting complete signing `A`.  The triangle inequality
and (8) give

\[
 \Phi(A)\le b\Phi(G_r)+\Phi(W)
 <b(2r^2+2)+t,
\]

which is (1).

## 4. One global clique retains every saddle gain

Let `T` be the union of the second half `T_r` in every copy of `G_r`, and
evaluate `A^(K_T)` at the all-one state.  Its within-block contribution is,
by (9),

\[
 b(4r^2-2r).                                               \tag{13}
\]

The inter-block contribution is one real number `Z(W)`.  Replacing every
inter-block sign by its negative preserves (12) and negates `Z(W)`.  Choose
between `W` and `-W` so that `Z(W)>=0`.  Then (13) proves (2).  Notice that
this sign choice is global; it neither changes an internal gadget nor
depends on a maximizing Boolean state.

This is also why completing the graph does not wash out the clique witness.
The filler has leading-order switching norm, but its value on the one
required witness can always be given the favorable sign.

## 5. Exact asymptotics

Set `r=kappa b`.  Dividing (2) by

\[
 n^{3/2}=(4\kappa b^2)^{3/2}=8\kappa^{3/2}b^3
\]

gives

\[
 {b(4r^2-2r)\over n^{3/2}}
 ={\sqrt\kappa\over2}-{1\over4\sqrt\kappa\,b}.           \tag{14}
\]

The internal part of (1), after the same normalization, is

\[
 {b(2r^2+2)\over n^{3/2}}
 ={\sqrt\kappa\over4}
  +{1\over4\kappa^{3/2}b^2},                              \tag{15}
\]

and its filler part is exactly

\[
 {t\over n^{3/2}}
 =\sqrt{\left(1-{1\over b}\right)
         \left(\log2+{1\over4\kappa b^2}\right)}.         \tag{16}
\]

Equations (14)--(16) prove (3)--(5).  For `kappa=256`, the right side of
(5) is positive, so no error `o(n^(3/2))` -- and hence no dyadic
Dini-summable error of that order -- can repair the coherent estimate on
the class `Phi(A)=O(n^(3/2))`.

## 6. Scope of the no-go

The theorem rules out the implication

\[
 \Phi(A)=O(n^{3/2})
 \quad\Longrightarrow\quad
 L_{\rm cl}(A)\le\sqrt2\Phi(A)+o(n^{3/2}).                 \tag{17}
\]

It also rules out any proof that uses only a scale bound on `Phi(A)`, even
if it couples arbitrarily many phase faces: the counterfamily is already a
complete signing and the gap is leading order.

It does **not** rule out

\[
 \Phi(A)=m_n
 \quad\Longrightarrow\quad
 L_{\rm cl}(A)\le\sqrt2m_n+o_{\rm Dini}(n^{3/2}).           \tag{18}
\]

Indeed the upper bound in (4) has leading constant
`sqrt(kappa)/4+sqrt(log 2)`, much larger than the best known global scale.
Any surviving coherent argument must therefore use genuine global
minimality, or at least a quantitatively near-minimal leading constant; the
bare `Theta(n^(3/2))` hypothesis is now rigorously insufficient.
