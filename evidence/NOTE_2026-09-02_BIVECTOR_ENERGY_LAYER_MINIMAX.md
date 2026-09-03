# Bivector energy-layer minimax and low-degree no-go

**Status:** exact all-orders reduction and quantified convex/moment
obstructions.
This note does not prove the original MO limit.  It identifies the remaining
orientation problem as an `A`-weighted covering problem on decomposable
Boolean bivectors and proves that the ordinary fractional relaxation, the
displayed covariance relaxation, and normalized single-row even-moment
certificates are subcritical in the stated range.  It does not exclude every
fixed-degree SOS or Pluecker-aware hierarchy.  No small-order census is used.

Let `A` be an order-`n` signing, put `M=Phi(A)`, and identify the
upper-triangular entries of a skew signing `R` with
`r in {+1,-1}^N`, where `N=binom(n,2)`.  For Boolean states `x,y`, define

\[
 b_{xy,ij}={x_i y_j-x_j y_i\over2},\qquad
 d_{xy}=M-{|Q_A(x)+Q_A(y)|\over2}.                          \tag{1}
\]

Then `b_xy,ij` belongs to `{0,+1,-1}` and

\[
 \langle b_{xy},r\rangle={1\over2}x^TRy.                   \tag{2}
\]

## 1. Exact recentered minimax

The skew completion functional from Proposition 6.5 is

\[
 K(A,R)={1\over2}\max_{x,y}
 \bigl(|Q_A(x)+Q_A(y)|+|x^TRy|\bigr).
\]

Equations (1)--(2) give exactly

\[
 \boxed{
 K(A,R)-M
 =\max_{x,y}\bigl(|\langle b_{xy},r\rangle|-d_{xy}\bigr).
 }                                                           \tag{3}
\]

Consequently

\[
 \min_RK(A,R)=M+min_{r\in\{\pm1\}^N}
 \max_{x,y}\bigl(|\langle b_{xy},r\rangle|-d_{xy}\bigr).  \tag{4}
\]

Put

\[
 \epsilon_A(x,y)=M-\max\{|Q_A(x)|,|Q_A(y)|\}.
\]

The endpoint identity `|a+b|+|a-b|=2max(|a|,|b|)` gives

\[
 d_{xy}={|Q_A(x)-Q_A(y)|\over2}+\epsilon_A(x,y),            \tag{5}
\]

and hence

\[
 \boxed{
 2(K(A,R)-M)=\max_{x,y}
 \left(|x^TRy|-|Q_A(x)-Q_A(y)|-2\epsilon_A(x,y)\right).
 }                                                           \tag{6}
\]

An excess `t` in (3) is achieved if and only if

\[
 |x^TRy|\le |Q_A(x)-Q_A(y)|+2\epsilon_A(x,y)+2t             \tag{7}
\]

for every pair.  At the multiplier-two threshold,
`t=(sqrt(2)-1)M` up to the Dini error.

## 2. Energy layers and the two directed halves

Write `e(x)=M-|Q_A(x)|`.  The free term in (7), excluding `2t`, is

\[
 |Q_A(x)-Q_A(y)|+2\min(e(x),e(y))
 =\begin{cases}
 e(x)+e(y),&Q_A(x)Q_A(y)\ge0,\\
 2M-|e(x)-e(y)|,&Q_A(x)Q_A(y)\le0.
 \end{cases}                                                \tag{8}
\]

The formulas agree when one energy is zero.  Thus the small-width rows are
precisely pairs in the same positive or negative near-maximizer layer.

There is also an exact directed-cut form.  Let
`U={i:x_i=-y_i}`, put `S=A circ R`, and split the
`A_ij y_i y_j` energy across the cut into `F` on arcs from `U` to `U^c`
and `G` on the reverse arcs.  Then

\[
 Q_A(y)-Q_A(x)=2(F+G),\qquad x^TRy=2(G-F).                   \tag{9}
\]

Therefore

\[
 |x^TRy|-|Q_A(x)-Q_A(y)|
 =\begin{cases}
 4\min(|F|,|G|),&FG<0,\\
 -4\min(|F|,|G|),&FG\ge0.
 \end{cases}                                                \tag{10}
\]

A pair violates (7) if and only if

\[
 \boxed{FG<0\quad\hbox{and}\quad
 2\min(|F|,|G|)>t+\epsilon_A(x,y).}                         \tag{11}
\]

Same-sign directed halves are automatically harmless.  This removes them
from the live problem exactly, rather than heuristically.

## 3. The linear relaxation is identically zero

For the cube relaxation define

\[
 \tau_{\rm frac}=\min_{r\in[-1,1]^N}\max_{x,y}
 \bigl(|\langle b_{xy},r\rangle|-d_{xy}\bigr).             \tag{12}
\]

Taking `r=0` proves `tau_frac<=0`.  Conversely, if `z` is an
`A`-maximizer, the row `(z,z)` has `b_zz=0=d_zz`, so every `r` has
objective at least zero.  Hence

\[
 \boxed{\tau_{\rm frac}=0.}                                \tag{13}
\]

The exact affine dual is

\[
 \tau_{\rm frac}
 =\max_{\lambda\in\Delta(\{(x,y,\sigma)\})}
 \left[-\sum\lambda_{xy\sigma}d_{xy}
 -\left\|\sum\lambda_{xy\sigma}\sigma b_{xy}\right\|_1\right]
 =0.                                                        \tag{14}
\]

Thus no distribution of signed affine constraints detects the orientation
integrality gap.  The dual optimum can sit entirely on a diagonal endpoint
row.

## 4. Exact nonlinear covering dual

For a fixed threshold `t`, define the zero--one bad-pair matrix

\[
 H_{(x,y),r}
 =\mathbf1\{|\langle b_{xy},r\rangle|>d_{xy}+t\}.           \tag{15}
\]

Finite minimax duality gives

\[
 \boxed{
 \eta_t
 =\max_{\lambda\in\Delta(\{(x,y)\})}
     \min_{r\in\{\pm1\}^N}\mathbb E_\lambda H_{(x,y),r}
 =\min_{p\in\Delta(\{\pm1\}^N)}
     \max_{x,y}\mathbb E_pH_{(x,y),r}.
 }                                                           \tag{16}
\]

This has exact certificate meaning:

- `eta_t=0` if and only if a successful orientation exists;
- `eta_t>0` if and only if the target fails;
- a failure certificate is a state-pair distribution whose bad sets
  fractionally cover every tournament with positive mass;
- if the right side is zero, every orientation in the support of `p` is an
  actual successful certificate.

The witness distribution can be symmetrized under exchanging `x,y` and
independently negating either state.

For a uniformly random orientation, Hoeffding and Proposition 5.2 give at
the critical `t=(sqrt(2)-1)M`

\[
 \eta_t\le
 2\exp\left(-{2(\sqrt2-1)^2\over\pi^2}(n-1)\right).         \tag{17}
\]

If failure occurs, the uniform distribution on all `4^n` ordered state
pairs gives the elementary lower bound `eta_t>=4^(-n)`.  A counterexample
cover is therefore allowed to be exponentially thin; proving only
`eta_t=o(1)` cannot establish success.

## 5. Decomposable bivector geometry

The constraint normals are not arbitrary:

\[
 b_{xy}={x\wedge y\over2}\in\Lambda^2\mathbb R^n.           \tag{18}
\]

They satisfy

\[
 \|b_{xy}\|_2^2={n^2-(x\cdot y)^2\over4}
 =d_H(x,y)(n-d_H(x,y)),                                    \tag{19}
\]

\[
 \langle b_{xy},b_{uv}\rangle
 ={(x\cdot u)(y\cdot v)-(x\cdot v)(y\cdot u)\over4},     \tag{20}
\]

and every Pluecker relation

\[
 b_{ij}b_{kl}-b_{ik}b_{jl}+b_{il}b_{jk}=0.                 \tag{21}
\]

Uniform independent `x,y` give the tight-frame identity

\[
 \mathbb E_{x,y}b_{xy}b_{xy}^T={1\over2}I_N.               \tag{22}
\]

This is the genuine extra geometry erased by an arbitrary-row discrepancy
relaxation.

## 6. Covariance and normalized single-row moments are subcritical

Put `w_xy=t+d_xy` and consider the elliptope relaxation

\[
 \rho_t=
 \min_{\substack{X\succeq0\\\operatorname{diag}X=1}}
 \max_{x,y}{b_{xy}^TXb_{xy}\over w_{xy}^2}.                 \tag{23}
\]

Its exact SDP dual is

\[
 \rho_t=\max_{\lambda\in\Delta}
 \max_{\substack{z\in\mathbb R^N\\
 H_\lambda-\operatorname{Diag}(z)\succeq0}}
 \sum_ez_e,\qquad
 H_\lambda=\sum_{x,y}\lambda_{xy}{b_{xy}b_{xy}^T\over w_{xy}^2}.
                                                                    \tag{24}
\]

At the critical threshold, Proposition 5.2 gives
`t>=(sqrt(2)-1)n sqrt(n-1)/pi`, while (19) is at most `n^2/4`.
The feasible point `X=I_N` therefore proves

\[
 \boxed{
 \rho_t\le {\pi^2\over4(\sqrt2-1)^2(n-1)}
 ={14.3810675\ldots\over n-1}<1\quad(n\ge16).
 }                                                           \tag{25}
\]

Thus no covariance/elliptope witness can certify failure at the critical
threshold from order 16 onward, even after recognizing that all rows obey
the Pluecker identities.

More generally, Khintchine's inequality for uniform random `r` gives, for
every integer `p>=1`,

\[
 \mathbb E_r\left(
 { |\langle b_{xy},r\rangle|\over w_{xy}}
 \right)^{2p}
 \le(2p-1)!!
 \left({\pi^2\over4(\sqrt2-1)^2(n-1)}\right)^p.             \tag{26}
\]

Averaging against any proposed row distribution shows that some orientation
obeys the same averaged bound.  Since `(2p-1)!!<=(2p)^p`, no certificate of
the natural form

\[
 \inf_r\mathbb E_\lambda
 \left({|\langle b,r\rangle|\over w}\right)^{2p}>1
\]

can exist for

\[
 p<{2(\sqrt2-1)^2\over\pi^2}(n-1)
 =0.0347679\ldots(n-1).                                    \tag{27}
\]

In particular every fixed-degree normalized single-row even-moment
certificate of this form is blind for all sufficiently large orders.  This
calculation does not rule out fixed-degree SOS certificates with cross-row
products or reductions by the Pluecker ideal.  A hierarchy continuing this
specific moment majorant would have to reach degree `Omega(n)` and couple
the `A`-dependent widths to bivector intersections.

## 7. Quantitative concentration forced on the energy layers

Every skew signing satisfies

\[
 G(R):=\max_{x,y}|x^TRy|\ge n\mu_{n-1}.                    \tag{28}
\]

For a successful orientation, every pair attaining `G(R)` must therefore
satisfy

\[
 |Q_A(x)-Q_A(y)|+2\epsilon_A(x,y)\ge G(R)-2t.              \tag{29}
\]

If such a pair has same-sign energies, (8) becomes

\[
 e(x)+e(y)\ge G(R)-2t.                                     \tag{30}
\]

Using `M<=(1/2+o(1))n^(3/2)` at the critical threshold makes the right
side at least

\[
 \left(\sqrt{2/\pi}-(\sqrt2-1)-o(1)\right)n^{3/2}
 =(0.383670998\ldots-o(1))n^{3/2}.                          \tag{31}
\]

A single chosen state is never the obstruction.  Given any Boolean `x_0`,
take a regular tournament matrix `T` for odd `n`, or a nearly regular one
with row sums `+/-1` for even `n`, and set
`R=D_(x_0) T D_(x_0)`.  Then

\[
 \max_y|x_0^TRy|=0\quad(n\text{ odd}),\qquad
 \max_y|x_0^TRy|\le n\quad(n\text{ even}).                 \tag{32}
\]

Thus any genuine obstruction must spread across many mutually incompatible
same-sign high-energy states.

## 8. Verdict

The surviving implication is exact: construct an orientation whose large
skew bilinear pairs avoid the simultaneous same-sign near-maximizer layers
of `Q_A`, equivalently rule out every positive fractional cover in (16).
The linear dual is identically zero, the displayed covariance relaxation is
subcritical by (25), and the normalized single-row moment certificate is
subcritical by (27).  Re-running those exact relaxations cannot advance the
proof.  A new theorem may instead use a richer SOS/Pluecker coupling, couple
the `A`-dependent defects `d_xy` to intersections of decomposable Boolean
bivectors at degree proportional to `n`, or directly round the nonlinear
cover.  The original MO limit remains open.
