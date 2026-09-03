# Signed-regular arcsine rigidity for the directed half-cut

**Status:** proved all-orders lower bound with a positive square correction.
The leading constant meets the outgoing-half target exactly at the universal
`1/pi` scale.  This proves analytic sharpness and gives a new rigidity
condition, but it does not construct the orientation needed for the upper
bound and therefore does not settle the original MO problem.

## 1. Signed regular graph lemma

Let `H` be a symmetric zero-diagonal matrix of order `N`, with entries in
`{0,+1,-1}`, whose nonzero support is `d`-regular, `d>=2`.  Define

\[
 \Phi(H)=\max_{z\in\{\pm1\}^N}
 \left|\sum_{i<j}H_{ij}z_iz_j\right|.
\]

Then

\[
 \boxed{
 \Phi(H)\ge {Nd\over\pi}\arcsin{1\over\sqrt d}
 +{1\over4\pi d^{5/2}(1-1/d)^{3/2}}
   \sum_{\substack{i<j\\H_{ij}\ne0}}(H^2)_{ij}^2.
 }                                                           \tag{1}
\]

This strengthens the base Gaussian arcsine lower bound whenever adjacent
rows have nonzero signed codegree.

### Proof

Let `g` be a standard Gaussian vector and put

\[
 z^\pm=(I\pm H/\sqrt d)g,\qquad x_i^\pm=\operatorname{sgn}(z_i^\pm).
\]

Every coordinate of `z^+` and `z^-` has variance two.  For a support edge
`e={i,j}`,

\[
 r_e^\pm=\operatorname{Corr}(z_i^\pm,z_j^\pm)
 ={(H^2)_{ij}\over2d}\pm{H_{ij}\over\sqrt d}.              \tag{2}
\]

Put

\[
 u_e={H_{ij}(H^2)_{ij}\over2d},\qquad v={1\over\sqrt d}.
\]

The Gaussian sign identity
`E[sgn Z sgn W]=(2/pi)arcsin(Corr(Z,W))` and oddness of arcsine give

\[
 H_{ij}\bigl(\arcsin r_e^+-\arcsin r_e^-\bigr)
 =\Delta(u_e,v),
\quad
 \Delta(u,v)=\arcsin(u+v)-\arcsin(u-v).                    \tag{3}
\]

Therefore

\[
 \mathbb E Q_H(x^+)-\mathbb E Q_H(x^-)
 ={2\over\pi}\sum_e\Delta(u_e,v)\le2\Phi(H),              \tag{4}
\]

where the sum is over the `Nd/2` support edges.

It remains to quantify the midpoint gain.  With
`f(t)=(1-t^2)^(-1/2)`,

\[
 \Delta(u,v)=\int_{u-v}^{u+v}f(t)\,dt.
\]

The function `Delta` is even in `u`, and

\[
 \Delta''(u)=\int_{u-v}^{u+v}f''(t)\,dt,\qquad
 f''(t)={1+2t^2\over(1-t^2)^{5/2}}.
\]

Since `f''` is even and increasing in `|t|`, an interval of fixed length
`2v` has its smallest integral when centered at zero.  Hence

\[
 \Delta''(u)\ge2f'(v),\qquad
 \Delta(u,v)\ge2\arcsin(v)+f'(v)u^2,                       \tag{5}
\]

with

\[
 f'(v)={v\over(1-v^2)^{3/2}}.
\]

Substituting `u_e`, `v`, and `|E|=Nd/2` in (4) proves (1).
Endpoint correlations follow by continuity.  This completes the proof.

## 2. Exact application to the outgoing directed half

Let `A` be an order-`n` signing, let `S` be a tournament, and put
`R=A circ S`, a skew signing.  Form the sparse order-`2n` signing

\[
 K_0(A,R)=\begin{pmatrix}A&R\\-R&-A\end{pmatrix}.            \tag{6}
\]

Its only zero off-diagonal entries are the `n` cross-block matching edges,
so its support is `(2n-2)`-regular.  Every state pair can be written
uniquely as `x=s,y=D_Ts`.  The four-label calculation gives

\[
 Q_{K_0}(x,y)
 =4\sum_{\substack{u\in T,\ v\notin T\\u\to v}}
 A_{uv}s_us_v.
\]

Consequently

\[
 \boxed{\Phi(K_0(A,R))=4D_{\to}(A,S).}                     \tag{7}
\]

Put `d=2n-2`.  Direct block multiplication yields

\[
 K_0^2=\begin{pmatrix}P&C\\C&P\end{pmatrix},\qquad
 P=A^2-R^2,\quad C=AR-RA.                                 \tag{8}
\]

Only support edges enter the correction in (1).  The two within-half copies
and the ordered cross edges give exactly

\[
 \Sigma(A,R)=\sum_{i\ne j}
 \left((A^2-R^2)_{ij}^2+(AR-RA)_{ij}^2\right).              \tag{9}
\]

Combining (1), (7), and (9) proves

\[
 \boxed{
 D_{\to}(A,S)\ge
 {n(n-1)\over\pi}\arcsin{1\over\sqrt{2n-2}}
 +{\Sigma(A,R)\over
 16\pi(2n-2)^{5/2}(1-1/(2n-2))^{3/2}}.
 }                                                           \tag{10}
\]

In particular,

\[
 D_{\to}(A,S)\ge
 \left({1\over\sqrt2\,\pi}+o(1)\right)n^{3/2}.            \tag{11}
\]

The desired upper orientation is
`D_to<=M/sqrt(2)+o_Dini(n^(3/2))`.  At the universal lower scale
`M~n^(3/2)/pi`, (11) has exactly the same leading constant.  Thus
`1/sqrt(2)` is analytically sharp in the directed-half geometry; there is no
fixed positive constant margin available.

There is also a recurrence sharpness statement that does not assume the
lower-floor regime.  Define

\[
 U_n=\min_{\Phi(A)=m_n}\min_{S\ {\rm tournament}}D_{\to}(A,S).
\]

Filling the missing perfect matching of `K_0` with arbitrary signs changes
the norm by at most `n`.  Hence

\[
 m_{2n}\le4U_n+n.                                          \tag{11a}
\]

If one had `U_n<=c m_n+o(n^(3/2))` with a fixed `c<1/sqrt(2)` on an
eventual dyadic tail, then
`alpha_(2n)<=sqrt(2)c alpha_n+o(1)` would force the normalized sequence to
zero on that tail, contradicting Proposition 5.2.  Therefore, for every
fixed `n_0>=2`,

\[
 \boxed{
 \limsup_{j\to\infty}{U_{2^jn_0}\over m_{2^jn_0}}
 \ge{1\over\sqrt2}.
 }                                                           \tag{11b}
\]

So the paving constant is forced independently by both the Gaussian
geometry at the lower floor and the global dyadic recurrence.

## 3. Conditional approximate-mate rigidity

Suppose an orientation satisfies

\[
 D_{\to}(A,S)\le {M\over\sqrt2}+r_n.
\]

Rearranging (10) gives the exact necessary condition

\[
 \Sigma(A,R)\le
 16\pi d^{5/2}(1-1/d)^{3/2}
 \left[{M\over\sqrt2}+r_n
 -{n(n-1)\over\pi}\arcsin{1\over\sqrt{2n-2}}\right].      \tag{12}
\]

Writing `M=alpha_n n^(3/2)` and `r_n=epsilon_n n^(3/2)`, this becomes

\[
 {\Sigma(A,R)\over n^4}
 \le64\pi\left(\alpha_n-{1\over\pi}\right)
 +16\pi2^{5/2}\epsilon_n+o(1).                             \tag{13}
\]

Hence along any subsequence with `alpha_n->1/pi` and `epsilon_n->0`, a
successful orientation necessarily obeys

\[
 \Sigma(A,R)=o(n^4).                                       \tag{14}
\]

In words, `R` must be an approximate equal-square commuting mate:
`A^2-R^2` and `AR-RA` have submaximal off-diagonal Frobenius mass.

This rigidity must not be oversold.  With only the currently known interval
`1/pi<=alpha_n<=1/2`, the right side of (13) at `alpha=1/2` is about
`36.54`, whereas the elementary upper bound for `Sigma/n^4` is smaller;
there the conclusion is vacuous.  Also typical `O(sqrt(n))` off-diagonal
entries contribute only `O(sqrt(n))` to (10), below the leading scale.
The rigidity is informative specifically near the universal lower floor or
under a stronger error budget.

## 4. Consequence for the proof search

The outgoing-half target is now squeezed at exactly the critical leading
constant when `M` is near its universal floor.  Any successful general proof
must do one of two genuinely sharp things:

1. construct an orientation attaining the signed-regular arcsine floor and
   control the positive correction; or
2. exploit that an actual optimal `M` sits far enough above `1/pi` to pay for
   the correction, while retaining a Dini-summable amplification error.

The lemma supplies a new theorem and a quantitative structural gate, not the
missing upper construction.  The multiplier-two ray and the original MO
limit remain open.
