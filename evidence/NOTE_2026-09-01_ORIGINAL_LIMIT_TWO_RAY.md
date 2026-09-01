# Original limit: two-ray convergence reduction

**Status:** proved theorem (conditional convergence criterion); the two
amplification estimates themselves remain open.

## Scope correction

MathOverflow 413935 asks whether

\[
\alpha_n=\frac{m_n}{n^{3/2}}
\]

converges.  Identifying the value is not part of the stated requirement.
Consequently the Paley E(1) program is one value-specific sufficient route,
not the unique acceptance gate for the original problem.

Throughout this note,
`Q_A(x)=sum_(i<j) a_ij x_i x_j` and `Phi(A)=max_x |Q_A(x)|`, so
`m_n=min_A Phi(A)`.

## The two-ray theorem

Put

\[
H(n)=m_n^{2/3},\qquad h(n)=\frac{H(n)}n=\alpha_n^{2/3}.
\]

The function `H` is nondecreasing because `m_n` is nondecreasing.  Let
`eta(n) >= 0`, and define

\[
\eta^*(N)=\sup_{u\ge N}\eta(u),\qquad
E(N)=\sum_{j\ge0}\eta^*(2^jN).
\]

Suppose `E(N) -> 0` and, for all sufficiently large `n`,

\[
H(2n)\le 2H(n)+2n\eta(n),\qquad
H(3n)\le 3H(n)+3n\eta(n).                    \tag{1}
\]

Then `h(n)`, and hence `alpha_n`, converges.

Indeed, each multiplication by `r in {2,3}` changes `h` by at most
`eta(current size)`.  Along any word in 2 and 3 the size before step `j` is
at least `2^j n`.  Thus, uniformly in `a,b >= 0`,

\[
h(2^a3^b n)\le h(n)+E(n).                       \tag{2}
\]

The sorted multiplicative semigroup

\[
\mathcal S=\{2^a3^b:a,b\ge0\}
\]

has consecutive ratios tending to one.  On taking logarithms this follows
from the irrationality of `log(2)/log(3)`: the irrational rotation by
`log(2)` modulo `log(3)` is dense, so a finite prefix is an arbitrarily fine
net; after adding a sufficiently large nonnegative multiple of `log(3)`,
every sufficiently large logarithmic interval of the prescribed length
contains `a log(2)+b log(3)`.

Fix `n` beyond that threshold.  For every sufficiently large `N`, choose
`s in S` with

\[
N\le sn=(1+o_N(1))N.
\]

Monotonicity and (2) give

\[
h(N)\le \frac{sn}{N}h(sn)
     \le (1+o_N(1))(h(n)+E(n)).
\]

Therefore

\[
\limsup_N h(N)\le h(n)+E(n).
\]

Let `n` tend to infinity along a liminf subsequence.  This proves
`limsup h <= liminf h`, and hence convergence.

Power-saving errors are sufficient, but unnecessarily strong.  For example,
it is enough in the original scale to prove, for either fixed `epsilon>0`,

\[
m_{rn}\le r^{3/2}m_n+
O\!\left(\frac{n^{3/2}}{(\log n)^{1+\epsilon}}\right),
\qquad r=2,3.                                      \tag{3}
\]

Indeed, `m_n=Theta(n^(3/2))` converts (3) into (1) with
`eta(n)=O((log n)^(-1-epsilon))`, whose dyadic tail is
`E(n)=O((log n)^(-epsilon))`.  The older power-saving target
`O(n^(3/2-delta))` is one special case.

It is also enough to prove a Dini-summable doubling estimate together with
the matching single `1:2` split

\[
H(3n)\le H(n)+H(2n)+\text{Dini-summable error},
\]

because the doubling estimate then supplies the tripling estimate.  Thus the
original problem no longer requires a composition theorem for every pair of
orders.

## Why one ray is not enough

For sufficiently small positive `epsilon`, the continuous function

\[
G(x)=x\left(1+\epsilon\sin(2\pi\log_2x)\right)
\]

is increasing, has bounded increments on the integers, and satisfies
`G(2x)=2G(x)` exactly.  Nevertheless `G(n)/n` has different subsequential
limits (take integers nearest `2^(k+theta)` for two different phases
`theta`).  Hence even exact doubling, monotonicity, and strong local
regularity do not force convergence.  A second multiplicatively independent
scale is load-bearing.

## Correction: ordinary two-block composition is still live

For a block signing `S=[[A,B],[B^T,C]]`, pairing `y` with `-y` gives the exact
identity

\[
\Phi(S)=\max_{x,y}\bigl(|Q_A(x)+Q_C(y)|+|x^TBy|\bigr). \tag{4}
\]

Thus `Phi(S)>=||B||_(infinity -> 1)>=(sqrt(2/pi)-o(1))n^(3/2)` for square
order-`n` cross blocks.  The old Section 10 incorrectly combined this lower
bound with a triangle upper bound to claim a universal additive cost at
least `0.282`.  A lower bound on a term in a nonsharp upper estimate cannot
prove that.  Equation (4) prevents pointwise cancellation, but the internal
and cross maxima need not occur at the same state.  Coupled two-block design
therefore remains live.

## Exact multiplier-two Hadamard state

For the natural orthogonal subclass, take a Hadamard `2x2` sign block on
every cloud pair and one arbitrary internal sign per cloud.  Write

\[
z_i=s_i(1,(-1)^{t_i}),\qquad s_i\in\{\pm1\},\quad t_i\in\mathbb F_2.
\]

Every oriented block has a unique form

\[
\frac12z_i^TB_{ij}z_j
=s_is_j(-1)^{t_it_j+\alpha_{ij}t_i+
                 \beta_{ij}t_j+\gamma_{ij}}.         \tag{5}
\]

The coefficient of `t_i t_j` is forced because the product of the four block
entries is `-1`; a vertex-factorable table has product `+1`.  For fixed `t`,
let `C_t` be the induced order-`n` signing in (5), and put
`K=max_t Phi(C_t)`.  The internal matching has magnitude at most `n`, so

\[
2K-n\le \Phi(S)\le2K+n.                             \tag{6}
\]

If `G_2(n)` minimizes `K` over all such frames, the exact Hadamard-doubling
target is

\[
G_2(n)\le\sqrt2\,m_n+n^{3/2}\omega(n),              \tag{7}
\]

where `omega>=0` has a vanishing dyadic Dini tail.  Orientation does not
remove the local four-state curvature; this does not rule out a nonlocal
coincidence among the induced signings.

There is a sharper necessary condition.  Let `G=C_0`, `J=C_1`, and write
`P(A)=max Q_A`, `N(A)=-min Q_A`.  For every `T subset [n]`, independently
maximize the `J[T]` and `G[T^c]` forms and then flip all spins in `T` to
reverse only the cross term.  This proves

\[
K\ge\max\{P(J[T])+P(G[T^c]),\ N(J[T])+N(G[T^c])\}.  \tag{8}
\]

Thus (7) requires one frame whose two endpoint signings obey a simultaneous
hereditary bound for every cut.  For that fixed frame, its mixed cross terms
cannot repair a cut that violates (8).  A different orientation changes the
endpoint signings, so finding a frame that passes all cuts remains open.

## New live gate

First attack (7)--(8), or bypass the Hadamard subclass with another coupled
doubling construction; then prove the `1:2` split.  Logarithmic
Dini-summable error is enough, so polynomial saving should not be imposed as
an unnecessary burden.  A finite Paley residue census does not address this
gate.
