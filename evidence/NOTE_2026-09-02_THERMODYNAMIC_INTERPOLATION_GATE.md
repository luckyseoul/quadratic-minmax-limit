# Critical-temperature interpolation: exact gate and graphon no-go

**Status:** proved finite-temperature reduction and proved counterexample to
ordinary graphon/right-convergence at the critical `n^(-1/2)` coupling scale.
This does **not** prove convergence of `m_n/n^(3/2)`.

This note tests a genuinely different route from the multiplier-two
orientation problem.  It does not use finite primes, a residue census, fixed
physical temperature, or a fixed-`c` linear free-energy lower bound.

Write

\[
 Q_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
 \Phi(A)=\max_x|Q_A(x)|,
 \qquad \alpha_n={m_n\over n^{3/2}}.
\]

## 1. The correctly scaled minimum pressure really would settle the problem

For `c>0`, define the critical symmetric pressure

\[
 s_n(c)={1\over n}\min_A\log
 \mathbb E_x\cosh\!\left({c\over\sqrt n}Q_A(x)\right).       \tag{1}
\]

For every `n` and `c>0`,

\[
 {s_n(c)\over c}\le \alpha_n
 \le {s_n(c)\over c}+{\log2\over c}.                         \tag{2}
\]

Indeed, for every signing `A`, the average in (1) is at most
`exp(c Phi(A)/sqrt(n))`.  Conversely, the two states `x_*,-x_*`
at which `|Q_A|=Phi(A)` contribute at least
`exp(c Phi(A)/sqrt(n))/2^n` to the normalized average.  Taking
logarithms and then minimizing proves (2).

Consequently, if `s_n(c)` converges as `n -> infinity` for every `c` in any
unbounded set, then `alpha_n` converges.  More generally,

\[
 \limsup_n\alpha_n-\liminf_n\alpha_n
 \le {\log2\over c}
   +{\limsup_n s_n(c)-\liminf_n s_n(c)\over c}.               \tag{3}
\]

Thus it would suffice that the pressure oscillation be `o(c)` along some
`c -> infinity`.  This is a valid thermodynamic reformulation; unlike the
disproved target in Proposition 6.9, it does not assert a particular value
for the pressure.

## 2. What block interpolation gives after the disorder minimum

It is useful first to keep the raw inverse temperature common.  Put

\[
 Z_A(\beta)=\mathbb E_x e^{\beta Q_A(x)},\qquad
 p_n(\beta)=\min_A\log Z_A(\beta).
\]

For all positive `n,m` and every `beta>=0`, one has the exact sandwich

\[
 p_n(\beta)+p_m(\beta)
 \le p_{n+m}(\beta)
 \le p_n(\beta)+p_m(\beta)+nm\log\cosh\beta.                 \tag{4}
\]

For the lower bound, restrict any order-`n+m` signing to diagonal blocks
`B,D`, with cross block `C`.  Pair `y` with `-y`:

\[
\begin{aligned}
 Z_A(\beta)
 &=\mathbb E_{x,y}e^{\beta(Q_B(x)+Q_D(y))}
                  \cosh(\beta x^TCy)\\
 &\ge Z_B(\beta)Z_D(\beta).
\end{aligned}                                                 \tag{5}
\]

For the upper bound, start with minimizers `B,D` and choose the entries of
`C` independently and uniformly.  For every fixed `x,y`,

\[
 \mathbb E_Ce^{\beta x^TCy}=\cosh(\beta)^{nm},
\]

so some integral `C` realizes the upper inequality in (4).

There is an analogous one-sided subadditivity statement for the symmetric
pressure.  If

\[
 c_A=\mathbb E\cosh(\beta Q_A),\qquad
 d_A=\mathbb E\sinh(\beta Q_A),
\]

and

\[
 q_n^{\rm sym}(\beta)=\min_A\log\mathbb E\cosh(\beta Q_A),
\]

replace one diagonal block by its negative so that `d_B d_D<=0`.  Random
cross completion then gives

\[
 q_{n+m}^{\rm sym}(\beta)
 \le q_n^{\rm sym}(\beta)+q_m^{\rm sym}(\beta)
      +nm\log\cosh\beta.                                    \tag{6}
\]

After subtracting `binom(n,2) log cosh(beta)`, (6) is genuinely subadditive
at a **common raw** `beta`.  This is the exact scope of the available Fekete
argument.

On the critical diagonal, however, set `N=2n` and
`beta=c/sqrt(2n)`.  If

\[
 f_n^+(c)={1\over n}p_n(c/\sqrt n),
\]

then (4) becomes

\[
 f_n^+(c/\sqrt2)\le f_{2n}^+(c)
 \le f_n^+(c/\sqrt2)
      +{n\over2}\log\cosh\!\left({c\over\sqrt{2n}}\right). \tag{7}
\]

The last term tends to `c^2/8`, not zero.  Moreover, the block pressure is at
`c/sqrt(2)`, not at `c`.  Thus common-temperature super/subadditivity leaves
an extensive `c^2 n/4+o(n)` total gap and does not yield the required
`o(n)` (let alone Dini-summable) cross-order control.

## 3. Ordinary dense graph convergence loses the critical pressure

The failure is not merely that a theorem has been applied with an awkward
normalization.  At the critical scale, the pressure is not continuous in the
ordinary signed cut metric.

Represent an order-`n` signing by its signed step kernel and write

\[
 \|A\|_\square={1\over n^2}\max_{S,T}
 \left|\sum_{i\in S,j\in T}a_{ij}\right|.                    \tag{8}
\]

Fix any `0<c<1`.  There are two deterministic signing sequences, both with
`||A_n||_square -> 0`, whose critical symmetric pressures are separated by
a fixed positive asymptotic gap.

### Conference sequence

Along an infinite Paley symmetric-conference sequence `C_n`,
`||C_n||_op=sqrt(n-1)`, hence

\[
 \|C_n\|_\square\le {\sqrt{n-1}\over n}\longrightarrow0.    \tag{9}
\]

Proposition 6.9 gives, for both signs,

\[
 \mathbb E_xe^{\mathord\pm cQ_{C_n}(x)/\sqrt n}
 \le\cosh(c\sqrt{1-1/n})^{n/2}.
\]

Therefore

\[
 \limsup_n {1\over n}\log\mathbb E_x
 \cosh(cQ_{C_n}(x)/\sqrt n)
 \le {1\over2}\log\cosh c.                                  \tag{10}
\]

### Quasirandom independent-sign sequence

Let `J` have independent uniform signs above the diagonal and put

\[
 X_J=\mathbb E_xe^{cQ_J(x)/\sqrt n},\qquad
 M_n=\mathbb E_JX_J=\cosh(c/\sqrt n)^{\binom n2}.             \tag{11}
\]

For two spin states put `z_i=x_i y_i`, `R=sum_i z_i`, and
`u=tanh^2(c/sqrt(n))`.  The numbers of edges on which the two edge
characters agree and disagree are

\[
 N_s={n^2+R^2-2n\over4},\qquad
 N_o={n^2-R^2\over4}.
\]

Independence of the edge signs gives the exact second-moment identity

\[
 {\mathbb E_JX_J^2\over M_n^2}
 =\mathbb E_z(1+u)^{N_s}(1-u)^{N_o}.                         \tag{12}
\]

This ratio is uniformly bounded for `c<1`.  Indeed,

\[
 (1+u)^{N_s}(1-u)^{N_o}
 \le \exp\!\left({c^2R^2\over2n}\right),
\]

and the Gaussian integral identity together with
`cosh(t)<=exp(t^2/2)` gives

\[
 \mathbb E_z\exp\!\left({c^2R^2\over2n}\right)
 \le {1\over\sqrt{1-c^2}}.                                  \tag{13}
\]

Paley--Zygmund now gives a fixed positive probability that
`X_J>=M_n/2`.  If
`Y_J=E_x cosh(cQ_J(x)/sqrt(n))`, then `E_J Y_J=M_n`, so Markov's
inequality also gives `Y_J<=K M_n` outside a set of arbitrarily small fixed
probability by choosing fixed `K` large enough.  Finally, a Hoeffding union
bound over the at most `4^n` pairs `S,T` shows that
`||J||_square=O(n^(-1/2))` with probability tending to one.  The three
events therefore intersect for every sufficiently large `n`.  Selecting one
such deterministic `J_n` gives

\[
 {M_n\over4}\le Y_{J_n}\le K M_n,
 \qquad \|J_n\|_\square=O(n^{-1/2}),
\]

and hence

\[
 \lim_n {1\over n}\log Y_{J_n}={c^2\over4}.                 \tag{14}
\]

The random-sequence limit and the conference-sequence limsup are separated
because

\[
 {c^2\over4}-{1\over2}\log\cosh c>0\qquad(c>0),             \tag{15}
\]

whose derivative is `(c-tanh(c))/2>0`.

Thus both sequences converge to the zero signed graphon (equivalently, to
edge density `1/2` in the unsigned encoding), but their critical pressures
are different.  The dense left/right-convergence theorem of Borgs--Chayes--
Lovasz--Sos--Vesztergombi concerns fixed `1/n`-scale interactions; it does
not extend to this moderate-deviation `1/sqrt(n)` functional.  Sparse
right-convergence does not repair the issue: these weighted complete systems
have diverging absolute weighted degree, and the signed cut-metric
counterexample remains.

## 4. Why Guerra--Toninelli does not survive the outer minimum

Guerra--Toninelli proves subadditivity after a quenched average over
independent centered disorder.  Here the disorder is the optimization
variable.  Their limit gives only

\[
 \min_A\log Z_A\le \mathbb E_J\log Z_J,                       \tag{16}
\]

with no reverse comparison.  The missing event is a lower tail over the
`2^(binom(n,2))` signings.  Ordinary bounded-difference concentration is at
speed only `n` for order-`n` deviations of `log Z`: changing one edge changes
`log Z` by at most `2c/sqrt(n)`.  It therefore cannot union-bound over a
disorder space of exponential size `exp(Theta(n^2))`.  This agrees with the
separate direct-limit audit: what would be needed is a joint Bernoulli lower-
tail principle at speed `n^2`, not the usual quenched thermodynamic limit.

Primary references for the scopes being separated here are Guerra--Toninelli,
[The Thermodynamic Limit in Mean Field Spin Glass Models](https://arxiv.org/abs/cond-mat/0204280),
and Borgs--Chayes--Lovasz--Sos--Vesztergombi,
[Convergent sequences of dense graphs II](https://doi.org/10.4007/annals.2012.176.1.2).
For sparse interpolation/right-convergence, see Gamarnik,
[Right-convergence of sparse random graphs](https://arxiv.org/abs/1202.3123),
and Borgs--Chayes--Gamarnik,
[Convergent sequences of sparse graphs: A large deviations approach](https://arxiv.org/abs/1302.4615).

## 5. Verdict

The minimum critical pressure (1) is a clean direct target and (2)--(3) show
exactly why it would solve the MO problem.  What is now ruled out is obtaining
its convergence from any of the following without a new premise:

1. common-raw-temperature Fekete interpolation;
2. the usual quenched Guerra--Toninelli limit;
3. ordinary dense graphon/left/right convergence; or
4. sparse right-convergence under its standard bounded-activity hypotheses.

A live thermodynamic continuation must control the optimized Bernoulli lower
tail at speed `n^2`, or introduce a stronger second-order limit object that
distinguishes conference signings from independent quasirandom signings.
