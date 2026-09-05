# Full row optimality does not bound tilted moments for arbitrary cavities

2026-09-05. Analytic counterexample; no computation is used.

## Scope

This is a counterexample for **arbitrary symmetric cavity measures**. It is
strictly positive on the entire cube, invariant under coordinate permutations
and global spin reversal, and satisfies every signed-row replacement
inequality. Nevertheless, its row-tilted second and fourth field moments
diverge at critical row scale.

No realization as the cavity law of the actual quadratic Ising host is proved.
In particular, this is **not** a counterexample signing, host, or pressure
comparison. Its consequence is that a proof of bounded tilted moments must
use additional structure of the actual cavity law or global host optimizer,
not only the full row-replacement inequalities.

## Proposition and construction

Fix `c > 0`. Let `d` range over even integers with `d >= max(2, c^2)`, and put

\[
u=\frac{c}{\sqrt d},\qquad H=ud=c\sqrt d,\qquad
p=\frac{u}{4\cosh H}.
\]

In particular, `0 < u <= 1` and `0 < p <= 1/4`. On `{ -1, 1 }^d`, let

- `B_d` be the uniform measure on the balanced shell `sum_i x_i = 0`;
- `E_d = (delta_{(1,...,1)} + delta_{(-1,...,-1)})/2`;
- `U_d` be the uniform measure on the entire cube.

Define

\[
\nu_0=(1-p)B_d+pE_d,\qquad
\nu=\frac12\nu_0+\frac12U_d.
\]

The measure `nu` is even, permutation invariant, and has strictly positive
mass at every cube point. For a row `b in { -1, 1 }^d`, write

\[
Z_b(u)=\int\cosh(u b\cdot x)\,d\nu(x).
\]

The row `a=(1,...,1)` is the unique global minimizer up to global sign. Define
its field and tilted measure by

\[
g(x)=u\sum_{i=1}^d x_i,\qquad
d\mu(x)=\frac{\cosh g(x)}{Z_a(u)}\,d\nu(x).
\]

With the finite constant

\[
K_c=\frac{5/4+e^{c^2/2}}2,
\]

one has

\[
\mathbb E_\mu g^2\ge\frac{c^3\sqrt d}{8K_c},\qquad
\mathbb E_\mu g^4\ge\frac{c^5d^{3/2}}{8K_c}.
\]

Thus neither moment is bounded by a constant depending only on `c`, even
under the complete family of row-replacement inequalities and full support.

## Proof of global row minimization

Since `Z_b=Z_{-b}`, it suffices to consider a row `b` with `k` negative
coordinates, where `0 <= k <= d/2`. On the balanced shell, if `I` is its
negative-coordinate set, then

\[
b\cdot x=-2\sum_{i\in I}x_i.
\]

Under `B_d`, distinct coordinates have covariance `-1/(d-1)` and each
coordinate has mean zero and variance one. Therefore

\[
\mathbb E_{B_d}(b\cdot x)^2
=\frac{4k(d-k)}{d-1}.
\]

Using `cosh y >= 1+y^2/2`, and `2(d-k)/(d-1) >= 1`, gives

\[
\mathbb E_{B_d}\cosh(u b\cdot x)-1
\ge\frac{2u^2k(d-k)}{d-1}\ge u^2k.
\]

The extreme-point component instead contributes `cosh(H-2uk)`, as compared
with `cosh H` for `a`. Since `0 <= H-2uk <= H`, the mean-value bound gives

\[
0\le\cosh H-\cosh(H-2uk)\le 2uk\cosh H.
\]

Its possible loss is consequently at most

\[
p\,[\cosh H-\cosh(H-2uk)]\le\frac{u^2k}{2}.
\]

Combining the balanced-shell gain with this loss yields

\[
\begin{aligned}
\mathbb E_{\nu_0}\cosh(u b\cdot x)
-\mathbb E_{\nu_0}\cosh(u a\cdot x)
&\ge (1-p)u^2k-\frac{u^2k}{2}\\
&\ge\frac{u^2k}{4}.
\end{aligned}
\]

The uniform-cube contribution is `(cosh u)^d`, independently of `b`, so

\[
Z_b(u)-Z_a(u)\ge\frac{u^2k}{8}.
\]

For an arbitrary row with `k` negative coordinates, the same assertion holds
with `k` replaced by `min(k,d-k)`. The gap is strictly positive unless
`b=a` or `b=-a`. This proves the asserted full row optimality.

## Tilted moment lower bounds

For the all-plus row,

\[
\mathbb E_{\nu_0}\cosh g=1-p+p\cosh H
=1-p+\frac u4\le\frac54.
\]

Also, `log cosh u <= u^2/2` implies

\[
\mathbb E_{U_d}\cosh g=(\cosh u)^d\le e^{c^2/2}.
\]

Hence `Z_a(u) <= K_c`. The rare extreme-point component of `nu` alone
contributes total unnormalized tilted mass

\[
\frac p2\cosh H=\frac u8
\]

at points where `|g|=H`. For every `q > 0`, this gives

\[
\mathbb E_\mu |g|^q
\ge\frac{uH^q}{8K_c}
=\frac{c^{q+1}d^{(q-1)/2}}{8K_c}.
\]

The cases `q=2` and `q=4` prove the proposition.

## The full noise hierarchy also holds

For clarity, define `Z_a(v)=integral cosh(v sum_i x_i) dnu(x)` at any real
amplitude `v`, with the same fixed measure `nu` constructed above. Let
`rho in [-1,1]`, and choose independent row signs `R_i` with mean `rho`.
Full row optimality gives

\[
Z_a(u)\le\mathbb E_R Z_R(u).
\]

Set

\[
v=\operatorname{arctanh}(\rho\tanh u),\qquad
C=\cosh u\sqrt{1-\rho^2\tanh^2u}.
\]

For `x_i in {-1,1}`, the identity

\[
\cosh u+\rho x_i\sinh u=C e^{v x_i}
\]

and its counterpart with `x_i` negated show

\[
\mathbb E_R\cosh\left(u\sum_i R_ix_i\right)
=C^d\cosh\left(v\sum_i x_i\right).
\]

Consequently this example satisfies, for every `rho in [-1,1]`,

\[
\boxed{\quad
Z_a(u)\le
\left[\cosh u\sqrt{1-\rho^2\tanh^2u}\right]^d
Z_a\!\left(\operatorname{arctanh}(\rho\tanh u)\right).
\quad}
\]

Indeed, every inequality obtained by averaging any collection of row
replacements holds, since each replacement holds separately. Thus the
complete replacement hierarchy, including this noise family, is compatible
with the divergent tilted moments above.

## Compatibility with the local radial bound

The example also retains the useful first radial bound. Flipping only row
coordinate `i` gives

\[
\frac{Z_{a^{(i)}}(u)}{Z_a(u)}
=\cosh(2u)-\sinh(2u)\,
\mathbb E_\mu[x_i\tanh g]\ge1.
\]

Hence `E_mu[x_i tanh g] <= tanh u`, and summing gives

\[
0\le\mathbb E_\mu[g\tanh g]
\le du\tanh u\le du^2=c^2.
\]

The bounded first radial quantity and the divergent second and fourth
moments therefore coexist in a full-support cavity with strict global row
optimality. Transferring a moment bound to the actual Ising problem requires
an additional actual-host structural input; that input is not supplied or
refuted here.
