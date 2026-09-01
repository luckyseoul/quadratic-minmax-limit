# Conference obstruction to every fixed-temperature free-energy target

**Status:** proved counterexample theorem (Proposition 6.9). No finite census
is used. The old statement that `c=3` remains viable is superseded.

Let

\[
 Q_C(x)=\sum_{i<j}C_{ij}x_ix_j={1\over2}x^TCx,
 \qquad x\in\{\mathord\pm1\}^n,
\]

where `C` is a symmetric conference signing:

\[
 C=C^T,\quad C_{ii}=0,\quad C_{ij}\in\{\mathord\pm1\},
 \quad C^2=(n-1)I.
\]

## 1. A sharp half-projection Laplace bound

Let `P` be an orthogonal projection on `R^n`, of rank `n/2`, with
`P_ii=1/2` for every `i`. For every `t>0`,

\[
 \boxed{\mathbb E_x e^{-t x^TPx}
 \le\left({1+e^{-2t}\over2}\right)^{n/2}.}                       \tag{1}
\]

To prove this, choose `U` with `UU^T=I`, `U^TU=P`, and write `u_i` for its
columns. The Gaussian Fourier identity gives

\[
 \mathbb E_xe^{-t x^TPx}
 =(4\pi t)^{-n/4}\int e^{-\|\xi\|^2/(4t)}
   \prod_i\cos\langle u_i,\xi\rangle\,d\xi.                    \tag{2}
\]

Take absolute values and set `v_i=sqrt(2)u_i`. Then every `v_i` is a unit
vector and

\[
 \sum_i{1\over2}v_iv_i^T=I.                                    \tag{3}
\]

Apply the geometric rank-one Brascamp--Lieb inequality with weights `1/2`
to

\[
 f(s)=e^{-s^2/(4t)}\cos^2(s/\sqrt2).
\]

The integral in (2), after absolute values, is at most

\[
 \left(\int_{\mathbb R}f(s)\,ds\right)^{n/2}.
\]

Since

\[
 \int_{\mathbb R}f(s)\,ds
 =\sqrt{4\pi t}\,{1+e^{-2t}\over2},                              \tag{4}
\]

the normalization in (2) cancels and proves (1). Equality holds for a
direct sum of rank-one projections onto `(1,1)/sqrt(2)`, so the bound is
sharp among projections with these parameters.

## 2. Conference free-energy upper bound

Put `lambda=sqrt(n-1)` and

\[
 P_\pm={1\over2}\left(I\pm{C\over\lambda}\right).
\]

Both are rank-`n/2` orthogonal projections with diagonal `1/2`, and

\[
 Q_C(x)=\lambda\left({n\over2}-x^TP_-x\right),\qquad
 -Q_C(x)=\lambda\left({n\over2}-x^TP_+x\right).                 \tag{5}
\]

Applying (1) separately to the two signs gives the exact domination

\[
 \boxed{\mathbb E_x e^{\mathord\pm\beta Q_C(x)}
 \le\cosh(\beta\sqrt{n-1})^{n/2}.}                              \tag{6}
\]

Therefore

\[
 \boxed{\mathbb E_x\cosh\left({c\over\sqrt n}Q_C(x)\right)
 \le\cosh\left(c\sqrt{1-{1\over n}}\right)^{n/2}.}             \tag{7}
\]

Along any sequence of symmetric conference signings `C_n` whose orders tend
to infinity, for every fixed `c>0`,

\[
 \limsup_{n\to\infty}{1\over n}
 \log\mathbb E_x\cosh\left({c\over\sqrt n}Q_{C_n}(x)\right)
 \le {1\over2}\log\cosh c < {c\over2}.                        \tag{8}
\]

At `c=3`, the coefficient on the right is

\[
 {1\over2}\log\cosh3=1.154664252\ldots<1.5,
\]

so the proposed lower bound fails by at least

\[
 {3-\log\cosh3\over2}\,n-o(n)
 =0.345335747\ldots n-o(n).                                    \tag{9}
\]

## 3. Even-Eulerian consequence

For the signed even-Eulerian polynomial,

\[
 \mathbb E_x\cosh(\beta Q_C(x))
 = (\cosh\beta)^{\binom n2}P_C(\tanh\beta).
\]

Since
`binom(n,2) log cosh(c/sqrt(n))=c^2 n/4+o(n)`, equation (8) gives

\[
 \limsup {1\over n}\log P_{C_n}(\tanh(c/\sqrt n))
 \le {1\over2}\log\cosh c-{c^2\over4}.                         \tag{10}
\]

This is strictly below the formerly proposed sufficient lower bound

\[
 \left({c\over2}-{c^2\over4}\right)n-o(n)
\]

for every fixed `c>0`. Symmetric Paley conference matrices exist at order
`q+1` for every prime `q=1 mod 4`, so this is an infinite counterfamily, not
an isolated finite-order exception.

The result does not rule out a temperature `c=c_n` tending to infinity. Put
`t_n=c_n sqrt(1-1/n)`. The exact conference loss relative to `c_n n/2` is

\[
 {n\over2}(c_n-\log\cosh t_n)
 ={\log2\over2}n+{c_n\over2(1+\sqrt{1-1/n})}+o(n)
 =O(n+c_n)=o(c_n n).
\]

A growing-temperature criterion with a suitably uniform error therefore
remains logically possible.

Primary analytic input: H. J. Brascamp and E. H. Lieb, *Best constants in
Young's inequality, its converse, and its generalization to more than three
functions*, Advances in Mathematics 20 (1976), 151--173,
doi:10.1016/0001-8708(76)90184-5.
