# Proposition 15.751: the generic fourth residual shell is empty

**Status:** proved theorem, with one fixed exhaustive four-cube certificate.

**Scope:** generic branch B at `k=4p+6`, `p=1 mod 4`, `p>=29`.
Together with Propositions 15.739--15.743 and the branch-A/C arithmetic in
15.735, this closes `k=4p+6` for every prime `p>=13`. It does not close any
later layer or residual (ii) globally.

## 1. Exact surviving cell

At `t=3`, branches A and C still have surplus below the number of opposite
directions and are excluded by the argument of 15.735. In branch B equality
forces every opposite direction to have `Q=3` and scaled mean `p+7`. For
`p>=17`, the phase-zero floor table and Proposition 15.688 force `b=0`.
Consequently the surviving local object would be a nonzero nonnegative
integer-valued quadratic

\[
 B:J\left(p,{p+1\over2}\right)\to\mathbb Z_{\ge0},
 \qquad 4p\,\mathbb EB=p+7.                         \tag{1}
\]

Write `H=max B`.

## 2. A dimension-free half-mean cube lemma

If a nonnegative integer-valued quadratic `g` on a Boolean cube has
`E g=1/2`, then

\[
                         \max g\le3.                \tag{2}
\]

Indeed, a nonzero degree-two cube polynomial has support density at least
`1/4`; equality for a nonnegative integral polynomial makes it Boolean.
Every restriction mean lies in `(1/4)Z`, because its multilinear
coefficients are integral finite differences.

Take a minimal-dimensional counterexample with maximum `M>=4` at the origin.
On every coordinate facet through the origin, the mean can be neither
`1/4` (support equality), nor `1/2` (minimality). If it is one, the opposite
facet vanishes and the polynomial factors as `(1-x_i)` times a nonnegative
integral affine function of mean one, whose maximum is at most two. Hence
the two facet means are `3/4` and `1/4`. Every vertex except the origin is
therefore Boolean.

Zero third differences now give, for distinct `i,j,k`,

\[
 M=g(ijk)-g(ij)-g(ik)-g(jk)+g(i)+g(j)+g(k)\le4.    \tag{3}
\]

Equality forces all singleton and triple values to be one and all pair
values to be zero. The constant, linear, and quadratic coefficients are
then `4,-3,2`; on any four-set the value is `4-12+12=4`, contradicting
Booleanity away from the origin. This proves (2). It is sharp:

\[
 g=3-2s+{s\choose2},\qquad s=x_1+x_2+x_3+x_4,
\]

has layer values `3,1,0,0,1` and mean `1/2`.

## 3. Heights at least two are impossible

Choose a maximizing middle set. Proposition 15.688's paired-cube operator
gives

\[
 TB(X)={H+p\mathbb EB\over p+1}
       ={4H+p+7\over4(p+1)}.                       \tag{4}
\]

If `H>=2`, every paired cube through `X` has quarter-integral mean at least
`1/2`; mean `1/4` would make it Boolean despite containing `H`. Thus

\[
                         H\ge {p-5\over4}.          \tag{5}
\]

Writing `p=4r+1`, the exact stabilizer bound in 15.688 is

\[
 H\le{(p+7)(p+3)\over4(p-1)},
 \qquad TB(X)\le {p+7\over2(p-1)}<{3\over4}.       \tag{6}
\]

Therefore some paired cube has mean exactly `1/2`, and (2) gives `H<=3`.
For `p>=29`, (5) gives `H>3`, a contradiction.

## 4. Height one gives a six-coordinate junta

It remains to exclude `H=1`. Complement the slice, put
`k=(p-1)/2`, and write `f:J(p,k)->{0,1}` with

\[
                         \mu=\mathbb Ef={p+7\over4p}. \tag{7}
\]

For a coordinate transposition define

\[
 I_{ij}={1\over4}\Pr[f(X)\ne f(X^{(ij)})].         \tag{8}
\]

For a relevant pair, condition on `x_i=1,x_j=0`. The difference is a
nonzero `{-1,0,1}`-valued affine function on
`J(p-2,(p-3)/2)`. Put `p-2=2r+1`, sort its coefficients, and subtract the
middle coefficient. Differences of coefficients are integral, and the
largest minus smallest `r`-subset sum is

\[
 \sum_{u\ne r+1}|c_u-c_{r+1}|\le2.                \tag{9}
\]

The possibilities are a constant, one unit deviation, one double
deviation, two same-sign unit deviations, or two opposite unit deviations.
Their support densities are respectively at least
`1,r/(2r+1),1,r/(2r+1),(r+1)/(2r+1)`. Hence

\[
 I_{ij}\ge{(p+1)(p-3)\over16p(p-2)}.              \tag{10}
\]

Zero influence is an equivalence relation, since transposition invariances
generate the third transposition. If the largest equivalence class has
complement size `L`, the relevance graph is complete multipartite and has
at least `pL/2` edges.

For the Johnson harmonic decomposition `f=f_0+f_1+f_2`, the transposition
Laplacian has eigenvalue `e(p+1-e)` on level `e`. The normalization in (8)
is essential:

\[
 \sum_{i<j}I_{ij}
 ={1\over2}\sum_{e=1}^2e(p+1-e)\|f_e\|_2^2
 \le(p-1)\mu(1-\mu).                              \tag{11}
\]

A dictator gives `sum I=(p/2) Var(f)`, independently checking the factor
`1/2`. Combining (10)--(11) yields

\[
 L\le {2(p-1)(p-2)(p+7)(3p-7)\over
              p^2(p+1)(p-3)}<7.                  \tag{12}
\]

For the strict inequality, after `p=x+29` the cleared numerator is

\[
 x^4+92x^3+3107x^2+45296x+237300>0.
\]

Thus `L<=6`.

Average a degree-two slice representative over the largest invariant
coordinate class. It has the form

\[
 A_0(x_J)+A_1(x_J)s_C+a_2 {s_C\choose2}.
\]

On the slice substitute `s_C=k-s_J`. Since `L<=6<min(k,p-k)`, every bit
pattern on `J` extends to the slice. The result is a Boolean quadratic on
the full `L`-cube. Every relevant cube coordinate has a nonzero affine
derivative and hence influence at least `1/2`, while Parseval gives total
influence at most `8 Var(f)<=2`. Therefore the cube quadratic depends on at
most four coordinates.

## 5. Fixed four-bit certificate

An exact Mobius transform of all 65,536 Boolean truth tables on four bits
leaves 222 degree-at-most-two tables and fourteen layer profiles. Their
density on `J(p,(p-1)/2)` belongs to exactly

\[
 \left\{0,1,{p-3\over4p},{p+1\over4p},{p-1\over2p},
 {p+1\over2p},{3p-1\over4p},{3(p+1)\over4p}\right\}. \tag{13}
\]

The target `(p+7)/(4p)` lies strictly between `(p+1)/(4p)` and
`(p-1)/(2p)` for `p>9`, so height one is impossible.

The CUDA kernels on Soulkiller's V100 and Orin, and OpenCL kernels on Nuka's
AMD `gfx1201` and Jellyfin's Arc A380, independently returned 222 tables,
the same fourteen-profile histogram, and the same digest

```text
63c9daf2b117b540a5199b1b007cb4e6997ba01704fbc6017efaaa9735859396
```

The scalar exact replay in `src/e1_gmin_m4_prop15751.py` reproduces all
three. The computation classifies only a fixed four-bit object; it does not
enumerate primes, slices, graphs, or residual cells.

## 6. Conclusion and exact remainder

Both `H>=2` and `H=1` are impossible, so generic branch B is empty for every
`p=1 mod 4`, `p>=29`. Existing exact certificates close `p=13,17`; branch
A/C arithmetic closes the other congruence class. Hence

\[
 \boxed{k=4p+6\text{ is impossible for every prime }p\ge13.}
\]

Residual (ii) remains open at critical `p=5,7`, at `p=11,k>=50`, at
`p=13,k=60,u=6` and later p13 layers, at every `p>=17,t>=4` layer, and in
the positive `p=7,z=7` branch. No global predicate is flipped.
