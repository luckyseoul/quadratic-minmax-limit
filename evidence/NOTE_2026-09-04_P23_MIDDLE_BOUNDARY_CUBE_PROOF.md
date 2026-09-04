# The p23 mean-46 middle-boundary swap-cube bridge

Classification: **proved theorem** for the local equality cells stated below.
This supplies the middle-boundary bridge in candidate Proposition 15.771;
it does not by itself close that endpoint or residual (ii).

Changed premise: the pending implementation checked parity-half ranks but
did not construct the slice cubes on which those ranks apply. The argument
here gives that construction through every slice point, including both
extreme intersection layers. It uses no graph, orbit, or five-set census.

## The local theorem

Let `B` be a subset of a 23-element set, with
`b=|B| in {6,8,10,12,14,16,18}`. There is no nonnegative integer-valued
function `A` on `J(23,12)` which is the restriction of a real polynomial of
degree at most two and satisfies both

```text
A(X) = |X intersect B| + 1 (mod 2),
E[A] = 1.
```

In the Proposition 15.632 normalization the second condition is precisely
scaled mean `a=2p E[A]=46`.

## Positive quadrature makes every even contact pointwise

Write `r=|X intersect B|`. Its feasible range is

```text
max(0,b-11) <= r <= min(b,12).
```

The conditional average `q(r)=E[A(X) | |X intersect B|=r]` is a polynomial
of degree at most two. Indeed, the averages of a coordinate in `B`, a
coordinate outside `B`, a pair in `B`, a pair outside `B`, and a cross pair
are respectively

```text
r/b,
(12-r)/(23-b),
r(r-1)/(b(b-1)),
(12-r)(11-r)/((23-b)(22-b)),
r(12-r)/(b(23-b)).
```

All displayed denominators are nonzero in the present range. The uniform
slice distribution has first three moments

```text
E[1]=1,  E[r]=12b/23,  E[r^2]=6b(b+1)/23.
```

The following strictly positive quadratures have exactly those moments.
Each vector in the last column is divided by the displayed denominator.

| b | nodes r | denominator | numerators |
|---|---|---:|---|
| 6 | 0,2,4,6 | 92 | 1,39,51,1 |
| 8 | 0,2,4,6,8 | 115 | 3,3,93,13,3 |
| 10 | 0,2,4,6,8,10 | 368 | 3,3,138,218,3,3 |
| 12 | 2,4,6,8,10,12 | 184 | 3,3,154,18,3,3 |
| 14 | 4,6,8,10,12 | 230 | 3,83,138,3,3 |
| 16 | 6,8,10,12 | 23 | 1,18,3,1 |
| 18 | 8,10,12 | 46 | 15,30,1 |

These are all feasible even intersection values, not a proper subset of
them. Since `A` is a nonnegative odd integer on each even layer, `q(r)>=1`
there. Moment matching and `E[A]=1` give

```text
1 = E[A] = sum_r w_r q(r) >= sum_r w_r = 1.
```

Every weight is positive. Consequently every even-layer average equals
one, and every individual value averaged is at least one. Hence

```text
A(X)=1 whenever |X intersect B| is even.                 (1)
```

This is a pointwise equality on the original slice; it does not assume that
`A` depends only on `r` or on the smaller parity side.

## A maximal cross-boundary swap cube through every slice point

Put `S=B` if `b<=11`, and `S=B^c` otherwise. Let `d=|S|=min(b,23-b)`
and let `T=S^c`. In the seven cases above, respectively,

```text
d = 6,8,10,11,9,7,5.
```

Fix any `X in J(23,12)` and put `s=|X intersect S|`. Pair every element
of `S` with a distinct element of `T`, always making the pair contain
exactly one element of `X`. Such a pairing exists by the following two
independent injections:

```text
X intersect S     --> T minus X,
S minus X        --> X intersect T.
```

Their required and available cardinalities are, respectively,

```text
s     <= 11-d+s,
d-s   <= 12-s.
```

Both inequalities hold for every `0<=s<=d`, because `d<=11`; the two
target sets are disjoint. Thus all `d` pairs can be chosen simultaneously.
This includes `s=0` and `s=d`, so no extreme intersection is omitted.

Write the pairs as `(s_i,t_i)`, with `s_i in S` and `t_i in T`. The
unpaired selected elements form a fixed set `R subset T` of size `12-d`.
For `x in {0,1}^d`, define

```text
X(x) = R union {s_i : x_i=1} union {t_i : x_i=0}.
```

All coordinates in this expression are distinct. Every `X(x)` has
cardinality `(12-d)+d=12`, and the chosen original `X` is one of these
vertices. This is a genuine `d`-dimensional slice swap cube: there is no
constraint on its bit patterns. Its incidence coordinates are constants,
`x_i`, or `1-x_i`; consequently `F(x)=A(X(x))-1` has cube degree at most
two after multilinearization.

The boundary intersection is

```text
|X(x) intersect B| = |x|       if S=B,
|X(x) intersect B| = 12-|x|    if S=B^c.
```

Because 12 is even, both cases have the same parity as `|x|`. In particular,
equation (1) says that `F` vanishes on the even half of this cube. Choosing
the bits to index the smaller side is important: indexing the larger-side
choices in the second case would produce a parity-shifted coordinate
formula instead.

## Parity-half injectivity and the contradiction

A degree-at-most-two function on `{0,1}^d` is a linear combination of the
Walsh characters

```text
chi_U(x) = (-1)^(sum_(i in U) x_i),       |U|<=2.
```

On the even subgroup `E={x:|x| even}`, two characters have the same
restriction precisely when their index sets are equal or complementary.
For `d>=5`, distinct sets of size at most two cannot be complementary.
Equivalently, for such index sets the exact character Gram matrix is

```text
sum_(x in E) chi_U(x) chi_V(x) = 2^(d-1) delta_(U,V).
```

Thus restriction to `E` is injective for degree-at-most-two functions.
Since `F` vanishes on `E`, it vanishes on the entire cube. At any odd vertex
this gives `A(X(x))=1`, whereas the prescribed phase-one parity requires
`A(X(x))` to be even. This contradiction proves the theorem.

One such cube would already suffice for the contradiction. The construction
above proves the stronger covering statement requested by the acceptance
gate: every point of `J(23,12)` lies in a suitable cube. No separate
small-side junta theorem or slice-ideal reduction is used for these middle
boundaries.

## Independent exact checks

An independent local rational/integer check on 2026-09-04 verified all seven
quadrature moment identities, all 63 `(b,s)` capacity cases, and the complete
integer Walsh Gram matrices for `d=5,...,11`. For each boundary the layer
cardinality sum

```text
sum_(s=0)^d binom(d,s) binom(23-d,12-s) = binom(23,12) = 1,352,078
```

also passed. This counts the covered domain by intersection sizes; it is not
a slice-function or graph enumeration. The exact Gram ranks, ordered by
dimension `d=5,...,11`, are `16,22,29,37,46,56,67`, with diagonal entries
`16,32,64,128,256,512,1024` and all off-diagonal entries zero. These checks
support the explicit proof above and are not substitutes for it.

The other endpoint obligations are proved and independently reviewed in
`NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md`: the `b=4,20`
general-slice classifications and the phase-zero mass-32 bridge complete
Proposition 15.771. The canonical p23 frontier is now `t>=12,k>=116`.
Residual (ii), E1, `L=1/2`, and the original limit remain open.
