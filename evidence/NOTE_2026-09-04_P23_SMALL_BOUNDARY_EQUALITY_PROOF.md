# The small-boundary equality bridge for Proposition 15.771

Date: 2026-09-04.

Status: **proved general-slice equality classification** for the phase-one
mean-46 cells with `b=4` and `b=20`. This note supplies one of the three
explicit bridges requested by
`NOTE_2026-09-04_P23_THIRD_POST_BAND_REVIEW.md`. It does not, by itself,
promote Proposition 15.771 to reviewed endpoint closure.

The argument applies to every degree-at-most-two function on `J(23,12)`
with the stated integrality, parity, nonnegativity, and mean. It does not
assume in advance that the function depends only on the boundary, and it
uses no graph, orbit, five-set, or full-slice enumeration.

## 1. Setup and pointwise contacts

Let `X` range over the twelve-subsets of a 23-element coordinate set. Write
`A(X)` for a nonnegative integer-valued function having a multilinear
polynomial representative of degree at most two. The mean-46 condition is

```text
E[A]=1,                  2*23 E[A]=46.
```

If `B` is the phase-one parity boundary, the parity condition is

```text
A(X) = 1+|X intersect B|  (mod 2).
```

For `b=4`, set `S=B` and `d=4`. For `b=20`, set `S=B^c` and `d=3`.
In the latter case, `|X intersect B|=12-|X intersect S|`, so the same
parity formula holds with `S` in place of `B`. The small-side coordinates
always mean the indicators of membership in the original twelve-set `X`;
we are not complementing the slice variables.

Put `r=|X intersect S|` and `O=S^c`. Averaging any degree-two representative
over permutations within `S` and within `O` gives a polynomial `a(r)` of
degree at most two: the averages of its terms are constants, multiples of
`r`, `12-r`, `r(r-1)`, `(12-r)(11-r)`, and `r(12-r)`.

The hypergeometric law of `r` has

```text
E[r]   = 12d/23,
E[r^2] = 6d(d+1)/23.
```

The following strictly positive probability weights have those same
first two moments:

| Small side | Nodes | Weights |
|---|---|---|
| `d=4` | `0,2,4` | `2/23,18/23,3/23` |
| `d=3` | `0,2` | `5/23,18/23` |

For `b=20`, the small-side nodes `0,2` are the original-boundary nodes
`12,10`, respectively. These are precisely the contact quadratures checked
in `mean_46_contact_quadratures()`.

At an even node, every value of `A` is a nonnegative odd integer and hence
at least one. Therefore

```text
1=E[A]=sum_r weight(r)*a(r) >= sum_r weight(r)=1.
```

Each contact weight is positive, so every contact average equals one.
Every summand of that average is at least one, so the conclusion is
pointwise:

```text
A(X)=1 whenever r=0,2       for d=3,
A(X)=1 whenever r=0,2,4     for d=4.                 (1)
```

## 2. An elementary fixed-weight quadratic kernel lemma

We need an exact identity for polynomial representatives, not only an
averaging assertion. Let

```text
q(y)=c+sum_j b_j y_j+sum_(i<j) b_ij y_i y_j,
```

where `y` has `N` Boolean coordinates and `2<=k<=N-2`. If `q` vanishes
on every `k`-subset, there are real numbers `lambda,mu_1,...,mu_N` such that
as functions on the entire Boolean cube,

```text
q(y)=(sum_j y_j-k)*(lambda+sum_j mu_j y_j).          (2)
```

Products in this identity are reduced using `y_j^2=y_j`.

To prove it, fix distinct `i,j`. For every `(k-1)`-subset `R` avoiding
`i,j`, subtraction of the two zero values gives

```text
0=q(R+i)-q(R+j)
 =b_i-b_j+sum_(a in R)(b_ia-b_ja).
```

Exchange one element of `R`. The inequalities on `k` allow this for any
two coordinates outside `i,j`; hence `b_ia-b_ja` is independent of such
`a`. Thus the off-diagonal pair array satisfies

```text
b_ia-b_ja=b_ib-b_jb          for four distinct i,j,a,b.       (3)
```

These equations imply `b_ij=mu_i+mu_j`. Explicitly, on three distinct
anchor coordinates take

```text
mu_1=(b_12+b_13-b_23)/2,
mu_2=(b_12+b_23-b_13)/2,
mu_3=(b_13+b_23-b_12)/2,
mu_i=b_1i-mu_1                 for i>=4.
```

Equation (3) first proves the claimed formula for pairs meeting anchors
`2,3` and then for every remaining pair. With `s=sum_j y_j`, the pair part
of `q` consequently equals `(s-1)*sum_j mu_j y_j`. On `s=k`, its linear
coefficients are `b_j+(k-1)mu_j`. A linear function constant on all
`k`-subsets has all coefficients equal, by a single-coordinate exchange.
Call their common value `lambda`. The constant value zero now gives
`c=-k lambda`, proving (2).

## 3. Removing all dependence on the outside coordinates

Write the small-side indicators as `x_i`, the outside indicators as `y_j`,
and set `F=A-1`. Start with an arbitrary degree-two representative

```text
F(x,y)=q(y)+sum_i a_i x_i+sum_(i<l) a_il x_i x_l
               +sum_(i,j) C_ij x_i y_j.                     (4)
```

There are `N=19` outside coordinates for `d=4`, and `N=20` for `d=3`.
The `r=0` part of (1) says that `q` vanishes on `J(N,12)`; both values
of `N` satisfy the kernel lemma. Applying (2), and then using the total
weight identity `sum_j y_j=12-sum_i x_i`, changes (4), on the full slice,
to a representative of the form

```text
F(x,y)=sum_i alpha_i x_i+sum_(i<l) alpha_il x_i x_l
                    +sum_(i,j) D_ij x_i y_j,                (5)
```

with no pure outside part. Explicitly,
`alpha_i=a_i-lambda`, `alpha_il=a_il`, and `D_ij=C_ij-mu_j`.

Fix any pair of distinct small-side indices `i,l`, and set precisely
those two small-side indicators to one. The `r=2` part of (1) says that

```text
alpha_i+alpha_l+alpha_il+sum_j(D_ij+D_lj)y_j=0
```

on every ten-subset of `O`. Since `0<10<N`, the same elementary linear
exchange argument implies that `D_ij+D_lj` is independent of `j`.

For two outside coordinates `j,j'`, put `v_i=D_ij-D_ij'`. We have
`v_i+v_l=0` for every distinct pair in the small side. There are at least
three small-side indices. For three of them, `v_1=-v_2`, `v_1=-v_3`, and
`v_2+v_3=0` imply all three are zero; pairing with one of these shows every
remaining `v_i=0`. Therefore every row `D_ij` is constant in `j`; write
its value as `t_i`.

The mixed term in (5) now equals

```text
(sum_i t_i x_i)*(12-sum_i x_i)
 =11 sum_i t_i x_i-sum_(i<l)(t_i+t_l)x_i x_l.
```

Thus (5) is a degree-two polynomial in `x` alone on every point of the
original slice. Every Boolean pattern on `S` extends to a twelve-subset,
because the required outside sizes range from `8` to `12` for `d=4`,
or from `9` to `12` for `d=3`, within the available `19` or `20` positions.
Consequently the resulting small-cube polynomial retains nonnegativity,
integrality, and parity on its entire cube.

This proves the previously implicit globalization. It uses only the
pointwise contacts `r=0,2`; the additional `r=4` contact will classify the
four-coordinate polynomial.

## 4. The complete small-cube equality types

Let `e_i` be the singleton value of `A` at small-side coordinate `i`.
By parity and nonnegativity, every `e_i` is a nonnegative even integer.
The values one at the empty pattern and at every pair determine the
quadratic polynomial uniquely:

```text
A(x)=1+sum_i(e_i-1)x_i
        +sum_(i<j)(2-e_i-e_j)x_i x_j.                       (6)
```

In particular, its value at any triple `T` is

```text
A(1_T)=4-sum_(i in T)e_i.                                  (7)
```

For `d=4`, its full-set value is `9-2 sum_i e_i`. The `r=4` contact forces
`sum_i e_i=4`. Equation (7) then says that each triple value equals the
singleton value at the omitted coordinate. The only possibilities are
the permutations of `4000` and `2200`: respectively four and six labeled
forms, totaling ten. Conversely, (6) for any of these assignments is
nonnegative, has the required parity on every pattern, and equals one at
every contact. The quadrature gives mean one.

For `d=3`, put `E=sum_i e_i`. Equation (7), nonnegativity, and parity give
`A(111)=4-E>=0`, so `E` is `0`, `2`, or `4`. The complete list is

| Singleton type | Triple value | Number of labeled forms |
|---|---:|---:|
| `000` | `4` | `1` |
| `200` | `2` | `3` |
| `220` | `0` | `3` |
| `400` | `0` | `3` |

Conversely, each of these ten forms satisfies the entire cube's
nonnegativity and parity conditions, and the contact quadrature gives
mean one. There are no additional forms with hidden outside dependence,
by Section 3.

## 5. Integral signed representatives and offsets

Use the original-slice signed coordinates `z_i=2x_i-1` and expand the
target `3+2A` from (6). Its pair coefficient at `{i,j}` is

```text
c_ij=1-(e_i+e_j)/2,
```

which is integral. For `d=4`, the constraint `E=4` gives

```text
3+2A=5+sum_(i<j)c_ij z_i z_j.
```

Every signed linear coefficient is zero, and the coefficient offset
`constant+sum(linear coefficients)` is exactly `5` for both equality
types.

For `d=3`, the signed expansion is

```text
3+2A=5+sum_i [1-(E-e_i)/2] z_i+sum_(i<j)c_ij z_i z_j.
```

Its constant, linear, and pair coefficients are all integers. Its offset
is `5+3-E=8-E`, giving the complete list

| Singleton type | Signed linear coefficients in that order | Offset |
|---|---|---:|
| `000` | `(1,1,1)` | `8` |
| `200` | `(1,0,0)` | `6` |
| `220` | `(0,0,-1)` | `4` |
| `400` | `(1,-1,-1)` | `4` |

These are explicit integral representatives on the original signed
slice, suitable for the already established isolated-chart coefficient
comparison. In particular, the offset-eight `000;4` endpoint is present;
it cannot be discarded by assuming that only the positive singleton
types occur.

The ten four-coordinate forms, ten three-coordinate forms, and offsets
`5` and `{4,6,8}` are exactly those encoded by
`mean_46_small_support_equality_catalog()` in
`src/e1_gmin_m4_prop15771.py`. The argument above proves their
exhaustiveness for general slice quadratics, not merely for the listed
small-cube tables.
