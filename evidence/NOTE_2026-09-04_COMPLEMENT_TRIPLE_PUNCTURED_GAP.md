# Complement-triple punctured gap and the first new equality

Date: 2026-09-04.

Classification: **proved all-prime local gap and equality theorem** in the
range below. This repairs the previously implicit complement-triple
floor-plus-two step of Proposition 15.770 and identifies a genuinely new
equality at floor plus four. It is not an endpoint or global residual-(ii)
closure by itself.

## 1. Exact theorem and the missing premise

Let `p>=29` be odd, put `m=(p+1)/2`, and fix a three-set `C`. Let `A` be
a nonnegative integer-valued quadratic function on `J(p,m)` with parity

```text
A(X) = r (mod 2),         r=|X intersect C|.
```

Suppose

```text
2p E[A] = 2p-6+delta,         0<=delta<=4.                 (1)
```

Then exactly the following possibilities occur:

* `delta=0`: `A=(r-2)^2`.
* `delta=4`: choose a pair `{i,j} subset C` and let `k` be its remaining
  point. Then

  ```text
  A=(1-x_i-x_j)^2+(1-x_k)
   =2-x_i-x_j-x_k+2x_i x_j.                              (2)
  ```

There are exactly three labeled forms in the second case. In particular,
scaled mean `2p-4` is impossible, whereas scaled mean `2p-2` is attained.

For `p=1 mod 4`, `m` is odd. A phase-one boundary `B=C^c` has parity
`|X intersect B|+1=m-r+1=r (mod2)`, so this is precisely the
`b=p-3` hard-cell theorem needed at the generic post-band frontier.

The baseline `A_0=(r-2)^2` is not the pointwise parity minimum at `r=0`:
it equals four there, while parity only requires a nonnegative even value.
Thus `(A-A_0)/2` can be negative on that layer. Proposition 15.688 cannot
be applied to that difference as a globally nonnegative lift without a
separate argument. The proof below supplies the omitted argument at excess
two and classifies, rather than incorrectly excludes, excess four.

## 2. A neighboring-slice lower bound without a prime hypothesis

Let `N>=6` be even and let `f` be a nonzero nonnegative integer-valued
quadratic on `J(N,N/2+1)` or `J(N,N/2-1)`. Then

```text
E[f] >= (N-2)(N-4) / (4N(N-1)).                          (3)
```

For the upper slice, restrict each coordinate to one. Every section is a
nonnegative integral quadratic on `J(N-1,N/2)`, an odd middle slice.
At most two of the `N` sections can vanish identically. Otherwise choose
three coordinates whose one-sections vanish and a support point `X` of
`f`. That point avoids all three coordinates. Exchange three distinct
selected points of `X` with the three forbidden coordinates to make a
three-dimensional swap cube. Every nonzero cube vertex contains a forbidden
coordinate and has value zero. A degree-two cube polynomial has zero third
difference, so its origin also has value zero, contradicting the choice
of `X`.

For the lower slice, instead restrict each coordinate to zero. If three
zero-sections vanished, a support point would contain all three
coordinates. Swapping them with three unselected points gives the same
contradiction. Each nonzero section is on `J(N-1,N/2-1)`; complementing
its variables makes it an odd middle slice of weight `N/2`.

Proposition 15.688 applies to every odd order `n>=5`, not only primes.
Consequently each nonzero section in either argument has mean at least

```text
(n-3)/(4n) = (N-4)/(4(N-1)),         n=N-1.
```

The average of the `N` section means equals `E[f]`. For example, on a
fixed-weight `k` slice,

```text
(1/N) sum_a E[f | x_a=1]
 = E[f*(sum_a x_a)/k]
 = E[f],
```

and the analogous identity uses `N-k` for zero-sections. At least `N-2`
sections are nonzero, proving (3).

No slice enumeration or support-classification theorem is used here.

## 3. Positive quadrature kills the singleton and triple differences

Put `L=(A-A_0)/2`, where `A_0=(r-2)^2`. This is an integer-valued
quadratic. On every point with `r=1,2,3`, it is nonnegative, since `A_0`
equals the pointwise parity minimum on those layers. On `r=0` we only
know `L>=-2`.

For a specified `S subset C`, write `mu_S` for the mean of `L` on the
outside slice obtained by setting `X intersect C=S`. The conditional mean
averaged uniformly over all small-side patterns `S` with `|S|=r`, and
then over their outside completions, is a quadratic in `r`. No assertion
that the individual `mu_S` depend only on `|S|` is being made. The exact
positive quadrature at nodes `1,2,3` has weights

```text
3(p-3)/(4p),       3/p,       (p-3)/(4p).
```

These match the first three hypergeometric moments, and the baseline has
scaled mean `2p-6`. Equation (1) therefore gives the exact identity

```text
delta = 4p E[L]
      = (p-3)*(sum_(|S|=1) mu_S + mu_C)
        +4*sum_(|S|=2) mu_S.                            (4)
```

All terms on the right are nonnegative. Set `N=p-3`. A singleton section
of `L` lives on `J(N,N/2+1)`, and the triple section lives on
`J(N,N/2-1)`. If any of these four sections were nonzero, (3) would make
its contribution to (4) at least

```text
(p-5)(p-7)/(4(p-4)) > 4.                                (5)
```

For `p>=29`, the strict gap is explicit:

```text
(p-5)(p-7)/(4(p-4)) -4
 = ((p-29)^2+30(p-29)+128)/(4(p-4)) >0.
```

At the endpoint the bound is `132/25`. Thus `delta<=4` forces

```text
L=0 pointwise on every r=1 and r=3 configuration.         (6)
```

The difference is still not assumed nonnegative on `r=0`.

## 4. Globalizing the remaining difference

Use the three complemented small-side bits `w_i=1-x_i` and outside bits
`y_j`. Write `s=sum_i w_i`, `k_0=m-3`, and `N=p-3`. The full slice now
has the identity

```text
sum_j y_j = k_0+s.
```

The contact layers (6) are `s=0,2`. Start with an arbitrary quadratic
representative

```text
L(w,y)=q(y)+sum_i a_i w_i+sum_(i<l) a_il w_i w_l
                    +sum_(i,j) C_ij w_i y_j.
```

At `s=0`, `q` vanishes on `J(N,k_0)`. The elementary fixed-weight
quadratic kernel lemma proved in Section 2 of
[the small-boundary note](NOTE_2026-09-04_P23_SMALL_BOUNDARY_EQUALITY_PROOF.md)
applies because `2<=k_0<=N-2`. It writes

```text
q(y)=(sum_j y_j-k_0)*(lambda+sum_j mu_j y_j).
```

Substituting `sum y-k_0=s` eliminates the pure outside part, changing
the linear coefficients to `a_i+lambda` and the mixed coefficients to
`D_ij=C_ij+mu_j`. Fix two small bits to one. The `s=2` contact says that
`D_ij+D_lj` is constant in `j`, by linear exchange on the outside
`(k_0+2)`-slice. For any outside indices `j,j'`, the three numbers
`D_ij-D_ij'` have all pairwise sums zero and therefore are all zero.
Thus `D_ij=t_i` is independent of `j`.

The mixed term is now

```text
(sum_i t_i w_i)*(k_0+s)
 =(k_0+1)sum_i t_i w_i+sum_(i<l)(t_i+t_l)w_i w_l.
```

Hence `L` depends only on the three small-side bits on the entire original
slice. Every bit pattern extends, since the outside weights run from
`m-3` to `m`. This proves global small-side dependence without assuming it.

## 5. The exact residual values and signed offsets

Return to the original `x` bits on `C`. Let `a_ij` be the value of `L`
on the pair `{i,j}`. These are nonnegative integers. The singleton and
triple values vanish by (6). The zero third difference on the three-cube
then gives

```text
L(000) = -sum_(i<j) a_ij.
```

Equation (4) reduces to

```text
delta = 4 sum_(i<j) a_ij.                               (7)
```

For `0<=delta<=4`, either every pair value is zero, yielding `L=0` and
the original baseline, or exactly one pair value is one. In the latter
case the eight values give exactly (2). Conversely (2) is a sum of two
Boolean nonnegative quadratic/linear terms, has parity `r`, and has mean
`(p-1)/p`, so it attains `delta=4` for every odd `p` in the stated range.

On the original signed slice `z_i=2x_i-1`, its target is

```text
3+2A = 5+z_i z_j-z_k.
```

This is an integral signed representative with coefficient offset
`constant+sum(linear coefficients)=4`. The `delta=0` complement-triple
baseline has target `5-sum_C z_i+sum_(i<j in C) z_i z_j` and offset two,
as in Proposition 15.768.

## 6. Scope at the next generic endpoint

For `p=1 mod4`, the prospective next layer is `t=q-1`, where
`q=(p-1)/2`. The `b=p-3` low-row means `2p-6`, `2p-4`, and `2p-2`
are now respectively the old offset-two equality, impossible, and the new
offset-four pair-plus-complement-literal equality. In particular, the
floor-plus-two exclusion needed in Proposition 15.770 is justified by this
theorem, not by subtracting the non-Boolean baseline and asserting global
nonnegativity.

At this target a quotient-one mean `2p` would require residue `u=q`.
That residue has quotient sum `q<m`, so it already contains a quotient-zero
baseline row. Classifying every possible high row of mean `2p` is not
needed if its common-row baseline ledger closes; this note does not claim
such a classification. Common-graph normalization and the remaining
opposite-cell contradictions are separate obligations.

## 7. The p=1 mod4 local lift mass p-1 is excluded

Let `p=1 mod4`, `p>=29`, and let `L` be a nonzero nonnegative integral
quadratic with `4p E[L]=p-1`. Suppose first that its maximum `H` is at
least two. The sharp paired-cube and stabilizer inequalities from
Proposition 15.688 give

```text
H >= (2(p+1)-(p-1))/4 = (p+3)/4,
H <= (p-1)(p+3)/(4(p-1)) = (p+3)/4.
```

Thus `H=(p+3)/4>=8`. Every cube through a maximizer has mean at least
one half: cube means are quarter-integral, and a nonnegative integral
quadratic of cube mean one quarter is Boolean. The average over the
maximizing cubes is

```text
(H+p E[L])/(p+1)
=((p+3)/4+(p-1)/4)/(p+1)=1/2.
```

Hence every such cube has mean exactly one half, and Proposition 15.751
bounds its maximum by three, a contradiction.

The remaining height-one case is Boolean, of density `(p-1)/(4p)`.
The established corrected Johnson influence bound gives at most

```text
J = 2(p-1)^2(p-2)(3p+1)/(p^2(p+1)(p-3))
```

coordinates outside the largest zero-influence class. In fact `J<6`
uniformly. If `D=p^2(p+1)(p-3)` and `x=p-29>=0`, its exact gap is

```text
6D-2(p-1)^2(p-2)(3p+1)
=10p^3-40p^2+2p+4
=10x^3+830x^2+22912x+210312 >0.
```

Thus at most five slice coordinates remain, fewer than either slice side;
all their patterns extend. Symmetrization supplies a degree-two cube
representative, and the existing cube influence bound leaves at most four
active coordinates. The fixed fourteen profiles of Proposition 15.751 have
exactly the following density values:

```text
0, 1, (p-3)/(4p), (p+1)/(4p), (p-1)/(2p), (p+1)/(2p),
(3p-1)/(4p), 3(p+1)/(4p).
```

The target `(p-1)/(4p)` lies strictly between the first two positive values
`(p-3)/(4p)` and `(p+1)/(4p)`, and is absent. This proves the local
mass-`p-1` exclusion. In particular, the equality of the older
height-at-least-two floor with `p-1` is not treated as an exclusion by
itself; the maximizing-cube equality argument above is essential.

## 8. The p=1 mod4 local lift mass p+11 is excluded

Now suppose `4p E[L]=p+11`, still with `p=1 mod4`, `p>=29`. If
`H>=2`, the same paired-cube inequality gives

```text
H >= (p-9)/4 >=5.
```

The stabilizer bound gives `H<=(p+11)(p+3)/(4(p-1))`, so the average
maximizing-cube mean is at most `(p+11)/(2(p-1))`. Its gap below three
quarters is exactly

```text
3/4-(p+11)/(2(p-1)) = (p-25)/(4(p-1)) >0.
```

Each such cube again has quarter-integral mean at least one half. Some
maximizing cube therefore has mean exactly one half, forcing `H<=3` by
Proposition 15.751, a contradiction.

For a Boolean lift, the corrected Johnson bound is now

```text
J = 2(p-1)(p-2)(p+11)(3p-11)/(p^2(p+1)(p-3)) <8.
```

The exact all-parameter gap, with the same `D` and `x`, is

```text
8D-2(p-1)(p-2)(p+11)(3p-11)
=2p^4-42p^3+338p^2-814p+484
=2x^4+190x^3+6776x^2+107936x+651360 >0.
```

At `p=29`, the bound is `76608/10933`, slightly larger than seven;
claiming a six-coordinate junta there would be incorrect. The sufficient
bound is at most seven coordinates, still fewer than `q>=14`, and the
same cube influence step leaves at most four active coordinates. Density
`(p+11)/(4p)` lies strictly between `(p+1)/(4p)` and `(p-1)/(2p)`, so
the same fixed catalog excludes it. No new catalog is constructed.

The executable package is `src/e1_gmin_m4_complement_triple_gap.py`; its
three public local certificates are `complement_triple_gap_certificate`,
`p1_p_minus_one_local_exclusion`, and `p1_p_plus_eleven_local_exclusion`.
Exact coefficientwise translation proves the generic positive polynomial
inequalities; endpoint instantiations are checks, not a sample-based proof.
