# The p=37 next all-finite endpoint is impossible

Date: 2026-08-28. This is Proposition 15.680. It closes the all-finite
boundary `p=37,s=30`; it does **not** close the same boundary at
`p=17,19,23,29,31,41`, later all-finite sizes, residual (ii), R1, global
QVAR, Type I, E(1), or the limit.

## 1. Exact boundary ledger

Put `p=37`, `P=p+1=38`, and `m=19`. The second even all-finite size above
`3(p-1)/4=27` is `s=30`. The exact phase-one quotient/floor replay has one
row:

```text
u_1=18,  quotient sum=1,  profile b=2 (18 times), b=30 (once),
D_1=504.
```

The pair-deficit budget is `s(s-1)=870`. The exact phase-zero replay leaves
only

```text
u_0   D_0   pair slack   minimizing b-profile
 2    328       38       0^10 2^1 30^8
 3    330       36       0^11     30^8
 4    358        8       0^11 2^1 30^7
 5    360        6       0^12     30^7.
```

Residue zero has deficit `388`, hence total deficit `892>870`; residue one
is infeasible. The interior minimum is strictly increasing after `u=2`, and
already `u=6` again has deficit `388`. The final four residues have deficit
at least `(m-4)s=450`. Thus the four displayed rows are exhaustive, not
just the minima of a selected family.

For each retained `u`, the quotient sum is `m-u<m`, so some phase-zero
direction has quotient zero. Its scaled mean is `2u<=10`. Every positive
even fibre count has phase-zero floor at least `P=38`, hence this direction
has `b=0` and

```text
A=2B,   2u=4p E[B],
```

where `B` is a nonzero nonnegative integer-valued quadratic on `J(37,19)`.
Proposition 15.642 gives nonzero lift cost at least ten and immediately
excludes `u=2,3,4`. It is sharp at the remaining `u=5` row, so a new
argument is required.

## 2. A paired-cube lemma for Boolean quadratics

The following statement holds for every odd `p>=5`:

> If `f` is a nonzero Boolean polynomial of degree at most two on
> `J(p,(p+1)/2)`, then
> `E[f] >= (p-3)/(4p)`.

Write `p=2m-1`. Fix a middle set `X`. Choose one point of `X` uniformly as
a leftover, pair the other `m-1` points bijectively with the complement,
and choose one endpoint of each pair. Together with the leftover, these
choices form a Boolean `(m-1)`-cube through `X`.

Let `T` average first over such paired cubes through `X` and then over the
cube. Directly on the monomial basis through degree two, with
`rho=1/(2m)=1/(p+1)`, one has

```text
T(1)=1,
T(x_i)=1/2+rho*x_i,
T(x_i*x_j)=1/4+rho*x_i*x_j.
```

The middle-slice means of `x_i` and `x_i*x_j` show equivalently that

```text
T(f)=rho*f+(1-rho)E[f]
```

for every quadratic `f`. If `f(X)=1`, its restriction to every paired cube
through `X` is nonzero. The elementary cube polynomial-distance lemma says
that a nonzero degree-two polynomial on a cube has support density at least
`1/4`. Since `f` is Boolean,

```text
rho+(1-rho)E[f] = T(f)(X) >= 1/4.
```

Solving gives

```text
E[f] >= (1/4-rho)/(1-rho) = (p-3)/(4p).
```

This restriction proof is finite and self-contained; it does not assume a
classification of Boolean slice quadratics or a conjectural exact minimum
support theorem.

## 3. The mass-ten lift cannot exist

Suppose now that the `u=5` row exists. Then

```text
E[B]=5/74.
```

At `p=37=4*9+1`, the stabilizer identity of Proposition 15.642 has endpoint
weight `9/370`. Applying it at any point where `B=h` gives

```text
5/74 >= (9/370)h,
```

so every value of `B` lies in `{0,1,2}`.

The exact degree-two slice-distance floor is

```text
Pr(B!=0) >= C(33,17)/C(37,19) = 171/2590.
```

Consequently the density of points where `B=2` is at most

```text
5/74 - 171/2590 = 2/1295.
```

If the value two occurred, `B(B-1)` would be a nonzero degree-four
polynomial supported exactly on those points. The same exact finite
slice-distance lemma gives

```text
Pr(B=2) >= C(29,15)/C(37,19)
        = 1938/441595
        = 2/1295 + 1256/441595,
```

a contradiction. Thus `B` is Boolean. The paired-cube lemma now gives

```text
E[B] >= 17/74 > 5/74,
```

the final contradiction. Therefore `u=5` is impossible, all four pair rows
are excluded, and the boundary `p=37,s=30` does not exist.

## 4. Literature and duplicate check

The exact degree-two and degree-four finite support inputs are both direct
instances of Lemma 2 in Amireddy--Behera--Srinivasan--Sudan [47].
Filmus's slice-junta threshold [51] and the paired constructions of
Kiermaier--Mannaert--Wassermann [52] are close context, but neither supplies
the paired-cube transition identity or this mass-ten exclusion. The latter
paper explicitly treats the general minimum-size problem as open and gives
constructions; no classification claim is imported here.

Exact OEIS API searches for the distinctive normalized blocks
`171,175,1938,1256` and `328,330,358,360,504`, plus the individual larger
values `441595` and `1194096750`, all returned `null`. This is a
duplicate/context check, not a sequence-submission or priority claim.

## 5. Reproduction

- theorem and exact arithmetic replay: `src/e1_gmin_m4_prop15680.py`;
- machine-readable record: `evidence/e1_gmin_m4_prop15680.json`;
- regression tests: `tests/test_prop15680.py`.
