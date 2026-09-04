# Sharp local small-mass spectrum

Date: 2026-09-04.

Classification: **proved infinite-family local theorem**, for every prime
`p>=29`, in both congruence classes. This supplies a local input only; it
does not by itself close residual (ii), E1, or the original limit problem.

The changed premise is the two-stage use of the already proved, genuinely
dimension-free cube height bounds at means `1/2` and `3/4`. It replaces
individual exclusions such as `p+5` with a sharp interval theorem. No prime,
graph, orbit, Johnson-cell, or Boolean-cube census is rerun. The existing
fixed four-bit certificate is read with a pinned artifact hash.

## 1. Statement and sharp endpoint

Let `m=(p+1)/2`, and let `C` be a nonzero nonnegative integer-valued
quadratic of degree at most two on `J(p,m)`. Set

```text
M=4p E[C],                 T_p=2p-10.
```

Then

```text
0<M<T_p  implies  C is Boolean and M belongs to {p-3,p+1}.       (1)
```

Both listed masses are attained: `(1-x_i)(1-x_j)` has mass `p-3`, and
`x_i*x_j` has mass `p+1`. Non-Booleanity first becomes possible exactly at
the strict endpoint. For a fixed set `R` of size four or five, put

```text
r=|X intersect R|,         C=3-2r+binom(r,2).                    (2)
```

The values are `3,1,0,0,1` for four coordinates and `3,1,0,0,1,3` for five.
Every pattern extends to the middle slice in our range. In both cases

```text
E[C]=(p-5)/(2p),           M=2p-10,          max C=3.
```

Thus the strict inequality in (1) cannot be replaced by a weak one. No
classification of all endpoint equalities is asserted.

## 2. Every smaller non-Boolean lift is impossible

Write `H=max C`, and suppose `H>=2`. Through a maximizing middle-slice
point `X`, Proposition 15.688 supplies paired Boolean cubes and the exact
averaged-cube operator

```text
T C(X)=(H+p E[C])/(p+1)=(4H+M)/(4(p+1)).                       (3)
```

An integer-valued degree-two polynomial restricted to a Boolean cube has
integer multilinear coefficients, so its mean lies in `(1/4)Z`. A nonzero
degree-two cube polynomial has support density at least `1/4`. At mean
`1/4`, nonnegative integral values are therefore all zero or one. Every
cube through `X` contains `H>=2`, so every such cube has mean at least
`1/2`. Equation (3) gives

```text
H >= [2(p+1)-M]/4 > [2(p+1)-(2p-10)]/4 = 3.                  (4)
```

Proposition 15.751's dimension-free half-mean theorem says that a
nonnegative integral cube quadratic of mean `1/2` has maximum at most
three. Consequently **every** cube through `X` has mean at least `3/4`.
Using (3) again,

```text
H >= [3(p+1)-M]/4 > (p+13)/4 > 6.                            (5)
```

The dimension-free three-quarter-mean theorem in
[Proposition 15.768, Section 4](NOTE_2026-09-04_P1_FIRST_POST_BAND_CLOSE.md)
says that mean `3/4` implies maximum at most six. Its statement and proof
are about arbitrary integral cube quadratics, with no prime-congruence
restriction. Thus every cube through `X` must now have mean at least one.

But the exact stabilizer inequalities from 15.688 give

```text
H <= M/4                             when p=3 mod4,
H <= M(p+3)/(4(p-1))                  when p=1 mod4.
```

Substituting in (3) yields the two cases

```text
T C(X) <= M/[2(p+1)] < 1             when p=3 mod4,
T C(X) <= M/[2(p-1)] < 1             when p=1 mod4.           (6)
```

This contradicts the average of cubes all having mean at least one.
Hence `H=1`, proving Booleanity. Notice that the strict hypothesis was
used in (4); the examples (2) have height exactly three and do not enter
the bootstrap.

## 3. Booleanity gives the exact two-mass spectrum

Complement the slice and write `f:J(p,q)->{0,1}`, where `q=(p-1)/2`.
The transposition-influence calculation of
[Proposition 15.751, Section 4](NOTE_2026-09-01_GENERIC_T3_INFLUENCE_CLOSE.md)
is independent of the special density and of `p mod4`. With normalization

```text
I_ij=(1/4)Pr[f(X)!=f(X^(ij))],
```

a relevant pair has

```text
I_ij >= (p+1)(p-3)/(16p(p-2)).                              (7)
```

Indeed, after conditioning on `x_i=1,x_j=0`, the derivative is a nonzero
`{-1,0,1}`-valued affine function on `J(p-2,(p-3)/2)`. Sort its
coefficients and subtract the median coefficient. The difference between
the largest and smallest fixed-size subset sums is the sum of absolute
deviations from that median, and is at most two. The five possibilities
are the constant case, one unit deviation, one double deviation, two
same-sign unit deviations, and two opposite unit deviations. Their support
densities are at least `1,r/(2r+1),1,r/(2r+1),(r+1)/(2r+1)`, respectively,
where `r=(p-3)/2`. Multiplying by the conditioning probability gives (7).

Zero influence is an equivalence relation. If the largest invariant
coordinate class has complement size `L`, the complete multipartite
relevance graph has at least `pL/2` edges. The Johnson Laplacian identity
gives

```text
sum_(i<j) I_ij <= (p-1)mu(1-mu),       mu=E[f].
```

Therefore

```text
L <= 32(p-1)(p-2)mu(1-mu)/[(p+1)(p-3)]
  <= 8(p-1)(p-2)/[(p+1)(p-3)] < 8.                         (8)
```

The gap in the final inequality is
`8(p-5)/[(p+1)(p-3)]>0`. Thus `L<=7<q`. Average a degree-two representative
over the largest invariant coordinate class `K`. The result has form
`A_0(x_J)+A_1(x_J)s_K+a_2 binom(s_K,2)`. Substitute `s_K=q-s_J`.
Every bit pattern on the at most seven kept coordinates extends to the
slice, so this gives a Boolean degree-at-most-two polynomial on the full
kept-coordinate cube. A relevant cube coordinate has influence at least
`1/2`, while the total cube influence is at most `8 Var(f)<=2`.
Consequently at most four coordinates are actually needed.

The reviewed fixed four-bit certificate has 222 tables and fourteen layer
profiles. Its exact middle-slice density list is

```text
0, 1, (p-3)/(4p), (p+1)/(4p), (p-1)/(2p), (p+1)/(2p),
(3p-1)/(4p), 3(p+1)/(4p).
```

Only `p-3` and `p+1` remain after multiplication by `4p` and restriction
to `0<M<2p-10`. This proves (1). The artifact
`e1_gmin_m4_prop15751.json` is pinned at SHA-256
`b25b0f8896fccfb01c48c92b6266724bf484a0b14f4841655229a647bc2a61a7`.
Its existing table/signature digest is
`63c9daf2b117b540a5199b1b007cb4e6997ba01704fbc6017efaaa9735859396`.
Evaluating its fourteen recorded profiles is not a new catalog.

## 4. Both affine-parity phases, every even boundary

Let `A>=0` be integer-valued and quadratic on `J(p,m)`, with

```text
A(x) = eta + sum_(i in B)x_i (mod2),    eta in {0,1}, |B|=b even.
a=2p E[A],                            0<=a<2p-10.
```

The exact even-boundary floor table in Proposition 15.734 permits only
`b=0,2,p-1` in phase zero and `b=2,p-1` in phase one below this threshold.
Every other boundary has floor at least `2p-6>2p-10`.

At `b=0,eta=0`, pointwise parity gives `A=2C`; then `a=4p E[C]`, so the
possible masses are `0,p-3,p+1` by (1). For the remaining boundaries use
the genuine zero/one pointwise parity minimum `A_0`:

| Boundary | Phase | Pointwise parity minimum | Scaled mass |
|---|---:|---|---|
| `b=2` | 0 | `(x_i-x_j)^2` | `p+1` |
| `b=2` | 1 | `1-(x_i-x_j)^2` | `p-1` |
| `b=p-1` | either | `x_j` if `eta+m` is even, otherwise `1-x_j` | `p+1` or `p-1` |

These are actual parity minima at every slice point: since `A` is
nonnegative and has that parity, `(A-A_0)/2` is globally nonnegative,
integer-valued, and quadratic. Its scaled mass is `a-2p E[A_0]`. If
positive, it is less than `p-9` or `p-11`, hence below the sharp `p-3`
lift floor of 15.688. Thus no positive lift exists. This subtraction does
not use the non-Boolean complement-triple baseline or an omitted layer.

The exact local sets are consequently:

| Prime class | Phase zero | Phase one |
|---|---|---|
| `p=1 mod4` | `{0,p-3,p-1,p+1}` | `{p-1,p+1}` |
| `p=3 mod4` | `{0,p-3,p+1}` | `{p-1}` |

Every listed mass has the displayed local polynomial or a twice-Boolean
example from Section 1. This makes no claim that the examples glue to a
single residual graph. In either prime class the two-phase union is
`{0,p-3,p-1,p+1}`.

## 5. Implementation and verification scope

The implementation is `src/e1_gmin_m4_small_mass_spectrum.py`. Its public
APIs are `small_mass_spectrum(p)`, `affine_parity_small_mass_spectrum(p)`,
and `local_mass_exclusion(p,mass)`. The last API certifies only the open
interval `p+1<mass<2p-10`; outside it both `proved` and `excluded` are
false. In particular it does not exclude either attained endpoint.

Cube, stabilizer, fixed-catalog, floor-table, and pointwise-baseline proof
dependencies are checked before a result can be marked proved. The focused
tests include failure injection, strict endpoint preservation, both prime
classes and parity phases, artifact equality, and a guard that forbids
calling the old four-bit enumerator. Tests and evidence generation are
offloaded; no controller test suite is used.

From the staged remote checkout:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. /home/nick/.venvs/mo-exact/bin/python src/e1_gmin_m4_small_mass_spectrum.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. /home/nick/.venvs/mo-exact/bin/python -m pytest -o addopts= -q -n 0 tests/test_small_mass_spectrum.py
```
