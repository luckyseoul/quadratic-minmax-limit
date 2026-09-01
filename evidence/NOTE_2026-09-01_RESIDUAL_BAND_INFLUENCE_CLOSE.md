# Proposition 15.752: influence rigidity closes a residual-II band

**Status:** proved infinite-family theorem. No finite prime, graph, orbit,
slice, or residual-cell census is used. The only exhaustive ingredient is
Proposition 15.751's fixed classification of 65,536 four-bit truth tables.

## 1. New local theorem

Let `p>=23` be prime and let `B` be a nonzero nonnegative integer-valued
quadratic on `J(p,(p+1)/2)`. Then

```text
4p E[B] != p+9.                                           (1)
```

Put `H=max B`. Paired cubes through a maximizer have average mean

```text
TB(X)=(4H+p+9)/(4(p+1)).                                  (2)
```

If `H>=2`, each cube mean belongs to `(1/4)Z` and is at least `1/2`.
Thus `H>=(p-7)/4`. The exact stabilizer bound gives

```text
H <= (p+9)/4                         (p=3 mod 4),
H <= (p+9)(p+3)/(4(p-1))             (p=1 mod 4).
```

Substitution in (2) makes `TB(X)<3/4` for `p>=23`. Some paired cube
therefore has mean exactly `1/2`. Proposition 15.751's dimension-free cube
theorem gives `H<=3`, whereas `(p-7)/4>=4`. Contradiction.

It remains that `H=1`, so `B` is Boolean with density

```text
mu=(p+9)/(4p).
```

The transposition derivative floor from 15.751 is unchanged. The Johnson
Laplacian identity gives, for the complement size `L` of the largest
zero-influence coordinate class,

```text
L <= 6(p-1)(p-2)(p+9)/(p^2(p+1)) < 7.                    (3)
```

The cleared numerator in (3) is

```text
p^3-29p^2+150p-108,
```

which is 168 at `p=23` and is strictly increasing thereafter. Hence
`L<=6`. Symmetrization extends the slice representative to a Boolean cube
quadratic, and cube influence reduces it to four actual coordinates. The
four-bit catalog has exactly the density list

```text
0, 1, (p-3)/(4p), (p+1)/(4p), (p-1)/(2p), (p+1)/(2p),
(3p-1)/(4p), 3(p+1)/(4p).
```

The target `(p+9)/(4p)` lies strictly between `(p+1)/(4p)` and
`(p-1)/(2p)`, proving (1).

The threshold cannot be lowered to `p=19` by this mechanism. If `R` is a
four-set, `r=|X intersect R|`, and

```text
B=3-2r+binom(r,2),
```

then the layer values are `3,1,0,0,1` and
`4p E[B]=2p-10=p+9=28` at `p=19`. This is only a local quadratic, not a
residual graph, so the `p=19` fifth shell remains a separate global target.

## 2. Exact shell arithmetic

Write

```text
k=4p+2t,  |H|=4p+2t+1,  q=(p-1)/2,  m=q+1.
```

The theorem covers

```text
p=1 mod 4:  p>=29,  4<=t<=q-4=(p-9)/2,
p=3 mod 4:  p>=23,  4<=t<=q-3=(p-7)/2.             (4)
```

Throughout (4), `p^2+1-2|H|>0`; signed transport gives an isolated chart
with `I=0` and every directional odd-fibre count even. In the phase-one
type,

```text
a_d=2u+(p+1)k_d,  sum k_d=m+t-u.                    (5)
```

The exact floor table and the sharp `p-3` integral-lift floor reduce (5) to
the same three baselines as 15.734--15.735:

```text
A: b=2,     P=4, Q_min=3, a_min=8;
B: b=p-1,   P=5, Q_min=2, a_min=6   (p=1 mod 4);
C: b=p-1,   P=3, Q_min=4, a_min=8   (p=3 mod 4).
```

The coefficient congruence and opposite-edge nonnegativity force the
displayed `P` values with quotient zero. The opposite `Q` totals are

```text
A: 4q+t,    B: 3q+t,    C: 5q+t+1.                 (6)
```

The minimum means 8,6,8 are below every nonzero phase-zero fibre floor and
below the `p-3` lift floor, so no `Q_min` direction exists. Raising every
direction once leaves only `t-4`, `t-3`, `t-4` surplus units respectively,
all below `m`. Consequently some direction has

```text
A/C: Q=Q_min+1 and scaled mean p+9,
B:   Q=3       and scaled mean p+7.                  (7)
```

At the means in (7), every nonzero odd-fibre option is an explicit baseline
plus fewer than `p-3` units. Thus the cell has `b=0`, equals `2B`, and
contradicts (1) in A/C or Proposition 15.751 in B.

Therefore every layer in (4) is empty for every boundary size. In
particular,

```text
k=4p+8 is impossible for every prime p>=23.          (8)
```

Together with Propositions 15.734, 15.735, and 15.751, this gives a
contiguous closed band from `t=0` through the endpoint in (4).

## 3. Exact remaining scope

Proposition 15.752 does not close residual II globally. The remaining fifth
shell endpoints are `p=13,k=60,u=6`, `p=17,k=76`, and `p=19,k=84`, in
addition to the prior `p<=11` gate. Later layers beyond (4), the critical
small primes, and the separate positive `p=7,z=7` branch also remain open.

Replay:

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15752.py
PYTHONPATH=src pytest -q tests/test_prop15752.py
```
