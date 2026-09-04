# Proposition 15.769: first `p=3 mod 4` layer beyond 15.752

**Status:** proved infinite-family theorem using only exact symbolic floor
identities and Proposition 15.751's fixed four-bit catalog.  No prime, graph,
orbit, slice, or residual-cell census is used.

## Statement

For every prime `p>=31` with `p=3 mod 4`, residual (ii) is empty at

```text
q=(p-1)/2,  m=q+1,  t=q-2=(p-5)/2,
k=4p+2t=5p-5,  |E(H)|=5p-4.
```

This is exactly the first layer beyond Proposition 15.752's
`p=3 mod 4` interval, which ends at `t=q-3`.  The parameterized local theorem
does not treat the same layer at `p=23`; Proposition 15.769's exceptional
companion now closes that one endpoint by the equality/moment argument in
`NOTE_2026-09-04_P23_POST_BAND_MOMENT_CLOSE.md`.  Later layers, the remaining
critical small primes, and residual (ii) globally remain open.

## 1. Exhaustive residue ledger

An isolated projective vertex exists because

```text
p^2+1-2(5p-4)=p^2-10p+9>0.
```

After signed transport it to infinity, write the hard phase-one means as

```text
a_L=2u+(p+1)k_L,       sum_L k_L=m+t-u.
```

At `t=m-3`, exact phase-one floors and Proposition 15.688's sharp `p-3`
lift floor leave only the residues `u=m-3` and `u=m-1`.

- For `u<m-3`, any quotient-one cell is a nonzero lift of mass below `p-3`.
- At `u=m-3`, every quotient is one and every hard mean is `2p-4`.  This is
  the one genuinely new branch.
- At `u=m-2`, a quotient-zero cell is forced below the phase-one floor.
- At `u=m-1`, the two endpoint equality types are the old XNOR and
  complementary-literal branches.  Their common-row ledgers again force a
  `p+9` opposite cell, so Proposition 15.752 already excludes them.

Only the `u=m-3` branch needs a new argument below.

## 2. Classifying the sharp `p-3` lift

At `u=m-3`, a hard cell's parity baseline is one of

```text
A_0=(1-x_i-x_j)^2,       3+2A_0=4+z_i z_j,
A_0=1-x_j,               3+2A_0=4-z_j.
```

Both baselines have scaled mean `p-1`.  Therefore
`B=(A-A_0)/2` is a nonzero nonnegative integral quadratic with

```text
4p E[B]=p-3.
```

Equality in Proposition 15.688 forces `B` to be Boolean.  The exact Johnson
transposition-influence calculation bounds the complement of the largest
zero-influence coordinate class by

```text
6(p-1)(p-2)/p^2 < 6.
```

Thus at most five slice coordinates remain.  Since `5<q` for `p>=31`, every
pattern extends to `J(p,q)`.  Symmetrization extends it to a Boolean
degree-at-most-two cube function, and cube influence leaves at most four
active coordinates.

The fixed four-bit catalog has exactly ten tables at density `(p-3)/(4p)`:
six selected-pair tables and four all-equal-triple tables.  Complementing
back to `J(p,m)` gives exactly the omitted-pair and all-equal-triple lifts.

## 3. Common-row rigidity and the `p+13` cell

Combining either lift with either parity baseline gives coefficient offsets

```text
2, 3, 4, 5.
```

For every hard direction, the local coefficient sum and the common graph row
give

```text
sum_a q_L(a)=pP_L-3p-(2p-4)=hT-P_L.                 (1)
```

Equation (1) first forces one common hard parallel count `P`.  The slice
kernel congruence says `q` divides `P-offset`.  Since the edge bound gives
`P<=9<q+2`, it follows that `P=offset`; different families cannot mix.
For each of the four cases,

```text
hT=(p+1)P-5p+4,
a(Q)=(p+1)Q+hT-3p.                                   (2)
```

The row `Q=8-P` has scaled mean 12.  It is below both the least nonzero
phase-zero odd-fibre floor `p+1` and the sharp nonzero lift floor `p-3`, so
it is impossible.  Exact opposite-edge accounting then gives

```text
sum_L (Q_L-(9-P))=m-9.
```

Consequently at least nine directions have `Q=9-P`, and (2) gives the same
scaled mean `p+13` in all four families.  At this mean, each nonzero
phase-zero baseline would leave a positive lift of mass 12, again below
`p-3`.  Hence the cell has `b=0`, is `A=2C`, and requires

```text
C:J(p,m)->Z_{>=0},  C nonzero,  4p E[C]=p+13.       (3)
```

## 4. Uniform local `p+13` theorem

Let `H=max C`.  If `H>=2`, paired-cube averaging and the `p=3 mod 4`
stabilizer bound give

```text
H >= (p-11)/4 > 3,
T C(X) <= (p+13)/(2(p+1)) < 3/4.
```

Quarter integrality therefore produces a paired cube of mean exactly `1/2`.
Proposition 15.751's dimension-free half-mean theorem says `H<=3`, a
contradiction.

If `H=1`, then `C` is Boolean of density `(p+13)/(4p)`.  The Johnson
influence bound for the complement `L` of the largest zero-influence class is

```text
L <= 2(p-1)(p-2)(p+13)(3p-13) /
       (p^2(p+1)(p-3)) < 8.                           (4)
```

After writing `p=x+31`, the cleared gap in (4), up to a positive factor, is

```text
x^4+99x^3+3670x^2+60728x+381824>0.
```

Thus `L<=7<q`; every pattern extends to the complementary middle slice, and
cube influence again reduces to four active coordinates.  The fixed catalog
has possible densities

```text
0, 1, (p-3)/(4p), (p+1)/(4p), (p-1)/(2p), (p+1)/(2p),
(3p-1)/(4p), 3(p+1)/(4p).
```

The target `(p+13)/(4p)` lies strictly between `(p+1)/(4p)` and
`(p-1)/(2p)` for `p>=31`, so it is absent.  This excludes (3) and closes the
only new branch.

## 5. Sharp threshold and scope

At `p=23`, if `r=|X intersect R|` for a fixed four-set `R`, then

```text
C=3-2r+binom(r,2)
```

has layer values `[3,1,0,0,1]` and `4p E[C]=36=p+13`.  This is only a local
quadratic witnessing sharpness of the parameterized local lemma.  It is
**not** a residual graph and does not establish that the `p=23` residual
layer survives.  The exceptional equality theorem classifies all such
height-three cells as `F4` or `F5`, globalizes them across paired cubes, and
leaves only the `F5` coefficient offset.  Twelve hard roots then force common
quartic and octic moment identities, while an exact scan of all
`binom(23,5)=33,649` five-sets finds no simultaneous opposite-row solution.
Thus `p=23,t=9,k=110` is closed separately.

Proposition 15.274's slope obstruction applies only to its dual-bad,
two-level subcase.  It does not apply to this multilevel isolated-chart
equality branch.

Therefore

```text
boxed: k=5p-5 is impossible for every prime p=3 mod 4, p>=31.
exceptional companion: p=23,t=9,k=110 is also impossible.
```

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python src/e1_gmin_m4_prop15769.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 tests/test_prop15769.py
```
