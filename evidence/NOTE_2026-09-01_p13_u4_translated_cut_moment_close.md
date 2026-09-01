# p=13,t=4,u=4: translated-cut moment close

**Proposition:** 15.749  
**Status:** proved branch theorem with an exact aggregate certificate

## Input left by Proposition 15.748

The omitted-pair `P=3` branch is already empty.  In the remaining
all-equal-triple `P=5` branch, exactly two opposite directions have `Q=3`
and are literal cells.  The other five have `Q=4`.  Proposition 15.748
records 336 exact `z=2` moment candidates for each hard sign.

The common signed total is `hT=9`.  For an opposite `Q=4` direction, let

```text
q_a=(-h) sum_(|L(delta)|=a) chi(delta)m_delta,   1<=a<=6.
```

Then

```text
sum_a q_a=-hT-Q=-13,       sum_a |q_a|<=61-Q=57.
```

The phase-zero cell satisfies

```text
(-h)S_H=Q+sum q_a-2 cut_W=-9-2 cut_W=3+2A,
```

so `A=-6-cut_W>=0`.  Summing over the thirteen additive translates of
any seven-set gives

```text
c dot q <= -78
```

for each of Proposition 15.740's already-certified 74 translated-cut
vectors `c`.

## Exact coordinate bounds

In the canonical cut-vector order, write `c_i` for row `i` and `1` for the
six-coordinate all-one vector.  Direct coefficient identities give

```text
 e_1 = (19/9)1 - c_0/18 - c_6/6 - c_34/18,
-e_1 = (29/15)1 - c_63/15 - c_69/30 - c_71/6 - c_73/30.
```

All cut coefficients on the right are nonpositive, so `sum q=-13` and
`c_i dot q<=-78` yield

```text
-52/9 <= q_1 <= 26/15.
```

Multiplication by `F_13^*/{+-1}` acts transitively on the six cyclic
distance bins and preserves the complete 74-vector catalog.  Hence the same
bounds hold in every coordinate.  Integrality sharpens them to

```text
-5 <= q_a <= 1.
```

## Moment intersection

Exact list recovery in `[-5,1]^6`, using only the row sum and all 74 cuts,
gives 522 rows and 492 distinct triples

```text
(N2,N4,N6)=(sum a^2 q_a, sum a^4 q_a, sum a^6 q_a) mod 13.
```

The row and moment-list hashes are

```text
2e4cf7f733ffd6d85a68a6b37ebd380d93d962e6119e9306673ddd3a1df8cb35
16e529a3ea7263b66f7af61cf6eaa59441747622c62d0dbc9a8267837c76f378
```

For each of Proposition 15.748's 336 records, reconstruct the homogeneous
forms

```text
M2=cR2,       M4=R2 Q2,       M6=R2 Q4.
```

The stored hard values uniquely recover `Q4` in the degree-four evaluation
code.  The records use hard normalization `hM_(2r)`, whereas an opposite
row has moment `(-h)M_(2r)`, so the required local triple is the negative of
the reconstructed form values.

For either hard sign the five nonroot evaluations range over the same
48-element alphabet.  Its intersection with the 492 admissible `Q=4`
triples is exactly

```text
(1,0,3), (2,0,1), (3,0,3), (4,0,10),
(5,0,1), (6,0,1), (7,0,12), (8,0,12),
(9,0,3), (10,0,10), (11,0,12), (12,0,10).
```

Thus every admissible intersection triple has `N4=0`.  Independently, the
per-survivor compatibility histogram is

```text
0 compatible Q4 directions: 252 survivors
1 compatible Q4 direction:    42 survivors
2 compatible Q4 directions:   42 survivors
```

No record supplies all five required directions.

## Contradiction

Each of the two literal directions is already a root of the common quartic
`M4`.  The moment intersection forces every one of the other five `Q=4`
directions to be another root.  Seven distinct projective roots force a
homogeneous binary quartic to vanish identically.  Proposition 15.748's
hard alphabet has no zero fourth moment, a contradiction.

Therefore the `P=5` branch is empty.  Proposition 15.747 already closes
`P=3`, so

```text
p=13,t=4,u=4 is empty.
```

The exact `p=13,k=60` remainder is now `u=6`.  Residual (ii) globally,
multi-level Type I, and the quadratic-minmax limit remain open.

## Artifacts

- `src/e1_gmin_m4_prop15749.py`
- `tests/test_prop15749.py`
- `evidence/e1_gmin_m4_prop15749.json`
