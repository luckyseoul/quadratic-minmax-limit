# First all-finite survivor closed from `p=19`

Date: 2026-08-28. This is Proposition 15.677. Together with Proposition
15.675, it excludes the first even all-finite boundary size strictly above
`3(p-1)/4` for every prime `p>=19`. It does not close the exceptional
smaller endpoint `p=17`, the next even boundary size, general residual (ii),
R1, global QVAR, Type I, or the limit.

## 1. The two residue classes left by 15.675

Put

```text
P=p+1=2m,    q=(p-1)/2,
```

and let `s` be the first even integer strictly above `3(p-1)/4`. Proposition
15.675 already closes `p=3,5 mod 8`. In the remaining classes,

```text
p=1 mod 8: s=(3p+5)/4,
p=7 mod 8: s=(3p+3)/4.
```

For either quadratic direction type, the exact same-type mean identity is

```text
a_d=2u+P k_d,    0<=u<m,    sum_d k_d=m-u.             (1)
```

Phase one has only `u_1=m-1`: exactly `m-1` directions are floor-attaining
`b=2` xnor baselines and the remaining direction has `b=s`. Its exact pair
deficit is

```text
D_1=(m-1)(s-2).                                        (2)
```

At phase-zero residue `u_0=2`, the exact minimum is

```text
D_0=((m+3)/2)s-2   if m is odd,
D_0=((m+2)/2)s     if m is even.                       (3)
```

After (2)--(3), the unused pair-deficit budget is `(p-1)/4` in the
`p=1 mod 8` class and `(p-7)/4` in the `p=7 mod 8` class. The exact quotient
increments then leave

- `u_0=2` in both classes;
- `u_0=3` additionally when `p=1 mod 8`.

Residues zero and one are excluded by the exact low-floor ledger; residues
four and higher already exceed the unused pair budget. The inequalities are
uniform from `p=23`. This lower endpoint suffices because the first prime
`p=7 mod 8` in scope is 23 and the first later prime `p=1 mod 8` is 41.

## 2. Independent pre-lift coefficient normal form

Let `j-1` be the phase-zero parallel-count offset and let `l` be the common
parallel count in the `m-1` xnor baseline directions. Exact coefficient and
edge counts give

```text
N_0=mj-u_0,    N_1=ml+1,
E=N_0+N_1,     I=4p+1-E,
q | I+l-4,     hence q | u_0-j.                        (4)
```

The support bound gives `j+l<=7`. Since `q>=11`, (4) forces `j=u_0`.
The exact inter-fibre `l1` lower bound excludes every `u_0=3` coefficient
row and every odd `l` row. Before using quadratic-lift mass, only two
arithmetic regimes remain:

```text
l=2: I=2p,   E=2p+1,
l=4: I=p-1,  E=3p+2.                                  (5)
```

This normal form is useful diagnostically, but the next argument eliminates
the branches earlier and uniformly.

## 3. The zero-quotient lift contradiction

Every retained phase-zero residue satisfies `u_0>0`. Equation (1) has
quotient sum `m-u_0<m`, so at least one phase-zero direction has `k_d=0`.
Its scaled mean is `a_d=2u_0`, equal to four or six.

The boundary is all-finite and even, hence `b_d` is even. The phase-zero
floor at every nonzero even `b_d` is at least `p+1`, while
`2u_0<=6<p+1`. Therefore this zero-quotient direction has `b_d=0`.
Its pointwise parity is zero, so its nonnegative slack factors as

```text
A_d=2B_d,
```

where `B_d` is a nonzero nonnegative integer-valued quadratic on the middle
slice. It is nonzero because `a_d>0`. Consequently

```text
2u_0 = 2p E[A_d] = 4p E[B_d].                          (6)
```

Proposition 15.642, using the exact degree-two polynomial-distance lemma on
the slice, gives

```text
4p E[B_d] >= nonbaseline_scaled_cost_floor(p) >= 8
```

for every `p>=23`. This contradicts (6) for both `u_0=2` and `u_0=3`.
Thus neither pre-lift regime in (5), nor any other retained residue row, can
come from a residual graph.

## 4. Combined theorem and the `p=17` endpoint

Proposition 15.675 handles primes `p>=19` congruent to 3 or 5 modulo 8.
The argument above handles primes `p>=23` congruent to 1 or 7 modulo 8.
There is no omitted prime between these ranges, so the first all-finite
survivor is excluded for every prime `p>=19`.

No claim is made for `p=17`. Its exact residue ledger has an additional
phase-zero `u_0=0` row, so the zero-quotient argument does not apply. Keeping
this endpoint explicit prevents the uniform theorem from silently assuming
a false residue reduction.

## 5. Near-perfect fibre profiles

For completeness, an independent `l1` classification applies to every one
of the `m-1` xnor baseline directions in (5). If `n_t` counts infinity
neighbours in fibre `t`, the only profiles not already over transverse-edge
capacity are

- `l=2`: all `n_t=2`, or one 1, one 3, and all remaining entries 2;
- `l=4`: `p-1` entries 1 and one 0, or one 2, `p-3` entries 1, and two 0s.

These are pre-lift profiles only; Section 3 proves that no in-scope graph
reaches them.

## 6. Literature, OEIS, and reproduction

The only imported ingredient in the decisive step is Lemma 2 of
Amireddy--Behera--Srinivasan--Sudan,
[A Near-Optimal Polynomial Distance Lemma over Boolean Slices](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2025.11),
already specialized and proved into the exact finite-parameter ledger of
Proposition 15.642. Targeted searches found no source combining that slice
support bound with Paley same-type mean residues and the zero-quotient
direction. Searches for the sample lift-cost values returned unrelated OEIS
entries; this proposition is a structural implication, not a new sequence
claim.

Reproduce the machine-readable arithmetic record with

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15677.py
PYTHONPATH=src pytest -q tests/test_prop15642.py tests/test_prop15669.py \
  tests/test_prop15675.py tests/test_prop15676.py tests/test_prop15677.py
```

The generated record is `evidence/e1_gmin_m4_prop15677.json`.
