# Complete `p=5` negative two-point exclusion

For `D={infinity,v}`, `c_H=-1`, and `p=5`, the parity baseline in every
direction is `1-x_j`, with scaled mean `a=4`.  Write `a_d=4+r_d`.  Each
quadratic type has three directions and exact excess budget

```text
sum_d r_d = 6.
```

Every `r_d` is a nonnegative even integer.  The exact same-type directional
mean identity gives

```text
r_d-r_e = 6(P_d-P_e).
```

Therefore all three excesses in one type are congruent modulo six.  There
are exactly two possibilities:

```text
unique:       (r_d) = (6,0,0),   P_exception=P_baseline+1;
distributed:  (r_d) = (2,2,2),   all three P_d equal.
```

If `e_+` and `e_-` indicate the unique profile in the positive and negative
types, and `x,y` are the corresponding baseline/common parallel counts,
then

```text
E = 3(x+y)+e_++e_-,
I = 21-E.
```

We impose the following exact necessary arithmetic:

- `I` is positive and odd because infinity is in the boundary;
- `3y+e_-` is odd because the finite negative-edge product is `-1`;
- `I-1<=2E` because an `E`-edge finite graph must realize the finite
  boundary `star symmetric-difference {v}`;
- a type with a zero-lift baseline has even baseline parallel count.

Exactly 24 arithmetic profiles survive.  Square multiplications in
`F_25^*` and Frobenius have two orbits on opposite-type exceptional pairs
(sizes six and three) and are transitive on each single direction type.
Consequently the 24 profiles have exactly 33 placement orbits:

```text
9 unique/unique profiles x 2 pair orbits = 18
5 unique/distributed profiles             =  5
5 distributed/unique profiles             =  5
5 distributed/distributed profiles         =  5
                                             --
                                             33
```

For each orbit, the finite CP-SAT model selects edges among the 300 finite
affine pairs and the 25 infinity-star leaves.  It enforces exact finite and
infinity edge counts, all six direction counts, the 25 boundary XORs,
negative edge-product parity, and all 60 affine score identities.  The lift
masses are three in a unique exceptional direction, zero in its two
baselines, and one in every distributed direction.

A fresh 33-way reproducibility run certified every representative
`INFEASIBLE`, with no feasible or unknown result.  The raw certificate hash
is

```text
2352ca1040a7989d5850730fafcad311adf1275dd06383578fb9e187698be69e
```

and the permanent archive is

```text
/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-26-negative-p5/
  p5_negative_two_point_certificate_2026-08-26.tar.gz
SHA256 c19ed8a0d50ffad3f7386d3d6100ce25213c3ba4d1e7cc49acea0275d3796a41
```

Thus `p=5` is excluded.  Combined with Propositions 15.647--15.649, the
negative-product infinity-plus-point branch is closed for every odd prime
`p>=5`.  Positive-product finite cases, other boundary profiles, residual
(ii) as a whole, and R1 remain open.
