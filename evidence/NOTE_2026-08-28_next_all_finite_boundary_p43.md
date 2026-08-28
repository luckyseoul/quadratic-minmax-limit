# The next all-finite boundary is impossible for every prime p>=43

Date: 2026-08-28. This is Proposition 15.679. Propositions 15.675,
15.677, and 15.678 exclude the first even all-finite boundary size above
`3(p-1)/4` for every prime `p>=17`. Here the same exact mean residue and
quadratic-lift machinery excludes the **next** even size for every prime
`p>=43`.

This is a uniform second-boundary theorem, not residual (ii). At this size
the seven smaller endpoints `p=17,19,23,29,31,37,41` remain open, as do
later all-finite sizes and the infinity-present remainder.

## 1. Parameters and exact phase-one type

Write

```text
P=p+1,     m=P/2.
```

The first survivor from Proposition 15.675 and the next even size `s` are

```text
p mod 8       first survivor        s=first+2
   1            (3p+5)/4             (3p+13)/4
   3            (3p-1)/4             (3p+7)/4
   5            (3p+1)/4             (3p+9)/4
   7            (3p+3)/4             (3p+11)/4.
```

For every in-scope residue class,

```text
s <= p-5.
```

Indeed, the four margins `p-5-s` are respectively
`(p-33)/4,(p-27)/4,(p-29)/4,(p-31)/4`.

For one quadratic direction type, all exact scaled means have the common
residue form

```text
a_d=2u+P k_d,      0<=u<m,      sum_d k_d=m-u.       (1)
```

In phase one, the smallest even-fibre floor is `P-2`, attained only at
`b=2`; all other even counts through `s` have floor `2P-2`. If `u=0`, the
apparent quotient-one `b=2` value has forbidden lift two, so every quotient
is at least two and (1) is impossible. If `1<=u<=m-2`, every quotient is at
least one, whereas their sum `m-u` is less than the `m` directions. Thus
only `u=m-1` remains. Exactly one quotient is one and the rest are zero;
maximizing its fibre count gives

```text
phase one:  b=2 on m-1 directions, b=s on one direction,
D_1=(m-1)(s-2).                                      (2)
```

## 2. Exact phase-zero residue reduction

For phase zero and `2<=u<=m-5`, the largest permissible fibre counts at
quotients zero, one, and two are respectively

```text
k=0: b=0,       k=1: b=2,       k=2: b=s.            (3)
```

At quotient two, `b=s` dominates `b=4`; quotients greater than two spend
more of the fixed sum without reducing the deficit further. Put

```text
t=m-u,  x=floor(t/2),  y=t mod 2,  z=m-x-y.
```

The exact phase-zero minimum is therefore

```text
x directions at b=s, y at b=2, z at b=0,
D_0(u)=z s+y(s-2).                                   (4)
```

As `u` increases, successive increments in (4) are `2` and `s-2`, so
`D_0(u)` is strictly increasing.

There are three exceptional ranges to check.

### Residue zero

At `u=0`, the largest fibre counts for quotients zero through three are
`0,2,4,s`; the quotient-two value at `b>=6` has forbidden lift two. Hence

```text
D_0(0)=(m-floor(m/3))s-2(m mod 3).                   (5)
```

For `p>=47`, (5), (2), and the pair budget give the lower gap

```text
D_0(0)+D_1-s(s-1)
 >= s(p-29)/12-p-3
 >= (3p^2-128p-347)/48 > 0.                         (6)
```

The derivative of the final numerator is positive in this range, and its
value at 47 is 264. At `p=43`, the exact values are

```text
D_0(0)=508, D_1=672, s(s-1)=1122, gap=58.
```

Thus residue zero is always over budget. Residue one is infeasible because
its quotient-zero lift is the forbidden value two, forcing all `m`
quotients positive although their sum is `m-1`.

### Interior residues at least eight

Since (4) is increasing, it is enough to evaluate `u=8`. Direct substitution
in (4) and (2) gives the following positive pair gaps:

```text
p mod 8     D_0(8)+D_1-s(s-1)
   1          2s-p-1 = (p+11)/2
   3          3s-p+1 = (5p+25)/4
   5          3s-p-1 = (5p+23)/4
   7          2s-p+1 = (p+13)/2.                    (7)
```

Therefore every interior residue `u>=8` is over budget.

### The final four residues

If `u>=m-4`, the quotient sum is at most four. Every nonzero fibre count
uses a positive quotient, so at least `m-4` directions have `b=0`. Hence

```text
D_0 >= (m-4)s.
```

After adding (2), the pair gap is at least

```text
s(p-s-3)-p+1
 >= (3p^2-84p-159)/16 > 0                           (8)
```

for `p>=43`. Equations (5)--(8) leave only

```text
2 <= u <= 7.                                         (9)
```

This is a constant-size residue window independent of `p`.

## 3. The zero-quotient quadratic-lift contradiction

For every residue in (9), equation (1) has quotient sum `m-u<m`, so some
phase-zero direction has `k_d=0` and scaled mean `2u<=14`. Every nonzero
even fibre count has floor at least `P`, while `14<P`. Thus this direction
has `b=0`.

Its parity target is zero. Its nonnegative pointwise slack factors as

```text
A_d=2B_d,
```

where `B_d` is a nonzero nonnegative integer-valued polynomial of degree at
most two on `J(p,m)`. It is nonzero because `u>0`, and

```text
2u = 2p E[A_d] = 4p E[B_d].                          (10)
```

Proposition 15.642, using Amireddy--Behera--Srinivasan--Sudan's exact
degree-two polynomial-distance lemma on the slice, gives

```text
4p E[B_d] >= (p^2-1)/(4(p-2)).                       (11)
```

For `p>=59`, the right side is strictly greater than 14 because

```text
p^2-56p+111>0,
```

and that polynomial is positive and increasing from `p=59`. Equations
(9)--(11) contradict each other.

The only primes from 43 through 58 are 43, 47, and 53. Their exact rows from
(2)--(4) are

```text
p    pair-surviving u values       max 2u     Prop. 15.642 floor
43          2,3,4                     8                12
47          2,3,4,5,6                12                14
53          2,3,4,5                  10                14.
```

Each row again contradicts (10). Therefore

```text
for every prime p>=43, the second even all-finite boundary size above
3(p-1)/4 is impossible.                              (12)
```

## 4. Independent arithmetic replay

The source contains a second dynamic program over every common residue,
every even `b<=s`, every quotient, the exact symbolic floor, and the
forbidden two-unit lift. It agrees with (2), (4), and (5) at
`p=43,47,53,59,73`. The theorem record also checks representative primes in
all four classes through 101.

This replay verifies the finite arithmetic. The uniform content of the
proof is the symbolic inequalities (6)--(8) and (11), not extrapolation
from those samples.

## 5. Scope

Proposition 15.679 closes exactly one additional all-finite size for every
prime `p>=43`. It does **not** close:

- this same size at `p=17,19,23,29,31,37,41`;
- any later all-finite size;
- the strict-deficit infinity-plus-`p` shell;
- residual (ii), R1, global QVAR, Type I, E(1), or the limit.

## 6. Literature and OEIS check

The load-bearing external input is the exact degree-two slice-distance
bound in Amireddy--Behera--Srinivasan--Sudan [47]. Ball--Csajbók [42]
studies few odd secants for projective sets of size `q+2`, not the present
roughly `3p/4` affine boundary with two Paley direction types. Targeted
searches of the direction and slice-polynomial literature found no theorem
combining (1), the pair budget, and the zero-quotient lift argument.

A fresh 2026 search also found Amireddy--Behera--Srinivasan--Sudan--
Willumsgaard [50], which develops low-degree **testing** over Boolean slices.
It does not replace or strengthen the exact finite support floor used in
(11), and it does not address the affine directional budget.

Exact OEIS API searches for the larger ledger blocks
`11130,7176,3922` and `21756,14016,7696` both returned `null`. This is a
duplicate/context check only; no sequence submission or priority claim is
made.

## 7. Reproduction

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15679.py
PYTHONPATH=src pytest -q tests/test_prop15679.py
```

The generated exact record is
`evidence/e1_gmin_m4_prop15679.json`.
