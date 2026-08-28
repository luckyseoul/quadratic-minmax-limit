# Quantized close of the first all-finite survivor

Date: 2026-08-28. This is Proposition 15.675. For every odd prime `p>=19`
with `p=3 or 5 mod 8`, it excludes the first even all-finite boundary size
strictly above Proposition 15.669's uniform range `3(p-1)/4`. The same
relaxation has explicit negative gaps for `p=1 or 7 mod 8`, so those classes
remain open.

## 1. The added information

Let

```text
P=p+1=2m
```

and let `s` be the first even integer strictly above `3(p-1)/4`. In one
quadratic direction type, the exact directional mean formula gives a common
even residue modulo `P`:

```text
a_d=2u+P k_d,    0<=u<m,    sum_d k_d=m-u.             (1)
```

Proposition 15.669 used the exact floors and pair deficit but not (1). The
new argument imposes (1) before minimizing the type deficit.

Throughout this first-survivor range, the relevant floors are

| `b` | phase zero | phase one |
|---:|---:|---:|
| `0` | `0` | `2P-2` |
| `2` | `P` | `P-2` |
| `4` | `2P-8` | `2P-2` |
| `6<=b<=s` | `2P-2` | `2P-2` |

The nonzero-lift theorem forbids a value exactly two above any floor.

## 2. Phase one is rigid

For `u<=m-2`, every phase-one direction has `k_d>=1`, but (1) requires
`sum k_d=m-u<m`. At `u=0`, the apparent first values over both floor classes
are forbidden two-unit lifts, making the contradiction stronger. Thus only
`u=m-1` survives.

There `sum k_d=1`. The minimum deficit has `m-1` floor-`P-2` directions at
`b=2` and one floor-`2P-2` direction at the largest possible count `b=s`:

```text
D_1=(m-1)(s-2).                                        (2)
```

## 3. Phase zero has residue four

For every `u>=2`, the deficit-optimal quotient weights at `b=0,2,s` are
respectively `0,1,2`. If `t` is the parity bit of `m+u`, the minimum is

```text
D_0(u)=s(m+u+t)/2 - 2t.
```

Increasing `u` alternately raises this expression by `s-2` and by `2`, so
the minimum is `u=2`. Consequently

```text
m even: D_0=s(m+2)/2,
m odd:  D_0=s(m+3)/2-2.                                (3)
```

At `u=0`, quotient weights `0,1,2,3` occur at `b=0,2,4,s`. A weight-three
middle direction saves `s` deficit units, while each remaining quotient
unit saves two. Writing `m=3h+t`, `0<=t<=2`, gives

```text
D_0(0)=(m-h)s-2t,
```

which is strictly larger than (3) for every `m>=10`. The six classes of
`m mod 6` give margins `(a-1)s`, `as-4`, `as-2`, `(a-1)s`,
`(a-1)s+2`, or `as-2`, with the applicable `a` positive. The case `u=1`
is impossible because the residue two is itself a forbidden lift over the
zero floor.

## 4. Exact modulo-eight gap

The first survivor is

| `p mod 8` | `s` |
|---:|---:|
| 1 | `(3p+5)/4` |
| 3 | `(3p-1)/4` |
| 5 | `(3p+1)/4` |
| 7 | `(3p+3)/4` |

Substituting these four values into (2)--(3) and subtracting the exact pair
budget `s(s-1)` gives

| `p mod 8` | quantized deficit gap |
|---:|---:|
| 1 | `-(p-1)/4` |
| 3 | `(p+1)/2` |
| 5 | `(p-1)/2` |
| 7 | `-(p-7)/4` |

The positive middle two rows prove the exclusion for `p=3 or 5 mod 8`.
The negative outer rows are an exact route boundary, not evidence of actual
graphs.

## 5. Verification, literature, and OEIS

`src/e1_gmin_m4_prop15675.py` contains both the symbolic formulas and an
independent dynamic program over every common residue, every allowed `b`,
and every relaxed lift of size at least four. The DP agrees on all sample
primes through 43. NUKA independently returned the same gaps

```text
p=19,23,29,31,37,41,43:  10,-4,14,-6,18,-10,22.
```

The adjacent odd-secant literature includes Ball--Csajbók,
[On sets of points with few odd secants](https://arxiv.org/abs/1711.10876),
but the searched sources do not impose Paley quadratic-type mean residues on
the odd-fibre deficit optimization. Targeted OEIS searches for the positive
gap values and the signed sample list returned unrelated arithmetic
progressions, Aurifeuillean numbers, and combinatorial arrays. The four gap
formulas are elementary residue-class evaluations, not a sequence claim.

## 6. Reproduction and scope

```bash
python src/e1_gmin_m4_prop15675.py
python -m pytest -q tests/test_prop15675.py
```

The generated record is `evidence/e1_gmin_m4_prop15675.json`. The first
survivor for `p=1,7 mod 8`, the next even size, the infinity-present
remainder, general residual (ii), R1, global QVAR, Type I, and the limit all
remain open.
