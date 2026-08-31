# Binary quadratic moment close of the first three `p=11` layers

**Date:** 2026-08-31
**Proposition:** 15.737
**Result status:** proved theorem, using Proposition 15.736's exhaustive finite
Boolean-quadratic catalog

## Result and exact scope

Residual (ii) at `p=11` is empty for `k=44,46,48`, equivalently
`|H|=45,47,49`. Proposition 15.736 supplies the sharp Boolean catalog. At
layer index `t=0,1,2`, its coefficient offsets exclude the hard-`b=2`
branch. The other branch has at least `5-t` hard low `b=10` baselines and at
least `4-t` opposite minimum directions with all-equal-triple targets. A
single binary quadratic moment makes those two requirements incompatible.

This does not treat critical `p=5,7`, even `k>=50` at `p=11`, or residual
(ii) as a whole.

## 1. Direct p=11 floor replay

For `t=0,1,2`, at most `2(45+2t)` of the 122 projective vertices are
nonisolated. The isolated gaps are `32,28,24`. Proposition 15.721's signed
PSL transport therefore gives `I=0`, an all-finite boundary, and even
directional odd-fibre counts.

The phase-one type budget is `12(6+t)`. Write

```text
a_d=2u+12k_d,    sum_d k_d=6+t-u.
```

For `0<=u<=t`, a `k=1` direction is forced at mean `12+2u`. The exact
phase-one even floors are

```text
b:      0  2  4  6  8 10
floor: 22 10 22 18 22 10.
```

Thus every compatible cell is a nonzero lift of excess `2,4,6`, below
Proposition 15.688's sharp floor eight. For `t<u<5`, `k=0` is below every
floor but the quotient sum is less than six. Hence only `u=5` survives,
with at least `5-t` exact mean-ten directions. Proposition 15.652's positive
quadrature makes their baselines pointwise:

```text
b=2:   A=(1-x_i-x_j)^2,  eps S_H=4+z_i z_j, offset 4;
b=10:  A=1-x_j,          eps S_H=4-z_j,     offset 3.
```

Equal means give equal parallel counts, so the offsets cannot mix modulo
`q=5`.

With `I=0`, the two branches force respectively `P=4` and `P=3`. Their
opposite parallel sums are `20+t` and `26+t`. A minimum opposite direction
has `Q=3` or `Q=4` and scaled mean eight. A nonzero `b` costs at least 12;
at `b=0`, Proposition 15.688 gives `A=2B` with `B` Boolean at equality.
Proposition 15.736's exhaustive catalog has offsets 2 and 4. Neither matches
`Q=3 mod 5`, excluding the hard-`b=2` branch. At `Q=4`, only the all-equal
triple remains. Since the surplus is `t+2<6`, at least `4-t` opposite
directions attain that minimum.

## 2. The cross-direction moment

For an `F_11`-linear fibre functional `L`, define

```text
M_H(L) = sum_({u,v} in H) chi(u-v) (L(u)-L(v))^2.
```

This is a homogeneous binary quadratic in the two coefficients of `L`. If
`K^L_st` denotes the signed selected-edge sum between fibres `s,t`, then

```text
M_H(L) = sum_(s<t) K^L_st (s-t)^2,
```

because the edges parallel to `ker L` contribute zero. The sign convention
is

```text
eps_L S_H = P_L + eps_L sum_(s<t) K^L_st z_s z_t.
```

## 3. At least three hard stars force `M_H=0`

In the remaining branch, each of at least `5-t>=3` hard low directions has
`P=3` and target `eps_L S_H=4-z_j`. Subtract the target and write the
slice-zero quadratic as the multilinearization of
`(sum z_s-1)(c+sum a_s z_s)`. Constant and linear comparison gives `c=0`,
`a_j=-1`, and every other `a_s=0`. Thus

```text
eps_L K^L_st = -1  when exactly one of s,t is j,
                  0  otherwise.
```

It follows that

```text
M_H(L) = -eps_L sum_(t!=j)(j-t)^2 = 0 in F_11,
```

since `sum_(u in F_11^*)u^2=385=0 mod 11`. These hard low directions are
distinct points of the projective line of linear forms. A nonzero homogeneous
binary quadratic has at most two projective zeros. Therefore `M_H` is
identically zero.

## 4. No all-equal triangle has zero moment

An opposite minimum direction has `P=4` and target

```text
eps_L S_H = 4+z_i z_j+z_i z_k+z_j z_k.
```

The same coefficient comparison has `c=a_s=0` and gives
`eps_L K^L_st=1` on the triangle `{i,j,k}` and zero off it. Thus `M_H=0`
would require

```text
(i-j)^2+(i-k)^2+(j-k)^2=0.
```

For three distinct fibres, translate and scale to `(i,j,k)=(0,1,r)`. The
left side becomes `2(r^2-r+1)`. Its discriminant is `-3=8 mod 11`, while the
nonzero squares modulo 11 are `{1,3,4,5,9}`. Hence the moment never vanishes.
The deterministic audit also checks all `C(11,3)=165` triples directly.

One all-equal target already contradicts `M_H=0`; the edge ledger forces at
least `4-t>=2`. This excludes the last branch and closes all three layers.

The argument stops exactly at `t=3`: the hard branch then guarantees only
two star zeros, and low-residue excess eight reaches the sharp equality
floor, introducing additional equality cells not classified here.

## Reproduction

- theorem and exact arithmetic: `src/e1_gmin_m4_prop15737.py`;
- regression tests: `tests/test_prop15737.py`;
- deterministic certificate: `evidence/e1_gmin_m4_prop15737.json`.
