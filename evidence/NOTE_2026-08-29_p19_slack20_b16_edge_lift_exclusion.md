# p=19 slack-20 b=16 edge-lift exclusion

**Date:** 2026-08-29
**Proposition:** 15.696
**Status:** proved computationally; one profile excluded, endpoint open

The remaining mixed slack-20 row has phase profiles

```text
phase zero: {0:5, 2:1, 16:4}
phase one:  {2:9, 16:1}.
```

The phase-one floors saturate the exact type budget:

```text
9*18 + 38 = 200.
```

The `b=16` floor polynomial has values

```text
q(7)=0, q(8)=1, q(9)=4/3, q(10)=1,
```

and contacts `t=7,8,10`. Equality forces the original integer slack to be
pointwise `0,1,1` on those three layers. The quadratic evaluation map on
these layers in `J(19,10)` has rank 169. A checked-in 169-row minor has rank
169 modulo two. Its rational kernel has the two-dimensional form

```text
(alpha dot y_C)(2-|y_C|),  sum(alpha)=0,
```

where `C` is the set of three empty fibres. On `t=9`, the three possible
values are nonnegative even integers with sum four, hence there are two
orbits, `{0,2,2}` and `{0,0,4}`. Their pointwise normal forms are

```text
A_022 = 1-z+u+v-2uv,
A_400 = 1+3z-u-v-2zu-2zv+2uv.
```

Comparing pair coefficients gives, for every fibre pair,

```text
w_st = T_st + g - h_s - h_t,
P_d  = 4 + 9g - I.
```

For a rigid `b=2` direction, `T` is `+1` on its odd pair and zero elsewhere.
For the `b=16` form, `T` is `-1` on every pair containing `z` and on `{u,v}`.
Summed signed-cross capacity, the phase-zero means, exact parallel-edge
accounting, even infinity degree, and a selected infinity-zero edge leave
exactly

```text
I = 2,8,10,12,18,20,28,30,38,48.
```

The exact affine edge-lift model selects all 77 edges from the 65,341 edge
variables, imposes the odd-degree boundary, edge-product sign, all twenty
direction profiles and means, 1,720 rigid coefficient identities, and the
single phase-zero elevation. Each infinity-degree shard for each of the two
normal-form orbits returned exact `INFEASIBLE` under OR-Tools CP-SAT
9.15.6755 on soulkiller or jellyfin. There are twenty logical
shape/degree shards. The hard `{0,2,2},I=28` shard is losslessly partitioned
by the role of the unique elevated phase-zero direction (`b=0,2,16`), and all
three subshards are `INFEASIBLE`. Thus the corrected archive has 22 raw JSON
files, all individually hashed and audited by Proposition 15.696.

## Finite-field sign correction

The original raw shard archive used integer subtraction on the encoded
values `a+b*p`. That is not subtraction in the chosen `F_{p^2}` basis when
the low component borrows into the high component. The affected model now
forms

```text
((a0-b0) mod p) + ((a1-b1) mod p)*p
```

before applying the quadratic character. A regression test compares every
one of the 65,341 model edge signs with
`minmax_quadratic.paley_conference_prime_power(19)`. The original shards are
retained only as an invalid historical record; Proposition 15.696 reads the
`*_correct.json` archive generated after this repair. No original shard,
timeout, or `UNKNOWN` result is used as evidence.

It is enough to solve `c_H=+1`. Multiplication by a nonsquare fixes infinity
and finite zero, flips every finite-finite conference sign and the direction
type, and therefore flips both `c_H` and `eps` while preserving their phase
product. Since infinity has even degree, the number of finite-finite selected
edges is odd. Thus any `c_H=-1` witness transfers to a `c_H=+1` witness.

The square torus fixes infinity and zero and acts regularly on the ten
directions of either quadratic type, justifying normalization of the unique
phase-one `b=16` direction to the first phase-one slot.

The p=19 remainder is now four profiles:

```text
{20:1, 24:1, 28:1, 32:1}.
```

## Reproduction

- `src/e1_gmin_m4_prop15696.py`
- `scripts/p19_slack20_b16_lift_cpsat.py`
- `tests/test_prop15696.py`
- `evidence/e1_gmin_m4_prop15696.json`
- `evidence/p19_slack20_b16_lift_shards/`

This proposition uses only completed `INFEASIBLE` shards. No timeout or
`UNKNOWN` result is treated as evidence.
