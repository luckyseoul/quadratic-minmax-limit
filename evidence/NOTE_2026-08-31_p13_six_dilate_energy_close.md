# Exact six-dilate energy close of the generic p=13 fourth shell

**Date:** 2026-08-31
**Proposition:** 15.742
**Result status:** exhaustive finite certificate
**Scope:** the last generic `p=13,t=3` partition `1^4 2^3`

## The implication closed

Proposition 15.741 reduced the last generic `p=13,k=58` branch to ten
nonexact six-bin difference rows of one common 59-edge graph: three elevated
hard rows and seven opposite rows.  If

`q_L(a)=epsilon_L sum_(L delta=+-a) chi(delta)m_delta`, `1<=a<=6`,

then its exact common-graph Parseval identity is

`sum_(L nonexact) ||q_L||_2^2 = 707+26C`,

where `C=sum_delta binom(m_delta,2)>=0`.  The unresolved implication was to
show that the required integral difference system has no solution.  The
binary midpoint lift is not needed.

## Six exact row inequalities

Take the interval seven-set `{0,1,...,6}` and its six nonzero multiplicative
dilates.  In the natural cyclic-distance order their translated-cut vectors
are

```
( 2,  4,  6,  8, 10, 12)
(12,  2, 10,  4,  8,  6)
( 8, 10,  2,  6, 12,  4)
( 6, 12,  8,  2,  4, 10)
(10,  6,  4, 12,  2,  8)
( 4,  8, 12, 10,  6,  2).
```

All six belong to Proposition 15.740's exact 74-vector catalog.  Proposition
15.741 gives every elevated row the necessary conditions

```
q in Z^6,  sum(q)=11,  ||q||_1<=53,  ||q||_2^2<=86,
sum_(a=1)^6 a^2 q(a)=0 (mod 13),  c*q<=91 for the six rows c above,
```

and every opposite row the necessary conditions

```
q in Z^6,  sum(q)=-20,  ||q||_1<=56,  ||q||_2^2<=106,
sum_(a=1)^6 a^2 q(a)=0 (mod 13),  c*q<=-130 for the six rows c above.
```

The sums are the signed nonparallel masses.  The `l1` bounds are simply the
numbers of nonparallel edges, `59-6=53` and `59-3=56`.  The congruence is the
already-proved global identity `M_2=0`.  The cut bounds are the sums over the
thirteen translates of one balanced seven-set.  The broad energy bounds
`86,106` are the exact six-dilate spectral bounds already proved in 15.741;
using them to make the following finite search bounded is not circular.

## Exhaustive row maxima

For each type, enumerate five integer coordinates and recover the sixth from
the fixed sum.  Energy alone gives `|q_a|<=9` in the elevated case and
`|q_a|<=10` in the opposite case, so the implementation's wider coordinate
intervals omit nothing.  Filtering by `l1`, energy, the degree-two
congruence, and the six displayed cut inequalities gives

| row type | before six cuts | after six cuts | exact maximum energy |
|---|---:|---:|---:|
| elevated | 5,844 | 30 | 31 |
| opposite | 1,704 | 24 | 82 |

The surviving-row SHA-256 digests are respectively

```
7bff1ebb77ac362b5089b46588f603be812ea4a96fdeaa2f2b52881803b486b5
5226c7ee0c44d3cf7e460db2a309a08368667686dfcb2cd45ec24ce932081c1a
```

There are six maximizers in each case.  The elevated representatives are
the six multiplicative-distance images of `(0,3,1,4,1,2)` and have energy
31.  The opposite representatives are the six multiplicative-distance
images of `(-6,-1,-4,-2,-4,-3)` and have energy 82.  Equivalently, these
are cyclic images after ordering the coordinates as `(1,2,4,5,3,6)`.

The primary certificate is the direct bounded integer enumeration.  A
separately encoded 19-variable, 22-constraint CP-SAT audit uses the broad
`l1` coordinate domains, deliberately omits the prior `86,106` energy caps,
and adds only the forbidden conditions `||q||_2^2>=32` and `>=83`; with one
worker, seed zero, and a safely exhaustive degree-two quotient domain
`[-200,200]`, both models return exact status `INFEASIBLE`.  Their canonical
model-proto hashes under OR-Tools 9.15.6755 are

```
557e700271596217961bf7f5a6db8107bc32dbdea718e2f62e0cbf4ad8765db3
72df1b51c8f369bce8d4133491a74c1a290cb6b878031e014ec7c2b3fc3b0603
```

Independent DFS reconstruction reproduced both row sets, counts, maximizers,
and hashes without importing the primary enumerator.  The same certificate
was replayed on `jellyfin`, `orin`, `nuka`, and `soulkiller`; all four hosts
returned `5844 -> 30, max 31, INFEASIBLE` and
`1704 -> 24, max 82, INFEASIBLE`.

## Contradiction and scope

The ten nonexact rows would have energy at most

`3*31 + 7*82 = 667`.

But the common-graph identity requires at least

`707+26C >= 707`.

The gap is 40, so the integral difference system is empty.  This closes the
generic four-exact `p=13,t=3` partition.  Together with Proposition 15.739's
exceptional close, it closes `p=13,k=58` completely.

This is not a classification of directional coefficient matrices.  It
enumerates only a six-coordinate necessary superset and becomes decisive
only when its sharp row-energy maxima are compared with the common-graph
Parseval identity.  It uses neither the quartic value code, a root-quartet
orbit split, nor midpoint variables.  Residual (ii), multi-level Type I, and
the quadratic-minmax limit remain open; the remaining residual ranges include
critical `p=5,7`, `p=11,k>=50`, `p=13,k>=60`, and the generic `p>=17`
later layers beginning at `k=4p+6`.

Canonical executable evidence is generated by
`src/e1_gmin_m4_prop15742.py` and checked by
`tests/test_prop15742.py`.
