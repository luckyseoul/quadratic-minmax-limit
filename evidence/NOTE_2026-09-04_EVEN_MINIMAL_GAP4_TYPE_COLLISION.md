# Even minimal-gap-four: the two-direction-type collision

**Status:** proved infinite-band theorem; not a global bridge or residual-(ii)
closure.  No prime, graph, orbit, or eigenshell census is used.

Let `C` be the Paley conference matrix of order `p^2+1`, `p>=11`, and let
`H` have even size

```text
h=4p+2t.
```

Assume the branch left open by Proposition 15.764:

```text
T_H^eps(y)>=4  for every y in E_eps and eps in {+1,-1}.       (1)
```

The theorem excludes (1) in the following ranges:

```text
p in {11,13,17}:       1<=t<=3;
p>=19, p=3 mod 4:      1<=t<=(p+5)/4;
p>=19, p=1 mod 4:      1<=t<=(p+7)/4.                       (2)
```

In particular `h=4p+2` is impossible for every prime `p>=11`.

## 1. Exact directional identities

Throughout (2),

```text
p^2+1-2h>0,
h < 4+(p^2-1)/2.                                           (3)
```

Thus `H` has an isolated vertex.  Proposition 15.721 transports it to
infinity without changing (1), `h`, or the relative flip set, so all edges
of `H` are finite.

Put

```text
M=p+1,  m=M/2,  tau=sum_(uv in H) C_uv.
```

For a projective `F_p`-direction `d`, let `eps` be the quadratic type of its
kernel and let `P_d` count the edges of `H` parallel to `d`.  On the affine
middle-slice chart define

```text
A_d=(eps S_H-4)/2,       a_d=2p E_d[A_d].                  (4)
```

By (1), `A_d` is a nonnegative integer-valued quadratic on
`J(p,(p+1)/2)`.  A parallel edge has signed expectation one.  A nonparallel
edge has signed expectation `-eps C_uv/p`.  Hence

```text
a_d=M P_d-eps tau-4p.                                     (5)
```

There are `m` directions of each type, and every finite edge is parallel to
the unique direction of type `C_uv`.  Therefore

```text
sum_(d:type eps) P_d=(h+eps tau)/2,
sum_(d:type eps) a_d=M t.                                 (6)
```

In particular every type has average `2t`.  Also

```text
a_d == 4-eps tau (mod M),                                 (7)
```

and every `a_d` is even.

## 2. A positive chart has a uniform mass floor

Let `D` be the odd-degree boundary of `H`, and in direction `d` let `B_d`
be the affine fibres containing an odd number of points of `D`.  The edge
product identity makes `A_d mod 2` an affine parity with active-coordinate
set `B_d`.

If this parity is nonconstant, Proposition 15.750 gives probability at least
`(p-1)/(2p)` to its odd class, so

```text
a_d>=p-1.
```

If it is constantly odd, `A_d>=1`, so `a_d>=2p`.  If it is constantly even
and nonzero, write `A_d=2B_d'`.  Proposition 15.681 gives

```text
a_d=4p E[B_d'] >= beta_p,

beta_p = (p+1)/2  if p=3 mod 4,
         (p-1)/2  if p=1 mod 4.                            (8)
```

Both values of `beta_p` are even.  Thus universally

```text
a_d=0 or a_d>=beta_p.                                     (9)
```

## 3. Rigidity of a zero direction

If `a_d=0`, nonnegativity makes `A_d` identically zero.  Write the signed
sum of the nonparallel edges between fibres `s,t` as `K_st`.  Then

```text
P_d + sum_(s<t) K_st z_s z_t = 4   on sum_s z_s=1.
```

The Johnson one-swap lemma used in Proposition 15.750 makes all `K_st`
equal to an integer `kappa`.  Since `sum_(s<t)z_s z_t=-(p-1)/2`,

```text
P_d-(p-1)kappa/2=4.                                       (10)
```

If `kappa<=-1`, then `P_d<=4-(p-1)/2<0`.  If `kappa>=1`, the absolute block
sums require at least

```text
P_d + binom(p,2) kappa
  = 4 + ((p^2-1)/2) kappa
```

edges, contrary to (3).  Hence

```text
kappa=0,  P_d=4,  eps tau=4.                              (11)
```

The last equality is (5) with `a_d=0`.

## 4. Collision of the two quadratic types

Suppose first that neither type contains a zero direction.  Let `r_+` and
`r_-` be the representatives in `[0,M)` of the residues (7).  In (2),
`2t<M`.  A residue below `beta_p` cannot be used: the next nonnegative value
in the same congruence class is already above `M`, hence above the type
average.  Consequently (9) gives

```text
beta_p<=r_+,r_-<=2t,
r_++r_- == 8 (mod M).                                     (12)
```

Since `2 beta_p>8` and `r_++r_-<2M`, (12) forces

```text
r_++r_-=M+8.
```

But (2) is chosen exactly so that `4t<M+8`.  Therefore at least one type
has a zero direction.  By (11), the opposite type has residue eight.

For `p>=19`, `beta_p>8`; hence a positive integer in the residue-eight
class is at least `M+8`, larger than the average `2t` in (2).  The opposite
type cannot have all positive entries, while residue eight forbids a zero.
This is the contradiction.  For `p=11,13,17` and `t<=3`, the least positive
residue-eight entry is eight and `2t<8`, giving the same contradiction.

## 5. Why this does not close residual (ii)

For an odd residual-(ii) separator `H=G union {e}` with
`h=4p+2t+1`, the baseline is three.  The same calculation gives

```text
a_d=M P_d-eps tau-3p,
sum_(d:type eps) a_d=m(M+2t),
a_d == 3-eps tau (mod M).                                 (13)
```

The type average is now `M+2t`, not `2t`, so (9) does not force a zero.
This failure is exact.  For every `t>=2`, take `tau=1`.  In the plus type
use entries

```text
a_d=M+2+M q_d,   P_d=4+q_d,   sum q_d=t-1;
```

and in the minus type use

```text
a_d=M+4+M q_d,   P_d=4+q_d,   sum q_d=t-2.
```

All entries are positive, both sums in (13) are exact, and the parallel
counts sum respectively to `(h+1)/2` and `(h-1)/2`.  This is not claimed to
come from a common Paley graph.  It proves precisely that the zero-direction
collision cannot be iterated into a residual-(ii) closure without a new
common-graph realizability constraint.

Replay:

```bash
PYTHONPATH=src python src/e1_gmin_m4_prop15766.py
PYTHONPATH=src pytest -q tests/test_prop15766.py
```
