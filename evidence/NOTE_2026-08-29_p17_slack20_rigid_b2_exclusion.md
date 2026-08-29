# Proposition 15.707: all 78 p17 slack-twenty profiles are impossible

The pair-slack-twenty block splits by its exact residue pair:

```text
(u_0,u_1)=(0,8): 69 profiles
(u_0,u_1)=(8,8):  9 profiles.
```

The first block closes by extending Proposition 15.706's global Paley-sign
identity. The second closes geometrically using two undetermined directions
and only the arc classifications already audited in Propositions
15.700--15.705.

## The 69 `(0,8)` profiles

For phase `a`, every direction has scaled mean

```text
2*u_a + 18*q_d,
```

where the nine nonnegative quotients sum to `9-u_a`. At phase one the
quotient sum is one and every profile has eight or nine `b=2` directions.
Therefore at least eight attain their floor 16 and are rigid.

At phase zero the least admissible quotients are

```text
b:      0  2  4  6  8  10 12 14 16
q_min:  0  1  2  3  3   3  3  3  2.
```

If `E=9-sum_d q_min(b_d)`, at most `E` directions can be raised above their
minimum. Every `b=0` direction at quotient zero and every `b=2` direction at
quotient one is rigid. The exact lower-bound histogram for retained rigid
`b in {0,2}` directions is

```text
retained:  3  4   5   6  7
profiles:  2 10  26  28  3.
```

Thus all 69 profiles retain rigid directions of both quadratic types. For a
phase-one rigid `b=2` direction the global finite-edge sign identity is

```text
S = c_H*(5 + 144*g_+ - 17*I).
```

For either phase-zero rigid floor, `b=0` with `(M,T)=(0,0)` or `b=2` with
`(M,T)=(18,-1)`, it is

```text
S = -c_H*(3 + 144*g_- - 17*I).
```

Equating them gives

```text
17*I = 4 + 72*(g_+ + g_-).
```

Hence `I=68`; the one remaining finite edge gives affine odd-boundary size
66, 68, or 70 rather than 16. All 69 profiles are impossible.

## The nine `(8,8)` profiles

Their undetermined-direction histogram is `{2:5,3:4}`. Choose two such
infinity points. If deleting the minimum `r<=5` boundary points repairs the
boundary to an arc `A`, adjoining the two infinity points gives an arc `K`
of size `18-r`. No deleted point gains a `K`-secant through either infinity
point, because that would create a boundary chord in an undetermined
direction.

- If `r<=3`, then `|K|>=15`, so Sticker's exhaustive classification puts
  `K` on a conic. If `h>=1` original boundary points are off the conic, the
  undetermined conic point makes two omitted conic points destroy only one
  secant through each off-conic point. Thus slack is at least `4h(7-h)`,
  whose values for `h=1,2,3` are 24, 40, and 48. If `h=0`, the boundary is
  an arc and has slack zero.
- If `r=4`, an incomplete 14-arc reaches the conic branch. The unique
  complete 14-arc has minimum outside secant index two, forcing slack at
  least `4*4*2=32`.
- If `r=5`, equality forces five deleted points of secant index one. The
  eight complete 13-arcs have at most three such outside points. An
  incomplete 13-arc reaches either a conic or the unique complete 14-arc;
  deleting one point from that 14-arc leaves at most four index-one outside
  points. Both bounds contradict the required five.

Therefore all nine profiles are impossible. No new classification, solver,
or numerical timeout is used.

An independent native-XOR affine-Radon model deduplicates the 78 arithmetic
rows to 69 boundary signatures and provides a lossless normalized audit path.
It was implemented but not run as a proof search because the analytical and
finite-geometric argument already closes the block. No solver output is used
by this proposition.

The full p17 ledger drops from 639 to 561 profiles, all of pair slack at
least 24. The endpoint and every top-level gate remain open.
