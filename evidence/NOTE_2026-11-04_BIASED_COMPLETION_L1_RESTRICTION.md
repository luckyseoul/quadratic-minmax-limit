# Biased-completion L1 rigidity for principal restriction

2026-11-04. Deterministic all-orders restriction theorem. The original
MathOverflow convergence problem remains OPEN. This note is not entered in
`src/original_mo_status.py` as a completion claim.

## Statement

Let a complete symmetric zero-diagonal signing be block-decomposed as

```text
A = [[B, C],
     [C^T, D]]
```

with principal blocks of orders n and t. Put

```text
F = Phi(A),   H = Phi(D).
```

For any Boolean state x on the B block, write

```text
q = Q_B(x),
L = ||C^T x||_1.
```

Then

```text
boxed:
F - |q| >= max_(0<=s<=1) { s L - s^2 H }.                  (1)
```

For H>0 this is equivalently

```text
F - |q| >=
    L^2/(4H),    if L <= 2H,
    L-H,         if L >= 2H.                               (2)
```

If H=0, (1) gives simply `F-|q| >= L`.

In particular, if x is an extremizer of B and

```text
delta = F - Phi(B) <= H,
```

then

```text
boxed:
||C^T x||_1 <= 2 sqrt(H delta).                            (3)
```

Thus a principal restriction that loses little Boolean norm forces every
chosen extremizer used in (3) to be an approximate Boolean L1-kernel vector
of the entire omitted cross block.

## Proof

Choose `sigma in {+-1}` so that `sigma q=|q|`. Let

```text
z_j = sign((C^T x)_j),
```

with an arbitrary sign when the coordinate vanishes. For `0<=s<=1`, let
the coordinates of a random Boolean vector Y on the D block be independent
with means

```text
E Y_j = sigma s z_j.
```

Then

```text
sigma E[x^T C Y]
  = sigma sum_j (C^T x)_j E Y_j
  = s ||C^T x||_1
  = sL.                                                     (4)
```

Independence gives, for every j<k,

```text
E[Y_j Y_k] = s^2 z_j z_k,
```

because the two factors of sigma cancel. Hence

```text
sigma E Q_D(Y)
  = sigma s^2 Q_D(z)
  >= -s^2 Phi(D)
  = -s^2 H.                                                 (5)
```

Therefore

```text
E[ sigma Q_A(x,Y) ]
  >= |q| + sL - s^2 H.                                     (6)
```

Every realization satisfies `sigma Q_A(x,Y)<=F`, so taking expectations
and rearranging proves (1). Maximizing the concave quadratic over `[0,1]`
gives (2).

For (3), take `|q|=Phi(B)`. If `L>2H`, (2) gives
`delta>=L-H>H`, contradicting `delta<=H`. Hence `L<=2H`, so
`delta>=L^2/(4H)`, which is exactly (3). QED.

## Relation to previous work

Repository deduplication was run against both
`luckyseoul/quadratic-minmax-limit` and
`Robby955/mo-413935-research` before this note was added.

The closest existing results are:

1. `NOTE_2026-09-24_ORTHOGONAL_COMPLETION_RESTRICTION.md`:
   ```text
   Q_B(x)^2 + binom(t,2) + ||C^T x||_2^2 <= F^2.
   ```
   This is an exact L2/second-moment completion identity.

2. `NOTE_2026-09-24_FOURTH_MOMENT_RESTRICTION_LOSS.md`:
   for a restriction extremizer,
   ```text
   Phi(B) <= F - c4 sqrt(binom(t,2)+||C^T x||_2^2).
   ```
   This yields L2 approximate-kernel rigidity.

3. `NOTE_2026-09-24_HYPERCONTRACTIVE_RESTRICTION_LOSS.md` gives the
   earlier weaker L2 one-sided completion estimate.

No equivalent L1 inequality (1)--(3), or the biased-completion proof above,
was found in either repository. The new statement is not a restatement of
the existing L2 estimates: the relation `||v||_1>=||v||_2` runs in the
wrong direction for deriving (3) from those upper bounds.

## Scope

This does not by itself prove convergence. To turn (3) into a cross-order
modulus one still needs a theorem excluding sufficiently good Boolean
L1-kernel extremizers for the cross blocks of exact minimizers at the
relevant scales, or an averaging/packing argument showing that such
extremizers cannot persist across the needed family of restrictions.

That is the exact remaining implication exposed by this note.
