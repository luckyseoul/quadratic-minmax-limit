# Outer-layer opposite-sign pair-cover defect

This is an all-orders consequence of the coefficient edge-flip witnesses.
It gives a quantitative reason that one positive and one negative outer-layer
state cannot supply the needed global-minimality cover.  It does not establish
the required growing multiplicity of such states, and does not close the RG2
rounding problem.

Let `A` be a complete signing with `M=Phi(A)`.  Suppose that Boolean states
`x,y` obey

```
Q_A(x) >= M-eps_plus,       Q_A(y) <= -M+eps_minus.
```

For an edge `e={i,j}`, write `s_x(e)=A_ij x_i x_j` and similarly for `y`.
The edge-flip condition uses a positive outer state only on edges with
`s_x(e)=-1`, and a negative outer state only on edges with `s_y(e)=+1`.

Put `r=x circ y`, let `d` be the Hamming distance of `x,y`, and let
`b=d(n-d)` be the number of cut edges of `r`.  On those cut edges set
`C_ij=A_ij x_i x_j`, and let

```
u = #{ {i,j}: r_i r_j=-1 and C_ij=+1 }.
```

Every one of these `u` edges is uncovered by *both* witnesses: it has
`s_x=+1`, while `s_y=C_ij r_i r_j=-1`.  Splitting the two energies into
the cut and its complement gives the exact identity

```
Q_A(y)-Q_A(x) = 2b-4u.                                    (1)
```

The outer-layer inequalities consequently imply

```
u >= b/2 + M/2 - (eps_plus+eps_minus)/4.                   (2)
```

Since `u<=b`, they also force the Hamming-cut separation

```
d(n-d) >= M - (eps_plus+eps_minus)/2.                      (3)
```

For the actual one-edge global-minimality witnesses, both errors are at most
two.  Thus every positive/negative outer-two pair has cut size at least
`M-2`, and its two witness sets leave at least `b/2+M/2-1` edges uncovered.
In particular, a two-state positive/negative cover is impossible whenever
`M>2`.

This strengthens the observation that one-edge minimality supplies only
witnesses rather than mass: the two opposite outer layers cannot collapse to
one mutually covering pair.  It still permits a larger, diffuse witness
family to cover all edges, so no anchor-capacity or RG2 conclusion follows.

## Same-sign pair defect

The defect is not peculiar to opposite signs.  Let `x,z` both be positive
outer states, with errors `eps_x,eps_z`, and use the same notation
`r=x circ z`, `b=d(n-d)`, and `a=binom(n,2)-b`.  On the `a` edges where
`r_i r_j=+1`, put `C_ij=A_ij x_i x_j` and let `u_+` count those with
`C_ij=+1`.  Such edges are uncovered by both positive witnesses.  The exact
identity is

```
Q_A(x)+Q_A(z) = 2(-a+2u_+),                                (4)
```

and therefore

```
u_+ >= a/2 + M/2 - (eps_x+eps_z)/4.                         (5)
```

For two negative outer states, define `u_-` to count the `a` agreement-side
edges with `C_ij=-1`.  Those are uncovered by both negative witnesses, and

```
Q_A(x)+Q_A(z) = 2(a-2u_-),                                  (6)
u_- >= a/2 + M/2 - (eps_x+eps_z)/4.                         (7)
```

Thus every pair of outer-two witness states has a certified uncovered set:
the disagreement cut for an opposite-sign pair and the agreement side for a
same-sign pair.  Pairwise defects alone do not rule out a cover by three or
more witnesses, so this remains a boundary lemma rather than a multiplicity
theorem.

## Exact sign regression

`tests/test_outer_layer_pair_cover_defect.py` exhausts every signing through
order four and every ordered state pair, checking (1) directly.  It passed on
the controller under pytest and independently on NUKA under Python 3 by
calling the same test function directly (NUKA does not have pytest installed).
