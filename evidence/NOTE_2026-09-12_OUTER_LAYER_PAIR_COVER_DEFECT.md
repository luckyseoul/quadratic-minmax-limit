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

## Exact sign regression

`tests/test_outer_layer_pair_cover_defect.py` exhausts every signing through
order four and every ordered state pair, checking (1) directly.  It passed on
the controller under pytest and independently on NUKA under Python 3 by
calling the same test function directly (NUKA does not have pytest installed).
