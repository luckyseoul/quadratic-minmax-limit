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

## Homogeneous three-witness obstruction

There is a first genuinely multi-state consequence.  Let `x,z,w` be three
positive states with respective errors `eps_x,eps_z,eps_w`.  Gauge by `x` and
put `r=x circ z`, `s=x circ w`.  The vertices split into four signature cells
according to `(r_i,s_i)`.  Let

```
N = sum_cells binom(|cell|, 2).
```

Thus `N` is the number of edges on which all three state-products agree.  On
such an edge write `C_ij=A_ij x_i x_j`.  If the three positive witness sets
cover every edge, then necessarily `C_ij=-1` on this whole common-agreement
region.  Indeed `C_ij=+1` makes all three positive witness conditions fail.

The exact signed-energy expansion is

```
Q_A(x)+Q_A(z)+Q_A(w) = sum_{i<j} C_ij (1+r_i r_j+s_i s_j).  (8)
```

The coefficient has magnitude three on the common-agreement region and
magnitude one elsewhere.  Its sign is forced against the positive maximum
there, so edgewise maximization in (8) gives

```
Q_A(x)+Q_A(z)+Q_A(w) <= binom(n,2)-4N.                     (9)
```

Convexity of `binom(t,2)` gives

```
N >= n^2/8-n/2,     hence     binom(n,2)-4N <= 3n/2.        (10)
```

Combining (9) with the three outer-layer lower bounds proves

```
eps_x+eps_z+eps_w >= 3M-3n/2.                              (11)
```

More generally, let `W_+` count the common-agreement edges with `C_ij=+1`.
They are precisely the edges jointly uncovered by the three positive
witnesses.  Keeping their actual sign in the edgewise estimate gives

```
W_+ >= (3M-(eps_x+eps_z+eps_w)-3n/2)/6.                   (11a)
```

For outer-two states this is `W_+ >= M/2-n/4-1`.  The identical statement
holds for three negative states after reversing signs.

The same proof after changing all signs applies to three negative states.
Consequently, three outer-two witnesses of one sign cannot cover every edge
as soon as `M>n/2+2`.  This is an all-orders obstruction to homogeneous
three-state covers.  By itself, (11) does not address mixed-sign triples or
provide the growing witness-family theorem needed for RG2.

## Mixed-sign triple obstruction

In fact the apparent mixed-sign exception is also impossible at the
outer-two scale.  Let `x` be positive outer and `y,z` negative outer, with
errors `eps_x,eps_y,eps_z`.  Gauge by `x`, write `r=x circ y`, `s=x circ z`,
and set `w=x circ y circ z`.  On an edge where both relative products are
negative, `r_i r_j=s_i s_j=-1`, the three witnesses all fail exactly when
`C_ij=A_ij x_i x_j=+1`.  Thus a three-state cover forces `C_ij=-1` on that
region.  If `N_{--}` is its number of edges, the four-state Walsh identity
is exact:

```
Q_A(x)-Q_A(y)-Q_A(z)+Q_A(w) = -4N_{--}.                  (12)
```

The three outer inequalities and `Q_A(w)>=-M` consequently give

```
eps_x+eps_y+eps_z >= 2M+4N_{--}.                          (13)
```

Without assuming a cover, let `W_mix` count the double-negative edges with
`C_ij=+1`.  These are jointly uncovered by the mixed triple.  The left side
of (12) is then `-4N_{--}+8W_mix`, so `Q_A(w)>=-M` gives

```
W_mix >= M/4 + N_{--}/2 - (eps_x+eps_y+eps_z)/8.           (13a)
```

Thus every mixed outer-two triple leaves at least `M/4-3/4` jointly
uncovered edges (and more when the double-negative region is large).

The same statement, with signs reversed, holds for two positive and one
negative state.  In particular a mixed-sign triple of outer-two witnesses
cannot cover every edge whenever `M>3`.  Together with (11), **no triple of
outer-two one-edge witnesses covers every edge once `M>n/2+2`**.  This is a
uniform all-orders cardinality lower bound of four on a putative witness
cover.  Equations (11a) and (13a) also give linear-in-`M` triple-defect mass
bounds.  They remain far short of the growing-family bound required by the
RG2/Banaszczyk route: a fourth witness can in principle cover a set of this
size.

## Exact sign regression

`tests/test_outer_layer_pair_cover_defect.py` exhausts every signing through
order four and every ordered state pair/triple, checking the displayed Walsh
identities directly.  It passed on the controller under pytest and
independently on NUKA under Python 3 by calling the same test functions
directly (NUKA does not have pytest installed).
