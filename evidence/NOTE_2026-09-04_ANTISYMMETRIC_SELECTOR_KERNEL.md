# The physical-half selector is not an anti-Radon invariant

Date: 2026-09-04

**Status:** the mod-two degree boundary of the full signed source chain is
determined by its edge--Radon target, but the degree boundary of one physical
half of a ternary antisymmetric lift is not.  An explicit Type-K circuit
changes both the full projective-kernel selector word and its aggregate
parity.  Thus the adaptive selector obstruction for localized Mobius halves
cannot be promoted to arbitrary antisymmetric lifts by target invariance.
For the ordinary unsigned anti-Radon map, a second exact construction shows
more: ternary kernel halves realize every unit selector word, so the induced
linear signature code is the full space and has minimum weight one.  This
full-code statement is **not** asserted after Paley source signing; the
normalized code and its minimum weight remain open.  The `p=31` Type-K
barrier itself is normalized-valid.  These are barrier results, not a
residual-(ii) witness or closure.

## 1. The full signed-chain boundary is invariant

Let `R` be the ordinary integral edge--Radon map on
`V=F_p^2`, for an odd prime `p`, and let

\[
 \partial_2:\mathbb Z^E\longrightarrow\mathbb F_2^V
\]

be edge boundary after coefficient reduction modulo two.  Every Type-P
ridge has zero `partial_2`: on a midpoint fibre the two endpoint maps
`a -> a-delta` and `a -> a+delta` are permutations of the same affine line.
Every Type-K ridge also has zero `partial_2`: the relation

\[
 L(\delta_1)^2=L(\delta_2)^2
\]

says that its two difference classes have the same unordered pair of
endpoint fibres, so their boundaries cancel modulo two.

The ridge theorem gives

\[
 p\ker_{\mathbb Z}R\subseteq K_{\rm ridge}.
\]

Therefore, for every `z in ker_Z R`,

\[
 \partial_2z=\partial_2(pz)=0,                         \tag{1}
\]

because `p` is odd.  In particular (1) holds on the antisymmetric kernel.
The full signed-chain boundary, and hence every fixed family of selector
pairings with that boundary, is target-determined.

This statement is not the invariant needed by the localized-Mobius proof.
That proof pairs selectors with the boundary of the physical edge chosen
from each active inversion orbit.  Passing from the signed anti-chain to
that positive half loses (1).

## 2. An all-odd-prime countercircuit

Take

\[
 L(x,y)=x,\qquad \delta_1=(1,0),\qquad\delta_2=(1,1),
 \qquad g=\mathbf1_{1}-\mathbf1_{-1}.
\]

The two difference classes have equal nonzero squared `L`-projection.  The
corresponding Type-K ridge is

\[
 z=u^K_{L,\delta_1,\delta_2,g}=\mathbf1_H-\mathbf1_{-H}, \tag{2}
\]

where

\[
 H=
 \{((1,y),[\delta_1]):y\in\mathbb F_p\}
 \ \cup\
 \{((-1,y),[\delta_2]):y\in\mathbb F_p\}.               \tag{3}
\]

Here `(a,[delta])` denotes the edge `{a-delta,a+delta}`.  Thus `H` has
`2p` edges, `H` and `-H` are disjoint, (2) has `4p` nonzero actual-edge
coordinates, is ternary and antisymmetric, and

\[
                         Rz=0.                            \tag{4}
\]

The first family in (3) has boundary equal to the two affine lines `x=0`
and `x=2`.  The second has boundary equal to `x=-2` and `x=0`.  Consequently

\[
              \partial_2H=\mathbf1_{\{x=2\}}
                           +\mathbf1_{\{x=-2\}}.          \tag{5}
\]

It is centrally symmetric.  For each projective functional `N`, let `S_N`
contain zero and one representative of every nonzero antipodal pair on
`ker N`.  If `N=L`, (5) misses `ker N`, so its selector bit is zero.  If
`N != L`, the two lines in (5) meet `ker N` in exactly one antipodal point
pair, so its selector bit is one.  In the canonical order beginning with
`L`, the word and its aggregate are therefore

\[
 \boxed{q(\partial_2H)=(0,1,1,\ldots,1),\quad
        \operatorname{wt}q=p,\quad
        \sum_Nq_N=p=1\pmod2.}                             \tag{6}
\]

The empty half and (3) have the same zero antisymmetric Radon target by
(4), but their selector words and aggregate parities differ.  Neither is
target-determined on arbitrary ternary antisymmetric lifts.

Equation (1) remains consistent: the full chain has boundary
`partial_2 H + partial_2(-H)=0`, since (5) is central.

## 3. The unsigned physical-half code is the full space

There is also a ternary anti-kernel half whose selector word cancels all but
one of the bits in (6).  Continue to call the two coordinate functionals
`L(x,y)=x` and `M(x,y)=y`, and put `c=1/2`.  Define `C` as the union of

\[
 \binom{\{(c,y):y\ne0\}}2
 \quad\hbox{and}\quad
 \{\{(c,0),(-c,y)\}:y\ne0\}.                              \tag{7}
\]

Thus `C` is a clique on the `p-1` nonzero-height points of the line `x=c`,
together with a star from the omitted point `(c,0)` to the nonzero-height
points of the opposite line.

The chain `1_C-1_{-C}` lies in the unsigned anti-Radon kernel.  This follows
row by row.

- In direction `L`, the clique is parallel and every star edge projects to
  the inversion-fixed pair `{c,-c}`.
- In direction `M`, the clique and star together project to the complete
  graph on `F_p`.
- Write any remaining functional as `N=aL+bM` with `ab != 0` and put
  `u=ac`.  The clique projects to the complete graph on
  `F_p minus {u}`.  The star supplies every pair from `u` except
  `{u,-u}`, and supplies one parallel edge.  Hence the off-diagonal image is
  every label pair except the inversion-fixed pair `{u,-u}`, together with
  one parallel edge.  This row is again invariant under label negation.

Therefore `R(C)=R(-C)`.  Its mod-two boundary is

\[
 \partial_2C=\{(c,y),(-c,y):y\ne0\}.                       \tag{8}
\]

The corresponding selector word is zero at `L` and `M`, and one at every
other direction.  The transverse half `H` from Section 2 has word zero at
`L` and one everywhere else.  The two halves use disjoint inversion orbits:
the clique has `L(delta)=0`, its star has midpoint `L(a)=0`, while the two
transverse families have nonzero `L(delta)` and midpoint values `+1,-1`.
Consequently their union is still ternary and

\[
 q\bigl(\partial_2(H\mathbin\cup C)\bigr)=e_M.             \tag{9}
\]

Its physical-half and anti-chain support sizes are

\[
 |H\mathbin\cup C|=2p+{p(p-1)\over2}={p(p+3)\over2},
 \qquad |\operatorname{supp}(1_{H\cup C}-1_{-(H\cup C)})|
 =p(p+3).                                                   \tag{10}
\]

Choose any desired projective direction for `M` and any independent
direction for `L`.  The same coordinate construction gives its unit word.
Thus actual ternary anti-kernel halves realize all `p+1` standard basis
vectors.  In particular, for the **unsigned** convention,

\[
 \boxed{\text{induced linear selector code}=\mathbb F_2^{p+1},
        \qquad d_{\min}=1.}                                \tag{11}
\]

This is an all-odd-prime row proof, not an extrapolation from the finite
checks at `p=3,7,11`.

The conclusion (11) must not be silently transferred to normalized physical
graph coordinates.  Multiplying an anti-kernel vector by the Paley column
sign preserves the normalized kernel equation, but can reverse which member
of individual inversion orbits lies in its positive physical half.  No
normalized full-code or minimum-weight theorem is proved here.

## 4. Exact p=31 replay and normalized signs

At `p=31`, the positive half has 62 edges and the full circuit has support
124.  Direct sparse application of all 32 edge--Radon blocks returns the
zero target.  Its half-boundary has weight 62, the selector word has weight
31 and aggregate one, while the full-chain boundary has weight zero.

For the normalized residual convention, the Paley column sign is constant
on a spatial difference class.  Both classes used above have sign `+1` at
`p=31`:

\[
 \eta(1^2+0^2)=1,\qquad \eta(1^2+1^2)=\eta(2)=1.
\]

Hence the normalized source signing does not change this `p=31`
countercircuit or its physical positive half.

Replay with

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
  tests/test_antisymmetric_selector_kernel.py
```

The executable record is
`src/e1_gmin_m4_antisymmetric_selector_kernel.py`.  It checks the exact
Radon image, inversion antisymmetry, both boundary weights, every selector
bit, the aggregate parity, and the `p=31` Paley signs.

## 5. Consequence

The adaptive product theorem remains valid for its proved localized-Mobius
family.  What fails is the proposed extension step saying that its selector
aggregate depends only on the full antisymmetric target.  Arbitrary
anti-Radon kernel freedom already changes that aggregate through the
elementary circuit (2).  Any wider obstruction must use additional source
structure (for example the localized shape or coupled symmetric data), not
the antisymmetric Radon target alone.
