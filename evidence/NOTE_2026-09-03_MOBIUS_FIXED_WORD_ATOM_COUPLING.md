# The branch-C `j=0` fixed-word/atom coupling

Date: 2026-09-03

Status: proved an exact target-sensitive necessary-and-sufficient criterion
for the **fixed-word and central odd-moment layer** of the all-active
branch-C `j=0` endpoint. It forces a new surplus of collisions in the
affine-block image of the actual Möbius support. Two distinct localized
halves share at most eight block types, but a direct pair can share three,
so the hoped-for bound one is false. The correct bounds leave room and do
not settle the even common moments, nonfixed target cells, the endpoint, or
residual (ii).

## 1. Use the actual target atoms

Put

\[
 p=4r+3,\qquad m={p+1\over2}=2r+2,\qquad h={p-1\over2}=m-1.
\]

In a hard direction `D`, let `e_D` be its balanced compact count and put

\[
                         n_D=e_D.
\]

In an opposite direction the target has exactly `r-1` all-equal
`p-3` atoms and `Q_D-r-2` compact `p+1` atoms. Thus its total triangle
count is

\[
                         n_D=Q_D-3.                       \tag{1}
\]

These are the atoms from Proposition 15.758, with their actual signs. After
reduction modulo two, both a compact signed triangle and an all-positive
triangle are the same three-cycle. The signs therefore disappear, but the
distinction between the two atom counts in (1) has not been discarded.

Writing `E=t+1`, the exact balanced totals give

\[
 \sum_Dn_D=E+\{E+5m-5-3m\}=2E+2m-5,                    \tag{2}
\]

which is odd. Throughout the branch-C interval every `n_D` is positive;
hard counts are at most `h-3` and opposite counts at most `h-2`.

## 2. Exact fixed-word coordinates

Let `M` be the square incidence matrix whose columns are the paired
nonorigin affine blocks. For the actual ternary Möbius support `U`, define

\[
 c_U(D,\beta)=\#\{O\in U:\Phi(O)={\cal B}(D,\beta)\}\pmod2. \tag{3}
\]

This vector is indexed by all `d h=(p+1)(p-1)/2` block types. Let
`ell_D` be the singleton coordinate `beta=j_D^2` in a hard row and zero in
an opposite row. If the one fixed endpoint edge is `[x]`, let `F=L_x` be
its parallel direction and put

\[
 s_x(D,\beta)=
 \begin{cases}
  1,&D\ne F,\ \beta=D(x)^2,\\
  0,&D=F.
 \end{cases}                                             \tag{4}
\]

For each target atom in row `D`, record whether its three distinct labels
contain an antipodal pair `{s,-s}`. A triangle contains at most one such
pair. Let `z_D(beta)` be the parity of the atoms whose antipodal pair has
square `beta`.

The explicit fixed-edge inverse is

\[
 a_{[v]}=g_{L_v}(0)+\sum_Lg_L(L(v)^2).
\]

The `P_(L_v)=g_(L_v)(0)` term occurs once outside and once inside the sum,
so it cancels pointwise. The hard literal contributes its affine block
`B(D,j_D^2)`, and the fixed cells of the triangle atoms contribute `M z`.
Consequently

\[
                         M^Ta_Y=\ell+z.                  \tag{5}
\]

Also `M^T Phi(U)=c_U` and `M^T e_[x]=s_x`. Since `M^TM=I`, the endpoint
equality is therefore exactly

\[
 \boxed{a_Y+\Phi(U)=e_{[x]}
        \quad\Longleftrightarrow\quad
        z=c_U+\ell+s_x.}                                \tag{6}
\]

This is the missing coupling to the actual `U`; neither centrality nor a
common-moment equation has been used as a substitute for it.

## 3. The atom-capacity syndrome

Centrality pairs every nonfixed label edge with its negative. The total
coefficient parity of `n_D` triangles is `n_D`, while all paired nonfixed
edges contribute evenly. Hence

\[
                     |z_D|\equiv n_D\pmod2.              \tag{7}
\]

Since one atom contains at most one antipodal pair, `|z_D|<=n_D`.
Conversely these two conditions are sufficient at this layer when the atom
labels are allowed to be chosen, subject to the prescribed type counts. For
every one in `z_D`, use a centered compact atom `K(s,-s;0)` or the centered
all-equal triangle `{s,-s,0}` of the prescribed type. Because
`n_D-|z_D|` is even, fill the remaining atom slots in equal pairs of square
classes. Every chosen atom is individually central and annihilates every
odd moment. This is an existential labeling statement with the exact
all-equal/compact count split; it is not sufficiency for an already fixed
list of atom labels.

Thus (6) gives the fully quantified row criterion

\[
 \boxed{
 |c_U(D,\cdot)+\ell_D+s_{x,D}|\le n_D,
 \qquad
 |c_U(D,\cdot)+\ell_D+s_{x,D}|\equiv n_D\pmod2
 }                                                        \tag{8}
\]

for every direction `D`. Under the exact `j=0` parallel slices, the parity
in (8) is automatic. Indeed, if `c_D` is the parity of the number of
auxiliaries in direction `D`, then the nonzero-block count in group `D`
has parity `c_D`: in one half every midpoint direction occurs twice, except
that the `t=0` zero-`Phi` orbit is removed together with its partner
`t=-2`, and that partner lies in the auxiliary direction. Hence precisely
the auxiliary group is toggled. Moreover,

\[
 c+1_{\rm hard}=u=P+e_F,
 \qquad P_D\equiv1+n_D,
 \qquad 1_{D\ne F}=1+e_F(D).
\]

Their sum is exactly `n_D`.

Equation (8) is necessary and sufficient only for the central fixed-word
and odd-moment layer. The centered construction need not satisfy the even
common moments or the prescribed nonfixed target coefficients.

## 4. A forced surplus of shadow block collisions

Every localized half has one zero-`Phi` orbit and `p-2` distinct nonzero
block types. Modulo two, deleting any physical cancellation changes two raw
occurrences and therefore does not change (3). Put

\[
 R_\Phi=m(p-2),\quad C=|c_U|,\quad
 \Lambda={R_\Phi-C\over2}.                              \tag{9}
\]

Let

\[
 q(x)=\#\{D\text{ hard}:D\ne F,\ D(x)^2=j_D^2\}.
\]

The words `ell` and `s_x` have weights `m` and `2m-1`, with exactly `q`
common coordinates. Summing (8), using (2), and applying the triangle
inequality gives

\[
 C\le(2E+2m-5)+(3m-1-2q).
\]

For `kappa_0=t_max-t+1`, direct substitution in (9) yields

\[
 \boxed{\Lambda\ge\kappa_0+m+q.}                        \tag{10}
\]

This has an exact support interpretation even under arbitrary triple and
higher overlaps. Let `kappa_1` and `kappa_z` be the physical cancellation
units on nonzero- and zero-`Phi` orbits. Then

\[
 \kappa_1+\kappa_z=\kappa_0,
 \qquad |U_{\rm np}|=R_\Phi-2\kappa_1.
\]

For each block `b`, the number of distinct surviving physical `U`-orbits
above `b` has parity `c_U(b)`. Therefore

\[
 \sigma={|U_{\rm np}|-C\over2}\in\mathbf Z_{\ge0},
 \qquad
 \Lambda=\kappa_1+\sigma
         =\kappa_0-\kappa_z+\sigma.                      \tag{11}
\]

Combining (10)--(11) gives the target-sensitive collision demand

\[
 \boxed{\sigma\ge\kappa_z+m+q.}                         \tag{12}
\]

Thus the `j=0` endpoint needs at least `m+q` pairs of **distinct surviving
orbits** which have the same nonzero `Phi` block, beyond the physical
cancellations themselves. Zero-`Phi` cancellations increase that demand.

## 5. The correct two-half block bound is eight

For one half with target/auxiliary basis `(L,M)`, center `j`, and
`z=t+1`, direct midpoint/difference calculation gives the normalized block
functional

\[
 N_z=-{z+1\over j}L+{z^2\over j(z-1)}M,
 \qquad z\ne0,1.                                       \tag{13}
\]

Its block is `{[v]:N_z(v)^2=1}`. If `e,f` is the primal basis dual to
`L,M`, set `A=je`, `B=j(e+f)`. Equation (13) is equivalently

\[
                         (N(A)+2)N(B)+1=0,               \tag{14}
\]

with the point `N(A)=N(B)=-1` omitted. This irreducible affine hyperbola
contains no pair `N,-N`; moreover `u=N(A)+2=1-z` recovers `z`. Hence the
`p-2` nonzero block types of one half are distinct.

Two halves share a block exactly when `N_i=+N_k` or `N_i=-N_k`. On the
first hyperbola use `u=N(A_i)+2`; any linear evaluation has the form

\[
                         N(V)=a(u-2)-{b\over u}.
\]

After multiplication by `u^2`, either signed equation for the second
hyperbola is a polynomial of degree at most four. It cannot vanish
identically: equality of the irreducible conics first identifies `B` from
the linear term and then `A` from the quadratic term, forcing the same
target direction because `L` annihilates `B-A`. Thus there are at most four
common blocks for each sign and at most eight total.

The stronger proposed bound one is false. Over `F_31`, the two halves

\[
 ((L,M),j)=(((1,1),(1,3)),2),\qquad
 (((0,1),(1,7)),3)
\]

have Paley-hard target directions, no common physical edge, and a ternary
sum, but their `Phi` blocks agree at the three parameter pairs

\[
                         (9,20),\ (19,12),\ (25,18).     \tag{15}
\]

This is a fail-when-wrong counterexample to a proof shortcut, not a census.

Charging pairs of surviving orbits to pairs of halves gives only

\[
                         \sigma\le8{m\choose2},          \tag{16}
\]

and the raw-occurrence bound is sharper in the relevant scalar comparison.
Even at `t=t_min` and the worst formal `q=m`, the remaining raw margin is

\[
 {m(2m-3)\over2}-\kappa_{0,\max}-m-q
 ={(m-1)(m-4)\over2}>0.                                \tag{17}
\]

Therefore (10)--(12) are a genuine new geometric gate, but no scalar
contradiction. A closure must control the global intersections of all `m`
hyperbolas together with the even target atoms and nonfixed cells.

## Reproduction

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_mobius_fixed_word_atom_coupling.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_mobius_fixed_word_atom_coupling.py
```
