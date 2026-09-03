# Prescribed centers on an equal-square common block

Date: 2026-09-03

Status: proved an exact Hall formulation for the fully doubled
equal-square cover, a deficiency-two obstruction that also excludes the
one-single incidence profile, and a one-way construction from any
prescribed anchor line whose link graph is a pseudoforest.  A direct
$p=31$ prescribed nonzero hard-center list realizes the obstruction for
one fixed common direction $K$.  It is not proved that some prescribed
anchor line always has the required pseudoforest link graph.  The results
also do **not** exclude a cover after changing $K$, using unequal square
classes, or imposing a different halved-code word.  Residual (ii) remains
open.

## 1. Equal squares turn the intrinsic criterion into anchors

Fix independent projective functionals (K,L_i), a nonzero (r), and the
common midpoint and difference blocks

\[
       X=B_{K,r^2},\qquad C=B_{K,r^2}.
\]

Represent every class of (X) by its unique point (a) satisfying
(K(a)=r).  For the prescribed nonzero center (j_i), let

\[
 L_i(z_i)=0,\qquad
 L_i(a_i^+)={j_i\over2},\qquad
 L_i(a_i^-)=-{j_i\over2}.                                \tag{1}
\]

These are three distinct points because (L_i) restricts to an affine
bijection on (K(a)=r).  Put

\[
 A_i=\{a_i^+,a_i^-\},\qquad
 T_i=\{z_i,a_i^+,a_i^-\}.                                \tag{2}
\]

The intrinsic criterion (22) in the affine-block-cover note is

\[
 l_1l_2(l_1+l_2)\ne0,qquad
 (r-\eta r)l_1-(r+\eta r)l_2=rj_i.                       \tag{3}
\]

For (eta=+1), equation (3) says (l_2=-j_i/2).  For
(eta=-1), it says (l_1=j_i/2).  The nonzero-product condition and
distinctness then say that the other point is outside (T_i).  Hence a
localized half with target ((L_i,j_i)) can supply two distinct midpoint
classes over (C) if and only if its unordered pair is

\[
             \boxed{\{x_i,y_i\},\quad
                    x_i\in A_i,\quad y_i\in X\setminus T_i.} \tag{4}
\]

For every pair in (4), the auxiliary is the unique functional determined
by

\[
              M_i(a)=L_i(a)-{j_i^2\over4L_i(a)}
\]

at its two points.  Thus (4) is both necessary and sufficient; it is not
only a counting relaxation.

## 2. Exact Hall criterion for the fully doubled profile

Let (m=(p+1)/2).  Give half (i) two left slots, with neighborhoods

\[
        N(i,\mathrm{anchor})=A_i,\qquad
        N(i,\mathrm{free})=X\setminus T_i.                \tag{5}
\]

If all (m) halves are doubled, they give (2m=p+1) incidences on the
(p)-point set (X).  They cover (X) exactly when every point is used
once except one point (d), which is used twice.  Duplicate (d) on the
right side of (5).  The cover is then exactly a perfect matching.

Consequently a fully doubled saturated cover exists if and only if there
is (d\in X) such that, for every (P,Q\subseteq\{1,\ldots,m\}),

\[
 \left|\left(\bigcup_{i\in P}A_i\right)
 \cup\left(X\setminus\bigcap_{i\in Q}T_i\right)\right|
 +\mathbf1_{\{d\text{ belongs to this union}\}}
 \ \ge\ |P|+|Q|.                                        \tag{6}
\]

When (Q=\varnothing), the second term in the union is omitted.  This is
ordinary Hall applied to (5), with the duplicate copy of (d), so (6) is
an exact necessary-and-sufficient criterion.

The anchor-only specialization gives the immediate necessary condition

\[
                  \left|\bigcup_{i\in P}A_i\right|
                         \ge |P|-1                       \tag{7}
\]

for every (P).

## 3. A deficiency-two obstruction covers both profiles

Suppose some subfamily (P) has

\[
                  \left|\bigcup_{i\in P}A_i\right|
                         \le |P|-2.                       \tag{8}
\]

If every half is doubled, the halves in (P) already lose at least two
of the (p+1) incidences to repeated anchors, leaving at most (p-1)
distinct midpoint classes.

The only other profile capable of covering (p) classes has exactly one
single half and (m-1) doubled halves.  At most one member of (P) can
be the single half.  Hence at least (|P|-1) doubled halves choose anchors
inside a set of size at most (|P|-2).  A repetition is forced, while
this profile has exactly (p) total incidences and therefore cannot
tolerate any repetition.  Thus (8) excludes both possible saturated-cover
profiles.

## 4. Direct (p=31) prescribed-center obstruction for one (K)

Work over (mathbf F_{31}).  Take

\[
                 K=(1,4),\qquad J=(1,0),                 \tag{9}
\]

and write (a(q)) for the point with (K(a(q))=1) and (J(a(q))=q).
The Paley sign of (K) is (-1), so all sixteen Paley-hard directions
are independent of (K).

Take the four block coordinates

\[
                         0,\ 2,\ 11,\ 14.                \tag{10}
\]

Their six pairwise midpoints are

\[
                         1,\ 7,\ 8,\ 21,\ 22,\ 28.       \tag{11}
\]

For an edge ({u,v}) in (10), let (L_{uv}) be the canonical
functional annihilating (a((u+v)/2)), and prescribe

\[
                         j_{uv}=2L_{uv}(a(u)).             \tag{12}
\]

Direct arithmetic gives:

| anchors ({u,v}) | zero (q) | hard (L_{uv}) | (j_{uv}) |
|---|---:|---:|---:|
| ({0,2}) | 1 | ((0,1)) | 16 |
| ({0,11}) | 21 | ((1,29)) | 30 |
| ({0,14}) | 7 | ((1,15)) | 23 |
| ({2,11}) | 22 | ((1,16)) | 27 |
| ({2,14}) | 8 | ((1,9)) | 15 |
| ({11,14}) | 28 | ((1,3)) | 7 |

Each displayed direction has Paley sign (+1), the six directions are
distinct, and every center is nonzero.  Set the centers of the other ten
hard directions equal to (1).  This is a complete prescribed nonzero
hard-center list.

The six displayed anchor sets are the edges of (K_4) on (10), so their
union has size four.  Equation (8) holds with equality:

\[
                             4=6-2.                       \tag{13}
\]

Therefore no choice of auxiliaries can give a saturated cover of
(B_{K,1}) over (C=B_{K,1}) for this center list.

This is a fixed-(K), equal-square obstruction.  It disproves automatic
saturation for arbitrary prescribed centers on a proposed common block.
It does not show that the same center list fails for every other common
direction (K'), nor does it address unequal square classes, mutual
ternarity, or containment of a full dual word.

## 5. A good prescribed anchor line gives a cover

There is a complementary one-way construction that exposes the exact
remaining quantifier.  Assume explicitly that every prescribed hard center
$j_i$ is nonzero.  Choose one hard target $(L_0,j_0)$ and put

\[
 K=L_0,\qquad \alpha=\beta=(j_0/2)^2,
 \qquad \ell=\{a:L_0(a)=j_0/2\}.                         \tag{14}
\]

The $L_0$-localized half is now dependent on the common direction.  At
$z=1$, equivalently $t=0$, its edge has $L_0,M_0$ coordinates

\[
                         \{(j_0,0),(0,0)\}.              \tag{15}
\]

Given any $a\in\ell$, choose $M_0$ to annihilate $a$.  Then (15) is the
physical edge $\{2a,0\}$, whose midpoint is $a$ and whose difference is
$-a$.  It therefore supplies an arbitrary single projective midpoint
class in $B_{K,\alpha}$ over the same difference block.  This common-block
parameter is unique: the midpoint condition gives $z^2=1$, while the
difference condition gives $(z-2)^2=1$.  Their only common solution in odd
characteristic is $z=1$.

For each of the other $h=(p-1)/2$ hard directions, define on $\ell$

\[
 \begin{aligned}
 z_i&:\ L_i(z_i)=0,\\
 A_i&=\{a:L_i(a)=j_i/2\text{ or }-j_i/2\},\\
 T_i&=A_i\cup\{z_i\}.
 \end{aligned}                                           \tag{16}
\]

The zero points are distinct.  Indeed, a common $z_i=z_k$ would be a
nonzero vector annihilated by two distinct projective functionals
$L_i,L_k$, which is impossible in dimension two.  Each $T_i$ is a
nondegenerate three-term affine progression on $\ell$, centered at $z_i$.
In characteristic other than three such a three-set has a unique
progression center, so the $T_i$ are also distinct.

Form the graph $G_\ell$ with vertex set $\ell$ and one edge $A_i$ for
every $i\ne0$.  It is simple: two equal anchor edges would have the same
midpoint $z_i$, contradicting the preceding paragraph.  Assigning distinct
anchors to the $h$ doubled halves is an edge-to-incident-vertex SDR.  Hall
says that this is possible exactly when every edge subfamily has at least
as many incident vertices.  For a graph, that is equivalent to

\[
              \boxed{\text{every component of }G_\ell
                     \text{ is a tree or is unicyclic}.} \tag{17}
\]

Thus the exact anchor condition is that $G_\ell$ be a pseudoforest.

Assume (17) and choose $h$ distinct anchors.  Their complement
$R\subset\ell$ has $h+1$ points.  The free slot of half $i$ has
neighborhood $R\setminus T_i$.  For a set $Q$ of free slots,

\[
 \left|\bigcup_{i\in Q}(R\setminus T_i)\right|
       =|R|-\left|R\cap\bigcap_{i\in Q}T_i\right|.        \tag{18}
\]

Hall follows in three exact ranges.

- If $|Q|\le h-2$, the intersection in (18) has size at most three, so
  the union has size at least $h-2\ge|Q|$.
- If $|Q|=h-1$, at least two distinct triples occur.  Their common
  intersection has size at most two, so the union has size at least
  $h-1$.
- If $|Q|=h$, the common intersection has size at most one.  To see this,
  fix two hypothetical common points.  There are at most three distinct
  nondegenerate three-term progressions containing that pair, one for each
  possible choice of the progression center.  Here $h\ge15>3$, and all
  $T_i$ are distinct.

The $h$ free slots therefore match into $h$ points of $R$.  The dependent
half (15) supplies the one remaining point.  We have proved the conditional
implication

\[
 \boxed{G_\ell\text{ pseudoforest}
        \ \Longrightarrow
        \text{a saturated equal-square common-block incidence cover}.} \tag{19}
\]

This is only a one-way construction.  The open step is to prove that some
one of the $p+1$ prescribed anchor lines has a pseudoforest link graph, or
to construct centers for which every base line fails.  Equation (19) does
not establish either assertion.  It also does not prove mutual ternarity
or containment of a full halved dual support.

## 6. Focused replay

The implementation verifies the displayed arithmetic, the dependent
$t=0$ singleton, representative pseudoforest and bicycle graphs, and the
theorem ledger.  It performs no search over primes or center lists.

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_prescribed_center_common_block.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_prescribed_center_common_block.py
