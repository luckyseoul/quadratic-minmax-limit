# Equianharmonic component packing and the exact-fiber threshold

## 1. Scope

Let \(p=4r+3\ge31\), and consider the constant tangent-conic atom fiber
after the normalization \(X=2x\).  In the only surviving nonstar case,

\[
 q^2+q+1=0,\qquad
 \Phi(X)=qX+1-q,\qquad
 \Psi(X)=qX+q-1.
 \tag{1}
\]

Thus \(p\equiv7\pmod {12}\), \(r\equiv1\pmod3\), and the target consists
of the \((p-1)/3\) nonfixed \(\Phi\)-triangles with the unique
self-antipodal edge deleted.  Its signed alignment score is \(p-2\).

The Reed--Solomon reduction leaves constants \(n=+1\) and \(n=-1\).
These are the same case for the present argument.  Apply \(x\mapsto-x\) to
every label of every atom, preserving the distinguished label in each
compact triple.  Every physical edge is then replaced by its antipode, so
in a fixed oriented edge-orbit basis the entire signed atom chain is
multiplied by \(-1\).  All-equal and compact atom types are preserved, while
supports, pairing components, and deficits are unchanged.  This is a
bijection between the \(-1\) and \(+1\) constant fibers.  We normalize to
\(n=+1\) below, and every conclusion covers both constant signs.

This note proves two statements.

1. An exact equianharmonic fiber containing \(r-1\) all-equal atoms and
   \(b\le r\) compact atoms necessarily satisfies

   \[
   b\ge {2r+7\over3}.                                  \tag{2}
   \]

2. Combining (2) with the already-proved line, two-line, boundary-cubic,
   and conic dichotomies proves zero-odd row centrality when
   \(3b\le2r+4\).  If \(p\equiv11\pmod {12}\), the equianharmonic case does
   not exist, and the same dependencies cover every \(0\le b\le r\).

These are odd-Radon statements.  They do not exclude nonzero global forms,
construct a common edge lift, solve the degree-six/eight conditions, or
close residual (ii).

## 2. Pairing graph and its exact component identity

For an atom, let \(A,R,O\) count aligned target, reverse target, and other
or self-antipodal signed edge occurrences.  The alignment score is \(A-R\).
Set

\[
 \delta_{AE}=3-(A-R)=2R+O,\qquad
 \delta_K=2-(A-R)=2R+O-1.                 \tag{3}
\]

The second deficit is nonnegative.  Indeed, a formal score-three compact
atom would satisfy

\[
 b=\Phi^i(a),\qquad c=\Psi^j(a)=\Psi^\ell(b),\qquad
 i,j,\ell\in\{1,-1\}.                       \tag{4}
\]

Solving the eight affine cases under \(q^2+q+1=0\) gives four repeated-label
solutions, two inconsistent systems, and two distinct-label solutions for
which the purported third aligned edge is the deleted self-antipodal edge.
Those last two atoms have score two, not three.  Hence no distinct-label
compact atom has score three.  This also excludes an otherwise isolated
\(K=1,Z=0\) negative-deficit component.

Pair every reverse target occurrence with an aligned occurrence at the same
target coordinate.  Exact unit target multiplicity leaves exactly one
unpaired aligned occurrence at each coordinate.  Pair nonself off-target
occurrences of opposite orbit sign.  Regard a self-antipodal occurrence as
a cap, and make atoms the vertices of the resulting pairing multigraph.
In particular, the unpaired aligned-coordinate supports of two distinct
components are disjoint: an overlap would give target multiplicity at least
two unless a reverse occurrence were paired to it, and that pairing would
join the components.

For one connected component, write \(K,E,m,\mu,Z\) for its compact vertices,
all-equal vertices, pairing edges, cycle rank, and caps.  If \(\delta\) is
the sum of (3), then

\[
\begin{split}
 \delta
 &=2R+O-K\\
 &=2m+Z-K\\
 &=K+2E-2+2\mu+Z.                           \tag{5}
\end{split}
\]

Here \(m=K+E-1+\mu\).  Formula (5) is the bookkeeping invariant used below;
it includes multiplicities and reverse-target cancellations.

## 3. Positive-excess components

Call \(K-2\delta\) the component's excess.  Formula (5) gives the exact
identity

\[
 K-2\delta=4-K-4E-4\mu-2Z.                    \tag{5a}
\]

The only nonnegative integer tuples \((K,E,\mu,Z)\) for which the right-hand
side is positive are

\[
 (1,0,0,0),\quad(1,0,0,1),\quad
 (2,0,0,0),\quad(3,0,0,0),                    \tag{5b}
\]

with excesses \(3,1,2,1\), respectively.  The first tuple would have
\(\delta=-1\), equivalently an isolated score-three compact atom, and was
excluded in Section 2.  Thus every component outside the remaining three
types satisfies

\[
 K\le2\delta.                                      \tag{6}
\]

### 3.1 HH component

In a capless all-high compact tree, every vertex has one outside port and no
reverse port, so connectedness forces two high atoms whose defects cancel
antipodally.  The valid discrete solutions have

\[
 (K,\delta,K-2\delta)=(2,0,2).                 \tag{7}
\]

The one-parameter solutions instead repeat aligned target coordinates and
cannot be components of an exact unit fiber.

### 3.2 One-low three-atom component

In an all-compact tree with exactly one low atom, positive excess requires
\(\delta=1\).  Equation (3) leaves exactly

\[
 (R,O)=(0,2)\quad\hbox{or}\quad(R,O)=(1,0).      \tag{8}
\]

The first is an O2 star: the two outside ports of the low atom pair with two
high leaves.  In the second case, the two high leaves form the HH component,
and the low reverse occurrence cancels one HH aligned edge.  Such a low atom
must be a full \(\Psi\)-triangle: its positive edge is reverse and its two
negative edges are aligned.  It therefore replaces the cancelled edge by
the other two edges of the same \(\Phi\)-triangle.  Both cases in (8) have

\[
 (K,\delta,K-2\delta)=(3,1,1).                 \tag{9}
\]

### 3.3 Capped singleton

A one-cap component has positive excess only when it is a single high atom:

\[
 (K,\delta,K-2\delta)=(1,0,1).                 \tag{10}
\]

All remaining equality components, including O3 stars and one-low trees of
deficit at least two, have nonpositive excess and require no classification.

## 4. Exact symbolic supports

Write an unordered normalized target coordinate as \([x,y]\), and work in
\(\mathbb Z[q]/(q^2+q+1)\).  First, the high-atom parameterization used below
does not silently omit a nonself positive-outside type.  If both negative
occurrences are aligned, write

\[
 c=\Psi_s(a)=\Psi_t(b),\qquad s,t\in\{q,q^2\}.
\]

For \(s=t\), injectivity forces \(a=b\).  For \(s\ne t\), the equality gives
\(b=\Psi_{s/t}(a)\), so the positive occurrence is reverse target, not
outside.  It becomes outside only when \([a,b]\) is the deleted
self-antipodal orbit.  The two ordered solutions have the same support
\({\cal P}\) in (13), so they are exactly the capped singleton already
classified, not a missing nonself high family.

The eight valid discrete HH assignments collapse to one four-coordinate
support

\[
\begin{split}
 {\cal H}=\{&[-1+2q,-1-4q],[-3-2q,3+4q],\\
             &[-3-4q,1+4q],[1-2q,3+2q]\}.
                                                        \tag{11}
\end{split}
\]

Call these four edges \(h_1,\ldots,h_4\) in the displayed order.  For each
\(j\), replace \(h_j\) by the other two edges of its \(\Phi\)-triangle.  The
four resulting five-coordinate supports are

\[
\begin{split}
{\cal F}_1=\{&[-1+2q,5+2q],[-1-4q,5+2q],h_2,h_3,h_4\},\\
{\cal F}_2=\{&h_1,[-3-2q,3-2q],h_3,h_4,[3+4q,3-2q]\},\\
{\cal F}_3=\{&h_1,h_2,[-3-4q,5],[1+4q,5],h_4\},\\
{\cal F}_4=\{&[-1,1-2q],[-1,3+2q],h_1,h_2,h_3\}.
                                                        \tag{12}
\end{split}
\]

An exhaustive affine calculation of the O2 stars gives exactly these same
four supports.  The calculation is all-orders, not a prime census.  Its 256
sign/endpoint branches split as follows:

- 176 rank-two systems, whose 48 valid assignments collapse to (12);
- 72 inconsistent determinant-zero systems;
- eight consistent rank-one systems, each with an identical aligned-edge
  collision.

The rank-two determinant norms have only prime factor \(3\).  Nonzero
augmented-minor witnesses for the 72 inconsistent branches have only prime
factors \(2,3\), so no inconsistent branch becomes consistent in the stated
range.

The R1 claim is also exhaustively checked rather than assumed from the
triangle picture.  For each of the four coordinates of \({\cal H}\), the
verifier tries both orders of both physical orbit representatives, both
choices of each \(q\)-multiplier, and all three possible locations of the
unique reverse occurrence: 192 branches, 64 per location.  Exactly eight
branches are consistent, all with the positive occurrence reverse; they
collapse to the four triangle flips (12).  Neither negative-occurrence
location has a solution.  The nonzero residuals in those discarded R1
equations have norm-prime factors only \(2,3,7\).

Solving the self-antipodal defect equation gives five atom types.  Two have
an internal aligned-coordinate duplication.  The other three all have the
same two-coordinate support

\[
 {\cal P}=\{[3,q^2-q],[3,q-q^2]\}.              \tag{13}
\]

The exact intersection table is

\[
 |{\cal H}\cap{\cal F}_j|=3,\qquad
 |{\cal F}_i\cap{\cal F}_j|=2\ (i\ne j),\qquad
 {\cal P}\cap{\cal H}={\cal P}\cap{\cal F}_j=\varnothing.
                                                        \tag{14}
\]

All coordinate-equality and antipodal support-collision reductions in
(11)--(14) have
possible exceptional primes only in

\[
 \{2,3,7,13\}.                                      \tag{15}
\]

The verifier separately audits every nonzero equality residual used to
discard or type-filter an O2, HH, cap, score-three, positive-outside-high,
or R1 branch.  Those possible relation-change primes are

\[
 \{2,3,5,7,13,19\}.                              \tag{15a}
\]

Together with the determinant and augmented-minor audit, (15)--(15a) prove
that the classification is uniform for every \(p\ge31\) in scope.
Reduction modulo such a prime can neither create a discarded branch nor
remove one of the exact intersections in (14); an additional collision
could only invalidate a proposed free support.

Because distinct pairing components have disjoint unpaired support, (14)
gives the exact weighted disjoint-packing maximum

\[
 \max\{\,2\# {\cal H}+\# {\cal F}+\#{\cal P}\,\}=3.       \tag{16}
\]

The maximum is attained by the disjoint pair \({\cal H},{\cal P}\).  One
cannot combine \({\cal H}\) with an \({\cal F}_j\), and at most one
\({\cal F}_j\) can occur.  Summing (6)--(10) therefore gives

\[
 b-2\Delta\le3,                                    \tag{17}
\]

where \(\Delta\) is the total deficit of the full pairing graph.

## 5. Pure equianharmonic threshold

Put

\[
 C={p-1\over3},\qquad E=r-1,\qquad
 L=C-E={r+5\over3}.                                \tag{18}
\]

The target score is \(3C-1\).  Exact equality and (3) give

\[
 3E+2b-\Delta=3C-1,\qquad
 \Delta=2b-3L+1.                                  \tag{19}
\]

If \(b\le2L-2\), then

\[
 b-2\Delta=6L-2-3b\ge4,                           \tag{20}
\]

contradicting (17).  Hence every exact equianharmonic fiber satisfies

\[
 b\ge2L-1={2r+7\over3}.                            \tag{21}
\]

This is a necessary threshold, not an existence assertion above it.

## 6. Exact threshold witnesses at 31 and 43

The threshold is attained in the two first in-scope characteristics.  The
canonical \(p=31,b=7,k=11\) replay, including its atom list, 29 target
coordinates, and all 105 odd channels, is in
`NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md`, Section 5.

There is also a self-contained \(p=43,b=9,k=13\) replay.  Here \(q=36\),
and the nine all-equal triples are

\[
\begin{split}
 &(2,31,33),(3,26,37),(5,6,12),(9,27,30),\
 &(10,20,36),(11,13,42),(14,17,35),(16,21,29),\
 &(32,38,39).
\end{split}                                             \tag{22}
\]

Writing \((\{a,b,c\};d)\) for a compact triple with distinguished label
\(d\), the nine compact atoms are

\[
\begin{split}
 &(\{0,24,39\};0),(\{1,25,36\};36),(\{2,36,40\};2),\
 &(\{3,7,42\};42),(\{3,18,41\};3),(\{4,9,19\};9),\
 &(\{4,9,35\};35),(\{8,24,39\};39),(\{15,20,23\};20).
\end{split}                                             \tag{23}
\]

Exact integer replay gives the constant tangent-conic target on all 41
nonantipodal coordinates, with \(\ell^1=41\), and zero in every one of the
210 odd channels of degrees \(3,5,\ldots,41\).  The even syndromes are

\[
 F_6=(37,19,8),\qquad F_8=(18,17,10,32),              \tag{24}
\]

so this is sharpness for the odd exact-fiber threshold only; it is not a
simultaneous zero-six/eight or global common-form witness.

## 7. Zero-odd centrality corollary

Assume a branch-C row has \(r-1\) all-equal atoms, \(b\le r\) compact atoms,
and every odd contraction through degree \(p-2\) vanishes.  Its associated
Cartesian dual word has support at most \(3(r+b-1)\le3h-6\).  The existing
Couvreur/Cayley--Bacharach reduction leaves the following possibilities.

- A support containing a maximal line reduces to one or two maximal lines;
  both coefficient patterns are already excluded for every \(b\le r\).
- The boundary cubic configuration is already excluded for every
  \(p\ge31\).
- A high-intersection conic is the tangent conic.  Its nonconstant dual
  coefficient is excluded by the integer \(\ell^1\) bound; the star and
  nonequianharmonic constant cases are excluded by parity and signed-score
  incidence.  The only remaining constant case is (1).

If \(p\equiv7\pmod {12}\) and \(3b\le2r+4\), then \(b<(2r+7)/3\), so (21)
excludes the last conic case.  If \(p\equiv11\pmod {12}\), a nontrivial
cube root \(q\) does not exist in \(\mathbb F_p\), so the same conclusion
holds throughout \(0\le b\le r\).  Thus the aggregate signed edge chain is
centrally symmetric in these ranges, under the zero-odd hypothesis only.

## 8. Replay

The exact certificate is implemented in
`src/e1_gmin_m4_equianharmonic_component_packing.py` and tested by
`tests/test_equianharmonic_component_packing.py`.  The public entry points
separate the three claims:

- `equianharmonic_component_packing_certificate()` checks the local table;
- `equianharmonic_exact_fiber_threshold(p,b)` checks (21);
- `p43_equianharmonic_threshold_witness_certificate()` replays (22)--(24);
- `p3_odd_radon_centrality_component_upgrade(p,b)` checks the dependency-
  gated centrality corollary.
