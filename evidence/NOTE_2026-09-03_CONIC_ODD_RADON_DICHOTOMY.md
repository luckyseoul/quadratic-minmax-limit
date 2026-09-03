# Irreducible-conic odd-Radon dichotomy

This note treats only a low-weight odd-Radon word whose support contains the
minimal irreducible-conic configuration.  It neither assumes nor proves that
the degree-six/eight global forms vanish.

Let

\[
 p=4r+3\ge31,\qquad h=(p-1)/2=2r+1,\qquad m=h-2=2r-1,
\]

and let \(H\subset\mathbf F_p^*\) be the nonzero squares and
\(\Omega=H\times H\).  A row with \(r-1\) positive all-equal atoms and
\(b\le r\) compact atoms has signed-occurrence budget

\[
 N=3(r+b-1)\le 6r-3=3m.
\]

Its odd-zero condition gives a word \(W\) on \(\Omega\), orthogonal to all
plane polynomials of total degree at most \(m\), with

\[
 W(U,D)=n_E D\sigma,\qquad U=\sigma^2,\qquad D=(s-t)^2,qquad
 \sum_E|n_E|\le N.                                      \tag{1}
\]

## 1. Conic peeling

Suppose the support of \(W\) contains at least \(2m+2=p-3\) points on an
irreducible conic \(Q=0\).  Put \(R=\operatorname{supp}(W)\setminus Q\).
For every polynomial \(G\) of degree at most \(m-2\), orthogonality applied
to \(QG\) says

\[
 \sum_{P\in R} W(P)Q(P)G(P)=0.                           \tag{2}
\]

If \(R\ne\varnothing\), the nonzero weights \(W(P)Q(P)\) therefore form a
word orthogonal through degree \(m-2\).  Such a word has support at least
\(m\): if it had at most \(m-1\) points, a product of at most \(m-2\) lines
could vanish at all but any prescribed support point.  On the other hand,

\[
 |R|\le3m-(2m+2)=m-2,
\]

a contradiction.  Thus the **entire** word is supported on \(Q\).  This is
the step needed to pass from “contains a conic circuit” to “conic-supported”;
no external low-weight classification is being silently used here.

## 2. Every high-intersection irreducible conic is triangle-tangent

The conic has at least \(p-3\) points in \(\Omega\).  It is consequently a
smooth rational conic; a geometrically reducible but \(\mathbf F_p\)-irreducible
conic has at most its conjugate lines' intersection as an \(\mathbf F_p\)
point.  Parameterize it by \(\mathbf P^1\), writing its three projective
coordinates as binary quadratics \(X,Y,Z\), so \(U=X/Z\), \(D=Y/Z\).

Let \(\mathcal E\) be the roots of \(XYZ\); \(|\mathcal E|\le6\).  Away
from \(\mathcal E\), expand the two square indicators.  If none of the
geometric square classes of \(U,D,UD\) is constant, each corresponding
quadratic cover of \(\mathbf P^1\) has at most four branch points and genus
at most one.  The three character sums are bounded by \(2\sqrt p\), and
removing \(\mathcal E\) costs at most six in each sum.  Hence

\[
 4|Q\cap\Omega|\le p+19+6\sqrt p<4(p-3).                 \tag{3}
\]

If exactly one of the three square classes is constant, the same calculation
gives

\[
 4|Q\cap\Omega|\le2p+14+4\sqrt p<4(p-3).                \tag{4}
\]

Both strict inequalities hold at \(p=31\) and their left-minus-right
margins increase thereafter.  Since the three classes satisfy
\([UD]=[U]+[D]\), the only remaining case is that all three are constant.
A nonsquare constant for either \(U\) or \(D\) would give no point of
\(\Omega\), so in fact \(U\) and \(D\) are squares in
\(\mathbf F_p(\mathbf P^1)\).

Because \(X,Z\) are quadratics and the nonconstant ratio \(X/Z\) is a
square, cancellation cannot leave a ratio of two distinct linear forms;
thus \(X=uL_1^2\), \(Z=L_0^2\).  Similarly \(Y=dL_2^2\).  The three linear
forms are pairwise nonproportional.  Sending their roots to
\(0,1,\infty\) yields the exact normal form

\[
 U=u z^2,qquad D=d(z-1)^2,qquad u,d\in H.              \tag{5}
\]

Consequently

\[
 Q\cap\Omega=\{z\in\mathbf F_p:z\ne0,1\},\qquad
 |Q\cap\Omega|=p-2.                                     \tag{6}
\]

These are precisely the conics tangent to the three sides of the coordinate
triangle.

## 3. The dual weights leave only one constant branch

The restriction of plane forms of degree \(m\) to a smooth conic is the
complete space \(H^0(\mathbf P^1,\mathcal O(2m))\).  Here
\(2m=p-5\).  For \(T=\mathbf F_p\setminus\{0,1\}\),

\[
 P_T(x)=\frac{x^p-x}{x(x-1)},\qquad P_T'(z)=-\frac1{z(z-1)}.
\]

The two-dimensional dual Reed--Solomon space is therefore

\[
 W(z)=c\,z(z-1)(Az+B).                                  \tag{7}
\]

Choose the edge-orbit orientation with \(\sigma=a z\), \(a^2=u\).  Dividing
(7) by \(D\sigma\) gives

\[
 n(z)=\alpha+\frac{\beta}{z-1}.                         \tag{8}
\]

If \(\beta\ne0\), the values in (8), for \(z\in T\), are every field
element except two.  The total least-absolute-residue mass is at least

\[
 h(h+1)-2h=h(h-1)>3m\ge N,                              \tag{9}
\]

impossible.  Hence \(\beta=0\).  If the remaining nonzero constant residue
had least absolute value at least two, its mass would be
\(2(p-2)>3m\); one nonminimal integer lift already gives the same
contradiction.  Thus every actual orbit difference is the same integer
\(+1\), or all are \(-1\), on all \(p-2\) conic points.  Numerically this
branch can occur only when

\[
 3(r+b-1)\ge p-2\quad\Longleftrightarrow\quad3b\ge r+4. \tag{10}
\]

Equations (5), (8), and (10) are the exact conic-supported reduction.

## 4. Physical graph and the star/nonstar split

Scale \(u=1\), write \(d=k^2\), and choose the physical edge

\[
 E_z(k)=\left\{
 \frac{(1+k)z-k}{2},\frac{(1-k)z+k}{2}
 \right\},\qquad z\ne0,1.                              \tag{11}
\]

For \(k=\pm1\), one endpoint is constant.  After projecting label vertices
by \(x\sim-x\), the graph in (11) has the two odd vertices
\([0]\) and \([1/2]\).  Every all-equal or compact atom becomes a triangle
boundary modulo two, so this star branch is impossible.

For \(k\ne\pm1\), put

\[
 q=\frac{1-k}{1+k}.
\]

After translating by \(1/2\), the full edge relation is \(t\mapsto qt\).
The omitted \(z=1\) edge is the loop at \(1/2\), and the omitted \(z=0\)
edge joins the antipodal pair \(-k/2,k/2\).  Adding just that invisible
self-antipodal edge makes (11) Eulerian, including after quotienting by
\(x\sim-x\).  Thus neither the integer mass bound nor the quotient-cycle
parity invariant excludes the nonstar constant branch.

The signed atom incidence nevertheless forces the equianharmonic case.
Normalize the constant word to \(+1\), which loses no generality after
globally negating all labels, and put \(X=2x\).  An aligned positive physical
edge lies in the graph of

\[
 \Phi(X)=qX+1-q,
\]

whereas a physical edge carrying coefficient \(-1\) is aligned precisely
when it lies in the graph of

\[
 \Psi(X)=qX+q-1=(-1)\circ\Phi\circ(-1)(X).
\]

Give each atom its alignment score: aligned occurrences count \(+1\),
reverse occurrences count \(-1\), and all other or self-antipodal
occurrences count zero.  If \(q^3\ne1\), an all-equal atom has score at most
two.  Indeed, score three would be a triangle in the disjoint-cycle graph of
\(\Phi\), forcing a three-cycle and hence \(q^3=1\).

A compact atom of score three has positive edge \(ab\) in the \(\Phi\)
graph and both negative edges \(ac,bc\) in the \(\Psi\) graph.  Thus, for
some \(i,j,\ell\in\{1,-1\}\),

\[
 b=\Phi^i(a),\qquad c=\Psi^j(a)=\Psi^\ell(b).           \tag{11a}
\]

Solving these eight linear cases is decisive.  Four give the degenerate
solution \(a=b=1\); two pairs merely exchange \(a,b\).  With
\(S=q^2+q+1\ne0\), the only two possible unordered compact atoms are

\[
\begin{split}
 &\left\{-\frac{q+2}{q},-(2q+1);-3\right\},\\
 &\left\{\frac{q^2-q-1}{S},
          -\frac{q^2+q-1}{S};
          -\frac{q^2-q+1}{S}\right\},
\end{split}                                             \tag{11b}
\]

where the entry after the semicolon is distinguished.  A missing target
edge or an accidental collision can only delete a candidate.

This constant-size list remains sufficient even when atoms may repeat.
Let \(T=(r-1)+b\le2r-1\) be the number of atoms, with multiplicity, and let
\(n_3\) count score-three compact occurrences.  Write \(A,R,O\) for the
total aligned, reverse, and other physical edge occurrences.  Exact equality
to the \(p-2=4r+1\) unit target coordinates gives

\[
 A-R=4r+1,\qquad 2R+O=3T-(4r+1).                       \tag{11c}
\]

Every non-score-three atom satisfies
\(2R_i+O_i\ge1+R_i\).  The \(n_3\) exceptional occurrences contribute
\(3n_3\) aligned occurrences on at most six target coordinates.  Since the
final coefficient on each coordinate is one, cancellation forces

\[
 R\ge\max(0,3n_3-6).
\]

Consequently \(2R+O\ge T-2\) when \(n_3\le2\), while for \(n_3\ge3\),

\[
 2R+O\ge T-n_3+3n_3-6\ge T.
\]

In all cases (11c) would give

\[
 4r+1=3T-(2R+O)\le2T+2\le4r,
\]

a contradiction.  Therefore every nonstar constant tangent-conic atom
fiber with \(b\le r\) necessarily has \(q^3=1\).

Since \(q\ne1\), this survivor is characterized exactly by

\[
 q^3=1, q\ne1
 \quad\Longleftrightarrow\quad
 q^2+q+1=0
 \quad\Longleftrightarrow\quad
 k^2=-3.                                                 \tag{12}
\]

In particular \(p\equiv7\pmod {12}\).  The full nonloop graph is then a
disjoint union of \((p-1)/3\) triangles.

## 5. Exact p=31 counterexample

For \(p=31\), take \(r=b=7\) and \(k=11\).  Here \(k^2=-3\),
\(q=25\), and \(q\) has order three.  Six positive all-equal atoms are

\[
 (1,3,13),(2,7,8),(4,18,26),(5,20,23),(6,14,28),(24,25,30).
\]

Seven compact atoms, displayed as \((\{a,b,c\};c_*)\) with distinguished
label \(c_*\), are

\[
\begin{split}
 &(\{0,2,12\};2),\ (\{0,12,19\};12),\ (\{4,19,22\};4),\\
 &(\{9,12,16\};16),\ (\{10,14,17\};14),\\
 &(\{11,19,22\};19),\ (\{12,16,20\};20).
\end{split}                                             \tag{13}
\]

Direct integer replay of all nonantipodal edge orbits gives exactly \(+1\)
on the 29 edges (11) and zero elsewhere.  Hence all 105 odd contractions
of degrees \(3,5,\ldots,29\) vanish, but the chain is not central.  The six
all-equal triples in (13) are six of the order-three \(q\)-cycle triangles.
The signed-occurrence audit is

\[
A=31,\qquad R=2,\qquad O=6,\qquad A-R=29,\qquad 2R+O=10,
\]

where \(A,R,O\) count aligned, reverse, and off-conic/self-antipodal signed
occurrences.  Its even syndromes begin

\[
 F_6=(11,19,10),\qquad F_8=(12,11,23,6),                \tag{14}
\]

so this particular row is **not** a zero-six/eight witness.  It is instead
an exact counterexample to any attempted extension of odd-zero centrality
through the full balanced range using only low weight, integer \(\ell^1\),
and quotient parity.

The deterministic certificate in
`src/e1_gmin_m4_conic_odd_radon.py`, replayed by
`tests/test_conic_odd_radon.py`, checks (13), every edge orbit, all 105 odd
channels, and the displayed degree-six/eight syndromes.  A later exact MITM
certificate rules out simultaneous zero degree-six/eight syndromes in this
specific (p=31,b=7,k=11) constant-conic fibre; see
`NOTE_2026-09-03_P31_EQUIANHARMONIC_ZERO68_MITM.md`.  That is a local
zero-form exclusion, not a residual-(ii) closure.

The remaining structural question is not “exclude every conic.”  It is:
classify the equianharmonic constant tangent-conic atom fibers, and decide
whether their nonzero degree-six/eight syndromes can be the evaluations of
common global forms across all directions.  No general parametric
equianharmonic construction or optimal threshold is asserted here.

## 6. A real, but limited, projective interpolation obstruction

There is one exact cross-direction consequence for the **unsigned**
scaled-copy ansatz (equivalently, for a model whose normalization sign is
constant on all directions).  It rules out that most obvious algebraic
globalization of (13), while not by itself addressing the signed Paley
normalization or direction-dependent choices from different local fibers.

Suppose one fixed admissible row has nonzero degree-six and degree-eight
syndromes \(v_6,v_8\), and in every projective direction \(L\) one uses a
nonzero scalar copy of that same row, with label scale \(\rho_L\in
\mathbf F_p^*\).  Homogeneity gives

\[
 M_6(L)=\rho_L^6v_6,\qquad M_8(L)=\rho_L^8v_8.           \tag{15}
\]

Choose one nonzero coordinate of each vector.  If these unsigned rows came
from common binary forms, there would be binary forms \(A\) and \(B\), of
degrees six and eight, satisfying

\[
 A(L)=\rho_L^6,qquad B(L)=\rho_L^8
\]

after fixed nonzero rescaling.  Hence \(A^4-B^3\) vanishes at all \(p+1\)
projective points.  Its degree is 24, so for \(p\ge31\) it is identically
zero.  Unique factorization and \(\gcd(3,4)=1\) then give

\[
 A=cQ^3,qquad B=dQ^4                                  \tag{16}
\]

for a binary quadratic \(Q\) and compatible constants \(c,d\).
Equation (15) makes \(Q(L)\) nonzero and of one fixed quadratic character
for every \(L\in\mathbf P^1(\mathbf F_p)\).  This is impossible.  If \(Q\)
has zero discriminant it has a projective root; if it has nonzero
discriminant, the standard quadratic character sum over the affine line,
including its value at infinity, is zero rather than \(\pm(p+1)\).

Thus the p=31 witness cannot simply be scaled independently by one nonzero
field scalar in every direction **under this unsigned common-form ansatz**.
At least one direction must use a different local fiber (or a zero-syndrome
row) in that ansatz.  This is not yet an exceptional-row theorem for the
actual signed Paley system, and more general mixtures of local conic
witnesses remain open.  Nor is it a contradiction to the seven-channel
dominance calculation: that calculation allows arbitrary extension-valued
local labels and does not impose the common finite-field projective
evaluation in (15).

For the actual Paley split, however, the same calculation gives a
constructive half-system rather than an obstruction.  If a projective
direction \(L=(R:S)\) has kernel generator \(d=S-R\alpha\in\mathbf F_{p^2}\),
then the repository normalization uses

\[
 \epsilon_L=\chi_{p^2}(d)=\chi_p(N(d)).                 \tag{17}
\]

Thus \(\Delta(L)=N(S-R\alpha)\) is an anisotropic binary quadratic whose
character is the direction type.  Choose a constant \(c\) so that
\(Q=c\Delta\) is a square on every opposite-type direction.  On every such
direction carrying the same conic atom profile, choose
\(\rho_L^2=Q(L)\).  Then the scaled copies satisfy the signed common-form
requirements with the explicit forms

\[
 F_{6,j}=\epsilon_{\operatorname{opp}}(v_6)_j Q^3,qquad
 F_{8,j}=\epsilon_{\operatorname{opp}}(v_8)_j Q^4.                 \tag{18}
\]

This is exact because \(\epsilon_L\) is constant on that projective half
and the moments are homogeneous.  It does **not** construct the other
Paley-type rows.  At the upper endpoint of the balanced branch-C ray there
are three opposite rows with \(b=r-1\) and the remaining \(2r-1\) opposite
rows have \(b=r\); for \(p=31\), that is three \(b=6\) rows and thirteen
\(b=7\) rows.  Thirteen matching evaluations already determine both a
binary sextic and a binary octic, so the exceptional opposite rows and all
hard rows would be forced to match (18).  The exact global gate is thereby
reduced, not solved.
