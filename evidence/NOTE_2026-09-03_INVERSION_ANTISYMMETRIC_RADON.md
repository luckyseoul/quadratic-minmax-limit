# Central inversion and the antisymmetric edge--Radon block

**Status:** exact all-prime rank, kernel, and integral-cokernel theorem;
exact hard-star target; and a constructive all-prime ternary lift of the
entire hard-star antisymmetric target.  This closes the antisymmetric box,
not the coupled symmetric half, and does not close residual (ii).

This note uses no finite configuration census.  It starts from the edge--Radon
map of Propositions 15.757--15.760 and the branch-C row theorem that, under
zero odd global forms, makes the certified opposite rows centrally symmetric.

## 1. The two involutions

Let \(V=\mathbf F_p^2\), where \(p=2h+1\) is odd, and write an edge uniquely
as

\[
 e=(a,[\delta])=\{a-\delta,a+\delta\},
 \qquad
 [\delta]\in\Delta=(V\setminus\{0\})/\{\pm1\}.
\]

Put \(d=p+1\).  Then

\[
 |\Delta|=dh={p^2-1\over2}.
\]

Central inversion acts on source edges by

\[
 J(a,[\delta])=(-a,[\delta]).
\]

Its fixed edges are exactly the \(dh\) antipodal edges with midpoint \(a=0\).
On the unsigned edge--Radon target, use midpoint/difference cells

\[
 P_L,\qquad K_L(\alpha,\beta),
 \qquad
 \alpha={s+t\over2},\quad
 \beta=\left({s-t\over2}\right)^2\in H.
\]

The corresponding target involution is

\[
 I(P_L)=P_L,\qquad
 I K_L(\alpha,\beta)=K_L(-\alpha,\beta).
\]

The edge--Radon map is equivariant:

\[
                         RJ=IR.                         \tag{1}
\]

The normalized residual map differs only by source-edge and direction-block
signs.  The edge sign is determined by the parallel direction and hence is
unchanged by \(e\mapsto-e\).  Those sign changes commute with \(J\) and \(I\),
so every statement below transports verbatim.

Define the integral antisymmetric lattices

\[
 E^-=\ker_{\mathbf Z}(1+J),\qquad
 {\cal A}^-=\ker_{\mathbf Z}(1+I)\cap{\cal A},
\]

where \({\cal A}\) is Proposition 15.760's ordinary compatible target
lattice.  Explicitly, a member of \({\cal A}^-\) has

\[
 P_L=0,\qquad K_L(0,\beta)=0,\qquad
 K_L(-\alpha,\beta)=-K_L(\alpha,\beta).
                                                               \tag{2}
\]

The ordinary equal-total and parallel-sum equations are automatic in (2).
Thus there are \(h^2\) free antisymmetric cells in every projective
direction.

## 2. Exact characteristic-zero rank and kernel

There are \(p^2dh\) source edges and \(dh\) fixed edges.  Hence

\[
 \operatorname {rank}E^-
 ={p^2dh-dh\over2}
 =(dh)^2=d^2h^2.                                           \tag{3}
\]

Proposition 15.760 says that \(R(E)\) has finite index in \({\cal A}\).
Over \(\mathbf Q\), apply the equivariant projector \((1-J)/2\) on the
source and \((1-I)/2\) on the target.  It follows immediately that

\[
 R(E^-\otimes\mathbf Q)={\cal A}^-\otimes\mathbf Q.
\]

Consequently

\[
 \boxed{
 \operatorname {rank}R^-=dh^2,\qquad
 \operatorname {rank}\ker R^-=dph^2.}                    \tag{4}
\]

There is a useful block proof of (4).  Extend scalars to a field containing
the \(p\)-th roots of unity and Fourier transform in the midpoint \(a\).
Antisymmetry pairs every nonzero frequency \(\xi\) with \(-\xi\); there are
\(dh\) such frequency pairs.  The projective direction \(L=[\xi]\) is
unique.  At this frequency the \(L\)-row maps the \(dh\) difference classes
to the \(h\) projected-square coordinates by

\[
 (c_\delta)_{\delta\in\Delta}
 \longmapsto
 \left(
   \sum_{[\delta]:\,L(\delta)^2=\beta}c_\delta
 \right)_{\beta\in H}.                                    \tag{5}
\]

Each of the \(h\) displayed fibers contains \(p\) difference classes and
the \(h\) classes parallel to \(L\) are invisible.  Therefore (5) has rank
\(h\) and kernel dimension \(ph\).  Summing over the \(dh\) frequency
pairs gives (4).

More generally, for any set \({\cal S}\) of projective directions,

\[
 \operatorname {rank}R^-_{\cal S}=|{\cal S}|h^2.          \tag{6}
\]

The frequency supports of distinct directions are disjoint.  In branch C
the hard and opposite direction sets each have \(d/2=h+1\) members.
Therefore centrality of every opposite row imposes exactly
\((h+1)h^2\) independent rational equations, and the remaining map onto
the hard antisymmetric rows is still surjective.  Its kernel is the full
kernel in (4).  There is no hidden characteristic-zero relation on the hard
target.

## 3. Exact integral cokernel: precisely the odd moments

The same Fourier-block argument works after reduction modulo two.  The
reduction of \(E^-\) consists of functions even in \(a\), zero at \(a=0\).
Over an algebraic closure of \(\mathbf F_2\), it has the basis

\[
 \psi_\xi(a)+\psi_{-\xi}(a),
 \qquad \xi\in(V^*\setminus\{0\})/\{\pm1\}.
\]

The reduction of \({\cal A}^-\) similarly has, in each projected-square
cell, the basis

\[
 \psi_j(\alpha)+\psi_{-j}(\alpha),
 \qquad j\in\mathbf F_p^*/\{\pm1\}.
\]

Because \(p\) is odd, the line-fiber factor \(p\) is one in
\(\mathbf F_2\), and the block is again exactly (5).  Thus

\[
                  R^-\pmod2\ \hbox{is surjective}.        \tag{7}
\]

Let

\[
 C^-={\cal A}^-/R(E^-),\qquad
 G={\cal A}/R(E).
\]

The natural map \(C^-\to G\) is injective.  Indeed, if
\(y\in{\cal A}^-\) and \(y=Rz\), then

\[
                   2y=R(z-Jz),\qquad z-Jz\in E^-.
\]

Its kernel is therefore killed by two, while (7) says that \(C^-\) has no
two-primary quotient and hence no two-torsion.

The image is exactly the \((-1)\)-eigenspace \(G^-\).  Containment is clear.
Conversely, if \(g\in G\) satisfies \(Ig=-g\), choose a representative
\(y\in{\cal A}\).  The class of \((1-I)y\in{\cal A}^-\) is \(2g\).
Proposition 15.760 gives \(pG=0\), and two is invertible modulo the odd prime
\(p\), so an integer multiple of \((1-I)y\) represents \(g\).

The moment basis of Proposition 15.759 diagonalizes \(I\).  A homogeneous
moment row of degree \(q\) changes by \((-1)^q\) under
\((s,t)\mapsto(-s,-t)\).  Hence \(G^-\) is exactly the span of the
odd-degree moment rows.  With \(q=2a+1\), their number is

\[
 \begin{split}
 S_-(p)
 &=\sum_{a=1}^{h-1}a\bigl(p-(2a+1)\bigr)\\
 &=2\sum_{a=1}^{h-1}a(h-a)
 ={h(h-1)(h+1)\over3}.
 \end{split}
\]

Therefore

\[
 \boxed{
 {\cal A}^-/R(E^-)
 \cong(\mathbf Z/p\mathbf Z)^{\,h(h-1)(h+1)/3}.}          \tag{8}
\]

Equivalently, an integral antisymmetric target has an integral
antisymmetric preimage if and only if it satisfies the odd-degree moment
congruences of Proposition 15.759.  There is no additional Smith, parity,
or even-moment obstruction in this block.

In the branch-C hypothesis of this note all global odd forms vanish.
Once the opposite rows are central, the remaining hard-supported
antisymmetric target therefore already has an unrestricted integral
antisymmetric lift by (8).

## 4. Exact hard-star target

For one fiber-label set \(\mathbf F_p\), let \(S_j\) be the unit star
containing all \(p-1\) unordered pairs incident with \(j\).  The
complement-literal term \(1-x_j\) contributes the signed coefficient chain
\(-S_j\).  Its inversion difference is

\[
                         A_j=S_{-j}-S_j.                  \tag{9}
\]

In midpoint/difference coordinates this is the explicit row

\[
 A_j(\alpha,\beta)
 =\mathbf1_{(\alpha+j)^2=\beta}
  -\mathbf1_{(\alpha-j)^2=\beta}.                         \tag{10}
\]

In particular \(A_0=0\).  Its signed vertex boundary is

\[
                 \partial A_j=(p-2)(e_{-j}-e_j).          \tag{11}
\]

For a compact atom

\[
 K(a,b;c)=\{a,b\}-\{a,c\}-\{b,c\},
\]

put \(K^-=K(a,b;c)-K(-a,-b;-c)\).  Then

\[
                 \partial K^-=2(e_{-c}-e_c).              \tag{12}
\]

Thus if a hard row \(L\) is

\[
 C_L=-S_{j_L}+\sum_{i=1}^{e_L}K(a_{Li},b_{Li};c_{Li}),
\]

its exact antisymmetric target, in the normalized convention, is

\[
 \Delta_L=\epsilon_L\left(
 A_{j_L}+\sum_{i=1}^{e_L}
 \bigl[
 K(a_{Li},b_{Li};c_{Li})
 -K(-a_{Li},-b_{Li};-c_{Li})
 \bigr]\right).                                          \tag{13}
\]

For a centrally symmetric opposite row, \(\Delta_L=0\).  Equations
(10)--(13) are the complete target forced by centrality; no source graph
has been assumed.

The separate hard-row compact-residual theorem now says that, on the full
balanced branch-C ray for \(p=4r+3\ge31\), zero odd global forms also force
the compact residual of every hard row to be central.  Since hard directions
have \(\epsilon_L=+1\), the antisymmetric target therefore reduces exactly to

\[
 \Delta_L=\begin{cases}
 A_{j_L},&L\text{ hard},\\
 0,&L\text{ opposite}.
 \end{cases}                                              \tag{13a}
\]

There are \((p+1)/2\) directions of each kind, and the centers \(j_L\) in
(13a) may be arbitrary.  No coherence such as \(j_L=L(v)\) is assumed.

Every unit star has zero contraction against every moment polynomial of
degree at most \(p-2\): for fixed \(j\), the sum over its other endpoint is
a field sum of a polynomial of degree at most \(p-2\), and the omitted
diagonal term is already zero.  Hence the hard-star part (9) contributes
nothing to the odd-moment cokernel (8).  Modulo two the compact boundaries
in (12) disappear, leaving the pair \(e_j+e_{-j}\), but (7) proves that this
also creates no binary obstruction.

## 5. The antisymmetric ternary box

Choose one representative \(e\) from each nonfixed source orbit
\(\{e,-e\}\), and let \(\Omega\) be the resulting set.  Define the
antisymmetric column map

\[
 D:\mathbf Z^\Omega\longrightarrow{\cal A}^-,
 \qquad D{\bf e}_{\{e,-e\}}=R{\bf e}_e-R{\bf e}_{-e},
                                                               \tag{14}
\]

with the normalized source signs included when using the residual map.
For a simple graph \(x\), set

\[
                         w_e=x_e-x_{-e}\in\{-1,0,1\}.     \tag{15}
\]

Then

\[
                         Dw=(1-I)y.                      \tag{16}
\]

Conversely, every \(w\in\{-1,0,1\}^\Omega\) is the inversion difference of
some binary source: use \((x_e,x_{-e})=(1,0)\), \((0,1)\), or either
\((0,0)\) and \((1,1)\), according as \(w_e=1,-1,0\).  Therefore

\[
 \boxed{
 \text{the antisymmetric target is Boolean-realizable}
 \Longleftrightarrow
 \exists\,w\in\{-1,0,1\}^{\Omega}:Dw=(1-I)y.}            \tag{17}
\]

This is necessary and sufficient for the antisymmetric half alone.  The
full graph must still solve the symmetric half with the coupled pair totals
\[
 s_e=x_e+x_{-e},\qquad
 (s_e,w_e)\in\{(0,0),(1,1),(1,-1),(2,0)\},
\]
and the fixed antipodal edges remain independent binary variables.

There is an exact nonnegative integer defect for (17):

\[
 \Gamma(w)={1\over2}\sum_{e\in\Omega}w_e^2(w_e^2-1),
 \qquad
 \Gamma_D(a)=\min_{Dw=a}\Gamma(w).                        \tag{18}
\]

For integral \(w\), every summand is nonnegative and vanishes exactly at
\(-1,0,1\).  Hence

\[
 \boxed{
 \Gamma_D((1-I)y)=0
 \Longleftrightarrow
 \text{the antisymmetric ternary fiber is nonempty}.}     \tag{19}
\]

If the fiber is empty, the first positive defect is at least six.  This is
a smaller exact nonlinear gate than the original Boolean fiber: opposite
centrality removes half the target blocks, inversion pairs the nonfixed
source edges, and the only antisymmetric box is ternary.

## 6. A direction-localized Mobius trade

The ternary gate (17) for the target (13a) always has an explicit solution.
The construction below works for every odd prime and does not use a finite
configuration census.

Fix one projective functional \(L\), a second functional \(M\) independent
of \(L\), and \(j\ne0\).  For

\[
             t\in T:=\mathbf F_p\setminus\{-1\},\qquad
             f(t)={t\over t+1},                             \tag{20}
\]

define two points by their \((L,M)\)-coordinates:

\[
             u_t=j(1,f(t)),\qquad v_t=j(t,t),               \tag{21}
\]

and put

\[
 E_{L,M,j}=\{\{u_t,v_t\}:t\in T\},\qquad
 z_{L,M,j}={\bf1}_{-E_{L,M,j}}-{\bf1}_{E_{L,M,j}}.          \tag{22}
\]

These are genuine distinct edges.  If \(u_t=v_t\), the first coordinate
forces \(t=1\), while the second would say \(1/2=1\).  The projected
\(L\)-cell \(\{j,jt\}\) determines \(t\), so different parameters give
different edges.  An edge could be fixed by inversion, or an edge of
\(E\) could lie in \(-E\), only if its \(L\)-cell were equal to its
negative; because \(j\ne0\), this forces \(t=-1\), which was excluded.
Thus (22) has exactly \(p-1\) nonfixed inversion orbits and is ternary.

In direction \(L\), the parameter \(t=1\) gives one parallel edge in each
half and cancels from \(P_L\).  The remaining parameters give

\[
 R_L z_{L,M,j}
 =\bigl(S_{-j}-\{j,-j\}\bigr)
  -\bigl(S_j-\{j,-j\}\bigr)=A_j.                           \tag{23}
\]

Every other direction cancels exactly.  After scaling its functional,
such a direction is either \(N=M\) or \(N=L+mM\) with \(m\ne0\).  Divide
the two projected endpoint labels by \(j\).  In the latter case they are

\[
 A_t=1+m f(t)={(1+m)t+1\over t+1},\qquad
 B_t=(1+m)t.                                               \tag{24}
\]

If \(m\ne-1\), put \(c=1+m\) and

\[
 \phi_m(t)=-{ct+1\over c(t+1)}.                            \tag{25}
\]

Here \(c\ne0,1\).  Consequently \(\phi_m(t)\ne-1\), and direct
substitution gives

\[
                 (A_{\phi_m(t)},B_{\phi_m(t)})
                    =(-B_t,-A_t).                          \tag{26}
\]

The injectivity of \(B_t=ct\) also gives \(\phi_m^2(t)=t\).  If a
parameter is fixed, its projected unordered cell is itself fixed by
central negation, so it contributes zero to (22); otherwise the two
parameters cancel as a pair.  For the two exceptional projective rows the
same calculation is

\[
 \begin{array}{c|c|c}
 N&\phi(t)&\text{projected ordered pair after }\phi\\ \hline
 L-M&-t-2&(-A_t,-B_t)\\
 M&-t/(t+1)&(-B_t,-A_t).
 \end{array}                                               \tag{27}
\]

Both maps preserve \(T\) and square to the identity.  Parallel projected
edges cancel in the common \(P_N\) coordinate as well.  Hence

\[
 \boxed{R_Lz_{L,M,j}=A_j,\qquad R_Nz_{L,M,j}=0\quad(N\ne L).}
                                                               \tag{28}
\]

The zero center needs no trade because \(A_0=0\).

### Disjoint simultaneous trades

It remains to keep the sum of one trade per hard direction inside the
ternary box.  This can always be done by choosing the auxiliary
functionals greedily.  For fixed \(L\), exactly \(p^2-p=p(p-1)\)
functionals \(M\) are independent of \(L\).

Fix one previously used inversion orbit and ask for how many \(M\) it can
occur in \(E_{L,M,j}\), accounting at once for intersection with either
the old \(E\) or the old \(-E\).  Orient the old orbit so that its candidate
endpoint has \(L\)-label \(j\).  This loses no second orientation: a
candidate orbit whose endpoint labels were \(\{j,-j\}\) would have
\(t=-1\), and that parameter is not in \(T\).

* At \(t=0\), the other endpoint in (21) is the origin.  One old origin
  orbit forbids at most the \(p-1\) nonzero functionals having the required
  kernel.
* At \(t\ne0,1\), the two endpoints in (21) are linearly independent and
  the endpoint of label \(j\) is unique.  Their two prescribed \(M\)-values
  determine at most one functional.
* At \(t=1\), both endpoints have \(L\)-label \(j\), so exchanging them can
  determine at most two functionals.

A previous trade has one origin orbit and \(p-2\) other orbits.  Moreover,
in the coordinates of (21),

\[
 \det(u_t-v_t,u_s-v_s)
 ={j^2(t-s)(t+s)\over(t+1)(s+1)}.                           \tag{29}
\]

Thus it has at most two edges in any spatial direction.  In particular at
most two of its non-origin orbits can incur the extra ordering in the
third bullet.  One previous trade therefore forbids at most

\[
                   (p-1)+(p-2)+2=2p-1                     \tag{30}
\]

auxiliary functionals.  Before the last of at most \((p+1)/2\) hard trades,
there are at most \((p-1)/2\) prior trades.  Even the union bound leaves

\[
 p(p-1)-{p-1\over2}(2p-1)={p-1\over2}>0                  \tag{31}
\]

choices.  This proves that all nonzero-orbit supports can be made pairwise
disjoint, for arbitrary centers \(j_L\).

Finally return to the normalized residual map.  It is
\(\widetilde R_N(w)=\epsilon_NR_N(\tau w)\), where
\(\tau_{-e}=\tau_e\in\{\pm1\}\).  Put \(w_e=\tau_ez_e\).  This preserves
ternarity and antisymmetry.  On a hard direction \(\epsilon_L=+1\), so
(28) is exactly (13a); all other rows remain zero.  If \(q\) hard centers
are nonzero, the resulting source uses exactly \(q(p-1)\) inversion orbits.
It is the inversion difference of a simple graph: on every used orbit,
select the unique edge on which \(w_e=1\).

Therefore

\[
 \boxed{\Gamma_D((1-I)y)=0\text{ for every branch-C target in (13a).}} \tag{32}
\]

## 7. Exact consequence and remaining obstruction

Central opposite rows plus hard-star antisymmetry do **not** yield a
contradiction.  The linear diagnostics already showed why:

1. the hard target is arbitrary over \(\mathbf Q\) after the opposite blocks
   vanish, by (6);
2. reduction modulo two is surjective, by (7);
3. zero odd global forms are exactly the remaining integral compatibility,
   by (8); and
4. the unit-star contribution itself is invisible to every odd moment.

The Mobius trade proves that the nonlinear value of \(\Gamma_D\) is zero for
the exact hard-star target, even when all centers are unrelated.  What
remains is only the coupled symmetric-half equation.  On every orbit used by
(22), its pair total is forced to be \(s_e=1\); on every unused nonfixed
orbit one may choose \(s_e=0\) or \(2\), and the fixed antipodal edges are
independent binary variables.  The auxiliary choices guaranteed by (31)
give genuine freedom, but no argument here matches those symmetric totals
to the prescribed full rows.

For orientation, at \(p=31\), \(h=15\), \(d=32\), the formulas give source
antisymmetric rank \(230400\), target rank \(7200\), full kernel rank
\(223200\), and rank \(3600\) for either half of the projective directions.
The construction uses at most \(16\cdot30=480\) nonfixed edge orbits for all
hard stars.  These numbers are evaluations of symbolic formulas, not a
finite search.

Residual (ii), E1, \(L=1/2\), and the original MathOverflow problem remain
OPEN.
