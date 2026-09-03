# Affine-block projections of Möbius-contained halved dual words

Date: 2026-09-03

Status: proved a structural cover theorem and the exact per-half resultant.
Any nonzero halved-code dual support contained in the actual branch-C
Möbius support either contains an already-excluded fixed-transverse
rectangle or forces an affine midpoint block over one difference block to
be covered by every half. Distinct target directions alone do not rule out
the latter: an explicit construction saturates the cover with freely
chosen centers. Prescribed centers, mutual ternarity, full dual-support
containment, punctured surjectivity, and residual (ii) remain open.

## 1. Cells and block-parity projections

Let \(p=2h+1\) be an odd prime and

\[
 \Delta=(\mathbf F_p^2\setminus\{0\})/\{\pm1\}.
\]

For a projective functional \(K\), put

\[
 A_K=\{[a]:K(a)=0\},\qquad
 B_{K,\alpha}=\{[a]:K(a)^2=\alpha\},
\]

where \(\alpha\) ranges over the \(h\) nonzero squares. These \(h+1\)
sets partition \(\Delta\), with

\[
 |A_K|=h,\qquad |B_{K,\alpha}|=p.                         \tag{1}
\]

The nonorigin blocks \(B_{K,\alpha}\), over all \(K,\alpha\), form the
paired-affine-line design from the fixed-edge-elimination note. Its square
incidence matrix \(M\) satisfies

\[
                         M^{\mathsf T}M=I
                         \quad\text{over }\mathbf F_2.    \tag{2}
\]

Let \(w\in\operatorname{Row}(D)\) be a word of the full halved row code.
For a difference block \(C=B_{K,\beta}\), define its midpoint parity
projection

\[
 r_C([a])=\sum_{[\delta]\in C}w([a],[\delta])\pmod2.      \tag{3}
\]

## 2. Exact projection theorem

For every \(C=B_{K,\beta}\), the word \(r_C\) is a disjoint union of
cells from

\[
                   A_K,\quad B_{K,\alpha}\quad(\alpha\ne0). \tag{4}
\]

Thus, for some \(\epsilon_C\in\{0,1\}\) and
\(q_C\in\{0,\ldots,h\}\),

\[
                         |r_C|=\epsilon_Ch+q_Cp.          \tag{5}
\]

It is enough to check the raw halved \(P,F,C\) row generators, because the
change to \((C,\Phi)\) is an invertible row operation.

- A paired row \(C_{L,\alpha,\gamma}\) projects to
  \(B_{L,\alpha}\) when \(C=B_{L,\gamma}\), and to zero otherwise.
- A fixed-transverse row \(F_{L,\gamma}\) projects to \(A_L\) when
  \(C=B_{L,\gamma}\), and to zero otherwise.
- A parallel row \(P_L\) projects to zero if \(L=K\), and to all of
  \(\Delta\) if \(L\ne K\). Indeed, \(A_L\) misses \(C\) in the first
  case and meets it in one antipodal class in the second.

The first two assertions use (2): two distinct affine blocks have even
intersection, while a block has odd self-intersection \(p\). The last
generator also lies in the span (4), because those cells partition
\(\Delta\). This proves (4)--(5).

The projections determine \(w\). For each fixed \([a]\), the vector
\((r_C([a]))_C\) is the difference fibre \(w([a],\cdot)\) multiplied by
the invertible matrix \(M\).

If \(q_C=0\) for every \(C\), then \(w\) is exactly an XOR of
fixed-transverse rectangles

\[
                         F_C=A_K\times C.                 \tag{6}
\]

Indeed, \(F_C\) has \(r_C=A_K\) and all its other block projections zero.
Subtract the indicated rectangles from \(w\); every projection then
vanishes, so invertibility of \(M\) kills the remainder. These rectangles
are pairwise disjoint: blocks of one direction are disjoint in the
difference factor, while \(A_K\cap A_L=\varnothing\) for \(K\ne L\)
inside \(\Delta\).

## 3. What one Möbius half can cover

Now specialize to \(p\equiv3\pmod4\). For a localized half
\(H(L,M,j)\), put \(z=t+1\in\mathbf F_p^*\). Its midpoint satisfies

\[
                         L(a)(L-M)(a)=j^2/4.              \tag{7}
\]

This nonzero conic meets an origin line \(A_K\) in at most one antipodal
class. It meets a paired nonorigin affine block \(B_{K,\alpha}\) in at
most two antipodal classes: intersect one of its affine lines with the
conic, then use central symmetry for the negative mate. The conic has no
affine-line component because the right side of (7) is nonzero.

If \(K=xL+yM\), direct substitution gives

\[
 K(\delta_z)={j\over2}
 \left((x+y)(z-2)+{y\over z}\right).                      \tag{8}
\]

For \(C=B_{K,\beta}\), choose \(s^2=\beta\). Each equation
\(K(\delta_z)=\pm s\), after multiplication by \(z\), is a nonzero
quadratic. Hence

\[
 \begin{aligned}
 |H\cap(\Delta\times C)|&\le4,\\
 |\operatorname{pr}_a(H)\cap A_K|&\le1,\\
 |\operatorname{pr}_a(H)\cap B_{K,\alpha}|&\le2.
 \end{aligned}                                           \tag{9}
\]

A half has exactly two columns above each of its midpoint classes, so

\[
                  |H\cap(A_K\times C)|\le2.              \tag{10}
\]

Let \(U\) be the actual support after all orientations and cancellations,
and assume only

\[
                         U\subseteq\bigcup_{i=1}^mH_i.    \tag{11}
\]

Neither ternarity nor an orientation convention is used in the cover
argument. They matter when constructing \(U\), but cancellations can only
shrink the right side of (11). Equations (9)--(10) imply

\[
 |U\cap(\Delta\times C)|\le4m,\qquad
 |U\cap(A_K\times C)|\le2m.                              \tag{12}
\]

## 4. The forced all-halves cover

Take the branch-C number of halves

\[
                         m={p+1\over2}=h+1.               \tag{13}
\]

Their raw number of columns is \(m(p-1)=|\Delta|\), so every dual word
contained in \(U\) automatically has weight at most \(|\Delta|\).

Suppose

\[
 0\ne w\in\operatorname{Row}(D),\qquad
 \operatorname{supp}(w)\subseteq U.
\]

If every \(q_C\) in (5) vanished, Section 2 would make \(w\) a disjoint
sum of rectangles (6). Any nonzero summand would have to lie in \(U\), but

\[
                  2m=p+1<ph=|A_K\times C|                \tag{14}
\]

for the branch primes \(p\ge31\), contradicting (12). Therefore some
\(r_C\) contains a full affine midpoint block \(B_{K,\alpha}\).

Moreover, (5) and (12) give

\[
                         \epsilon_Ch+q_Cp\le4m=2p+2.      \tag{15}
\]

For \(p\ge7\), this leaves only

\[
 A_K,\quad B_{K,\alpha},\quad
 A_K\cup B_{K,\alpha},\quad
 B_{K,\alpha_1}\cup B_{K,\alpha_2}.                      \tag{16}
\]

In particular \(q_C\le2\), and \(q_C=2\) forces
\(\epsilon_C=0\). The first alternative cannot be the only nonzero type
because of (14).

To cover the \(p=2h+1\) midpoint classes of an affine block in (16), while
each half supplies at most two by (9), requires all \(h+1\) halves. At
least \(h\) of them must supply two distinct classes. If two affine blocks
occur, every half must meet both, and the required \(2p\) incidences are
within two of the absolute capacity \(4m=2p+2\).

Thus any contained nonzero dual word forces an all-halves common affine
block cover over one difference block. This is stronger than the
single-rectangle exclusion, but it is not yet a contradiction.

## 5. Exact per-half resultant

Fix one half \(H(L,M,j)\), write

\[
                         K=xL+yM,\qquad A=x+y,             \tag{17}
\]

and choose \(r^2=\alpha\), \(s^2=\beta\). Equations (7)--(8) give

\[
 \begin{aligned}
 K(a_z)&={j\over2}\left(Az-{y\over z}\right),\\
 K(\delta_z)&={j\over2}\left(A(z-2)+{y\over z}\right).
 \end{aligned}                                           \tag{18}
\]

Put \(K(a_z)=\epsilon r\) and
\(K(\delta_z)=\eta s\), where \(\epsilon,\eta\in\{\pm1\}\).
Eliminating \(z\) and \(\epsilon\) yields

\[
 \boxed{\alpha-\beta-Axj^2=2\eta Ajs},\qquad
 \boxed{(\alpha-\beta-Axj^2)^2=4A^2j^2\beta}.             \tag{19}
\]

The second equation is the sign-free resultant. When the signed equation
holds, the two parameter candidates are

\[
                 z_\epsilon={Aj+\eta s+\epsilon r\over Aj}. \tag{20}
\]

They give two valid distinct midpoint classes precisely on the
nondegenerate locus

\[
                         Ay(Aj+\eta s)\ne0.               \tag{21}
\]

If \(A=0\) or \(y=0\), the midpoint conic meets the block in at most one
class. If \(Aj+\eta s=0\), the two candidates in (20) are antipodal and
again give one class.

There is an intrinsic form. Assume \(K,L\) are independent and represent
two distinct classes of \(B_{K,\alpha}\) by \(a_1,a_2\) with
\(K(a_1)=K(a_2)=r\). Put \(l_k=L(a_k)\). Some auxiliary \(M\) makes the
half supply both classes over \(C=B_{K,\beta}\) if and only if, after
possibly ordering the classes, for one \(\eta\in\{\pm1\}\),

\[
 l_1l_2(l_1+l_2)\ne0,\qquad
 (r-\eta s)l_1-(r+\eta s)l_2=rj.                         \tag{22}
\]

The auxiliary is then unique:

\[
                   M(a_k)=l_k-{j^2\over4l_k}\quad(k=1,2). \tag{23}
\]

The points are independent, so (23) defines a functional. It is not
proportional to \(L\): otherwise (7) would force \(l_1^2=l_2^2\), contrary
to distinctness and the nonzero sum in (22). The two relevant parameters
are \(z_1=2l_1/j\) and \(z_2=-2l_2/j\). If \(K=xL+yM\), subtracting the
two instances of (23) gives \(A=r/(l_1+l_2)\), and direct substitution
shows

\[
                   K(\delta_{z_1})=K(\delta_{z_2})=\eta s.
\]

This proves necessity and sufficiency of (22).

## 6. Distinct directions do not obstruct saturation

Condition (22) gives an explicit counterexample to any obstruction using
only distinctness of the target directions. The construction works for
every \(p>7\), hence throughout branch C.

Choose a functional \(J\) independent of \(K\). Identify the \(p\) classes
of \(B_{K,\alpha}\) with the points

\[
              a(q): K(a(q))=r,\quad J(a(q))=q,
              \qquad q\in\mathbf F_p.                    \tag{24}
\]

Choose \(h+1\) ordered pairs of distinct field elements whose union is all
of \(\mathbf F_p\): take \(h\) disjoint pairs covering \(p-1\) elements,
then pair the remaining element with one already used element. For the
\(i\)-th pair \((q_i,q_i')\), choose \(u_i\in\mathbf F_p\), all distinct,
and put

\[
 L_i=J+u_iK,\quad l_i=q_i+u_ir,\quad l_i'=q_i'+u_ir,      \tag{25}
\]

\[
 j_i={(r-s)l_i-(r+s)l_i'\over r}.                        \tag{26}
\]

For each pair, exclude the values of \(u_i\) making \(l_i\), \(l_i'\),
\(l_i+l_i'\), or \(j_i\) zero, and exclude the previously used values.
There are at most four pair-dependent exclusions. At the last step there
are at most \(h+4<p\) forbidden values when \(p>7\), so the distinct
\(u_i\) may be chosen greedily. This proves availability throughout the
stated range; it is not an empirical search.

Define \(M_i\) by (23). Equation (22) with \(\eta=+1\) holds by
construction. Thus every half \(H(L_i,M_i,j_i)\) supplies its two assigned
midpoint classes over the same difference block \(C\). The \(L_i\) are
pairwise distinct, and the \(p+1\) incidences cover all \(p\) midpoint
classes of \(B_{K,\alpha}\), with one repeated incidence.

This does not settle the actual branch-C problem. It chooses the centers
\(j_i\), so it need not respect the prescribed hard-center list. It also
does not prove that the full signed sum of halves is ternary, or that the
support of a complete nonzero dual word is contained in the resulting used
set. It proves exactly that distinct target directions alone supply no
product or counting obstruction to the saturated cover. Any next
obstruction must use the prescribed pairs \((L_i,j_i)\), mutual overlaps
and ternarity, or the rest of the dual-support equations.

## 7. Fail-when-wrong replay

The implementation checks the block intersections, generator projections,
one normalized midpoint conic, and the bounds (9) only at \(p=3,7\). It
also constructs the explicit saturated cover at \(p=31\). It does not
enumerate choices of halves or dual words.

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_symmetric_halved_mobius_cover.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_symmetric_halved_mobius_cover.py
