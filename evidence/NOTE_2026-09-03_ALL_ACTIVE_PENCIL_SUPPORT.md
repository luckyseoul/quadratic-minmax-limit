# The all-active hard-star equality case has \(c\ge p-1\)

**Status:** exact all-prime exclusion of the \(c=p-2\) equality case when
every hard center is nonzero, for balanced branch C at \(p=4r+3\ge31\).
This is a support sharpening, not a common Boolean lift and not a closure of
residual (ii).

## 1. Equality input

Assume zero odd global forms.  The hard compact residual and every opposite
row are central, so the antisymmetric target consists of the hard unit-star
differences.  Let \(c\) be the number of nonfixed inversion orbits occupied
on exactly one side by a putative graph.

The hard-star support theorem gives \(c\ge p-2\).  Suppose equality holds
and every one of the

\[
                         m={p+1\over2}=2r+2               \tag{1}
\]

hard centers is nonzero.  Since \(m\ge16>8\), the equality-pencil theorem
applies: all \(p-2\) single source-edge orbits have a common projective
endpoint \([P]\), and

\[
                         j_L^2=L(P)^2                     \tag{2}
\]

for every hard functional \(L\).

Choose in each inversion orbit the representative whose common endpoint is
\(P\), and let \(Q\) be the set of its other endpoints.  No member of
\(Q\) is \(P\), because source loops are absent, or \(-P\), because these
are nonfixed inversion orbits.  Thus \(|Q|=p-2\).

## 2. Exact two-point completion

For a hard row, equality means that its \(p-2\) selected source orbits map
bijectively to the \(p-2\) central cell-orbits in \(A_{j_L}\).  Put
\(c_L=L(P)\ne0\), where nonvanishing follows from (2) and the active-center
hypothesis.  The central cell-orbits incident with \(c_L\) are

\[
 [\{c_L,a\}],\qquad a\in\mathbf F_p\setminus\{c_L,-c_L\}. \tag{3}
\]

The map \(a\mapsto[\{c_L,a\}]\) is injective on this domain.  Indeed, an
equality with the negative cell would force \(a=-c_L\), which (3) excludes.
The sign choice \(j_L=\pm c_L\), the Paley column sign, and the choice of
which physical edge is selected from an inversion orbit affect coefficients
and orientations, but not this cell support.  Equality puts exactly one
single source orbit in every cell orbit in (3).  Consequently hard-row
bijectivity says exactly

\[
             \boxed{L(Q)=\mathbf F_p\setminus\{L(P),-L(P)\}.} \tag{4}
\]

Adjoin the two missing points:

\[
                            S=Q\cup\{P,-P\}.              \tag{5}
\]

For every hard \(L\), equation (4) makes \(L:S\to\mathbf F_p\) a
bijection.  Hence none of the \(m\) spatial directions annihilated by the
hard functionals is determined by a pair of points of \(S\).  The
\(p\)-point set \(S\) therefore determines at most

\[
                         (p+1)-m={p+1\over2}              \tag{6}
\]

affine directions.

The prime-order Rédei--Megyesi direction theorem says that, for odd prime
\(p\), a noncollinear
\(p\)-point subset of \(\operatorname{AG}(2,p)\) determines at least
\((p+3)/2\) directions.  Equations (5)--(6) force \(S\) to be collinear.
This avoids any unproved near-maximum Paley-clique stability statement.
We use the form recorded in T. Szőnyi, *On the number of directions
determined by a set of points in an affine Galois plane*, JCTA 74 (1996),
141--146.

Because \(L(P)\ne0\) for every hard \(L\), we have \(P\ne0\), so \(P\) and
\(-P\) are distinct.  Their affine line contains their midpoint, the origin;
it is the linear line \(\langle P\rangle\).  Since it already
contains \(p\) points, \(S=\langle P\rangle\), and

\[
                       Q=\langle P\rangle\setminus\{P,-P\}. \tag{7}
\]

## 3. The opposite parallel-row contradiction

Every one of the \(p-2\) single graph edges is either a representative
\(\{P,q\}\), \(q\in Q\), or its antipode \(\{-P,-q\}\).  Both are parallel
to \(\langle P\rangle\).  Let \(M\) annihilate \(P\).  It is an opposite
direction: if it were hard, (2) and the active-center hypothesis would give
\(M(P)\ne0\).

Thus the actual parallel count in row \(M\) is at least \(p-2=4r+1\).
On the full balanced branch-C ray, every opposite quota obeys

\[
                             Q_M\le2r+2.                  \tag{8}
\]

The bound follows already at the upper endpoint: the opposite total is

\[
 10r+6+(4r^2-2r-5)=4r^2+8r+1,
\]

whose ceiling after division among \(2r+2\) balanced rows is \(2r+2\).
Since \(4r+1>2r+2\), equations (7)--(8) contradict the target parallel
row.  Additional double or fixed source orbits only increase that
nonnegative actual parallel count.

Therefore equality is impossible and

\[
 \boxed{c\ge p-1\quad\text{when all hard star centers are nonzero}.} \tag{9}
\]

This theorem says nothing when some hard center is zero, and \(c\ge p-1\)
alone does not solve the coupled symmetric Boolean fibre.  Residual (ii),
E1, \(L=1/2\), and the original MathOverflow problem remain OPEN.

## Reproduction

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
      tests/test_all_active_pencil_support.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_all_active_pencil_support.py
