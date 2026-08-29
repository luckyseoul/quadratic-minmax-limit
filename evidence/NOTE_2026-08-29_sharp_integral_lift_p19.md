# Sharp integral quadratic-lift floor and the p=19 residue-zero reduction

**Date:** 2026-08-29
**Proposition:** 15.688
**Status:** proved; advances but does not close the \(p=19\) endpoint,
residual (ii), or \(L\)

## 1. Sharp lift theorem

Let \(p\ge5\) be odd, \(m=(p+1)/2\), and let \(B\not\equiv0\) be a
nonnegative, integer-valued polynomial of degree at most two on the middle
slice \(J(p,m)\). Then

\[
\boxed{4p\,\mathbb E B\ge p-3.}
\tag{1}
\]

The inequality is sharp: for two coordinates \(i\ne j\),

\[
B(x)=(1-x_i)(1-x_j)
\]

is Boolean and quadratic, and

\[
\mathbb E B
=\frac{\binom{p-2}{m}}{\binom pm}
=\frac{p-3}{4p}.
\tag{2}
\]

### Paired-cube input

Put \(H=\max B\) and choose \(X\) with \(B(X)=H\). Use the paired cube
through \(X\) from Propositions 15.680--15.681. Its averaging operator obeys

\[
TB(X)=\frac{H+p\mathbb E B}{p+1}.
\tag{3}
\]

The restriction of \(B\) to this cube is a nonzero degree-two polynomial,
so its support density is at least \(1/4\). It is also a nonnegative
integer-valued quadratic on a Boolean cube. Its unique multilinear monomial
coefficients are integers, by finite differences of its integer values.
Consequently its mean lies in \(\frac14\mathbb Z\): constant, linear, and
quadratic monomials have means \(1,1/2,1/4\).

If \(H=1\), the cube mean is at least \(1/4\), and (3) gives

\[
4p\mathbb E B\ge p-3.
\tag{4}
\]

If \(H\ge2\), the cube mean cannot equal \(1/4\). Equality would require
support density exactly \(1/4\) and value one at every support point, but
the cube contains \(X\) with value at least two. Quarter-integrality therefore
raises the cube mean to at least \(1/2\), and (3) gives

\[
4p\mathbb E B\ge2(p+1)-4H.
\tag{5}
\]

### Stabilizer input

Proposition 15.642's exact stabilizer identities retain the endpoint value
\(B(X)=H\), so

\[
4p\mathbb E B\ge
\begin{cases}
4H,&p\equiv3\pmod4,\\[1mm]
\dfrac{4r}{r+1}H,&p=4r+1.
\end{cases}
\tag{6}
\]

For \(p\equiv3\pmod4\), the maximum of the two affine lower bounds (5) and
(6) is minimized where \(H=(p+1)/4\), at value \(p+1\). For \(p=4r+1\),
it is minimized where \(H=r+1\), at value \(4r=p-1\). Hence every
\(H\ge2\) costs strictly more than \(p-3\), while the \(H=1\) branch gives
(4). This proves (1).

Equality in (1) forces \(H=1\), so \(B\) is Boolean. Equation (3) then
forces every paired cube through every support point to attain the minimum
degree-two support density \(1/4\).

## 2. The live p=19 consequence

At \(p=19\), the second even all-finite boundary has size \(s=16\) and pair
budget

\[
s(s-1)=240.
\]

The exact quotient/floor recurrence of Proposition 15.681 leaves one
phase-one row

\[
u_1=9,qquad D_1=126,qquad 9[b=2]+[b=16],
\]

and the pair-surviving phase-zero residues

\[
u_0\in\{0,2,3,4,6\}.
\]

For every positive \(u_0\), its quotient sum \(m-u_0\) is strictly below
the ten directions, so one direction has quotient zero. Its scaled mean is

\[
2u_0\in\{4,6,8,12\}.
\]

The least positive even fibre floor is \(20\), so that direction has
\(b=0\) and pointwise slack \(A=2B\), where \(B\) is a nonzero
nonnegative integer-valued quadratic. But (1) requires

\[
4p\mathbb E Bge p-3=16,
\]

excluding all four positive residues.

The sole arithmetic survivor is therefore

\[
\boxed{
u_0=0,quad u_1=9,quad
5[b=0]+5[b=16],quad 9[b=2]+[b=16].}
\tag{7}
\]

Its deficits are \(80\) and \(126\), leaving pair slack \(34\). The
quadratic-lift theorem does not exclude (7), so the \(p=19\) endpoint and
all top-level gates remain open.

## 3. Exact next lemma

Close the residue-zero incidence configuration (7): prove that no 16-point
finite boundary realizes those two directional profiles together with the
exact type budgets of Proposition 15.632. This one finite-geometric lemma
would close the entire \(p=19\) second endpoint.

## Reproduction

- symbolic theorem and exact p=19 row replay:
  `src/e1_gmin_m4_prop15688.py`;
- machine-readable record: `evidence/e1_gmin_m4_prop15688.json`;
- regression tests: `tests/test_prop15688.py`.
