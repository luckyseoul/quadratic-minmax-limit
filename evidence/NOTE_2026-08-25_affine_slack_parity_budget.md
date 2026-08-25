# Affine slack budget and quadratic parity lifts

Date: 2026-08-25.  This is Proposition 15.632.  It gives an all-prime
necessary condition for the non-Walsh multi-level residual and excludes its
Eulerian-boundary branch.  It does **not** close residual (ii), Type I, R1,
global QVAR, or the limit.

## 1. Setup

Let `H` be an odd set of signed Paley edges and write

\[
 S_H(y)=\sum_{e\in H}C_e y_e,\qquad h=|H|.
\]

For a projective \(\mathbb F_p\)-direction \(d\), let \(\epsilon_d\) be the
quadratic type of its kernel.  There are \(m=(p+1)/2\) directions of each
type.  The affine halfspaces in direction \(d\) are Boolean
\(\epsilon_dp\)-eigenvectors and are indexed by the middle slice
\(J(p,m)\).

If `H` separates both affine eigenshells with margin three, then

\[
 A_d(y)={\epsilon_d S_H(y)-3\over2}
\]

is a nonnegative integer-valued quadratic on that middle slice.  Put

\[
 a_d=2p\,\mathbb E_d A_d.
\]

## 2. Exact directional budgets

Multiply the directional mean contribution of one edge by \(p\epsilon_d\).
The result is

| edge relative to the direction | contribution |
|---|---:|
| incident with infinity | \(1\) |
| finite and parallel | \(p\) |
| finite and transverse | \(-\epsilon_d C_e\) |

Thus \(a_d\) is an even integer: every edge contribution is odd modulo two,
while \(3p\) is odd and \(h\) is odd.  Nonnegativity follows from the
separator hypothesis.

Summing over all directions, each edge contributes \(p+1\).  More strongly,
each edge contributes \(m\) inside *each* quadratic-type half.  For an
infinity edge this is immediate.  For a finite edge of sign \(c\), the
type-\(c\) half contributes

\[
 p-(m-1)=m,
\]

and the type-\(-c\) half contributes \(m\) transverse copies of `+1`.
Consequently

\[
 \boxed{\sum_{d:\epsilon_d=\tau}a_d
       =m(h-3p)\quad(\tau=\pm1)},                         \tag{1}
\]

and the total is \((p+1)(h-3p)\).  At residual size \(h=4p+1\), each
quadratic-type half has budget \((p+1)^2/2\), and the total is \((p+1)^2\).

## 3. Boundary parity

Let \(D\) be the odd-degree boundary of `H` and
\(c_H=\prod_{e\in H}C_e\).  Since

\[
 \prod_{e\in H}C_e y_e=c_H\prod_{v\in D}y_v,
\]

direct reduction modulo two gives

\[
 (-1)^{A_d}
 =\epsilon_d(-1)^{(h-3)/2}c_H\prod_{v\in D}y_v.           \tag{2}
\]

Let \(B_d\) be the affine fibres in direction \(d\) containing an odd
number of finite points of \(D\), and put \(b_d=|B_d|\).  If \(x_s\) is the
middle-slice membership bit of fibre \(s\), then

\[
 \prod_{v\in D}y_v
 =\epsilon_d^{\mathbf1_{\infty\in D}}
  (-1)^{b_d+\sum_{s\in B_d}x_s}.
\]

These two product formulas therefore give an explicit phase \(\eta_d\) such that

\[
 A_d(x)\equiv\sum_{s\in B_d}x_s+\eta_d\pmod2.             \tag{3}
\]

The implementation verifies (3) point by point, not only in expectation.

## 4. Exact quadratic majorant

Average \(A_d\) under
`Sym(B_d) x Sym(B_d^c)`.  The result is a quadratic \(q(t)\), where
\(t=|X\cap B_d|\).  Every value entering the average is a nonnegative
integer with parity \(t+\eta_d\), hence

\[
 q(t)\ge (t+\eta_d\bmod2)                                 \tag{4}
\]

on the hypergeometric support

\[
 \Pr(t)=\frac{\binom{b_d}{t}\binom{p-b_d}{m-t}}
               {\binom pm}.
\]

Define \(M(p,b,\eta)\) as the minimum hypergeometric expectation of a
quadratic satisfying (4).  This is a three-variable linear program.  With
at least three support points its vertices have three active evaluation
constraints, so interpolation over all triples computes it exactly over
`Fraction`; one- and two-point supports are immediate.

Since \(a_d\) is even,

\[
 a_d\ge2\lceil pM(p,b_d,\eta_d)\rceil.                    \tag{5}
\]

Combining (1) and (5) gives the type-split necessary conditions

\[
 \boxed{\sum_{d:\epsilon_d=\tau}
 2\lceil pM(p,b_d,\eta_d)\rceil
 \le {p+1\over2}(h-3p),\qquad\tau=\pm1.}                 \tag{6}
\]

This uses integer *magnitude*, not only Walsh parity.

## 5. An all-prime branch kill

At \(h=4p+1\), suppose `H` is Eulerian, so \(D=\varnothing\).  Then
\(b_d=0\) in every direction and the residual phase is
\(-\epsilon_dc_H\).  One quadratic-type half therefore has phase one.
For `b=0, phase=1`, the parity function is constantly one, so
\(M(p,0,1)=1\) and every direction in that half costs \(2p\).  Its required
cost is

\[
 p(p+1),
\]

but its exact budget is only \((p+1)^2/2\).  The contradiction gap is

\[
 p(p+1)-{(p+1)^2\over2}={p^2-1\over2}>0.                 \tag{7}
\]

Hence the Eulerian-boundary residual branch is empty for every odd prime.

## 6. Why this is not the full close

The corrected affine model at \(p=5\) has a genuine integral solution.  For
the distinguished edge `(0,1)`, take the twenty-edge `G` stored as
`AFFINE_P5_G` in the source and set `H=G union {(0,1)}`.  Direct enumeration
of all six affine middle-slice families gives

\[
 (a_d)=(12,4,0,6,10,4),\qquad \sum a_d=36=(p+1)^2,
\]

with slack supports

\[
 \{0,1,2,3,4\},\{0,2\},\{0\},\{0,2\},\{0,2,4\},\{0,2\}.
\]

Its boundary is infinity together with one affine line.  Thus (6) is a real
necessary condition and kills a general branch, but it is not secretly a
claim that the affine system—or the full residual—is empty.

## 7. Literature check

Searches on 2026-08-25 for the exact parity-majorant/type-split mechanism
found adjacent, not duplicate, frameworks:

- low-degree analysis on Boolean slices / Johnson schemes:
  [Filmus--Kindler--Mossel--Wimmer](https://arxiv.org/abs/1504.01689) and
  [Kiermaier--Mannaert--Wassermann](https://arxiv.org/abs/2405.07572);
- nonnegative quadratic polynomials on finite Boolean sets:
  [Blekherman--Gouveia--Pfeiffer](https://arxiv.org/abs/1402.4199);
- finite-geometry incidence codes:
  [Sin--Sorci--Xiang](https://arxiv.org/abs/1908.06824).

No searched source stated the combined edge-boundary parity lift, exact
hypergeometric quadratic-majorant LP, or square/nonsquare directional budget
(6).  This is a search record, not a novelty claim.  There is no nontrivial
new integer sequence here requiring an OEIS identification; the explicit
Eulerian gap `(p^2-1)/2` is already a closed polynomial.

## Reproduction

- theorem and exact tables: `src/e1_gmin_m4_prop15632.py`;
- machine-readable evidence: `evidence/e1_gmin_m4_prop15632.json`;
- corrected finite affine/full-shell solver:
  `scripts/residual_affine_johnson_milp.py`;
- pointwise, random-edge, symmetry, and branch-kill tests:
  `tests/test_prop15632.py`.
