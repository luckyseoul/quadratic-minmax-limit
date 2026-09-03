# Grouped uncertainty from the canonical square-root remainder

Date: 2026-09-03

Status: the grouped uncertainty inequality is proved for every support and
every odd prime.  Consequences for the two-dimensional halved row code and
the structured Möbius puncture are deliberately separated below; this note
does not by itself classify minimum words, prove puncture robustness, solve
the symmetric Boolean box, close residual (ii), close E1, or solve the
original MathOverflow problem.

## 1. Statement

Let

\[
 \Delta=(\mathbb F_p^2\setminus\{0\})/\{+1,-1\}
\]

for an odd prime \(p\).  The paired non-origin affine-block incidence matrix
is denoted by \(M\), and its \(h=(p-1)/2\) output coordinates in projective
direction \(A\) form the group \(B_A\).  Then every nonzero binary point word
\(f\) satisfies

\[
 \boxed{\operatorname{wt}(f)+
 \#\{A:(M^Tf)|_{B_A}\ne0\}\ge p+1.}                 \tag{1}
\]

Equivalently, if \(S=\operatorname{supp}f\), \(s=|S|\), and \(z\) is the
number of silent direction groups, then

\[
 z\le s.                                                   \tag{2}
\]

For odd \(s\), (2) is the existing radial-parity theorem.  The rest of this
note proves it for even \(s\).

## 2. A canonical square-root remainder

Write \(s=2n\), choose distinct antipodal representatives \(v_1,\ldots,v_{2n}\),
and put \(\ell_i(u)=u(v_i)\) for a projective functional parameter
\(u=(U,V)\).  Introduce a second variable \(X\), assigned homogeneous weight
two, and form

\[
 P_u(X)=\prod_{i=1}^{2n}\bigl(X-\ell_i(u)^2\bigr)
       =X^{2n}+c_1(u)X^{2n-1}+\cdots+c_{2n}(u),             \tag{3}
\]

where \(c_k\) is a homogeneous binary form of degree \(2k\).

Because 2 is invertible in \(\mathbb F_p\), there is a unique monic

\[
 Q_u(X)=X^n+q_1(u)X^{n-1}+\cdots+q_n(u)                    \tag{4}
\]

whose square agrees with the top \(n\) non-leading coefficients of \(P\).
Indeed, recursively,

\[
 q_k=\frac12\left(c_k-
       \sum_{\substack{i+j=k\\1\le i,j<k}}q_iq_j\right),
       \qquad 1\le k\le n.                                \tag{5}
\]

Each \(q_k\) is homogeneous of degree \(2k\).  Thus

\[
 P_u(X)-Q_u(X)^2
   =\sum_{k=n+1}^{2n}R_k(u)X^{2n-k},                       \tag{6}
\]

with \(R_k\) homogeneous of degree \(2k\le4n=2s\).

At least one \(R_k\) is nonzero.  Otherwise (3) would be a square in the
unique factorization domain \(\mathbb F_p[U,V,X]\).  Its monic linear factors
in \(X\) are distinct: equality
\(\ell_i^2=\ell_j^2\) as binary forms gives
\((\ell_i-\ell_j)(\ell_i+\ell_j)=0\), hence
\(v_i=\pm v_j\), contrary to distinctness in \(\Delta\).  Therefore (3) is
squarefree as a polynomial in \(X\), not a square.

## 3. Silent directions are double zeros

Fix a silent projective direction \(u_0\).  The nonzero squared-projection
fibres \(\ell_i(u_0)^2\) have even multiplicity by definition.  Their total
size and \(s\) are even, so the zero (radial) fibre also has even
multiplicity.  Hence every root of \(P_{u_0}(X)\) occurs evenly.

Use a local parameter \(t\) at \(u_0\).  Pair indices inside each equal-root
fibre.  If one pair has

\[
 \ell_i(t)^2=a+\alpha t+O(t^2),\qquad
 \ell_j(t)^2=a+\beta t+O(t^2),
\]

then

\[
\begin{aligned}
 &(X-\ell_i(t)^2)(X-\ell_j(t)^2)\\
 &\qquad=\left(X-a-\frac{\alpha+\beta}{2}t\right)^2
          \pmod {t^2}.                                     \tag{7}
\end{aligned}
\]

Multiplying (7) over all pairs shows that \(P_t(X)\) is the square of a
monic degree-\(n\) polynomial modulo \(t^2\).  The recursive construction
(5) is the unique monic top-half square root over the local ring
\(\mathbb F_p[t]/(t^2)\).  Consequently its remainder (6) is zero modulo
\(t^2\):

\[
 \operatorname{ord}_{u_0}R_k\ge2
 \quad\text{for every }k.                                  \tag{8}
\]

This local argument also covers a radial root \(a=0\).  The projective point
at infinity is covered by taking the other affine chart.

## 4. Degree count

Choose a nonzero coefficient \(R_k\) from (6).  It is a nonzero homogeneous
binary form of degree \(2k\), and each of the \(z\) distinct silent
projective directions is a root of multiplicity at least two.  Therefore

\[
 2z\le 2k\le4n=2s,
\]

which proves \(z\le s\) for even support.  Combining this with radial parity
for odd support proves (1) for every support and every odd prime.  No
congruence assumption on \(p\) is needed.

## 5. Scope of the immediate consequences

The earlier row-code note proves the conditional implication

\[
 \text{(1)}\quad\Longrightarrow\quad d(\operatorname{Row}(D))\ge ph.
\]

The fixed-transverse rectangles already have weight \(ph\), so after that
conditional derivation is independently rechecked, (1) gives equality of
the minimum distance.  It does not automatically classify all minimum words.

Likewise, a puncture of size strictly below \(ph\) cannot contain a nonzero
dual support once the row-code consequence is imported.  An actual Möbius
puncture may be larger than this, so (1) alone does not establish its
surjectivity or solve the divided Boolean fibre.

## 6. Executable transcription checks

- `src/e1_gmin_m4_grouped_uncertainty_square.py`
- `tests/test_grouped_uncertainty_square.py`

The replay constructs (3)--(6) for small explicit supports and verifies
double projective vanishing at their silent directions.  These checks are
fail-when-wrong tests, not computational evidence for the theorem.
