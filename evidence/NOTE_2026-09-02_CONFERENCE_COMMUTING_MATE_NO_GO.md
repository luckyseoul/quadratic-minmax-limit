# Conference commuting mates: exact parity obstruction and literature audit

**Status:** exact all-orders obstruction for the proposed conference-mate
route, plus a live primary-source audit.  This does **not** close residual
`(ii)`.  It proves that a symmetric conference signing cannot commute with
any skew signing at all, not merely that symmetric and skew conference
matrices fail to coexist at the same order.

Let `A` be a symmetric zero-diagonal matrix with signs off the diagonal and
let `R` be a skew zero-diagonal matrix with signs off the diagonal.  The
opposite-diagonal block used in Proposition 6.5e/6.5i is

\[
 K_0(A,R)=\begin{pmatrix}A&R\\-R&-A\end{pmatrix},
\]

and direct multiplication gives

\[
 K_0(A,R)^2=
 \begin{pmatrix}
 A^2-R^2&AR-RA\\
 AR-RA&A^2-R^2
 \end{pmatrix}.                                      \tag{1}
\]

Thus a symmetric conference `A`, a skew conference `R`, and `AR=RA`
would indeed make (1) equal to `2(n-1)I`.  The commutation premise is
impossible.

## 1. The one-line diagonal obstruction

For every vertex `i`, symmetry of `A` and skewness of `R` give

\[
\begin{aligned}
 (AR-RA)_{ii}
 &=\sum_{j\ne i}\bigl(A_{ij}R_{ji}-R_{ij}A_{ji}\bigr)\\
 &=-2\sum_{j\ne i}A_{ij}R_{ij}.                       \tag{2}
\end{aligned}
\]

If `n` is even, the last sum contains the odd number `n-1` of signs.
Consequently

\[
 (AR-RA)_{ii}\equiv2\pmod4,
 \qquad |(AR-RA)_{ii}|\ge2                            \tag{3}
\]

for every `i`.  In particular,

\[
 \boxed{AR\ne RA},\qquad
 \boxed{\|AR-RA\|_F^2\ge4n}.                         \tag{4}
\]

Every real symmetric conference matrix has even order (in fact order
`2 mod 4`), so (4) applies to every possible order.  This proof uses no
conference identity for `R`: **there is no skew orientation signing that
commutes with a symmetric conference signing.**

For a primary reference on the order restriction, Mathon's opening
definition states that a symmetric conference matrix has order
`n congruent to 2 modulo 4`:

* Rudolf Mathon, *Symmetric Conference Matrices of Order \(pq^2+1\)*,
  Canadian Journal of Mathematics 30 (1978), 321--331,
  [DOI 10.4153/CJM-1978-029-1](https://doi.org/10.4153/CJM-1978-029-1).

## 2. Independent mod-two/eigenspace proof

There is a second proof which is useful for detecting mistaken relaxations.
Modulo two, every even-order skew signing has the same reduction

\[
 \bar R=J-I=J+I.
\]

Since `J^2=nJ=0` over `F_2` for even `n`,

\[
 (J+I)^2=I.                                           \tag{5}
\]

Hence every even-order skew signing has odd determinant and is invertible
over the reals.

If `A` is symmetric conference, `A^2=(n-1)I` and `tr(A)=0`.  Its two
eigenspaces for `+sqrt(n-1)` and `-sqrt(n-1)` therefore both have dimension
`n/2`.  For `n=2 mod 4`, these dimensions are odd.  If `AR=RA`, then `R`
preserves both eigenspaces.  Each restriction is a real skew-symmetric
operator on an odd-dimensional space, hence singular.  That contradicts
(5).  This independently proves the same stronger no-go.

## 3. The terminology trap in orthogonal-design papers

The standard orthogonal-design convention is

\[
 X,Y\text{ amicable}\iff XY^T=YX^T,
 \qquad
 X,Y\text{ anti-amicable}\iff XY^T=-YX^T.             \tag{6}
\]

See the definitions and constructions in Ebrahimian, Kharaghani, and
Suda, *Some Constructions for Amicable Orthogonal Designs*,
[arXiv:1509.03627](https://arxiv.org/abs/1509.03627).  When `A^T=A` and
`R^T=-R`, (6) reverses the word one might guess:

\[
 \begin{array}{c|c}
 \text{OD terminology}&\text{ordinary product relation}\\ \hline
 A,R\text{ amicable}&AR+RA=0,\\
 A,R\text{ anti-amicable}&AR-RA=0.
 \end{array}                                           \tag{7}
\]

The requested commuting pair is therefore an **anti-amicable** pair in
that literature.  Searches only for "amicable conference matrices" mainly
return the anticommuting relation and cannot support (1).

The newest close-looking result makes the mismatch explicit.  Glazyrin
defines an amicable Hadamard pair `(X,Y)` with `X` symmetric and
`Y=I+C` skew Hadamard.  Theorem 2 derives

\[
 XC+CX=0.                                              \tag{8}
\]

Corollary 3 supplies such Hadamard pairs for orders `2^ell`, `q+1` for
prime powers `q=3 mod 4`, certain core powers, and products of these
orders.  But `X` is a symmetric **Hadamard** matrix, not a zero-diagonal
symmetric conference matrix, and (8) is anticommutation, not the
commutation in (1).  In the symmetric-Paley case, Lemma 3.1 says exactly:

> If `q=1 mod 4` is a prime power and `Y` is the symmetric Paley
> conference matrix of order `q+1`, there exists a Hermitian complex
> Hadamard matrix `X` of order `q+1` such that `XY+YX=0`.

Source: Alexey Glazyrin, *New constructions of optimal arrangements of
\(2d\) lines in \(\mathbb C^d\)*, Theorem 2, Corollary 3, and Lemma 3.1,
[arXiv:2608.16116v1](https://arxiv.org/html/2608.16116v1).
This is a genuine infinite family, but both the matrix class and product
relation are wrong for the proposed `K_0` square.

Therefore the set of real orders carrying the requested symmetric-
conference/skew-conference commuting pair is **empty**, not merely sparse
or not known to be dense.

## 4. Near-conference and complex variants do not repair the bridge

### 4.1 Skew EW matrices

The most relevant real near-conference objects are skew EW matrices of
order `n=4t+2`.  Greaves--Suda prove that for the associated skew Seidel
matrix `S`, the rank-two square defect has, after switching, the exact form

\[
 (n+1)I+S^2
 =2\,\operatorname{diag}(J_{n/2},J_{n/2}).             \tag{9}
\]

See Lemma 4.8 of Gary Greaves and Sho Suda, *Symmetric and skew-symmetric
\(\{0,\pm1\}\)-matrices with large determinants*,
[arXiv:1601.02769](https://arxiv.org/abs/1601.02769).  Equation (9) is the
natural rank-two relaxation of a skew conference square in precisely the
right residue class.  It still cannot commute with `A`: `S` is a skew
signing of even order, so (2)--(5) apply.  In particular, passing from
conference matrices to the determinant-optimal EW class does not evade the
obstruction.

### 4.2 Complex conference matrices

Complex phases evade the real residue restrictions.  Et-Taoui constructs
parametric complex conference families, and Et-Taoui--Makhlouf classify
complex symmetric, skew-symmetric, and Hermitian conference matrices
through order six:

* Boumediene Et-Taoui, *Complex conference matrices, complex Hadamard
  matrices and complex equiangular tight frames*,
  [arXiv:1409.5720](https://arxiv.org/abs/1409.5720).
* Boumediene Et-Taoui and Abdenacer Makhlouf, *Complex skew-symmetric
  conference matrices*, Linear and Multilinear Algebra 70 (2022),
  6648--6663,
  [DOI 10.1080/03081087.2021.1967848](https://doi.org/10.1080/03081087.2021.1967848).
* Boumediene Et-Taoui, *A note on complex conference matrices*, Finite
  Fields and Their Applications 91 (2023), 102251,
  [DOI 10.1016/j.ffa.2023.102251](https://doi.org/10.1016/j.ffa.2023.102251).

These matrices have arbitrary unimodular off-diagonal entries.  They are
not real tournament signings.  Realification replaces a phase by a `2 by
2` rotation block: general phases produce entries outside `{0,+-1}`, and
even fourth-root phases produce monomial blocks with additional zeros on a
quadratic number of positions.  Neither operation yields the complete
real skew signing required in Proposition 6.5.  Filling those positions is
not a lower-order border operation.  Thus the complex results give no
valid orientation or real-signing upper bound.

## 5. Broader live literature search for the simultaneous orientation

The conference route was checked alongside the requested directed-cut,
structured-discrepancy, chaining, Gaussian-comparison, and coding-theory
search.  No theorem found has the needed quantifiers and sharp constant.
The closest primary results and their exact mismatches are as follows.

1. **Offset-preserving Gaussian comparison (closest analytic statement).**
   Van Handel's Theorem 1.3 says that if a finite centered process `X_t`
   has every increment dominated by the corresponding increment of a
   centered Gaussian process `G_t`, then, for arbitrary deterministic
   offsets `m_t`,

   \[
    E\sup_t(X_t+m_t)\le E\sup_t(cG_t+m_t)
   \]

   for a universal constant `c`.  This unusually does preserve the exact
   `A`-dependent defects by taking
   `t=(x,y,sigma)` and `m_t=-d_A(x,y)`.  However, every such scalar convex
   comparison already needs `c>=sqrt(pi/2)` (test a Rademacher variable
   with `f(u)=|u|`).  For a standard Gaussian skew matrix `G` and fixed
   Boolean `x`,

   \[
    E\max_y x^TGy=E\|G^Tx\|_1
      =n\sqrt{2(n-1)/\pi}.
   \]

   At the critical floor `M/n^(3/2) -> 1/pi`, the resulting Gaussian
   certificate is at least

   \[
    \left({1\over2}-{1\over\pi}+o(1)\right)n^{3/2},
   \]

   while the required width is only
   `((sqrt(2)-1)/pi+o(1)) n^(3/2)`.  Numerically these coefficients are
   `0.18169...` and `0.13185...`.  Hence even the sharp scalar comparison
   constant is too costly for this route.  Source: Ramon van Handel,
   *On the subgaussian comparison theorem*, Theorems 1.1 and 1.3,
   [arXiv:2512.18588](https://arxiv.org/abs/2512.18588).

2. **Bernoulli-process decomposition and chaining.**  Bednorz--Latała
   prove that for `T subset ell_2(I)` with
   `b(T)=E sup_(t in T) sum_i t_i epsilon_i<infinity`, one can decompose
   `T subset T_1+T_2` with
   `sup_(t in T_1)||t||_1 <= L b(T)` and
   `E sup_(t in T_2) sum_i t_i g_i <= L b(T)`, for a universal `L`.
   Liu--Zadik's 2026 comparison gives universal-constant equivalences
   among their Bernoulli functional, an integrated randomized-Dudley
   functional, and a decomposition functional.  Both are expectation-of-
   supremum theorems for independent signs; neither produces one signing
   lying in all of the nonuniform slabs
   `|<b_xy,r>| <= (sqrt(2)-1)M+d_A(x,y)`, and neither supplies the sharp
   constant.  Sources:
   [Bednorz--Latała, arXiv:1305.4292](https://arxiv.org/abs/1305.4292) and
   [Liu--Zadik, arXiv:2608.11031](https://arxiv.org/abs/2608.11031).

3. **Constructive partial coloring.**  Lovett--Meka's Theorem 4 applies to
   vectors `v_j in R^N` and thresholds `c_j>=0` only under

   \[
    \sum_j e^{-c_j^2/16}\le N/16.                     \tag{10}
   \]

   It then finds, with probability at least `0.1`, a point respecting
   `|<x-x_0,v_j>| <= c_j||v_j||_2` while making at least half its
   coordinates nearly integral.  In the Hamming-central layer here there
   are `exp((log 4-o(1))n)` ordered state pairs, while the universal upper
   scale `M <= (1/2+o(1))n^(3/2)` gives
   `c_j^2/16 <= (1/8+o(1))n`.  The left side of (10) therefore remains at
   least `exp((log 4-1/8-o(1))n)`, whereas `N/16=Theta(n^2)`.  The theorem
   is exponentially outside its hypothesis before the first partial-
   coloring step.  Source: Shachar Lovett and Raghu Meka,
   *Constructive Discrepancy Minimization by Walking on the Edges*,
   [arXiv:1203.5747](https://arxiv.org/abs/1203.5747).

4. **Gaussian discrepancy.**  Chewi--Gerber--Rigollet--Turner define
   `Gdisc(B)=min E||Bg||_infinity` over centered Gaussian `g` with covariance
   diagonal one, and prove `Gdisc(B) <= sqrt(2/pi) disc(B)`.  Their online
   theorem is for Gaussian rank at least two; the paper explicitly leaves
   rounding Gaussian discrepancy to a rank-one sign vector as an open
   problem.  This is the missing quantifier in the present application,
   not a bridge.  Source:
   [arXiv:2109.08280](https://arxiv.org/abs/2109.08280).

5. **Spherically punctured Reed--Muller codes.**  The augmented cut code is
   the degree-one Reed--Muller code evaluated on the weight-two slice, so
   the unsigned minimax value is its covering radius.  Dumer--Kapralova
   and Gini--Meaux develop this exact code family, but the general covering
   radius is itself unknown, and a covering-radius theorem would still not
   choose an orientation for the exponentially many `A`-dependent signed
   bivector slabs.  Sources:
   [Dumer--Kapralova thesis](https://escholarship.org/uc/item/36d0023f) and
   [Gini--Meaux, IACR ePrint 2022/408](https://eprint.iacr.org/2022/408.pdf).

Targeted GitHub searches for the exact spherically-punctured covering
radius, tournament cut-norm skew signing, and directed-half-cut phrases
returned no implementation or theorem matching the target.  The targeted
Zenodo search returned no conference/amicable bridge and one tournament-
discrepancy manuscript about maximum acyclic subgraphs, not signed cuts:
Vladimir Riabov, *A Spectral Barrier for Paley Tournament Discrepancy*,
[Zenodo 20670429](https://zenodo.org/records/20670429).  The original
MathOverflow question still has no posted solution.

**Live-search conclusion:** no serious theorem located on arXiv,
MathOverflow, journal sites, GitHub, or Zenodo proves the simultaneous
signed directed-half-cut orientation inequality.  The newest algebraic
lead (Glazyrin) is rigorously in the wrong relation/class, and the closest
analytic lead (van Handel) preserves the right offsets but already loses
too much at the best possible scalar comparison constant.

## 6. Exact replay

The checker verifies (2)--(5), the terminology translation (7), Paley
conference squares at small prime orders, and all `2^15=32768` skew
signings against the order-six Paley conference matrix.  The finite
exhaustion finds zero commuting mates; its smallest diagonal commutator
Frobenius square is `24=4n`.  It is a regression check, not the proof.

```bash
PYTHONPATH=src pytest -q tests/test_conference_commuting_mate_no_go.py
PYTHONPATH=src python src/conference_commuting_mate_no_go.py
```

Fresh result: `25 passed`.  The full order-six audit also reports minimum
full commutator Frobenius square `104`.  Residual `(ii)` remains open; only
the conference commuting-mate route is closed.
