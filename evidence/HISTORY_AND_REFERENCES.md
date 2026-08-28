# History of the problem, and older references

**Date:** 2026-08-22  
**Sources:** live fetches (MathOverflow API, X thread, Paata’s 2019 blog and CV, arXiv, Goethals–Seidel PDF, Wikipedia/Paley construction).  
**Not a close:** nothing here proves that \(\lim\alpha_n\) exists. Do not treat sandwich, Paley \(\rho=1\), or conference matrices as settlement.

Related in-repo notes: `evidence/MO_THREAD_REAUDIT.md` (thread-only time-savers), `LONG_HORIZON_GOAL.md`, `solution.md`.

---

## 1. The question

Paata Ivanisvili, MathOverflow [413935](https://mathoverflow.net/questions/413935) (16 Jan 2022; 0 answers as of the 2026-08-22 scrape):

\[
\lim_{n\to\infty}\,
n^{-3/2}
\min_{a_{ij}=\pm1}
\max_{x_j=\pm1}
\Bigl|\sum_{1\le i<j\le n} a_{ij}x_ix_j\Bigr|.
\]

He wrote that there is **no significant motivation** — “pure curiosity” — and that he is **not** interested in numerics or in bounds unless they give significant evidence that the limit **does not** exist. If it exists, he wants the proof.

Tags used: combinatorics, discrete geometry, Fourier analysis, eigenvalues, hypercube. In comments he wished he could also tag spin glasses, polynomials, max-cut, and **Littlewood’s 4/3 inequality** (tag limit \(<6\)).

This repo writes \(m_n\) for the inner min-max and \(\alpha_n=m_n/n^{3/2}\). The long-horizon goal is existence of \(L=\lim\alpha_n\) (or a proof that it fails), not a better sandwich.

---

## 2. What the MO comments actually say

Fetched via StackExchange API (`questions/413935/comments`, 11 comments, `has_more=false`) plus a Firecrawl scrape of the HTML.

| Who | Date | Load-bearing content |
|-----|------|----------------------|
| **domotorp** | 16 Jan 2022 | Seidel rewrite \(\frac12\lim\min_A\max_x n^{-3/2}x^\top Ax\). Retracted dropping the absolute value. “More eigenvalues than discrepancy.” |
| **Paata** | 17 Jan 2022 | Agrees limsup finite. Liminf \(\ge C>0\); \(C=2^{-5/2}\) if the calculus is right. Those arguments **do not** prove or disprove existence. |
| **Shannon Starr** | 17 Jan 2022 | Limsup finite by Talagrand (typical random \(A\)); first guessed the limit is 0, then retracted. Pointed at Cris Moore, Lenka Zdeborová, Aukosh Jagannath / Ben Arous large deviations, and Garry Bowlin (bipartition of \(K_n\)). |
| **Paata** | 18 Jan 2022 | Liminf via Defant–Mastyło–Pérez on a 2-homogeneous polynomial; link to the 2019 blog. Thanks for names. |
| **Andrei Z.** | 25 Jul 2026 | Put the sequence in OEIS (same day as the X prize). |
| **Rafael Hipólito** | 25 Jul 2026 | Antisymmetric \(A\) would give 0 — then noticed the sum is only the upper triangle. |

No later mathematical answer on MO. Starr’s Bowlin pointer is **maximum frustration in bipartite signed graphs** (Electron. J. Combin. 19 (2012), #P10); Starr himself said it probably does not have good enough results. Do not reopen multipartite comparison as existence (`solution.md` §9–§10).

---

## 3. The X prize

[PI010101/status/2081070728422752329](https://x.com/PI010101/status/2081070728422752329) (25 Jul 2026):

- The question arose **almost accidentally** while he was working on something else.
- About **five hours** of his own time, no solution.
- Least-upvoted of his MO questions.
- Posted before AI was good at mathematics; the challenge is whether AI can beat those five hours.
- Rule: full AI transcript, **no substantial human mathematical hints**.
- Prize via X Money; he will check the proof himself.

Replies on the original post are not a proof. A later remark that a Grok 4.5 writeup of \(\lim=1/2\) had a gap is **not** in the original thread body.

Author constraint, already in the MO post: sandwich (liminf \(>0\), limsup \(<\infty\)) does **not** settle existence. He said so in 2022.

---

## 4. How Paata learned the object (education, not Paley 1933)

Live CV and bios (Princeton PDF, Simons, UCI ScholarConnect, Google Scholar):

| Years | Fact |
|-------|------|
| 2006–2011 | B.S./M.S. mathematics, St. Petersburg State University, diploma with honor |
| 2011–2013 | Research engineer, Chebyshev Laboratory (Smirnov) |
| 2011–2015 | Dual Ph.D.: SPbU + Michigan State, **Alexander Volberg** |
| First paper (2010, undergraduate) | with **S. V. Kislyakov**: sparse spectrum / uniformly convergent Fourier series |

That is the St Petersburg analysis line (Havin–Kislyakov–Vasyunin–Stolyarov–Zatitskiy–Osipov), then Volberg’s Bellman / Hamming-cube school. Volberg’s INRIA Bellman lectures are [arXiv:1106.3899](https://arxiv.org/abs/1106.3899) (June 2011), the summer Paata starts the Ph.D.

**Textbooks and papers that sit in that pipeline**, not a Paley-design course:

- Havin / Nikolski harmonic analysis; Kislyakov
- Koosis on \(H^p\)
- Kolmogorov–Fomin, Natanson (SPbU undergraduate analysis)
- Stein (singular integrals) as the international companion
- Volberg: Bellman functions, Burkholder, hypercontractivity on the cube
- Later: Fourier analysis on \(\{\pm1\}^n\), Bonami–Beckner, Littlewood 4/3, Blei

Paley conference matrices are **not** a 2006–2011 SPbU analysis staple. They enter this problem from the combinatorial side once one asks which \(A\) realises a small \(\max_x|x^\top Ax|\).

The 6 Dec 2019 blog
[A little bit of Fourier analysis](https://extremal010101.wordpress.com/2019/12/06/a-little-bit-of-fourier-analysis/)
is the explicit context clue. He frames the **dual** of the prize: Dean’s problem / max-cut / SK Hamiltonian

\[
\mathbb{E}_A\max_{x\in\{\pm1\}^n}\sum_{i<j}a_{ij}x_ix_j\sim c\,n^{3/2},
\]

Parisi–Talagrand free energy, then the **deterministic** coefficient inequality that gives his liminf: **Defant–Mastyło–Pérez** [arXiv:1706.03670](https://arxiv.org/abs/1706.03670) (Boolean Fourier spectrum of degree-\(d\) functions). For \(d=2\) that is a 4/3-type bound on \(\sum|a_{ij}|^{4/3}\) versus \(\|f\|_\infty\), which is how he gets \(2^{-5/2}\).

He traces the **ideas** of DMP to:

1. **Littlewood 1930**, *On bounded bilinear forms in an infinite number of variables*, Quart. J. Math. (Oxford) **1**, 164–174 — the 4/3 inequality, forerunner of Grothendieck.
2. **Ron Blei**, *Analysis in Integer and Fractional Dimensions*, CUP 2001, **Theorem 34** — existence of a dimension-free constant \(B(d)\) for homogeneous Walsh polynomials.
3. Bayart–Pellegrino–Seoane-Sepúlveda (and Núñez-Alarcón–Pellegrino–Seoane-Sepúlveda): the constant down to \(C^{\sqrt{d\ln d}}\).

Hypercontractivity / Harris-type lemmas are in the same blog proof, not Paley graphs.

So: he learned the min-max as a **Fourier / polynomial / SK dual**, not as a 1933 Hadamard construction. The Paley matrix is a later candidate for the **min**, not the source of the question.

---

## 5. Pre-internet lineage of the *object*

The bilinear form on \(\{\pm1\}^n\) with coefficients \(\pm1\) is old. The **existence of this specific limit** is not.

| Year | Source | What it actually gives |
|------|--------|------------------------|
| 1923 | Khintchine | \(L^p\) norms of Rademacher sums — typical size, not \(\min_A\max_x\) |
| **1930** | **Littlewood 4/3** | Deterministic coefficient vs \(\|f\|_\infty\) for bilinear forms. Ancestor of Paata’s liminf. |
| 1932 | Paley–Zygmund | Random series; same Paley, different paper |
| **1933** | **R. E. A. C. Paley, *On Orthogonal Matrices*, J. Math. Phys. 12, 311–320** | Conference / Hadamard matrices from quadratic residues. First wording of the Hadamard conjecture. Motivated by **Coxeter polytopes**, not by min-max of \(x^\top Ax\). Same volume as Todd/Coxeter. Paley used “orthogonal” for mutually orthogonal \(\pm1\) rows. |
| 1950s | Belevitch | Names “conference matrices” (telephony networks) |
| 1966–67 | van Lint–Seidel; **Goethals–Seidel**, *Orthogonal matrices with zero diagonal*, Canad. J. Math. **19** (1967), 1001–1010 | \(C\)-matrices; Paley equivalent to a negacyclic \(C\)-matrix; \(CC^\top=qI\). Motivated by Coxeter polytopes (again) and equilateral point sets in elliptic geometry. Dedicated to Coxeter. |
| 1971 | Delsarte–Goethals–Seidel, *Orthogonal matrices with zero diagonal II* | Same family |
| 1985 | Spencer, linear discrepancy | **Different problem.** domotorp already pushed away from discrepancy tags. |
| 2001 | Blei, Thm 34 | Dimension-free \(B(d)\) for Walsh polynomials |
| 2006 | Talagrand (Parisi formula) | **Typical** \(A\): \(\mathbb{E}\max_x\sim c\,n^{3/2}\). Gives Paata’s limsup, not the min. |
| 2017 | DMP [1706.03670](https://arxiv.org/abs/1706.03670) | Paata’s published liminf route |
| 2017 | Jones, [arXiv:1702.00285](https://arxiv.org/abs/1702.00285) | History of Paley graphs; Paley’s 1933 paper is **only** Hadamard construction |
| 2021–22 | Goryainov–Lin [2104.08839](https://arxiv.org/abs/2104.08839), Goryainov–Shalaginov–Yip [2203.16081](https://arxiv.org/abs/2203.16081) | Eigenfunctions / canonical cliques of Paley of **square order**. Already cited in-repo. They describe eigenspaces of \(P(p^2)\), not a proof that \(\lim\alpha_n\) exists. |

Paley 1933 constructs a family with \(C^\top C=qI\), so \(\|Cy\|_2=\sqrt{q}\|y\|_2\). For \(n=q+1=p^2+1\) that is the conference graph used in this repo. That gives a **sequence** along which \(\max_x|x^\top Ax|\) is order \(n^{3/2}\), hence a candidate for \(\limsup\alpha_n\le 1/2\) once Paley is shown near-optimal in the \(\rho=1\) sense. Paley never asked whether \(\lim\alpha_n\) exists.

Goethals–Seidel 1967, opening paragraph (Cambridge PDF): questions in the theory of polytopes, posed by Coxeter, led Paley to Hadamard matrices in which he used \(C\)-matrices; van Lint–Seidel discussed \(C\)-matrices for equilateral point sets; Belevitch initiated the study in conference telephony.

---

## 6. Two problems that share a formula and are not the same

The 2019 blog and the 2022 MO question are **accidental duals**:

| | Typical \(A\) (SK / Dean / blog) | Worst \(A\) (MO / this repo) |
|--|----------------------------------|------------------------------|
| Object | \(\mathbb{E}_A\max_x\sum a_{ij}x_ix_j\) | \(\min_A\max_x\bigl|\sum a_{ij}x_ix_j\bigr|\) |
| Scaling | \(n^{3/2}\) | \(n^{3/2}\) |
| Status | Parisi–Talagrand: limit exists, \(\approx 0.76\) after normalisation | Existence of the limit is the prize |
| Method | Free energy / Guerra–Talagrand | Deterministic Fourier on the cube + combinatorics of Paley |

Littlewood 4/3 and DMP bound **coefficients of a fixed polynomial** versus \(\|f\|_\infty\). Applied to a 2-homogeneous Walsh polynomial they give Paata’s liminf. They do not identify \(\min_A\).

Talagrand’s typical-\(A\) theorem gives limsup finite because \(\min\le\) typical. It does not force the min to track the typical value, and it does not prove the limit exists.

---

## 7. What this repo uses from the older combinatorics

On Paley orders \(n=p^2+1\):

- Paley conference \(C\), \(CC^\top=qI\), \(n=q+1\).
- Max+ \(=\{y\in\{\pm1\}^n:Cy=py\}\).
- Paley graph of **square order** \(P(p^2)\) as the block graph of the orthogonal array \(\mathrm{OA}(m,p)\) with \(m=(p+1)/2\) quadratic slopes (Goryainov–Lin; Godsil–Royle).
- Canonical cliques = affine \(\mathbb F_p\)-lines in square directions; balanced clique indicators span an eigenspace.

That is the right language for leftover 1 / GLOBAL QVAR. It is **not** a 1933 existence proof of \(\lim\alpha_n\).

---

## 7b. A blind-spot check (2026-08-23): searching the object, not Paata's reading list

§§4–7 trace the problem through **Paata's own path** (SPbU analysis → Volberg
Bellman/hypercube school → Littlewood/Blei/DMP/Talagrand). That search cannot
find results Paata himself never read. This section is a deliberately
*different* search: the exact matrix class \(G(n)\) (symmetric, zero
diagonal, \(\pm1\) off-diagonal — a Seidel matrix) hunted across communities
disconnected from that lineage: classical algebraic graph theory, TCS
(Grothendieck/cut-norm), and statistical physics (worst-case spin glass).

**Found, and not previously in this file:**

> E. Spence, *Eigenvalues of a class of \((0,\pm1)\) symmetric matrices*,
> Linear Algebra Appl. **166** (1992).

Spence studies \(G(n)\) directly and defines
\(g(n):=\min_{C\in G(n)}\lambda_{\max}(C^2)\), the minimum achievable squared
spectral radius. When a conference matrix exists at that order (Cameron–
Delsarte–Goethals' term for the orthogonal case), \(g(n)=n-1\) exactly —
with \(n=p^2+1\) that is \(p^2\), matching this repo's \(Cy=py\) / \(\rho=1\)
setup exactly. When none exists, \(g(n)=n+3\) (a near-miss, not a gap).

**What this is:** independent confirmation, from a completely disconnected
1990s classical-algebraic-graph-theory paper (not cited by, and almost
certainly unknown to, anyone in Paata's analysis lineage), that
\(\min_A(\text{spectral radius of }A)\sim\sqrt n\) holds at **essentially
every** \(n\), not only Paley orders. Combined with the elementary bound
\(\max_{x\in\{\pm1\}^n}|x^\top Ax|\le n\cdot(\text{spectral radius of }A)\),
that reproduces the \(\limsup\le\tfrac12\) direction of the sandwich
genuinely independently of the Paley/\(\rho=1\) route already in this repo.

**What this is NOT:** a resolution of anything. It bounds only the
*continuous* relaxation of the max; it says nothing about whether the
spectral-radius-minimizing matrix's top eigenvector is achievable by a
\(\pm1\) vector (exactly the \(\rho=1\) question this repo already tracks),
and it does not touch E(1), any leftover, or existence of \(\lim\alpha_n\).
Confirms and cross-validates one already-known bound from an independent
source; does not advance the open content. **Do not reopen as progress on
E(1) or the limit.**

Method note for future passes: search the *object* (exact matrix class,
exact quantity) across communities the problem-poser is unlikely to have
read, rather than re-tracing citation lineage from a known starting point —
that is the mechanism, not literature breadth per se.

---

## 7c. GQR codes, inversive circles, and Soto--Andrade operators (2026-08-24)

The W2 investigation exposed a second classical line that was absent from
the original history: generalized quadratic-residue codes of length
(p^2+1).

J. H. van Lint's 1979 paper proves that the two PSL orbits of circles in the
Miquelian inversive plane are exactly the minimum-weight supports in the two
extended GQR codes (Theorem 4.4), and that each orbit spans its code.  In the
repo's notation these are the square/nonsquare circle modules already
identified as (H_0^perp) and (H_0).  This is not merely an analogy: it
names the exact binary code containing all Max-minus differences.

The new note `evidence/NOTE_2026-08-24_w2_gqr_circle_route.md` uses tangent
pencils plus an explicit norm-one-torus reduction to Katz's `t=-2`
Soto--Andrade character estimate to prove that nonsquare
circles meeting the fixed edge ({0,infty}) evenly span
(H_0capker(e_0+e_infty)) for every odd prime.  The nonlinear step also
closes explicitly.  Every subset T of F_p of size (p+1)/2 gives an affine
halfspace +p eigenvector; after nonsquare dilation, every T containing zero
is a completion of one fixed nonsquare circle.  Membership of two transverse
coordinates can be prescribed to put any outside edge in U.  Thus every
eligible circle is a difference inside U, so W1, W2, and Walsh 15.406 E hold
for every odd prime.  The explicit p=19 affine witness supersedes the earlier
HiGHS/CP-SAT UNKNOWN result.  This closes the Walsh slice, not the unrelated
5+-level or other E1 leftovers.

Katz's 1993 estimate gives the uniform (2sqrt p) control for degree-two
Soto--Andrade sums (apart from two explicit exceptional cases).  Kable's
2002 paper supplies the PGL harmonic-analysis and orthogonality dictionary.
Imin Chen's modular-curve paper is a separate useful model: exact
double-coset identities force nonvanishing of Legendre/Soto eigenvalues by
an algebraic-integer ramification argument.  No direct W2 operator identity
has yet been extracted from Chen, so do not cite it as a close.

---

## 7d. The profile-glued eigenspace lattice (2026-08-24)

Proposition 15.629 identifies the saturated integral lattice

\[
L=\ker_{\mathbb Z}(C-pI)
\]

as an odd-index overlattice of the square-circle lattice, with every glue
coordinate supplied by the same polynomial coefficient kernels that control
the Max+ profiles.  It proves

\[
[L:A]=p^{\binom{m-1}{2}},\quad
\det L=2p^{m^2},\quad
L^*=P\mathbb Z^n,\quad
L^*/L\cong\mathbb Z/2\oplus(\mathbb Z/p)^{m^2},\quad
\operatorname{level}(L)=4p.
\]

Böttcher et al.'s rational equiangular-tight-frame lattice is exactly the
nearby dual object (P\mathbb Z^n=L^*).  Their work establishes the correct
ETF/lattice setting but does not identify this integral kernel, its profile
glue, or the odd-coset first-shell moment.  Chapman's conference-matrix
lattices use a different skew-conference/unimodular construction.  Rains--
Sloane shadow theory explains what a modular theta-series route would need,
but it does not apply verbatim: the present even lattice has level (4p),
growing with (p), and R1 concerns a coset shell rather than the ordinary
minimum shell.

OEIS searches of the exact index (25937424601=11^{10}) and the determinant/
index formulas returned generic prime-power occurrences only.  No direct
catalogue entry or prior profile-glue theorem was found.  This is a genuine
structural advance toward R1, not a proof of R1 or of the global limit.

---

## 7e. The exact dual minimum shell (2026-08-25)

Proposition 15.630 proves for every odd prime that

\[
\min(L^*)={1\over2},\qquad
\operatorname{Min}(L^*)=\{\pm Pe_i\},\qquad
|\operatorname{Min}(L^*)|=2(p^2+1).
\]

The proof is arithmetic rather than an extrapolated lattice census.  Circle
pairings obey an exact Parseval identity.  A nonzero common profile sum is
handled by sharp integer balancing.  At common sum zero, the profile-glue
annihilator is a family of projective Reed--Solomon MDS moment codes; Newton
identities turn their support constraints into the sharp positive-mass floor
`(p-1)/2`.

After scaling by `sqrt(2)`, this is the rational ETF lattice studied by
Böttcher--Fukshansky--Garcia--Maharaj--Needell.  Their 2016 paper proved the
exact signed-frame shell computationally for `(5,10)` and `(13,26)` and left
the `(25,50)` shell as `?, ?` in Table 1.  Proposition 15.630 recovers those
finite cases and settles the first (minimum-shell) question for the standard
Paley family at all odd primes; it does not separately settle whether a
minimal-vector subset is a lattice basis.  OEIS
searches found no object entry for the resulting signed-minimum counts.

This does not prove R1: the remaining harmonic coefficient belongs to the
first shell of the odd coset `y0+2L`, not the ordinary minimum shell of `L*`.

---

## 7f. The odd coset has a radial dual phase (2026-08-25)

Proposition 15.631 proves that for every Max+ odd vector `y0 in L` and every
`u=Pz in L*`,

\[
 \langle u,y_0\rangle\equiv 2p\|u\|^2\pmod2.
\]

Consequently the degree-four Poisson transform of `y0+2L` has the scalar
twist `(-1)^(2p||u||^2)`: it has no unidentified glue-class phase and is
independent of the chosen Max+ vector.  Together with Proposition 15.630,
this gives an exact first dual gap from `1/2` to at least `(p-1)/p` and an
explicit positive phased first-shell coefficient
`||W||_F^2/[8(d+2)]` for the R1 harmonic.

This resembles the characteristic-vector/shadow machinery used by Elkies
and by Rains--Sloane, but it is not an application of their unimodular or
strongly modular theorems.  The Paley lattice is non-unimodular with growing
level `4p`, and R1 asks for a degree-four harmonic coefficient of one odd
coset.  Higher dual shells still have uncontrolled signed harmonic mass, so
the result sharpens the exact R1 target without closing it.

---

## 7g. Affine slack parity lifts and split budgets (2026-08-25)

Proposition 15.632 attacks the non-Walsh residual rather than R1.  An odd
candidate edge set produces, in every affine direction, a nonnegative
integer quadratic on the middle Johnson slice.  Its odd-degree graph
boundary fixes the pointwise parity of that quadratic.  Symmetrization over
the odd and even fibres turns the minimum possible mean into an exact
three-variable hypergeometric quadratic-majorant LP.

The directional mean has an additional Paley split: square and nonsquare
directions each receive exactly half of the total slack budget.  Combining
that split with the parity majorant excludes an Eulerian residual boundary
for every odd prime, with gap `(p^2-1)/2`.  A corrected p=5 affine witness
with nonempty line boundary shows why this is a branch kill rather than a
full residual close.

The closest searched literature concerns low-degree functions on Johnson
slices, nonnegative polynomials on finite Boolean sets, and finite-geometry
incidence codes (references 31--33 below).  No searched source gave the
combined edge-boundary parity lift and square/nonsquare budget.  This is a
search record, not a claim of priority.

---

## 7h. Complete second dual shell and the square-circle operator (2026-08-25)

Proposition 15.633 classifies every vector in the next Paley-dual shell.
For `p>=5` it is exactly the disjoint union, with both signs, of projected
signed point-pairs and the signed complements of square Miquelian circles.
Its signed count is

\[
 p(p+1)(p^2+1),
\]

with the exceptional count `30` at `p=3`.  The proof combines the integral
circle profiles from Props. 15.629--15.630, Newton/MDS equality cases, and a
half-conic rigidity lemma whose anisotropic case is forced by Hasse's bound
for a genus-one double cover.

Proposition 15.634 then diagonalizes the associated square-circle
two-secant graph and the projected signed-complement tensor Gram matrix.
The resulting complete degree-four harmonic second shell has three explicit
eigenvalues, all strictly negative for every `p>=11`.  Thus the previously
positive first dual shell is not an anti-cancellation theorem: the very next
shell cancels it in every channel, and a valid R1 proof must control the
later tail or use a genuinely multi-scale theta identity.

Exact OEIS searches for the values
`780,2800,16104,30940,88740,137560,292560,732540` and the polynomial formula
returned no matching sequence. Searches of conference-ETF lattices,
Miquelian-circle designs, modular-lattice shadows, and related literature
found adjacent ingredients but no duplicate of the shell classification or
operator identity. This records the search performed; it is not an
unqualified priority claim for the bare polynomial sequence.

---

## 7i. Third dual norm and exact p=11 shell (2026-08-25)

Proposition 15.635 sharpens the shell gap without claiming an all-prime
third-shell classification.  For every `p>=11`, the third dual norm is
`(p+1)/p`; every nonminimal odd-phase vector has scaled norm at least
`3p-6`.  The signed point-pair orbit at the third norm has a negative scalar
degree-four harmonic operator.  Exact saturated-dual enumeration proves
that this orbit is the complete third shell at `p=11`.

The count `p^2(p^2+1)` is not a new sequence.  It is the elementary family
OEIS A071253 and occurs in A069187; individual searches correctly returned
`14762`, `28730`, and `83810`.  The mathematical content here is the Paley
lattice norm gap, harmonic operator, and exact shell identification, not the
bare integer formula.

---

## 7j. Complete third shell via a finite-field coefficient gap (2026-08-25)

Proposition 15.636 excludes the sole equality profile left by Proposition
15.635. Its positive and negative root polynomials would differ by a
constant and cover all but two field elements, with one repeated root.
Affine normalization turns their product into a geometric-series polynomial.
If that polynomial plus a constant were a square, the reversed square root
would force a long coefficient gap in

\[
(1+cy)^{(p+1)/2}(1-y)^{(p-1)/2}.
\]

High Hasse derivatives evaluated at the two roots force simultaneously
`c=1/3` and `c^2=1/5`, hence `p|4`. Thus the complete third shell is the
signed projected point-pair orbit for every `p>=11`.

This sits near ideal Prouhet--Tarry--Escott configurations, Rédei's theory
of fully reducible lacunary polynomials, and Biró's theorem on low-degree
polynomials taking two values on the multiplicative group. The searched
statements do not directly contain the near-partition multiplicity pattern,
the two Hasse equations, or the Paley-lattice shell result. This is a search
record, not a claim that every underlying polynomial idea is new.

---

## 7k. Zero-common-sum gap at the next even energy (2026-08-25)

Proposition 15.637 applies the same finite-field polynomial viewpoint one
step beyond the complete third shell. At profile energy `p+3`, the MDS mass
bound leaves three possible active-profile counts. The one-profile branch
has only two integer multiplicity patterns. One reduces by Frobenius to two
root polynomials sharing their repeated root. For the other, the reversed
square root `S=sqrt(N/D)` has a coefficient gap of length `(p-3)/2`.
Its first-order differential equation descends through that gap and forces
`N=D`, again identifying repeated and omitted roots.

The two dense active-count branches are then excluded by universal
low-degree moment recurrences for a signed pair. Enough ordinary pair
directions force the relevant binary quartic or cubic to vanish
identically. Their factored defects exclude energy-four profiles and the
doubled energy-six pattern; in the remaining six-unit pattern, four equal
power sums and Newton identities contradict disjointness. Thus the entire
zero-common-sum channel is absent at this candidate energy. Nonzero common
sums remain.

Rédei's coefficient-gap framework remains the closest historical setting.
The specific square-root differential descent and moment-recurrence closure
are elementary and are recorded here as the mechanisms used, without a
broad novelty claim.

---

## 7l. Empty first post-third even candidate shell (2026-08-25)

Proposition 15.638 closes the three nonzero common sums left by Proposition
15.637. Their tiny profile-energy defects reduce every profile to one of a
finite set of integer multiplicity patterns. Binary cubic and quartic
moment recurrences exclude the nondegenerate exceptions, and Newton
identities exclude the remaining unequal root multisets.

The only residual configuration would make a binary quadratic square on
one quadratic-character half of the projective line while having two roots
in that same half. Multiplying by the anisotropic norm form produces a
smooth genus-one character sum of absolute value \(p-3\), contradicting
Hasse's \(2\sqrt p\) bound for \(p\ge11\). Therefore the complete scaled
dual shell \(2(p+3)\) is empty.

Searches located general conference-ETF lattice work and standard
Hasse-bound applications, but not this Paley-dual shell gap. Searching the
candidate values \(28,32,40,44,52,64,68,\ldots\) in OEIS returns unrelated
sequences; no integer-sequence novelty claim is made.

---

## 7m. Complete first nonminimal odd Paley-dual shell (2026-08-25)

Proposition 15.639 classifies the complete first possible nonminimal odd
shell, at scaled norm \(3p-6\), for every prime \(p\ge11\). A parity and
short-coordinate argument reduces every vector to equality in the
common-sum-one profile bound. The sparse equality branch is the known
square-circle family. In the dense branch, divisibility of the degree-three
moment difference by its two-root degree-two difference adds one minimum
vector and reduces the problem to the already complete second and third
shells and the empty shell from Proposition 15.638.

The result is a disjoint union of negative signed conference triangles and
incident point--square-circle vectors, with signed count

\[
 \frac{p^2(p-1)(p+7)(p^2+1)}6.
\]

Exact saturated-dual enumeration at \(p=11\) independently returns the
predicted \(442,860\) vectors at scaled norm \(27\). Targeted searches for
the formula, its \(p=11\) value, and the Paley/conference geometric
description found no matching shell theorem or OEIS entry. The elementary
half-count of negative triangles follows immediately from
\(\operatorname{tr}(C^3)=0\); novelty is not claimed for that identity or
for the bare polynomial sequence.

---

## 7n. Through-point circle frame and quartic saddle (2026-08-25)

Proposition 15.640 computes the complete degree-four harmonic operator of
the shell classified in Proposition 15.639. The key new incidence identity
is that signed square-circle complements through a fixed point form a tight
frame for the codimension-one slice of the positive Paley eigenspace. It
follows directly by splitting the circles through infinity into orthogonal
parallel classes with Gram matrix \(p(pI-J)\).

The negative-triangle family contributes a scalar fourth-moment operator.
The point--circle family contributes another scalar plus the same
square-circle tensor operator diagonalized in Proposition 15.634. Their
complete harmonic sum is indefinite for every \(p\ge11\): one
circle-kernel eigenvalue is negative and two circle-image eigenvalues are
positive. The radial odd-coset phase reverses all three signs in the Poisson
shadow.

Exact searches for the three rational eigenvalue formulas, the \(p=11\)
spectrum `-582/7,258/7,426/7`, and the through-point frame identity found no
matching Paley-lattice theorem. Individual OEIS searches for larger raw
numerators likewise returned no entries. General harmonic ETF and
Miquelian-plane sources provide context, not this shell operator. This is a
search record and does not assert unreviewed priority.

---

## 7o. Mod-seven closure of the p=7 unsaturated four-point branch (2026-08-26)

Proposition 15.655 turns the complete Johnson slack catalogs into right
sides of one common finite-field system. The 280 exact affine score rows,
edge count, and distinguished edge form a `282 x 1225` matrix of rank 147
over `F_7`. Its 135 left-null dependencies reject all 1,716,742,440 catalog
tuples in the 2,408 fixed elevation cases covering the 518 unsaturated
boundary orbits. A separate implementation reconstructs the matrix,
dependencies, target coefficients, and coverage and also obtains zero
survivors. Together with the nonsquare anti-isometry of Proposition 15.654,
this closes both product signs and hence every `p=7` size-four boundary.

The closest literature found after the result is the theory of finite-field
codes generated by Paley graph incidence matrices, including
Ghinelli--Key (2011), and the broader Johnson/finite-incidence-code context.
Those papers concern different incidence matrices and do not state this
affine-score rank, catalog syndrome system, or boundary exclusion. Exact
searches for the rank description produced no matching theorem. Direct
OEIS searches returned no result for the total `1716742440`, the largest
single-case product `3939012`, the five pattern counts, or the rank tuple.
These are duplicate/context checks, not an integer-sequence novelty claim.

---

## 7p. Full-shell closure of the p=5 four-point branch (2026-08-26)

Proposition 15.656 uses the complete `p=5` eigenshell rather than only its
affine subfamily. After antipodal quotient there are 130 score rows, each
normalized edge column sums to 26, and every putative 21-edge graph has
total shell slack 78. Boundary/product parity reduces each slack to a fixed
bit plus twice a bounded lift. The resulting `132 x 325` shell matrix has
rank 67 over both `F_5` and `F_7` and 65 left dependencies.

Complete finite scans close 712 orbit cases modulo five. The sole timeout is
independently infeasible modulo seven. An exact nonsquare anti-isometry and
a fresh 489-orbit bijection transfer the no-infinity negative-sign result to
the positive sign. Thus all 1,202 floor-surviving orbit/sign cases, covering
26,450 boundary/sign cases, are excluded. With Proposition 15.632 this
closes every `p=5` size-four boundary; with 15.652--15.655, every size-four
boundary is closed for every odd `p>=5`. Boundary size at least six remains.

Targeted searches found ordinary Paley incidence-code ranks and Johnson
scheme methods, but not this full-eigenshell matrix or bounded-syndrome
exclusion. Individual counts occur in unrelated OEIS sequences; the rank
and coverage tuples did not match a relevant entry. No sequence novelty
claim is made.

The requested read of Ivanisvili--Stolyarov--Vasyunin--Zatitskii,
arXiv:2305.03523, found a potentially useful Bellman-function analogy for
compressing R1/QVAR to two moments, but no theorem directly applicable to
the finite Paley residual system. It did not enter Proposition 15.656.

---

## 7q. Pair-deficit closure of every six-point boundary for `p>=11` (2026-08-26)

Proposition 15.657 extends the exact positive quadrature of Proposition
15.652 through six odd fibres. If a projective direction has finite
boundary-fibre multiplicities `n_i` and `b_d` odd fibres, then
`s-b_d <= 2 sum_i binom(n_i,2)`. Every finite pair belongs to exactly one
projective direction, so `sum_d(s-b_d) <= s(s-1)`: the deficit budget is 30
for six finite points and 20 for infinity plus five finite points.

For infinity plus five points, the resulting phase-independent floor bound
exceeds the total affine slack budget by `p^2-9p+10`, positive from `p=11`.
For six finite points and `p>=13`, the excess is `p^2-12p+7`, positive from
`p=13`. At `p=11`, the two quadratic types have opposite phases and their
separate budgets require deficits at least 20 and 18, contradicting 30.
Thus every six-point boundary is closed for every odd prime `p>=11`;
Propositions 15.658--15.661 subsequently close both infinity-present signs
and six finite points at `p=7`, as well as every `p=5` size-six case.

Targeted searches of Johnson-slice, finite-incidence, and Paley
degree-parity sources found adjacent methods but no duplicate of the
directional slack/pair-deficit/type-budget argument. Li--Zhou studies
parity-uniform induced Paley subgraphs and MDS-code applications, not the
edge-set boundaries or affine slack system here. The result proposes no new
integer sequence, so no OEIS submission search was made. This is a
duplicate/context search record, not an unqualified priority claim.

---

## 7r. Complete `p=5` size-six exclusion (2026-08-27)

Proposition 15.660 rebuilds the four exact `p=5` size-six boundary catalogs
for both product signs and infinity bits. The no-infinity catalogs contain
159,050 survivors and 6,766 square-semilinear orbits per sign, exchanged by
a nonsquare multiplier. Complete coarse full-shell SCIP batches leave seven
classes before cross-infinity symmetry; infinity orbit 1144 maps to
no-infinity orbit 881, leaving exactly six signed-symmetry classes.

Independent layered audits close classes `0`, `881`, `2529`, `3032`,
`4731`, and `4939`. They reconstruct every pattern, degree, and recursive
crossing-edge quotient and verify every recorded SCIP infeasibility status.
The from-definitions global audit is true, and all 97 class-881 artifact
hashes were rechecked in the permanent `/mnt/storage/` archive. Thus every
`p=5` size-six boundary is excluded. Proposition 15.661 subsequently closes
the six-finite `p=7` branch; larger boundaries and the top-level theorem
remain open.

Targeted arXiv searches found no duplicate Paley residual classification.
The count 159,050 occurs in OEIS A063533 for an unrelated construction from
twin-prime Pythagorean triples, so the match is coincidental and supplies no
structural input or sequence claim.

---

## 7s. Complete `p=7` six-finite exclusion (2026-08-27)

Proposition 15.661 evaluates the exact six-finite floor budget on all
`C(49,6)=13,983,816` boundaries. It leaves 3,856,300 boundaries and 80,704
square-semilinear orbits. For 80,519 ordinary orbits, every type-floor sum
is 24 or 32, so at most one non-singleton catalog occurs per type.
Simultaneously joined left-null signatures over `F_3` and `F_7` reject all
160,745 exact elevation cases.

The other 185 orbits require larger mean allocations. Compact exact models
retain all 35 slack values per direction, their parity, the 14 primitive
integer degree-two relations, the exact type means, and both modular edge
systems. They close 92 orbits immediately. The 93 timeouts split into all
930 exact mean allocations: 810 close directly, and the remaining 120
low-catalog leaves close by exact two- or three-table signature joins.

An independent NUKA NumPy sweep reproduces the V100 survivor hash and full
histogram. Its serial quotient reproduces the ordered orbit catalog and
profile histogram canonically, and its serial ordinary sweep reproduces all
160,745 rejections. A nonsquare signed anti-isometry transfers the opposite
product sign. Together with Propositions 15.657--15.660, every size-six
boundary is therefore closed for every odd prime `p>=5`. Boundaries of size
at least eight and the top-level analytic remainders remain open.

Targeted literature searches found adjacent generalized-Paley and
Johnson-slice methods, but no duplicate of this finite residual
classification. Exact OEIS searches find 80,704 in A060716, A133751, and
A133756, and 160,745 in A254067, all for unrelated constructions;
3,856,300 has no exact hit. No sequence claim is made from these counts.

---

## 7t. Exact p=11 profile and quartic-trace reconstruction (2026-08-28)

Proposition 15.667 reduces all \(11^{10}\) words of the p=11 glue-dual code
to 21,437,340 translation/scalar representatives and 2,558,543 weighted
six-profile tuples. Exact one-profile dynamic programs collapse 604 quartic
value distributions to 13 affine types. Five 31-bit prime moduli admitting
primitive eleventh roots have product
larger than unrestricted integer bounds for the ordinary coefficient and
the common-sum second and fourth moments through scaled norm 120, so CRT
recovers all three sequences exactly. Every one of the 51 nonempty shells
passes the independent tight-frame second-moment identity.

The common-sum fourth moment gives the trace of the positive raw quartic
shell operator directly. This exposed a normalization error in the first
version of Proposition 15.665: the scaled-norm 20 and 24 rows used
\(H(u/2)\) from the Poisson-shadow calculations as though it were \(H(u)\).
Multiplication by the required quartic factor 16 changes their raw traces to
\(89792/11\) and \(7076\). The general positive-operator and conservation
theorem is unchanged.

The exact ordinary prefix through exponent 88 has full rank 41 and uniquely
determines its scalar modular form. The exact quartic-trace prefix first has
full affine rank 32 at exponent 92. The resulting reconstructions through
800 predict all 32 and 28 held-out profile coefficients through 120,
respectively, and every reconstructed ordinary coefficient and raw trace has
the required nonnegative integrality or rational positivity.

Exact QSopt_ex primal/dual certificates then impose the reconstructed raw
trace as a conserved mass at every shell through 800. Seven of eight
component target endpoints are already identical at truncations 120 and
800; the circle-low-Weil maximum contracts from about 880.0044 to 874.9202.
All four final intervals remain broad and two-sided. Thus scalar
trace-conservation is a genuine tightening but does not close R1. The next
live refinement is to compute channel-resolved moments
\(\operatorname{tr}(R_sO)\) and \(\operatorname{tr}(R_sO^2)\) for the
square-circle tensor operator, or finer PSL-twisted traces.

Targeted searches for the corrected rational traces, the exact half-cusp
trace numerator, profile/common-coordinate fourth-moment terminology, and
Paley quartic theta combinations found only adjacent weighted-theta,
shell-design, and complete-weight-enumerator literature. Individual OEIS
searches of five large reconstructed coefficients returned no entries. This
is a duplicate/context check and makes no sequence or unreviewed priority
claim.

---

## 8. What is not in the older literature

- No 1930s–1970s theorem that \(\lim n^{-3/2}\min_A\max_x|x^\top Ax|\) exists.
- Sandwich (liminf \(\ge 2^{-5/2}\) or \(1/\pi\), limsup \(\le 1/2\) or Talagrand’s typical constant) is exactly what Paata already said does **not** settle existence.
- Bowlin bipartition, Spencer discrepancy, Seidel switching, and Paley eigenfunctions are nearby objects. None of them is the prize.
- Paata’s own cube papers (Poincaré 3/2, square functions, Jackson on the hypercube, KKL, hypercontractivity) sit in the right neighbourhood. None of them is this limit.

Do not reopen as existence proofs: BH / DMP as a substitute for E(1); Talagrand typical-\(A\) as \(\lim=c\); Bowlin; “Paley \(\Rightarrow\lim=1/2\)” without E(1).

---

## 9. Bibliography (primary)

1. J. E. Littlewood, *On bounded bilinear forms in an infinite number of variables*, Quart. J. Math. (Oxford) **1** (1930), 164–174.
2. R. E. A. C. Paley, *On Orthogonal Matrices*, J. Math. Phys. **12** (1933), 311–320.
3. R. E. A. C. Paley and A. Zygmund, *A note on analytic functions in the unit circle*, Proc. Camb. Phil. Soc. **28** (1932) (Paley–Zygmund; random series).
4. A. Khintchine, *Über dyadische Brüche*, Math. Z. **18** (1923), 109–116.
5. V. Belevitch, conference-matrix papers, 1950s (telephony); see Goethals–Seidel for the citation trail.
6. J. H. van Lint and J. J. Seidel, equiangular lines / \(C\)-matrices (1966).
7. J. M. Goethals and J. J. Seidel, *Orthogonal matrices with zero diagonal*, Canad. J. Math. **19** (1967), 1001–1010. [Cambridge PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/28A22B5B4DBB2C7DA9DF99DD89013A7A/S0008414X00054997a.pdf/div-class-title-orthogonal-matrices-with-zero-diagonal-div.pdf).
8. P. Delsarte, J. M. Goethals and J. J. Seidel, *Orthogonal matrices with zero diagonal II*, Canad. J. Math. **23** (1971), 816–832.
9. R. Blei, *Analysis in Integer and Fractional Dimensions*, Cambridge Univ. Press, 2001, Theorem 34.
10. M. Talagrand, *The Parisi formula*, Ann. of Math. **163** (2006), 221–263.
11. A. Defant, M. Mastyło, A. Pérez, *On the Fourier spectrum of functions on Boolean cubes*, [arXiv:1706.03670](https://arxiv.org/abs/1706.03670) (2017).
12. G. A. Jones, *Paley and the Paley graphs*, [arXiv:1702.00285](https://arxiv.org/abs/1702.00285) (2017).
13. S. Goryainov and H. Lin, *On balanced characteristic functions of canonical cliques in Paley graphs of square order*, [arXiv:2104.08839](https://arxiv.org/abs/2104.08839) (2021).
14. S. Goryainov, L. Shalaginov, C. H. Yip, *On eigenfunctions and maximal cliques of generalised Paley graphs of square order*, [arXiv:2203.16081](https://arxiv.org/abs/2203.16081) (2022).
15. P. Ivanisvili, *A little bit of Fourier analysis*, 6 Dec 2019, https://extremal010101.wordpress.com/2019/12/06/a-little-bit-of-fourier-analysis/
16. P. Ivanisvili, MathOverflow 413935 (16 Jan 2022), https://mathoverflow.net/questions/413935
17. P. Ivanisvili, X, 25 Jul 2026, https://x.com/PI010101/status/2081070728422752329
18. G. Bowlin, *Maximum frustration of bipartite signed graphs*, Electron. J. Combin. **19** (2012), #P10 (Starr pointer; not an existence proof).
19. A. Volberg, *Bellman function technique in Harmonic Analysis*, [arXiv:1106.3899](https://arxiv.org/abs/1106.3899) (2011).
20. P. Ivanisvili, CV (Princeton, c. 2018), https://web.math.princeton.edu/~paatai/CV-Paata.pdf
21. E. Spence, *Eigenvalues of a class of \((0,\pm1)\) symmetric matrices*, Linear Algebra Appl. **166** (1992), https://www.sciencedirect.com/science/article/pii/002437959290280N (§7b; found by blind object-search, cross-validates \(\limsup\le\tfrac12\) independently, does not advance the open content).
22. J. H. van Lint, *Generalized Quadratic-Residue Codes*, IEEE Trans. Inform. Theory / Eindhoven report (1979), https://pure.tue.nl/ws/files/4233047/593589.pdf (Theorem 4.4: Miquelian circles are the minimum-word supports of the extended GQR codes).
23. N. M. Katz, *Estimates for Soto-Andrade sums*, J. reine angew. Math. **438** (1993), 143--162, https://web.math.princeton.edu/~nmk/old/sotosums.pdf
24. A. C. Kable, *Legendre sums, Soto-Andrade sums, and Kloosterman sums*, Pacific J. Math. **206** (2002), 139--157, https://msp.org/pjm/2002/206-1/pjm-v206-n1-p09-s.pdf
25. I. Chen, *On relations between Jacobians of certain modular curves*, [arXiv:math/9809209](https://arxiv.org/abs/math/9809209) (double-coset/nonvanishing model only; not a W2 proof).
26. A. Böttcher, L. Fukshansky, S. R. Garcia, H. Maharaj, D. Needell, *Lattices from tight equiangular frames*, Linear Algebra Appl. **510** (2016), 395--420, https://www1.cmc.edu/pages/faculty/lenny/papers/lattices_frames.pdf
27. R. Chapman, *Conference matrices and unimodular lattices*, [arXiv:math/0007116](https://arxiv.org/abs/math/0007116) (nearby conference-lattice construction; not the integral (+p)-kernel here).
28. E. M. Rains and N. J. A. Sloane, *The Shadow Theory of Modular and Unimodular Lattices*, J. Number Theory **73** (1998), 359--389, [arXiv:math/0207294](https://arxiv.org/abs/math/0207294) (methodological comparison; its strongly modular hypotheses do not close R1).
29. N. D. Elkies, *A characterization of the `Z^n` lattice*, Math. Res. Lett. **2** (1995), 321--326, [arXiv:math/9906019](https://arxiv.org/abs/math/9906019) (characteristic-vector/shadow precedent for unimodular lattices; not directly applicable to the level-`4p` Paley lattice).
30. N. D. Elkies, *Lattices and codes with long shadows*, Math. Res. Lett. **7** (2000), 151--164, [arXiv:math/9906086](https://arxiv.org/abs/math/9906086) (modular-form use of odd cosets in the unimodular setting; methodological comparison only).
31. Y. Filmus, G. Kindler, E. Mossel, and K. Wimmer, *Invariance principle on the slice*, [arXiv:1504.01689](https://arxiv.org/abs/1504.01689) (Johnson-slice low-degree framework; no parity-budget theorem).
32. G. Blekherman, J. Gouveia, and J. Pfeiffer, *Sums of squares on the hypercube*, [arXiv:1402.4199](https://arxiv.org/abs/1402.4199) (nonnegative polynomial framework; adjacent, not the exact hypergeometric LP used here).
33. P. Sin, J. Sorci, and Q. Xiang, *Linear representations of finite geometries and associated LDPC codes*, [arXiv:1908.06824](https://arxiv.org/abs/1908.06824) (finite incidence-code context).
34. L. Rédei, *Lacunary Polynomials over Finite Fields*, North-Holland, 1973 (fully reducible coefficient-gap framework; adjacent to 15.636).
35. A. Biró, *On Polynomials over Prime Fields Taking Only Two Values on the Multiplicative Group*, Finite Fields Appl. **6** (2000), 302--308, https://www.renyi.hu/~biroand/pdfs/TwoValues.pdf (two-value degree bounds; adjacent, not a direct proof of 15.636).
36. D. Ghinelli and J. D. Key, *Codes from incidence matrices and line graphs of Paley graphs*, Adv. Math. Commun. **5** (2011), 93--108, https://doi.org/10.3934/amc.2011.5.93 (Paley incidence-code ranks; nearby finite-field coding context, not the 15.655 affine-score system).
37. P. Ivanisvili, D. Stolyarov, V. Vasyunin, P. Zatitskii, *Bellman functions on simple non-convex domains in the plane*, [arXiv:2305.03523](https://arxiv.org/abs/2305.03523) (minimal locally concave Bellman construction; possible R1/QVAR analogy, not a Paley residual theorem).
38. Q. Li and Y. Zhou, *On induced subgraphs with degree parity conditions in Paley graphs and Paley tournaments*, [arXiv:2512.19312](https://arxiv.org/abs/2512.19312) (degree-parity context for induced Paley subgraphs and MDS codes; not the directional affine-slack boundary problem of 15.657).

Secondary: Paley construction (Wikipedia, quoting the 1933 Hadamard-conjecture sentence); conference matrix (Wikipedia, van Lint–Seidel sum-of-two-squares obstruction).
