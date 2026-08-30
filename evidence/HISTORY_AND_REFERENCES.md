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
trace-conservation is a genuine tightening but does not close R1. Proposition
15.668 subsequently computes the channel-resolved contractions
\(\operatorname{tr}(R_sO)\) and \(\operatorname{tr}(R_sO^2)\); even that
strict refinement does not close the cone, leaving finer PSL-twisted or
uniform transport inequalities as the live route.

Targeted searches for the corrected rational traces, the exact half-cusp
trace numerator, profile/common-coordinate fourth-moment terminology, and
Paley quartic theta combinations found only adjacent weighted-theta,
shell-design, and complete-weight-enumerator literature. Individual OEIS
searches of five large reconstructed coefficients returned no entries. This
is a duplicate/context check and makes no sequence or unreviewed priority
claim.

---

## 7u. Exact p=11 broad-channel reconstruction and finite R1 (2026-08-28)

Proposition 15.668 retains the marked Legendre-convolution statistic

\[
 U_4(a)=\sum_c\left(\sum_s\eta(s-c)a_s\right)^4
\]

in the quartic glue-profile census. Input-affine reduction leaves 1,007
profile types and output-affine reduction leaves 20 canonical dynamic
programs. The complete glue code reduces from \(11^{10}\) words to
21,437,340 translation/nonzero-scalar representatives and 2,584,901 weighted
rich-profile tuples. Five-modulus CRT is exact because its modulus product
exceeds independent unrestricted bounds for every recovered integer through
exponent 120.

The scalar trace, \(z^Tz\), and \(z^TA_2z\) recover the raw shell mass on the
square-circle kernel, low, and high spaces, of dimensions 1220, 305, and 244.
All three affine modular spaces have rank 32 by exponent 92. Solving on that
prefix and holding out exponents 93--120 reproduces all 28 withheld
coefficients in each channel exactly. Their continuations through exponent
800 are nonnegative, sum coefficientwise to the aggregate trace, and satisfy
the exact dimension-weighted transformed-target identity.

Eight separate QSopt_ex endpoint problems impose channelwise shell-mass and
target conservation. Every rational primal constraint and dual stationarity
equation was independently checked. The certified target intervals are much
narrower than the aggregate intervals from 15.667, but each still contains a
value mapping to \(\Phi<6\). Thus broad square-circle channel conservation is
an exact failed proof route; it is not a counterexample to R1.

The independent full \(p=11\) Max+ census gives

\[
 \|\delta\|^2={1382747375360\over583792784981}
 < {61\over6}={n\over12}
 < {22143\over1682},
\]

with strong margin
\(27314875631681/3502756709886>0\). Strong R1, and hence R1, is therefore an
exact finite theorem at \(p=11\). This does not supply an all-prime proof.
General R1, global QVAR, the non-Walsh remainder, Type I, and the limit remain
open.

Tasaka's weighted-theta survey and Ozeki's association-scheme/Siegel-theta
work confirm nearby standard machinery but contain no channel transport
inequality for these Paley shells. Targeted searches of the large exact
values and dimension triple found no relevant OEIS match; no novelty claim is
drawn from that negative search. The complete 33-file archive is under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1-broad-channel/`
and its manifest hash is
`d1ef69b9af7007c0d2f09a3a5ea8a014cde62d9ed6109175cf4a6496d06b3f07`.

---

## 7v. Uniform non-Walsh boundary-range exclusion (2026-08-28)

Proposition 15.669 evaluates the exact hypergeometric parity-majorant floor
through the entire middle range. For every odd \(p\ge17\), every
\(5\le b\le p-5\), and either parity phase, the constant quadratic one is
optimal. The dual certificate is constructive: complement to
\(b\le(p-1)/2\), bracket the hypergeometric mean by adjacent parity-one
contacts for the minimum variance, use the two endpoint contacts for the
maximum variance, and mix the two rational measures to match the exact
variance. The four upper and four lower envelope margins are positive by
elementary factorizations and two concave-quadratic endpoint checks.

At residual size \(4p+1\), combining these floors with the separate
quadratic-type budgets and

\[
 \sum_d(s-b_d)\le s(s-1)
\]

gives a uniform geometric range theorem. Every all-finite even boundary with
\(6\le s\le3(p-1)/4\) is impossible, because four times the deficit gap is
\(s(3p+3-4s)-4p>0\). Every infinity-present boundary with odd finite
\(5\le s\le p-4\) is also impossible. The infinity proof includes the
complemented endpoint exception: in phase one and \(p\equiv1\pmod4\),
\(b=p-4\) can save six units, but allowing two such non-\(b=1\) directions
still misses the type saving target by \(p-11>0\).

Exact rational count-profile programs extend the exclusion to infinity plus
seven points at \(p=11\), and at \(p=13\) to eight finite points and infinity
plus seven or nine. The first profiles beyond the proved ranges survive the
floor-and-pair relaxation. They are not affine boundaries or residual graphs,
so residual (ii), R1, global QVAR, Type I, and the limit remain open.

Targeted duplicate searches found adjacent work on low-degree analysis of
Johnson slices, nonnegative polynomials on the cube, finite incidence codes,
and the non-binary Johnson scheme of Bernard--Crampé--Vinet [41], but no
source stating the full-middle parity quadrature or its Paley type-budget
consequence. Individual OEIS searches for larger endpoint values found only
unrelated sequences (and no entry for 750002). This is context checking, not
an unqualified priority or sequence claim.

---

## 7w. Finite p=11 size-eight boundary exclusion (2026-08-28)

Proposition 15.670 closes the first finite survivor of Proposition 15.669.
Every eight-point set in \(AG(2,11)\) has an affine-similarity image
containing field points zero and one. The exact pointed-set identity

\[
 \binom{121}{8}\,8\cdot7
 =\binom{119}{6}\,121\cdot120
\]

reduces exclusion to all \(3,470,108,187\) normalized sets. Multiplication
by a nonsquare swaps the two quadratic direction types and simultaneously
transfers the boundary phase sign; testing both signs is therefore lossless.
The source checks this action on all scalars, translations, and directions.

Exact direct-rank CUDA and HIP censuses independently test every normalized
set for both signs. Both complete cost-pair histograms agree, both signs have
zero survivors under the exact type budget 72, and the exact minimum larger
type cost is 76. A separate CPU combinations traversal matches every
histogram entry in a 100,000-set prefix. Thus every finite \(p=11\)
size-eight boundary is impossible. Infinity plus nine and larger boundaries
remain; residual (ii) is not closed.

The adjacent literature found in the post-result search studies odd secants
of \(q+2\)-point sets [42] and prescribed-direction integral point sets with
Paley connections [43], not this eight-point split-floor exclusion. OEIS
contains the two raw counts in unrelated binomial tables: 3,470,108,187 in
A004379 and 899,749,479,915 in A126450. Those matches are context only and
carry no sequence or priority claim.

The permanent machine record is archived under
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-p11-size-eight-boundary/`.

## 7x. Rigid-sign near-line exclusion (2026-08-28)

Proposition 15.671 attacks an actual realization of the first general
infinity-present floor survivor: infinity together with `p-2` collinear
finite points. In one product-sign branch, the `b=1` line floor and the
complementary `b=p-2` transverse floors use the complete split budget and
force every quadratic slack pointwise. The resulting targets are linear in
the line direction and a single two-fibre xor in every transverse direction.

Coefficient comparison gives parallel-count congruences whose sum is
`I=3+sigma mod (p-1)/2`, where `I` is the infinity-edge count. For
`p=1 mod 4,c_H=-1`, this contradicts the required odd parity of `I`. For
`p=3 mod 4,c_H=+1`, the signed inter-fibre `l1` capacities and total edge
count force an even quotient that must be odd, once `p>=19`. Thus one sign
of every collinear infinity-plus-`(p-2)` branch is excluded uniformly. The
opposite sign and noncollinear boundaries remain open. Post-result searches
found no matching Paley odd-fibre/coefficient theorem; the prime lists are
only the standard modulo-four prime sequences (OEIS A002144 and A045326),
not a new sequence.

## 7y. Complete collinear near-line exclusion (2026-08-28)

Proposition 15.672 handles the sign opposite Proposition 15.671. The exact
directional mean `a_d=I+(p+1)P_d-eps_d*T-3p` quantizes same-type excesses in
units of `p+1`. The split budgets and the degree-two lift support floor leave
exactly one exceptional direction per type. Writing `x,y` for the baseline
parallel counts gives `I=4p-1-(p+1)(x+y)/2` and `x+y<=7`.

Every type retains a transverse xnor baseline. Its coefficient congruence,
after substituting the global edge count, says `(p-1)/2` divides both `x+1`
and `y+1`. This forces `x+y>=p-3>=8` in the applicable range, a
contradiction. The opposite sign is therefore empty from `p=11` or `p=13`
according to the residue class. Combined with Proposition 15.671, both signs
of every collinear infinity-plus-`(p-2)` boundary are excluded for every
prime `p>=13`. Noncollinear boundaries remain open.

## 7z. Complete endpoint-only near-line exclusion (2026-08-28)

Proposition 15.673 removes collinearity while retaining the endpoint
hypothesis `b_d in {1,p-2}`. Exact same-type means are congruent modulo
`p+1`; together with the four-unit minimum lift cost, this leaves only two
baseline parallel counts and at most one exception per type. The endpoint
pair-deficit equality case is an arc with exactly three undetermined
directions. Adjoining any two gives a `p`-arc, and Segre's odd-order `p`-arc
theorem [44] puts it on a conic. Two choices force one conic through all
three collinear infinity points, a contradiction.

The remaining coefficient arithmetic splits into four residue/phase rows.
Three are excluded by `(p-1)/2` divisibility and the boundary support bound
`I<=p-2+2E`. The only endpoint is `p=17,(x,y)=(0,7)`, where a complementary
baseline direction prescribes an inter-fibre matrix of exact `l1` norm at
least 75 but has only 57 transverse edges. Thus every endpoint-only
infinity-plus-`(p-2)` boundary is excluded for both signs and every prime
`p>=17`. Non-endpoint profiles remain open. Literature and OEIS searches
found no duplicate of the combined Paley mean/coefficient argument and no
relevant integer-sequence interpretation.

## 7aa. Complete infinity-plus-`(p-2)` shell exclusion (2026-08-28)

Proposition 15.674 removes 15.673's endpoint hypothesis. Across every odd
fibre count, all directional floors are at least `p-1`, and every
intermediate floor is strictly above `p+1`. The exact same-type sum then
permits only residues zero and `p-1`: an intermediate count can occur only
as the unique mean-`2p` exception of its type. Pair deficit excludes two
`b=1` baseline types, while two complementary baseline types determine at
most two directions and are collinear. The forced mixed pair has exactly
15.673's four coefficient-congruence rows and exact `p=17` norm obstruction.
Thus the full infinity-plus-`(p-2)` shell is excluded for both signs and all
primes `p>=17`; the next shell has `p` finite points.

The determined-direction step is elementary but adjacent to the Rédei
direction literature, including Lev [45]. Targeted literature searches found
no source combining those ideas with Paley type sums and parity-majorant
floors. OEIS searches for the two arithmetic-candidate value lists and the
pair `75,57` produced unrelated matches only; no sequence claim is made.

## 7ab. Quantized first all-finite survivor half-close (2026-08-28)

Proposition 15.675 returns to the first even all-finite size above
`3(p-1)/4` and adds the exact same-type mean residue omitted by 15.669's
floor-only optimization. Phase one becomes rigid; phase zero has unique
minimum residue four. The exact deficit gaps in the four classes
`p=1,3,5,7 mod 8` are

```text
-(p-1)/4, (p+1)/2, (p-1)/2, -(p-7)/4.
```

The positive middle rows close the first survivor for every prime `p>=19`
congruent to `3` or `5 mod 8`; the negative outer rows remain open. An
independent NUKA replay agrees with the local symbolic and dynamic-program
ledgers. Ball--Csajbók [42] is adjacent odd-secant literature but does not
state this Paley type-residue optimization. OEIS searches returned unrelated
matches to the linear gap values and no sequence interpretation is claimed.

## 7ac. Infinity-plus-`p` arc-equality close (2026-08-28)

Proposition 15.676 treats the pair-deficit equality branch of the next
infinity-present shell. Equality makes the `p` finite points a `p`-arc, so
Segre's theorem [44] puts them on a conic. Classifying the line at infinity
leaves two affine profiles: a tangent conic has `p` directions with `b=1`
and one with `b=p`; an external conic with one point deleted has `m+1`
directions with `b=1` and `m-1` with `b=3`.

The external profile exceeds the exact Paley type-floor capacity in both
phases. For the tangent profile, same-type mean quantization retains a
baseline `b=1` direction in each type; its coefficient congruence forces
both baseline parallel counts to vanish, leaving zero or two finite edges
and violating boundary support. Thus pair-deficit equality is closed for
both signs and all primes `p>=17`. Strict deficit and the full shell remain
open. Van de Voorde [46] is adjacent external-conic/determined-direction
literature but does not supply the Paley floor or coefficient argument.
Targeted literature and OEIS checks found no duplicate or relevant sequence
interpretation.

## 7ad. First all-finite survivor closed from `p=19` (2026-08-28)

Proposition 15.677 resolves the two modulo-eight classes left by 15.675.
For `p=1,7 mod 8` from `p=23`, exact quotient/deficit arithmetic leaves
phase-zero residue `u_0=2`, and additionally `u_0=3` in the first class.
Because `sum k_d=m-u_0<m`, one phase-zero direction has quotient zero. Its
mean is four or six, below the phase-zero floor at every nonzero even fibre
count, so it has `b=0` and pointwise form `A_d=2B_d` with `B_d` a nonzero
nonnegative integer-valued quadratic. Proposition 15.642's exact slice
support bound gives `4p E[B_d]>=8`, contradicting that mean.

Together with 15.675, this closes the first even all-finite size above
`3(p-1)/4` for every prime `p>=19`. The exact audit caught and retains a
smaller exception: at `p=17`, residue `u_0=0` also survives the uniform
reduction, so that endpoint is not claimed. Amireddy--Behera--Srinivasan--
Sudan [47] supplies the degree-two slice-distance lemma already certified in
15.642; targeted searches found no source combining it with the Paley
same-type residue and zero-quotient argument. OEIS searches for the lift-cost
samples returned unrelated entries, and no sequence claim is made.

## 7ae. Exceptional `p=17` first-survivor close (2026-08-28)

Proposition 15.678 resolves the additional `u_0=0` row left explicitly by
15.677. The exact `p=17,s=14` ledger first excludes `u_0=2` by the six-unit
lift floor and `u_0=3` by the coefficient `l1` bound. Pair-slack
divisibility then leaves two `u_0=0` profiles. Both attain pair equality and
are 14-arcs with six seven-secants, eight six-secants, one one-secant, and
three undetermined directions.

Adjoining two undetermined infinity points gives a 16-arc. Sticker's full
classification [48] records exactly one PGL-equivalence class of 16-arcs in
`PG(2,17)`. Since conic-minus-two is one such arc, every 16-arc is
conic-contained. The third undetermined infinity point is off that conic;
an external conic point has eight conic secants and an internal point nine,
so deleting the four points outside the original 14-set leaves at least four
secants through it. This contradiction closes the endpoint. Keri [49]
provides the earlier MDS/superregular classification in the same small-field
range; Sticker reports agreement.

Together with 15.675/15.677, the first all-finite survivor is now excluded
for every prime `p>=17`. Later all-finite sizes remain open. An exact OEIS
search for `17633,21064,6814,629`, the distinctive larger block of the
`PG(2,17)` arc-class row, returned no result; this is a duplicate/context
check, not a sequence or priority claim.

## 7af. Next all-finite boundary closed from `p=43` (2026-08-28)

Proposition 15.679 advances one full size beyond the first survivor. At the
second even size above `3(p-1)/4`, phase one is rigid and exact phase-zero
quotient arithmetic leaves only residues `2<=u<=7`. Every retained row
forces a quotient-zero `b=0` direction with mean at most 14 and pointwise
form `A=2B`. The degree-two slice-distance floor [47] exceeds 14 from
`p=59`; exact pair/lift ledgers close the only smaller in-scope primes
`43,47,53`.

Thus this boundary is impossible for every prime `p>=43`. This proposition
leaves the same size at `p=17,19,23,29,31,37,41` open; Propositions
15.680--15.683 subsequently close `p=37,29,31,41`. Ball--Csajbók [42] is adjacent
odd-secant literature but has a different `q+2`-point scope. The newly
posted low-degree testing paper [50] concerns robust testing rather than an
exact support floor or affine type budget. Exact OEIS searches for
`11130,7176,3922` and `21756,14016,7696` returned no result; no sequence or
priority claim is made.

## 7ag. The `p=37,s=30` next-boundary endpoint closes (2026-08-28)

Proposition 15.680 resolves the only row where Proposition 15.642's
degree-two lift floor was attained exactly. The exact pair ledger at
`p=37,s=30` leaves `u=2,3,4,5`; all force a quotient-zero `b=0`
direction. The old lift floor excludes the first three.

At `u=5`, the lift has mean `5/74`. Stabilizer averaging forces its values
into `{0,1,2}`. The degree-two distance floor bounds the density of value
two by `2/1295`, while applying the degree-four instance of [47] to
`B(B-1)` requires at least `1938/441595`; hence the lift is Boolean. A new
paired-cube restriction identity gives the explicit all-prime bound

```text
nonzero Boolean quadratic on J(p,(p+1)/2): density >= (p-3)/(4p).
```

Its proof averages over Boolean cubes obtained by leaving one middle-set
point fixed and pairing every other point with the complement. On all
quadratics the transition is `Tf=rho*f+(1-rho)E[f]`, with
`rho=1/(p+1)`; the elementary cube distance `1/4` gives the formula. At
`p=37`, `17/74>5/74`, closing the endpoint.

Filmus's exact slice-junta threshold [51] and the paired Boolean
constructions of Kiermaier--Mannaert--Wassermann [52] are adjacent context,
not inputs to the proof. Targeted searches found no prior statement of the
paired-cube transition bound or this mass-ten use. OEIS API searches for
`171,175,1938,1256`, `328,330,358,360,504`, and the individual large
values `441595` and `1194096750` all returned `null`; no sequence or
priority claim is made.

## 7ah. Integral paired cubes and the `p=29` endpoint (2026-08-28)

Proposition 15.681 observes that 15.680's paired-cube argument applies
before Booleanization. If a nonzero nonnegative integral quadratic has
value `h` at one middle-slice point and scaled mass `c=4p E[B]`, cube
distance gives `c>=p+1-4h`. Combining this with 15.642's stabilizer
inequality eliminates `h` and gives

```text
c >= (p+1)/2 for p=3 mod 4,
c >= (p-1)/2 for p=1 mod 4.
```

The floors `14,16,18,20` remove every positive-residue row at the
`p=29,31,37,41` second all-finite endpoints. At `p=29,s=24`, the sole
residue-zero row has five exact phase-labelled profiles. Globally they are
either 24-arcs with at least four undetermined directions or one-triple
near-arcs with six. In the latter case one deletion gives a 23-arc.

Coolsaet--Sticker [53, Table 5] classify all 25- and 26-arcs in
`PG(2,29)` into 10 and 5 projective classes. Exact `PGL(2,29)` enumeration
gives the same 10 and 5 orbit counts for five- and four-point complements
on a conic, proving that all classified arcs of those sizes are
conic-contained. This agrees with Chao--Kaneta's earlier maximum
nonclassical size 24 [54]. Extending through three undetermined infinity
directions would force one conic through three collinear points, closing
`p=29`. The geometry/classification facts are prior work; the new step is
their use with the Paley endpoint ledger and the integral paired-cube
bound.

OEIS searches identified `27405=C(30,4)` and `142506=C(30,5)` in standard
binomial tables. Searches for the endpoint deficit block, arc-class tail,
and PGL orbit-size block found no relevant sequence. No sequence or
priority claim is made.

## 7ai. The `p=31,s=26` next-boundary endpoint closes (2026-08-28)

Proposition 15.682 applies 15.681's scaled lift floor 16 to remove every
positive phase-zero residue at `p=31`. Exact enumeration of the residue-zero
row leaves fourteen phase-labelled profiles. Eleven are 26-arcs with at
least three undetermined directions. The other three have exactly one
3-secant and at least five undetermined directions; deleting a point of the
triple gives a 25-arc.

Coolsaet's exhaustive classification [55] has no complete arc in
`PG(2,31)` of size 23 through 31. Since every finite arc extends greedily to
a complete arc, every 27- and 28-arc extends to a complete 32-arc and hence,
by Segre's theorem [44], lies on a nondegenerate conic. Adjoin two
undetermined infinity points to a 26-arc, or delete the triple point and
adjoin two to the resulting 25-arc. Repeating with two pairs chosen from
three undetermined infinity points gives conics sharing at least 25 affine
points, so they coincide and would contain three collinear infinity points.
This contradiction closes `p=31,s=26`.

OEIS A000509 records the equivalent complete-arc parameter `m'(31)=22` and
explicitly notes the absence of complete sizes 23 through 31; it was a
useful pointer to [55], not a new-sequence claim. Searches for the endpoint
deficit block `284,286,288,290,360,362,364,366` and the seven global secant
profiles found no additional relevant sequence. Proposition 15.683
subsequently closes `p=41`; the same boundary remains open at
`p=17,19,23`.

## 7aj. The `p=41,s=34` next-boundary endpoint closes (2026-08-28)

Proposition 15.683 applies 15.681's integral lift floor to remove every
positive phase-zero residue. Exact completion-bounded enumeration of the
residue-zero row leaves seven 34-arcs and two one-triple near-arcs. In every
profile, fourteen directions are perfect matchings, twenty miss two points,
and the remaining eight directions contain only three floor-secant units.

For an arc profile, each exceptional direction has at least 28 tangents.
Ball--Lavrauw's polynomial tangent envelope [44, Theorem 11] has degree 18
and restricts to the square of the tangent polynomial on every point
pencil. Root counting therefore makes all eight exceptional direction
pencils double components, leaving a conic. The three exceptional secant
edges touch at least three points; at each such point two remaining tangents
force its point-pencil line into that conic, a degree contradiction.

For a near-arc profile, delete a triple point while preserving the unique
exceptional floor-secant as a pair. The resulting 33-arc has exceptional
tangent counts `33^7,31`. Their eight double components leave a quartic
from the degree-20 envelope. The surviving pair endpoints force two
point-pencil factors; each of the other 31 points then forces its pencil
into the residual conic. This again contradicts degree and closes
`p=41,s=34`. Only `p=17,19,23` remain at this boundary.

Targeted literature searches found no previous use of eight high-tangent
direction pencils with these exact profiles. OEIS searches for
`14,7,20,1`, `18,16,2`, and `34,41,561` returned unrelated results; no
sequence or priority claim is made.

## 7ak. The `p=23,s=20` endpoint reduces to 203 profiles (2026-08-28)

Proposition 15.684 first closes every positive phase-zero residue. The
universal paired-cube floor from 15.681 removes residues 2 through 5. Sharp
cube value floors remove scaled mass 12 and all but the height-four case at
mass 16. In that final case, stabilizer equality makes the quadratic vanish
on the middle shell of `J(23,12)`. The shell-restriction kernel has exact
dimension 23 and is `(t-6)V_1`; one- and two-replacement evaluations of the
affine factor give an immediate modulo-three integrality contradiction.

The exact residue-zero ledger has 1,247 phase-labelled profiles and 485
global shapes. Ball--Lavrauw's tangent envelope [44] excludes all 363 arc
profiles. Coolsaet--Sticker's exhaustive complete-arc classification of
`PG(2,23)` [56] says that complete arcs occur only at sizes
`10,12,13,14,15,16,17,24`, with the unique 24-arc a conic. A line of
occupancy `n` can be repaired by at most one quarter of its pair-slack
contribution. Once this produces an arc of size at least 18, the
classification gives a conic core; an off-conic point count then forces
slack at least 24. This excludes 1,044 profiles in total and leaves exactly
203, so the endpoint is reduced but not closed.

Exact OEIS searches for `363,264,189,136,94,68`, `1247,1044,203`, and
`112449,4341514,1828196` returned no sequence. Individual terms occur in
unrelated entries; no sequence or priority claim is made.

## 7al. The unique `p=23` slack-12 profile closes (2026-08-28)

Proposition 15.685 sharpens the conic-core repair from 15.684. A
slack-twelve realization repairs after at most three deletions. Fewer than
three, or an incomplete 17-arc after three, already gives an 18-point conic
core and contradicts the slack-24 off-conic floor. The remaining branch is
a complete 17-arc `A` plus three outside points.

If `mu_A(x)` is the number of arc secants through an outside point, then
the exact line-slack values at occupancies three, four, and five give
`slack(S)>=4 sum mu_A(x)`. Completeness and total slack twelve force all
three added points to have `mu_A=1`.

Coolsaet--Sticker [56] classify exactly five complete-17-arc classes in
`PG(2,23)`. Five explicit representatives were verified point by point.
Their full outside-point secant-multiplicity histograms are pairwise
distinct, so they are inequivalent and exhaust the five classes. Their
numbers of multiplicity-one points are `0,0,1,0,0`, which excludes the
profile and reduces the exact `p=23` remainder from 203 to 202.

Secant index is standard arc terminology, but targeted literature searches
found no prior table of these five exact histograms. OEIS searches for the
high-count blocks `68,172,190,86` and `69,171,196,78` found no relevant
sequence. No novelty or sequence-submission claim is made.

## 7am. The unique `p=23` slack-16 profile closes (2026-08-28)

Proposition 15.686 reuses the complete-17-arc certificate with the unique
undetermined direction in the slack-sixteen row. At most three repair
deletions, followed by adjoining that infinity point, produce an 18-arc and
the conic-core contradiction. In the four-deletion branch, the repaired
16-arc plus the infinity point is a 17-arc and must be complete.

The undetermined direction ensures that no secant through a deleted point
uses the added infinity point. Completeness and the exact line-slack
inequality force all four deleted points to have secant multiplicity one.
The five classified classes have at most one such point. This excludes the
row and leaves exactly 201 `p=23` profiles, all of slack at least 20.

Targeted literature and OEIS searches found no previous statement matching
this four-point repair obstruction or the exact remaining histogram. No
novelty or sequence-submission claim is made.

## 7an. All 68 `p=23` slack-20 profiles close (2026-08-28)

Proposition 15.687 extends the conic-core count through five deleted points:
for `1<=h<=5` off-conic points, the slack floor `4h(7-h)` has minimum
24. The 68 exact rows have two, three, or four undetermined directions with
counts `2,36,30`. For the 66 rows with at least three directions, two
overlapping infinity-point pairs either extend to the same conic, forcing
three collinear conic points, or produce a complete 17-arc that would need
five multiplicity-one outside points.

For either two-direction row, the only hard branch uses all five repair
deletions. Adjoining both infinity points gives a complete 17-arc, and the
line-slack equality forces all five deleted points to have secant
multiplicity one. The exhaustive five-class certificate permits at most
one. Thus the full block closes, leaving 133 `p=23` profiles of slack at
least 24.

Targeted literature and OEIS searches found no matching statement or
sequence for the exact obstruction or remainder. No novelty or
sequence-submission claim is made.

## 7ao. Cold theorem audit, global route kill, and sharp p=19 lift (2026-08-29)

The original limit was reread independently of the accumulated exact gate.
The weakest sufficient statement is only the ratio-dense Paley-tail deficit

\[
\Phi(C_p)-m_{p^2+1}=o(p^3),
\]

now recorded as Proposition 15.20e. The every-prime exact gap-two four-unit
architecture is strictly stronger.

A general probabilistic counter-mechanism kills local Paley stability:
random signings can be chosen with both \(\Phi(A)\) and
\(\Phi(A\circ C)\) of order \(n^{3/2}\), and greedy edge descent preserves
distance \(\binom n2/2-O(n^{3/2})\) from the signed Paley orbit. The product
frame's complete second moments are also independent of \(A\). The surviving
stability input is instead the all-subsets witness hierarchy supplied by a
closest **global** minimizer.

The augmented cut-code formulation was reduced exactly to the signed
even-Eulerian high-temperature identity

\[
\mathbb E\cosh(\beta Q_a)
=(\cosh\beta)^{\binom n2}P_a(\tanh\beta).
\]

The initially proposed \(\beta=2/\sqrt n\) lower bound
\(\log P_a\ge-o(n)\) would prove the target constant, but a
fractional-moment construction now disproves it by a linear margin. The
corrected fixed-\(c\) target is
\(\log P_a(\tanh(c/\sqrt n))\ge(c/2-c^2/4)n-o(n)\); the construction rules
out every \(c<2.0843108\ldots\), while \(c=3\) remains a clean viable
target. Classical Delsarte and fixed-\(L^q\) norms still discard the phase
consistency needed at this scale.

The conference-class audit corrected another shortcut. Spectral defect zero
identifies all conference classes, not Paley. Craigen's Lemma 7 proves only
that a regular conference matrix has square \(n-1\); it does not prove the
converse. Orders 2, 10, and 26 are regularizable, while order 50 is the first
unresolved square case found in the classification literature. Separately,
Mathon's order-\(5r^2+1\) family has irrational conference eigenvalues and
hence \(\rho<1\) at each order, but no uniform gap is known.
Momihara--Suda's excess bound can approach one at the
\(\Theta(n^{-3/2})\) normalized scale, while Mathon's visible quotient only
gives \(\liminf\rho\ge8/(5\sqrt5)\). The exact missing theorem is
switched-row variance \(\Omega(r^2)\) for every Boolean switching.

On the live residual front, Proposition 15.688 combines paired-cube
quarter-integrality with Proposition 15.642's exact stabilizer weights to
prove the sharp theorem

\[
4p\mathbb E B\ge p-3
\]

for every nonzero nonnegative integral quadratic on the middle slice. The
Boolean quadratic \((1-x_i)(1-x_j)\) attains equality. At the \(p=19,s=16\)
second boundary, this removes every positive residue. The minimum
residue-zero pair has impossible slack 34 but is not the full row: exact
completion leaves 143 profiles. Proposition 15.689 uses the complete-arc
spectrum and conic-core counting to exclude all 129 profiles of slack at
most twelve, leaving fourteen at slack 16--32. The endpoint remains open.

The full-Max+ dilation energy was also normalized exactly. If \(S_K\) is
the square-torus energy in the cold strategy, then
\[
S_K=12(q-1)\|\delta\|^2/n.
\]
Thus \(S_K\le q-1\) is precisely strong R1, not an auxiliary lemma.
Representation/character data and additive-autocorrelation PSD do not imply
it: explicit abstract spectra and PSD autocorrelations violate the target.
Any proof must use the Boolean rank-one identity and exact cancellation
among all Max+ orbit types.

The degree-four truncated-moment import of
Infusino--Kuna--Lebowitz--Speer and Prékopa was also evaluated. It proves
that the constant quartic majorant is already optimal in the large-prime
bulk, so quartics cannot improve the existing \(2p\) floor there. At
\((p,b,\eta)=(19,9,0)\), an explicit negative quartic expectation prevents
blind extension to the small endpoint.

Full audit and attack map:
`evidence/STRATEGY_2026-08-29_COLD_REVIEW.md` and
`evidence/NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md`.

## 7ap. Binary affine-Radon normal form for p=19 (2026-08-29)

Proposition 15.692 rewrites the fourteen profiles left by 15.689 without a
solver relaxation. If `A` is the binary affine line-point incidence matrix,
then `A^T A=I+J` for odd order. It follows by equal dimensions that `A` is
an isomorphism from even point words to tuples of even directional line
words, with inverse `x=A^T r`. Thus the p=19 parity-profile problem has no
unseen linear compatibility condition: its exact remaining equation is
`wt(A^T r)=16`.

The inverse-weight congruence modulo four accepts every survivor. The fixed
first two moments of the selected-stripe count also admit exact nonnegative
distributions supported on the even values `4,6,8`, so pairwise moments
cannot force positive odd density. This is a rigorous method barrier and
normal form, not endpoint closure. Searches on binary affine incidence
codes, odd secants, and Rédei direction methods found adjacent literature
but no theorem supplying the missing fixed-block inverse-weight bound.

## 7aq. The p=19 slack-sixteen block closes (2026-08-29)

Proposition 15.693 uses the classified secant-index distribution of
14-arcs in `PG(2,19)`. Each of the seven slack-sixteen profiles has at least
three undetermined directions and repair depth at most four. Depth at most
three gives a conic contradiction after adjoining two infinity points. At
depth four, the repaired 12-arc plus two infinity points is a complete
14-arc `K`; otherwise the complete-arc spectrum again extends it to a conic.

The four deleted points have secant index one with respect to `K`, and each
unused undetermined infinity point contributes another index-one outside
point. Thus `c1(K)>=5`. Al-Zangana's exhaustive list of all 83 projective
14-arc classes gives `c1<=4`, a contradiction. The exact p=19 remainder
drops from fourteen to seven profiles with slack histogram
`{20:4,24:1,28:1,32:1}`. No top-level implication changes.

## 7ar. The p=19 slack-twenty equality case becomes rigid (2026-08-29)

Proposition 15.694 combines the exact five-deletion depth from 15.693 with
equality in the line-secant slack bound. A putative witness is the disjoint
union of an 11-arc and a 5-arc; every deleted point lies on exactly one core
secant. Only eight per-line `(core,deleted)` occupancy pairs survive, and
the bad lines have one of three patterns: five triples, one quadruple plus
three triples, or two quadruples plus one triple.

Adding any two undetermined infinity points gives a 13-arc with at least
seven or eight outside points of secant index one. Al-Zangana's exhaustive
13-arc classification permits at most nine, so this filters but does not
exclude the four profiles. Bounded exact SAT and CP-SAT trials returned
UNKNOWN and were not promoted as evidence.

## 7as. Three p=19 slack-twenty rows sharpen (2026-08-29)

Proposition 15.695 excludes both rows whose exceptional phase-one direction
has `b=14`: equality forces the directional quadratic to equal one on three
intersection layers, and a fixed pair-inclusion minor has full rank 171
modulo 101. Proposition 15.696 handles the mixed `b=16` row. Its equality
layers have rank 169 and two integral kernel orbits; all twenty corrected
infinity-degree edge-lift shards are infeasible. A subsequent audit found and
repaired integer subtraction on encoded `F_{19^2}` elements. Every logical
shard was rerun with componentwise field subtraction and a full canonical
conference-sign regression. The hard `022/I=28` shard was partitioned by all
three possible elevated phase-zero roles, so the twenty logical shards have
22 exact raw certificates. Only this corrected archive supports 15.696.
The p=19 endpoint is thereby reduced to four profiles with slack histogram
`{20:1,24:1,28:1,32:1}`.

Proposition 15.697 gives a structural reduction of the all-`b=2` slack-20
profile. Stabilizer equality and an exact rank-152 layer kernel reduce a
hypothetical maximum-five lift to an additive ten-by-nine cross-difference
matrix. Exhaustion of all `2^18` first-row/first-column assignments leaves
21 labelled matrices in four symmetry types, each with an explicit negative
layer; hence the elevated lift is Boolean. Exact coefficient `l1` bounds for
the five phase-zero zero-lifts reduce the infinity degree to `0,20,38`.

Filmus's publication page states that a short note with Antoine Vinciguerra
proves the Boolean restriction threshold equals the sharp junta threshold
`2d`. Conditional on that external statement, a four-variable cube audit
classifies 3,420 possible Boolean lifts. The note's SharePoint PDF returned
HTTP 403 during this audit, so the catalog is not used as unconditional
evidence. Bounded edge-lift trials returned `UNKNOWN`; no profile or
top-level gate closes in 15.697.

## 7at. The p=19 slack-twenty block closes (2026-08-29)

Proposition 15.698 imports 15.694's forced 11-arc plus five-deletion repair
into the exact binary affine-Radon model for the final all-`b=2` profile.
The model has 1,184,892 clauses, 741 native XOR equations, and 776 exact
cardinality constraints. CryptoMiniSat 5.11.21 returned `UNSATISFIABLE` on
nuka and independently on soulkiller's registered-ECC CPU. The lossless
normalization chooses a retained core point and its partner on a phase-zero
`b=0` line and sends them to zero and one by a square affine similarity.

Thus the boundary itself is impossible; this is not an edge-lift timeout.
Nonsquare dilation transfers the two product signs. All p=19 slack-twenty
profiles are closed, reducing the endpoint from four profiles to three with
histogram `{24:1,28:1,32:1}`. The endpoint and every top-level gate remain
open.

## 7au. The p=19 second all-finite endpoint closes (2026-08-29)

Proposition 15.699 applies the exact affine-Radon inverse model directly to
the three profiles left by 15.698. No edge-lift variables, floor relaxation,
or repair hypothesis is used. Slack 24 is UNSAT on nuka and soulkiller ECC;
slack 28 is UNSAT on soulkiller ECC; slack 32 is UNSAT on jellyfin and
soulkiller ECC. Each model has 361 point variables, 380 affine line parity
variables, 741 native XOR equations, and the exact directional cardinality
histogram. Nonsquare dilation transfers the two signs.

Consequently the p=19 second all-finite endpoint is closed. The remaining
next-boundary endpoints are p=17 and p=23; later boundary sizes and all
top-level gates remain open.

## 7av. The p=17 slack-zero block reduces to two conic profiles (2026-08-29)

Proposition 15.700 derives the exceptional `p=17,s=16` arithmetic ledger
directly. The sharp `p-3=14` integral-lift floor excludes phase-zero residues
2 through 6; exact completion, the pair budget, and slack divisibility leave
1,575 phase-labelled profiles. Exactly 247 have pair slack zero.

Slack zero makes the boundary a 16-arc. Sticker's exhaustive classification
[48] has one 16-arc class in `PG(2,17)`, represented by conic-minus-two.
Fixing one conic and enumerating every line at infinity and eligible deleted
pair gives 21,267 affine cases: 20,808 external, 306 tangent, and 153 secant.
Their exact directional Paley census has 53 profiles after including the
nonsquare phase swap. Only two occur in the arithmetic list, both tangent at
infinity. Therefore 245 rows are excluded and the exact p=17 remainder falls
to 1,330 profiles, including two zero-slack rows.

The follow-up fixed-boundary model adds all coefficient identities forced by
floor equality. Five-minute Soulkiller and Nuka runs returned `UNKNOWN`, so
they have no evidentiary role. The p=17 endpoint and all top-level gates
remain open. Targeted literature and OEIS searches found no prior listing of
the 247-to-2 directional profile reduction.

## 7aw. The p=17 low-positive-slack conic-core reduction (2026-08-29)

Proposition 15.701 uses the same exhaustive p17 arc classification [48], now
at size fifteen. There is one PGL class of 15-arcs, and conic-minus-three is
a representative, so every 15-arc is conic-contained.

Pair slack `4r` repairs to an arc after at most `r` deletions. The exact
profile census shows that the repaired arc reaches size fifteen directly for
all 227 slack-four profiles, after adjoining one undetermined infinity point
for 128 of 195 slack-eight profiles, and after adjoining two for 43 of 155
slack-twelve profiles. If `h<=3` original points lie off the resulting conic,
retained conic secants force slack at least `4h(6-h)>=20`; if `h=0`, the
original boundary is itself an arc. Thus all 398 qualifying profiles are
impossible.

The exact p17 remainder falls from 1,330 to 932 rows: two at slack zero, 67
at slack eight, 112 at slack twelve, and 751 at slack at least sixteen. The
endpoint and every top-level gate remain open. Targeted literature searches
found no prior statement of this profile reduction. OEIS has no matching
count sequence; the isolated values 398, 932, and 1330 occur together in
A137707 but its Secondary Wythoff Array construction has no identified
connection and is not evidence.

## 7ax. The unique complete p17 14-arc excludes 146 more rows (2026-08-29)

Proposition 15.702 uses Sticker's complete-arc table [48], which has one
complete 14-arc class in `PG(2,17)`. An explicit representative was checked
over all 307 projective lines and points. Its line occupancy census is
`{0:146,1:70,2:91}`, and its outside secant-index census is
`{2:4,3:4,4:76,5:128,6:75,7:6}`. In particular, the unique class has no
outside point of secant index one.

Slack-eight equality would require a two-point repair to that complete arc
with both deleted points of index one, so all remaining 67 slack-eight rows
are impossible. For slack twelve, one undetermined direction extends the
three-deletion repaired 13-arc to a 14-arc. If it is incomplete, the conic
core argument applies; if complete, the undetermined direction ensures the
three deleted points each retain secant index at least two, forcing slack at
least 24. This excludes the 79 one-direction rows.

The p17 remainder drops from 932 to 786 exact profiles: two at slack zero,
33 at slack twelve, and 751 at slack at least sixteen. The endpoint remains
open. Targeted literature and OEIS searches found no prior statement of the
secant-index census/profile exclusion.

## 7ay. Complete p17 13-arcs close the slack-twelve block (2026-08-29)

Proposition 15.703 treats the 33 remaining slack-twelve profiles, all with no
undetermined direction. Equality repairs any realization by three deletions
to a 13-arc. A normalized PGL generator fixes a quadrangle, enumerates and
blocks each representative's complete projective orbit, and produces eight
pairwise inequivalent complete 13-arcs. Their stabilizer-order fingerprint
`1,2,2,2,2,3,4,6` matches Sticker's exhaustive table [48], whose published
count of eight makes the generated list complete.

The eight outside secant-index histograms have index-one-point counts
`0,0,0,0,0,0,2,3`. Slack-twelve equality requires three such points. The
only candidate triple reconstructs slack sixteen, so no complete-core case
survives. If the repaired 13-arc is incomplete, extend it to a 14-arc. The
conic branch was excluded by Proposition 15.701. Deleting each point from
the unique complete 14-arc gives index-one-count histogram `{0:4,1:8,4:2}`;
all eight candidate triples reconstruct slack twenty. If the original
deleted set contains the extension point, the other two points each have
complete-14-arc secant index at least two and force slack at least sixteen.

Thus all 33 slack-twelve profiles are impossible, conditionally only on the
published complete-class count. The p17 remainder is 753 exact profiles:
two tangent-conic slack-zero rows and 751 rows of slack at least sixteen.
The endpoint and every top-level gate remain open. Targeted literature and
OEIS searches found no prior statement of the eight outside histograms, the
index-one vector, or this profile exclusion.

## 7az. Free directions exclude 99 p17 slack-sixteen rows (2026-08-29)

Proposition 15.704 splits the 112 slack-sixteen profiles by undetermined
directions as `{0:13,1:47,2:47,3:5}`. Repair takes at most four deletions.
With two directions, their infinity points produce a 14-arc: the complete
class forces slack at least 32 and the incomplete branch reaches a conic,
whose positive off-conic slack floor is at least 20. With three directions,
two overlapping infinity-point pairs must extend to the same conic, which
cannot contain three collinear infinity points. This excludes 52 rows.

For one direction and repair depth four, the repaired arc plus its infinity
point is a 13-arc. A complete core would need four outside index-one points,
but the eight complete classes have counts `0,0,0,0,0,0,2,3`. An incomplete
core extends either to a conic or to the unique complete 14-arc. In the latter
case, the 14 possible one-point deletions supply two index-one quadruples and
26 raw infinity placements with slack histogram `{16:2,28:16,32:8}`. The two
slack-sixteen placements fail the exact undetermined-direction condition;
the eight valid placements all have slack 32. Thus all 47 one-direction rows
also fail.

The exact p17 remainder drops from 753 to 654: two slack-zero profiles,
thirteen zero-direction slack-sixteen profiles, and 639 profiles of slack at
least twenty. The endpoint and every top-level gate remain open. Targeted
literature and OEIS searches found no prior statement of the infinity-
placement census or this 99-profile exclusion; isolated numerical matches
are unrelated and are not evidence.

## 7ba. The p17 slack-sixteen block closes (2026-08-29)

Proposition 15.705 exhausts all 629 PGL classes of twelve-arcs in
`PG(2,17)` and all 97,122 four-point extensions within the exact
slack-sixteen core-secant charge. Only 47 extensions have an allowed line
pattern. Across every affine chart and both phase labellings, all 6,345
resulting cases miss the thirteen target profiles left by 15.704. Therefore
all thirteen zero-direction slack-sixteen rows are impossible. The p17
remainder drops from 654 to 641 profiles: two at slack zero and 639 of slack
at least twenty. The endpoint and all top-level gates remain open.

## 7bb. A global Paley-sign identity closes p17 slack zero (2026-08-29)

Proposition 15.706 excludes the two slack-zero profiles left by 15.700.
Every possible mean allocation retains a rigid `b=2` direction of each
quadratic type. For such a direction, summing all 136 inter-fibre coefficient
identities and comparing with its exact mean gives

```text
P_d = 4 + 8g_d - I.
```

Writing `S` for the signed sum of every selected finite edge gives a second
expression for the same directional cross sum. Comparing one direction of
each type eliminates `S` and forces

```text
17I = 4 + 72(g_+ + g_-).
```

Hence `I=68 (mod 72)`, and the full range `0<=I<=69` leaves only `I=68`.
There is then one finite edge, so the affine odd boundary has size 66, 68,
or 70 rather than 16. Both profiles are impossible without a solver or any
new classification input. The exact p17 remainder is 639 profiles, all of
pair slack at least twenty; the endpoint and every top-level gate remain
open.

Targeted GitHub code and MathOverflow searches for the distinctive
`17I=4 (mod 72)` identity, the 639-profile p17 remainder, and global
directional Paley-sign sums found no prior statement of this exclusion.

The same audit found that three exploratory edge-lift scripts had subtracted
encoded `F_{p^2}` integers instead of subtracting their two field components.
The scripts now use componentwise subtraction and are checked against the
canonical Paley conference matrix by a full-edge regression test. This
invalidated the original raw CP-SAT archive used by Proposition 15.696; its
corrected twenty-shard rerun is recorded separately rather than silently
reusing the old files.

## 7bc. The p17 slack-twenty block closes (2026-08-29)

Proposition 15.707 first extends 15.706's global-sign argument to all 69
slack-twenty profiles with residue pair `(u_0,u_1)=(0,8)`. Exact quotient
minima force between three and seven rigid phase-zero directions with
`b in {0,2}`, while phase one retains at least eight rigid `b=2` directions.
The `b=0` floor `(M,T)=(0,0)` and `b=2` floor `(18,-1)` give the same
phase-zero global-sign constant. Comparing with phase one again forces
`I=68`, whose one finite edge cannot have affine odd boundary 16.

The other nine profiles have residue pair `(8,8)` and two or three
undetermined directions. Adjoining two infinity points after a minimum arc
repair handles every depth. Depth at most three reaches a conic and gives
the improved positive slack floor `4h(7-h)>=24`. Depth four reaches the
unique complete 14-arc, whose minimum outside secant index two forces slack
at least 32. At depth five, equality needs five index-one outside points;
the eight complete 13-arcs permit at most three, and every 13-subarc of the
unique complete 14-arc permits at most four. These are exactly the
classification certificates already audited in 15.701--15.703.

Thus all 78 slack-twenty profiles are impossible. The p17 remainder drops
from 639 to 561 profiles, all of pair slack at least twenty-four. The endpoint
and every top-level gate remain open. No new classification or solver is
used.

Targeted GitHub and MathOverflow searches found no prior statement of this
69-plus-9 exclusion or the improved `4h(7-h)` undetermined-conic bound. An
OEIS search on the full reduction ledger found only isolated numerical
matches in unrelated constructions; none supplies mathematical evidence.

## 7bd. A unique-even-fibre identity closes p17 slack twenty-four (2026-08-29)

Proposition 15.708 excludes all 54 slack-twenty-four profiles. The 45 rows
with residue pair `(0,8)` retain rigid phase-zero `b=0` and phase-one `b=2`
directions. Their global-sign comparison forces `I=68` and gauge sum 16,
whereas nonnegative parallel counts require the two gauges to sum to at least
17.

The nine `(8,8)` rows retain rigid phase-zero `b=16` directions. Comparing
with phase-one `b=2` forces `I=4`, phase-zero gauge one, 64 phase-zero finite
edges, and one phase-one finite edge. For the unique even fibre `j` of the
canonical floor `1-x_j`, summing the exact coefficient cells incident with
`j` gives

```text
N_j-delta_j = -15z_j-I.
```

Here `N_j` is a nonnegative phase-zero crossing-edge count, `delta_j<=1`,
and `z_j>=0`; hence `N_j<=-3`. This closes the block without a solver or a
new arc classification. The p17 remainder falls from 561 to 507 profiles,
all of pair slack at least twenty-eight. The endpoint and all top-level gates
remain open.

Targeted GitHub code and MathOverflow searches for the distinctive
`17I=-4 (mod 72)`, `N_j-delta_j=-15z_j-I`, and unique-even-fibre Paley
identities found no prior statement of this exclusion. OEIS searches on the
reduction counts returned only unrelated isolated numerical matches and are
not evidence.

## 7be. The rigid-anchor identities remove every p17 `u_1=8` row (2026-08-29)

Proposition 15.709 applies 15.708 uniformly across the full live ledger. All
280 profiles with phase-one residue eight retain a rigid phase-one `b=2`
core. Of these, 66 have phase-zero residue zero and retain rigid `b=0`; the
other 214 have phase-zero residue eight and retain rigid `b=16`. The global
gauge and unique-even-fibre contradictions respectively exclude the two
blocks without reference to their pair slack.

The p17 remainder falls from 507 to 227 profiles, all with phase-one residue
zero and pair slack at least 96. Their residue split is `(0,0):181`,
`(7,0):9`, `(8,0):37`. The endpoint and all top-level gates remain open.

Targeted GitHub code and MathOverflow searches found no prior statement of
the 280-row rigid-anchor sweep. OEIS searches on `507,280,227` and the
surviving residue counts returned only unrelated isolated matches; they are
not mathematical evidence.

## 7bf. Complementary phase-one `b=16` identities leave nineteen p17 rows (2026-08-29)

Proposition 15.710 uses the nine rigid phase-one `b=16` directions present in
every row left by 15.709. For 176 rows, an actual floor-rigid phase-zero
`b=0` anchor forces
infinity degree 60 and gauge sum 14, while nonnegative parallel counts force
gauge sum at least 15. Thirty-two rows have `b=16` anchors in both phases;
their comparison forces infinity degree 68 and gauge sum 16, while
nonnegativity forces at least 17.

The p17 remainder falls from 227 to nineteen profiles, with residue split
`(0,0):5,(7,0):9,(8,0):5`. The endpoint and all top-level gates remain open.

Targeted GitHub code and MathOverflow searches found no prior statement of
the complementary `17I=12 (mod 72)` gauge contradiction. OEIS returned only
unrelated arithmetic sequences; those numerical matches have no bearing on
the proof.

## 7bg. Uniform directional means exclude the five p17 residue-zero rows (2026-08-29)

Proposition 15.711 handles the five `(u_0,u_1)=(0,0)` profiles left by
15.710. Avoiding its rigid phase-zero `b=0` anchor consumes every free
quotient increment, so every direction in both phases has mean 18. The two
global mean identities reduce the infinity degree to `6,24,42,60`. In each
case the rigid phase-one `b=16` identity forces all finite edges into phase
one and fixes the gauge `g=1,3,5,7`. Nonnegative cross-cell counts then give
`I<=g+1+15 floor(g/2)`, with respective upper bounds `2,19,36,53`, excluding
all four candidates.

The p17 remainder falls from nineteen to fourteen profiles, with residue
split `(7,0):9,(8,0):5`. The endpoint and all top-level gates remain open.

Targeted GitHub-code and MathOverflow searches found no prior occurrence of
the exact uniform-mean or fibre-capacity identities. Related Paley graph,
character-sum, spectral, and equitable-partition literature contains no
matching exclusion. OEIS matches for `6,24,42,60` and `2,19,36,53` are
unrelated numerical coincidences and have no evidentiary role.

## 7bh. Szőnyi's direction bound closes the p17 endpoint (2026-08-29)

Proposition 15.712 observes that all fourteen rows left by 15.711 have
phase-one profile `{16:9}`. For a sixteen-point affine set, `b_d=16` forces
sixteen singleton fibres, so none of those nine directions is determined.
The boundary therefore determines at most nine directions. Szőnyi's theorem
requires a noncollinear `k<=p` subset of `AG(2,p)` to determine at least
`(k+3)/2` directions, which is ten at `(k,p)=(16,17)`. The boundary must be
collinear, but its resulting profile `{0:1,16:8}/{16:9}` is absent.

All fourteen profiles are excluded and the `p=17,s=16` endpoint is closed.
Residual (ii) and every top-level gate remain open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found the
general theorem but no prior occurrence of the exact nine-undetermined-
direction, phase-labelled `AG(2,17)` application. Numerical OEIS hits are
unrelated and have no evidentiary role.

## 7bi. Szőnyi cuts the positive p7 infinity-plus-seven envelope (2026-08-29)

Proposition 15.713 projects the positive-product branch to the two labelled
multisets of four odd-fibre counts. The phase-zero floors and pair-deficit
budget leave 1,217 ordered projected `b`-profile pairs. If at least four of
the eight directions have `b_d=7`, the seven-point affine boundary determines
at most four directions, so Szőnyi's theorem forces it to be a complete
affine line. Exactly two labelled line profiles survive. The other 208 pairs
are excluded, leaving a 1,009-profile projected outer envelope. This does
not count residue/quotient-labelled states and does not close the branch.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found the
general theorem but no prior occurrence of this exact `p=7` application,
histogram, or reduction. Numerical OEIS matches are unrelated and have no
evidentiary role.

## 7bj. Complete mod-seven exclusion of positive p7 infinity-seven z0 (2026-08-29)

Proposition 15.714 uses the unique phase-zero mean-eight Johnson catalogs to
turn each actual seven-point affine boundary with no undetermined direction
into a fixed right side of the 282-by-1,225 edge system. Complete V100
combinadic scans with two launch geometries test all 85,900,584 boundaries
against the 135 audited mod-seven left dependencies. Exactly 79,447,032 have
`z=0`, and none survives. The positive branch falls to 6,453,552 actual
boundaries and a 792-profile projected outer envelope; it remains open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found no
prior occurrence of the exact boundary census or mod-seven exclusion. The
generic `C(49,7)` total is catalogued, but the histogram and zero-survivor
result are not; numerical matches have no evidentiary role.

## 7bk. Complete mod-seven exclusion of positive p7 infinity-seven z1 (2026-08-29)

Proposition 15.715 observes that each actual positive `z=1` boundary has
exactly four mean allocations.  For each variable direction, 23 exact
mod-seven dependencies give a losslessly packed prefilter.  Complete V100
scans with two launch geometries agree on all 6,324,528 `z=1` boundaries and
the same 1,326 projected boundary candidates.  Host reconstruction tests all
four allocations against the complete catalogs on all 135 dependencies;
none survives.  The positive branch falls to 129,024 actual boundaries and
a 492-profile projected outer envelope at `z=2,3,7`; it remains open.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found only
the surrounding finite-direction, Paley-49, and modular-incidence theory,
not this exact census or exclusion.  Numerical matches have no evidentiary
role.

## 7bl. Pair-transversal and modular closure of positive p7 z2 (2026-08-29)

Proposition 15.716 enumerates the 141,120 pair-transversal incidences of the
positive infinity-plus-seven remainder and decomposes its 129,024 actual
boundaries into 104 orbits.  The `z=2` branch has 123,480 boundaries in 92
orbits.  Its complete mean ledger has 1,232 leaves: an independent audit
caught 48 residue-four leaves omitted by the first residue-zero count.

Using the unpointed 281-row translation-equivariant edge system, 112 exact
dependencies annihilate each one-high catalog, while full 135-coordinate
hash joins test every two- and four-catalog leaf.  Nuka rejects all 1,232
leaves modulo seven.  The positive remainder falls to 5,544 actual
boundaries in twelve orbits and a 212-profile projected envelope at `z=3,7`.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found no
prior occurrence of this orbit/catalog exclusion.  OEIS contains 123,480 in
an unrelated injective-function triangle; that numerical match has no
evidentiary role.

## 7bm. Same-tuple multimodular closure of positive p7 z3 (2026-08-30)

Proposition 15.717 reconstructs the ten exact affine-semilinear `z=3`
boundary orbits, covering 5,488 actual boundaries, and all 400 corrected mean
leaves.  Exact joins on the full 135-dimensional mod-seven dependency space
reject 398 leaves.  The other two leaves have three complete 1,764-row
catalogs and exactly four matching catalog-row triples apiece.

The follow-up retains those exact row indices, rebuilds each 281-entry integer
edge-system right side, and tests the same right side against complete
left-nullspace bases modulo `3,5,7,11`.  All eight mod-seven tuples fail
modulo three.  Therefore the full `z=3` branch is excluded.  The positive
remainder is 56 line boundaries in two `z=7` orbits and two projected line
profiles.

Targeted GitHub-code, MathOverflow, literature, and OEIS searches found no
prior occurrence of this exact orbit/catalog and same-tuple multimodular
exclusion.  Searches on 5,488 and 225,792 produced unrelated neural-network,
divisor, and OEIS matches; those numerical coincidences have no evidentiary
role.

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
39. K. Tasaka, *Harmonic strength of shells of lattices and weighted theta series*, [arXiv:2308.14309](https://arxiv.org/abs/2308.14309) (survey of the shell-harmonic/modular-form bridge; no Paley channel inequality).
40. M. Ozeki, *On the Siegel theta series of extremal lattices and their association schemes*, Kyushu J. Math. **68** (2014), 53--73, [doi:10.2206/kyushujm.68.053](https://doi.org/10.2206/kyushujm.68.053) (adjacent orbital/theta machinery; not a proof of R1).
41. P.-A. Bernard, N. Crampé, and L. Vinet, *A bivariate Q-polynomial structure for the non-binary Johnson scheme*, J. Combin. Theory Ser. A **202** (2024), 105829, [arXiv:2306.01882](https://arxiv.org/abs/2306.01882) (adjacent Johnson-scheme machinery; not the parity-majorant or affine pair-deficit theorem of 15.669).
42. S. Ball and B. Csajbók, *On sets of points with few odd secants*, [arXiv:1711.10876](https://arxiv.org/abs/1711.10876) (odd-secant bounds for \(q+2\)-point projective sets; not the eight-point split-type census of 15.670).
43. M. Kiermaier and S. Kurz, *Maximal integral point sets in affine planes over finite fields*, [arXiv:1401.2825](https://arxiv.org/abs/1401.2825) (prescribed directions and Paley-clique context; not the boundary parity-floor budget of 15.670).
44. S. Ball and M. Lavrauw, *Planar arcs*, J. Combin. Theory Ser. A **160** (2018), 261--287, [doi:10.1016/j.jcta.2018.06.015](https://doi.org/10.1016/j.jcta.2018.06.015), [arXiv:1705.10940](https://arxiv.org/abs/1705.10940) (modern primary restatement of Segre's odd-order arc theorems; the `q`-arc conic theorem is used in 15.673, and Theorem 11's polynomial tangent envelope is used in 15.683).
45. V. F. Lev, *Point distribution and perfect directions in* \(F_p^2\), [arXiv:1903.01518](https://arxiv.org/abs/1903.01518) (Rédei-type direction bounds adjacent to, but not implying, the determined-direction/type-budget argument of 15.674).
46. G. Van de Voorde, *On sets without tangents and exterior sets of a conic*, [arXiv:1201.0484](https://arxiv.org/abs/1201.0484) (adjacent external-conic and determined-direction geometry; not the Paley phase-budget exclusion of 15.676).
47. P. Amireddy, S. Behera, S. Srinivasan, and M. Sudan, *A Near-Optimal Polynomial Distance Lemma over Boolean Slices*, ICALP 2025, [doi:10.4230/LIPIcs.ICALP.2025.11](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2025.11) (Lemma 2 supplies the exact finite degree-two support floor used in 15.642 and 15.677).
48. H. Sticker, *Classification of Arcs in Small Desarguesian Projective Planes*, PhD thesis, Ghent University, 2012, [PDF](https://cage.ugent.be/geometry/Theses/57/PhDHeideSticker.pdf) (Section 5.3 gives the exhaustive PGL class counts for arcs not necessarily complete; the unique 16-arc class in `PG(2,17)` is used in 15.678).
49. G. Keri, *Types of superregular matrices and the number of n-arcs and complete n-arcs in PG(r,q)*, J. Combin. Des. **14** (2006), 363--390, [doi:10.1002/jcd.20091](https://doi.org/10.1002/jcd.20091) (earlier MDS/superregular classification in the same small-field range; contextual independent cross-check for [48]).
50. P. Amireddy, A. R. Behera, S. Srinivasan, M. Sudan, and S. V. Willumsgaard, *Low-Degree Testing Over Boolean Slices*, [arXiv:2608.21730](https://arxiv.org/abs/2608.21730) (robust low-degree testing on slices; not an exact distance-floor or Paley affine-boundary theorem).
51. Y. Filmus, *Junta threshold for low degree Boolean functions on the slice*, [arXiv:2203.04760](https://arxiv.org/abs/2203.04760) (proves the sharp `k>=2d` threshold for bounded-size juntas; adjacent context, not the explicit paired-cube density bound in 15.680).
52. M. Kiermaier, J. Mannaert, and A. Wassermann, *The paired construction for Boolean functions on the slice*, [arXiv:2510.02804](https://arxiv.org/abs/2510.02804) (small-support constructions and exact degrees; does not classify all Boolean quadratics or imply the mass-ten exclusion in 15.680).
53. K. Coolsaet and H. Sticker, *The complete k-arcs of PG(2,27) and PG(2,29)*, J. Combin. Des. **19** (2011), 111--130, [doi:10.1002/jcd.20261](https://doi.org/10.1002/jcd.20261), [open PDF](https://backoffice.biblio.ugent.be/download/1247338/1247417) (Table 5 gives the exhaustive projective class counts 10 and 5 for all 25- and 26-arcs in `PG(2,29)`, used in 15.681).
54. J. M. Chao and H. Kaneta, *Classical arcs in PG(r,q) for 23 <= q <= 29*, Discrete Math. **226** (2001), 377--385 (records maximum nonclassical plane-arc size 24 at `q=29`; independent prior confirmation of 15.681's classification consequence).
55. K. Coolsaet, *The Complete Arcs of PG(2,31)*, J. Combin. Des. **23** (2015), 522--533, [doi:10.1002/jcd.21410](https://doi.org/10.1002/jcd.21410) (exhaustive complete-arc classification used in 15.682; in particular, there are no complete arcs of sizes 23 through 31).
56. K. Coolsaet and H. Sticker, *A full classification of the complete k-arcs of PG(2,23) and PG(2,25)*, J. Combin. Des. **17** (2009), 459--477, [doi:10.1002/jcd.20211](https://doi.org/10.1002/jcd.20211) (complete-arc spectrum and class counts in `PG(2,23)`, used in 15.684 to force every arc of size at least 18 onto a conic).
57. R. Mathon, *Symmetric conference matrices of order \(pq^2+1\)*, Canad. J. Math. **30** (1978), 321--331, [doi:10.4153/CJM-1978-029-1](https://doi.org/10.4153/CJM-1978-029-1) (ratio-dense non-square-eigenvalue conference family; a uniform Boolean-radius gap is not known).
58. R. Craigen, *Trace, Symmetry and Orthogonality*, Canad. Math. Bull. **37** (1994), 461--467, [doi:10.4153/CMB-1994-067-1](https://doi.org/10.4153/CMB-1994-067-1) (Lemma 7 proves regular conference \(\Rightarrow n=k^2+1\), not the converse).
59. A. J. L. Paulus, *Conference matrices and graphs of order 26*, Eindhoven report, 1973, [PDF](https://pure.tue.nl/ws/files/4424819/252855.pdf) (the four order-26 conference classes have regular representatives).
60. M. Infusino, T. Kuna, J. L. Lebowitz, and E. R. Speer, *The truncated moment problem on \(\mathbb N_0\)*, J. Math. Anal. Appl. **452** (2017), 443--468, [arXiv:1504.02989](https://arxiv.org/abs/1504.02989) (degree-four contact-measure criterion used to audit quartic parity majorants).
61. A. Prékopa, *The discrete moment problem and linear programming*, Discrete Appl. Math. **27** (1990), 235--254, [doi:10.1016/0166-218X(90)90068-N](https://doi.org/10.1016/0166-218X(90)90068-N) (finite-support moment LP and endpoint constraints).
62. D. Bryant, D. Horsley, and W. Pettersson, *Cycle decompositions V: complete graphs into cycles of arbitrary lengths*, Proc. Lond. Math. Soc. **108** (2014), 1153--1192, [arXiv:1204.3709](https://arxiv.org/abs/1204.3709) (prescribed complete-graph cycle decompositions; contextual support for the broad Eulerian dual spectrum).
63. G. Faina, S. Marcugini, A. Milani, and F. Pambianco, *The spectrum of values k for complete k-arcs in PG(2,q) for q<=23*, Ars Combinatoria **47** (1997), 3--11, [open article](https://combinatorialpress.com/ars-articles/volume-047-ars-articles/the-spectrum-of-values-k-for-complete-k-arcs-in-pg2q-for-q-leq-23/) (complete-arc spectrum of `PG(2,19)` used in 15.689).
64. K. Momihara and S. Suda, *Conference matrices with maximum excess and two-intersection sets*, [arXiv:1611.01305](https://arxiv.org/abs/1611.01305) (Proposition 1.1 gives the arithmetic maximum-excess bound used to delimit the Mathon fixed-gap route).
65. E. B. Al-Zangana, *The Geometry of the Plane of Order Nineteen and its Application to Error-Correcting Codes*, PhD thesis, University of Sussex, 2011, Chapter 4, Section 4.22, p. 105 (exhaustive `PG(2,19)` 14-arc class counts and the `c1<=4` bound used in 15.693).
66. Y. Filmus and A. Vinciguerra, short note on the restriction threshold for bounded-degree functions on slices, linked from [Filmus's publication page](https://yuvalfilmus.cs.technion.ac.il/publications/papers/) (the page states the arithmetic-progression, hence Boolean, restriction-threshold result; the linked PDF was inaccessible during the 15.697 audit, so the result is used only conditionally).
67. T. Szőnyi, *On the number of directions determined by a set of points in an affine Galois plane*, J. Combin. Theory Ser. A **74** (1996), 141--146, [doi:10.1006/jcta.1996.0042](https://doi.org/10.1006/jcta.1996.0042) (the `k<=p` direction bound used in 15.712).
68. G. Somlai, *A new proof of Rédei's theorem on the number of directions*, Arch. Math. **122** (2024), 575--580, [doi:10.1007/s00013-024-01979-x](https://doi.org/10.1007/s00013-024-01979-x) (modern explicit restatement of Szőnyi's bound and proof context).

Secondary: Paley construction (Wikipedia, quoting the 1933 Hadamard-conjecture sentence); conference matrix (Wikipedia, van Lint–Seidel sum-of-two-squares obstruction).
