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

Secondary: Paley construction (Wikipedia, quoting the 1933 Hadamard-conjecture sentence); conference matrix (Wikipedia, van Lint–Seidel sum-of-two-squares obstruction).
