# Directed half-cut literature check (2026-09-02)

**Verdict.**  I found no published theorem, preprint, or implementation whose
hypotheses imply

\[
 \min_{S\ \mathrm{tournament}}\max_{U\subseteq[n]}
 \Phi(A^{E_S^+(U)})
 \leq \sqrt2\,\Phi(A)+o(n^{3/2})                         \tag{L}
\]

for an arbitrary signed complete graph `A` (or even for a globally
`Phi`-minimal one).  The searches below were deliberately restricted to the
new directed-half-cut formulation.  This is a negative search result, not a
claim that (L) is false or that no such theorem exists under different
terminology.

The exact quantifier pattern that was screened was

\[
 \exists S\quad\forall U\quad\forall y\in\{\pm1\}^n:
 \left|Q_{A^{E_S^+(U)}}(y)\right|
 \leq \sqrt2\,\Phi(A)+o(n^{3/2}).                       \tag{Q}
\]

No finite-order computation was run.

## 1. Closest sources and the decisive hypothesis mismatches

1. **Tournament switching is full-cut reversal.**  Babai and Cameron define
   switching a tournament at `Y` by reversing *all* arcs between `Y` and its
   complement; their results characterize automorphism groups and enumerate
   switching classes.  The target `E_S^+(U)` contains only the arcs that
   originally leave `U`, so it is a strict half of a mixed-direction cut in
   general.  There are no edge weights, Boolean test vector, or norm estimate
   in their theorem.

   L. Babai and P. J. Cameron, “Automorphisms and Enumeration of Switching
   Classes of Tournaments,” *Electronic Journal of Combinatorics* **7**
   (2000), R38. DOI: [10.37236/1516](https://doi.org/10.37236/1516).

2. **Seidel/two-graph switching is also full-cut switching.**  In matrix
   language it is `A -> D A D`, hence it preserves `Phi(A)` exactly.  A
   directed half-cut neighbor generally changes only a proper subset of
   `delta(U)` and leaves the Seidel switching class.  The current 2026 paper
   on identity switches studies switches producing an isomorphic graph; it
   supplies no estimate for partial flips.

   S. V. Gervacio, “On identity Seidel switches,”
   [arXiv:2601.04530](https://arxiv.org/abs/2601.04530) (2026).

3. **Cycle/cocycle reversal reverses an already directed cut in full.**
   Aichholzer et al. study flips that reverse every edge of a minimal directed
   cut (and cycle reversals for prescribed-outdegree orientations), with flip
   distance and distributive-lattice conclusions.  Backman's generalized
   cycle--cocycle system likewise concerns divisor equivalence and directed
   path/cut reversals.  Neither theorem first selects the outward arcs of an
   arbitrary mixed cut, and neither has signed edge weights or a cube norm.

   O. Aichholzer, J. Cardinal, T. Huynh, K. Knauer, T. Mütze, R. Steiner,
   and B. Vogtenhuber, “Flip Distances Between Graph Orientations,”
   *Algorithmica* **83** (2021), 116–143. DOI:
   [10.1007/s00453-020-00751-1](https://doi.org/10.1007/s00453-020-00751-1);
   [arXiv:1902.06103](https://arxiv.org/abs/1902.06103).

   S. Backman, “Riemann--Roch Theory for Graph Orientations,”
   [arXiv:1401.3309](https://arxiv.org/abs/1401.3309) (2014).

4. **Balanced-orientation theorems control connectivity or unsigned
   cardinality, not signed energies.**  Nash-Williams' theorem gives an
   orientation with
   `lambda_D(u,v) >= floor(lambda_G(u,v)/2)` (and vertex imbalance at most
   one).  Király--Szigeti simultaneously obtain this property for pairwise
   edge-disjoint subgraphs.  The newer `k`-cut-balanced condition is

   \[
     |\delta_D^+(U)|\geq |\delta_G(U)|/k\quad\hbox{for every }U.
   \]

   These conclusions count arcs and are invariant under neither the weights
   `A_ij y_i y_j` nor the simultaneous family of all `y`.  In fact, weighted
   and extra-constrained well-balanced-orientation variants are generally
   NP-hard; that is not a counterexample to (L), but it blocks a direct
   algorithmic transfer.

   C. St. J. A. Nash-Williams, “On Orientations, Connectivity and
   Odd-Vertex-Pairings in Finite Graphs,” *Canadian Journal of Mathematics*
   **12** (1960), 555–567. DOI:
   [10.4153/CJM-1960-049-6](https://doi.org/10.4153/CJM-1960-049-6).

   Z. Király and Z. Szigeti, “Simultaneous well-balanced orientations of
   graphs,” *Journal of Combinatorial Theory, Series B* **96** (2006),
   684–692. DOI:
   [10.1016/j.jctb.2006.01.002](https://doi.org/10.1016/j.jctb.2006.01.002).
   Its simultaneous theorem assumes the designated subgraphs are pairwise
   edge-disjoint.

   K. Chandrasekaran, S. Liu, and R. Ravi, “Minimum Cost Nowhere-Zero Flows
   and Cut-Balanced Orientations,” *ICALP 2025*, LIPIcs 334, Article 46.
   DOI: [10.4230/LIPIcs.ICALP.2025.46](https://doi.org/10.4230/LIPIcs.ICALP.2025.46);
   [arXiv:2504.18767](https://arxiv.org/abs/2504.18767).

   A. Bernáth, S. Iwata, T. Király, Z. Király, and Z. Szigeti, “Recent
   results on well-balanced orientations,” *Discrete Optimization* **5**
   (2008), 663–676. DOI:
   [10.1016/j.disopt.2008.03.001](https://doi.org/10.1016/j.disopt.2008.03.001).

5. **“Oriented discrepancy” currently means a selected spanning subgraph.**
   Gishboliner--Krivelevich--Michaeli prove that, under minimum-degree or
   random-graph hypotheses, every orientation contains a Hamilton cycle with
   many edges pointing in one cyclic direction.  It is an
   `forall orientation, exists cycle` cardinality theorem.  (Q) instead asks
   for one orientation controlling every cut and every signed Boolean test.

   L. Gishboliner, M. Krivelevich, and P. Michaeli, “Oriented discrepancy of
   Hamilton cycles,” *Journal of Graph Theory* **103** (2023), 780–792. DOI:
   [10.1002/jgt.22947](https://doi.org/10.1002/jgt.22947);
   [arXiv:2203.07148](https://arxiv.org/abs/2203.07148).

6. **Gain-graph/four-way switching matches the conjugation, but not the
   norm being bounded.**  The exact local identity already proved in this
   repository is

   \[
     A^{E_S^+(U)}=\operatorname{Re}(Z_U^*(A+iR)Z_U),
     \qquad Z_U\in\operatorname{diag}\{1,-i\}^n.        \tag{H}
   \]

   Mohar's four-way switching and general unit-gain switching are diagonal
   unitary similarities, so they preserve the spectrum of the *full*
   Hermitian matrix.  They do not bound
   `max_y |y^T Re(Z^* H Z)y|`, do not optimize the imaginary completion `R`,
   and do not compare that maximum with `Phi(Re H)`.  Thus (H) is a valid
   dictionary, but the published switching theorems stop immediately before
   the projection/cube-norm step needed for (L).

   B. Mohar, “Hermitian adjacency spectrum and switching equivalence of
   mixed graphs,” *Linear Algebra and its Applications* **489** (2016),
   324–340. DOI:
   [10.1016/j.laa.2015.10.018](https://doi.org/10.1016/j.laa.2015.10.018);
   [arXiv:1505.03373](https://arxiv.org/abs/1505.03373).

   M. Kadyan and B. Bhattacharjya, “Switching equivalence of Hermitian
   adjacency matrices of mixed graphs,”
   [arXiv:2103.13632](https://arxiv.org/abs/2103.13632) (2021).  Its
   hypothesis is `B=D^{-1}AD` for a diagonal unit-gain matrix `D`, and its
   conclusions concern switching classes and cospectrality.

## 2. What generic discrepancy gives (and why it does not close the gap)

For fixed `(U,y)`, the orientation-dependent part of
`Q_{A^{E_S^+(U)}}(y)` is a linear form in the
`N=binom(n,2)` independent orientation signs.  There are at most `4^n`
such `(U,y)` constraints.  More explicitly, if `I_{U,y}` is the `y`-energy
of edges not crossing `U`, then

\[
 Q_{A^{E_S^+(U)}}(y)=I_{U,y}
       +\sum_{e\in\delta(U)}c_{U,y,e}s_e,
 \qquad c_{U,y,e}\in\{\pm1\},                         \tag{D1}
\]

for a fixed convention for the independent tournament signs `s_e`.
Moreover

\[
 I_{U,y}={Q_A(y)+Q_A(y^U)\over2},\qquad |I_{U,y}|\leq\Phi(A), \tag{D2}
\]

where `y^U` flips the coordinates in `U`.  Applying the general Spencer
bound to the linear terms in (D1),

\[
 \operatorname{disc}(m,N)=O\!\left(
 \sqrt{N\log(2m/N)}\right),\qquad m\geq N,
\]

gives at best

\[
 \min_S\max_{U,y}|Q_{A^{E_S^+(U)}}(y)|
 \leq \Phi(A)+O(n^{3/2}).                              \tag{D3}
\]

This neither gives `o(n^(3/2))` nor the relative sharp constant `sqrt(2)`
in (L).  Banaszczyk's general vector-balancing theorem has the same issue
when specialized to an `l_infinity` body with exponentially many
coordinates.

J. Spencer, “Six Standard Deviations Suffice,” *Transactions of the American
Mathematical Society* **289** (1985), 679–706. DOI:
[10.2307/2000258](https://doi.org/10.2307/2000258).

W. Banaszczyk, “Balancing vectors and Gaussian measures of n-dimensional
convex bodies,” *Random Structures & Algorithms* **12** (1998), 351–360.
DOI:
[10.1002/(SICI)1098-2418(199807)12:4%3C351::AID-RSA3%3E3.0.CO;2-S](https://doi.org/10.1002/(SICI)1098-2418(199807)12:4%3C351::AID-RSA3%3E3.0.CO;2-S).

There is also a simple `Omega(n^(3/2))` lower bound for the ordinary
discrepancy of the system of all cuts of `K_n`: split the vertices into two
equal parts, choose a random subset on one side, retain linearly many columns
with signed sum of order `sqrt(n)`, and combine one of three cuts.  This does
**not** refute (L), whose constraints have extra oriented-incidence structure
and whose right side is itself of order `n^(3/2)`; it does show that an
unstructured “balance all cuts with lower-order error” lemma cannot be the
missing ingredient.

## 3. Repository searches

As of 2026-09-02, exact GitHub code searches returned zero results for
`"directed half-cut" graph`, `"one-sided cut flip" graph`, and
`"cut reversal" tournament discrepancy`.  The query
`"Seidel switching" tournament` returned seven code hits; the only
mathematical-research hit relevant to (L) was this repository itself, while
the others were unrelated notes/metadata.

Exact Zenodo API phrase searches returned zero records for
`"directed half-cut"`, `"partial cut reversal"`, `"one-sided cut flip"`,
`"Seidel switching"`, `"cut-balanced orientation"`, and
`"directed cut reversal"`.  The sole hit for `"orientation discrepancy"`
was unrelated to graph theory.  These database checks support only the
terminology-level negative result; they are not an exhaustive proof about
all deposited mathematics.

## 4. Actionable conclusion

There is no citation-ready black-box theorem among these literatures that
closes (L).  The only genuinely exact bridge is the Hermitian gain-switching
identity (H).  A successful proof still needs a new **completion inequality**:
choose a skew signing `R` for a given symmetric signing `A` so that the
fourth-phase cube norm of `A+iR` is at most `sqrt(2)` times the real cube norm
of `A`, up to Dini-lower-order error.  Neither Seidel switching, tournament
pushes, cut-balanced orientations, cycle/cocycle reversals, nor current
oriented-discrepancy theorems supply that inequality.
