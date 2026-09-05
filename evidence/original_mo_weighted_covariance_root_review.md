# Root review: weighted covariance, actual shells, and route constraints

2026-09-05. Verdict: PASS for the four complete analytic proofs below.
This is a full mathematical review with disclosed collaboration, not a
claim that arithmetic checks establish the all-orders statements.

## Reviewed sources and provenance

- Cubic alignment: final 338 lines, SHA256
  `60037f67234fbca8c17ee90bf52c7f4346b24e5f18eb5f2c922ebbd2d9382c2a`.
  Root read the complete initial source, checked the Section 6 extension,
  and then read every line of the complete final source. The exact worker
  authored the new cubic theorem; root contributed its connection to the
  diagnostic barrier and requested the explicit variable-u/limit-order
  derivation. The docs worker supplied a separate independent full review.
- Fixed-cap deflation: complete 342 lines, SHA256
  `22febfa722afb3e18878f23f8e140895da90a3eb41fe0179356b08232d44f27a`.
  The proof worker authored the construction following root's rate
  question. Root read all 342 lines; the exact worker supplied an
  independent complete 109-line review, also read fully by root.
- Diagonal-majorizer cross covariance: complete initial 381 lines,
  then the full revised Section 5 including both changed locations;
  final 384 lines, SHA256
  `0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
  Root proposed the full-D route and the separable padded covariance.
  The proof worker completed the all-Hermite analysis and reduction.
  The docs worker supplied the independent full review. The author's
  required epsilon<=1/2 series-domain clarification is incorporated;
  smaller bounded orders use the already separate direct argument.
- Weighted-shell upper: complete final 381 lines, SHA256
  `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
  The exact worker authored the field and metric proof; root contributed
  the real-parameter binning direction and completely checked the
  resulting proof. Root's read included the final covariance-dependency
  hash. The docs worker supplied the separate independent full review.

Exact imported bytes are checked against each frozen source.
No canonical proof is edited during import.

## Cubic threshold and diagnostic scope

Matched scalar-optimal frames give both CV=qU and C^TU=qV, their common
frame operator, and the simultaneous top-frame completion inequality.
Root checked the four-tensor flattening: each of the three fourth-moment
arrays has second-moment sum p^2 s, so the squared norm is at least
4p^2 s^2. The contraction upper on the same quantity gives
j_3>=-1+2s^2/mu, with no actual-rank inference.

The cubic clipping coefficient has square 2phi^2/3. Keeping that positive
correction alongside the exact degree-two cancellation proves the stated
completion lower. The quartic penalty maximization and strict rational
margin prove q>=12/5 implies Gamma>283/200, only under tau(C)=p q.
Gamma remains an upper certificate in its actual-sign application.
Section 6 separately stipulates leading energy sqrt(2); its general
trace formula permits variable u, and its repaired-field argument fixes
eta before the n limit. This clears that specified weak-Dirac diagnostic,
not arbitrary actual-source or nonsymmetric attainability.

Root completely read the 74-line checker and wrapper before its single
remote execution, then all five exact-fraction results and the execution
receipt. All five passed. The earlier 28 clipping comparisons were reused
and not rerun. These are arithmetic supplements to the analytic proof.

## Actual fixed-cap repair lower

Root checked the exact finite Hadamard modification, localized top-vector
overlap and untouched bulk subspaces. The SDP amplification bound controls
the ACTUAL signing cap despite its operator outlier. Both tensor signs
contain the stated commuting spectral subcovariances. Gaussian convexity
gives a lower, not upper, norm comparison. The sparse cube witness, both
symmetric covariance copies, factor sqrt(2), and diagonal-removal loss
all have the required directions and constants.

The amplitude is fixed before the threshold, then the order grows.
Every cap C>1/2 is covered. The conclusion excludes uniform exponents
greater than 1/2 for this repair under that fixed cap; it does not cover
C=1/2, exact minimizers, adaptive slack, or the full coupled cross repair.

## Same-source weighted covariance and upper

Root checked SDP symmetrization, attainment, strictly positive D, the
row-square lower bounds, symmetric-edge compression and independent
padding. The covariance is PSD before a Gaussian is used. The weighted
four-cycle decomposition includes every even degree with its correct
series domain, exact diagonal cancellation, and all retained profiles.
Their actual Boolean norm cost is O(n), not inferred from low rank.
Operator-error comparison, shifted-sign universality, and removal of
padding preserve the original source, drift and conditional objective.

For the field theorem root independently expanded the general midpoint
increment identity, including -Delta p Delta q+(Delta c)^2. An actual
representative covariance is PSD by tensor compression. The independent
state noise restores increment domination with the stated constants;
weighted bins and original integer triples have at most (2n^2+1)^6
nonempty cells. Cell concentration needs no independence among suprema.

The same D majorizes the original cross matrix H_B. Its metric radius
therefore uses original c, while cyclic trace identities retain actual
weighted p_D,q_D,c_D. Root checked the complete Boolean completion-square
argument, its affine consequence, the r_L r_R feedback coefficient, and
the explicit Delta_B and internal Delta_A discrepancies. Those terms
are not discarded or asserted small. Equation (4.6) is a valid but not
sharply evaluated upper on the actual cells.

## Scope of publication

The original MO problem is OPEN. Neither convergence, nonconvergence,
the value 1/2, nor the sharp conditional dyadic comparison is proved.
The current next target is an actual-coupled-cell upper evaluation using
source/conditional optimality or a proved compatibility estimate.
Unfinished follow-on canonical-gap work is outside this milestone.
No src modules, tests or global proof predicates are changed.
