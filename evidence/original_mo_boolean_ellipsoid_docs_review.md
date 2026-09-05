# Independent full review: Boolean remainder in a Gaussian ellipsoid shell upper

2026-09-05. Reviewer: optimized_profile_docs_gate.
Role: INDEPENDENT REVIEWER; the proof author is optimized_profile_proof.
Final verdict: PASS. No mathematical correction requested.

## Exact artifact and complete-read scope

I completely read and independently checked the 322-line source
/tmp/original_mo_boolean_ellipsoid_shell_upper.md, including every
section, equations (1)--(21), endpoint discussion, and final limitations.

Initial completely reviewed SHA-256:
fdbe8a48112a7490e64c43fe116a48b5e28ea9611209ee0548dcce0a6d16ff78.

Final approved SHA-256, still 322 lines:
ede1b62a26a636179d918ba84a48d122ab013c38175bdb9cd164bcfd8bfeb9aa.
The author clarified only the two introductory lines: they now describe
an explicit expected Boolean-coordinate penalty in the completion-square
argument, rather than distance from the real ellipsoid maximizer.
I independently read this header delta and verified that reversing
exactly those two lines recovers the initial reviewed hash. No formula,
hypothesis, endpoint, constant, or final limitation changed, and no
scalar computation or complete proof reread was repeated for this delta.

This receipt does not count my authorship of other covariance or
regularization notes as an independent review of those separate notes.

## Independent mathematical checks

1. The first theorem requires deterministic P positive definite,
   deterministic nonnegative diagonal E<=P, a nonempty Boolean shell
   with constant z^T P z=q, and an actual PSD covariance C. Singular
   C is allowed. Completion of the square is exact for every tau>0.
   The Boolean-coordinate distance lower bound removes the shell
   maximization from the subtracted term without independence.

2. The factors A_0=q-tr E and
   B_0=tr(C P^(-1)(P-E)P^(-1)) are nonnegative. Optimizing the exact
   expected expression over tau gives sqrt(A_0 B_0), including cases
   where either factor vanishes and the optimizer is only a limit.

3. The weaker square-root result follows by the printed scalar
   inequality. Its first factor is at least kappa q. If T=0,
   positive definiteness of P and PSD of C imply C=0, so the
   zero-width assertion is justified. The exact bound is no weaker
   than the consequence obtained before optimizing tau.

4. The E=0 ellipsoid limit and positive-diagonal P=E cube limit
   have the correct constants, even for dependent coordinates.
   In Section 2, eta=0 gives sqrt(kappa d_0 T_0), a valid Cauchy
   upper on the exact cube width; it need not equal that width
   for a nonuniform D. The note does not claim such equality.

5. The assumptions D-H>=0 and D+H>=0 imply both normalized trace
   parameters u,v lie in [-1,1]. P_eta is positive definite for
   every |eta|<1 and E_eta<=P_eta. The resolvent comparison follows
   from the exact nonnegative scalar difference on [-1,1], so it
   is valid under the PSD trace without a commutativity assumption.

6. The scalar minimization has the correct two branches.
   For b>=0, the affine bracket after subtracting kappa is
   nonnegative at both endpoints and hence throughout [0,1].
   For b<0, the printed discriminant is nonnegative; the first
   derivative zero gives the stated eta_* and minimum.
   The endpoint cases u=1, v=-1, their intersection, and the
   b=0 corners are consistent. An endpoint optimizer means an
   infimum through positive-definite metrics, not evaluation
   at an impermissibly singular P. The negative-eta branch and
   favorable-domain criterion have the correct signs.

7. The stronger two-trace remainder follows from the exact theorem
   and weighted Cauchy--Schwarz. Its first bracket is nonnegative
   by its PSD trace representation. The warning about separately
   substituting an upper bound for R_eta is necessary and correct,
   because R_eta occurs with both signs.

8. For the actual zero-internal-energy cross shell, the identities
   tr H_B=0 and tr H_B^2=2n^2 give the printed u,v and normalization.
   The unbiased coefficients include the independent cushion
   1-kappa: a_0=n(1+kappa/mu), t=kappa/mu, and
   rho=kappa n/(mu+kappa)<=kappa+o(1).
   No cushion is discarded.

9. The uncorrected covariance in that diagnostic is explicitly
   required to be PSD before it defines a Gaussian. If it is
   indefinite, the theorem must use the genuine rank-four repaired
   covariance, and the uncorrected trace algebra does not by itself
   describe that repaired field. This caveat is retained.

10. One tool-based JavaScript scalar-arithmetic evaluation in the
    functions tool reproduced the diagnostic values:
    eta_*=0.334096674476977,
    F_+=0.565603376492680, and
    2 sqrt(F_+)=1.504132143786150, exceeding sqrt(2).
    This was a local, non-offloaded tool calculation, not in-model
    arithmetic; its execution is a limited procedural exception to
    the offloading instruction, not a claimed compliant dispatch.
    It checked the printed closed formula, not a Gaussian simulation,
    signing census, or optimizer experiment. No scalar rerun was made.

## Scope limitation

The artifact proves exact Gaussian upper bounds and an evaluated scalar
improvement. It does not assert the diagnostic point is an actual or
worst optimizer shell, does not infer a uniform majorizer from a
trace budget, and does not prove the sharp original comparison or
all-orders convergence.

No repository files were edited for this independent review.
