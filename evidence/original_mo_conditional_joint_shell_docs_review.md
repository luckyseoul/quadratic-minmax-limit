# Independent full review: conditional cross joint-shell Gaussian upper

2026-09-05. Reviewer: optimized_profile_docs_gate.
Role: INDEPENDENT REVIEWER of this conditional joint-shell proof.
Final verdict: PASS.

## Frozen artifact and complete-read scope

Source: /tmp/original_mo_conditional_cross_joint_shell_upper.md.
Final complete source: 451 lines, equations (1)--(26).
Final reviewed SHA-256:
64d68bb2feaa59a8049d6bcc42f3ab94c845249c3088fa618916522412d0a68a.

I first read and independently checked the complete 358-line core,
Sections 1--6, at SHA-256:
40334a9436906c3f798b7b09524a53bc73f5a926da464784f6648a17f421cc92.
I then read and checked the complete added Section 7 and its updated
Section 6 transition in the final 451-line source. This complete
review covers both the joint-shell upper and the masked-cross floor,
not merely a summary or a selected formula. The final source hash
was checked again when preparing this receipt.

I am not the author of this conditional joint-shell note. My author
receipt for the separate intrinsic-mu repair is not counted as an
independent review of that separate proof.

## Independent mathematical checks

1. The symmetric cross-edge compression has the stated kernel
   -A tensor A+S_B, with S_B(X)=B X^T B. Its normalization has
   diagonal one and 0<=R<=3I, including the specified n>=2 range.

2. The joint field covariance is PSD with all printed constants.
   The 2-by-2 matrix H=[[q,-c],[-c,-p]] is minus a coordinate-swap
   conjugate of [[p,c],[c,-q]], not a raw adjugate. Its tensor
   compression supplies the required n L^2 bound. Adding vn I
   retains the actual independent Gaussian cushion.

3. The actual cross covariance contains the mixed exchange term
   -(x B y')(x' B y), not a product of same-state energies.
   Direct subtraction gives the nonnegative R quadratic form and
   (d-e)^2 correction in (10). The independent part contributes
   v(n-r_x)(n-r_y). Thus the displayed Gaussian comparison has the
   correct upper direction.

4. All three arcsine reference means have the stated factors:
   internal sums n(n-1), cross sum n^2, and negative cross sign.
   The denominator in the actual cushioned correlation is
   n(k L^2+vD). No independence assertion is needed.

5. The joint Hamming argument retains selection dependence and
   the factor 2n in the spectral-diameter normalization. The
   quadratic polarization and concavity step prove 0<=r<=1.
   The half-normal threshold bound, its exact exponential form,
   and the weaker quadratic deficit all have the printed constants.

6. The two-block feasible constraints use the exact bilinear
   polarization and joint concavity of geometric mean. Their
   flip-cost minimization is a valid improvement of the joint
   one-parameter bound, not an assumed optimizer relationship.

7. The Gaussian maximum bound accounts for both absolute-value
   phases, at most J=(2n^2+1)^3 shells, coefficient norm n, and
   covariance norm at most 3k+v. The conditional shell diamond is
   valid for an actual conditional minimizer with fixed A,-A.

8. Section 7 applies the arbitrary-observable shifted-sign theorem
   explicitly to cross coordinates and the augmented prior with
   arbitrary deterministic internal energy. The mask is stated,
   rather than silently inferred from a whole-observable theorem.

9. The whole-edge Hermite covariance decomposition restricts
   correctly to the principal cross coordinates. The retained
   rank-one term costs at most 2n^2/(sqrt(pi)D)=O(n); the residual
   operator error costs O(n^(5/4)). Combined with the shifted-sign
   comparison this proves the uniform O(n^(16/11)) expected-norm
   approximation.

10. Every thresholded cross matrix remains an admissible competitor
    with the same diagonal blocks. This yields the one-sided
    conditional floor (26), without claiming K is a full order-2n
    optimum or treating pressure optimality as norm optimality.

## Scope limitation and review procedure

No correction was requested. This is a valid finite-dimensional
self-consistency upper and conditional Gaussian floor. It does not
evaluate the upper against 2 sqrt(2) Phi(A), establish the original
all-orders convergence theorem, or apply the unrepaired field at
a different covariance normalization.

This receipt records analytic verification only. No numerical
experiment, finite-order census, source/test mutation, or unchanged
mathematical rerun was used to produce it.

