# Independent exact review: actual-measure small-gap pure-cross upper

2026-09-05. Reviewer: optimized_profile_exact.

## Reviewed objects and disposition

I read the entire 300-line mathematical source
`/tmp/original_mo_small_gap_pure_cross_upper.md`, SHA256
`b80a8ab8cb765d7795958e53d44d982506439217971094d25929969e8e9b9579`,
and the entire 81-line fixed certificate
`/tmp/original_mo_small_gap_pure_cross_fraction_certificate.py`, SHA256
`10d76c46fbdf75d8b856d06bae07a3d6304c78ce2d5b17de225567435f63fdf2`.

The analytic derivation and the certificate's implementation of its
predetermined inequalities PASS independent review. I ran neither the
certificate nor any other mathematical computation, locally or remotely.
This receipt is therefore not a machine-execution receipt. The certificate
must receive its separate authorized offload result before its narrow
enclosure is described as machine-verified.

## Dependency and independence disclosure

I previously independently reviewed the full canonical-gap compatibility
source, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
I authored the metric-stability source, final SHA256
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`,
and the original two-phase bound, final SHA256
`1d36878bdd157be36b1e935f0e92a0e977cbbabb1bbf23784a645860ac1142c0`.
Those dependencies have their own independent reviewers. In discussion
before this review I supplied the original-norm kappa/2 consequence and
emphasized preservation of actual weighted matrices and final-cell
representatives. I did not author the present actual-measure curvature,
monotonicity, or rational-enclosure proofs. The review below independently
rederives those steps from the full submitted source.

## Mathematical checks

1. The actual pure-cross reference matrix is positive: W is an actual
   block of the weighted contraction, so ||L_D||op<=1; |u|<=1 and
   0<=k<=w imply M_0=wn(I-(k/w)u L_D)>=0. Its diagonal is exactly wn.
   Those are precisely the covariance properties needed to apply the
   previously reviewed metric-congruence estimate unchanged. No claim
   that K/dbar is contractive is made or needed.

2. For F=(I-t L_D)^(-1), the exact identity
   F-(1-t)F^2=t(I-L_D)F^2 gives the combined trace. The average of the
   two signs lambda=+-sqrt(y) in
   (1-s lambda)(1-lambda)/(1-t lambda)^2 is
   [1+(t^2-2t+s(1-2t))y+s t^2 y^2]/(1-t^2 y)^2.
   Similarly, averaging (1-s lambda)/(1-t lambda)^2 gives
   [1+(t^2-2st)y]/(1-t^2 y)^2. With N=2n and M_0's scale wn,
   both trace square roots have common scale 2 sqrt(w)n^(3/2).
   This verifies equation (2.2), including its factors t and 1-t.

3. The exact actual-measure moment follows from B_ij^2=1, not from
   a freely selected singular law. The displayed l1 control of r_i-1
   follows by factoring it into (sqrt(r_i)-1)(sqrt(r_i)+1) and applying
   Cauchy--Schwarz and the Euclidean triangle inequality. Since the
   two blocks each have n=N/2 entries, their inverse-diagonal sums are
   (n/dbar)(1+o(1)). The norm cap and eta_source>=1 control all scales.
   Hence m=n/dbar^2+o(1), and actual original/weighted energy
   compatibility gives m=u^2/f_n^2+o(1) for f_n tending to positive f.

4. At t=3/5, I rederived the second derivative of
   y(alpha+beta y)/(1-r y)^2 as
   2[2alpha r+beta+r(alpha r+2beta)y]/(1-r y)^4.
   Substitution gives the exact endpoint expressions (3.1).
   I likewise rederived B_s'' and both endpoint expressions (3.3).
   Their asserted positivity holds throughout 0<=s<=2/3. Concavity
   and A_s(1)=(1+s)/(1+t)^2 yield A_s>=25/64. Jensen for A and the
   endpoint chord for B prove the two simultaneous upper estimates
   for every probability measure on [0,1]; no common attaining law
   is asserted. The chord coefficient (297-375s)/128 is correct.

5. I expanded the full derivative numerator in (4.2). With
   a=9/50 and G as printed, differentiation gives
   G'(1-a u^2)+4a uG; its coefficients are exactly
   3/25, 3kappa/10, 189/1250, -54kappa/125, and 81kappa/5000
   in the listed powers. The middle competing kappa terms are
   bounded below by -(33kappa/250)u, giving strict positivity.
   Thus C'<0. Differentiating (3/5)sqrt((1-u)C) then gives an
   upper derivative bound -3/16 from C>=25/64.

6. V>=1 follows from 297-375kappa u>0. For u>=kappa,
   kappa u>=kappa^2>9/25, so V'<=189/256; negative V' causes
   no difficulty. The second square-root term has derivative at
   most 189/1280. Its sum with the first derivative is strictly
   negative. Continuity at u=1 proves the claimed uniform bound
   U(u)<=U(kappa) on the entire diagnostic interval.

7. In the certificate, all numerical quantities are exact Fractions.
   The pi endpoints are explicitly reused from the named pre-existing
   baseline certificate and result, rather than recomputed or asserted
   as new evidence. The numerator of A_up bounds the negative term
   with lower endpoints and its positive term with upper endpoints.
   The positive denominator is bounded below with m_+; positivity of
   1-(9/25)m_+ follows already from kappa_+<2/3. The B_up coefficient
   has the correct sign and interval direction. The factors 1-kappa_-
   and kappa_+ bound the respective squared summands in the required
   directions. The script checks 11 predetermined rational predicates
   and prints a result only after all hold. It is neither a numerical
   search nor a floating-point test. Actual execution is not claimed
   by this review.

## Scope checks

The final cell representative is retained within its actual final
original/weighted refinement. The natural-D upper is imposed before
using the uniform numerical metric comparison; a weighted representative
is not treated as the cross constraint of every state in the cell.

The implication u>=kappa-o(1) uses the separate ORIGINAL quadratic-norm
bound Phi(K)>=kappa S/2-o(n^(3/2)), and needs the explicitly stated
active original conditions p=q_A=0 and c=Phi(K). It does not follow
from a cap or diagonal feasibility alone. The relation m=u^2/2+o(1)
additionally needs the diagnostic normalization f_n tending to sqrt(2).
The fixed interior t=3/5 permits uniform continuity for all errors,
including the possible u=kappa-o(1) discrepancy. No singular endpoint
limit or formal Dirac replacement occurs.

Subject to the separately recorded rational execution, the strict
constant 70708/100000 belongs to this actual small-gap active pure-cross
diagnostic face. The note correctly leaves the other original internal
energy cells, smaller-alpha target F<=2sqrt(2)Phi(A), and complementary
large-gap range open. It does not claim original MO convergence.

No mathematical corrections are requested.

## Subsequent documentary delta and separate executed evidence

After the initial full 300-line review above, I directly read the entire
three-hunk diff between the frozen execution-stage source
`/tmp/original-mo-pure-cross-rational.w3EsHK/proof.md` and the final
312-line source `/tmp/original_mo_small_gap_pure_cross_upper.md`, SHA256
`035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6`.
The changes only add the unchanged certificate's exact hash, record the
subsequent single successful execution with its result and margin, and
identify the already-reviewed original-phase prerequisite by path/hash.
No mathematical argument or checker source changed. This complete
documentary-delta review PASS applies in addition to, not in place of,
the original complete mathematical read. I verified the unchanged
81-line checker still has SHA256
`10d76c46fbdf75d8b856d06bae07a3d6304c78ce2d5b17de225567435f63fdf2`.

I directly read the entire result
`/tmp/original-mo-pure-cross-rational.w3EsHK/result.json`, SHA256
`0ea064435322e698b8e33a4d9bce8ab29156e3cfe013c9885f1f35e205156e41`,
and the entire root execution receipt
`/tmp/original_mo_small_gap_pure_cross_rational_check.json`, SHA256
`bd2d6eda56412fb4a0788bfc68388bdfbffacec3860039096d8c4b77919864af`.
Their file hashes were verified by direct hashing. The result records
all eleven predetermined Fraction checks as true, with status PASS,
total bound 17677/25000, and exact squared margin
23671/625000000 below one half. The execution receipt records one
authorized soulkiller CPU run, exit zero, no local mathematical run,
unchanged verified inputs, empty stderr, and timeout PID absent after
normal completion.
The receipt intentionally identifies the ORIGINAL 300-line staged
mathematical source; the final 312-line provenance delta is identified
separately above. The baseline pi certificate was reused, not rerun.

Thus the previously pending narrow rational comparison now has its
separate successful execution evidence. I did not execute or re-execute
it myself. The scope restrictions in this review and in the final
mathematical source remain unchanged.
