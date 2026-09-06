# Root contributing review: actual original-source strict-gain package

2026-09-06. PASS for the stated conditional theorem and same-source
transfer, not for global convergence or an all-profile assertion.

## Complete proof reads and frozen versions

I read the entire initial 611-line main proof at SHA256
0759910b88aa9586e2552de64c2dac96873893f110fa4509bbeffbdc85f3c663.
The sole final change explicitly says that the paired corollary's D
is diagonal. I then read that final corollary and the entire postlude
and verified the final 612-line source hash:
7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298.
This records the actual full-source-plus-final-delta review pattern,
not a claimed second end-to-end reread of all 612 lines.

I also completely read both supporting proofs:

- Scalar local-update bound, 209 lines, SHA256
  7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155.
- Actual internal-law transfer, 141 lines, SHA256
  f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8.

The full 411-line Gaussian predecessor was reread in this turn at
b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20.
The main proof rederives its new normalization, variance, distinguished-
coordinate, and Boolean-update arguments explicitly. The previously
verified pi enclosure is reused, not rerun.

All seven final author/reviewer receipts for these three proofs were
read completely: main author111/exact227/docs280, scalar author86/
exact104, and transfer author50/docs146. Their source hashes match.
No reviewer requests a correction. Their contribution boundaries are
retained, not replaced by a blanket claim of complete independence.

## Main analytic checks

The actual operator cap and weak three-atom law imply the two spectral-
projector Frobenius approximations. Complete sign squares fix diag(A^2)
at n-1, and zero diagonal fixes the signed-projector diagonal error.
Good-coordinate normalization of the positive spectral projector plus
independent unit-variance bad coordinates gives an actual correlation
R on ALL original coordinates, with bounded operator norm and
R-R0=o_F(sqrt(n)). The original matrix and energy are not trimmed.

Schur powers remain uniformly bounded in operator norm. For every odd
q>=3, complete off-diagonal entries bound tr(A R^{circ q}) by O(n).
This yields the original phase baseline 5kappa/8. The matrix identity
A^2/lambda^2=rho R-A/lambda+o_F(sqrt(n)) then retains the even-power
diagonal contribution in every higher chaos. The errors are uniform in
q, so summing the infinite nonnegative Hermite coefficients gives mean
higher variance at least 1-kappa-o(1), not half that value.

I specifically checked the first-chaos variance via the trace of
(A-lambda I)^2 R. Frobenius approximation of AR alone would not
justify it for singular R. The trace-of-square comparison does, and
implies the required empirical row covariance convergence.

The generalized coefficient estimates permit the single zero diagonal
entry. The distinguished-coordinate contraction is O(1/n) in squared
norm, and the earlier delocalized contractions control all remaining
nonconstant mixed terms. The finite-chaos characteristic equation,
uniform L2 tail, and compact covariance subsequences prove a joint
limit for ONE Gaussian input coordinate and ONE local field. The sign-
disagreement boundary has zero limiting probability because both
Gaussian marginals have positive variance, even at correlation one.
There is no growing-dimensional CLT or covariance assumption on sign(AX).

The actual independent Bernoulli update has exactly the epsilon linear
term and epsilon^2 quadratic term: diagonal exceptions vanish since
A_ii=0. A mismatch contributes four to the squared update displacement,
giving the penalty 2epsilon^2 C_n p_n. The ACTUAL cap 5/3 is used,
not the limiting nonzero atom 5/4. Every updated state is Boolean on
the same source; no bilinear norm or polarization factor is substituted.

With fixed admissible epsilon=1/10, the heterogeneous-variance chord
and mismatch bounds give gain at least 16/3125 in the limit. The
reused kappa interval gives
5kappa/8+16/3125>2/5+3/1100. The unconstrained scalar optimizer is
not needed. The optional scalar calculations are correctly qualified.

## Same-source transfer and consequence

The common good index labels leave q/n tending to one and retain a
complete principal source A_J on both halves. Principal interlacing
and congruence give the normalized second-moment error 2b+9epsilon.
Literal complete sign squares recover dbar/sqrt(q) tending to 5/3
from the FULL weighted internal law, without a separate trace cap.

The same comparison transfers the whole empirical law to atoms zero
and plus or minus 5/4 and gives limsup ||A_J||/sqrt(q)<=5/3.
Unbiased extension gives Phi(A_J)<=Phi(A); feasibility and the recovered
scale justify the normalized o(1) comparison. This is only an original-
norm lower-bound device, never a replacement paired field or covariance.

The new lower therefore excludes the specified ACTUAL near-scalar
internal-law regime at alpha tending to 2/5. The older FORMAL certificate
boundary remains true for its listed relaxation: it never supplied an
actual signing or active-state realization. No other profile is silently
excluded and no optimizer is assumed to have small diagonal dispersion.

## Contribution disclosure and stopping scope

I proposed the positive-projector/local-update strategy, the corrected
operator cap, the fixed-probability step, and the trace-of-square audit.
The proof worker derived and authored the robust main proof. The docs
worker derived scalar bounds, and the exact worker authored the transfer.
My review is therefore contributing. The exact worker independently
checked the NEW main-core steps with its older Gaussian-prerequisite and
transfer contributions disclosed; docs independently checked the main
non-scalar links while disclosing scalar authorship. The scalar and
transfer proofs each also have an independent new-derivation review.

The companion exact-rank barrier was read separately in full. Its
complete replicated-row family shows why exact parity rank obstruction
has no generic asymptotic stability upgrade, but its source norm tends
to 1/sqrt(2), not 2/5. It is not a counterexample to the source theorem
or a separate claimed milestone; its role is to prevent that false route
from being mistaken for the actual robust argument above.

No mathematical computation, scan, checker, numerical integration,
solver, or signing search was run. Original convergence, the limit
value, and the remaining all-profile/all-cell implications remain OPEN.
