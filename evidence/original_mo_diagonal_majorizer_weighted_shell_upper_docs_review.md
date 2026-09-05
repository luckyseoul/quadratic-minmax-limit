# Independent full review: weighted fields and original cross-energy metric

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact sources and independence

The reviewer directly read every line of the final 381-line source
`/tmp/original_mo_diagonal_majorizer_weighted_shell_upper.md`, SHA-256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.
This was a full read of the version with the corrected final covariance
dependency hash, not a review of the earlier dependency reference.
The reviewer did not contribute to this theorem, its binning argument,
or its derivation, and made no source edits.

The covariance prerequisite had just received the reviewer's complete
381-line initial read and full final Section 5 delta review. Its final
384-line SHA-256 is
`0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
The reviewer additionally read all 322 lines of the current
`NOTE_2026-09-05_BOOLEAN_ELLIPSOID_SHELL_UPPER.md`, SHA-256
`ede1b62a26a636179d918ba84a48d122ab013c38175bdb9cd164bcfd8bfeb9aa`.
The exact Boolean remainder and affine/two-trace hypotheses were
checked directly against the new field and metric.

## Covariance, energies, and the exact positive field

The source keeps the two triples separate: p,q,c are the original
integer quadratic and cross energies, while p_D,q_D,c_D use UAU,
VAV, and UBV. Compression of the normalized contraction bounds all
three weighted quantities by n in absolute value. The two covariance
identities (1.1)--(1.2) have the correct cross orientation and signs.

The padding cost (1.5) follows from the variance of each Boolean
pairing and the finite-state Gaussian maximum. The inverse-trace and
two-part product bounds give its O(n) rate. With fixed offsets the
same cost follows by objective Lipschitzness; no padding is removed
from the actual unit-variance shifted-sign law.

At an attained weighted triple, compressing T along the two normalized
Boolean block vectors gives the stated two-by-two matrix with norm
at most n. Negative swap conjugation gives H_theta with the same norm.
The specified compression of H_theta tensor T produces exactly the
nonconstant block in (2.1), including the positive p_D A_R block.
Its norm is at most n. Hence the field eigenvalues lie between
(w-k)n and (w+k)n, proving 0<=M_theta<=2wnI since 0<=k<=w.
This proof uses actual attainment of the representative parameters.

## Complete increment identity and real-parameter binning

The reviewer independently expanded the two increment variances.
For unequal weighted parameters, writing the two-state averages as
bar p, bar q, bar c gives, after dividing the excess by two,

    w(n-r_x)(n-r_y)
    +k[(bar p-a)(bar q-b)-bar c^2+bar c(d+e)-de
                                 -Delta p Delta q/4+(Delta c)^2/4],

where a=x-transpose A_L x-prime and b=y-transpose A_R y-prime.
Expanding the R_0 quadratic form on delta x tensor delta y yields
exactly (3.1). In particular the exchange term is minus de and the
parameter corrections are minus Delta p Delta q plus (Delta c)^2,
with the displayed factor one-quarter. On an exact weighted shell
these corrections vanish, giving (2.4) and valid increment domination.

Gaussian increment comparison applies with identical deterministic
offsets and permits singular covariance. The optional positive padding
of M_theta contributes exactly the increment needed to replace R_0
by R_D in the same identity. It is not used to justify an indefinite
unperturbed field; the base field was already positive semidefinite.

Inside a real weighted cell, each parameter difference has magnitude
at most delta. Thus the midpoint identity is bounded below by
-k delta-squared/4. Replacing its affine covariance by that of an
actual representative changes the operator by at most 2k delta.
Because the Boolean increment has squared norm at most 8n, its
half-variance cost is at most 8kn delta. This proves the exact
variance choice for the auxiliary independent state noises.

For delta<=4n that variance is at most 9kn delta. There are at most
2^(2n) states, giving exactly 6n sqrt(k delta log(2)) in (3.3).
Equal-state increments are zero, so adding independent noises creates
no exceptional diagonal case. The midpoint covariance is used only
algebraically; the comparison field is always the attained positive
representative field.

Taking delta=1/n and refining by the original integer triple gives
the stated bound (2n-squared+1)^6 on nonempty actual cells. Exact
weighted values are not assumed to have polynomial multiplicity.
The representatives in the final bound belong to the actual cells;
the original triple and weighted representative are not independently
chosen values from a relaxed box.

Each original cell supremum is Lipschitz in an underlying standard
Gaussian with constant at most n sqrt(2k+v). Gaussian concentration
and the exponential maximum bound prove (3.5) without independence
of cell suprema. Arbitrary deterministic cell offsets do not change
the Lipschitz constant. Augmenting both signs for an absolute maximum
changes log(m) to log(2m). The resulting O(n sqrt(log n)) selection
cost and O(sqrt(n)) within-cell cost are both o(n^(3/2)).

## Same-diagonal original cross constraint and the two field traces

Conjugation by J preserves the diagonal D and reverses only the cross
block. Averaging the two stated PSD constraints therefore proves
D-H_B>=0 and D+H_B>=0 with the ORIGINAL unweighted H_B.
Consequently L_D is a contraction, P_eta is positive definite, and
E_eta=(1-|eta|)D satisfies 0<=E_eta<=P_eta for every |eta|<1.

The Boolean radius is exactly d_0-2eta c. It involves the original
cross energy c, not the weighted c_D. The same constraints give
|2c/d_0|<=1. Constancy of the original c on every refined cell is
precisely the shell hypothesis required by the Boolean remainder
theorem, even though the weighted values are only binned.

Congruencing the field covariance gives the displayed Mhat_theta:
its nonconstant term is kD-inverse cal H_theta D-inverse. Cyclic
trace identities give T_eta=tr(M_theta P_eta-inverse) and
R_eta=tr(D P_eta-inverse M_theta P_eta-inverse) exactly, with
F_eta-squared in the second trace. No commutation or simultaneous
diagonalization assumption is made.

The completion-square proof quoted in Section 4 retains the correct
coordinatewise Boolean penalty. The radius-minus-trace term becomes
d_0(|eta|-eta u), and the second nonnegative factor becomes
T_eta-(1-|eta|)R_eta. Cauchy--Schwarz in the diagonal remainder gives
exactly (4.5), including the full factor (1-|eta|) outside sqrt(R_eta).
Both zero-factor cases are valid by the one-sided optimization limit.
Each metric is used at a fixed finite eta, so its infimum requires
no uncontrolled endpoint interchange.

For a fixed cell the centered Gaussian symmetry equates the expected
positive and negative field suprema. This justifies the |a_j| term
after absolute-value augmentation. Combining the cell comparison,
concentration, and field width bound yields precisely (4.6), with
the displayed constants. The original drift is (p_j-q_j)/2+s c_j;
its internal source has not been replaced by weighted energies.

## First feedback and the explicit unresolved compatibility

The field trace at eta zero is T_0=wn tr D-inverse. Its internal
contributions vanish because A has zero diagonal. In the trace against
L_D only the two off-diagonal blocks survive. Their products contain
B_ij-squared=1 and give exactly
-2k c_D (sum d_left^(-3/2))(sum d_right^(-3/2)), proving (5.1).
Positivity of Mhat_theta and contraction of L_D give |v_0|<=1.

The affine ellipsoid bound (5.2) follows from the prior theorem's
resolvent upper inequality. It is correctly described as weaker than
the exact two-trace bound, not an equality replacing that bound.
The w=0 case is separated before dividing by T_0.

The three compatibility differences in (5.3) give the exact identities
(5.4). Substitution into (5.1), with u=2c/d_0, yields (5.5) with the
stated lambda_D and the correct sign and coefficient of Delta_B.
This equation does not identify c_D from the original cross shell.

The diagonal blocks of the resolvent are the two stated inverses of
I-eta-squared W_D W_D-transpose and its reversed product. The internal
contributions to T_eta are precisely (5.6); those for R_eta use the
corresponding blocks of F_eta-squared. Their zero at eta=0 follows
from A's zero diagonal, not from an assumed sign or cancellation at
other eta. In particular original p=q=0 does not imply p_D=q_D=0.

The trace cap is not asserted to control any of these compatibility
terms. The final theorem is a genuine, but still unevaluated, upper
bound on actual coupled cells. It supplies neither an impossibility
theorem for their evaluation nor original MO convergence. This scope
is consistent with the covariance reduction and conditional floor.

No mathematical computation, certificate execution, signing census,
simulation, optimization, solver, or new test was run by this reviewer.
