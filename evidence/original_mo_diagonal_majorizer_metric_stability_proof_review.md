# Proof-worker review: actual-diagonal metric stability

2026-09-05. Full analytic PASS on the final 252-line source. No source
correction, mathematical computation, or numerical evaluation was needed.

Reviewed source:
`/tmp/original_mo_diagonal_majorizer_metric_stability.md`, SHA256
`ab473024c6ec7f2c87377c48bdf58a159236dea954f68df30dd6a32716875c1a`.

The reviewer read the entire final source and rechecked the relevant
field/cell/ellipsoid dependencies in the frozen weighted-shell theorem,
SHA256
`9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.

## Dispersion and actual energies

Both pointwise square-root inequalities follow by writing t=r^2:
their ratios to (t-1)^2/t are respectively r^2/(r+1)^2 and
1/(r+1)^2. Thus neither requires a maximum-diagonal bound.
The symmetric error factorization and contraction of T give the
uniform full-cube bound exactly as claimed. The error is zero-diagonal.
Extending a principal-block state by zero, then averaging over unused
Boolean coordinates, justifies the principal-block bound. Opposite
second-block signs justify the sharper cross bound without needing the
two diagonal error blocks to be negatives of one another.

Consequently the three energy inequalities and |u_D-u|<=2sqrt(delta)
are valid for actual states. The two radii separately belong to [-1,1]
by D+-H_B>=0 and the actual contraction L_D. Their validity does not
depend on small dispersion or on p=q=0.

## Noncommutative covariance congruence

The constant diagonal of M gives all three Frobenius identities used
in Section 3. The expansion of V_0MV_0-M has the correct order of
factors. Applying the nuclear/Frobenius product bound to those factors
gives exactly E_delta, with no unproved commutation involving D or M.

Cyclic trace identities give dbar T_eta and dbar R_eta with respectively
F_eta and F_eta^2 against V_0MV_0. The resolvent norm bounds therefore
give the two separate trace errors. For the cancellation trace the
spectral function is (t-b)/t^2 for t>=b, which is nonnegative and at
most 1/(4b). This proves the stronger combined error and the stated
upper bound on its reference trace.

## Complete upper-expression stability

Rewriting the valid natural-D upper with outer factor sqrt(N) leaves
dbar times the natural traces inside. The coefficient difference is
at most 2a sqrt(delta), while each coefficient is at most two.
Separating the product error using the reference combined trace gives
E_delta/(2b)+a wN^2 sqrt(delta)/(4b), as written.

For delta<=1 this is at most wN^2 sqrt(delta)/b, since its scalar
factor is [1+sqrt(1+delta)+a]/4<1. Taking square-root differences
gives the first term with the recorded constant one. The other
square-root difference is at most sqrt(kappa E_delta). Their sum is
strictly below the stated safe constant three after allowing b<=1.
The w=0 and delta=0 cases are correctly treated by direct equality.

This derivation retains the positive actual M and contractive actual
L_D throughout. It never introduces an indefinite scalar-rescaled
source covariance. It covers all attained internal energies.

## Cells, optimization, and scope

The representative is explicitly chosen within the final refined
cell. Its original c is consequently the exact original cross value
throughout that cell, so B_D is a genuine field upper there. The note
correctly treats B_flat as a compared numerical functional rather than
as an exact scalar-I ellipsoid constraint for every state of a weighted
bin. This distinction is essential and is preserved in (5.1).

The uniform pointwise stability bound over |eta|<=1-b_0 survives taking
the infimum and maximizing over the actual finite cells. The two
previously proved bin and selection errors have the same constants
and the same count m<=(2n^2+1)^6. The absolute-value augmentation uses
log(2m), and the original drift remains unchanged.

The order of limits is explicit: first delta_N tends to zero at fixed
b_0, then any endpoint analysis needs separate justification. No
uniform endpoint estimate or unsupported interchange is claimed.

The final section uses the correctly identified, independently proved
small-canonical-gap implication. Its trace-optimality requirement is
not silently imported into Sections 1--5, where D can be any actual
majorizer. It does not assert that actual optimizers have a small gap,
and it does not evaluate the remaining reference trace supremum.

## Review provenance

The exact worker authored the source. Earlier messages supplied this
reviewer with the proposed statements and constants; the reviewer
mentally checked their factorization and safe constants before this
full-file review, but did not author the metric-stability derivation
or edit its source. The reviewer did author the separately cited
canonical-gap compatibility theorem. This receipt records a complete
analytic check of the final metric-stability source with that prior
discussion and dependency involvement disclosed.
