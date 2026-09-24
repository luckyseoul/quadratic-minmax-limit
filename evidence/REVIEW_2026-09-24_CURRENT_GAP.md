# Review of the current convergence gap and one positive next step

Reviewed baseline: main at
`889a4c88c8aa0634693e9b7db54de6db4fa7c18d`.
The remote main reference was checked against this hash during the review.
Primary sources are the latest HANDOFF.md and its linked proof notes.
This is an analytic author review, not an independent referee report.

## 1. The eight requested commits

| Commit | Verified mathematical scope |
| --- | --- |
| `99e997e1a309514241232111167fab8407133869` | Handoff for the optimized paired-field constant; no additional theorem. |
| `80c6877d247c1bf63d33b33ed11087d7c81e30a3` | Status propagation of that constant. |
| `98225de03ca5bf2ac80d3dca499e15a426de4355` | CORE propagation of the same theorem. |
| `7de631a18a90210d8fb1489b672a218bd12445f3` | Exact optimization of the already-proved fixed-p inequality; the maximizing p is admissible. It proves the exact lower constant whose numerical value is 0.3258474004377944..., not convergence. |
| `c54cc30538be3bd76a33d6078963489cf346977d` | Withdraws the impossible partition-mass specialization and replaces it with an exact Rademacher-tail sufficient condition. It does not prove that a minimizing source satisfies the condition. |
| `f74525f83f0c896a24982239055cf4919460cf25` | Valid partition extension lemma. Its first scalar specialization must be read with the later correction. |
| `dcf7034c639544e2ceacf89c846edb94dc27d68a` | Finite one-sided product bound, interior deficit, and the additional p^2 eta/(24 alpha) mean-update term. All apply without a spectral cap; none gives the post-update parameters. |
| `d9705fced9dbaad096b4b2a09599ae54105bf9af` | Coherent-excess limsup floor and local Paley-gap modulus. Neither is the missing upper comparison. The relevant summable error is the excess over the neutral increment. |

The supporting notes are
[optimized field bound](NOTE_2026-09-23_OPTIMIZED_PAIRED_FIELD_LOWER.md),
[partition and exact tails](NOTE_2026-09-23_PARTITION_EXTENSION.md),
[interior gap](NOTE_2026-09-23_BOOLEAN_INTERIOR_GAP.md),
[coherent excess](NOTE_2026-09-23_PSI_MINIMIZER_EXCESS.md), and
[Paley modulus](NOTE_2026-09-23_PALEY_GAP_MODULUS.md).

## 2. Later results that supersede the quoted handoff

The [tilted paired-field theorem](NOTE_2026-09-23_TILTED_PAIRED_FIELD_LOWER.md)
uses t=993/1000 and raises the baseline to
`B_tilt=0.3258530333538241...`. It reuses the fixed-cap Gaussianization
and removes the cap only after taking n to infinity. A small numerical
optimization alone would not justify that limit passage; the proof does.

The [pointwise disagreement envelope](NOTE_2026-09-23_POINTWISE_DEFICIT_DISAGREEMENT.md)
is valid without a Gaussian law or a cap. For M=A/sqrt(n), with normalized
deficit D=alpha*n-Q_(sM)(x), it gives

    R <= D+sqrt(2alpha*n*D).

Its phase-averaged version recovers the already-optimized mean-update
bound; it is not an extra independent gain at that averaged level.
The pointwise information is stronger than retaining only the means.
Keep units explicit: a raw Q_A deficit O(sqrt(n)) is D=O(1) here and
therefore gives R=O(sqrt(n)). The note's D=O(sqrt(n)) statement uses
normalized energy and corresponds to a raw deficit O(n).

The [stable-state reduction](NOTE_2026-09-23_STABLE_STATE_EXTENSION_REDUCTION.md)
exactly restricts the one-vertex objective to oriented one-flip-stable
states. Integrality is essential: an improving flip gains at least two
in raw energy, while the row correlation can lose at most two.
The [stable-state geometry](NOTE_2026-09-23_STABLE_SKELETON_GEOMETRY.md)
then confines a pair's disagreement to light local-field coordinates.
This is a local counting bound, not a bound on the total number of stable
states or the weighted tail sum.

Thus the new work materially improves a universal lower bound and the
description of the dangerous state family. It has not supplied an estimate
forcing limsup alpha_n to meet liminf alpha_n.

## 3. The missing arguments are still specific

For repeated updates, all-law inequalities constrain permissible
(e,f,eta), but do not describe their dynamics. For example, they allow
e<alpha together with f=2e and eta=0, as for a law on stable states.
They do not force energy to approach the global Boolean norm. A successor
sign law is not Gaussian merely because its input law was a Gaussian sign
law. Reusing the source field estimate for that successor is unsupported.

For a one-vertex comparison, let A_n be an exact minimizer and set

    r_n=m_n*((1+1/n)^(3/2)-1).

It would suffice to find, at each sufficiently large order, an incident
row whose extension increment is at most r_n+rho_n, where rho_n>=0 and

    sum_n rho_n/(n+1)^(3/2) < infinity.

Then the positive increments of alpha_n are summably bounded, and adding
the remaining error tail makes a bounded monotone sequence. This proves
convergence without selecting a limit value.

The stable-state binomial-tail criterion is a sufficient way to obtain
such rows, not a proved property of minimizers and not a necessary
condition for extension. Neither the interior gap nor local stable-state
packing proves its hypothesis. Only raw deficits up to n-r contribute
to that exact-tail test; the global near-edge population in that window
is still uncontrolled.

Likewise, the Paley-gap modulus controls local rises only. A tail bound
uniform in the later prime, or another independent cross-order estimate,
is missing. The coherent-excess LOWER limsup bound cannot be substituted
for an upper bound, and doubling alone is not a convergence theorem.

## 4. Completed positive next step

The smallest additional implication attempted here was to quantify eta
for the INITIAL paired source, so the new interior term is provably used.
The [new proof](NOTE_2026-09-24_PAIRED_DISAGREEMENT_FLOOR.md) establishes

    (1/(2n)) sum_s E k_s >= 1/10-o_L(1),
    eta_n >= 1/1000-o_L(1),

with constants independent of each fixed operator cap L.
It uses the already-proved joint limit with one distinguished input,
the exact input/field covariance, and the paired third/fourth moment
inequality. It needs no new census or law for an updated field.

Let B=B_tilt and keep its existing optimizing p. The interior correction
and the same cap-removal order of limits then give the exact strict gain

    liminf alpha_n >= B_int,
    B_int=[B+sqrt(B^2+p^2/(6000*(1+p^2)))]/2 > B.

This resolves a source-specific part of eta control. The repeated-update
and cross-order conclusions in Section 3 remain unproved. The production
completion registry stays empty, so the original problem stays OPEN.

No test or numerical evaluation was run: AGENTS.md requires offloaded
mathematical checks and the NUKA hostname was not resolvable here. The new
result is an author-reviewed analytic proof with exact scalar comparisons;
no independent verification is claimed. Existing computational receipts
were read as corroboration within their stated finite/algebraic scopes.

Two clerical inconsistencies in the baseline are corrected with this
review: README's pre-optimization lower bound and STATUS's repeated new
number in the sentence describing the old t=1 bound. CORE's stray comma
in the tilted covariance formula is also corrected. No historical proof,
counterexample, or retired route is removed or restarted.
