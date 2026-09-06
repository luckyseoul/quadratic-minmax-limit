# Author full-read receipt: original-source near-flat strict gain

2026-09-06. This is a CONTRIBUTING-AUTHOR receipt, not an independent
review of a derivation in which its author participated.

## Frozen source read completely

    /tmp/original_mo_original_source_near_flat_strict_gain.md
    612 lines
    SHA256 7726b89e1c39429cde75ff887b981cbd3cf831adb17b04f20193a3c6dbb35298

I wrote the full source, then read every line of the 611-line version,
and read every line again after the sole final change explicitly stated
that the paired corollary's D is diagonal. The final source above passes
my complete analytic reread without an outstanding correction.
The earlier 611-line hash 0759910b88aa9586e2552de64c2dac96873893f110fa4509bbeffbdc85f3c663
is superseded only by that explicit-hypothesis wording change.

## Prerequisite reads and hash verification

The following sources were read completely during this proof chain.
The scalar support and internal-law transfer were also both read in
their entirety immediately before writing the main source; the entire
411-line Gaussian predecessor was reread in this final author turn.
Their current exact hashes were verified alongside the main source.

    original_mo_complete_cross_flat_spectral_gain.md
    411 lines
    b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20

    original_mo_original_source_local_update_scalar_gain.md
    209 lines
    7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155

    original_mo_near_scalar_internal_flat_law_transfer.md
    141 lines
    f65ce2200fd926ba969c9bc5bbaf8ecec8a79b8d228e0f17865fc56c9d9775a8

    original_mo_original_phase_spectral_moment.md
    262 lines
    7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600

    original_mo_source_cross_nuclear_trace_boundary.md
    106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556

The last source's already verified pi enclosure is reused, not rerun.
No new mathematical checker or numerical verification was executed.

## Core links checked

1. The bounded empirical flat law yields both Frobenius projector
   approximations. Completeness fixes diag(A^2)=n-1, and diag(A)=0
   fixes the signed-projector diagonal error. Normalizing only good
   coordinates and independently filling the o(n) bad coordinates
   gives an actual full correlation R with bounded operator norm and
   ||R-R0||_F=o(sqrt(n)). No exact finite-n flatness is assumed.

2. Uniform Schur-power operator bounds and the actual entries imply
   |tr(A R^{circ q})|=O(n) for every odd q>=3. The original positive
   phase has mean 5kappa/8+o(1) after normalization by n^(3/2).

3. The identity A^2/lambda^2-rho R+A/lambda=o_F(sqrt(n)) pairs with
   every higher odd Schur power uniformly. Its even entry powers
   contribute at least the diagonal n, giving mean higher-chaos
   variance >=1-kappa-o(1), not merely half that quantity. Individual
   variances lie in [0,2(1-kappa)/rho+o(1)].

4. First-chaos alignment is established by the actual trace of
   (A-lambda I)^2 R. The difference from R0 is controlled using the
   Frobenius norm of the square and requires no inverse of R.

5. The Gaussian lemma allows zero or smaller row coefficients. Its
   distinguished-coordinate contraction is O(1/n) in squared norm.
   The finite-chaos characteristic equation, uniform L2 tail, compact
   covariance subsequences, and the Gaussian sign-disagreement
   continuity set justify BOTH local-field absolute moments and the
   sign-mismatch probabilities, including degenerate Gaussian limits.

6. The independent Bernoulli update is Boolean on the original source.
   Zero diagonal removes the diagonal second-order term exactly.
   Its penalty is 2epsilon^2 C_n p_n with the ACTUAL cap
   C_n<=5/3+o(1), not the limiting atom 5/4. The fixed epsilon=1/10
   is admissible without optimizing over unconstrained probabilities.

7. The heterogeneous-variance chord and probability bounds yield
   liminf improvement >=16/3125>1/200. The original quadratic norm
   therefore has lower bound >=5kappa/8+16/3125>2/5+3/1100.
   There is no bilinear-norm substitution or lost factor of two.

8. The separately authored transfer lemma applies to one common
   original principal source and has precisely the spectral and
   operator premises needed above. Its composition excludes only the
   stated actual near-scalar internal-law regime at objective 2/5.

## Contribution and remaining scope

Root supplied the positive-projector/local-update strategy, corrected
the relevant operator cap, proposed the live fixed probability, and
emphasized the trace-of-square variance argument. I supplied the robust
correlation construction, the stronger mean, the distinguished-coordinate
Gaussianization extension, and the full proof write-up. The docs-gate
worker supplied the scalar support; the exact worker supplied the
separate transfer. The earlier 411-line argument itself has disclosed
root/exact/proof contributions. These are not independent-author claims.

No canonical repository file was edited by me for this result. Only
/tmp proof and receipt files were written, and read-only metadata
checks were used. No mathematical computation was run.
The all-profile/all-active-cell implication and the original MO target
remain OPEN. Publication, repository gates, commits, and backups are
root's separate workflow, not assertions made by this receipt.
