# Long-horizon goal (binding)

**Done means exactly one thing:** MathOverflow [413935](https://mathoverflow.net/questions/413935) is settled.

That is: prove that

\[
L=\lim_{n\to\infty}\alpha_n
\]

exists (whether or not the proof identifies its value), **or** prove that the
limit does not exist.  The original MathOverflow question asks for existence;
requiring the value as well is a strictly stronger, optional objective.

A session, plan, prop, test suite, or handoff update is **not** completion of this goal.

## The only acceptable terminal states

| Terminal state | What must be true |
|----------------|-------------------|
| **Existence CLOSED; value unidentified** | A proved convergence theorem applies to the actual sequence.  For example, both Dini-summable amplification rays in Proposition 6.3 are proved, so \(H(n)=m_n^{2/3}\) satisfies the two-ray criterion.  A reduction with an open hypothesis is not enough. |
| **\(L=1/2\) CLOSED** | E(1) on the Paley family \(n=p^2+1\) is proved Max+-free for all primes \(p\ge5\), denseness (Prop 6.1–6.2) is applied, and `e1_closed_general` is True via **real imports** from the hinge modules (not a handwritten `return True`). `solution.md` Main Theorem states \(L=1/2\). |
| **\(L=c\) CLOSED for a specific \(c\neq 1/2\)** | A proved value with the same wiring standard. |
| **Non-existence CLOSED** | A rigorous proof that \(\alpha_n\) does not converge.  Two ratio-dense subsequences with unequal proved limits are one sufficient mechanism, not an extra requirement on every proof. |

Anything else — including “honest OPEN,” “structure shipped,” “census holds at \(p=5,7\),” “AI-test later,” “good increment” — is **not done**.

## What is already proved (do not re-derive)

- Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le 1/2\).
- Denseness: the limit exists globally iff it exists along Paley orders \(n_k=p_k^2+1\).
- Two-ray convergence reduction (Prop 6.3): Dini-summable normalized defects
  at multipliers 2 and 3 suffice. In particular,
  errors `O(n/(log n)^(1+epsilon))` in `H`, or power-saving errors, force
  convergence. Propositions 6.4--6.5 identify the exact four-state
  equal-endpoint diamond for the all-Hadamard multiplier-2 construction;
  the hereditary endpoint conditions are automatic. Proposition 6.5a
  identifies that diamond exactly with the maximum norm of the directed
  half-cut flips of `A` under a chosen tournament orientation. This is a
  simultaneous upper-neighbor problem; global optimality alone gives only
  the reverse lower bound. Proposition 6.5b shows that the `sqrt(2)`
  multiplier is forced on every dyadic ray: a fixed smaller multiplier with
  vanishing normalized error would contract the sequence to zero. Random
  vertex orders currently stop at an exact variance gate requiring new
  global-minimizer control, and fixed-real Hermitian interlacing stops at one
  spectral edge. Proposition 6.5c gives a distinct opposite-diagonal diamond
  with arbitrary cross block. Its symmetric specialization is an exact
  hybrid-slice/QPSK problem; the coherent holomorphic choice already violates
  the zero-loss bound at a global order-four optimizer. Its general
  four-label form isolates the all-directed target: one tournament must
  `1/sqrt(2)`-pave every signed cut submatrix simultaneously. Proposition
  6.5d expresses the same skew gate as a nonlinear cover by decomposable
  Boolean bivectors and proves that its affine relaxation, displayed
  elliptope, and normalized single-row even-moment certificates are
  subcritical in their stated ranges. Richer SOS/Pluecker couplings remain
  open. Proposition 6.5e proves a signed-regular
  arcsine lower bound meeting the outgoing-half target at exactly the
  `1/sqrt(2)` constant, plus conditional approximate-commuting rigidity near
  the universal `1/pi` floor. Proposition 6.5f shields anchor-incident
  constraints for every fixed or low-signature family, conditional on an
  open vertex-cover condition. Proposition 6.5g constructs the required approximate
  square/commutator mate near that floor but proves its generic spectral
  conversion loses `pi/2`. Proposition 6.5h rules out the exact independent
  first-moment certificate, while Proposition 6.5i proves Gaussian
  saturation forces the Hamming-central two-half saddle near `alpha=1/pi`.
  The latter is necessary structure, not an orientation. Only a noncoherent, `A`-dependent cross
  construction or a direct high-degree cover theorem remains live there.
  Proposition 6.6 proves
  the diamond outside the explicit Hamming-central/joint-energy residue
  (6.20), but does not close the ray. Central Paley conference maximizers
  prove that those coarse residue data alone cannot close it. Proposition 6.7 gives an exact
  equal-endpoint tetrahedral frame for multiplier 3, with only a \(3n\)
  internal-edge error. Its single-skew form is an exact three-state diamond;
  two spectral shields are proved, but the unshielded complement remains
  open. Proposition 6.8 gives an independent `1:2` composition: its
  bi-balanced Hadamard cross block proves the exact two-state diamond when
  `k_A k_B<=n^2/100`, leaving only (6.42)--(6.43). Fixed finite anchor
  refinement changes only the `O(n)` border, but alternating paired states
  prove that this refinement cannot empty the residual using the current
  spectral bound. Proposition 6.9 proves that the
  former signed-Eulerian fallback fails for every fixed temperature `c>0`,
  including `c=3`; only a growing-temperature formulation remains logically
  possible.
- \(\rho=1\) on that Paley family, so \(\Phi(C_n)=\frac12 n\sqrt{n-1}\) and \(\alpha_{n_k}\to 1/2\) **if** E(1) holds there.
- 15.167's majorization algebra is conditionally valid, but its final bi-tight implication is **retracted** because `ker(G-(n/2)P1)` contains `ker G`. Proposition 15.720 instead excludes the required bi-tight levels 2 and 3 for every prime `p>=5` by a degree congruence. Its bi-tight level-4 corollary does not exclude one-sided tight level 4.
- Propositions 15.726--15.727 historically narrowed the first general
  residual-(ii) endpoint. Proposition 15.733 closes the former
  `p=31,R=10` case, and Proposition 15.734 supersedes that endpoint route by
  closing the whole `k=4p` shell for every boundary and every `p>=13`.
  Propositions 15.735--15.737 extend the first-three-shell close to every
  `p>=11`; Proposition 15.751 closes the fourth shell `k=4p+6` for every
  `p>=13`; Proposition 15.752 closes the fifth shell for every `p>=23`
  plus its stated contiguous higher band; and Proposition 15.753 closes its
  exceptional p17/p19 endpoints by exact common-energy certificates.
  Proposition 15.754 closes the last `p=13,k=60,u=6` residue by an exact
  finite aggregate/common-form certificate. Thus the fifth shell is closed
  for every prime `p>=13`. Propositions 15.768--15.770 close the first two
  generic post-band layers in each congruence class. They also close the
  exceptional `p=23,t=9,k=110` and `p=23,t=10,k=112` endpoints; the latter
  has eleven low roots and reuses the fixed 33,649-five-set quartic/octic
  certificate. Proposition 15.771 additionally closes `p=23,t=11,k=114`
  for every boundary size. The mean-46 equality classification and common
  row ledger force mass-32 opposite cells, which pointwise parity
  subtraction and the local `p+9` theorem exclude. See the
  [endpoint proof](evidence/NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md)
  and [four-node replay](evidence/p23_third_post_band_mesh_replay.json).
  Proposition 15.772 closes the third generic p1 post-band layer
  `p=1 mod 4,p>=29,t=q-1,k=5p-3`. Its punctured complement-triple
  theorem repairs the gap-two premise of 15.770 and retains the new
  gap-four equality of offset four; common-row normalization then forces
  forbidden opposite masses, including the newly excluded `p+11`.
  See the [generic proof](evidence/NOTE_2026-09-04_P1_THIRD_POST_BAND_CLOSE.md)
  and [local gap theorem](evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md).
  Proposition 15.773 then closes `t=q,k=5p-1` for every prime `p>=29`.
  All old low-row branches carry; the new flat mean-`2p` branch directly
  forces five forbidden mass-`p+9` rows, without a new equality
  classification. See the [joint proof](evidence/NOTE_2026-09-04_JOINT_5P_MINUS_ONE_CLOSE.md).
  Proposition 15.774 adds `t=q+1,q+2,k=5p+1,5p+3` for every
  prime `p>=29`. On `J(p,(p+1)/2)`, its strict spectrum
  `0<4p E[C]<2p-10` for nonnegative integral quadratics permits
  only Boolean masses `p-3,p+1`; the affine-parity union is
  `{0,p-3,p-1,p+1}`. Coupled type quotas exclude signed shell floors
  `r=3,4,5` through `5p+4,6p+4,7p+6` at `p>=37`, with uniform
  bounds `5p-12,6p-12,7p-12` at every `p>=29`. Residual carries
  cover p29/p31 as well. See the
  [local proof](evidence/NOTE_2026-09-04_SHARP_SMALL_MASS_SPECTRUM.md) and
  [two-type bridge](evidence/NOTE_2026-09-04_SMALL_MASS_TWO_TYPE_BRIDGE.md).
  The generic frontier is `t>=q+3,k>=5p+5`.
  493 technical tests passed on soulkiller; the final documentation-gate result is recorded in `evidence/small_mass_two_type_regression.json`.
  Provenance is recorded in `evidence/small_mass_two_type_regression.json`.
  Exact scalar survivors at the first uncovered eventual layer are
  not graph witnesses. No all-size witness localization is proved.
  Eventual E1 for all sufficiently large primes would suffice for
  `L=1/2` by denseness, but bounded-size exclusions are not eventual
  E1 and give no `o(p^3)` deficit estimate. All global gates remain open.
  Propositions 15.755--15.756 prove a sharp
  full-cube defect dichotomy and rule out arbitrary-boundary Weil/Parseval as
  a global close; both are reductions/method barriers, not residual closure.
- Exact Paley optimality is **false** (\(m_{10}=13<15\)). E(1) is asymptotic.

There are now two honest routes.  The direct route to the original question
is the two-ray amplification theorem of Proposition 6.3 (or any other proof
of convergence/non-existence).  The value-specific route to \(L=1/2\) is
**E(1)** on \(n=p^2+1\).  On that optional Paley route, the live gate is the
corrected dictionary returned by
`src/e1_main_chain_status.py`; as of 2026-09-04 through Proposition 15.774 it
reads:

1. **Required bi-tight levels 2 and 3:** TRUE by 15.720. The spectral floor,
   global mixed-\(k\) QVAR, and principal R1 are no longer acceptance gates.
2. **Residual (ii), even \(k\ge4p\):** OPEN. Propositions 15.734--15.737
   close `k in {4p,4p+2,4p+4}` for every boundary and every prime `p>=11`.
   Propositions 15.738--15.742 close `p=13,k=58`, and Proposition 15.743
   closes `p=17,k=74`. At `p=13,k=60`, the complete residue sieve leaves
   `u in {0,3,4,6}`; Propositions 15.744--15.745 close `u=3,0`, and
   Propositions 15.746--15.749 close `u=4`. Proposition 15.754 closes the
   remaining `u=6` by exact joint common-form and collision-energy
   certificates, completing `p=13,k=60`. Propositions 15.768--15.770 close
   the next two generic post-band layers and the exceptional
   `p=23,t=9,10` endpoints. Proposition 15.771 closes `p=23,t=11,k=114`.
   Proposition 15.772 closes `p=1 mod 4,p>=29,t=q-1,k=5p-3`.
   Proposition 15.773 closes `p>=29,t=q,k=5p-1` in both congruence classes.
   Proposition 15.774 closes `p>=29,t=q+1,q+2,k=5p+1,5p+3`.
   With `q=(p-1)/2`, the exact live frontier is
   critical `p=5,7`; `p=11,t>=3` (`k>=50`); `p=13,17,19,t>=5`;
   `p=23,t>=12` (`k>=116`); `p=1 mod 4,p>=29,t>=q+3`;
   `p=3 mod 4,p>=31,t>=q+3`; and the positive `p=7,z=7` subbranch.
   Proposition 15.751 closes generic branch B at `p>=29,t=3`, hence the
   entire fourth shell for `p>=13`; Proposition 15.752 closes the displayed
   higher band, and together with Propositions 15.753--15.754 closes the
   fifth shell for every `p>=13`.
   Historical
   endpoint/profile artifacts are not live gates; only the lemmas and
   certificates explicitly retained by the dedup audit remain valid evidence.
   Read `evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md` before reopening one.
3. **Type I, multi-level Max−:** TRUE by Proposition 15.750. Isolated-chart
   square rigidity and parity halving close `p>=11`; tracked exact integer
   Farkas identities close `p=5,7`. The old `3A+B` route remains incomplete
   but is no longer a global gate.
4. **Lemma D:** TRUE.

The older statement that residual (ii) was closed by 15.179/15.236/15.237
covered an earlier, narrower split and is not the live multi-level predicate.
Likewise, “residual (i) is the only leftover” is retired shorthand. The
authoritative open function is `residual_ii_k_ge_4p_ND_closed()`;
`type_I_multilevel_bad_case_ND_closed()` is now True. The valid bi-tight predicate is
`required_bitight_levels_empty_all_primes()`.

Path-C / 16N / 15.193 exhaustiveness is **not** required.

Closing the sole open residual-(ii) unit completes the remaining E1 ledger,
but it is not necessary for a different proof of convergence and is not goal
completion. The direct convergence front is to close the multiplier-two
diamond on the exact residue (6.20), then close Proposition 6.7's unshielded
tetrahedral diamond or close Proposition 6.8's two-state residual
(6.42)--(6.43). Endpoint
selection is retired on both explicit frames; independent skew-norm budgets
and an all-pairs composition theorem are not the live target. The disk
surrogate is strictly stronger, not equivalent: its asymptotic form would
prove a new `1/sqrt(2pi)` lower bound, and only its zero-error form is
currently disproved.

## Forbidden translations (these are the small-failure mode)

Do **not** replace this goal with any of:

- Ship a Max+-free identity / dead-path / layer rewrite whose `residual_i` / `e1` / \(L\) flags stay False, then treat the session as done.
- “Continue residual-(i)” as a plan whose acceptance criteria are “a new increment exists.”
- Recertifying that unsigned \(\sum|\mathrm{per}|\) exceeds \(B\) (15.231 already did this).
- Updating STATUS / HANDOFF / session notes as the deliverable.
- Census \(p\le 7\), SA, or Lipschitz scale-counting as E(1).
- Treating the all-prime, gap-2 Paley architecture as the only way to answer
  the original existence question.
- Soft-close: sandwich + denseness + \(\rho=1\) \(\Rightarrow L=1/2\) without E(1).
- Flipping predicates without a general Max+-free hinge imported for real.

If a proof attempt fails, **stop**. Record the failed mechanism in one paragraph. Do not wrap the failure as Prop 15.xxx with `proved: True` on an identity that does not flip a leftover.

## Session rule

A session may end only if one of these happened:

- one of the two fixed amplification rays in Proposition 6.3 is actually
  proved with a Dini-summable error, and the other ray remains the named direct
  convergence gate, or
- the remaining open E(1) predicates, residual (ii) and the minimal-four-gap
  implication bridge, **actually closed** (True via real imports), or
- a live route was **killed as a path** by a general counter-mechanism (not a
  small-\(p\) census), and the writeup names the replacement route, or
- A terminal state in the table was reached, or
- The human redirected the goal.

“We shipped structure and left \(L\) OPEN on purpose” is **not** an allowed end condition for this goal.

## After a terminal state

The writeup must stand on its own.  If closure uses the Paley route, an
independent cold read of `evidence/share/denseness_path_package.md` must call
that argument essentially correct; a different convergence proof needs its
own self-contained package. Channel is X + GitHub, not a MathOverflow answer.
