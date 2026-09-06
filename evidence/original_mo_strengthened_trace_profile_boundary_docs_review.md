# Independent review: strengthened formal profile and all-threshold boundary

2026-09-06. Verdict: PASS. No mathematical correction requested.

## 1. Final source, complete reads, and contribution boundary

I directly read the complete 227-line source initially supplied at
SHA256 `488ba3d4142e0b2d5f49e1eafad6e5e61a1fb9494c3d97e6531f33414c8ff846`.
After the two documentary publication-wording changes, I directly
reread all 227 lines of the final source:

    /tmp/original_mo_strengthened_trace_profile_all_threshold_boundary.md
    903ac72c78c60706fbcfef09e50abeda0a18fe05976e3efab89d65becdbfccf1

The equations, profile, and claims did not change in that final delta.
This receipt applies to the final hash above.

I had no role in selecting or deriving this new parameter profile,
checking it before the source was written, or deriving its supporting
line. Those new arguments are independently reviewed here. I previously
contributed the directional AM-GM refinement in Section 2.4 of the older
444-line actual-coupling source. That prerequisite contribution is
disclosed; it is not authorship of this new profile or supporting-line
argument, and I do not claim independence from my own earlier step.

I directly reread these complete prerequisites for this review:

    original_mo_source_cross_nuclear_trace_boundary.md, 444 lines
    106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556

    original_mo_small_gap_pure_cross_upper.md, 312 lines
    035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6

    original_mo_original_phase_spectral_moment.md, 262 lines
    7108222bd693fd65b11e552a7a4138654dd96d032bda24dc5c61d7abc92dc600

    original_mo_cross_singular_moment_rounding.md, 168 lines
    6d5129a1572842c76c8f11a008b0093cb3c340684a40219b7db8828fdeeaf756

The other two direct prerequisites were already completely read in
this working interval and their hashes were reverified:

    original_mo_near_scalar_diagonal_spectral_normalization.md, 280 lines
    c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2

    original_mo_near_scalar_cross_spectral_gain.md, 364 lines
    ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9

The latter has my separate full independent 190-line review, SHA256
`c74ca051421336eabd5380bbfc4537b24f495d3a7459a8ba0638c03f0055b6b0`.
The earlier coarse kappa enclosure is reused as established in the
444-line source. No baseline checker or new rational checker was run.

## 2. Formal model and the exact retained constraints

The profile has alpha=2/5, f=4/3, u=4/5, m=9/25, r=18/25,
and a=3/4. It satisfies `u=f sqrt(m)`, `r=2m`, and
`a^2=m/(1-m)` exactly. With P a projection and H annihilating its
range, the stated commuting block matrix has square
`diag(H^2+P,H^2+P)<=I`. Its internal second moment is
`(1-m)a^2=m`, while its full and cross laws are exactly those printed.
The prior weighted-trace or limiting rank interpretation remains formal;
it supplies no complete-signing or Boolean realization.

The full absolute moments are `21/25`, `18/25`, and `63/100`.
Substituting them into the two retained full-normalization bounds gives
exactly `301kappa/400` and `6kappa/7`, each strictly below 4/5.
The formal scalar moment-gap label is `1-(63/100)/(18/25)=1/8`.
It is not an actual realized SDP-gap assertion.

For the internal law, the normalized nonzero absolute eigenvalue of
the putative A/sqrt(n) is `a/sqrt(m)=5/4`, on mass 16/25.
Thus its nuclear moment is 4/5 and its cubic moment is 5/4.
The retained nuclear and common-zero-odd-diagonal cubic conditions
both require precisely `alpha>=5kappa/8`. The latter is the
zero-odd-diagonal formal specialization described in the older block
construction, not a generic inference from spectral symmetry alone.
Since kappa<16/25, the proposed alpha=2/5 passes strictly.

The source/cross coupling has left side 16/25 and right side
`3kappa/4<12/25`. The earlier cross cubic condition reduces to
`u>=kappa`, since the endpoint law has moment `integral y^(3/2)=m`.
For the new full-weighted cross-entry condition, v_2=m and zero
dispersion give `u>=kappa+(sqrt(kappa)-kappa)m`.
Its right side is `(16/25)kappa+(9/25)sqrt(kappa)`.
Using kappa<16/25 and sqrt(kappa)<4/5 gives exactly the upper
436/625, which is strictly below 4/5. Thus the new condition is
actually checked, not omitted after it eliminated the earlier profile.

These are precisely checks of the listed necessary inequalities.
They do not assign actual Boolean norms alpha and f to a real matrix,
or produce an original active state with the proposed u.

## 3. Same-law functional and global supporting line

For finite h, z is the standard Gaussian probability of the symmetric
interval [-|h|,|h|]. The product probability z^2 is bounded by the
Gaussian disk probability `1-exp(-h^2)`, giving
`w>=exp(-h^2)` and `0<=s=(k/w)u<=kappa u<64/125`.
This bound covers both signs of h. The scalar factor sqrt(w) is
separated without changing the completion coefficient sqrt(kappa).

Evaluating both terms of the general reference formula at the SAME
endpoint law gives (3.2). It is not a combination of two unrelated
extremizing laws. Absorbing 1-t into the second square root gives
exactly (3.3), so its extension to t=1 is continuous and finite.

The first-term bound uses `1-u=1/5`, `1-m=16/25`,
`1+s>=1`, and `(1+t)^2<=4`. Its squared coefficient is
73/500=365/2500, strictly above `(19/50)^2=361/2500`.
The displayed weak inequality also holds at t=0.

The identity for the second-term ratio is exact. Since
`t/(1+t)^2<=1/4`, it is at least `(1-s)/2>=61/250`.
Multiplying by kappa>63/100 yields the constants

    A=(63/100)(16/25)=252/625,
    B=(63/100)(9/25)(61/250)=34587/625000.

I checked the supporting-line arithmetic analytically. In particular
`(19/50)^2/A=361/1008`, so the remaining factor is 647/1008.
The product with B is `549*647/10000000=355203/10000000`.
Meanwhile `(47/250)^2=353440/10000000`, proving the strict
inequality in (4.4). The comparison vector therefore has norm
strictly less than one. Cauchy--Schwarz gives a STRICT lower supporting
line for every x in [0,1], since `(sqrt(A)x,sqrt(B))` never vanishes.

Adding that line with x=1-t to the first-term bound cancels t and
gives `U_s(t)>71/125`, uniformly for the entire stated rectangle of
t and s, including the continuous t=1 endpoint. The final exact
square comparison is `5041/15625>5000/15625`, with squared
margin 41/15625. This is a global analytic inequality, not a scan,
numerical minimum, or finite set of sampled thresholds or metrics.

## 4. Negative metrics, full original drift, and all endpoints

For the reflected metric, the first zero-atom coefficient increases.
The unit-atom squared numerator increases by `2(u-s)>0`, since
s<=kappa u<u. The second unit-atom numerator also increases,
by 4st. All relevant denominators are unchanged. Thus the negative
metric cannot reduce this same-law reference functional, and its
continuous endpoint is covered by the same canceled expression.

The original pure-cross formal drift is `|s_h|c`; after the stated
normalization it is exactly zf/2, not a substituted weighted drift.
The inequality `sqrt(1-z^2)>=1-z` gives the full certificate lower

    z(2/3)+(1-z)(71/125)>=71/125>sqrt(2)alpha.

The target remains `sqrt(2)alpha=2sqrt(2)/5`. The quantity f/2
is only the drift scale and is not substituted for the target. At
infinite positive or negative threshold, the noise vanishes and the
drift is 2/3, still above target. Thus both metric signs, both metric
endpoints, all finite thresholds, and both threshold infinities are
covered. No uncontrolled actual-field error is taken to a metric
endpoint; the result is about the exact formal reference functional.

## 5. Exact scope and disposition

The changed formal profile passes the stated strengthened trace
relaxation while its drift-plus-ellipsoid reference certificate fails
the original target uniformly. This explains why the new actual-entry
gain alone does not close that relaxation. It does not invalidate the
gain or its full-weighted near-scalar transfer, and does not retract
the earlier, separately scoped relaxation statement.

A lower bound on this UPPER certificate is not a lower bound on
actual Gaussian width or an original Boolean norm. The profile is
not an actual-signing or conditional-optimizer counterexample. Further
source-entry, frame, active-state, or optimality information can still
exclude it; those conditions are not supplied by the listed moment
tests. The stopping point concerns unchanged trace-only use of these
constraints, not all spectral methods or all convergence proofs.

No mathematical computation, checker, parameter search, metric scan,
test, canonical edit, or publication was performed by me. This receipt
is the only new file I wrote for this review. Final-source analytic
and scope verdict at SHA256 903ac72c...: PASS, with no requested changes.
