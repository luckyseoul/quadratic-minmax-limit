# Independent review: strengthened formal profile and all-threshold boundary

2026-09-06. Reviewer: optimized-profile proof worker.

Verdict: PASS for the precisely listed strengthened trace relaxation
and its stated drift-plus-ellipsoid reference certificate. No correction
is requested. This is not an actual-signing counterexample or a lower
bound on actual Gaussian width.

## 1. Complete final source and independence

I read all 227 lines of the initial source, then reread all 227 lines
of the final source after its two publication-status wording changes:
`/tmp/original_mo_strengthened_trace_profile_all_threshold_boundary.md`,
final SHA256
`903ac72c78c60706fbcfef09e50abeda0a18fe05976e3efab89d65becdbfccf1`.
The earlier SHA256
`488ba3d4142e0b2d5f49e1eafad6e5e61a1fb9494c3d97e6531f33414c8ff846`
is superseded. No formula changed between those versions.

The root supplied this profile and the exact worker derived its
supporting-line argument. I had no role in selecting these parameters,
deriving that supporting line, or authoring this diagnostic. This is
an independent review of the new profile and its full analytic argument.

I authored the earlier general pure-cross reference evaluation and
actual complete-cross gain/near-scalar transfer prerequisites. That
history is disclosed; this receipt does not claim independence from
the development of every prior theorem. I previously read the entire
444-line boundary, 312-line reference, 280-line normalization, and
364-line gain-transfer sources, and refreshed their exact hashes for
this review. I also directly rechecked Section 2 of the 312-line
source, especially the general-u formula (2.2); I did not mistakenly
reuse a formula specialized to u=kappa.

The refreshed prerequisite hashes match the final source:

    444-line original boundary
    106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556
    312-line general reference
    035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6
    280-line full normalization
    c679c9155845aa2b51c55e72b781a72f7122f27cb4b2d7c8be69fec178172fd2
    364-line near-scalar actual gain
    ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9

## 2. Formal model and retained necessary inequalities

I checked u=f sqrt(m)=4/5, r=2m=18/25, and
a=sqrt(m/(1-m))=3/4. The projection/complement construction has
HP=PH=0 and H^2=a^2(I-P), so its block contraction and the three
listed spectral laws are algebraically consistent. Its internal and
cross second moments both equal m. Nothing asserts sign-sized entries
or assigns actual Boolean norms to these real matrices.

The full moments are exactly mu_1=21/25, mu_2=18/25,
mu_3=63/100. The full cubic normalization becomes 301kappa/400
and the full nuclear normalization becomes 6kappa/7. Both lie below
u=4/5 under the reused coarse kappa enclosure. The moment-gap label
is 1-mu_3/r=1/8, with the needed formal-only qualification.

For the internal law, the source nuclear and balanced common-zero-
odd-diagonal cubic conditions both require alpha>=5kappa/8. One
can check the latter independently: the unweighted source absolute
law has nonzero value (3/4)/sqrt(m)=5/4 with mass 16/25;
its normalized cubic moment is 5/4, so its phase lower bound is
5kappa/8. Thus the source condition passes strictly at alpha=2/5.

The source/cross nuclear coupling has left side 16/25 and right
side 3kappa/4<12/25. The old cross cubic floor is kappa. The
new actual-entry condition is the convex combination
(16/25)kappa+(9/25)sqrt(kappa)<436/625<4/5. All the listed
constraints therefore pass, including the new gain rather than only
the earlier cubic floor.

These checks are not a realization theorem for actual source/cross
entries or the original norms alpha and f. They do not supply an
actual pure-cross state attaining u. The source explicitly retains
that distinction, which is essential to the diagnostic conclusion.

## 3. Same-law functional, thresholds, and the rational lower bound

The independent-Gaussian square-to-disk probability comparison gives
w>=exp(-h^2), hence s=(k/w)u<=kappa u<64/125. Factoring
sqrt(w) leaves the original completion constant sqrt(kappa), not
sqrt(k). Substituting the endpoint law into both expectations of
the general reference gives exactly (3.2). The first factor is
sqrt(1-u)=sqrt(1/5), as appropriate for this changed profile.

The cancellation giving T_2 in (3.3) is valid for t<1 and provides
its continuous t=1 limit. The same law is used in both terms;
neither separate extremal measures nor separate Jensen envelopes
have been substituted for those expectations.

The first term satisfies T_1>=(19/50)t because
73/500=365/2500>361/2500. In the second term,
t/(1+t)^2<=1/4 and s<64/125 give the stated ratio lower bound
61/250. Multiplying by the lower kappa enclosure yields exactly
A=252/625 and B=34587/625000.

I checked the supporting-line comparison by exact algebra:

    (19/50)^2/A=361/1008,
    B(1-361/1008)=355203/10000000,
    (47/250)^2=353440/10000000.

Thus the proposed two-coordinate vector has norm strictly below
one. Its Cauchy--Schwarz comparison is strict because B>0, and
gives sqrt(Ax^2+B)>(19/50)x+47/250 for every x in [0,1].
Taking x=1-t cancels the t dependence against the first-term lower
bound. Therefore U_s(t)>71/125 uniformly over the whole permitted
s interval and every t in [0,1], including both endpoints.

The target is the original sqrt(2)alpha=2sqrt(2)/5. Its squared
comparison has numerator difference 5041-5000=41 over 15625.
The rational supporting-line computation and that positive margin
are analytic identities, not outputs of a metric scan or new checker.

## 4. Negative metrics and original drift

For the negative metric, the zero-atom first coefficient increases.
At the unit atom the first numerator difference is 2(u-s)>0;
the second unit-atom numerator also increases. Thus the negative-
metric reference is no smaller than its positive counterpart at the
same magnitude. The cancellation already checked handles its limit.

The actual convention for the ORIGINAL formal pure-cross drift is
zf/2 after division by 2n^(3/2). Adding it to the noise reference
and using sqrt(1-z^2)>=1-z gives a lower bound at least 71/125,
since f/2=2/3 is larger. Infinite thresholds have zero noise and
drift 2/3, so they cannot evade the conclusion. The target has not
been changed to f/2 at any stage.

## 5. Exact stopping point

The enlarged LISTED trace constraints remain insufficient to make
this reference certificate at most the desired target, even after
all shifted Gaussian sign thresholds and all signed metrics in this
one-cross-shell family are allowed. Reoptimizing those parameters
under only these constraints cannot close this route.

This does not prove insufficiency of every conceivable moment
constraint or spectral method. It also does not lower-bound the
actual Gaussian width: a lower bound on the value of an upper
certificate has a different meaning. Actual entry/Boolean-state/
conditional-optimality/frame information may still exclude the
profile or permit a sharper width bound outside this certificate.
The new actual-cross theorem and its near-scalar transfer remain
valid and are not retracted. The original global target remains open.

No mathematical computation, new enclosure, checker, solver, metric
evaluation, signing search, or parameter search was run. This
review used source reads, hash checks, and analytic reasoning only.
No canonical repository file or the reviewed source was edited.
