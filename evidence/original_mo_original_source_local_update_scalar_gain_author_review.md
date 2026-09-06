# Author receipt: explicit-cap scalar local-update gain

2026-09-06. Author checkpoint, not an independent proof review.

Final source: `/tmp/original_mo_original_source_local_update_scalar_gain.md`,
complete 209 lines, SHA256
`7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`.

I directly reread the complete final written source after its revision
and checked every coefficient, rational comparison, probability bound,
and scope statement. The hash has been reverified for this receipt.
Author verdict: PASS for the scalar implications expressly stated there.

## Live fixed-update statement

The scalar setting is rho=16/25, kappa=2/pi,
a^2=kappa/rho, 0<=v<=V=2(1-kappa)/rho and mu=E v.
The source retains the proposed G and p formulas and the explicit
operator-penalty coefficient C in

    J_C(epsilon)=epsilon G-2C epsilon^2 p.

For mu>=1-kappa and C<=5/3, the admissible FIXED epsilon=1/10
gives J_C>16/3125>1/200. The baseline 5kappa/8 is greater
than 35/88, so baseline plus this gain exceeds 2/5+3/1100.
Under the more permissive C<=3, fixed epsilon=1/20 gives
J_C>9/3125>1/400, with baseline margin 1/4400 above 2/5.

Both thresholds remain eventually valid under an o(1) loss in the
strong mean premise and an o(1) enlargement of the respective cap.
The strict gaps also absorb o(1) errors in a separately justified
objective lower. No homogeneity or convergence of individual variances
is assumed, and no finite-order convergence rate is asserted.

## Scalar mechanisms and exact checks

The square-root chord gives G>=c mu with
c=(4/5)sqrt(kappa)/(sqrt(kappa)+sqrt(2-kappa))>8/25.
Arctan(x)<=x and Cauchy--Schwarz give
p<=sqrt(mu)/(pi a)<(8/25)sqrt(mu).
These bounds hold for every probability law on the specified interval.

At C<=5/3, the fixed-update lower is
(4/375)(3r^2-r), r=sqrt(mu)>3/5. This is increasing on that
range, and its boundary value is 16/3125. At C<=3, the analogous
polynomial is (10r^2-3r)/625, with boundary value 9/3125.
The rational comparisons 3200>3125 and 3600>3125 verify the
two gain thresholds. The baseline deficits use
2/5-35/88=1/440, followed by
1/200-1/440=3/1100 and 1/400-1/440=1/4400.

The interval 7/11<kappa<16/25 is obtained from the existing
verified pi enclosure, not a new computation. In particular
7*31415927=219911489<220000000 proves its upper endpoint
is below 22/7. All subsequent rational comparisons were checked
analytically; no old or new certificate was executed.

## Optional calculations and contribution limits

The earlier 162-line version with SHA256
`d2f36d20dfedc254114892b5bc1767dd383c5b08f21b6ee035d041a3cee2992c`
is superseded. Its half-mean calculation at C=5/4 remains only an
exact-scale specialization in Section 6. The limiting nonzero spectral
value 5/4 must not replace the live operator cap 5/3.

The optional quantity G^2/(8Cp) is explicitly the UNCONSTRAINED
quadratic maximum. Implementing its optimizer would require the extra
check G/(4Cp)<=1. Neither fixed-probability corollary needs that check.

Root supplied the scalar setting, corrected the relevant operator cap,
and relayed the proof worker's stronger mean premise. I derived the
distribution-free chord/probability bounds and exact-scale margins.
Root proposed the cap-5/3 extension and fixed 10% update; I checked
its constants and stability and supplied the cap-3, fixed-5% fallback.
Those are contributions, not an independent review of the combined lemma.

This scalar note does not prove the actual Gaussian local-field limit,
the actual variance constraints, or the claimed original quadratic-energy
effect of independent updates. Those premises require their separate
actual-source proof and complete review. No original-norm theorem is
being inferred here merely from a variance lower bound.

No mathematical run, checker, spectral scan, optimization, signing
search, or canonical repository edit was performed by me. The work
consisted of analytic reasoning, scoped /tmp writes, complete source
reads, and provenance checks. This receipt does not publish the result.
