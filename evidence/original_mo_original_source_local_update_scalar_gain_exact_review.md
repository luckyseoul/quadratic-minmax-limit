# Independent review: scalar original-source local-update gain

2026-09-06. Reviewer: optimized_profile_exact.

## Complete reads and independence

I directly read the complete 209-line source
`/tmp/original_mo_original_source_local_update_scalar_gain.md`, SHA256
`7de99c4bbf997fc25eafa2742cb55c220dc13fdf29d0b1ae535358ea8c73f155`.
I also fully read its 86-line author receipt, SHA256
`74693d17bc646355a5fdac4be9aace460c8f85fd3f983bf22d45491fd2939fb1`,
at `/tmp/original_mo_original_source_local_update_scalar_gain_author_receipt.md`.

I refreshed and fully read the complete 444-line named interval source
`/tmp/original_mo_source_cross_nuclear_trace_boundary.md`, SHA256
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.
All three hashes were checked directly. Its already verified pi enclosure
is reused; no old certificate was rerun. Only that enclosure, not the
source's separate coupling or relaxation claims, is used by this lemma.

I supplied no derivation, parameter choice, correction, or proof step to
the new scalar lemma. Root and the docs-gate worker's contributions are
expressly recorded in its source and receipt. I authored the separate
near-scalar internal transfer, which is NOT a premise of this scalar
lemma. I also contributed to the older 444-line interval source, so this
receipt claims independence for the NEW scalar derivation, not for every
older result in its provenance chain.

## Live implication

The definitions are rho=16/25, kappa=2/pi,
a^2=kappa/rho, V=2(1-kappa)/rho, and an arbitrary probability law
on 0<=v<=V with mean mu>0. The proposed gain and mismatch formulas are

    G=sqrt(kappa) E[sqrt(a^2+v)-a],
    p=(1/pi) E arctan(sqrt(v)/a),
    J_C(epsilon)=epsilon G-2C epsilon^2 p.

The scalar lemma does not establish these formulas for an actual signing.
It correctly keeps the operator cap C separate from the limiting nonzero
spectral value 5/4, and uses an admissible fixed update probability.

The old pi enclosure gives 7/11<kappa<16/25; the displayed integer
comparison 219911489<220000000 correctly proves the needed pi<22/7.
The endpoint chord for the concave square-root difference yields

    G>=c mu,
    c=(4/5)sqrt(kappa)/(sqrt(kappa)+sqrt(2-kappa))>8/25.

The strict inequality is equivalent to 13kappa>8, which follows from
kappa>7/11. All square roots and squared sides are positive.
Arctan(x)<=x and Cauchy--Schwarz give
p<=sqrt(mu)/(pi a)<(8/25)sqrt(mu), using
pi a=5/(2sqrt(kappa))>25/8.

For C<=5/3 and epsilon=1/10, these bounds give

    J_C(1/10)>(4/375)(3r^2-r),       r=sqrt(mu).

The strong mean premise mu>=1-kappa implies r>3/5. The displayed
polynomial is increasing there, so the gain exceeds 16/3125>1/200.
The baseline 5kappa/8 exceeds 35/88, and
2/5-35/88=1/440, 1/200-1/440=3/1100. Thus the claimed live margin
above 2/5 is correct. Neither homogeneity of v nor an optimized update
probability was assumed.

## Remaining displayed claims

I checked all optional calculations, not just the live conclusion.
For C<=3 and epsilon=1/20 the lower is
(10r^2-3r)/625, increasing for r>=3/5, whose endpoint value is
9/3125>1/400. Its baseline improvement is 1/4400.

Since mu>0 implies p>0, the unconstrained quadratic maximum is
I_C=G^2/(8Cp), with maximizer G/(4Cp). The same scalar bounds give
I_C>mu^(3/2)/(25C). The strong-mean lower constants 81/15625 at
C=5/3 and 9/3125 at C=3, and their stated comparisons, are correct.
The source correctly does NOT infer implementability without separately
checking that the unconstrained maximizer is at most one.

The optional half-mean calculation at C=5/4 gives
I_C>27sqrt(2)/15625>1/420, hence margin 1/9240 above 2/5 after
adding the baseline. The comparisons using sqrt(2)>7/5 and
79380>78125 are valid. This is explicitly an exact-scale specialization,
not a replacement for the live operator cap or strong mean hypothesis.

With mu>=1-kappa-o(1) and C<=5/3+o(1), bounded p and fixed epsilon
make the extra penalty o(1). The strict gap 16/3125>1/200 absorbs it.
The analogous C<=3 stability is valid. No variance convergence,
homogeneity, finite-order rate, or unconstrained optimizer is required.

## Verdict and scope

PASS for the entire frozen scalar source, with no required correction.
This is an independent review of its conditional scalar implications.
Actual projector normalization, the local Gaussian approximation and
variance constraints, joint sign mismatch, and the original quadratic
energy update still require their separate complete proof and review.
This receipt alone asserts no original-source theorem or profile exclusion.

No mathematical computation, checker, numerical integral, optimization,
signing construction, or search was run on any host. Tools were used only
for complete reads, provenance checks, and writing this receipt in /tmp.
No canonical repository file was edited and no publication was performed.
