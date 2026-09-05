# Root review: actual nuclear coupling and formal trace boundary

2026-09-05. Analytic/source review PASS. No mathematical execution.

## Exact source and review coverage

Source: `/tmp/original_mo_source_cross_nuclear_trace_boundary.md`.
Final complete length: 444 lines.
SHA-256:
`106cc8ae8bb4e2d7f4024f18ffc8114e123299a276005b7ce31ebab3ab74e556`.

Root read the complete initial 370-line note, then the complete final
444-line note, not only the added threshold section. All formulas and
scope statements passed. No checker, solver, census or optimization ran.

Root independently checked the finite source/cross derivation and
centered-law formal example. Root AUTHORED the shifted-threshold extension
in final Section 9. Therefore this receipt is not an independent review
of that extension. The exact worker independently checked and integrated
it; the proof worker independently reviewed the complete mathematical
argument. Separate final-source receipts document their coverage.

## Finite actual theorem

For any positive feasible same-diagonal D, the literal off-diagonal
sign squares and Jensen applied to the reciprocal square roots give

    ||K/dbar-T||_F^2 <= (N^2/dbar^2)(2delta+delta^2).

The internal nuclear-norm loss divided by n is consequently at most
sqrt[(2N/dbar^2)(2delta+delta^2)]. The n, N=2n, and dbar factors
in the nuclear prerequisite were checked directly. The cross identity
m=(n/dbar^2)ell h and ell h<=(1+delta)^2 give the correct LOWER
bound on sqrt(n/dbar^2). Operator-monotone square roots transfer
A_L^2<=I-WW^T to the actual trace of sqrt(I-WW^T), without a
commutation assumption. The same holds for A_R.

Feasibility implies dbar^2>=N-1; the original nuclear bound gives
alpha>=kappa sqrt(1-1/n)/2. Thus the claimed absolute
O(sqrt(delta)+1/n) error for 0<=delta<=1 needs no trace cap.
No trace optimality, small canonical gap, trimming of the actual field,
or assertion that K/dbar is a contraction is used.

## Formal data and retained constraints

At alpha=2/5, f=4/3, u=kappa and m=9kappa^2/16, the full and
cross moments were derived from the displayed commuting block model.
The projection construction preserves the stated second moments and
zero internal diagonals in rational approximations. It does NOT preserve
complete-entry magnitudes or identify alpha, f and u as actual Boolean
norms and an active cross value.

Root checked both full-normalization terms, the actual nuclear coupling
at zero error, and the internal nuclear and cubic leading bounds. The
balanced flat nonzero internal spectrum makes its zero-odd-diagonal
cubic bound equal kappa/[2sqrt(1-m)], exactly the nuclear condition.
The actual cross cubic necessary bound is also retained and saturated.
No template Gamma conclusion is silently imposed on an actual norm.

## Exact certificate failure and the threshold extension

The centered-law evaluation uses the SAME endpoint measure in both
traces. The Euclidean triangle inequality, affine-ratio endpoints and
every displayed rational comparison yield the uniform squared lower
40501/125000, exceeding the original target square 8/25 by
501/125000. Signed metrics and continuous endpoints are included.

For each finite shifted threshold, the planar Gaussian square is
contained in the stated disk. This proves k_h/(1-s_h^2)<=kappa.
The endpoint function C(s) is decreasing for 0<=s<=kappa^2,
so the entire earlier lower bound survives; the individual A_*>7/10
bound is correctly NOT reused at smaller s. Negative metrics cannot
lower it. The original unweighted pure-cross drift is |s_h| f/2
after division by 2n^(3/2), and is explicitly added to the noise.
The elementary inequality sqrt(1-z^2)>=1-z gives the lower bound
min{f/2,sqrt(40501/125000)}>sqrt(2)alpha for all z in [0,1].

The six mathematical prerequisites cited by the source were read fully
earlier in this research session. Root reread the complete original-phase
and cross-moment prerequisites and the relevant exact field/threshold
formulas during this review; their hashes match the source. The earlier
pi enclosure and certificate were reused without another execution.

## Status

PASS: finite actual necessary coupling and a precisely specified formal
drift-plus-ellipsoid certificate obstruction. Failure of an UPPER
certificate is not a lower bound on an actual Gaussian width, an actual
signing counterexample, or impossibility of other spectral methods.
Entrywise realizability, active Boolean states and their additional frame
constraints remain outside this relaxation. Original all-orders
convergence and its possible limit value remain OPEN.
