# Independent exact-worker review: full-SDP gap and compatibility

2026-09-05. Complete analytic PASS; no corrections requested.
The proof worker authored this source. The exact worker independently
read all 303 lines and checked the derivations, including the pure-cross
field corollary. No computation or source edit was performed.

Reviewed source:
`/tmp/original_mo_full_sdp_gap_weighted_compatibility.md`, SHA256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.

## Complete residual and diagonal-spread argument

The canonical primal has unit rows because the actual complete
zero-diagonal signing satisfies diag(K^2)=N-1. Its objective is
tr|K|^3/(N-1), so the stated gap is nonnegative. The same-diagonal
dual whitening gives Qcal Dcal-inverse Qcal<=2Qcal and the exact
4(N-1)g bound on the sum of the two weighted residual squares.

The symmetric orthogonal polar factor exists also on a zero eigenspace
and can be chosen to commute with K and |K|. Thus
`OR_2^T=KD-K|K|` and the commutator identity are exact. The treatment
of its second residual does not move O through D: it bounds the left
inverse square root in operator norm and uses orthogonality to retain
the correctly oriented right weighted residual. The squared triangle
inequality therefore gives precisely 8(N-1)||D-inverse||g.

Because all off-diagonal entries of K have squared magnitude one,
the commutator's squared Frobenius norm is exactly
`2[S tr(D-inverse)-N^2]`. This is the essential use of literal complete
sign entries. The local row-square lower bound
`d_i>=(N-1)^2/S` gives the stated cap-free dispersion bound
`delta<=4Sg/[(N-1)N^2]`.

Cauchy--Schwarz and tr(D-inverse)<=S/(N-1) justify
S>=N sqrt(N-1), eta>=1, and the exact normalized expression
`delta<=4 eta^2 gamma`. No small gap is inferred from a norm cap.

## Uniform original/weighted energy bounds

For t_i=d_i/(S/N), the identity
`E[(t-1)^2/t]=delta` and the pointwise square-root inequality prove
the claimed sum of squared square-root deviations. Factoring
`K=D^(1/2)TD^(1/2)` and using contraction of T gives
`Phi(K-dbar T)<=S sqrt(delta)` uniformly on the actual cube.

The principal-block quadratic estimate correctly uses averaging or
zero-diagonal multilinearity, and the cross-block estimate correctly
uses the two full quadratic half-values after a block sign flip.
Thus the internal energy errors are at most 2N sqrt(delta) and the
cross error at most N sqrt(delta), after dividing by dbar.
The principal-block mixed bilinear estimate has the valid additional
polarization factor, giving 4N sqrt(delta). The two diagonal error
blocks need not be negatives of one another, and the proof does not
assume that they are.

## Actual zero-original-source field corollary

The representative is explicitly selected within the FINAL refined
cell with original p=q_A=0. Hence the just-proved internal error
bounds really apply to its p_D and q_D. The actual field covariance
is positive by the prior weighted-field theorem.

The pure-cross covariance is independently positive: its off-diagonal
block norm is at most k|c_D| ||W_D||<=kn<=wn. Therefore no Gaussian
with an indefinite covariance enters the comparison. Removing the
two diagonal source blocks changes operator norm by at most
2kN sqrt(delta).

There are at most 2^N coefficient states, each of squared norm N.
The finite Gaussian-maximum comparison consequently gives exactly
`2sqrt(k log(2))N^(3/2)delta^(1/4)` for matching fixed offsets.
The safe absolute-value augmentation constant and the conversion to
`2sqrt(2k eta log(2))N^(3/2)gamma^(1/4)` are both correct.
The original cell-selection and bin-noise errors remain separate,
as the source explicitly states.

## Scope and review provenance

The theorem applies to every attained trace-optimal same-diagonal
majorizer, including nonunique optima and singular K. The quantitative
small-gap regime is an additional actual condition, not a statement
that conditional minimizers automatically satisfy it.
The pure-cross field still uses the actual weighted W_D and c_D;
its width is not evaluated, and K/dbar is not declared a contraction.
No original convergence claim is made.

During development this reviewer independently checked the weighted
commutator and uniform energy constants and supplied their square-root
spread proof in discussion. The complete written theorem and the
zero-original-source Gaussian corollary were authored by the proof
worker; this receipt records the subsequent full-source review, not
an assertion of having had no prior mathematical discussion.
