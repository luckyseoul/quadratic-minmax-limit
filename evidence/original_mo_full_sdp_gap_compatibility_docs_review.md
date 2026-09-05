# Independent full review: canonical SDP gap and weighted compatibility

2026-09-05. Reviewer: optimized_profile_docs_gate.
Verdict: PASS; no corrections requested.

## Exact source and independent scope

The reviewer directly read all 303 lines of
`/tmp/original_mo_full_sdp_gap_weighted_compatibility.md`, SHA-256
`3a1367bab1fe73aa24c0edbdb1bb583546e28ae82148f4cf5af749e49b9778f0`.
The reviewer did not contribute to this theorem or its derivation and
made no source edits. This new note and receipt remain separate from
the preceding four-note publication package.

The complete covariance and weighted-field prerequisites had also
received the reviewer's full independent reads. Their final hashes
match the exact references in this source:

- Covariance, final 384 lines:
  `0b3921d43d88424457ad2ed777ee158e8ac34c6751c995f0c3b86aee870e95ff`.
- Weighted field, final 381 lines:
  `9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f`.

The same-D existence proof uses symmetry and the complete-signing
structure, so its SDP symmetrization and strict positivity argument
also apply at the present general order N>=3. It does not require
the special paired form until Sections 5--6.

## Actual canonical primal and inverse-weighted residual

The rows of K/sqrt(q) and |K|/sqrt(q) are genuinely unit vectors,
because |K|-squared=K-squared and every diagonal entry is q=N-1.
Their objective is tr(|K|-cubed)/q, proving nonnegativity of the
stated FULL bipartite SDP gap, not a Boolean or cross-only gap.

The doubled dual covariance Qcal is between zero and 2Dcal. Its
quadratic trace on the actual canonical frame Z is exactly 2g,
with the correct unit-row and factor-two normalizations. Whitening
and squaring therefore gives Qcal Dcal-inverse Qcal<=2Qcal.
Substituting the two residual blocks proves precisely (2.1), with
right side 4qg. No maximum-diagonal substitution enters this step.

## Commutator and the exact diagonal-spread statistic

An orthogonal symmetric polar factor can be completed by signs on
the zero eigenspace. Its stated products with K and |K| remain exact
there. In particular O R_2-transpose=KD-K|K|, which gives (3.1)
with the required transpose and matrix order.

For the first normalized residual the RIGHT inverse-square-root
factor is bounded in operator norm. For the second, the LEFT factor
is bounded, followed by orthogonality of O and transpose invariance
of Frobenius norm. This proves both displayed inequalities without
commuting O through D. The squared triangle inequality then gives
the exact 8q b g bound in (3.2).

Completeness makes every off-diagonal K_ij-squared equal to one.
Expanding the normalized commutator norm therefore gives exactly
2[S tr(D-inverse)-N-squared], including the vanishing diagonal terms.
Combining it with (3.2) proves the first bound in (1.1).

The row-square contraction bound yields d_i>=q-squared/S and
tr(D-inverse)<=S/q. Thus b<=S/q-squared gives the second bound
in (1.1). Cauchy--Schwarz establishes delta>=0 and, with the inverse
trace bound, S>=N sqrt(q). The normalization eta>=1 and the identity
delta<=4eta-squared gamma in (1.3) have the correct powers and factors.
Since the canonical primal objective is positive, gamma lies in [0,1).

## Uniform cube rescaling and the three actual energy bounds

With t_i=d_i/dbar, the three averaging identities in Section 4 are
exact. The scalar inequality comparing the squared square-root
deviation with (t-1)-squared/t holds for every positive t. It gives
V<=S delta without an upper bound on any individual diagonal entry.

For every real cube vector, both rescaled vectors have norm at most
sqrt(S), and their difference has norm at most sqrt(V). Contraction
of T gives (4.1), hence Phi(K-dbar T)<=S sqrt(delta) with the correct
one-half convention. This is a pointwise uniform cube bound, not an
SDP-frame average or a bound restricted to a rounding law.

The two principal-block estimates follow by averaging the other
Boolean block; its quadratic term has mean zero because its diagonal
is zero. The cross estimate follows by flipping just one block and
subtracting the two quadratic energies. It has constant Phi(E),
whereas each principal quadratic energy has constant 2Phi(E).
Dividing by dbar=S/N gives exactly all three constants in (1.4).
The argument does not assume opposite diagonal blocks for E.

Ordinary polarization on the real cube gives the mixed-source bounds
4N sqrt(delta) in (5.1). The cross estimate already allows arbitrary
Boolean vectors from the two blocks. Substitution of (1.3) gives the
stated 4eta N sqrt(gamma) source and 2eta N sqrt(gamma) cross bounds.

## Positive pure-cross field stability

The representative is expressly chosen inside the FINAL refined cell
with original p=q_A=0. This preserves the original-zero premise when
applying the uniform bounds and is compatible with the earlier binning
argument. No representative from a larger unrelated weighted bin is used.

The actual M_theta is PSD by the weighted-field theorem. Its pure-cross
counterpart M_0 is independently PSD: the off-diagonal block has norm
at most k|c_D| ||W_D||<=kn<=wn. Thus neither Gaussian requires an
indefinite covariance or a scalar-D hypothesis.

The two removed diagonal blocks have norm at most
k max(|p_D|,|q_D|)<=2kN sqrt(delta), proving (6.3).
For at most 2^N coefficient states of squared norm N, the Gaussian
operator-error maximum comparison gives exactly
2sqrt(k log(2)) N^(3/2) delta^(1/4). Deterministic offsets are valid
and need not be constant. Augmenting both signs and using N+1<=2N
gives the stated safe absolute-value constant.

Taking the fourth root of delta<=4eta-squared gamma proves exactly
(6.5). Zero k or a zero gap creates no singular case in these bounds.
The earlier within-cell approximation and cell-selection errors remain
explicitly present when the corollary is used on binned cells.

## Conditional conclusion and remaining work

The gap is that of the literal COMPLETE matrix's vector SDP. Every
trace-optimal same-D solution is covered, including nonunique solutions;
singular K is covered without inverting K or |K|.

A fixed original-norm cap bounds eta but is not claimed to force
gamma to vanish. The theorem therefore identifies a quantified actual
SMALL-GAP regime, not an established property of original or conditional
minimizers. Its negligible internal-field conclusion still leaves the
actual weighted W_D and c_D in the pure-cross field and does not evaluate
that field's width. It gives no conference-scale operator cap, no license
to use an indefinite unweighted tensor covariance, and no original MO
convergence claim. All these limits are stated correctly in the source.

No mathematical computation, certificate execution, signing construction,
census, simulation, numerical SDP, search, solver, optimization, or new
test was run by this reviewer.
