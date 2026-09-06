# Author checkpoint: near-scalar actual-cross transfer

2026-09-05. Author: optimized-profile proof worker.

Final source: `/tmp/original_mo_near_scalar_cross_spectral_gain.md`,
complete 364 lines, SHA256
`ec911854e59788fabbb4e189d47849acedff15a1c80dbd9225a373a49e62d1f9`.

I reread the complete final written source and checked its mathematics
end to end. Author verdict: PASS for the actual full-cross-moment
transfer and its explicitly stated active-state consequence.
This author checkpoint is not an independent review or publication.

The initial 364-line draft had SHA256
`3a81d89fcf4345085a842640f122ccf9b51165aecfe9faf814b32d31135efc55`.
The sole change was to make the interlacing test function unambiguous:
Section 4 now says `(max(x,0))^(2k)`. No formula, theorem hypothesis,
constant, consequence, or source dependency changed. Root and both
reviewers were informed before their final receipts were requested.

## Exact new transfer

For an actual complete paired signing, feasible positive diagonal D,
S<=C(2n)^(3/2), and 0<=delta<=1/512, write W for the full actual
weighted cross block, m for its second singular moment and v_2 for
its fourth singular moment. The source proves

    beta(B)/(n dbar)>=kappa v_2/m+(sqrt(kappa)-kappa)m
           -[25kappa C^2+6(sqrt(kappa)-kappa)]delta^(1/3)-R_C(n),

where R_C(n) tends to zero uniformly. No finite-n rate is asserted.
An actual original pure-cross active state gives the same lower bound
for its actual u_D after the explicit further loss 2sqrt(delta).

The balanced complete block has q/n=a>=1/2, scalar scale
d'=(1+eta)dbar, and m'=q/d'^2>=1/(9C^2). It provides only the
exact original norm lower beta(B)>=beta(B_J). Its first and fourth
singular moments are compared back to FULL W using dilation
interlacing and contraction congruence, with error
2theta+4k a eta for k=1,2.

The positive-part finite lower retains both normalization factors
a and 1+eta. The one-sided comparison m'<=m_0<=m then avoids
amplifying a subtracted error by the larger coefficient. The literal
cross-entry identity m_0<=m<=m_0(1+delta)^2 controls the retained
gain term. These yield exactly the printed dispersion constant.
The error envelope over all q>=n/2 supplies one uniform R_C(n).
The delta=0 case is handled directly, without dividing by zero.

I also checked the compatibility argument afresh: Q=diag(sqrt(t_i))
has ||Qz||=||z||=sqrt(N) on full Boolean vectors, and its rescaling
error gives the cross loss 2sqrt(delta) without trace optimality.
The pure-cross active substitution separately uses c=beta(B)=Phi(K).

## Dependencies, contributions, and non-claims

The complete 411-line actual-cross theorem is the substantive imported
inequality. Its fixed hash is
`b30903b22c0b602464a864b78b59be6827bb0c110e6cc382c753f3ea0a16fb20`.
I reread the full 280-line near-scalar normalization and full 303-line
weighted compatibility sources at their hashes stated in Section 9;
their relevant arguments are rederived in this transfer.

Root proposed balanced good-coordinate trimming and returning the
first/fourth moments to full W. I derived the finite coefficients and
errors and authored the note. The exact worker independently checked
the candidate formulas; its optional stronger inverse-mean refinement
is not needed or used. Their prior prerequisite contributions are
to be disclosed in their own review receipts. The docs-gate worker's
full-source review is separate and independent of this derivation.

No mathematical computation, checker, solver, signing construction,
spectral scan, or numerical run was used. All artifacts are in /tmp;
no canonical repository source was edited by the proof worker.

The theorem removes the extra unweighted operator-cap hypothesis
inside the near-scalar fixed-trace-cap branch. It does not prove that
every optimizer has delta tending to zero, that all active cells are
pure-cross, or that every remaining cross measure permits the desired
ellipsoid upper. The original global inequality remains open.
