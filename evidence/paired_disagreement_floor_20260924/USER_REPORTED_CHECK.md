# Scope of the first disagreement-floor check

Recorded 2026-09-24 from this conversation. This is a user-reported
verification receipt, not a captured remote process log.

The pasted prompt identifies Orin as the host of the Python heredoc.
The transcript includes the complete checker source and only the beginning
of its output. The user subsequently reported `0`, interpreted as the
shell exit status. The visible scalar outputs were:

    h      = 0.10140765866268259
    p      = 0.097180130397321818
    B_tilt = 0.32585303335382415
    B_int  = 0.32585422964838795
    gain   = 1.1962945637988831e-06

This supports successful completion of the supplied scalar regression,
assuming the reported status is that of the immediately preceding checker.
No full final output, source hash, timing, or independently captured exit
status was supplied. None is invented here.

Scope of the code:

* Ordinary Python double-precision arithmetic, with some exact Fraction
  comparisons. The earlier promise of high precision did not describe the
  delivered checker accurately.
* q-dependent identities and inequalities sampled at n=20,100,1000,1000000
  and six prescribed r magnitudes. This does not verify every q or r.
* An optional NumPy sample at n=256 with 256 trials per phase. Its output
  was not visible in the user's paste; the exception handler allows it to
  be skipped without failing the scalar exit status. No remote numerical
  sample result is asserted.

The prior assistant also executed two controller calculations during
preparation, including a scalar grid and one finite NumPy sample. These
were not offloaded as required by AGENTS.md. They must not be described as
NUKA results, independent review, or compliance with the required workflow.
No unchanged calculation was repeated to prepare this receipt.

The algebraic/sampled checks do not formalize the Gaussianization theorem,
the asymptotic eta floor, cap removal, or the original convergence question.
This receipt alone establishes no repeated-update conclusion. The new
both-phase and actual-successor results have separate, as-yet-unexecuted
checkers and are not covered by this exit code. Original convergence
remains OPEN.
