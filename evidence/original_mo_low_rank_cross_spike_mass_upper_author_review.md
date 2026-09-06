# Author check: low-rank cross-spike mass upper

2026-09-06. AUTHOR / CONTRIBUTING check, not an independent whole-proof
review. The exact worker authored the new source and supplied new links.

Frozen source fully self-read, all 347 lines, without truncation:

    /tmp/original_mo_low_rank_cross_spike_mass_upper.md
    SHA256 30347140ecf9fb2458444fb152490c601fe81d8a1733e90f31be692126ecdf1c

Result: PASS of the displayed conditional theorem and corollary.
No mathematical correction was identified during the full self-read.
Before freezing, two clarifying lines were added to say that truncated
spectral traces retain the ORIGINAL dimension as denominator. The whole
347-line final source was then reread. Independent frozen-source review
remains separate and is required before publication.

Analytic checks covered:

- Actual weighted spike mass t and cross allocation v obey |v|<=t
  and |u-v|<=r0(1-t). The PSD metric uses e Pi on the spike space;
  its radius retains ORIGINAL u=2c/S exactly.
- Projection loss is bounded without a small-state-mass assumption by
  2n sqrt(w[d+N sqrt(delta(1+delta))]). It is o(n^(3/2)).
- The scalar translated-Gaussian bound follows by rotating independent
  copies and centering the even decreasing positive-part integrand.
- Weighted Fenchel/Cauchy--Schwarz retains the mass constraint inside
  the Boolean remainder, giving e sqrt(kappa S R)sqrt(1-t).
  Correlated and unequal projected Gaussian marginals are permitted.
- The transformed field couples uniformly to an op-controlled field
  at cost N sqrt(wn)sqrt(delta)(1+sqrt(1+delta)). Two such costs
  and a low-rank sphere net justify adaptive choice of direction.
- Finite mass bins require no separate v bins. Their selection is
  controlled using the ACTUAL g, and the mesh/net limits follow n.
- Fixed r0 above the limiting bulk edge gives rank(Pi)=o(n). The
  cutoff decreases only after n; the metric inverse gap remains 1/10.
- The pure-cross trace reduction imports the fully read 279-line
  theorem. Its actual operator-radius convergence premise is not used
  as a surrogate for weak-law convergence.
- The final envelope F0 sqrt(1+t/rho)+G0 sqrt(1-t) is concave, with
  F0/rho<7/25<G0. Thus every allocation is bounded by its t=0 value,
  which the frozen prerequisite bounds by 14/25<2sqrt(2)/5.

Fully read prerequisites and rechecked frozen hashes:

    381-line weighted-shell source:
    9aec82a5e808837ea626f2fd85f526cda1fffe883929711dfc2c6f396392f15f
    279-line actual-radius source:
    44fa3e7361e2142b20dce58d2dde727458db786529690f15e752390b8081725f
    312-line pure-cross source:
    035c8e9d042fe8b54773784988356d16ed7c1257f35c470c5c64aa68dd65cfa6

Root supplied the translated-Gaussian inequality and mass-inside-
remainder idea. The exact worker supplied allocation, weighted finite
bounds, conditioning and source authorship. The proof worker checked
those links and independently checked the derivative conclusion before
writing. These roles do not qualify either contributor as an independent
whole-new-source reviewer.

No mathematical program or test ran on any host. Only the new /tmp
source and this receipt were written; earlier sources remain frozen.
No canonical file, publication, or backup was changed or performed.
The result does not close the original MO convergence question.
