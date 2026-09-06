# Author check: actual two-cross-moment nuclear transfer

2026-09-06. This is an AUTHOR / CONTRIBUTING check, not an independent
whole-new-proof review.

Frozen source fully self-read, all 255 lines, without truncation:

    /tmp/original_mo_near_scalar_two_cross_moment_source_nuclear_transfer.md
    SHA256 32c1e47608c1dc06037ababb6c3c34934fdea5546df17527e94f6509562e6525

Result: PASS of the displayed conditional statements; no correction
was identified during the full self-read. Independent full-source
review remains a separate requirement before publication.

Checks completed analytically:

- The finite Schatten factorization has Frobenius factors M-z, 1-v,
  z, v. Trace positivity gives z<=d; concavity gives z_*=min(d,vM).
  No commutation of the full actual H with WW^T is assumed.
- Exact complete-cross sign squares and reciprocal half-means force
  n/dbar^2->m from positive first cross moment and delta->0 alone.
- One common label set has q/n->1. Congruence and interlacing give
  |tr H^2/n-(q-1)/dbar^2|<=2b+9epsilon, hence tr H^2/n->m<=1/2.
- Full nuclear control transfers by unnormalized trace-norm compression.
  Unbiased Boolean extension preserves the ORIGINAL source comparison;
  no paired covariance, cross matrix, active state, or optimizer changes.
- S(m,m,Delta)/sqrt(m) is exactly the stated continuous piecewise L.
  The angle inversion gives Delta_crit=m(4sqrt(m)-3sqrt(1-m))^2/25.
  Its exactness concerns the certified cap rectangle, not necessity
  for actual source exclusion. The interior rectangle uses 1917<1922.
- partial_w Psi<=1/2 gives F_C'(u)<0 for C>=1, 0<u<=1.
  Nuclear moments of complete capped sources stay bounded away from 0;
  subsequences, fixed C_m+eta, and continuity justify the composition.
- Compactness yields the no-moment-convergence version when every
  actual moment accumulation point lies in the stated closed region.

The 230-line endpoint transfer and 553-line all-law gain were fully
read in this work sequence and their frozen hashes were rechecked:

    230: 6a486df0fd46aa76259e3f02e3734eb2529162500f98f89af58e90562e6a2187
    553: 0a7c553e29d4e3ac1572edb0e3fc795bc4d252d090061181365f01764c500a51

Root contributed the finite Schatten bound. The exact worker checked
that bound, supplied threshold inversion, monotone composition and
source transfer, and wrote the full source. Root also checked the new
threshold link. Docs reported pre-writing checks without contributing
a new derivation; a full frozen-source docs review is still separate.

No mathematical programs or tests ran on any host. Only the new /tmp
source and this receipt were written. Earlier sources remain frozen;
no canonical file, commit, push, or backup was changed or performed.
This is conditional original-source progress, not global MO closure.
