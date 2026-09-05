# Author receipt: cubic-frame alignment threshold

2026-09-05. Author self-review PASS, with the completed five-comparison
exact-rational result read and matched to its recorded source hash.
This is an author receipt, not an independent-review claim.

## Reviewed artifacts

- Complete final 338-line proof:
  `/tmp/original_mo_scalar_template_cubic_alignment.md`, SHA256
  `60037f67234fbca8c17ee90bf52c7f4346b24e5f18eb5f2c922ebbd2d9382c2a`.
- Complete 74-line five-new-comparison checker:
  `/tmp/original_mo_scalar_template_cubic_rational_certificate.py`, SHA256
  `bc6b5fc08120a2ed645d16ce5a4762919776853bd94d264bf5ce26d6531979ee`.
- Complete result:
  `/tmp/original-mo-cubic-rational.PilWFP/result.json`, SHA256
  `4bc6760b06927a05c104123ff858546ddac1729306c9f7c94a8ad490bc91ad27`.
- Complete run metadata:
  `/tmp/original_mo_cubic_alignment_rational_check.json`, SHA256
  `d865982cd37f28d2138eeccebb76e57ae63a947761d240c32f4094f0f5c67842`.

The source and result hashes were checked in this review. The result
contains exactly five passing Fraction comparisons and the matching
checker hash. The metadata records the root's single bounded remote
run on soulkiller, start and end `2026-09-05T21:31:01Z`, exit zero,
and an absent process at the post-run check. This author did not
execute or repeat the checker and does not claim an independent
live verification of the remote-process metadata.

The Gaussian enclosures, positive D slope, and squared sqrt(2)
endpoint are reused from the previously verified 28-check baseline,
not rerun by the new checker. Its proof, checker, and result hashes
are explicitly retained in the new result and metadata.

## Mathematical scope checked

The source uses actual unit-row optimal frames with `tau(C)=pq`.
The common frame matrix yields
`j_3>=-1+2s^2/mu`; the clipped Gaussian cubic coefficient has square
`2phi^2/3`. These give

    Gamma(C)-1 >= D(q)+4q phi^2 s^2/(3mu)-sqrt(Rs/2).

Combining the actual-frame inequality
`mu<=(Gamma(C)-1)/(q-1)` with the five new rational comparisons
proves `q>=12/5 => Gamma(C)>283/200>sqrt(2)`.
The one-variable penalty maximization is proved analytically in the
source; the checker verifies only its required rational comparison.

The variable-u weak-Dirac diagnostic in Section 6 is derived explicitly
from the general trace expression. It takes n to infinity at each
fixed eta before eta approaches one, retaining the bounded-q repair
trace estimates. Leading energy `f=sqrt(2)` is separately stipulated,
not inferred from the Gamma certificate cap.

Gamma is an upper completion certificate in the constructed-template
application. A lower bound for Gamma does not itself lower-bound an
actual Boolean norm. No actual-rank bound, arbitrary scalar-diagonal
reduction, source compatibility, or original-MO closure is asserted.

## Provenance

The exact worker authored the cubic-frame correction, threshold proof,
new checker, and this receipt. The proof worker's reviewed baseline
completion theorem and exact Gaussian enclosures are retained as
dependencies. The root independently checked the frame correction,
identified and checked the target barrier connection, and requested
the explicit variable-u and order-of-limits extension. The docs worker
independently read the original 293-line proof and the full final
Section 6 extension. The proof and docs workers independently read
the new checker without executing it. Their independent reviews are
separate from this author self-review.
