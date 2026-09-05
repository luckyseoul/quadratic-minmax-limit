# Author receipt: fixed-internal-block cross-only regularization

2026-09-05. Author: optimized_profile_docs_gate.
This is an AUTHOR record, not an independent review of my own proof.

## Frozen source and exact delta provenance

Final source: /tmp/original_mo_cross_only_operator_regularization.md.
Length: 264 lines.
Final SHA-256:
27d9ab77768e8b7afa2d48d041cf3fe6bf3b66e8b16e481ca12abcf906a28d4f.

The initial complete source, also 264 lines, had SHA-256:
fcb20c448f022caf374c43ae1fca610e62164862a6d63d75b25c305bf8c34a71.
Both root and optimized_profile_proof independently read that complete
source and returned mathematical PASS. Root also completely read the
two exact prerequisite files.

Root requested one scope clarification in Section 5: explicitly state
epsilon_B>0 before selecting the cross cap. The final source differs
only by changing
"Applying (4) to an actual conditional minimizer gives, for"
to
"For epsilon_B>0, applying (4) to an actual conditional minimizer gives, for".
Reversing that line in a read-only stream reproduces the initial
SHA-256 exactly. The original final blank line is preserved.
The proof reviewer was asked to verify this exact delta independently.

An intermediate version with SHA-256 a4a7b6d8... omitted the final blank
line and is NOT the final artifact. It is superseded by the full final
hash above.

## Author checks and genuine restricted scope

The imported prerequisite is precisely the finite rectangular tensor-
rounding inequality (6) in Section 2 of
NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md. The proof checks its
hypotheses for the rectangular cross matrix B and derives the single
diagonal SDP majorizer of [[0,B],[B^T,0]] with trace at most
Lambda beta(B), Lambda=pi/log(1+sqrt(2)). Bipartite conjugation
provides both majorizer signs using the same diagonal.

The row/column trimming sets have the literal stated cardinality.
Only cross entries are refilled. The rectangular net estimate and
the Boolean maximum estimate have joint positive success probability
already at n=2, so one filler realizes both required bounds.

The essential fixed-internal-block step is the independent whole-set
spin-flip average on the two exceptional sets. It removes the
incident cross entries but also removes the two INTERNAL cuts.
Restoring those unchanged internal blocks incurs exactly the two
printed rectangular cut norms. The theorem does not silently assume
monotonicity of the conditional objective under cross-strip deletion.

The bounded-A corollary, cap (K+8)sqrt(n), and loss
(2+sqrt(2)K_A)sqrt(Lambda C/K)n^(3/2) follow with the stated constants.
A and -A remain literally fixed throughout this construction.

Applying the theorem to a conditional minimum produces a near-
conditional signing, not a new exact conditional optimizer. Its
epsilon_B n^(3/2) slack is retained in any Gaussian-floor use.

The optional source-selection paragraph uses the actual prior
same-order original-norm regularization theorem. It records
K_A=O(epsilon_A^(-2)) and cross K=O(epsilon_A^(-4)epsilon_B^(-2)).
It neither transfers exact optimality to the source selected that
way nor compares its conditional minimum with that of a prescribed
exact original minimizer.

The arbitrary-exact-source limitation and competition with the
evaluated Gaussian cap dependence are explicit. The theorem does
not establish the sharp cross upper or original all-orders convergence.

No repository source or tests were changed, and no numerical
experiment, finite-order census, or unchanged mathematical rerun was
used in writing this proof or author record.

