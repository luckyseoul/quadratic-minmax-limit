# Result: elementary knife-edge certificate for the order-16 doubling straddle

## Prior art (not claimed as new)

The research archive
[luckyseoul/quadratic-minmax-limit](https://github.com/luckyseoul/quadratic-minmax-limit)
already identifies the **m₁₆ straddle**: with \(F(8)=10\),

\[
2\sqrt{2}\,F(8)=2^{3/2}\cdot10\approx28.284,
\]

and \(F(16)\in\{28,30\}\) (now tightened from an older \([27,30]\) bracket)
sits on either side of that threshold. See e.g.

- `evidence/NOTE_2026-09-11_GROWTH_CHAIN_M16.md`
- `evidence/NOTE_2026-09-12_COMPLETION_DISCREPANCY_FRAMEWORK.md`
(“first open doubling instance”)

So the **qualitative** statement “\(F(16)=28\) meets the doubling target;
\(F(16)=30\) misses it” is **not new**.

## What this note adds

An **elementary integer certificate** that the CORE §7 / \(H=F^{2/3}\) form
of the same doubling inequality is equivalent to that dichotomy, with a
one-line proof that needs no floating point:

## Theorem

Assume \(F(8)=10\) and \(F(16)\in\{28,30\}\). Let \(H(n)=F(n)^{2/3}\). Then

\[
H(16)\le 2\,H(8)
\quad\Longleftrightarrow\quad
F(16)=28.
\]

### Proof

\(H(16)\le 2H(8)\) is exactly \(F(16)\le 2^{3/2}F(8)\). Cubing the comparison
\(H(28)<2H(8)<H(30)\) gives the integer inequalities

\[
28^2 < 8\cdot10^2 < 30^2
\qquad\text{i.e.}\qquad
784 < 800 < 900,
\]

which are immediate. Hence \(H(28)<2H(8)<H(30)\), and the dichotomy follows. □

### Defect if \(F(16)=30\)

\[
\eta_8=\frac{H(30)-2H(8)}{16}\approx0.0232322611769.
\]

## Machine checks in this repo

| Check | Result | Artifact |
| --- | --- | --- |
| Integer knife-edge | `784<800<900` verified | `results/experiment7_knife_edge.json` |
| Framed Hadamard | independent `F(16)≤30` PASS | upstream verifier re-run |
| Scalar \(B_2\) locator | \(T\approx0.32586697\), \(B_2>B_{\mathrm{opt}}\) | `results/experiment5_B2_ladder.json` |
| Search for \(\Phi\le28\) | no witness (best \(30\)–\(32\)) | `results/experiment5_*` |

## Status

The original limit remains open; \(F(16)\) remains open. This repo’s
contribution relative to upstream is the cubed-integer packaging plus the
reproducible experiment harness, not a new asymptotic theorem.
