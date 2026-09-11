# Exact Paley structural two-block lifts

For a Paley conference signing `C` of order `n=q+1`, this artifact fixes the
complete cross block `B=C+I` and exactly evaluates

\[
K=\begin{pmatrix}C&B\\B^T&-C\end{pmatrix}.
\]

This is a defined algebraic family, distinct from the seeded local-repair
search.  The off-diagonal cross block is complete because the zero diagonal
of `C` becomes `+1` in `B`.

| q | source order | source Phi(C) | lifted order | exact Phi(K) | verification |
| ---: | ---: | ---: | ---: | ---: | --- |
| 5 | 6 | 5 | 12 | 18 | direct local enumeration |
| 13 | 14 | 21 | 28 | 70 | NUKA HIP and 16-thread CPU replay |
| 17 | 18 | 33 | 36 | 108 | NUKA HIP and separate Soulkiller 88-thread CPU replay |

The order-28 entry gives the rigorous finite construction bound `m_28<=70`.
The order-36 entry similarly gives `m_36<=108`.  The raw q=13 and q=17
cross blocks are retained here. Their SHA-256 values are respectively
`87d02b925ee8ae1634332ac76a637f309b080321472ca649082d620b969205d0`
and `82f5517bc9bf5a012207958d1ef418d15941fff97032673e943714878139728b`.

The generic CPU and HIP scorers have SHA-256
`b1b6c95f582e86f807b181aa9b56599f8a56cb331f1f6be8edc28a86659d8b3c`
and `15e81fdeabfd4bab1de1b2620b242845ff1cd1aa80cf792dacc90af346dcc763`.

## Cone-repair neighborhood at q=13

The multi-active Navier cone driver was started from the exact `q=13`
structural block, rather than a random block. Its 16 ranked coordinates were
drawn from the exact tied-maximizer band. No single or pair flip lowered 70.
The expanded one-stage test then exhaustively scored all 696 ranked one-,
two-, and three-entry proposals, again with no strict improvement. This is a
finite local-basin observation only; it neither proves a cross-block optimum
nor excludes a larger or differently ranked repair. The retained result is
`cone_triple_result.json` (SHA-256
`4dd38233e2561d4470ccc0bc98119bf5f617046320a62c858b83fc696e4dc004`),
and the driver SHA-256 is
`94046b35ea4009e4889a1d79363b601845838f9ea791a61e9c17e59b34ed3b66`.

## Interpretation boundary

Three finite orders do not prove a formula for this family, establish an
asymptotic upper bound, or solve the original convergence question.  In
particular, the normalized values of these examples are not asserted to be
monotone.  Their value is to provide exact, independently replayed data for
a specific scalable algebraic cross-block mechanism.
