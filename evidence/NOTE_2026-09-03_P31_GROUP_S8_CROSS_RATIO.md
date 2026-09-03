# The `p=31`, weight-eight grouped branch: exact fourth-direction reduction

## Scope

This note freezes a symmetry reduction and an exact-certificate plan.  It does
**not** report the large enumeration, exclude the weight-eight case, prove the
group-support lemma, determine the halved row-code minimum distance, close
residual (ii), or solve the original MathOverflow problem.

Let `S` be a support of eight antipodal point classes.  A counterexample to

\[
 |S|+\#\{\text{active direction groups}\}\ge 32
\]

would have at least nine silent projective directions.

## Four generic cross-ratio cases suffice

Move any three silent directions by `PGL(2,31)` to the unordered triple
\(\{0,1,\infty\}\).  Its setwise stabilizer is `S_3`; on a fourth direction
\(\lambda\) it gives the six anharmonic transforms

\[
 \lambda,\quad 1-\lambda,\quad \lambda^{-1},\quad
 (1-\lambda)^{-1},\quad \frac{\lambda}{\lambda-1},\quad
 \frac{\lambda-1}{\lambda}.
\]

Their orbits on \(\mathbb F_{31}\setminus\{0,1\}\) are exactly

\[
\begin{aligned}
 &(2,16,30),\\
 &(3,11,15,17,21,29),\\
 &(4,8,10,22,24,28),\\
 &(5,7,9,23,25,27),\\
 &(6,26),\\
 &(12,13,14,18,19,20).
\end{aligned}
\]

The size-three orbit is harmonic and the size-two orbit is equianharmonic.
Their union contains only five directions.  After fixing the silent triple, a
putative weight-eight counterexample still has at least six silent directions.
Those six cannot all lie in that five-point exceptional union.  Consequently
one further silent direction lies in a generic orbit.  It is therefore enough
to treat the four representatives

\[
 \boxed{\lambda=3,4,5,12}.
\]

This is a complete reduction from 29 possible fourth directions to four exact
cases.  It does not assume that a prescribed one of the six remaining silent
directions is generic.

## Exact `4+4` certificate plan

For one representative \(\lambda\), attach to every one of the 480 antipodal
point classes the 60-bit incidence signature consisting of the fifteen
nonzero squared-projection blocks in each of
\(0,1,\infty,\lambda\).  Enumerate all

\[
 \binom{480}{4}=2,184,297,480
\]

four-subsets, sort their signatures, and inspect equal-signature groups for
disjoint pairs.

If two disjoint four-subsets have equal signatures, their eight-point union
has even parity in all sixty recorded blocks.  Conversely, if an eight-point
support is silent in the four fixed directions, every one of its
\(\binom84/2=35\) unordered `4+4` partitions produces such an equal-signature
pair.  The omitted zero fibre causes no gap: the total support weight is even,
so even parity in all fifteen nonzero fibres forces even parity in the zero
fibre as well.

Every candidate union must then be replayed directly against all 32 direction
groups and retained only if at least nine are silent.  That replay, together
with the four generic cross-ratio cases, would be an exact weight-eight
certificate.

A simple aligned record containing a 64-bit signature and packed indices uses
16 bytes, or about 32.55 GiB for one case.  Alternatively, a 64-bit key array
and 32-bit combinadic-rank array use about 24.41 GiB for one bank.  These are
resource estimates only; no enumeration was launched for this note.

## Executable theorem record

- `src/e1_gmin_m4_p31_group_s8_cross_ratio.py`
- `tests/test_p31_group_s8_cross_ratio.py`

The test checks the exact orbit partition, the `6>5` coverage argument, the
`4+4` counts, and explicit non-closure flags.  It is a transcription guard,
not evidence that the proposed large certificate has been run.
