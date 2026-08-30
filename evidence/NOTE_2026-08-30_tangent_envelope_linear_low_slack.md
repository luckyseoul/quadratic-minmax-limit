# Tangent-envelope linear low-slack exclusion

**Date:** 2026-08-30
**Status:** proved narrowing of the first `|D|=p+1` shell; not a close of
that shell, residual (ii), Type I, or the limit.

## Result

For every prime `p>=17`, no normalized `p+1`-point boundary can have outside
pair slack

`1 <= R <= floor((p-4)/3)`.

Together with Proposition 15.724's closure of `R=0`, any positive survivor
must therefore satisfy

`R >= floor((p-1)/3)`.

This supersedes Proposition 15.722's square-root low-slack cutoff in the
current frontier.

## Proof certificate

Write the line-slack contribution as

`h(2r)=r(r-1)` and `h(2r+1)=r^2`.

Choose an inclusion-minimal deletion set `T` such that `A=D\T` is an arc,
and put `t=|T|`. The standard rich-line deletion gives `1<=t<=R`.
Minimality makes every `z in T` incident with at least one `A`-secant.
If `s_z` counts those secants and `I=sum_z s_z`, then linewise

`I <= R`,

because a line with two `A`-points and `u` deleted points contributes `u`
to `I`, while `h(2+u)>=u`.

The arc has size `|A|=p+1-t=p+2-(t+1)`. The hypothesis
`3t<=3R<=p-4` is stronger than the size condition in the odd-order
tangent-envelope theorem. That theorem gives a nonzero dual polynomial of
degree `2(t+1)` whose value in the pencil at `P in A` is the square of the
product of the tangent forms at `P`.

For `z in T`, exactly `|A|-2s_z` lines through `z` are `A`-tangents. More
than `2(t+1)` of them would make the envelope vanish identically on the dual
line `z*`. An `A`-secant through `z`, guaranteed by minimality, gives a point
of `z*` where the envelope is nonzero. Hence

`s_z >= (p-1-3t)/2`

and therefore

`I >= F(t)=t(p-1-3t)/2`.

The quadratic `F` is concave on `1<=t<=R`, so its minimum is at an endpoint.
From `3R<=p-4`,

`F(1)=(p-4)/2 >= 3R/2 > R`

and

`F(R)=R(p-1-3R)/2 >= 3R/2 > R`.

Thus `I>R`, contradicting `I<=R`.

## Scope and source

The proof is symbolic and all-prime in its stated range; it does not infer a
general theorem from a finite census. The small table emitted by the evidence
module checks only arithmetic consequences and sample cutoffs.

The finite-geometry input is the odd-order polynomial tangent envelope in
S. Ball and M. Lavrauw, [*Planar arcs*](https://arxiv.org/abs/1705.10940),
J. Combin. Theory Ser. A **160** (2018), 261--287. It is Theorem 11 in
arXiv v4; the strict inequality available here also meets the hypothesis in
the authors' current renumbered manuscript.

## Artifacts

- `src/e1_gmin_m4_prop15726.py`
- `tests/test_prop15726.py`
- `evidence/e1_gmin_m4_prop15726.json`
