# Baseline fibre profiles in the negative branch

Date: 2026-08-25. This is Proposition 15.645. It is an exact geometric
refinement inside Proposition 15.644's sufficiently-large-prime normal form.
It does not itself exclude that branch or close any top-level gate;
Proposition 15.646 subsequently bypasses the simultaneous geometry and
excludes the normal form by summing the signed inter-fibre identities.

Let `S` be the `2p-1` finite endpoints of the infinity edges. Fix one of the
`p-1` baseline directions, let `n_s=|S intersect fibre_s|`, and let `j` be
the fibre containing the boundary point `v`. Proposition 15.644 gives

\[
\sum_sn_s=2p-1,qquad P_d=2,qquad E-P_d=2p.              \tag{1}
\]

The baseline target is `4-z_j`. With

\[
w_s=n_s+\mathbf1_{s=j},qquad a_s=w_s-2,                \tag{2}
\]

we have `sum a_s=0`, and the additive signed inter-fibre matrix satisfies

\[
|K_{st}|=|a_s+a_t|.                                     \tag{3}
\]

Its entrywise `l1` norm is at most the number `2p` of transverse selected
edges:

\[
L(a):=\sum_{s<t}|a_s+a_t|\le2p.                         \tag{4}
\]

Put

\[
A=\sum_{a_s>0}a_s=-\sum_{a_s<0}a_s.                    \tag{5}
\]

If there are `u` positive, `v` negative, and `z` zero entries, direct
separation of positive-positive, negative-negative, zero, and mixed pairs
gives

\[
L=(u+v-2+2z)A+\sum_{i=1}^u\sum_{j=1}^v|x_i-y_j|,       \tag{6}
\]

where the positive magnitudes `x_i` and negative magnitudes `y_j` both sum
to `A`. Since all magnitudes are positive integers, `u+v<=min(2A,p)`. Thus

\[
L\ge
\begin{cases}
2A(p-A-1),&2A\le p,\\
(p-2)A,&2A>p.
\end{cases}                                             \tag{7}
\]

For `p>=7`, every `A>=2` makes (7) strictly larger than `2p`. Consequently
`A=0` or `A=1`. In the first case every `a_s=0`. In the second, integrality
forces one entry `+1`, one entry `-1`, and all others zero. Therefore

\[
\boxed{w=(2,\ldots,2)\quad\text{or}\quad
w\text{ is obtained from it by one unit transfer}.}     \tag{8}
\]

Equivalently, every baseline direction sees the ideal intersection profile
`n_j=1`, `n_s=2` for `s!=j`, possibly modified by one transfer after adding
the boundary-fibre unit in (2). The ideal profile is exactly what a union of
two nonparallel affine lines through `v` produces in every other direction.
Proving that the simultaneous one-transfer alternatives collapse to those
two lines remains the next finite-geometry step.

The formulas are verified by `src/e1_gmin_m4_prop15645.py`, with generated
evidence in `evidence/e1_gmin_m4_prop15645.json` and tests in
`tests/test_prop15645.py`.
