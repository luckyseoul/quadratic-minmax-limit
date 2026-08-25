# Negative-product infinity-point boundary exclusion

Date: 2026-08-25. This is Proposition 15.646. It excludes the
`D={infinity,v}`, `c_H=-1` branch for every sufficiently large odd prime by
combining Proposition 15.644's normal form with one signed-sum identity. It
does not give an explicit threshold, handle the positive branch at
`p=5,7,11,13`, classify other boundaries, or close residual (ii), Type I,
R1, global QVAR, or the limit.

## 1. Input normal form

Proposition 15.644 proves that, for every sufficiently large odd prime,

\[
 I=2p-1,\qquad E=2p+2,\qquad P_d=2
\]

in each of the `p-1` baseline directions. There is one exceptional direction
of each quadratic type. If `U,V` denote their positive- and negative-type
parallel counts, respectively, then

\[
 (U,V)=(3,1)\quad\hbox{or}\quad(1,3).                 \tag{1}
\]

There are `q=(p-1)/2` baseline directions of each type.

## 2. Every baseline transverse sum is zero

Fix a baseline direction `d` of type `epsilon_d`. With the notation of
Propositions 15.643--15.645, coefficient comparison gives

\[
 K_{st}=-\epsilon_d(a_s+a_t),\qquad \sum_s a_s=0.      \tag{2}
\]

Here `K_st` is the signed Paley sum of selected finite edges joining fibres
`s,t`. Summing (2) over unordered fibre pairs gives

\[
 \sum_{s<t}K_{st}
 =-\epsilon_d(p-1)\sum_sa_s=0.                         \tag{3}
\]

The left side is exactly the signed sum of finite selected edges transverse
to `d`. Thus every baseline direction requires transverse signed sum zero.
No classification of the individual fibre profile is needed.

## 3. Global sign bookkeeping contradicts zero

Every finite edge parallel to a direction of type `epsilon` has Paley sign
`epsilon`: its difference is an `F_p^*` multiple of a kernel generator, and
every element of `F_p^*` is a square in `F_{p^2}`. The two positive and two
negative edges in each pair of baseline directions cancel in the global
finite-edge signed sum. Consequently

\[
 W:=\sum_{e\text{ finite}}C_e=U-V.                    \tag{4}
\]

If `(U,V)=(3,1)`, then `W=2`. In any negative-type baseline direction its
two parallel edges have total sign `-2`, so the transverse signed sum is

\[
 W-(-2)=4,                                              \tag{5}
\]

contradicting (3). If `(U,V)=(1,3)`, then `W=-2`. In any
positive-type baseline direction its two parallel edges have total sign
`+2`, so the transverse signed sum is

\[
 W-(+2)=-4,                                             \tag{6}
\]

again contradicting (3). Since `q>=1`, the required baseline type exists in
both cases. This excludes the complete normal form and hence the
negative-product infinity-point boundary for all sufficiently large primes.

## 4. Computational discovery and independent checks

An exploratory CP-SAT model first showed that the ideal/one-transfer fibre
marginals alone admit non-two-line point sets at `p=5,7`, while the full
signed inter-fibre system is infeasible in every normalized case. Removing
the boundary and edge-product constraints preserved infeasibility; removing
the exact `K` identities did not. Summing those identities exposed (3)--(6).
The computation is discovery evidence only; the proof above is elementary
and independent of CP-SAT.

## 5. Literature and OEIS check

Searches after the finding located standard Paley graph/conference-matrix
sources and work connecting affine directions to square-order Paley graphs,
but no prior statement matching this residual signed transverse-sum
obstruction. OEIS searches for the elementary counts `2p-1`, `2p+2`, and
the split `(3,1)` returned unrelated linear sequences; this proposition does
not define a new integer sequence.

## Reproduction

- arithmetic verifier: `src/e1_gmin_m4_prop15646.py`;
- generated evidence: `evidence/e1_gmin_m4_prop15646.json`;
- tests: `tests/test_prop15646.py`;
- exploratory finite scouts: `scripts/residual_two_line_fibre_cpsat.py` and
  `scripts/residual_negative_full_cpsat.py`.
