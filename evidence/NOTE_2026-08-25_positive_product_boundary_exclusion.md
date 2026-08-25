# Positive-product infinity-point boundary exclusion

Date: 2026-08-25. This is Proposition 15.643. Under the residual affine
hypotheses, it excludes `D={infinity,v}`, `c_H=+1` for every odd `p>=17`.
The four smaller primes and the negative-product branch remain open, so this
does not close residual (ii), Type I, R1, global QVAR, or the limit.

## 1. Pointwise baseline

Let `H` have `h=4p+1` edges and odd-degree boundary
`D={infinity,v}`. Proposition 15.642 proves that `c_H=+1` consumes the
complete type-split budget with the parity baseline. Thus in every
projective direction `d`,

\[
A_d(x)=x_{s_d(v)},\qquad \epsilon_d S_H(x)=3+2x_{s_d(v)}. \tag{1}
\]

Write `z_s=2x_s-1`; on the middle slice, `sum_s z_s=1`, and (1) becomes

\[
\epsilon_dS_H(z)=4+z_j,                                 \tag{2}
\]

where `j=s_d(v)`.

## 2. The additive inter-fibre matrix

Fix `d`. Let

- `I` be the number of infinity edges of `H` (independent of `d`);
- `n_s` be the number of those edges entering fibre `s`, so `sum n_s=I`;
- `P_d` be the number of finite `H`-edges parallel to `d`;
- `K_{st}` be the signed Paley sum of finite `H`-edges between distinct
  fibres `s,t`.

Parallel finite edges have Paley sign `epsilon_d`. Therefore

\[
\epsilon_dS_H(z)=P_d+\sum_sn_sz_s
 +\epsilon_d\sum_{s<t}K_{st}z_sz_t.                     \tag{3}
\]

A quadratic vanishing on the middle slice is a multiple of
`sum x_s-m` by a linear polynomial. Comparing the cross coefficients in
(2)--(3), equivalently taking all four-point second differences, shows that

\[
K_{st}=u_s+u_t.                                          \tag{4}
\]

Using `sum z_s=1` and `z_s^2=1`,

\[
\sum_{s<t}(u_s+u_t)z_sz_t=\sum_su_s(z_s-1).             \tag{5}
\]

Comparison with (2) now gives a scalar `c_d` with

\[
u_s=\epsilon_d(c_d+\mathbf1_{s=j}-n_s),\qquad
(p-1)c_d=I+P_d-5.                                       \tag{6}
\]

Since every `K_{st}` is an integer, (6) gives the key divisibility

\[
{2(I+P_d-5)\over p-1}\in\mathbb Z.                     \tag{7}
\]

## 3. Global parallel-count arithmetic

Put `q=(p-1)/2`. Equation (7) says all `P_d` have one residue `r` modulo
`q`. Write

\[
P_d=r+qk_d,\quad 0\le r<q,\qquad I=5-r+qk_0.            \tag{8}
\]

Every finite edge is parallel in exactly one of the `p+1` projective
directions, so `sum_d P_d=4p+1-I`. Substitution into (8) gives

\[
q(8-k_0-\sum_dk_d)=pr.                                  \tag{9}
\]

Because `gcd(p,q)=1`, (9) implies `q` divides `r`; the chosen range forces
`r=0`. Hence

\[
\boxed{I=5+qk_0,\qquad P_d=qk_d,\qquad
       \sum_dk_d=8-k_0.}                                \tag{10}
\]

For `p>=17`, nonnegativity of `I` and the last identity give
`0<=k_0<=8` and `k_d>=0`.

## 4. Directional l1 obstruction

From (6) and `c_d=(k_0+k_d)/2`, put
`w_s=n_s-1_{s=j}`. Then

\[
|K_{st}|=|k_0+k_d-w_s-w_t|,qquad
\sum_sw_s=I-1=4+qk_0.                                  \tag{11}
\]

The number of transverse selected edges in direction `d` is

\[
4p+1-I-P_d=q(8-k_0-k_d).                                \tag{12}
\]

Since an integer signed sum has absolute value at most the number of terms,
(11)--(12) imply

\[
\begin{aligned}
q(8-k_0-k_d)
&\ge\sum_{s<t}|K_{st}|\\
&\ge\left|\sum_{s<t}(k_0+k_d-w_s-w_t)\right|\\
&=q|k_0+pk_d-8|.
\end{aligned}                                           \tag{13}
\]

If `k_d>=1` and `p>=17`, the quantity inside the absolute value is positive,
and (13) yields

\[
2k_0+(p+1)k_d\le16,                                     \tag{14}
\]

which is impossible. Thus every `k_d=0`. Equation (10) then forces
`k_0=8`, `I=4p+1`, and no finite edges. The only arithmetic endpoint is an
all-infinity star with `4p+1` leaves. Its odd-degree boundary has size
`4p+2`, not two. This contradiction proves the branch empty.

## 5. Boundary of the result

The inequality permits populated directions at `p=5,7,11,13`; they are not
claimed empty. The `c_H=-1` branch has at most three exceptional directions
of each quadratic type by Proposition 15.642, but its baseline equation is
`4-z_j` and requires a separate global count. Other nonempty boundary
profiles also remain.

## 6. Literature and OEIS check

Searches for this repeated additive inter-fibre/rank-two obstruction found
general Boolean-slice harmonic and polynomial-distance frameworks, but no
matching finite-affine edge-count theorem. The only external input remains
the slice-distance lemma already cited in Proposition 15.642; Proposition
15.643 itself is the elementary coefficient, divisibility, and `l1`
argument above. The resulting quantities are closed linear formulas, not an
unidentified sequence requiring an OEIS claim.

## Reproduction

- arithmetic certificate: `src/e1_gmin_m4_prop15643.py`;
- generated evidence: `evidence/e1_gmin_m4_prop15643.json`;
- tests through odd `p<=201`: `tests/test_prop15643.py`.
