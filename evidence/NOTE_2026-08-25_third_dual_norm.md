# The third dual norm and the exact p=11 third shell

Date: 2026-08-25. Proposition 15.635. This proves the next dual norm for
every `p>=11` and classifies its complete shell at `p=11`. It does not claim
the all-prime shell classification, control later shells, or prove R1.

Write `s=2p||x||^2`, and retain the integral circle profiles

\[
 a_{j,r}=\langle x,v_{j,r}\rangle,\qquad
 t=\sum_r a_{j,r},\qquad
 \sum_{j,r}a_{j,r}^2={s+t^2\over2}.                \tag{1}
\]

For `|t|=ap+b`, `0<=b<p`, balancing one length-`p` profile gives

\[
 g_p(t)=(p+1)f_p(t)-t^2
       =pa^2+2ab+b(p+1-b).                         \tag{2}
\]

After the equality values `|t|=1,p` giving `s=p`, the next odd values of
(2) are `|t|=3,p-2`, giving `3p-6`. It remains to exclude a low-energy
excitation of the `t=1,p` equality cases.

For `t=1`, the degree-one glue makes
`mu_j=sum_r r a_{j,r}` the evaluations of one linear form. Let `u` be the
corresponding point and subtract the minimum vector `Pe_u`. The resulting
zero-sum profiles `b_j=a_j-delta_{mu_j}` are active in, say, `h` directions.
Put

\[
 \Delta={1\over2}\left(\sum_{j,r}a_{j,r}^2-{p+1\over2}\right),
 \qquad s=p+4\Delta.
\]

If `M_j` is the positive mass of `b_j`, then

\[
 \sum_r a_{j,r}^2-1
 =\sum_r b_{j,r}^2+2b_{j,\mu_j}\ge2M_j-2.
\]

Indeed `sum b_j^2>=sum |b_j|=2M_j`; if the distinguished entry is
`-c<0`, its square supplies the extra `c(c-1)`, and
`c(c-1)-2c>=-2`. Thus `Delta_j>=M_j-1`, while every active original
profile costs at least two extra squares, so `Delta>=h`. The MDS/Newton
argument of Proposition 15.630 gives

\[
 \sum_jM_j\ge h(R-h),\qquad R={p+1\over2},\quad h<R.
\]

Therefore

\[
 \Delta\ge\max\{h,h(R-h-1)\}\ge R-2,
\]

since `h(R-h-1)-(R-2)=(h-1)(R-h-2)`. The `h=R` case is larger. For
`t=p`, subtracting `Pe_infinity` gives the same conclusion directly.
Consequently every odd-phase vector outside the minimum shell satisfies

\[
 \boxed{s\ge3p-6.}                                  \tag{3}
\]

For even `t`, (1) and the zero-sum parity show that the second norm
`2(p-1)` can increase only in steps of four; (2) has no intervening even
value. Since `3p-6>2(p+1)` for `p>=11`, the next possible norm is

\[
 \boxed{s_3=2(p+1),\qquad ||x||^2={p+1\over p}.}     \tag{4}
\]

It is attained by

\[
 x=\pm P(e_i+C_{ij}e_j),\qquad i<j,                 \tag{5}
\]

giving `p^2(p^2+1)` signed vectors. These are distinct for `p>=5`: a
collision would give a nonzero `-p` eigenvector of `l1` norm at most four,
whereas the coordinate maximum bound forces `l1>=p+1`. For admissible `W`,
summing the harmonic
polynomial over (5) uses

\[
 \sum_{i<j}(x^TWx)^2=2||W||_F^2,\qquad
 \sum_{i<j}x^TW^2x=p(p+1)||W||_F^2,
\]

and yields the negative scalar

\[
 \boxed{-{p^2+4p-3\over4(p^2+5)}||W||_F^2.}         \tag{6}
\]

At `p=11`, exact PARI/GP `qfminim` on the saturated dual Gram form returned

```
P=11
BOUND=24
SIGNED_COUNT=31110
MAXNORM=24
ELAPSED_MS=215744
```

The proved first and second signed counts are `244` and `16104`; the
remainder is `14762=11^2(11^2+1)`, exactly (5). Thus (5) is the complete
third shell at `p=11`. For `p>11`, (4) and the orbit (5) are uniform
theorems, but completeness of that shell is not asserted.

The bare count is not a new integer sequence: OEIS A071253 records
`n^2(n^2+1)`, and A069187 contains the values `90,650,2450,14762,28730,...`
in an equivalent squarefree-core description. Targeted literature searches
for the Paley ETF lattice and next-shell terminology returned only generic
theta-series and conference-lattice references, not this norm-gap,
harmonic-scalar, or exact shell statement. This is a search record, not an
unqualified novelty claim.
