# Degree-congruence obstruction for the required bi-tight levels

**Status:** required bi-tight levels 2 and 3 proved empty for every prime
`p>=5`; bi-tight level 4 is an additional corollary. No new finite-prime
computation.

Put `q=p^2`, `n=q+1`, and let `C` be the symmetric conference matrix,
`C^2=qI`.  Let `H` have edge indicator `h`, size `|H|=sp`, degrees `d_i`,
and density

\[
a=\frac{|H|}{\binom n2}=\frac{2s}{np}.
\]

Write `kappa=h-a1`, `A=C odot kappa`, and `B=C odot h=A+aC`.  Proposition
15.272 together with 15.207 proves

\[
\ker(G_++G_-)=\mathrm{scheme}\oplus\mathrm{cross}.
\]

If `H` is bi-tight of level `s`, then its centered indicator lies in this
kernel: subtracting `a1` changes the full-edge scores `+pn/2` and `-pn/2`
by exactly `+s` and `-s`.  Hence

\[
A=D_fC+CD_f+X,\qquad \sum_i f_i=0,\qquad CX+XC=0.
\]

Since `XC` is skew-symmetric, its diagonal is zero.  Comparing the diagonal
of `AC` gives

\[
f_i=\frac{d_i-aq}{q-1}.
\]

Absorb `aC` by putting `g_i=f_i+a/2`; then

\[
B=D_gC+CD_g+X,
\qquad
g_i=\frac{d_i-s/p}{q-1}.
\]

The projection onto matrices commuting with `C` is

\[
\operatorname{Comm}(M)=\frac12\left(M+\frac1q CMC\right).
\]

It kills `X` and fixes `D_gC+CD_g`.  Thus, for `i!=j`,

\[
h_{ij}+\frac{C_{ij}(CBC)_{ij}}q=2(g_i+g_j).
\]

The quantity `C_ij(CBC)_ij` is an integer.  Rearranging therefore forces

\[
q-1\mid 2(d_i+d_j)-4ps,
\]

or equivalently

\[
\boxed{d_i+d_j\equiv2ps\pmod{(p^2-1)/2}}.
\]

Subtracting two such equations with one common vertex shows that all degrees
have one common residue modulo

\[
M=\frac{p^2-1}{2}.
\]

For `s=2`, `M>2p` for every `p>=5`.  Since every degree is at most `2p`, all
degrees are equal.  The handshake identity gives

\[
d=\frac{4p}{p^2+1}\in(0,1),
\]

impossible for an integer degree.

For `s=3` and `p>=7`, `M>3p`, so the same argument gives the impossible
integer degree `6p/(p^2+1)`.  At the sole exceptional size comparison
`p=5`, one has `M=12`, `n=26`, and total degree `30`.  Writing every degree
as `r+12m_i` forces `r` to be `0` or `1`, but then respectively `30` or `4`
would have to be divisible by `12`.  Both fail.

Therefore the level-2 and level-3 bi-tight alternatives used by the E(1)
no-descent reductions are empty for every prime `p>=5`.

For `s=4`, the same equal-degree argument works from `p>=11`. At `p=5`,
`M=12`, `n=26`, and total degree is `40`; common residue `0` or `1` leaves
remainder `4` or `2` modulo `12`. At `p=7`, `M=24`, `n=50`, and total degree
is `56`; residue `0` or `1` leaves remainder `8` or `6` modulo `24`. Thus
bi-tight level 4 is also empty for every prime `p>=5`.

This does **not**
claim that every Max-plus-tight cover is empty, nor that every possible
bi-tight level is empty. The E(1) implication chain requires bi-tight levels
2 and 3. In particular, the one-sided Max-plus/Max-minus-tight level-4 cases
inside residual (ii) are **not** excluded by this theorem. The spectral
floor, global QVAR, and principal R1 are not needed for the required
bi-tight gate.
