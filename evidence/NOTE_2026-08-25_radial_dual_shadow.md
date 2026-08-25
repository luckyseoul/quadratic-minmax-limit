# Radial dual-shadow transform of the Max+ odd coset

## Result

Retain the notation of Propositions 15.629--15.630:

\[
L=\ker_{\mathbb Z}(C-pI),\qquad
P={1\over2}(I+C/p),\qquad L^*=P\mathbb Z^n,
\]

and let `d=(p^2+1)/2`.  If `y0` is any Max+ vector, then `y0` is an
odd integral vector in `L`.  For `u=Pz` in `L*`,

\[
 \langle u,y_0\rangle=\langle z,y_0\rangle
 \equiv\sum_i z_i\pmod2.
\]

On the other hand,

\[
 2p\|u\|^2=p\|z\|^2+z^TCz\equiv\sum_i z_i\pmod2,
\]

because `p` is odd, `z_i^2 congruent to z_i (mod 2)`, and the symmetric
zero-diagonal matrix `C` has even quadratic form `z^TCz`.  Therefore

\[
 e^{\pi i\langle u,y_0\rangle}=(-1)^{2p\|u\|^2}.             \tag{1}
\]

The dual phase of the odd coset is independent of the chosen Max+ vector
and of every glue label except norm parity.

For every homogeneous harmonic polynomial `H` of degree four and `t>0`,
Poisson summation on `y0+2L` now reads

\[
 \sum_{y\in y_0+2L}H(y)e^{-\pi t\|y\|^2}
 ={t^{-d/2-4}\over\operatorname{vol}(2L)}
 \sum_{u\in L^*}(-1)^{2p\|u\|^2}H(u/2)
 e^{-\pi\|u\|^2/(4t)}.                                    \tag{2}
\]

This is an exact scalar norm-parity twist, rather than an unknown
vector-valued glue-class phase.

## Exact first dual gap

Proposition 15.630 gives

\[
 \|u\|^2={1\over2}\quad\Longleftrightarrow\quad
 u\in\{\pm Pe_i\}.
\]

Every other nonzero dual vector satisfies

\[
 \|u\|^2\ge {p-1\over p}.                                  \tag{3}
\]

For common circle-profile sum zero this is the MDS/Newton bound in
Proposition 15.630.  For nonzero common sum, write `|t|=ap+b` with
`0<=b<p`.  The sharp balancing expression

\[
 (p+1)f_p(t)-t^2=pa^2+2ab+b(p+1-b)
\]

equals `p` only for `|t| in {1,p}`; away from those minimum vectors it is
at least `2(p-1)`.

## Degree-four first-shell coefficient

Let `W` be symmetric on `V+`, with

\[
 PWP=W,\qquad \operatorname{diag}W=0,\qquad\operatorname{tr}W=0,
\]

and define

\[
 H_W(x)=(x^TWx)^2-{4\over d+4}\|x\|^2x^TW^2x
 +{2\operatorname{tr}(W^2)\over(d+2)(d+4)}\|x\|^4.
\]

Direct differentiation shows that `H_W` is harmonic.  Since the `Pe_i`
form a tight frame, have squared norm `1/2`, and satisfy
`(Pe_i)^TW(Pe_i)=W_ii=0`,

\[
 \sum_{u\in\{\pm Pe_i\}}H_W(u)
 =-{2\over d+2}\|W\|_F^2.                                 \tag{4}
\]

The minimum shell has phase `(-1)^p=-1`, while homogeneity gives
`H_W(u/2)=H_W(u)/16`.  Its contribution to the transformed sum in (2) is
therefore

\[
 {1\over8(d+2)}\|W\|_F^2.                                  \tag{5}
\]

## What this does and does not close

Equations (1)--(5) remove the unknown phase and determine the first dual
shell, its gap, and its degree-four coefficient exactly.  They do **not**
yet prove R1: degree-four harmonic sums on the higher dual shells have no
proved sign, and a single Gaussian parameter cannot simultaneously isolate
the first primal odd-coset shell and the first dual shell.

The mechanism is analogous to characteristic-vector and shadow arguments
for odd unimodular lattices, especially Elkies and Rains--Sloane.  Those
results do not apply verbatim here: `L` is non-unimodular of growing level
`4p`, and the needed coefficient is a degree-four harmonic coefficient of a
specific odd coset.  The exact norm-parity identity above is the Paley
lattice input that the general shadow formalism does not supply.

No R1, global-QVAR, or final-limit flag is flipped by this proposition.
