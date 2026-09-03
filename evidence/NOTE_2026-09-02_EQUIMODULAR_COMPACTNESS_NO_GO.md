# Equimodular compactness and block-blow-up no-go

Date: 2026-09-02

Status: **proved obstruction to the natural weighted-compactness argument;
the limit of `K_n` remains open.**  This note uses no finite census.

Write

\[
 Q_c(x)=\sum_{i<j}c_{ij}x_ix_j,\qquad
 \Phi(c)=\max_{x\in\{\pm1\}^n}|Q_c(x)|,
\]

and let

\[
 I_1(Q_c)=\sum_i {\bf E}\left|\sum_{j\ne i}c_{ij}x_j\right|.
\]

For complete equimodular forms the sharp constant is

\[
 K_n=\max_{|c_{ij}|=1}{I_1(Q_c)\over\Phi(c)}
     ={n\mu_{n-1}\over m_n},
 \qquad \mu_r={\bf E}|\epsilon_1+\cdots+\epsilon_r|.
\]

The question tested here is whether one can put these forms inside a nested
weighted compact class, obtain a limit there, and recover `K_n` by padding or
by filling the zero coefficients with pseudorandom sign blocks.

## 1. The nested weighted relaxation has a limit

Define

\[
 U_n=\sup_{0\ne c\in\mathbb R^{E(K_n)}}
             {I_1(Q_c)\over\Phi(c)}.                         \tag{1}
\]

After imposing `Phi(c)=1`, this is a maximum over a compact
finite-dimensional set.  Extending a coefficient array by one zero row and
column preserves both `I_1` and `Phi`, so

\[
 U_n\le U_{n+1}.                                             \tag{2}
\]

It is useful to record a self-contained uniform bound.  For fixed `x`, put
`h_i=sum_(j ne i)c_ij x_j` and choose `y_i=sgn(h_i)`.  If
`u=(x+y)/2` and `v=(x-y)/2`, then `u,v` belong to `{0,+-1}^n` and

\[
 \sum_i|h_i|=y^TCx=2Q_c(u)-2Q_c(v),                          \tag{3}
\]

where `C` is the symmetric zero-diagonal matrix with off-diagonal entries
`c_ij`.  Every value of `Q_c` on `{0,+-1}^n` is a conditional expectation of
values on the Boolean cube, and therefore has absolute value at most
`Phi(c)`.  Averaging (3) gives

\[
 I_1(Q_c)\le4\Phi(c).                                       \tag{4}
\]

Consequently `U_n` increases to a finite limit `U_infty in [2,4]`; the lower
bound comes from the single monomial `x_1x_2`, whose ratio is exactly two.
Thus weighted relaxation really does manufacture a convergent sequence.
What fails is the return from that relaxation to the complete flat slice.

## 2. Compactness loses every normalized flat extremizer

Let `A_n` be any complete sign form and normalize it by
`F_n=Q_(A_n)/Phi(A_n)`.  The coefficient `ell_2` norm is

\[
 \|\widehat F_n\|_2={\sqrt{\binom n2}\over\Phi(A_n)}.
\]

The universal bound already proved in the repository,

\[
 \Phi(A_n)\ge {n\sqrt{n-1}\over\pi},
\]

therefore gives the exact decay

\[
 \|\widehat F_n\|_2\le {\pi\over\sqrt{2n}}\longrightarrow0. \tag{5}
\]

In particular, after zero-padding into one infinite coefficient space, every
normalized flat sequence converges strongly in coefficient `ell_2` (and
coordinatewise) to zero.  On the other hand, if `A_n` is optimal then

\[
 I_1(F_n)=K_n
 \ge \sqrt{{n\over2(n+2)\log2}},                             \tag{6}
\]

using the elementary random-sign upper bound for `m_n` and
`mu_(n-1)>=sqrt((n-1)/2)`.  The right side has a positive limit.

Equations (5)--(6) are a direct discontinuity certificate: total `L^1`
influence is not continuous on the coefficient compactification in which the
flat forms acquire a limit.  A topology strong enough to retain `I_1` cannot
make this escaping sequence compact without adding extra fluctuation data.

This is also a literal failure of density in the usual function norms.  Embed
the weighted extremizer `g=x_1x_2` in every dimension.  Walsh orthogonality
and (5) show, for every normalized complete flat form `F_n`,

\[
 \|F_n-g\|_2^2
 =1+\|F_n\|_2^2-2\widehat F_n(\{1,2\})\longrightarrow1.       \tag{7}
\]

Hence normalized complete flat forms do not approximate even this fixed
weighted form in `L^2`, and therefore not in uniform norm.  Equality of the
*extremal values* `K_n` and `U_n` could still occur accidentally, but it
cannot follow from density or ordinary compactness.

## 3. A zero-mean sign block is a leading-order object

For every `r by s` sign matrix `B`, define its switching norm by

\[
 \|B\|_{\infty\to1}
 =\max_{u\in\{\pm1\}^r,v\in\{\pm1\}^s}|u^TBv|.
\]

For random `v`, choose `u_i=sgn((Bv)_i)`.  Each row dot product is a sum of
`s` Rademacher variables, independently of the signs in that row.  Therefore

\[
 \boxed{\quad
 \|B\|_{\infty\to1}\ge r\mu_s,
 \qquad \|B\|_{\infty\to1}\ge s\mu_r.\quad}                 \tag{8}
\]

This lower bound is completely unaffected by requiring zero total sum, zero
row and column sums, Hadamard orthogonality, or pseudorandomness.

If `B` is the cross block of a complete quadratic form `F` between disjoint
vertex sets, then

\[
 \Phi(F)\ge\|B\|_{\infty\to1}.                              \tag{9}
\]

Indeed, fix `u,v`, average every outside coordinate, and compare the two
assignments `(u,v)` and `(u,-v)`.  The within-block terms are identical, the
outside terms average to zero, and half the difference is `u^TBv`.

Thus a nominally zero-weight macro-edge filled by a `t by t` sign block
already forces

\[
 \Phi(F)\ge t\mu_t
          =(\sqrt{2/\pi}+o(1))t^{3/2}.                       \tag{10}
\]

For a blow-up of a fixed weighted seed this is a fixed nonzero fraction of
the critical `(nt)^(3/2)` scale; it is not `o(t^(3/2))`.  More generally, a
balanced `R by S` rectangle of filler signs forces
`Phi(F)>=max(R mu_S,S mu_R)=Theta((R+S)^(3/2))`.

There is a second, weighted version of the same obstruction.  Khintchine's
inequality followed by (4) gives every symmetric zero-diagonal coefficient
array `R` the bound

\[
 \boxed{\quad
 \Phi(R)\ge {1\over4\sqrt2}
   \sum_i\left(\sum_{j\ne i}r_{ij}^2\right)^{1/2}.\quad}      \tag{11}
\]

In the usual biased-sign encoding of a macro coefficient, block means must
be `O(t^(-1/2))` to keep the coherent part on the critical `t^(3/2)` scale.
After subtracting those means, the residual coefficients still have
absolute value `1-o(1)` on a dense set.  Equation (11) says that this
residual itself has norm `Theta(N^(3/2))`.  Consequently it cannot be put in
a triangle-inequality error term.  Correlating it with the coherent part may
help, but designing precisely that cancellation is again the original flat
minimax problem, not a consequence of the weighted seed.

The alternative scaling fails even more directly.  If a sign block has a
nonzero constant mean, its all-one bilinear value is `Theta(t^2)`, so (9)
makes the normalized total influence of a fixed-size blow-up tend to zero.
Hence block homogenization faces an exact dichotomy:

* constant block means create a fatal `t^2` coherent value;
* critical `t^(-1/2)` means leave `t^(3/2)` filler discrepancy.

There is no perturbative regime between them that preserves the sharp
influence ratio.

## 4. Precise conclusion

The unrestricted constants `U_n` are a legitimate dimension-monotone compact
relaxation and they converge.  None of the three desired transfer steps is
valid:

1. zero-padding leaves the equimodular class;
2. normalized flat forms escape to zero while `I_1` stays positive, so the
   relevant functional is discontinuous in the natural compactification;
3. pseudorandom zero-mean completion has unavoidable leading-order switching
   norm, so a block blow-up cannot inherit the weighted ratio with a vanishing
   error.

This does **not** prove that `K_n` diverges, nor does it rule out a tailored
nonperturbative comparison theorem.  It proves that weighted compactness,
ordinary padding, and negligible-filler blow-up do not provide such a
theorem.  A successful compact route must retain a second-order fluctuation
object and prove an all-orders comparison inside the flat class itself.
