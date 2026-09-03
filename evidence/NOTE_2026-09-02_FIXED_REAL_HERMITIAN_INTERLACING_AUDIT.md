# Fixed-real-part Hermitian interlacing audit

**Status:** exact positive adaptation for the largest eigenvalue; exact
obstruction to the spectral-radius step.  This is an all-order argument and
does not use a coefficient-signing census.  It does not prove the multiplier-
two estimate.

The primary source is Gary Greaves, Bojan Mohar, and Suil O,
[“Interlacing families and the Hermitian spectral norm of digraphs,”
arXiv:1602.06274](https://arxiv.org/abs/1602.06274).  The relevant ingredients
are its Lemma 2.1 (the matching-polynomial expectation), Lemma 3.3 (real
rootedness for mixtures of rank-one Hermitian updates), Theorem 3.4
(interlacing), and the symmetry step immediately before Corollary 1.2.

## 1. The fixed real part does not break one-sided interlacing

Let `A` be a real symmetric signing with zero diagonal on a graph `G`, and
choose independent skew signs `R_uv=-R_vu in {+1,-1}` on its edges.  Put

\[
 H_R=A+iR,
 \qquad f_R(x)=\det(xI-H_R).
\]

For every edge `e={u,v}`, choose once and for all an ordering `(u,v)`.  If
`sigma=R_uv`, define

\[
 w_{e,\sigma}
 =2^{1/4}\left(e_u+{-A_{uv}+i\sigma\over\sqrt2}e_v\right).
\]

Since `A_uv^2=1`, direct multiplication gives

\[
 (w_{e,\sigma}w_{e,\sigma}^*)_{uv}=-A_{uv}-i\sigma,
 \qquad
 (w_{e,\sigma}w_{e,\sigma}^*)_{uu}
 =(w_{e,\sigma}w_{e,\sigma}^*)_{vv}=\sqrt2.
\]

Consequently, with `D_G` the degree matrix,

\[
 \sum_{e\in E(G)}w_{e,R_e}w_{e,R_e}^*=\sqrt2D_G-H_R,              \tag{1}
\]

and hence

\[
 f_R(x)=\det\left(xI-\sqrt2D_G+
              \sum_e w_{e,R_e}w_{e,R_e}^*\right).                \tag{2}
\]

This is exactly the rank-one form used in Greaves--Mohar--O.  More
explicitly, set `Delta=sqrt(2) max_v deg(v)`.  Every convex combination at
every node of the sign-decision tree is, up to a positive scalar and the
translation `y=x-Delta`, a polynomial of the form

\[
 \sum_{S\subseteq E}\prod_{e\in S}p_e\prod_{e\notin S}(1-p_e)
 \det\left(yI+(\Delta I-\sqrt2D_G)
       +\sum_{e\in S}w_{e,+}w_{e,+}^*
       +\sum_{e\notin S}w_{e,-}w_{e,-}^*\right).                 \tag{3}
\]

Here fixed decisions have `p_e` equal to zero or one, the next decision has
arbitrary `p_e in [0,1]`, and undecided signs have `p_e=1/2`.  The fixed
matrix `Delta I-sqrt(2)D_G` is positive semidefinite, so Lemma 3.3 of the
paper makes (3) real-rooted.  The common-interlacing criterion used in its
Theorem 3.4 therefore applies verbatim.

Thus

\[
 \boxed{\{\det(xI-A-iR):R\text{ skew signing}\}
        \text{ is an interlacing family}.}                       \tag{4}
\]

In particular, if

\[
 p_A(x):=\mathbb E_R\det(xI-A-iR),
\]

then `p_A` is real-rooted and some `R` satisfies

\[
 \lambda_{\max}(A+iR)\le \lambda_{\maxroot}(p_A).                \tag{5}
\]

The adaptation is therefore not blocked at the interlacing lemma.  What is
lost is the simple expected polynomial and, more importantly, the passage
from one spectral edge to both spectral edges.

## 2. Exact expected characteristic polynomial

For a vertex set `W`, write `A[W]` for the corresponding principal
submatrix.  Expanding the determinant by permutations gives

\[
 \boxed{
 p_A(x)=\sum_{M\text{ a matching in }G}(-1)^{|M|}
 \det\bigl(xI-A[V\setminus V(M)]\bigr).
 }                                                               \tag{6}
\]

The empty determinant is one.  Here is a direct proof of all factors.  A
fixed point of a permutation contributes `x`.  On a transposition `{u,v}`,

\[
 (-H_{uv})(-H_{vu})=|A_{uv}+iR_{uv}|^2=2,
\]

instead of the value `A_uv^2=1` in `det(xI-A)`.  In every permutation cycle
of length at least three, each unoriented edge variable occurs only once, so
independence and `E R_uv=0` replace every factor `A_uv+iR_uv` by `A_uv`.
Thus a permutation with transposition set `T` has twice its `A`-determinant
weight for each edge of `T`.  Expanding `2^{|T|}` as the number of subsets
`M subseteq T` gives (6): first prescribe the extra transpositions `M`, then
take the determinant on the unused vertices.

When `A=0`, (6) is the ordinary matching polynomial, exactly Lemma 2.1 of
the paper.  For nonzero `A`, cycles of length at least three survive.  Hence
`p_A` depends on the signed cycle products of `A`; it is not the matching
polynomial and the matching-polynomial/universal-cover root estimate in the
paper does not transfer.

## 3. Exact failure of the spectral-radius conclusion

Greaves--Mohar--O first choose an orientation with a controlled largest
eigenvalue.  They then obtain spectral-radius control because the spectrum
of the pure matrix `iR` is symmetric about zero.  The fixed real part destroys
that symmetry: in general `A+iR` is not cospectral with `-(A+iR)`.

This is a real obstruction, not just a missing argument.  Take `G=K_3` and
`A=J-I`.  Formula (6) gives

\[
 p_A(x)=\det(xI-A)-3x=x^3-6x-2.                                 \tag{7}
\]

For any skew signing, writing
`(R_12,R_13,R_23)=(r,s,t)`, direct expansion gives

\[
 \det(A+iR)=2(1-rt+rs+st)\in\{4,-4\}.
\]

Also `tr((A+iR)^2)=12`, so every leaf polynomial is one of

\[
 x^3-6x-4=(x+2)(x^2-2x-2),
 \qquad
 x^3-6x+4=(x-2)(x^2+2x-2).                                     \tag{8}
\]

Every orientation therefore has

\[
 \|A+iR\|_{\mathrm{op}}=1+\sqrt3.                               \tag{9}
\]

On the other hand, for `q(x)=x^3-6x-2`,

\[
 q(-(1+\sqrt3))=-6<0<q(-\sqrt2),
 \qquad
 q(\sqrt2)<0<q(1+\sqrt3)=2.
\]

Using the monotonicity intervals cut out by `q'(x)=3x^2-6` shows that every
root of `p_A` lies strictly inside `(-(1+sqrt(3)),1+sqrt(3))`.  Thus

\[
 \min_R\|A+iR\|_{\mathrm{op}}>\rho(p_A).                         \tag{10}
\]

So there is no fixed-real-part analogue of Corollary 1.2 obtained by
replacing “largest root” with “largest absolute root.”  Applying the
one-sided theorem separately to `A` and `-A` merely produces two potentially
different signings and cannot control both ends at one leaf.  Example (10)
shows that the desired simultaneous selection is false for the expected
characteristic polynomial itself.

## 4. Why a full spectral bound is stronger than the cube target

For the fourth-phase vector `w=(x+iy)/sqrt(2)`, one has `||w||_2^2=n` and

\[
 B(A,R)=\max_w|w^*(A+iR)w|\le n\|A+iR\|_{\mathrm{op}}.            \tag{11}
\]

Therefore the zero-leading-error directed-half-cut target would follow from

\[
 \|A+iR\|_{\mathrm{op}}\le {2\sqrt2\,M\over n},
 \qquad M=\Phi(A).                                               \tag{12}
\]

But every such completion has

\[
 \|A+iR\|_F^2=2n(n-1),
 \qquad
 \|A+iR\|_{\mathrm{op}}\ge\sqrt{2(n-1)}.                        \tag{13}
\]

Combining (12)--(13), including an `o(sqrt(n))` Dini error, forces

\[
 M\ge {n\sqrt{n-1}\over2}-o(n^{3/2}).                            \tag{14}
\]

For an optimal `A`, the known upper bound is
`M<=(1/2+o(1))n^(3/2)`.  Hence a full operator-norm theorem strong enough for
the half-cut target would already force the sharper conclusion
`m_n/n^(3/2) -> 1/2`, rather than merely prove convergence.  It is
Frobenius-critical: at the upper scale, (12) asks for
`(sqrt(2)+o(1))sqrt(n)`, the smallest possible leading constant in (13).
Even the original pure-orientation universal-cover scale for `K_n` is
`2sqrt(n-2)`, already too large, and the fixed-real-part expectation has no
matching-polynomial estimate at all.

The usable conclusion is therefore exactly one-sided: the fixed-real-part
characteristic polynomials interlace and (6) is their real-rooted average.
Neither result controls the two spectral edges needed for (11), and replacing
the fourth-phase cube by the full sphere imposes a substantially stronger,
value-determining demand.
