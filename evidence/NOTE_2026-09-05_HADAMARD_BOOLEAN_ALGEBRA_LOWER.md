# An invariant Boolean algebra gives actual lower bounds in a modified embedding

2026-09-05. Additive analytic note. This modifies the regular Hadamard
background in the finite-template sparse-flip construction. The new
background has an explicit algebra of actual Boolean vectors fixed by
its normalized operator. Their block means have a dense exact grid.
This proves a genuine Boolean lower bound, not just a lower bound on
the magnetization completion upper Gamma.

The dense paired-mean direction, completion functional, and symmetric
dilation were developed collaboratively with the parent. The exact
worker supplied the invariant-algebra construction and its integration
with the sparse flips. The clipped-Gaussian constant D(q) is reused
from the proof worker's companion finite-template theorem.

## 1. Companion construction and the one change of background

The companion *Actual sign matrices with a finite template on a
Hadamard background*, source
`original_mo_hadamard_sparse_flip_template.md`, proves the sparse-flip
existence estimate and the Boolean upper Gamma. Its frozen SHA256 is
`0d2355f94734b4c1e950c1e05c6df75df38b5ce181ba7fce550a4245e11328ed`.

Here we retain its notation for beta, vector SDP tau, and the
finite-template completion

\[
 \Gamma(C)=\max_{a,b\in[-1,1]^p}\left[
 {a^TCb\over p}
 +\sqrt{\left(1-{\|a\|^2\over p}\right)
        \left(1-{\|b\|^2\over p}\right)}\right].            \tag{1}
\]

We do NOT simply reuse its background `O_p tensor O_m` for the new
lower. That operator mixes blocks, so independently constructing a
paired Boolean vector in each block would not suffice.

Let F=J_4-2I_4, let P_4 be the permutation swapping 1 with 3 and 2
with 4, and put

\[
 H'_4=F P_4,\qquad O'_4=H'_4/2,\qquad
                         w=(1,1,-1,-1)^T.
\]

The matrices F and P_4 commute, P_4 is a symmetric involution, and
F^2=4I_4. Thus O'_4 is symmetric orthogonal and H'_4 has sign entries.
Moreover

\[
             O'_4 1_4=1_4,\qquad O'_4 w=w.                 \tag{2}
\]

Indeed Fw=-2w and P_4w=-w. For n=4^k define

\[
                H'_n=(H'_4)^{\otimes k},\qquad
                O'_n=H'_n/\sqrt n=(O'_4)^{\otimes k}.       \tag{3}
\]

This is again a complete regular real Hadamard sign matrix.

## 2. Fixed Boolean algebra, equal-size blocks, and the exact mean grid

Index coordinates by omega in {1,2,3,4}^k and write
`s_j(omega)=w_(omega_j)` in {-1,1}. Each binary pattern
`(s_1,...,s_k)` has exactly 2^k preimages. Products of the s_j are
fixed by O'_n, by (2), and they span every real function of these k
binary variables. Therefore

\[
 O'_n f(s_1,\ldots,s_k)=f(s_1,\ldots,s_k)
 \quad\hbox{for every real function }f\hbox{ on }\{-1,1\}^k.
                                                               \tag{4}
\]

In particular this includes every Boolean-valued f, so (4) is an
exact sign-to-sign identity, not a Gaussian or sphere relaxation.

Fix p=4^ell, set r=2ell so p=2^r, and let k>=r+1. Partition the n
coordinates by their first r bits. There are p blocks, each of size
m=n/p. Their normalized indicator vectors define an isometry J_n
from R^p to R^n. Since block-constant functions belong to the algebra
(4),

\[
                  O'_n J_n=J_n,\qquad J_n^T O'_nJ_n=I_p. \tag{5}
\]

Inside each block, the remaining k-r bits are uniform on a set of
size L=2^(k-r). A Boolean function on them has precisely the mean grid

\[
                   -1+{2j\over L},\qquad 0\le j\le L.    \tag{6}
\]

These means can be chosen independently in all p blocks. The grid
becomes dense in [-1,1] because k tends to infinity with p fixed.

More generally, choose any two prescribed grid means a_i,b_i in a
block. Order its L binary patterns, and use nested initial segments
as the plus-sets of two Boolean functions f_i,g_i. Their exact overlap
is

\[
                         E_i[f_i g_i]=1-|a_i-b_i|.        \tag{7}
\]

Here E_i is the uniform average in that block; the equal multiplicity
of binary patterns preserves it. Combining these choices across the
p blocks gives actual global sign vectors x,y fixed by O'_n, with
block means a,b, and

\[
                    {x^TO'_n y\over n}
                    ={x^Ty\over n}
                    =1-{\|a-b\|_1\over p}.               \tag{8}
\]

## 3. Actual sparse-flip sign matrices for the modified model

For any fixed real p by p template C, define

\[
                       T'_n=O'_n+J_n(C-I_p)J_n^T.         \tag{9}
\]

This acts as C on block constants and as the orthogonal O'_n on
their perpendicular complement. Thus the spectral decomposition in
the companion construction is unchanged, except that its compressed
background O_p has been replaced by I_p.

For clarity, the exact sparse flip rule also changes accordingly:
put `M=p(C-I_p)`, `gamma=max(1,max_ij|M_ij|)`, and independently flip
the background entry H'_n[a,b], whose block indices are i,j, with
probability

\[
               {\gamma-H'_n[a,b]M_{ij}\over2\sqrt n}.     \tag{10}
\]

For n>=gamma^2 these are valid probabilities. Entrywise expectation
gives `E B_n/sqrt(n)=T'_n-gamma O'_n/sqrt(n)`. Centered entries have
norm at most two and maximum row/column variance sum at most
`4 gamma sqrt(n)`. The SAME rectangular Bernstein and Markov proof
in the companion therefore supplies deterministic actual complete
sign matrices, with O_C(n^(3/2)) flips, satisfying

\[
 \left\|{B_n\over\sqrt n}-T'_n\right\|_{\rm op}
             \le\epsilon_n=n^{-1/8}+\gamma n^{-1/2}
                              \longrightarrow0.          \tag{11}
\]

All statements below concern these actual sign realizations. As in
the companion, if `q=||C||op>=1` and `tau(C)=pq`, then

\[
 {\|B_n\|_{\rm op}\over\sqrt n}\longrightarrow q,
 \qquad {\tau(B_n)\over n^{3/2}}\longrightarrow q,
 \qquad {\tau(B_n)\over n\|B_n\|_{\rm op}}\longrightarrow1. \tag{12}
\]

Their empirical singular-value law after division by sqrt(n) tends
to delta_1, with only the fixed template's non-bulk limiting outliers.
Exact finite-n scalar optimality is not asserted.

## 4. A genuine Boolean lower and an exact PSD identification

For the actual paired signs (8), equation (9) gives

\[
 {x^TT'_n y\over n}
 ={a^TCb\over p}
       +1-{\|a-b\|_1+a^Tb\over p}.                       \tag{13}
\]

Define

\[
 \Lambda_I(C)=\max_{a,b\in[-1,1]^p}
       \left[{a^TCb\over p}
              +1-{\|a-b\|_1+a^Tb\over p}\right].         \tag{14}
\]

The objective is continuous on a compact cube. Its exact grid maxima
converge to this maximum. Equations (11), (13), and the unchanged
orthogonal-complement upper from the companion therefore prove

\[
 \boxed{\displaystyle
 \Lambda_I(C)\le\liminf_n{\beta(B_n)\over n^{3/2}}
 \le\limsup_n{\beta(B_n)\over n^{3/2}}\le\Gamma(C).}       \tag{15}
\]

Unlike a lower bound on Gamma alone, the first inequality in (15)
is witnessed by actual sign vectors. In particular, taking a=b gives

\[
 \Lambda_I(C)\ge
            1+{1\over p}\max_{a\in[-1,1]^p}a^T(C-I_p)a. \tag{16}
\]

There is an exact identification when C is symmetric positive
semidefinite:

\[
 \boxed{\displaystyle
 \Lambda_I(C)=\Gamma(C)
     =1+{1\over p}\max_{a\in[-1,1]^p}a^T(C-I_p)a,
 \qquad {\beta(B_n)\over n^{3/2}}\longrightarrow\Gamma(C).} \tag{17}
\]

Indeed PSD Cauchy--Schwarz gives
`a^TCb <= (a^TCa+b^TCb)/2`; arithmetic-geometric mean gives the
analogous upper for the residual square root in (1). Adding shows
that Gamma is at most the maximum on the right of (17), while
choosing a=b gives the reverse inequality.

For completeness, Lambda_I<=Gamma also holds without PSD. The
Bernoulli overlap in (7) has covariance
`1-|a_i-b_i|-a_i b_i`, at most
`sqrt((1-a_i^2)(1-b_i^2))` by Cauchy--Schwarz. Average over i and
apply Cauchy--Schwarz again to bound it by the square root in (1).
Together with (16) this completes the PSD equality and, by (15),
its actual Boolean-limit assertion.

## 5. A positive quadratic top frame gives an actual clipping lower

Suppose C is symmetric, `q=||C||op>1`, and there exists a p by d
matrix W with unit rows and

\[
                              C W=q W.                   \tag{18}
\]

This is a positive quadratic SDP top frame. It implies tau(C)=pq,
but we do not infer (18) from the bipartite equality alone.

Let G be a standard real Gaussian and let f:R->[-1,1] be measurable
and odd. Put

\[
 c_f=E[Gf(G)],\qquad v_f=E[f(G)^2],\qquad
                         D_f(q)=q(2c_f^2-v_f)-v_f.
\]

Then this modified actual sign family satisfies

\[
 \boxed{\displaystyle
 \liminf_n{\beta(B_n)\over n^{3/2}}
                              \ge1+D_f(q).}               \tag{19}
\]

Here is the full rounding proof. Let g be standard Gaussian in R^d,
let X=Wg, and take a_i=f(X_i), an admissible random vector in (16).
Writing `r(x)=f(x)-c_f x`, Gaussian conditional expectation cancels
the two mixed terms because E[G r(G)]=0. Thus

\[
 E[a^TCa]=c_f^2\operatorname{tr}(W^TCW)
                     +E[r(X)^TCr(X)]
              \ge p q c_f^2-p q(v_f-c_f^2),
\]

using (18), `||C||op=q`, and `E||r(X)||^2=p(v_f-c_f^2)`.
Also `E||a||^2=pv_f`. Taking the expected objective in (16) proves
(19). Indefiniteness of C causes no problem: the residual was bounded
below by -q times its squared norm, not assumed nonnegative.

For `f(x)=clip(x,-1,1)`, write

\[
 P=\Pr\{|G|\le1\},\qquad \phi={e^{-1/2}\over\sqrt{2\pi}}.
\]

Gaussian integration by parts gives `c_f=P`, `v_f=1-2phi` and

\[
             D_f(q)=q(2P^2-1+2\phi)-1+2\phi.             \tag{20}
\]

Sections 5--6 of the companion *A quantitative completion bound for
scalar-optimal finite templates*, source
`original_mo_scalar_template_gamma_bound.md`, establish the exact
elementary enclosures
`0.68268<P<0.68270`, `0.24197<phi<0.24198` and hence

\[
 2P^2-1+2\phi>0,\qquad
 D_f(5/2)=5P^2-7/2+7\phi>0.524049912.                     \tag{21}
\]

The general actual theorem (19) does not depend on decimal bounds.
Reusing those same rational enclosures, rather than a new numerical
calculation, yields the concrete consequence

\[
 \boxed{\displaystyle
 q\ge5/2\quad\Longrightarrow\quad
 \liminf_n{\beta(B_n)\over n^{3/2}}
                           >1.524049912>\sqrt2,}          \tag{22}
\]

under the additional positive-frame condition (18).

## 6. Symmetric dilation supplies that frame, but changes the model

Let C_0 be any real p by p template with p=4^ell, `||C_0||op=q>1`,
and `tau(C_0)=pq`. Choose optimal unit-row frames U,V. Equality in
the operator/Cauchy--Schwarz bound gives

\[
                   C_0V=qU,\qquad C_0^TU=qV.             \tag{23}
\]

Define the symmetric, generally INDEFINITE template

\[
 \widetilde C=\begin{pmatrix}0&C_0\\C_0^T&0\end{pmatrix},
 \qquad W=\begin{pmatrix}U\\V\end{pmatrix}.
\]

Then `||Ctilde||op=q` and `Ctilde W=qW`, with unit rows in W.
Its dimension 2p is not a power of four, so to fit the explicit
background construction use

\[
                         C_* =\widetilde C\otimes I_2,
\]

of dimension 4p=4^(ell+1). This is permutation-equivalent to two
copies of Ctilde. Repeating W accordingly gives a unit-row positive
top frame for C_*, while `||C_*||op=q` and `tau(C_*)=4pq`.

The finite-template Boolean/SDP ratio is preserved by this dilation
and duplication: `beta(Ctilde)=2 beta(C_0)`, since its two bilinear
terms optimize independently, and another factor two comes from the
two copies. Thus `beta(C_*)/(4pq)=beta(C_0)/(pq)`.

Nevertheless, when q>=5/2, the actual sign family implanted from C_*
in the MODIFIED invariant-algebra background has Boolean norm above
sqrt(2)n^(3/2) asymptotically, by (22). This is a genuine actual-norm
exclusion for this specified construction and symmetry class.

It does not say that implanting the original nonsymmetric C_0 has
the same actual Boolean norm, or that all complete sign matrices
with similar spectral bulk are of this form. The symmetric dilation
changes the template and the large sign family; it is not an allowed
identification of arbitrary actual matrices.

## 7. Remaining scope

The invariant Boolean algebra makes the mean grid and overlaps exact,
and (15)--(22) give actual lower bounds for the modified family. The
PSD case even has an identified actual limit. For general nonsymmetric
templates, Lambda_I and Gamma may differ, and no equality is asserted.

No numerical sign construction or optimization was run. No original
internal signing, conditional minimizer, joint-shell realization, or
intrinsic source compatibility is supplied. This additive result does
not prove the original conditional upper or the original MO limit.
