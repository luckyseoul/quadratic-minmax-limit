# The exact second shell of the Paley dual lattice

Date: 2026-08-25. Proposition 15.633. This is a proved lattice theorem and
sharpens the radial-shadow attack on R1. It does not control all higher dual
shells and therefore does not by itself prove R1 or global QVAR.

Let

\[
L=\ker_{\mathbb Z}(C-pI),\qquad
P={1\over2}(I+C/p),\qquad L^*=P\mathbb Z^n,
\]

where `C` is the symmetric Paley conference matrix on
`P^1(F_{p^2})`, `n=p^2+1`, and `d=n/2`. For every odd prime `p>=5`,

\[
\boxed{\{u\in L^*: \|u\|^2=(p-1)/p\}
=\mathcal P\mathbin{\dot\cup}\mathcal C,}
\]

where

\[
\mathcal P=
 \{\mathord\pm P(e_i-C_{ij}e_j):i<j\}
\]

and

\[
\mathcal C=
 \{\mathord\pm w_S/p:S\text{ a square }\mathbb F_p\text{-subline}\}.
\]

Here `w_S` is zero on the Miquelian circle `S`, is `+/-1` off `S`,
and satisfies `Cw_S=pw_S`. Consequently

\[
|\mathcal P|=p^2(p^2+1),\qquad
|\mathcal C|=p(p^2+1),
\]

and the signed second-shell count is

\[
\boxed{N_2(p)=p(p+1)(p^2+1).}
\]

The two displayed cardinalities count both signs. At `p=3` the two
candidate descriptions coincide; the signed second shell has size `30`.

## Circle-profile identity and possible common sums

Use the `m=(p+1)/2` square-direction circle words `v_{j,s}` from
Propositions 15.629--15.630. For `x in L*`, set

\[
a_{j,s}=\langle x,v_{j,s}\rangle\in\mathbb Z,
\qquad t=\sum_s a_{j,s}=2p x_\infty.
\]

The sum is independent of `j`, and

\[
\sum_{j,s}a_{j,s}^2=p\|x\|^2+{t^2\over2}.       \tag{1}
\]

For `||x||^2=(p-1)/p`, the right side is `p-1+t^2/2`. If `f_p(t)`
is the minimum squared energy of an integral length-`p` profile with sum
`t`, balancing gives, for `|t|=ap+b`, `0<=b<p`,

\[
(p+1)f_p(t)-t^2=pa^2+2ab+b(p+1-b).              \tag{2}
\]

Equations (1)--(2) imply that the left side of (2) is at most
`2(p-1)`. Equation (1) also forces `t` to be even. Solving these two
conditions gives, up to replacing `x` by `-x`,

\[
\boxed{t\in\{0,2,p-1\}.}                         \tag{3}
\]

For `t=2` or `t=p-1`, equality holds in the balancing bound for every
profile. Thus a `t=2` profile consists of two `+1` entries, while a
`t=p-1` profile consists of one zero and `p-1` entries equal to `+1`.

## The `t=p-1` case

Write the zero of profile `j` as `s_j`. The degree-one profile-glue
condition says that `(s_j)_j` is the evaluation of one linear form at the
projective direction forms `t_j`. Hence there is one point `u` such that
`s_j=t_j(u)` for every `j`. The profiles are exactly

\[
a_{j,s}=1-1_{s=t_j(u)}
=\langle e_\infty-e_u,v_{j,s}\rangle.
\]

Since all circle words span the rational `+p` eigenspace, they determine
`x`, and

\[
x=P(e_\infty-e_u).
\]

This accounts for the `p^2` point pairs containing infinity.

## The `t=0` equality cases

Let `h` be the number of active profiles and let `M_j` be the positive
mass in an active profile. The MDS/Newton argument of Proposition 15.630
gives

\[
M_j\ge m-h,\qquad
\sum_jM_j\ge h(m-h)\ge m-1.
\]

Equality in (1) forces every nonzero profile entry to be `+/-1`, total
positive mass `m-1`, and

\[
h\in\{1,m-1\}.                                   \tag{4}
\]

If `h=1`, the active profile has one zero and `(p-1)/2` entries of each
sign. Its positive and negative root multisets have equal power sums
through degree `(p-3)/2`. Newton identities show that their monic root
polynomials have all coefficients except their constants in common. After
translating the zero to zero, their product is `X^(p-1)-1`. Comparing
coefficients in this product forces the common nonconstant part to be
`X^((p-1)/2)`. Therefore the profile is, up to sign,

\[
a_s=\eta(s-s_0),                                 \tag{5}
\]

where `eta` is the Legendre character of `F_p`, extended by zero.

For the standard square circle `S=P^1(F_p)`, write a field element as
`a+b omega` and define

\[
w_\infty=0,\qquad w_{a+b\omega}=\eta(b).
\]

For fixed transverse coordinates `b,c`, the quadratic-character sum over
`a` is `p-1` when `b=c` and `-1` otherwise. It follows immediately that

\[
Cw=pw.                                             \tag{6}
\]

Moreover, every `z in L` has square-direction profiles modulo `p` of
degree at most `(p-3)/2`. Orthogonality of `eta` to those powers gives
`<w,z>=0 (mod p)`, so `w/p in L*`. Signed PSL transport sends (6) to
every square circle. There are `p(p^2+1)/2` such circles.

If `h=m-1`, every active profile is `delta_alpha-delta_beta`. The
degree-one glue says that the differences `alpha_j-beta_j` are evaluations
of one linear form, whose zero is the unique inactive direction. At
`p=5`, the two active pencils already determine two points. For `p>=7`,
the degree-two glue says that
`alpha_j^2-beta_j^2` is a binary quadratic; its zero at the inactive
direction makes it divisible by the degree-one difference form, so the
sums `alpha_j+beta_j` are evaluations of a second linear form. Thus the
positive and negative lines concur at two points `u,v`, joined in the
inactive square direction, and

\[
x=P(e_u-e_v).
\]

This accounts for the finite square-separated point pairs.

## The half-conic rigidity lemma and `t=2`

A `t=2` profile is `delta_alpha+delta_beta`. Degree-one glue makes
`alpha_j+beta_j` the values of a binary linear form. Degree-two glue makes
`alpha_j beta_j` the values of a binary quadratic form; for `p=5`, the
three direction values interpolate that quadratic automatically. Hence

\[
D(t_j)=(\alpha_j-\beta_j)^2                    \tag{7}
\]

for one binary quadratic `D`. The square direction forms are one
quadratic-character half of `P^1(F_p)` for an anisotropic norm form `N`.
Equation (7) says that `D` is a nonzero square on every point of that half.

The required rigidity statement is:

> If `p>=5`, `N` is an anisotropic binary quadratic, and `D` is square
> and nonzero on all points where `eta(N)=epsilon`, then either
> `D=ell^2`, with the root of `ell` in the opposite half, or `D=cN`,
> with `eta(c)=epsilon`.

To prove it, first classify `D` by its discriminant. A split nondegenerate
binary quadratic has only `(p-1)/2` positive projective values and cannot
cover a half of size `(p+1)/2`. A rank-one form must be a square `ell^2`,
and its one root must lie outside the selected half. If `D` is anisotropic,
it has exactly `(p+1)/2` positive values, so its positive set must equal the
selected half. Thus `eta(DN)` is constant on all `p+1` rational projective
points.

If `D` is not proportional to `N`, the double cover

\[
Y^2=D(X,Z)N(X,Z)
\]

is a smooth genus-one curve. Constant character `+1` would give it
`2(p+1)` rational points; constant character `-1` would give none. Both
contradict the Hasse interval

\[
p+1-2\sqrt p\le\#E(\mathbb F_p)\le p+1+2\sqrt p
\]

for every `p>=5`. Therefore `D=cN`. This proves the lemma.

The rank-one case of the lemma gives two finite points whose difference
has nonsquare direction, hence `P(e_u+e_v)`. The `cN` case is realized by
the signed complements of square circles not containing infinity. Their
counts per translation centre are respectively

\[
{p^2-1\over4},\qquad {p-1\over2},
\]

so the exhibited families exhaust every `t=2` equality case.

Combining all three values of `t` proves the shell classification and the
count `p(p+1)(p^2+1)`.

## Exact degree-four harmonic decomposition

Let `W` satisfy

\[
PWP=W,\qquad \operatorname{diag}W=0,\qquad
\operatorname{tr}W=0,
\]

and let `H_W` be the degree-four harmonic polynomial from Proposition
15.631. Put `F=||W||_F^2`. On the point-pair orbit,

\[
u=P(e_i-C_{ij}e_j),\quad
(u^TWu)^2=4W_{ij}^2.
\]

Using `CW=pW`,

\[
\sum_{i<j}(u^TWu)^2=2F,
\qquad
\sum_{i<j}u^TW^2u=p(p-1)F.
\]

After including both signs and evaluating at `u/2`, the complete pair
orbit contributes

\[
\boxed{{1\over4}\left(1-{(p-1)^2\over d+2}\right)F.}       \tag{8}
\]

For the square-circle representatives `u=w_S/p`, signed PSL irreducibility
and the trace determine

\[
\sum_S u u^T=(p-1)P.
\]

The circle contribution is therefore

\[
\boxed{
 {1\over8p^4}\sum_S(w_S^TWw_S)^2
 -{(p-1)^2\over4p(d+2)}F.}                                \tag{9}
\]

Equation (9) identifies the first anisotropic dual-shell term exactly: it
is a positive-semidefinite square-circle evaluation operator plus an
explicit scalar. This is the usable R1 advance. It does not assign a sign
to all later shells.

## Independent exact audits and searches

PARI/GP `qfminim` on the saturated dual Gram form gives the following
complete antipodal counts through scaled norm `2(p-1)`:

| p | first shell | second shell | formula for second |
|---:|---:|---:|---:|
| 3 | 10 | 15 | exceptional |
| 5 | 26 | 390 | 390 |
| 7 | 50 | 1400 | 1400 |
| 11 | 122 | 8052 | 8052 |

The CUDA sparse scan supplies exactly the point-pair orbit. Exact numerator
comparison leaves `65,175,671` vectors at `p=5,7,11`; in every case those
vectors have `p+1` zeros, all other numerators equal to `+/-2`, and their
zero sets are precisely the square circles.

Exact OEIS searches for `780,2800,16104,30940,88740,137560,292560,732540`
and for the polynomial formula returned no matching sequence. Web and
literature searches found the adjacent conference-ETF lattice papers, but
no second-shell classification or count. The bare polynomial sequence is
not asserted to be new independently of this lattice interpretation.
