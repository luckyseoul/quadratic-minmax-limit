# Actual diagonal mixed moments for the cross-rectangle row function

Status: proved exact all-orders mixed-moment identities (2026-09-05).
The proofs are analytic; finite regression checks are supplementary.
This extends the interface in
`evidence/NOTE_2026-09-03_CROSS_RECTANGLE_FOURIER_STABILITY.md` by retaining
the actual diagonal blocks. It does not prove its equation (33).

All expectations below are normalized uniform expectations. Write
`mu_4={1,i,-1,-i}`. Matrices denoted Hermitian have their usual complex
adjoints. All diagonal blocks in this note have zero diagonal.

## 1. The row function as a Boolean absolute-value average

Let `w` be uniform on `mu_4^r`, and write, independently in each coordinate,

\[
 w_j=\frac{\alpha_j+i\beta_j}{1+i},
 \qquad \alpha_j,\beta_j\in\{\pm1\}.
\]

All `2r` Boolean signs are independent. With

\[
 A_0=\sum_j\alpha_j,\qquad B_0=\sum_j\beta_j,
\]

the unphased row function is

\[
 h_r(w)=\max(|A_0|,|B_0|)
 =\frac{|A_0+B_0|+|A_0-B_0|}{2}.                 \tag{D1}
\]

Thus, if `f(alpha,beta)=|sum alpha + sum beta|`, then

\[
 h_r=\frac{f(\alpha,\beta)+f(\alpha,-\beta)}2.
\]

In particular `E h_r=eta_r=E|S_(2r)|` exactly, not just asymptotically.

For `gamma in (Z/4Z)^r`, choose each exponent in `{0,1,2,3}`, put

\[
 \chi_\gamma(w)=\prod_jw_j^{\gamma_j},
 \qquad d(\gamma)=\#\{j:\gamma_j=1\text{ or }3\}
                  +2\#\{j:\gamma_j=2\},
\]

and use the Fourier convention
`hhat(gamma)=E[h_r conjugate(chi_gamma)]`.
The number `d(gamma)` is its Lee weight.

The needed coordinate expansions are

\[
 \bar w_j=\frac{\alpha_j-i\beta_j}{1-i},
 \quad \bar w_j^2=w_j^2=\alpha_j\beta_j,
 \quad \bar w_j^3=w_j=\frac{\alpha_j+i\beta_j}{1+i}. \tag{D2}
\]

Every monomial in the expansion of `conjugate(chi_gamma)` therefore uses
exactly `d(gamma)` DISTINCT Boolean variables. This includes exponent `2`:
it uses both `alpha_j` and `beta_j`, rather than a repeated Boolean variable.

For `N=2r`, set

\[
 c_{N,d}=\mathbb E\left[|S_N|\prod_{j=1}^d\varepsilon_j\right].
\]

Exchangeability says that every degree-`d` square-free Boolean monomial has
this coefficient against `f`. The sum of the coefficients in (D2)'s full
product is its value at all `alpha=beta=1`, namely `1`. After replacing every
`beta` by `-beta`, that sum is instead

\[
 \prod_j\overline{\left(\frac{1-i}{1+i}\right)}^{\gamma_j}
 =i^{\sum_j\gamma_j}.
\]

Consequently the exact coefficient identity is

\[
 \boxed{\widehat h_r(\gamma)
 =c_{2r,d(\gamma)}\frac{1+i^{\sum_j\gamma_j}}2.}  \tag{D3}
\]

If `d` is odd, `c_(2r,d)=0` by simultaneous reversal of all Boolean signs.
For even `d`, the remaining factor is `1` when `sum gamma=0 mod4`, and `0`
when `sum gamma=2 mod4`. Thus every balanced character of a given even Lee
weight has the same real coefficient. Here “balanced” means total exponent
zero modulo four.

## 2. The exact Lee-weight-four coefficient, including r=2 and r=3

For `N=2r>=4`, conditioning on the first two signs and using the step-two
second difference of `|t|` on the even lattice gives

\[
 c_{N,2}=\mathbb P(S_{N-2}=0)=a_r.
\]

For degree four the same conditioning gives

\[
 c_{N,4}
 =\mathbb E\left[\varepsilon_3\varepsilon_4
       1_{\{\varepsilon_3+\cdots+\varepsilon_N=0\}}\right].
\]

Put `T=epsilon_5+...+epsilon_N`, including the deterministic empty sum `T=0`
when `N=4`. Writing `p_t=P(T=t)`, the four choices of signs 3 and 4 yield

\[
 c_{N,4}=\frac{p_{-2}+p_2-2p_0}{4}
         =\frac{p_2-p_0}{2}.                     \tag{D4}
\]

With `m=N-4`, the central-binomial ratios are

\[
 p_2=\frac{m}{m+2}p_0,
 \qquad a_r=\frac{m+1}{m+2}p_0.
\]

These formulas also hold for `m=0`, where `p_2=0,p_0=1`. Hence

\[
 \boxed{c_{2r,4}=-\frac{a_r}{2r-3}.}             \tag{D5}
\]

Write `d_r=a_r/(2r-3)>0`. Equations (D3)--(D5) prove that every balanced
Lee-weight-four character has coefficient `-d_r`.

The small-index cases require no separate limiting argument:

- `r=2`: `a_2=1/2`, `d_2=1/2`. The character `w_1^2 w_2^2` uses the four
  independent Boolean variables `alpha_1,beta_1,alpha_2,beta_2`; its coefficient
  is `-1/2`.
- `r=3`: `a_3=3/8`, `d_3=1/8`. Both the two-coordinate characters
  `w_i^2 w_j^2` and the three-coordinate fork characters
  `w_i^2 w_j w_k`, `w_i^2 conjugate(w_j) conjugate(w_k)` have coefficient
  `-1/8`. No four-distinct-coordinate term is needed or assumed.

## 3. Exact Fourier decomposition of a diagonal-block energy squared

Let `B` be an arbitrary zero-diagonal Hermitian matrix of order `r`, and put

\[
 b(w)=w^*Bw=\sum_{j\ne k}B_{jk}\bar w_jw_k.
\]

Squaring shows that its Fourier support has only Lee weights `0,2,4`:

1. A directed edge followed by its reverse contributes the constant
   `p_0=tr(B^2)`.
2. A directed path `j->k->m`, with `j!=m`, contributes to weight two. The
   two orders of its factors give coefficient `2(B^2)_(jm)`.
3. Repeating one directed edge gives weight four on TWO coordinates. For
   the character `w_j^2 w_k^2`, the coefficient is
   `B_(jk)^2+B_(kj)^2=2 Re(B_(jk)^2)`.
4. Two edges with a common source or common target give weight four on
   THREE coordinates. Their coefficients are respectively
   `2B_(ij)B_(ik)` for `w_i^2 w_j w_k` and its conjugate for
   `w_i^2 conjugate(w_j) conjugate(w_k)`.
5. For four distinct coordinates, the coefficient of
   `conjugate(w_i) conjugate(w_j) w_k w_l` is
   `2(B_(ik)B_(jl)+B_(il)B_(jk))`.

This is exhaustive: zero diagonal excludes loops, and two directed edges
are reverse, a cancelling path, repeated, a source/target fork, or disjoint.
In particular, case 3 must NOT be discarded as if the fourth phases were
Gaussian; it survives at `r=2`. Cases 3 and 4 exhaust the weight-four part
at `r=3`.

Let `p_2,p_4` denote the weight-two and weight-four components of `b^2`.
For every fourth-phase vector `v`,

\[
 \begin{aligned}
 p_0&=\operatorname{tr}B^2,\\
 p_2(v)&=2v^*(B^2-\operatorname{diag}(B^2))v
       =2v^*B^2v-2\operatorname{tr}B^2,\\
 p_4(v)&=b(v)^2-2v^*B^2v+\operatorname{tr}B^2.
 \end{aligned}                                      \tag{D6}
\]

The diagonal subtraction uses `|v_j|=1`; no large-order assumption occurs.

## 4. The row-function/actual-energy covariance

Let `Q in mu_4^(ell x r)`, write its rows as `q_a`, and set

\[
 F_Q(w)=\sum_{a=1}^{\ell}h_r(q_a\circ w),
 \quad v_a=\overline{q_a}^{\,T},
 \quad T=Q^*Q=\sum_av_av_a^*.
\]

Changing variables `u=q_a circ w` in one row, then using (D3)--(D6), gives

\[
 \mathbb E[h_r(q_a\circ w)b(w)^2]
 =\eta_r p_0+a_r p_2(v_a)-d_r p_4(v_a).           \tag{D7}
\]

This also fixes the conjugation convention: the diagonal energies are
evaluated at `v_a=conjugate(q_a)`, not at `q_a`.
Summing (D7), and using
`E F_Q=ell eta_r`, `E b^2=tr B^2`, proves

\[
 \boxed{\begin{aligned}
 \operatorname{Cov}(F_Q,b^2)
 ={}&2(a_r+d_r)\operatorname{tr}\bigl(B^2(T-\ell I)\bigr)\\
 &-d_r\left(\sum_a(v_a^*Bv_a)^2
                     -\ell\operatorname{tr}B^2\right).
 \end{aligned}}                                      \tag{D8}
\]

For comparison, the weight-two identity alone gives

\[
 \mathbb E[F_Q b]=a_r\operatorname{tr}(BT),          \tag{D9}
\]

because `E b=tr B=0`.
The covariance (D8), unlike (D9), detects the actual squared diagonal
energies at the conjugated row phases.

If `T=ell I`, the exact surviving term is

\[
 \operatorname{Cov}(F_Q,b^2)
 =d_r\left(\ell\operatorname{tr}B^2
                   -\sum_a(v_a^*Bv_a)^2\right).    \tag{D10}
\]

As a direct symbolic check at `r=2,ell=1,q_1=(1,1)`, take
`B=[[0,u+iv],[u-iv,0]]`. Its square is `(u^2+v^2)I`, and (D8) reduces to
`Cov(F_Q,b^2)=v^2-u^2`. This uses precisely the repeated-edge term from
case 3 above.

## 5. Joint left/right contraction, with the zero-diagonal corrections

Let

\[
 L=H_L,\quad B=H_R,\quad G=C+iR=(1+i)Q,
\]

and take independent uniform fourth-phase vectors `z,w`. Write

\[
 a=z^*Lz,\quad b=w^*Bw,\quad s=z^*Gw,\quad D=a+b.
\]

For independent fourth phases `u_i`, the exact fourth-moment rule is

\[
 \mathbb E[\bar u_i u_j\bar u_k u_l]
 =\delta_{ij}\delta_{kl}+\delta_{il}\delta_{kj}
   -1_{\{i=j=k=l\}}.                              \tag{D11}
\]

There is no `i=k,j=l` contraction when the two indices differ, because
`E u_i^2=0`. The last term in (D11) removes the double count when all four
indices coincide. Consequently, for arbitrary matrices `M,N`,

\[
 \mathbb E[(u^*Mu)(u^*Nu)]
 =\operatorname{tr}M\operatorname{tr}N+
   \operatorname{tr}(MN)-\sum_iM_{ii}N_{ii}.        \tag{D12}
\]

In particular, if `N` has zero diagonal then its trace and the last term
vanish. Applying (D12) first in `z` and then in `w` gives

\[
 \boxed{\mathbb E[|s|^2ab]
       =\operatorname{tr}(LGBG^*).}                \tag{D13}
\]

More explicitly, with `K(w)=Gww^*G^*`,

\[
 \mathbb E_z[|s|^2a]=\operatorname{tr}(K(w)L)
                   =w^*G^*LGw;
\]

the zero diagonal of `L` removes the trace and diagonal-contraction terms.
Multiplication by `b` and averaging `w` removes those same two terms using
the zero diagonal of `B`, leaving (D13).

For the two marginal squared energies, use the weight-zero/two components
of `a^2` from (D6). Against the quadratic `z^*Kz`, its weight-four part is
orthogonal. Hence, for ANY matrix `K`,

\[
 \mathbb E_z[(z^*Kz)a^2]
 =\operatorname{tr}K\operatorname{tr}L^2
   +2\operatorname{tr}((K-\operatorname{diag}K)L^2).\tag{D14}
\]

Here the subtraction of `diag K` is essential; replacing this formula by
the Gaussian expression would be incorrect. Since

\[
 \mathbb E_w K(w)=GG^*,\qquad
 \operatorname{diag}(GG^*)=2rI_\ell,
\]

(D14) yields

\[
 \mathbb E[|s|^2a^2]
 =\|G\|_F^2\operatorname{tr}L^2
  +2\operatorname{tr}(L^2(GG^*-2rI)).              \tag{D15}
\]

Likewise,

\[
 \mathbb E[|s|^2b^2]
 =\|G\|_F^2\operatorname{tr}B^2
  +2\operatorname{tr}(B^2(G^*G-2\ell I)).          \tag{D16}
\]

Finally `E|s|^2=||G||_F^2` and
`E D^2=tr L^2+tr B^2`, since the two diagonal energies are independent and
mean zero. Combining (D13)--(D16) proves

\[
 \boxed{\begin{aligned}
 \operatorname{Cov}(|s|^2,D^2)
 ={}&2\operatorname{tr}(L^2(GG^*-2rI))\\
 &+2\operatorname{tr}(B^2(G^*G-2\ell I))\\
 &+2\operatorname{tr}(LGBG^*).
 \end{aligned}}                                    \tag{D17}
\]

For an arbitrary rectangular `G`, the same proof replaces `2rI` and
`2ell I` by `diag(GG*)` and `diag(G*G)`. The displayed constants use exactly
the sign-pair entry condition `|G_(ij)|^2=2`.

When `ell=r` and `G=sqrt(2r)U` with `U` unitary, the marginal Gram defects
vanish, but the actual diagonal interaction remains:

\[
 \begin{aligned}
 \operatorname{Cov}(|s|^2,D^2)
 &=4r\operatorname{tr}(LUBU^*)\\
 &=2r\left(\|L+UBU^*\|_F^2-\|L\|_F^2-\|B\|_F^2\right).
 \end{aligned}                                      \tag{D18}
\]

Thus the extra datum is opposition of the two actual diagonal blocks in
the cross-block basis. Gram perfection alone does not eliminate that datum.

## 6. Exact remaining quantifier; no pointwise claim

No global-minimizer assumption was used in (D8) or (D17). These identities
therefore identify mixed quantities that a construction or a global-
minimality theorem could control, but do not supply that control.
Even an established negative covariance would be an averaged statement;
equation (33) is a uniform maximum statement.

For the notation of the existing note, set

\[
 \delta_z(w)=F_Q(w)-\rho(z^*Gw)\ge0.
\]

For a fixed selected globally minimal real signing `A`, the missing
assertion is the existence of ONE compatible choice of all skew blocks,
with a Dini-admissible error `e_n=o(n^(3/2))`, such that

\[
 \boxed{|a(z)+b(w)|\le
   2(\sqrt2M-F_Q(w))+2\delta_z(w)+e_n
       \quad\text{for EVERY }z,w.}                \tag{D19}
\]

Taking the maximum over `z` shows that (D19) is exactly equation (33), not
a stronger surrogate. Restricting it only to cross-maximizers (`delta=0`)
is necessary but not sufficient: near-maximizers must also satisfy their
own diagonal-payment allowance.

For the multiplier-two route it suffices, at each sufficiently large
order, to SELECT a globally minimizing `A` and one such skew completion.
Requiring (D19) for every global minimizer would be a stronger, unnecessary
quantifier. The sufficient quantifier is `for every order, there exist a
global minimizer A and a compatible S, such that every z,w satisfies (D19)`.

Conditioning Parseval on an energy stratum introduces the stratum's full
Fourier correlation matrix and provides no uniform bound on its maximizing
rows from the identities proved here. Moreover, a common fourth-phase
rotation of an entire block leaves its diagonal energy, `F_Q`, and `rho`
unchanged; that symmetry was already exhausted in equation (29).

The exact outstanding quantifier is therefore a SAME-choice-of-skew-blocks,
ALL-phase-states diagonal-payment bound. Nothing in this note claims to
derive it, or to close multiplier two or the original MathOverflow limit.

The exact additive compatibility restriction, including the quadratic/sign
realization gap, is proved separately in
`NOTE_2026-09-05_DIAGONAL_PAYMENT_COMPATIBILITY.md`.
`src/original_mo_diagonal_moments.py` and
`tests/test_original_mo_diagonal_moments.py` give exact-arithmetic reference
formulas and independent finite expectation checks, including the small
orders where repeated-edge and three-coordinate fourth-phase terms matter.
