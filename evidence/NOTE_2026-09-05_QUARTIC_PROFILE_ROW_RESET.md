# Actual whole-row and multirow reset budgets for the quartic profile

2026-09-05. **Analytic necessary conditions at an actual global optimum.**
These are finite multi-edge variations, not a sum of independently
reoptimized edge flips. They do not determine the signed internal/cross
gap integral or prove order transport. No computation is used.

## 1. Setup and the deleted Gibbs law

Use the exact balanced profile and symmetric pressure from
`NOTE_2026-09-05_QUARTIC_PENALIZED_PROFILE_IDENTITY.md`. Thus N=2n,

\[
 M=\frac{\sqrt{2-t}A_I+\sqrt t A_C}{\sqrt N},\qquad
 F(M)=\log\mathbb E_x\cosh\left(\frac c2x^TMx\right),
 \quad V(M)=\operatorname{tr}M^4,
\]

where c,lambda>0 and 0<=t<=1. The complete signing A is an ACTUAL
global minimizer of F(M)+lambda V(M) over all signings at these fixed
parameters. All rows have squared norm

\[
 d=(M^2)_{ii}=1-\frac{2-t}{N}\le1.
\]

Fix a vertex i. In block order (i,T), write

\[
 M=\begin{pmatrix}0&w^T\\w&B\end{pmatrix},
 \qquad \|w\|_2^2=d.
\]

Here B is the principal matrix on T with the ORIGINAL normalization
and weights. It is not renormalized to order N-1, and it need not
minimize any smaller-order problem. Set

\[
 F(B)=\log\mathbb E_{x_T}\cosh(cQ_B(x_T)),\qquad
 \mu_B(\sigma,x_T)=\frac{\exp(\sigma cQ_B(x_T))}
                         {2^N\exp(F(B))},
 \quad Q_B(x_T)=\tfrac12x_T^TBx_T.
\]

This is the ACTUAL augmented symmetric Gibbs law of the deleted
matrix: sigma is not replaced by a fixed phase or assigned uniform
weight after conditioning on x_T. With h(x_T)=c w^Tx_T, summation over
the removed spin gives exactly

\[
 p_i:=F(M)-F(B)=\log\mathbb E_{\mu_B}\cosh h\ge0.       \tag{1}
\]

## 2. Exact random-row comparison

Replace all signs in row i, symmetrically in column i, by independent
fair signs, preserving each magnitude |w_j| and leaving B untouched.
Let the resulting vector be w' and matrix M'. The zero cross weights
at t=0 cause no difficulty. Global optimality and Jensen's inequality
give

\[
 F(M)+\lambda V(M)
 \le\mathbb E_{w'}[F(M')+\lambda V(M')]
 \le F(B)+\sum_j\log\cosh(c|w_j|)
                         +\lambda\mathbb E_{w'}V(M').     \tag{2}
\]

Indeed, for the normalized symmetric partition function Z=exp(F),
independent sign averaging gives the exact identity

\[
 \mathbb E_{w'} Z(M')=Z(B)\prod_j\cosh(c|w_j|).
\]

The quartic trace has the exact block expansion

\[
 V(M)-V(B)=4w^TB^2w+2d^2
         =4(M^4)_{ii}-2d^2.                              \tag{3}
\]

Since (B^2)_{jj}=d-w_j^2 and the w'_j are independent with mean zero,

\[
 \mathbb E_{w'}[V(M')-V(B)]
 =4\sum_jw_j^2(d-w_j^2)+2d^2
 =6d^2-4\sum_jw_j^4.                                    \tag{4}
\]

Combining (1)--(4) proves the whole-row inequality

\[
 \boxed{\quad
 p_i+4\lambda(M^4)_{ii}
 \le\sum_j\log\cosh(c|w_j|)
                  +8\lambda d^2-4\lambda\sum_jw_j^4
 \le\frac{c^2d}{2}+8\lambda d^2-4\lambda\sum_jw_j^4.
 \quad}                                                   \tag{5}
\]

The negative quartic terms on the right of the equivalent bound for
p_i are retained. In particular, (5) recovers
`(M^4)ii <= 2d^2+c^2 d/(8 lambda)` and also controls a full finite
row reset rather than only its single-edge linearization.

Since `(M^4)ii>=d^2`, a convenient weaker consequence is

\[
 0\le p_i\le\frac{c^2d}{2}+4\lambda d^2\le\frac{c^2}{2}+4\lambda.
                                                               \tag{6}
\]

## 3. Cavity tails and the actual-law distinction

Equation (1) is a logarithmic moment bound under the DELETED Gibbs law.
For every s>=0, (5) therefore gives

\[
 \mu_B\{|h|\ge s\}
 \le 2\exp\!\left[
 \frac{c^2d}{2}+8\lambda d^2-4\lambda\sum_jw_j^4
                      -4\lambda(M^4)_{ii}-s\right].      \tag{7}
\]

The remaining-spin marginal of the full Gibbs law is instead
`d(mu_M)_(sigma,x_T)/d mu_B = exp(-p_i) cosh h`.
Consequently (7) must NOT automatically be asserted for the full
Gibbs law: its density contains the same exponential tilt that
produced the cavity moment.

One can separately bound actual relative entropy using the previously
proved row optimality inequality. On the full space (sigma,x), let
nu_i be mu_B times an independent uniform removed spin. Then exactly

\[
 D(\mu_M\Vert\nu_i)=\mathcal E_i-p_i,\qquad
 \mathcal E_i=\mathbb E_{\mu_M}[h\tanh h].                 \tag{8}
\]

The row inequality from the quartic-profile note is

\[
 \mathcal E_i+8\lambda(M^4)_{ii}
                   +8\lambda\sum_jw_j^4
 \le c^2d+16\lambda d^2.
\]

Thus the actual-law entropy obeys the uniform bound

\[
 0\le D(\mu_M\Vert\nu_i)
 \le c^2d+16\lambda d^2-8\lambda(M^4)_{ii}
                         -8\lambda\sum_jw_j^4-p_i
 \le c^2+8\lambda.                                      \tag{9}
\]

This compares actual Gibbs measures on the same state space. It does
not identify mu_B with a minimizing smaller-order Gibbs measure.

## 4. A finite reset of any set of rows

More generally, let S have k vertices, let T be its complement, and
retain B=M_T. Replace independently all signs of edges incident to S,
preserving the original edge magnitudes and B. For j in T put

\[
 b_j=\sum_{i\in S}m_{ij}^2,\qquad (B^2)_{jj}=d-b_j,
\]

and let E_S be the unordered set of edges incident to S. Write

\[
 W_S=\sum_{e\in E_S}m_e^2\le kd,
 \qquad H_S=\sum_{e\in E_S}m_e^4.
\]

Exactly the same partition-function averaging and global-optimality
argument gives

\[
 F(M)-F(B)+\lambda[V(M)-V(B)]
 \le\sum_{e\in E_S}\log\cosh(c|m_e|)
                             +\lambda\mathcal V_S,       \tag{10}
\]

where the EXPECTED quartic increment of the random refill is exactly

\[
 \boxed{\quad
 \mathcal V_S
 =2kd^2+4d\sum_{j\in T}b_j-2\sum_{j\in T}b_j^2-2H_S
 \le6kd^2.
 \quad}                                                   \tag{11}
\]

To verify (11), expand the trace over closed walks of length four.
The walks with a repeated edge contribute
`2 sum_i d_i^2 - sum_(i,j) m_ij^4`. A simple four-cycle survives sign
averaging precisely when all its vertices lie in T; those terms are
exactly the same in V(B). Subtracting the retained repeated-edge terms
gives
`2[k d^2 + sum_T(d^2-(d-b_j)^2)] - 2H_S`, which is (11).
The last inequality uses `sum_T b_j<=kd`.

Both terms on the left of (10) are nonnegative. Pressure monotonicity
follows by Jensen averaging the spins in S in each fixed phase.
Quartic monotonicity follows directly from the block identity, writing
`M=[[C,R],[R^T,B]]`:

\[
 \begin{aligned}
 V(M)-V(B)={}&\operatorname{tr}(C^2+RR^T)^2
       +2\operatorname{tr}(B^2R^TR)
       +\operatorname{tr}(R^TR)^2\\
       &+2\|CR+RB\|_F^2\ge0.
 \end{aligned}
\]

Therefore a genuine finite multirow budget is

\[
 \boxed{\quad
 0\le F(M)-F(M_T)+\lambda[\operatorname{tr}M^4-
                                      \operatorname{tr}M_T^4]
 \le\frac{c^2}{2}W_S+\lambda\mathcal V_S
 \le\left(\frac{c^2d}{2}+6\lambda d^2\right)|S|.
 \quad}                                                   \tag{12}
\]

For S={i}, (10)--(11) give exactly (5), not a different normalization.
For S empty the assertion is zero; for S the whole vertex set take
F(M_T)=V(M_T)=0. No smaller-order minimality is needed in any case.

## Scope

The conclusions hold for every active global optimizer, every profile
parameter including both endpoints, and every chosen row set. They
bound actual pressure and quartic cost jointly under deletion and
random refill. An O(|S|) cavity budget is not an o(N) transport bound
when |S| is a positive fraction of N. No sign of the mixed
internal/cross reset-gap integral follows here.
