# Two-sided quartic-force bounds and the actual weighted Gibbs kernel

2026-09-05. **Proved matrix and actual-optimizer inequalities.**
The first inequality is algebraic. The second uses the actual augmented
Gibbs law and unrestricted penalized edge optimality. No unsigned
covariance substitution, computation, or order-transport assertion is used.

## 1. A two-sided diagonal majorizer for the quartic force

Let M be any real symmetric matrix. Define the entrywise product and
the diagonal matrix

\[
 H_4=M\circ M^3,\qquad D_4=\operatorname{diag}((M^4)_{ii}).
\]

Then

\[
 \boxed{-D_4\preceq H_4\preceq D_4.}                       \tag{1}
\]

Here the middle product is entrywise, whereas M^3 and M^4 are ordinary
matrix powers. In particular,

\[
 \|H_4\|_{\rm op}\le\max_i(M^4)_{ii}.                     \tag{2}
\]

Proof. For any real vector z put Z=diag(z). Direct expansion gives

\[
 z^TD_4z=\operatorname{tr}(Z^2M^4),\qquad
 z^TH_4z=\operatorname{tr}(ZMZM^3).
\]

Diagonalize M orthogonally, with eigenvalues alpha_a, and write B for
Z in this eigenbasis. Since B is real symmetric, the preceding two
forms are respectively

\[
 \frac12\sum_{a,b}(\alpha_a^4+\alpha_b^4)B_{ab}^2,
 \qquad
 \frac12\sum_{a,b}\alpha_a\alpha_b
                    (\alpha_a^2+\alpha_b^2)B_{ab}^2.
\]

Both their sum and their difference are nonnegative, because for real
p,q,

\[
 p^4+q^4-pq(p^2+q^2)=(p-q)^2(p^2+pq+q^2)\ge0,
\]
\[
 p^4+q^4+pq(p^2+q^2)=(p+q)^2(p^2-pq+q^2)\ge0.
\]

This proves (1), and D_4 is bounded above by its largest diagonal
entry times the identity, proving (2).

## 2. Actual penalized balanced-profile optimizer

Let N=2n, c,lambda>0, and 0<=t<=1. For a complete symmetric
zero-diagonal signing A define

\[
 M=\frac{\sqrt{2-t}A_I+\sqrt t A_C}{\sqrt N},\qquad
 F(M)=\log\mathbb E_x\cosh(cQ_M(x)),\qquad
 Q_M(x)=\tfrac12x^TMx.
\]

Assume A actually minimizes `F(M)+lambda tr(M^4)` over all signings
at these fixed parameters. Every row has the same squared norm

\[
 d=(M^2)_{ii}=1-\frac{2-t}{N}\le1.
\]

Use its actual Gibbs law on (sigma,x), with uniform reference measure
and density proportional to `exp(sigma cQ_M(x))`. Write

\[
 \Gamma_{ij}=\langle\sigma x_ix_j\rangle,\qquad
 H_\Gamma=M\circ\Gamma.
\]

The following bound is uniform in N and t:

\[
 \boxed{\quad\|H_\Gamma\|_{\rm op}
                 \le3c+\frac{40\lambda}{c}.\quad}         \tag{3}
\]

In particular it is uniform for 0<lambda<=1 at fixed c. It is a bound
on the weighted SIGNED kernel, not on either unsigned phase covariance.

## 3. Exact edge gaps and actual row budgets

For an unordered edge e={i,j}, put m=M_ij, u=c|m|, and
`r=A_ij Gamma_ij`. Flipping A_ij has pressure increment

\[
 \phi_e=\log[\cosh(2u)-r\sinh(2u)].
\]

Let

\[
 q_e=\phi_e+2ur.
\]

This is the logarithmic MGF of the centered variable
`-2u(tau-r)`, where `tau=A_ij sigma x_i x_j` takes values +/-1.
Jensen and the variance bound for a two-point variable give EXACTLY

\[
 0\le q_e\le2u^2.                                        \tag{4}
\]

No Taylor truncation is involved. Direct expansion with the edge
perturbation `-2m(e_i e_j^T+e_j e_i^T)` gives the finite quartic change

\[
 \Delta_eV=-16m(M^3)_{ij}+R_e,\qquad
 R_e=32dm^2-16m^4\ge0.                                   \tag{5}
\]

The actual penalized gap `g_e=phi_e+lambda Delta_eV` is nonnegative
by the assumed global sign optimality. Define

\[
 \mathcal E_i=\sum_{j\ne i}cM_{ij}\Gamma_{ij}
     =\mathbb E[h_i\tanh h_i]\ge0,\qquad
 h_i=c\sum_j M_{ij}x_j.
\]

The equality follows by conditioning on sigma and all spins except
x_i. Summing (4)--(5) over the incident edges gives

\[
 0\le\sum_{j\ne i}g_{ij}
 \le-2\mathcal E_i+2c^2d
     -16\lambda(M^4)_{ii}+32\lambda d^2
                              -16\lambda\sum_jm_{ij}^4.
                                                               \tag{6}
\]

In particular,

\[
 (M^4)_{ii}\le D_0:=2+\frac{c^2}{8\lambda}.                \tag{7}
\]

Also `(M^4)ii=sum_j ((M^2)ij)^2>=d^2`. Substituting this and
`E_i>=0` in (6) proves the sharper nonnegative gap-row budget

\[
 \sum_{j\ne i}g_{ij}\le2c^2+16\lambda.                    \tag{8}
\]

Zero-weight edges at t=0 have zero increments and require no division.

## 4. The exact force-kernel decomposition

Make q,R,g into symmetric matrices with zero diagonals and put
`B=q+lambda R`. Equations (4)--(5) give the EXACT matrix identity

\[
 2cH_\Gamma=B-g-16\lambda H_4.                            \tag{9}
\]

The matrices B and g are entrywise nonnegative. Their row budgets are

\[
 \sum_j B_{ij}\le2c^2d+32\lambda d^2\le2c^2+32\lambda,
 \qquad \sum_jg_{ij}\le2c^2+16\lambda.
\]

For a symmetric entrywise-nonnegative matrix the operator norm is at
most the maximum row sum. By (1)--(2) and (7), also
`||H_4||_op<=D_0`. Applying the triangle inequality to (9) gives

\[
 \begin{aligned}
 2c\|H_\Gamma\|_{\rm op}
 &\le(2c^2+32\lambda)+(2c^2+16\lambda)+16\lambda D_0\\
 &=6c^2+80\lambda,
 \end{aligned}
\]

which is (3). The proof concerns the actual optimized signing and its
actual signed Gibbs correlations. These operator bounds alone do not
assert that a specified balanced cut has discrepancy o(N).
