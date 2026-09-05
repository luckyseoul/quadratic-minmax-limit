# Structural rigidity of leading half-product near-minimizers

2026-09-05. **All-orders strengthening of the exact-minimum results.**
The operator-norm and sparse-energy conclusions remain true with a
vanishing normalized optimality gap. In particular, they do not
distinguish exact minima from leading half-product near-minima.

## 1. Statement and approximate deletion budget

Fix c>0. Let A_N be complete signed hosts with beta_N=c/sqrt(N) and

\[
 0\le a_{A_N}(\beta_N)-R_N(\beta_N)=\delta_N N,
 \qquad \delta_N\longrightarrow0.
\]

Then

\[
 \boxed{\quad \|A_N\|_{\rm op}=o(N^{3/4}).\quad}              \tag{1}
\]

Also, for EVERY sequence k_N=o(N), uniformly over vertex sets S of
size at most k_N and all spin configurations x, with T=[N] minus S,

\[
 \boxed{\quad
 \max_{|S|\le k_N}\max_x
 |Q_{A_N}(x)-Q_{(A_N)_T}(x_T)|=o(N^{3/2}).
 \quad}                                                       \tag{2}
\]

All definitions are those of the exact-minimum spectral note. Repeating
its edge-randomized completion comparison, now with the allowed gap,
gives

\[
 0\le a_A(\beta)-a_{A_T}(\beta)
       \le\frac{c^2|S|}{2}+\delta_N N.                       \tag{3}
\]

The universal block-flip lower bound `a_A>=a_(A_S)+a_(A_T)` is unchanged.
With C=c+(log 2)/c, the same extreme-configuration argument gives

\[
 \boxed{\quad
 \Phi(A_S)\le C\sqrt N\,|S|+
                          \frac{2\delta_N}{c}N^{3/2}
            \quad(S\subseteq[N]).\quad}                     \tag{4}
\]

In particular Phi(A)=O_c(N^(3/2)). The reviewed elementary norm
comparison `||M||_op^2<=16 Phi(M)` for complete signings M (Section 4
of `evidence/NOTE_2026-09-05_GAUSSIAN_SIGN_INFORMATION_SCALE.md`)
therefore gives a constant L_0 depending only on c with

\[
 \|A_N\|_{\rm op}\le L_0N^{3/4}                             \tag{5}
\]

for all sufficiently large N. The same comparison applies to every
principal signing A_S. It follows by polarization bounding the real
bilinear cube norm by 4 Phi(M), decomposing complex vectors into real
and imaginary parts, and interpolating the complex infinity-to-one
norm against the one-to-infinity norm, which is one. Order zero or
one principal restrictions have zero operator norm and cause no issue.

## 2. Qualitative delocalization without fractional rounding errors

Suppose (1) fails along a subsequence, with unit eigenvectors v_N and
eigenvalues satisfying `|lambda_N|>=epsilon N^(3/4)` for a fixed
epsilon>0. We claim that ||v_N||_1/sqrt(N) is bounded below by a
positive constant on that subsequence, after discarding finitely many
orders.

If not, pass to a further subsequence on which
`d_N=||v_N||_1/sqrt(N)` tends to zero. Truncate v_N to

\[
 S_N=\{i:|v_{N,i}|\ge\sqrt{d_N}/\sqrt N\}.
\]

Then `|S_N|<=sqrt(d_N)N=o(N)` and

\[
 \|v_N-(v_N)_{S_N}\|_2^2
 \le\frac{\sqrt{d_N}}{\sqrt N}\|v_N\|_1=d_N^{3/2}=o(1).
\]

By (5), deleting this small Euclidean tail changes the Rayleigh
quadratic form by o(N^(3/4)). Thus
`||(A_N)_(S_N)||_op>=epsilon N^(3/4)-o(N^(3/4))`.
But (4) makes Phi((A_N)_(S_N))=o(N^(3/2)); the same norm comparison
then makes its operator norm o(N^(3/4)), a contradiction.

We have proved the asserted positive lower bound on ||v_N||_1/sqrt(N).
The elementary first/second-moment argument of Section 3 of the exact
spectral note now supplies fixed delta>0 and eta>0 and sets I_N with
|I_N|>=eta N on which

\[
 \delta/2\le\sqrt N\,|v_{N,i}|\le4/\delta.                   \tag{6}
\]

This avoids applying the exact note's fractional hereditary rounding
formula with a nonuniform additive near-minimizer error.

## 3. The actual sparse-pinning contradiction survives the gap

Use the exact spectral proof's pinning construction with
`r=(epsilon/2)N^(1/4)` and independent ternary signs Z_i of expectation
r v_i. The row norm bound makes all pin probabilities at most one-half.
Equations (9)--(11) there give a realization with support size
O(N^(3/4)) and squared field error O_c(N^(3/4)). Here those estimates
use only the complete-sign column norms and the unit eigenvector,
and so are unchanged.

The lower eigenvalue bound, (5), and (6) then give a fixed positive
density of complement field coordinates in a fixed interval [a,H],
with a>0. Their complement has q=N-o(N) vertices, and (4) gives

\[
 \Phi(\pm\beta A_T)\le cCq+2\delta_N N\le(cC+1)q
\]

eventually. The reviewed ACTUAL-Gibbs moderate-field theorem therefore
supplies a positive linear log-MGF response in both phases. Opposite
pinning still cancels the internal pinned energy, giving
`a_A-a_(A_T)>=Kq-O(N^(3/4))` for a fixed K>0. The approximate
upper budget (3) is `O_c(N^(3/4))+delta_N N=o(N)`.
This contradiction proves (1).

## 4. Uniform sparse-energy robustness also survives the gap

This part does not require the spectral conclusion. Fix k_N=o(N),
any |S|<=k_N, and signs z on S. Set w=beta A_(T,S)z and b=cC.
The opposite-pinning inequality is now

\[
 \frac12\{\Psi_{\beta A_T}(w)+\Psi_{-\beta A_T}(w)\}
 \le\left(\frac{c^2}{2}+\log2\right)|S|+\delta_N N=o(N),       \tag{7}
\]

uniformly in those choices. For fixed H>b put J={i in T: |w_i|>H}.
Choose the signs on J to realize the absolute cross-field sum and
apply (4) to S union J. Reversing the S spin block gives

\[
 H|J|\le\sum_{i\in J}|w_i|
              \le b(|S|+|J|)+2\delta_N N.                   \tag{8}
\]

Hence |J|=o(N) and the sum in (8) is o(N), uniformly. For any fixed
0<a<=H, a positive density of field coordinates with magnitude at
least a would, after removing J, leave a positive density in [a,H].
The actual moderate-field response with the norm cap cC+1 contradicts
(7). Therefore the number of such coordinates is o(N) uniformly.

Split the field sum below a, between a and H, and above H. Equations
(7)--(8) and the last conclusion give `||w||_1<=aN+o(N)`. First
take the limsup divided by N and then let a decrease to zero. Thus
`||w||_1=o(N)` uniformly over S,z. The cross energy is consequently
o(N^(3/2)), and (4) bounds the internal S energy by o(N^(3/2)).
This proves (2).

## Scope

These statements concern leading near-minimality of the half-product
pressure at each fixed c>0. They do not assert leading near-minimality
of Phi, O(sqrt(N)) spectral flatness, an o(k_N sqrt(N)) deletion rate,
or a fixed-positive-fraction order comparison. They are compatible
with the sparse-module full-strength canonical failure family, whose
module eigenvalues are themselves o(N^(3/4)). The original convergence
problem remains open.
