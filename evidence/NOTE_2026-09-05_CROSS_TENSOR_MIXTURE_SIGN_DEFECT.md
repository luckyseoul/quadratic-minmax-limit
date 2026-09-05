# Cross SDP sign defects and one tensor-mixture rounding theorem

2026-09-05. Analytic draft. This combines the actual canonical spectral
Gram with an optimal SDP Gram. The result is a necessary relation between
the Boolean/SDP ratio and the canonical-primal gap. It does not evaluate
the conditional Gaussian upper or prove original convergence.

The odd-tensor mixture construction and the cap-free residual estimate
were developed and cross-checked collaboratively. All arguments below are
analytic; no signing search or numerical mathematical calculation is used.

## 1. Actual matrix, actual SDP, and canonical spectral Gram

Let B be an n by n real sign matrix, n>=1. Put

\[
 \kappa=2/\pi,\qquad c_0=\operatorname{arsinh}1
                         =\log(1+\sqrt2),
\]
\[
 \beta=\max_{x,y\in\{-1,1\}^n}|x^TBy|,
 \qquad
 \tau=\max_{\|u_i\|=\|v_j\|=1}
                         \sum_{ij}B_{ij}\langle u_i,v_j\rangle.
                                                               \tag{1}
\]

Write |B|=(B^TB)^(1/2), let sigma_j be its eigenvalues including zeros,
and define

\[
 H_B=\begin{pmatrix}0&B\\B^T&0\end{pmatrix},\quad
 Z={1\over\sqrt n}\begin{pmatrix}B\\|B|\end{pmatrix},\quad
 W={B|B|\over n},\quad
 S={\operatorname{tr}|B|^3\over n},\quad g=\tau-S.
                                                               \tag{2}
\]

Every row of Z has norm one, since BB^T and B^TB have diagonal n.
Thus W is the cross Gram of two actual unit-vector families. In
particular |W_ij|<=1 and

\[
 \langle B,W\rangle=S,\qquad
 \|W\|_F^2={\operatorname{tr}|B|^4\over n^2},\qquad
                  n^{3/2}\le S\le\tau.                       \tag{3}
\]

The first lower bound follows from sum sigma_j^2=n^2 and the power
mean inequality. The nonnegative number g is the actual gap of this
canonical feasible primal, not an assumed consequence of signing
optimality.

Fix any actual optimal diagonal SDP dual

\[
 D=\operatorname{diag}(D_r,D_c)
   =\operatorname{diag}(d_1,\ldots,d_n,e_1,\ldots,e_n),
 \qquad D-H_B\succeq0,\quad\operatorname{tr}D=2\tau.       \tag{4}
\]

Finite-dimensional strict feasibility and SDP duality give such a D.
The entries d_i,e_j are positive by the principal minor d_i e_j>=1.
Reciprocal block rescaling of (4) and optimality imply
sum d_i=sum e_j=tau. Bipartite sign conjugation gives D+H_B>=0.

## 2. Cap-free negative-mass estimate

Define the actual sign-defect mass

\[
                       N_- =\sum_{ij}(B_{ij}W_{ij})_-.
\]

Then, for every complete cross sign matrix,

\[
                         \boxed{N_-\le g/2.}               \tag{5}
\]

Here is a full residual proof, including the column orientation. Put
Q=D-H_B. Since 0<=Q<=2D,

\[
 QD^{-1}Q\preceq2Q,\qquad
 \operatorname{tr}(Z^TQZ)=2g,
 \qquad\|D^{-1/2}QZ\|_F^2\le4g.                          \tag{6}
\]

Define

\[
 R_1=D_rB-B|B|,\quad R_2=D_c|B|-|B|^2,\quad
                       R_c=BD_c-B|B|.
\]

Thus QZ=n^(-1/2)(R_1,R_2)^T. An orthogonal polar factor U with
B=U|B| exists even if B is singular, and R_c=U R_2^T. Consequently

\[
 \|D_r^{-1/2}R_1\|_F^2+
 \|R_cD_c^{-1/2}\|_F^2
 =\|D_r^{-1/2}R_1\|_F^2+
 \|D_c^{-1/2}R_2\|_F^2\le4ng.                            \tag{7}
\]

For every real z and a>0, z_-<=(z-a)^2/(4a). Apply this first with
z=B_ij W_ij, a=d_i/n and then with a=e_j/n. The identities

\[
 B_{ij}W_{ij}={d_i\over n}-{B_{ij}(R_1)_{ij}\over n}
             ={e_j\over n}-{B_{ij}(R_c)_{ij}\over n}
\]

give separately

\[
 N_-\le{\|D_r^{-1/2}R_1\|_F^2\over4n},\qquad
 N_-\le{\|R_cD_c^{-1/2}\|_F^2\over4n}.
\]

Adding these two inequalities and using (7) proves (5). No operator
bound, maximum-diagonal bound, or small-gap hypothesis was used.

There is a stronger special case, which is not presumed in general:
if an optimal D is scalar, D=dI_(2n), then tau=nd, d>=||B||op and

\[
                         \boxed{N_-\le g/4.}               \tag{8}
\]

Indeed let E=W-(d/n)B. Then

\[
 \|E\|_F^2={1\over n^2}\sum_j\sigma_j^2(d-\sigma_j)^2
 \le {d\over n^2}\sum_j\sigma_j^2(d-\sigma_j)
 ={d\over n}g.
\]

The same scalar inequality with a=d/n gives
N_-<=n||E||_F^2/(4d)<=g/4. Scalar dual optimality does not force g=0.

## 3. One exact finite-order tensor-mixture bound

Fix a mixing parameter 0<=t<1, and set

\[
 a=\operatorname{arsinh}(1-t),\qquad r=\sin a+t<1,
 \qquad M_t={r\over(1-r^2)^{3/2}},
 \qquad\lambda(t)={1+\sec a\over2}.
                                                               \tag{9}
\]

The strict inequality for r follows from sin a<a<1-t for t<1.
Then every B satisfies

\[
 \boxed{\displaystyle
 \beta\ge\kappa\big[(a+t)\tau-t\lambda(t)g\big]
       -{\kappa t^2 M_t\over2}
                         {\operatorname{tr}|B|^4\over n^2}.}
                                                               \tag{10}
\]

If an optimal diagonal is scalar, the stronger version replaces
lambda(t) by

\[
                         \lambda_s(t)={3+\sec a\over4}.   \tag{11}
\]

To prove (10), choose optimal SDP unit vectors u_i,v_j, with cross
Gram V_ij=<u_i,v_j> and <B,V>=tau. Lift them into odd tensor powers,
with coefficients sqrt(a^(2k+1)/(2k+1)!), and put the sign (-1)^k
on the second family. The squared norm of either lift is sinh a=1-t;
the lifted cross Gram is sin(aV_ij).

Take the orthogonal direct sum of those lifts with sqrt(t) times the
two canonical unit families in Z. The resulting families are unit and
have cross Gram

\[
                         T_{ij}=\sin(aV_{ij})+tW_{ij}.
                                                               \tag{12}
\]

Although tensor notation uses a countable Hilbert direct sum, the Gram
of these 2n vectors has a finite-dimensional realization. Ordinary
finite-dimensional Gaussian sign rounding therefore gives

\[
                 \beta\ge\kappa\sum_{ij}B_{ij}\arcsin T_{ij}.
                                                               \tag{13}
\]

Every point between sin(aV_ij) and T_ij lies in [-r,r]. The second
derivative of arcsin has absolute value at most M_t there. Taylor's
theorem and arcsin(sin(aV_ij))=aV_ij yield

\[
 \arcsin T_{ij}=aV_{ij}+tW_{ij}\sec(aV_{ij})+E_{ij},
 \qquad |E_{ij}|\le{M_t t^2\over2}W_{ij}^2.              \tag{14}
\]

Since 1<=sec(aV_ij)<=sec a,

\[
 \sum_{ij}B_{ij}W_{ij}\sec(aV_{ij})
               \ge S-(\sec a-1)N_-.
\]

Equations (13)--(14), followed by (5), give

\[
 \beta\ge\kappa[a\tau+tS-t(\sec a-1)N_-]
              -{\kappa t^2M_t\over2}\|W\|_F^2,
\]

which is (10), since S=tau-g. Using (8) gives (11).

At t=0 the error vanishes and (10) is exactly

\[
                       \beta\ge\kappa c_0\tau,
 \qquad \tau\le K_G\beta,\quad K_G=(\kappa c_0)^{-1}.     \tag{15}
\]

No claim is made by substituting t=1 into (10): M_t diverges as t
approaches 1. Limits in the next section keep t fixed first.

## 4. Uniform error under only a Boolean cap

For every real entry-bounded matrix T, interpolation gives

\[
             \|T\|_{\rm op}^2\le2\beta_{\mathbb R}(T),
 \quad \beta_{\mathbb R}(T)=
       \max_{\|x\|_\infty,\|y\|_\infty\le1}|x^TTy|.     \tag{16}
\]

For completeness, its complex infinity-to-one norm is at most twice
the real cube norm: rotate u^*Tv to be real, write u=p+iq,v=r+is,
and bound its real part p^TTr+q^TTs by twice beta_R(T). Each of the
four real vectors belongs to the real unit cube. The complex
one-to-infinity norm is at most one. Complex Riesz--Thorin
interpolation halfway gives (16); the Euclidean operator norm of a
real matrix is the same over real or complex vectors.

For B, the real cube norm equals beta by multilinearity. With
L=||B||op, (3) and (16) therefore give the exact estimate

\[
 {\operatorname{tr}|B|^4\over n^2\tau}
 \le{L\operatorname{tr}|B|^3\over n^2\tau}
 \le{L\over n}\le {\sqrt{2\beta}\over n}.               \tag{17}
\]

In particular, if beta<=C n^(3/2), then the right side of (17) is
at most sqrt(2C)n^(-1/4). Thus the error in (10), divided by kappa
tau, is at most

\[
                \epsilon_{n,t}={t^2M_t\over2}
                                      \sqrt{2C}\,n^{-1/4}.   \tag{18}
\]

Since (15) gives tau<=K_G C n^(3/2), its absolute error is
O_(C,t)(n^(5/4)). Only the Boolean cap is assumed; no additional
operator cap or original-optimizer structure is hidden here.

## 5. Normalized necessary curve and its endpoint meaning

Set

\[
 u={\beta\over\tau},\qquad U={u\over\kappa},\qquad
 \rho={S\over\tau},\qquad \delta={g\over\tau}=1-\rho.
\]

For every fixed 0<t<1, (10) and (18) imply

\[
 U\ge a+t-t\lambda(t)\delta-\epsilon_{n,t},\qquad
 \delta\ge{a+t-U-\epsilon_{n,t}\over t\lambda(t)}.       \tag{19}
\]

Consider any sequence n->infinity with beta<=C n^(3/2), and any
subsequence on which U and rho converge to U_* and rho_*.
Then c_0<=U_*<=1/kappa and

\[
 \boxed{\displaystyle
 1-\rho_*\ge G(U_*),\qquad
 G(U)=\max\left(0,\sup_{0<t<1}
       {\operatorname{arsinh}(1-t)+t-U\over
          t[1+\sec(\operatorname{arsinh}(1-t))]/2}\right).}
                                                               \tag{20}
\]

If each matrix has a scalar optimal diagonal, replace the denominator
coefficient by lambda_s(t), giving a stronger curve G_s(U).
Equation (20) follows for each fixed t before taking the supremum;
it does not require a uniform Taylor estimate as t approaches 1.

These explicit one-dimensional envelopes genuinely strengthen the
cubic-only limiting constraint rho_*<=U_* whenever c_0<=U_*<1.
To see this analytically, reparameterize t=1-sinh a for 0<a<c_0.
The quotient in (20) becomes

\[
 F_U(a)={a+1-\sinh a-U\over
              (1-\sinh a)(1+\sec a)/2}.
\]

Its continuous endpoint value is F_U(0)=1-U and its right derivative
there is F_U'(0)=1-U. Thus G(U)>1-U for every c_0<=U<1.
The same derivative statement holds with lambda_s. For U>=1,
a+t<1 makes every numerator negative, so both curves equal zero.

At the former tensor-rounding endpoint u_*=1/K_G, equivalently
U_*=c_0, taking t down to zero in (20) also gives the exact bounds

\[
 1-\rho_*\ge {2(1-1/\sqrt2)\over1+\sec c_0},\qquad
 1-\rho_*\ge {4(1-1/\sqrt2)\over3+\sec c_0}
                    \quad\hbox{if the optimal diagonal is scalar}.
                                                               \tag{21}
\]

Here a=c_0-t/sqrt(2)+O(t^2), and this final t limit is taken only
after the n limit. Independently of how much (21) improves the
envelope, the strict inequality G(c_0)>1-c_0 proved above excludes
the old joint relaxation profile U_*=rho_*=c_0 for actual B.

## 6. Scope for conditional shells

The ratio u in (19)--(21) is beta(B)/tau(B). On an active cross
shell with |c|=beta(B), it also equals |c|/tau(B). If |c|<beta(B),
one must not replace beta by |c| in these constraints.

For an actual conditional cross minimizer with beta(B)<=F_A^* and
F_A^*=O(n^(3/2)), the Boolean cap required above is available. All
spectral and SDP quantities still belong to its actual B. No claim
is made that an optimal dual is scalar, that its canonical gap is
small, or that every abstract moment profile satisfying (20) is
attainable by a complete signing.

The theorem supplies an additional necessary constraint for the
Boolean-sensitive shell resolvent. The sharp evaluation over actual
nonflat duals, joint shells, and the full original covariance remains
unresolved. No conclusion about the original MO limit follows yet.
