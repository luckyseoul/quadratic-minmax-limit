# Quartic spectral regularization of the actual optimized balanced profile

2026-09-05. **All-orders regularized identity; the signed gap integral is OPEN.**
The quartic penalty permits every discrete sign variation and has an exact,
signing-independent balanced-profile remainder. This is not a proof that the
remaining integral is nonnegative or subextensive. No computation is used.

## 1. Actual pressure, profile, and uniform fourth-moment bound

Fix c>0 and lambda>0. Let N=2n and let A range over ALL complete symmetric
zero-diagonal signings. Write A_I,A_C for its internal and cross edge
matrices with respect to two blocks of n vertices. Put

\[
 a(t)=\sqrt{2-t},\quad b(t)=\sqrt t,\quad
 M_A(t)=\frac{a(t)A_I+b(t)A_C}{\sqrt N},\quad 0\le t\le1,
\]
\[
 F_A(t)=\log\mathbb E_x\cosh\left(\frac c2x^TM_A(t)x\right),
 \quad V_A(t)=\operatorname{tr}M_A(t)^4,\quad
 G_N(t)=\min_A[F_A(t)+\lambda V_A(t)].                       \tag{1}
\]

The expectation is uniform over Boolean spins. No independent-Gaussian
replacement is made. For an active minimizing signing use its ACTUAL
augmented Gibbs law on (sigma,x), proportional to
`exp(sigma c x^T M_A(t) x/2)`, and write

\[
 \Gamma_{ij}=\langle\sigma x_ix_j\rangle,\quad
 r_e=A_e\Gamma_e,\quad
 u_I=ca/\sqrt N,\quad u_C=cb/\sqrt N.
\]

Every row of the profile has the same squared norm:

\[
 d(t)=(M_A(t)^2)_{ii}=1-\frac{2-t}{N},\qquad
 \operatorname{tr}M_A(t)^2=Nd(t)\le N.                      \tag{2}
\]

The elementary random-sign quarter-net estimate supplies a complete
signing B with `||B||_op<=8 sqrt(N)` at every order N>=2. If D is the
diagonal block-sign matrix, then

\[
 M_B(t)=\frac{(a+b)B+(a-b)DBD}{2\sqrt N},\quad
 \|M_B(t)\|_{\rm op}\le8a\le8\sqrt2.
\]

Consequently `F_B(t)<=4 sqrt(2)c N` and, by (2), `V_B(t)<=128N`.
Since F_A(t)>=0, every active minimizer satisfies the uniform bound

\[
 \boxed{\quad V_A(t)\le B_{c,\lambda}N,\qquad
 B_{c,\lambda}=128+\frac{4\sqrt2 c}{\lambda}.\quad}         \tag{3}
\]

This controls the physical profile M_A(t), not the unweighted cross
signs at t=0. It supplies no unsigned covariance diffuseness claim.

## 2. Exact finite edge flip of the quartic penalty

Fix an edge e={i,j}, put m=(M_A(t))_(ij), and let
`S_e=e_i e_j^T+e_j e_i^T`. Flipping A_e replaces M by M+H with
H=-2m S_e. Cyclicity of trace gives the EXACT expansion

\[
 \begin{split}
 \operatorname{tr}(M+H)^4-\operatorname{tr}M^4
 ={}&4\operatorname{tr}(M^3H)+4\operatorname{tr}(M^2H^2)
    +2\operatorname{tr}(MHMH)\\
   &+4\operatorname{tr}(MH^3)+\operatorname{tr}H^4.
 \end{split}
\]

Use M_ii=M_jj=0, `(M^2)ii=(M^2)jj=d`, and `S_e^2=P_e`, the
projection onto the two edge coordinates. The result is

\[
 \boxed{\quad
 \Delta_e V=-16m(M^3)_{ij}+R_e,\qquad
 R_e=32d m^2-16m^4\ge0.
 \quad}                                                       \tag{4}
\]

The last inequality follows from d>=m^2. This is a finite-flip formula,
not a Hessian evaluated only at the original matrix.

The actual pressure flip is also exact:

\[
 \phi_e=F_{A^e}-F_A
       =\log[\cosh(2u_e)-r_e\sinh(2u_e)].
\]

Every flip is admissible in (1). Thus the penalized gaps satisfy

\[
 \boxed{\qquad g_e:=\phi_e+\lambda\Delta_e V\ge0.\qquad}   \tag{5}
\]

It is g_e, not necessarily phi_e, that is nonnegative. No hard-cap
boundary variation or larger-cap optimality is being assumed.

## 3. Signed variation control for the actual penalized optimizer

For u_e>0, exponentiating (5) and using `exp(-z)>=1-z` gives

\[
 r_e\le\tanh u_e+\frac{\lambda\Delta_e V}{\sinh(2u_e)}.
\]

Since `u_e=c|m|`, (4) and `sinh(2u)>=2u` imply

\[
 (r_e)_+\le\tanh u_e+
       \frac{8\lambda}{c}|(M^3)_{ij}|+
       \frac{16\lambda}{c}d|m|.                            \tag{6}
\]

Pairing one entire spin block with its negative gives
`E cosh(u I+v C)=E[cosh(u I)cosh(v C)]`. Therefore both physical
group radial derivatives are nonnegative, for EVERY signing:
`sum_(e in g) r_e>=0` separately for g=I,C. Hence
`sum_g |r_e|<=2 sum_g (r_e)_+`.

The Schatten moment inequality and (3) give

\[
 \|M^3\|_F=(\operatorname{tr}|M|^6)^{1/2}
 \le(\operatorname{tr}M^4)^{3/4}
 \le B_{c,\lambda}^{3/4}N^{3/4},
\]
\[
 \sum_{i<j}|(M^3)_{ij}|\le\frac N2\|M^3\|_F,
 \qquad \sum_{i<j}|m_{ij}|\le\frac{\sqrt2}{2}N^{3/2}.
\]

Summing (6), using d<=1 and `tanh u<=u`, yields the explicit uniform
actual signed-correlation bound

\[
 \boxed{\quad L_N(t):=\sum_e|\Gamma_e|
 \le C_0N^{3/2}+C_1N^{7/4}=:\mathcal L_N,\quad}             \tag{7}
\]
\[
 C_0=\sqrt2\left(c+\frac{16\lambda}{c}\right),\qquad
 C_1=\frac{8\lambda}{c}B_{c,\lambda}^{3/4}.
\]

At t=0 cross correlations vanish by block-flip symmetry, so (7) holds
there without division by zero. In particular the signed matrix has
`||Gamma||_F^2<=N+2 mathcal L_N=o(N^2)` at fixed c,lambda. No bound on
either unsigned phase covariance has been used or concluded.

### 3.1 Actual row optimality and the sharp signed-sum bound

The following retains negative terms that the preceding global estimate
discarded. For each vertex put

\[
 \mathcal E_i=\sum_{j\ne i}u_{ij}r_{ij}
     =\mathbb E[g_i\tanh g_i]\ge0,\qquad
 g_i=c\sum_j M_{ij}x_j.
\]

The equality follows by conditioning on sigma and all spins except x_i.
Since the second derivative of `log E exp(s tau_e)` is at most one,
the actual edge gap obeys `phi_e<=-2u_er_e+2u_e^2`. Summing (4)--(5)
over the edges incident to i therefore gives

\[
 \boxed{\quad
 \mathcal E_i+8\lambda(M^4)_{ii}
                 +8\lambda\sum_j m_{ij}^4
 \le c^2d+16\lambda d^2.
 \quad}                                                     \tag{7a}
\]

In particular, EVERY row, not only its average, satisfies

\[
 (M^4)_{ii}\le D_0:=2+\frac{c^2}{8\lambda},\qquad
 \sum_{j\ne i}g_{ij}\le2c^2d+16\lambda d^2
                         \le2c^2+16\lambda.                \tag{7b}
\]

For the second bound use `(M^4)ii>=((M^2)ii)^2=d^2` in the same
row expansion. The two uses of g are distinguished by indices:
g_i is a local field, while g_ij is the nonnegative penalized edge gap.

The source comparison in Section 1 and the extreme-state bound give

\[
 \Phi(M)\le C_\Phi N,\qquad
 C_\Phi=\frac{4\sqrt2c+128\lambda+2\log2}{c}.              \tag{7c}
\]

Indeed, `c Phi(M)-(N+1)log2<=F_A<=G_N`. We now use the elementary
real vector-rounding inequality with `kappa=log(1+sqrt(2))` and
`K_0=pi/(2kappa)`:
for unit vectors v_i,w_k,
`|sum M_ik <v_i,w_k>|<=K_0 max_(x,y signs)|x^TMy|`.
This is the proved tensor-lift lemma: odd tensor powers with weights
`sqrt(kappa^(2l+1)/(2l+1)!)`, and alternating signs in the second lift,
give unit norms and inner product `sin(kappa<v_i,w_k>)`; Gaussian
sign rounding then has expectation `(2kappa/pi)<v_i,w_k>`.
Finite Gram realization suffices. Vectors of norm at most one can
be padded in mutually orthogonal extra coordinates without changing
the cross inner products. Thus no external Grothendieck bound is needed.

Let s_i be row i of `sign(M^3)`, with zero signs chosen as +1. Then
`||s_i||_2=sqrt(N)`, while (7b) gives
`||(M^2)_(k,.)||_2<=sqrt(D_0)`. Consequently

\[
 \begin{aligned}
 \sum_{i,j}|(M^3)_{ij}|
 &=\sum_{i,k}M_{ik}\langle s_i,(M^2)_{k,.}\rangle\\
 &\le K_0\sqrt{ND_0}\max_{x,y}|x^TMy|\\
 &\le4K_0\sqrt{ND_0}\Phi(M)
 \le4K_0C_\Phi\sqrt{D_0}\,N^{3/2}.
 \end{aligned}                                              \tag{7d}
\]

The real bilinear bound in the last line is zero-diagonal cube
interpolation and polarization. The unordered off-diagonal sum is at
most half the full entry sum. Repeating (6) with (7d) therefore upgrades
(7) to

\[
 \boxed{\quad L_N(t)\le C_*N^{3/2}=:\mathcal L_N^*,\qquad
 C_*=\sqrt2(c+16\lambda/c)
       +\frac{32\lambda K_0C_\Phi\sqrt{D_0}}{c}.
 \quad}                                                     \tag{7e}
\]

This is actual signed Gibbs control under unrestricted penalized edge
optimality. The constants C_* are also uniformly bounded at fixed c
over 0<lambda<=1, since `lambda sqrt(D_0)` is bounded there.

Equations (4) and (7b) also give `|M^3_ij|<=sqrt(dD_0)`, hence
`max_e g_e=O_(c,lambda)(N^(-1/2))`. Together with the bounded gap row
sums, this is a diffuse bounded-row gap graph. It does not by itself
exclude an extensive internal/cross cut imbalance.

## 4. Exact balanced cancellation and the integrated identity

For 0<t<=1 define the penalized reset-gap imbalance

\[
 \mathcal J_N^{\lambda}(t)=
 \frac1{4t}\sum_{e\in C}g_e
 -\frac1{4(2-t)}\sum_{e\in I}g_e
 =\sum_e\frac{u_e'}{2u_e}g_e.                               \tag{8}
\]

The leading term of (4) has the exact profile chain rule

\[
 \sum_e\frac{u_e'}{2u_e}[-16m_e(M^3)_{ij}]=-V_A'(t).
\]

More importantly, the remaining term is INDEPENDENT OF THE SIGNING:

\[
 \boxed{\quad
 \sum_e\frac{u_e'}{2u_e}R_e
 =6-2t-\frac{12-6t}{N}=:\mathcal R_N(t),\qquad
 \int_0^1\mathcal R_N(t)dt=5-\frac9N.
 \quad}                                                       \tag{9}
\]

To check this directly, use `m_e^2=a^2/N` or `b^2/N`, group sizes
`k_I=n(n-1), k_C=n^2`, and `(u_e'/u_e)m_e^2=+/-1/(2N)`.
The sum of the 32d m_e^2 terms in (9) is 4d. The sum of the
-16m_e^4 terms is `2-2t-(4-2t)/N`. Substituting (2) proves (9),
including N=2. No extensive penalty remainder is discarded.

For clarity, the pressure calculation can avoid an unspecified Taylor
error. If `b_r(s)=log(cosh s+r sinh s)`, then `|b_r''''(s)|<=2`.
Expanding `phi_e=b_(r_e)(-2u_e)` gives

\[
 \phi_e=-2u_er_e+2u_e^2(1-r_e^2)
       +\frac83u_e^3r_e(1-r_e^2)+\rho_e,
 \qquad |\rho_e|\le\frac43u_e^4.                           \tag{10}
\]

Let `J_F=sum_e u_e' phi_e/(2u_e)`. Solving (10) for r_e and summing
u_e' r_e gives

\[
 F_A'=\sum_eu_eu_e'(1-r_e^2)-J_F+E_3,
\]
\[
 |E_3|\le\frac43\sum_e|u_e'|u_e^2|r_e|
                  +\frac23\sum_e|u_e'|u_e^3.
\]

Equations (4), (8), and (9) say
`J_N^lambda=J_F-lambda V_A'+lambda R_N(t)`.
Each finite branch of (1) is absolutely continuous; its finite minimum
has the active-branch derivative almost everywhere. Therefore

\[
 G_N'(t)=\frac{c^2}{4}-\mathcal J_N^\lambda(t)
                     +\lambda\mathcal R_N(t)+E_N(t),       \tag{11}
\]

where `sum_e u_eu_e'=c^2/4` exactly. Using
`|u_e'|u_e=c^2/(2N)`, `max u_e<=c sqrt(2/N)`, and (7e),

\[
 |E_N(t)|\le
 \left(\frac{c^2}{2N}+\frac{2\sqrt2c^3}{3N^{3/2}}\right)
                 \mathcal L_N^*+\frac{c^4}{6}
 =O_{c,\lambda}(\sqrt N).                                  \tag{12}
\]

These bounds are uniform in the interior. They and the finite-branch
absolute-continuity bound also prove integrability of (8) at t=0.
Integrating (11) proves

\[
 \boxed{\quad
 G_N(1)-G_N(0)=\frac{c^2}{4}+\lambda\left(5-\frac9N\right)
       -\int_0^1\mathcal J_N^\lambda(t)dt
       +O_{c,\lambda}(\sqrt N).
 \quad}                                                       \tag{13}
\]

At fixed c the error constant in (13) is uniform over 0<lambda<=1,
by the corresponding uniform bound on C_* following (7e).

The nonnegative edge gaps in (5) have MIXED weights in (8). Neither
(3), (7), nor (9) establishes a sign or small-oh bound for their
integrated imbalance. That is the remaining unproved step.

## 5. Endpoint compatibility and approximation of the original pressure

Define the uniformly weighted penalized minimum

\[
 P_N^\lambda(c)=\min_A\left[
  \log\mathbb E_x\cosh(cQ_A(x)/\sqrt N)
                +\lambda\operatorname{tr}(A/\sqrt N)^4\right].
\]

At order one the unique zero matrix has zero pressure and penalty.

Then `G_N(1)=P_N^lambda(c)`. At t=0 the matrix M_A is block diagonal,
with blocks A_L/sqrt(n), A_R/sqrt(n), so the penalty is EXACTLY
additive with the correct smaller-order normalization. Writing
`a_A=(log Z_A^+ +log Z_A^-)/2`, phase pairing gives

\[
 G_{2n}(0)=2\min_{A\in\mathcal S_n}
       [a_A(c/\sqrt n)+\lambda\operatorname{tr}(A/\sqrt n)^4]
       \le2P_n^\lambda(c).                                \tag{14}
\]

The equality uses A,-A, whose quartic penalties are equal. It does
not identify the symmetric and half-product pressure minima.

The regularized pressures approximate the original symmetric minima.
Put `P_N(c)=min_A log E cosh(cQ_A/sqrt(N))` and
`C_c=c/4+2 log(2)/c`. Independent-edge averaging gives
`P_N(c)<=c^2(N-1)/4`; the extreme-state lower bound then shows that
an actual P_N(c) minimizer satisfies `Phi(A)<=C_c N^(3/2)`.

Use the proved same-order diagonal-majorizer trimming with threshold K,
where `gamma_0=4pi/log(1+sqrt(2))`. It deletes at most
`gamma_0 C_c N/K` vertices and admits a SAME incident-edge filler with
operator norm at most 8 sqrt(N). At the fixed c, independent filler
signs have exact phase partition-function ratio expectation
`(cosh(c/sqrt(N)))^e`. Two Markov bounds with threshold log(8) have
total failure probability at most 1/4; the existing joint operator and
Boolean failures total less than 1/2. Hence one filler satisfies all
three requirements. Principal restriction lowers each phase pressure
by Jensen. The resulting signing A' therefore has

\[
 \|A'\|_{\rm op}\le(K+8)\sqrt N,\qquad
 F_c(A')\le P_N(c)+\frac{\gamma_0C_cc^2}{2K}N+\log8.
\]

Completeness supplies the useful sharper penalty estimate
`tr(A'/sqrt(N))^4 <= (K+8)^2 (N-1)`, since its second trace is N-1.
Consequently, for every N>=2 and K>0,

\[
 \boxed{\quad
 0\le\frac{P_N^\lambda(c)-P_N(c)}N
 \le\frac{\gamma_0C_cc^2}{2K}
       +\lambda(K+8)^2+\frac{\log8}{N}.
 \quad}                                                       \tag{15}
\]

Taking K=lambda^(-1/3), for 0<lambda<=1, gives a limiting normalized
error `O_c(lambda^(1/3))`. Thus the original normalized-pressure
convergence problem is equivalent to vanishing oscillation of these
regularized sequences as lambda decreases to zero. Convergence at every
fixed lambda would suffice, but is not established by (13).

## 6. A bounded fifth spectral moment

The same attained diagonal-SDP majorizer used in same-order regularization
supplies a nonnegative diagonal D with `D+/-M>=0` and
`tr D<=gamma_0 Phi(M)`, where `gamma_0=4pi/log(1+sqrt(2))`.
Let M_+,M_- be the positive and negative spectral parts. Taking PSD
traces separately against their fourth powers gives

\[
 \begin{aligned}
 \operatorname{tr}|M|^5
 &=\operatorname{tr}(M M_+^4)
                         +\operatorname{tr}((-M)M_-^4)\\
 &\le\operatorname{tr}(D M_+^4)+\operatorname{tr}(D M_-^4)\\
 &=\operatorname{tr}(D M^4)
 \le D_0\operatorname{tr}D
 \le\gamma_0C_\Phi D_0 N.
 \end{aligned}                                               \tag{16}
\]

Thus actual quartically penalized profile minima have a uniform fifth
spectral moment as well, and `||M||_op=O_(c,lambda)(N^(1/5))`.
The diagonal nature of D is essential in the penultimate inequality.
This is additional regularity, not an integrated gap sign or transport
theorem.

## Scope

The new ingredients are unrestricted actual penalized sign variations,
uniform row-fourth and fifth-moment control, an exact signing-independent
balanced quartic remainder, an O(sqrt(N)) integrated error outside the
unresolved gap, and the compatible actual pressure endpoints and
approximation. The pressure construction is the fixed-temperature
specialization of the independently reviewed same-filler profile note.
No cap-boundary optimality, Gaussian universality, unsigned covariance
diffuseness, or favorable sign of the reset imbalance is assumed.
Even a dyadic consequence of a future gap estimate would still require
the separate all-orders transport argument. The original MO limit is OPEN.
