# Source-preserving whole-edge rounding and a universal simple Gaussian model

2026-09-05.  Analytic theorem; independent complete proof reads passed.
This note proves a distributional reduction on the SAME original order.
It does not infer a doubled-order upper bound from source optimality.

## 1. Symmetric edge compression, normalization, and order two

Let `K` be any complete symmetric zero-diagonal signing of order `n>=2`,
and let `L=||K||op`.  On the n^2 ordered matrix coordinates define

\[
 R_0=I-K\otimes K/L^2.
\]

The tensor eigenvalues lie in `[-L^2,L^2]`, so
`0<=R_0<=2I`, and its diagonal is one.  For each unordered edge
`e={i,j}`, `i<j`, let

\[
 F_e=(E_{ij}+E_{ji})/\sqrt2.
\]

These are orthonormal for the Frobenius inner product.  Let
`T` denote the compression of `K tensor K` to their span.  Thus

\[
 T_{ij,kl}=K_{ik}K_{jl}+K_{il}K_{jk},\qquad T_{e,e}=1.
                                                               \tag{1}
\]

The compressed covariance is `I-T/L^2` and its diagonal is
`1-1/L^2`, not one.  Put `D=L^2-1`.  For `n>=3`, its correctly
normalized correlation matrix is

\[
 \boxed{\displaystyle R={L^2I-T\over D}.}                  \tag{2}
\]

To check the denominator and uniform bound, `tr(K^2)=n(n-1)` gives
`L^2>=n-1`.  For n=3, every triangle is switching-equivalent to
`+(J-I)` or `-(J-I)`, so L=2.  Consequently `L^2>=3` for every
n>=3.  Thus D>=2 and

\[
 \operatorname{diag}R=1,\qquad
 0\preceq R\preceq{2L^2\over L^2-1}I\preceq3I.            \tag{3}
\]

There is a genuine degeneracy at n=2: L=1 and the sole symmetric
off-diagonal coordinate has compressed variance zero.  Formula (2)
is undefined and cannot be normalized.  At that order ONLY, use an
independent standard Gaussian on the sole unordered edge, so R=[1].

If `G~N(0,R)` on unordered edges, its associated symmetric matrix
has entries G_e above and below the diagonal, and zero diagonal.
There is no further factor sqrt(2) in these actual matrix entries.

## 2. The model and the proved norm comparison

Fix any real h, deterministically before the disorder.  On unordered
edges define

\[
 B_{h,e}=\operatorname{sign}(G_e+hK_e),\quad
 s=2\Phi_{\rm Gauss}(h)-1,\quad
 k=4\phi(h)^2,\quad v=1-s^2-k\ge0.
\]

Copy these entries symmetrically and set the diagonal to zero.
Thus B_h is an actual complete signing of order n and `E B_h=sK`.
Let W be an independent Gaussian matrix with independent standard
UNORDERED edge entries, copied symmetrically, and zero diagonal.
The simple mean Gaussian model is

\[
 \boxed{\displaystyle Z_h=sK+\sqrt{k}\,G+\sqrt v\,W.}      \tag{4}
\]

For ANY fixed deterministic real internal energy I_n(x) and
`|theta|<=1`, write

\[
 \mathcal M_I(H)=\max_{x\in\{-1,1\}^n}
       \left|I_n(x)+\theta\sum_{i<j}H_{ij}x_ix_j\right|.
\]

There is an absolute constant C, independent of K, h, I_n and theta,
such that

\[
 \boxed{\displaystyle
 |\mathbb E\mathcal M_I(B_h)-\mathbb E\mathcal M_I(Z_h)|
                 \le C n^{16/11}.}                        \tag{5}
\]

In particular, with I_n=0 and theta=1, this compares the ORIGINAL
whole-order quadratic norms `Phi(B_h)` and `Phi(Z_h)`.  No assumption
on either `Phi(K)` or `||K||/sqrt(n)` is needed.  The proof below
retains the nonvanishing even-Hermite term explicitly.

## 3. Exact edge correlations and the even correction

Assume n>=3 for now.  Let m=n(n-1)/2.  For distinct edges, (2) gives

\[
 R_{ij,kl}=-{K_{ik}K_{jl}+K_{il}K_{jk}\over D}.             \tag{6}
\]

Adjacent edges have correlations of magnitude 1/D; disjoint edges
have correlations in `{0,+2/D,-2/D}`.  All these arguments lie in
[-1,1], since D>=2.

Expand the centered scalar activation
`f_h(t)=sign(t+h)-s` in normalized Hermites, with coefficients c_j.
The already proved shifted-threshold identities give

\[
 c_1^2=k,\quad\sum_{j\ge2}c_j^2=v,\quad
 c_j(-h)=(-1)^{j+1}c_j(h).
\]

Define

\[
 o(t)=\sum_{\substack{j\ge3\\j\ {
 \rm odd}}}c_j^2t^j,\qquad
 e(t)=\sum_{\substack{j\ge2\\j\ {
 \rm even}}}c_j^2t^j.
\]

For |t|<=1, `|o(t)|<=|t|^3` and `0<=e(t)<=t^2`.
Set `d_e=K_e`, let `D_d=diag(d)`, and let `C_h` be the exact
centered covariance of B_h.  Its base covariance is `C_0=kR+vI`.
For distinct edges, the covariance remainder is

\[
 (C_h-C_0)_{e,f}=o(R_{e,f})+d_e d_f e(R_{e,f}),
\]

and its diagonal is zero.  The diagonal Hermite mass is already in vI.

Let E be the adjacency matrix of the line graph of the complete
graph: E_{e,f}=1 exactly when distinct edges share a vertex.  It
has constant row sum and operator norm `2(n-2)`.  Define the
four-cycle matrix on unordered edges by

\[
 Q_{ij,kl}=K_{ik}K_{jk}K_{il}K_{jl}.                        \tag{7}
\]

It is zero whenever the two edges share a vertex, including on
the diagonal.  For disjoint edges, the two summands in the numerator
of (6) have the same sign exactly when Q_{e,f}=1.  Therefore, putting

\[
 u=e(1/D),\qquad w={1\over2}e(2/D),
\]

the ENTIRE even-Hermite correction, not just its degree-two term,
is exactly

\[
 \boxed{\displaystyle
 R_{\rm even}
   =wdd^T+D_d\{wQ+(u-w)E-wI\}D_d.}                       \tag{8}
\]

Indeed, on a disjoint pair it equals `d_e d_f w(1+Q_{e,f})`;
on an adjacent pair it equals `d_e d_f u`; its diagonal is zero.
The first term in (8) is positive semidefinite of rank at most one.
Also `u<=1/D^2`, `w<=2/D^2`, and `|u-w|<=3/D^2`.

## 4. A universal operator bound for the four-cycle matrix

Introduce the full ordered-pair matrix

\[
 F_{ij,kl}=K_{ik}K_{jk}K_{il}K_{jl},
\]

including pairs i=j.  It is real symmetric.  Direct multiplication
gives the important nonnegative-entry identity

\[
 (F^2)_{ij,pq}
       =\left(\sum_r K_{ir}K_{jr}K_{pr}K_{qr}\right)^2.
                                                               \tag{9}
\]

Let V have rows indexed by ordered pairs and columns by vertices,
`V_{ij,r}=K_{ir}K_{jr}`.  Then

\[
 V^TV=(K^2)\circ(K^2)\preceq L^2(n-1)I.                  \tag{10}
\]

For (10), apply the positive Schur map for the PSD matrix K^2 to
`K^2<=L^2I`, and use `diag(K^2)=(n-1)1`.  Each row of V has squared
norm at most n-1.  The row sum in (9) is consequently

\[
 \sum_{p,q}(F^2)_{ij,pq}
       =V_{ij,\cdot}(V^TV)V_{ij,\cdot}^T
       \le L^2(n-1)^2.
\]

Since F^2 has nonnegative entries, its spectral radius is at most
its largest row sum.  Thus `||F||<=L(n-1)`.  Compression of F to
the symmetric off-diagonal orthonormal basis has entries 2Q_{e,f}.
Therefore

\[
 \boxed{\displaystyle \|Q\|_{\rm op}\le{L(n-1)\over2}.} \tag{11}
\]

This does not assume that K has bounded normalized operator norm.

## 5. Remainder size and its actual expected norm cost

There are 2(n-2) adjacent edges and (n-2)(n-3)/2 disjoint edges
in each covariance row.  Hence the odd-Hermite remainder satisfies

\[
 \|R_{\rm odd}\|_{\rm op}
        \le{2(n-2)(2n-5)\over D^3}.                       \tag{12}
\]

Combining (8), (11), and (12) gives

\[
 C_h=(C_0+wdd^T)+\mathcal E,
\]
\[
 \boxed{\displaystyle
 \|\mathcal E\|_{\rm op}
 \le {2(n-2)(2n-5)\over D^3}
       +{L(n-1)+6(n-2)+2\over D^2}
 \le {40\over\sqrt n}.}                                  \tag{13}
\]

For the final absolute bound when n>=4, use `D>=n/2` and
`D>=L^2/2`.  The odd term is at most 32/n; the last two terms
in the second numerator cost at most 24/n; and its first term
is at most `4/sqrt(n-1)<=5/sqrt(n)`.  Thus their sum is at most
`33/sqrt(n)` for n>=4.  For n=3, L=2 and D=3 verify (13)
directly.  Both C_h and C_0+wdd^T are PSD; the error need not have
a sign.

The rank-one Gaussian retained here is simply `sqrt(w) xi K`,
where xi is an independent scalar standard Gaussian.  Its expected
whole quadratic norm is

\[
 \sqrt w\,\mathbb E|\xi|\,\Phi(K)
       \le{2\over\sqrt\pi D}\Phi(K)=O(n).                 \tag{14}
\]

Even the trivial bound `Phi(K)<=n(n-1)/2` suffices, together with
`D>=n-2` and the exact n=3 denominator.  No original-norm cap is
needed for (14).

The previously proved Gaussian finite-maximum comparison with
arbitrary deterministic offsets says that covariance operator error
delta in m dimensions costs at most `sqrt(2 delta m log N)` for
N states whose coefficient vectors have norm at most sqrt(m).
It follows from Gaussian convex order, adding independent delta-I
noise, and the finite Gaussian exponential maximum bound.  Here
`m=n(n-1)/2`, `N<=2^(n+1)`, and delta=40/sqrt(n), so (13)
costs `O(n^(5/4))` in the expected actual maximum.  Removing the
independent scalar Gaussian in (14) then costs O(n), by the
pointwise whole-norm Lipschitz bound and conditional Jensen.

Thus, if Y_h is the Gaussian with mean sK and exact covariance C_h,

\[
 |\mathbb E\mathcal M_I(Y_h)-\mathbb E\mathcal M_I(Z_h)|
                    \le C_1 n^{5/4}.                     \tag{15}
\]

## 6. Shifted-sign universality and the exact original norm

The independently proved
`NOTE_2026-09-05_SHIFTED_SIGN_GAUSSIAN_UNIVERSALITY.md`
applies with covariance operator constant three by (3), m<=n^2,
and at most 2^(n+1) augmented states.  It gives

\[
 |\mathbb E\mathcal M_I(B_h)-\mathbb E\mathcal M_I(Y_h)|
                       \le C_2 n^{16/11},                 \tag{16}
\]

uniformly in h and all deterministic internal energies.  Since
5/4<16/11, (15)--(16) prove (5) for n>=3.  There are no diagonal
matching terms: the entire construction uses unordered edges and
has zero matrix diagonal.

For n=2 use the independent-edge fallback specified in Section 1.
There is only one non-Gaussian coordinate, with mean sK and variance
1-s^2; its matched Gaussian is exactly the simple model (4).
The expectation of their absolute difference under an independent
coupling is at most two by centering and Cauchy-Schwarz.  The
maximum objective is 1-Lipschitz in that coordinate, so its error
is at most two, also covered by (5).

If K is an exact ORIGINAL norm minimizer at order n, every outcome
of B_h is an admissible signing of that same order.  Therefore (5)
gives the genuine, uniform source-optimality consequence

\[
 \mathbb E\Phi(Z_h)\ge\Phi(K)-C n^{16/11}.                 \tag{17}
\]

## 7. Full symmetric lift and an integrated actual-optimizer constraint

Assume n>=3 in this section, so D>0.  Use the Frobenius-orthonormal
basis of ALL real symmetric matrices, including the diagonal, and let
H_R be a centered symmetric Gaussian with covariance operator

\[
 {2\over D}(L^2\operatorname{Id}-\mathcal K),\qquad
                   \mathcal K(X)=KXK.                    \tag{18}
\]

This operator is PSD.  Its off-diagonal MATRIX ENTRIES, rather than
their orthonormal-basis coefficients, have exactly covariance R:
each off-diagonal entry is its basis coefficient divided by sqrt(2).
Consequently H_R is a Gaussian extension of the edge matrix G above.
For Boolean x,y,

\[
 \operatorname{Cov}(Q_{H_R}(x),Q_{H_R}(y))
       ={L^2(x^Ty)^2-(x^TKy)^2\over2D}.                  \tag{19}
\]

Here `Q_H(x)=x^THx/2`, also when H has a diagonal.  This equals
the off-diagonal quadratic form plus the SINGLE common scalar
`tau_R=tr(H_R)/2`.  Its variance is

\[
 \operatorname{Var}(\tau_R)
       ={nL^2-\operatorname{tr}(K^2)\over2D}
       ={n[L^2-(n-1)]\over2D}\le n/2.                    \tag{20}
\]

The numerator is nonnegative by L^2>=n-1, and the last inequality
uses n>=2.  No independence of this scalar and the off-diagonal
matrix is asserted or needed.

Similarly extend W to H_W with covariance operator 2 Id on the
symmetric matrix space.  Its off-diagonal entries are independent
standard Gaussians, its diagonal entries have variance two, and
`Var(tr(H_W)/2)=n/2`.  Take H_R,H_W independent, and set

\[
 \widetilde Z_h=sK+\sqrt{k}\,H_R+\sqrt v\,H_W.
\]

Its off-diagonal matrix is the model Z_h of (4).  The added common
scalar `tau_h=tr(tilde Z_h)/2` has variance at most
`(k+v)n/2=(1-s^2)n/2<=n/2`, and therefore

\[
 |\mathbb E\Phi(\widetilde Z_h)-\mathbb E\Phi(Z_h)|
                   \le\mathbb E|\tau_h|\le\sqrt{n/\pi}.   \tag{21}
\]

The analogous absolute pressure error is at most beta times this
quantity.  For a single phase WITHOUT an absolute-value augmentation,
its expected pressure is exactly unchanged, since the common scalar
factors out of the partition and has mean zero.  The augmented phase
is retained below; it is not cancelled.

For a fixed beta>0 define the actual Gaussian cosh pressure

\[
 P(h)=\mathbb E\log\mathbb E_x
       \cosh\!\left(\beta Q_{\widetilde Z_h}(x)\right).
\]

At each realized disorder, use its CURRENT posterior on (sigma,x),
where sigma is the sign augmenting the cosh and the energy is
`beta sigma Q_(tilde Z_h)(x)`.  Define

\[
 \Gamma=\langle\sigma xx^T\rangle,\qquad
 U=\|\Gamma\|_F^2
     =\langle\sigma\sigma'(x^Tx')^2\rangle,\qquad
 T_K=\operatorname{tr}(K\Gamma K\Gamma)
     =\langle\sigma\sigma'(x^TKx')^2\rangle.              \tag{22}
\]

The prime denotes an independent replica of that SAME current
posterior.  Gamma is symmetric but need not be PSD.  In particular,
the factors sigma sigma' in (22) cannot be replaced by one.

Let V_R,V_W be the contractions of the actual posterior covariance
of the symmetric-matrix observable `sigma xx^T/2` with the Gaussian
covariances of H_R,H_W respectively.  Formula (19) gives exactly

\[
 V_R={L^2(n^2-U)-4\langle Q_K(x)^2\rangle+T_K\over2D},
 \qquad V_W={n^2-U\over2}.                               \tag{23}
\]

Both quantities are nonnegative.  These formulas retain the exact
replica K-contraction and the entire augmented posterior.  The
`x^T K^2 x` terms of an off-diagonal-only computation have been
accounted for by the explicit common-scalar extension (20)-(21),
not silently dropped from an augmented covariance calculation.

There is a useful universal upper bound even for signed Gamma:

\[
 |\operatorname{tr}(K\Gamma K\Gamma)|
       \le\|\Gamma\|_F\|K\Gamma K\|_F
       \le L^2\|\Gamma\|_F^2.                            \tag{24}
\]

The first step is Frobenius Cauchy--Schwarz; the second uses the
operator bound on both factors K.  No positivity of Gamma is used.
Combining (23)-(24), the full noise contraction therefore obeys the
pointwise CURRENT-posterior upper

\[
 \boxed{\displaystyle
 kV_R+vV_W\le{1\over2}
  \left[(1-s^2+k/D)n^2
           -{4k\over D}\langle Q_K(x)^2\rangle-vU\right].}\tag{25}
\]

Thus a negative actual energy-square term and a posterior-overlap
subtraction remain.  This variance upper does not, by itself, give
an upper bound for every interpolation derivative: a negative
coefficient multiplying a variance reverses its inequality direction.

Write `a(h)=s'(h)=2 phi(h)`, so `k=a^2`,
`k'=-2hk`, and `v'=2hk-2sa`.  Gaussian covariance differentiation,
with K and R held FIXED, gives the exact identity

\[
 \begin{split}
 {P'(h)\over\beta}
 ={}&a\,\mathbb E\langle\sigma Q_K(x)\rangle
             +{\beta\over2}\mathbb E[k'V_R+v'V_W]\\
 ={}&a\,\mathbb E\langle\sigma Q_K(x)\rangle
   -{\beta\over2}(sa+hk/D)\mathbb E[n^2-U]\\
 &\hspace{12mm}+{\beta hk\over2D}
                \mathbb E[4\langle Q_K(x)^2\rangle-T_K].
 \end{split}                                                \tag{26}
\]

Denote the right side by J_beta(h).  Every moment in (26) is taken
after the Gaussian noise and the deterministic mean are present.

For completeness the derivative is integrable from any finite h to
infinity.  At fixed n,beta there are finitely many spin states;
their energy and covariance derivatives have uniform finite bounds.
The scalar functions a, k', and v' are absolutely integrable on
every such interval.  As h tends to infinity, s tends to one and
k,v tend to zero, so Gaussian finite-maximum bounds and the pressure
Lipschitz bound show

\[
 P(h)\longrightarrow P_K:=\log\mathbb E_x\cosh(\beta Q_K(x)),
 \qquad
 \int_h^\infty J_\beta(u)\,du={P_K-P(h)\over\beta}.       \tag{27}
\]

Now assume K is an EXACT ORIGINAL norm minimizer, with M=Phi(K).
Equation (17) and the trace restoration (21) give
`E Phi(tilde Z_h)>=M-C n^(16/11)-sqrt(n/pi)`, uniformly in h.
The finite maximum-term bound and `P_K<=beta M` therefore prove

\[
 \boxed{\displaystyle
 \int_h^\infty J_\beta(u)\,du
 \le C n^{16/11}+\sqrt{n/\pi}
                      +{(n+1)\log2\over\beta}.}          \tag{28}
\]

In particular, choose `beta=n^(-5/11)`.  The actual source-optimality
constraint has the quantitative, uniform form

\[
 \boxed{\displaystyle
 \int_h^\infty J_{n^{-5/11}}(u)\,du\le C' n^{16/11},
 \qquad h\in\mathbb R,}                                  \tag{29}
\]

with an absolute C'.  This uses original NORM optimality, not
pressure optimality at the auxiliary beta.  It is a bound on the
integrated, fully coupled posterior expression (26), not a pointwise
sign assertion.

## Scope

Equations (17) and (28)-(29) are actual SAME-order source-optimality
constraints.  Turning them into a doubled-order upper comparison
requires a further valid mapping or argument, not a reversal of
these inequalities.  No closure of the original MO convergence
question follows from this note alone.
