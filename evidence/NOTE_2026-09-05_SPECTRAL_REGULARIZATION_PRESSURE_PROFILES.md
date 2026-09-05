# Pressure-profile extension of same-order spectral regularization

2026-09-05. **Proved same-order profile-preserving construction.**
This extends the reviewed norm-regularization theorem using the SAME
incident-edge filler. It does not prove order transport or convergence.

## 1. Statement

Use the notation of
`evidence/NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md`:

\[
 \Phi(A)=\max_{x\in\{-1,1\}^N}|Q_A(x)|,\qquad
 \Gamma=\frac{4\pi}{\log(1+\sqrt2)}.
\]

Suppose N>=2, A is a complete zero-diagonal signing, and
`Phi(A)<=C N^(3/2)`. Fix K>0 and a finite temperature cutoff U>=0.
Define

\[
 C'=C+2\sqrt{\Gamma C/K},\qquad
 q=\lfloor UN\rfloor+1,\qquad
 R_N=\log(8q)+C+C'.                                      \tag{1}
\]

There is ONE complete signing A' of the same order, agreeing with A
away from edges incident to a set S, such that

\[
 |S|\le\Gamma C N/K,\qquad
 \|A'\|_{\rm op}\le(K+8)\sqrt N,\qquad
 \Phi(A')\le C'N^{3/2},                                  \tag{2}
\]

and, simultaneously for every c in [0,U] and both sigma in {+1,-1},

\[
 \boxed{\quad
 \log Z_{A'}^\sigma(c/\sqrt N)
 \le \log Z_A^\sigma(c/\sqrt N)
       +\frac{\Gamma C c^2}{2K}N+R_N,
 \quad}                                                  \tag{3}
\]

where `Z_A^sigma(beta)=E_x exp(sigma beta Q_A(x))` uses the uniform
spin measure. More precisely, writing T=[N] minus S and
`e=binom(N,2)-binom(|T|,2)`, the added term in (3) can be replaced by

\[
 e\log\cosh(c/\sqrt N)+R_N.                              \tag{4}
\]

The identical upper comparison holds for the half-product pressure
`a_A=(log Z_A^+ +log Z_A^-)/2` and symmetric pressure
`F_A=log((Z_A^+ +Z_A^-)/2)`. Thus for fixed C,K,U the normalized
profile loss is at most `Gamma C c^2/(2K)+O(log(N)/N)`, uniformly in c.

## 2. The same filler has enough simultaneous probability slack

Take the trace-majorizer exceptional set S from the regularization
theorem. If S is empty, choose A'=A. Otherwise, independently assign
fair signs to its e incident edges, leaving A_T fixed. The proved
operator and Boolean failure bounds for this filler F are

\[
 p_N=2e^{-(4-\log9)N}+2e^{-(2-\log2)N}.                    \tag{5}
\]

They decrease with N, and for N>=2 satisfy

\[
 p_N\le162e^{-8}+8e^{-4}<\frac12.
\]

For example, `exp(4)>32` gives the upper bound
`162/1024+8/32<1/2`. On their simultaneous success event, (2) holds,
and the stronger original one-sided bound
`Phi(A')<=Phi(A)+2 sqrt(Gamma C/K) N^(3/2)` also holds.

For each fixed beta>=0 and each sigma, independence of the new edges
gives the EXACT identity

\[
 \mathbb E_F\frac{Z_{A'}^\sigma(\beta)}{Z_{A_T}^\sigma(\beta)}
       =(\cosh\beta)^e.                                  \tag{6}
\]

Indeed, at fixed spins every new edge contributes the factor
`E exp(sigma beta xi_ij x_i x_j)=cosh(beta)`; averaging the remaining
old energy over spins gives Z_(A_T)^sigma. Thus Markov's inequality
gives

\[
 \Pr_F\left\{\log Z_{A'}^\sigma(\beta)
       -\log Z_{A_T}^\sigma(\beta)
       >e\log\cosh\beta+s\right\}\le e^{-s}.              \tag{7}
\]

Apply this at the q grid points `c_j=j/N`, `0<=j<=floor(UN)`, with
`beta_j=c_j/sqrt(N)`, for both phases, and put `s=log(8q)`.
The union of these 2q pressure failures has probability at most 1/4.
Together with (5), all required successes therefore have probability
at least `1-p_N-1/4>1/4`. Select ONE such filler. No independence
between success events, and no sampling computation, is required.

For each fixed T-spin vector, the old full energy has conditional
S-spin mean Q_(A_T). Jensen's inequality consequently gives

\[
 \log Z_{A_T}^\sigma(\beta)\le\log Z_A^\sigma(\beta).      \tag{8}
\]

Hence the selected filler satisfies the desired phase comparison at
every grid point, with error `e log cosh(beta_j)+log(8q)`.

## 3. Uniform interpolation between temperatures

For any fixed signing B,

\[
 \left|\frac{d}{dc}\log Z_B^\sigma(c/\sqrt N)\right|
       \le\Phi(B)/\sqrt N.                               \tag{9}
\]

The simultaneously selected Boolean bound in (2) therefore makes the
two phase pressures for A and A' Lipschitz in c with constants CN and
C'N, respectively. Given c in [0,U], choose `c_j=floor(Nc)/N`.
Then `0<=c-c_j<1/N`, and the grid comparison implies

\[
 \log Z_{A'}^\sigma(c/\sqrt N)-\log Z_A^\sigma(c/\sqrt N)
 \le e\log\cosh(c_j/\sqrt N)+\log(8q)+C+C'
 \le e\log\cosh(c/\sqrt N)+R_N.                           \tag{10}
\]

The last step uses c_j<=c and monotonicity of log cosh on [0,infinity).
Finally, `log cosh u<=u^2/2` and `e<=N|S|<=Gamma C N^2/K` prove (3).

If both phase partition functions increase by factors at most exp(D),
their geometric mean and their arithmetic mean increase by factors
at most exp(D). This proves the half-product and symmetric-pressure
claims directly; no additional probabilistic constraints are needed.

## Scope

The same signing A' works for the entire prescribed compact
temperature interval and retains both norm and operator conclusions.
The assertions compare actual pressures of actual integral signings.
They do not bound the discarded old energy or `Phi(A'-A)`, identify
the two pressure optima, imply phase balance, or provide a selected
order-extension theorem. Original convergence remains OPEN.
