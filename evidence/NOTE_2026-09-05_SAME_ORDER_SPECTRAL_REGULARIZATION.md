# Same-order bounded-operator approximation for the original Boolean minimum

2026-09-05. **All-orders reduction for the ORIGINAL norm objective.**
Any norm-bounded complete signing can be replaced, at the same order,
by a complete signing with bounded normalized operator norm and an
arbitrarily small additive normalized Boolean-norm loss. This applies
to exact minimizers of m_N, not only to a pressure surrogate.

This is a one-sided objective approximation. It is NOT a bound on
`Phi(A'-A)` or on the uniform difference of the two energy functions.
The original convergence problem remains open.

## 1. The theorem

Let S_N be the complete symmetric zero-diagonal signings of order N,
and put

\[
 Q_A(x)=\frac12x^TAx,\qquad
 \Phi(A)=\max_{x\in\{-1,1\}^N}|Q_A(x)|,
 \qquad m_N=\min_{A\in\mathcal S_N}\Phi(A).
\]

Set

\[
 \kappa=\log(1+\sqrt2),\qquad \Gamma=\frac{4\pi}{\kappa}.
\]

**Theorem.** Suppose N>=2, A in S_N, and Phi(A)<=C N^(3/2).
For every K>0 there are a vertex set S and a signing A' in S_N such
that A' agrees with A on all edges inside T=[N] minus S and

\[
 |S|\le\frac{\Gamma C}{K}N,\qquad
 \|A'\|_{\rm op}\le(K+8)\sqrt N,                              \tag{1}
\]
\[
 \boxed{\quad
 \Phi(A')\le\Phi(A)+2\sqrt{\frac{\Gamma C}{K}}N^{3/2}.
 \quad}                                                       \tag{2}
\]

In particular, for any 0<epsilon<=1 one may take
`K=4 Gamma C/epsilon^2`, obtaining an objective loss at most
epsilon N^(3/2), changing only edges incident to at most
epsilon^2 N/4 vertices, and retaining

\[
 \|A'\|_{\rm op}
 \le\left(8+\frac{4\Gamma C}{\varepsilon^2}\right)\sqrt N.
                                                               \tag{3}
\]

The integer cardinality inequality is literal; if its right side is
less than one then the selected exceptional set is empty.

## 2. A trace-controlled diagonal majorizer

The elementary tensor-rounding and SDP argument in Section 2 of
`evidence/NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md` gives diagonal
matrices D^+,D^- with nonnegative entries satisfying

\[
 D^+-A\succeq0,\qquad D^-+A\succeq0,
 \qquad \operatorname{tr}D^\pm\le\frac{2\pi}{\kappa}\Phi(A).
                                                               \tag{4}
\]

For precision, that argument first proves the finite vector inequality
with constant pi/(2 kappa) by tensor lifts having inner products
`sin(kappa <u_i,v_j>)`; Gaussian sign rounding recovers exactly
`(2 kappa/pi)<u_i,v_j>`. The real bilinear cube norm of a symmetric
zero-diagonal A is at most 4 Phi(A), by polarization and multiaffinity.
Thus each semidefinite maximum

\[
 \max\{\operatorname{tr}(\pm A X):X\succeq0,
                                   \operatorname{diag}X=\mathbf1\}
\]

is at most `(2pi/kappa) Phi(A)`. Its diagonal dual is (4).
Strict primal and dual feasibility give strong duality; objective
sublevels of the dual are compact since their diagonal entries are
nonnegative. Hence the dual minima are attained. This is the proved
diagonal-majorizer result, not an assumed spectral-flatness claim.

Put D=D^++D^-. Then

\[
 D\pm A\succeq0,\qquad
 \operatorname{tr}D\le\Gamma\Phi(A)\le\Gamma C N^{3/2}.        \tag{5}
\]

No identity shift is needed here: positive semidefiniteness, rather
than the field-response proof's positive definiteness, is enough.
Choose

\[
 S=\{i:D_{ii}>K\sqrt N\},\qquad T=[N]\setminus S,
 \qquad k=|S|.
\]

The trace bound proves `k<=Gamma C N/K`. Taking principal restrictions
in (5), and using `D_T<=K sqrt(N) I`, gives

\[
 -K\sqrt N\,I\preceq A_T\preceq K\sqrt N\,I,
 \qquad \|A_T\|_{\rm op}\le K\sqrt N.                        \tag{6}
\]

If k=0, take A'=A and the theorem is already proved. Otherwise
recomplete every edge incident to S as follows.

## 3. One filler satisfies BOTH required bounds

Independently assign a fair sign to each of the

\[
 e=\binom N2-\binom{N-k}2
   =Nk-\frac{k(k+1)}2\le Nk
\]

edges incident to S. Let F be their symmetric zero-diagonal matrix,
with zero entries on T times T. Set A'=A_T extended by zeros, plus F.
Then A' is a complete signing of exactly the original order N.

For any fixed real unit vector v,

\[
 v^TFv=2\sum_{\{i,j\}\text{ incident to }S}\xi_{ij}v_i v_j.
\]

The squared sum of these independent sign coefficients is at most

\[
 4\sum_{i<j}v_i^2v_j^2
      =2\left(1-\sum_i v_i^4\right)\le2.
\]

The elementary estimate `cosh u<=exp(u^2/2)` and exponential Markov
therefore give

\[
 \Pr\{|v^TFv|>u\}\le2e^{-u^2/4}.                             \tag{7}
\]

A 1/4-net of the real unit sphere has at most 9^N points, by the
disjoint-ball volume argument. For a symmetric matrix its net maximum
controls the operator norm within a factor two: replacing a unit
maximizer by a net point changes the quadratic form by at most half
the operator norm. Substituting u=4 sqrt(N) in (7) gives

\[
 \Pr\{\|F\|_{\rm op}>8\sqrt N\}
       \le2e^{-(4-\log9)N}.                                  \tag{8}
\]

For a fixed Boolean spin vector x, the quantity Q_F(x) is a sum of
e independent fair signs, each with coefficient of magnitude one.
Thus

\[
 \Pr\{|Q_F(x)|>u\}\le2e^{-u^2/(2e)}.
\]

A union bound over the 2^N spin vectors, with u=2 sqrt(N e), proves

\[
 \Pr\{\Phi(F)>2\sqrt{Ne}\}
       \le2e^{-(2-\log2)N}.                                  \tag{9}
\]

The two failure bounds sum to less than one for N>=2. For example,
at N=2 they sum to `162 e^(-8)+8 e^(-4)<1`, and both terms decrease
as N increases. Consequently a SINGLE filler realization satisfies

\[
 \|F\|_{\rm op}\le8\sqrt N,
 \qquad \Phi(F)\le2\sqrt{Ne}\le2N\sqrt k.                    \tag{10}
\]

No independence of these two success events is needed. The construction
is a probabilistic existence proof; it requires no sampling run.

## 4. The objective comparison is one-sided and stays at order N

For every T-spin vector y, uniform averaging over the deleted S spins
gives `E_(x_S) Q_A(x_S,y)=Q_(A_T)(y)`, because A has zero diagonal.
Therefore

\[
 \Phi(A_T)\le\Phi(A).                                       \tag{11}
\]

Combining (6), (10), and (11) yields

\[
 \|A'\|_{\rm op}\le(K+8)\sqrt N,
\]
\[
 \Phi(A')\le\Phi(A_T)+\Phi(F)
    \le\Phi(A)+2N\sqrt k
    \le\Phi(A)+2\sqrt{\Gamma C/K}\,N^{3/2}.
\]

These prove (1)--(2). The original incident edges have been discarded;
their OLD energy is not bounded by (10). Thus this argument provides
no estimate on Phi(A'-A), and none is used.

## 5. Uniform approximation of the original normalized minima

For L>=8 define the operator-constrained normalized optimum

\[
 \alpha_N^{(L)}=N^{-3/2}
   \min\{\Phi(A):A\in\mathcal S_N,\ \|A\|_{\rm op}\le L\sqrt N\},
 \qquad \alpha_N=m_N/N^{3/2}.
\]

The constrained class is nonempty: the all-edge version of (8)--(10)
already supplies a complete signing of operator norm at most 8 sqrt(N)
and Boolean norm at most 2 N^(3/2). In particular m_N<=2 N^(3/2)
for every N>=2. Applying the theorem to an exact norm minimizer with
this uniform C=2 gives, for EVERY N>=2 and K>0,

\[
 \boxed{\quad
 0\le\alpha_N^{(K+8)}-\alpha_N
          \le2\sqrt{\frac{2\Gamma}{K}}.
 \quad}                                                       \tag{12}
\]

One may replace C=2 by any better proved uniform or eventual bound
on alpha_N; the displayed version needs no asymptotic construction
input at all.

For any prescribed sequence K_N>0 tending to infinity, however slowly,
(12) supplies actual complete signings A_N' satisfying

\[
 \Phi(A_N')=m_N+o(N^{3/2}),\qquad
 \|A_N'\|_{\rm op}\le(K_N+8)\sqrt N.                         \tag{13}
\]

Thus the original norm minimum is approximable at leading order by
signings whose normalized operator norm grows arbitrarily slowly.
For fixed K, (12) instead gives an approximation error uniform in N
inside the genuinely bounded-operator class.

## 6. Exact remaining convergence reduction

For a bounded real sequence f_N write
`osc_N f_N=limsup_N f_N-liminf_N f_N`. Equation (12) implies

\[
 \left|\operatorname{osc}_N\alpha_N^{(K+8)}
              -\operatorname{osc}_N\alpha_N\right|
       \le2\sqrt{\frac{2\Gamma}{K}}.                         \tag{14}
\]

Consequently the original convergence statement is equivalent to

\[
 \boxed{\quad
 \lim_{K\to\infty}
       \operatorname{osc}_N\alpha_N^{(K+8)}=0.
 \quad}                                                       \tag{15}
\]

In particular, proving convergence of the operator-constrained
normalized optimum for every sufficiently large FIXED K would prove
the original limit. Neither that fixed-K convergence nor (15) is
established here. Letting K grow with N in (13) alone does not prove it.

## Duplication and scope audit

The exact-object search covered solution.md, CORE.md, the proposition
audit, and the evidence notes. Proposition 15.12 is a conditional
pointwise spectral bound under a cube/sphere-ratio hypothesis; it is
not this unconditional same-order replacement statement. Earlier
covariance spectral truncations concern U,V, not the original signing.
No prior same-order trace-majorizer trimming and recompletion theorem
was found in those sources.

The only imported local proof is the established elementary tensor/SDP
majorizer in (4). The filler estimates, joint success probability,
objective comparison, and convergence reduction are proved above.
No signing census, solver, numerical matrix experiment, or Gaussian
simulation was used.
