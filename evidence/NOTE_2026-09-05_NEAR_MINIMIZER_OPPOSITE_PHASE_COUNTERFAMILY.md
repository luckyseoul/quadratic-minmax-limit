# Near-minimizing complete signings with a superquadratic opposite-phase trace

## Scope and theorem

Fix \(c>0\). For every sufficiently large integer \(N\), there exist actual
complete symmetric zero-diagonal signings, at \(\beta_N=c/\sqrt N\), with
\[
T_A:=\operatorname{tr}(A U_A A V_A)=\Omega_c(N^{9/4}),            \tag{1}
\]
where \(U_A,V_A\) are the actual zero-field Gibbs covariances at
\(+\beta_N,-\beta_N\), and with either of the following near-minimality
properties:
\[
0\le\Phi(A)-m_N=O_c(N^{11/8}),                                  \tag{2}
\]
or
\[
0\le a_A(\beta_N)-R_N(\beta_N)=O_c(N^{7/8}).                     \tag{3}
\]
Here
\[
\begin{gathered}
Q_A(x)=\tfrac12x^TAx,\quad
\Phi(A)=\max_{x\in\{-1,1\}^N}|Q_A(x)|,\quad m_N=\min_A\Phi(A),\\
Z_A(t)=2^{-N}\sum_xe^{tQ_A(x)},\quad
a_A(t)=\tfrac12(\log Z_A(t)+\log Z_A(-t)),\quad
R_N(t)=\min_A a_A(t).
\end{gathered}
\]

Statements (2) and (3) are obtained by two potentially different choices
of the old base signing. No assertion is made that the same family
satisfies both. Their normalized errors are \(O_c(N^{-1/8})\), so these
are leading-order near-minimizers, not merely hosts with the correct
norm or pressure scale.

Both families also satisfy \(\Phi(A)=O_c(N^{3/2})\),
\(a_A(\beta_N)=O_c(N)\), and
\[
\frac{T_A}{2a_A'(\beta_N)/\beta_N}\longrightarrow\infty.         \tag{4}
\]
Nevertheless, the constructed signings are provably not even edge-local
half-product minima for large \(N\). Neither conclusion asserts exact
norm minimality. The construction therefore rules out replacing actual
minimality in a trace argument by generic leading-order near-minimality.
It does **not** disprove a trace bound for exact global minima, for
actual edge-local minima, or for quantitatively tighter approximation
classes. It does not settle the MO limit or a cross-order comparison.

The proof is probabilistic existence for complete signings at every
sufficiently large order. No arbitrary cavity laws, finite computations,
censuses, simulations, or external random-matrix theorems are used.

## 1. Parameters and a universal planted extension

Put
\[
A_c=\frac{4\log2}{c}+8c,\qquad
K=\max\{8,\lceil1024A_c\rceil\}.                                \tag{5}
\]
For each sufficiently large integer \(N\), define
\[
\ell=2\left\lceil\frac{K\sqrt N}{2}\right\rceil,\quad
s=\lfloor N^{1/4}\rfloor,\quad m=2s\ell,\quad L=N-m.             \tag{6}
\]
Then \(m<N\), \(\ell\ge8\), \(m=\Theta_c(N^{3/4})\), and
\[
K\sqrt N\le\ell\le(K+2)\sqrt N,\qquad
m\le2(K+2)N^{3/4}.                                             \tag{7}
\]
Fix any complete signing \(B\) on the \(L\) old vertices. It remains
unchanged. The argument below applies uniformly to every such \(B\).

The \(m\) new vertices form \(s\) modules of \(2\ell\) vertices each.
In each module there are \(\ell\) twin pairs, half in each of two
communities. Each community thus contains \(\ell\) spins. The internal
signs are the constant community pattern
\[
J=\begin{pmatrix}1&1\\1&-1\end{pmatrix},
\]
with zero diagonal.

There are \(t=m/2=s\ell\) new twin pairs in total. Independently choose:

- a uniform sign \(R_{pi}\) for every new pair \(p\) and old vertex \(i\);
- a uniform sign \(W_{pq}=W_{qp}\) for every pair of new twin pairs in
  different modules.

Set \(W_{pq}=0\) inside modules, including the diagonal. The old-to-pair
edge signs are \(R_{pi}(1,-1)\); the cross-module twin-pair block is
\[
W_{pq}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.                    \tag{8}
\]
All edges are now specified. The resulting matrix \(A\) is symmetric,
zero-diagonal, and has a sign on every off-diagonal entry.

For new spins write
\[
u_p=x_{p,+}-x_{p,-}\in\{-2,0,2\}.
\]
The filler energy, consisting of all edges joining distinct modules or
joining new and old vertices, is exactly
\[
Q_F(z,x)=\sum_{p,i}R_{pi}u_pz_i+
         \sum_{\substack{p<q\\\text{different modules}}}
                  W_{pq}u_pu_q.                              \tag{9}
\]
Here \(z\) is the old spin vector. Thus
\(Q_A(z,x)=Q_B(z)+Q_{\rm int}(x)+Q_F(z,x)\).

## 2. An elementary high-probability filler norm bound

Use orthonormal coordinates on the old subspace and on the pair-odd
subspace of the new vertices. The filler vanishes on every pair-even
vector. On old plus pair-odd coordinates its matrix is
\[
\mathcal F=
\begin{pmatrix}
0&\sqrt2R^T\\
\sqrt2R&2W
\end{pmatrix}.                                                  \tag{10}
\]
Its dimension is \(d=L+t\le N\).

For a fixed unit vector \(v=(a,b)\), the independent Rademacher
coefficients in \(v^T\mathcal Fv\) have squared sum at most
\[
8\sum_{p,i}b_p^2a_i^2+
16\sum_{p<q}b_p^2b_q^2
\le8\|a\|^2\|b\|^2+8\|b\|^4\le8.                              \tag{11}
\]
For completeness, independent signs \(\varepsilon_j\) obey
\[
\mathbb E e^{h\sum_jq_j\varepsilon_j}
 =\prod_j\cosh(hq_j)
 \le e^{h^2\sum_jq_j^2/2}.
\]
Exponential Markov, applied also to the negative sum, consequently gives
\[
\Pr\{|v^T\mathcal Fv|>u\}\le2e^{-u^2/16}.                      \tag{12}
\]

A maximal \(1/4\)-separated subset of the unit sphere is a \(1/4\)-net
of size at most \(9^d\), by disjoint-ball volume comparison. For a
symmetric matrix \(M\), such a net satisfies
\[
\|M\|_{\rm op}\le2\max_{v\text{ in the net}}|v^TMv|,
\]
because replacing a maximizing unit vector by a net point changes the
quadratic form by at most \(2(1/4)\|M\|_{\rm op}\).
Taking \(u=8\sqrt N\) in (12) and using the union bound proves
\[
\Pr\{\|\mathcal F\|_{\rm op}>16\sqrt N\}
 \le2e^{-(4-\log9)N}=o(1).                                    \tag{13}
\]
Call the complementary event \(\mathcal G_{\rm op}\).

On \(\mathcal G_{\rm op}\), each block of (10) has operator norm at most
\(16\sqrt N\). The old spin vector has norm \(\sqrt L\); the pair-odd
projection of the new spins has norm at most \(\sqrt m\). Since the
old-old filler block is zero, this yields the sharper Boolean bound
\[
|Q_F(z,x)|\le16N\sqrt m+8m\sqrt N.                              \tag{14}
\]
Using only \(\|Q_F\|_\infty\le N\|\mathcal F\|_{\rm op}/2\) would lose
the smallness supplied by \(m=o(N)\).

For each module its exact internal energy at community magnetizations
\(S_1,S_2\in[-\ell,\ell]\) is
\[
q(S_1,S_2)=\tfrac12(S_1^2+2S_1S_2-S_2^2),\qquad |q|\le\ell^2.
                                                                    \tag{15}
\]
The diagonal corrections cancel between the two equal-sized communities.
Therefore, on \(\mathcal G_{\rm op}\), uniformly in all spins,
\[
|Q_A(z,x)-Q_B(z)|\le
E_N:=\frac{m\ell}{2}+16N\sqrt m+8m\sqrt N
\le C_K N^{11/8},                                              \tag{16}
\]
where one may take
\[
C_K=(K+2)^2+16(K+2)+16\sqrt{2(K+2)}.
\]

## 3. Conditional entropy uses only the new vertices

Define the module deficits
\[
\begin{aligned}
\delta^+(S)&=\ell^2-q(S)
 =\ell^2-S_1^2+\tfrac12(S_2-S_1)^2,\\
\delta^-(S)&=\ell^2+q(S)
 =\ell^2-S_2^2+\tfrac12(S_1+S_2)^2.                             \tag{17}
\end{aligned}
\]
They are nonnegative. Let \(P=s\ell^2=m\ell/2\).

Condition on any old spin configuration \(z\). The conditional new-spin
Hamiltonian is \(H_z=Q_{\rm int}+Q_F(z,\cdot)\). Configurations with
every module at either corner \(\pm(\ell,\ell)\) have all twins agreeing,
so \(Q_F=0\) and \(H_z=P\). Corners \(\pm(\ell,-\ell)\) similarly give
\(H_z=-P\). These witnesses are valid for every \(z,R,W\).

For a Gibbs law on \(2^m\) states at inverse temperature \(\beta>0\),
\[
\mathbb E H\ge\max H-\frac{m\log2}{\beta}.
\]
Indeed, the normalized log partition is at least
\(\beta\max H-m\log2\), and \(\beta\mathbb EH\) is that log partition
plus nonnegative relative entropy. Apply this to \(H_z\) and \(-H_z\),
then average over the actual phase-dependent old-spin marginal. With
\(\beta=\beta_N\), this proves
\[
\begin{aligned}
\sum_b\mathbb E_+\delta_b^+
 &\le\frac{m\log2}{\beta}+\mathbb E_+Q_F,\\
\sum_b\mathbb E_-\delta_b^-
 &\le\frac{m\log2}{\beta}-\mathbb E_-Q_F.                        \tag{18}
\end{aligned}
\]
The signs in the second line matter. No independence of the old marginal
from the filler is assumed: the conditional inequalities hold pointwise
in the old spins.

## 4. Filler-sign averaging controls the thermal deficits

Here all Gibbs expectations are for the actual full Hamiltonian, including
the fixed old base \(B\). Consider a filler sign \(\xi\) whose coefficient
in \(Q_F\) is the spin observable \(Y\), with \(|Y|\le b\). Condition
on every other filler sign and interpolate its coefficient continuously
over \(u\in[-1,1]\). Then
\[
\frac d{du}\langle Y\rangle_+=\beta\operatorname{Var}_+(Y),
\qquad
\frac d{du}\langle Y\rangle_-=-\beta\operatorname{Var}_-(Y).
\]
In particular,
\[
0\le\mathbb E_\xi[\xi\langle Y\rangle_+]\le\beta b^2,\qquad
0\le-\mathbb E_\xi[\xi\langle Y\rangle_-]\le\beta b^2.            \tag{19}
\]
This follows by writing each sign expectation as half the difference
between its two endpoint values and using \(\operatorname{Var}(Y)\le b^2\).

In (9), the \(tL\) pair-old signs have \(b=2\); the cross-module pair
signs have \(b=4\), and number at most \(t^2/2\). Consequently,
\[
\sum b^2\le4tL+8t^2
          =2mL+2m^2=2mN.                                     \tag{20}
\]
Let
\[
d_b^+=\mathbb E_+\delta_b^+,\quad
d_b^-=\mathbb E_-\delta_b^-,\quad
\Delta=\sum_b(d_b^++d_b^-).
\]
Taking filler-sign expectations in (18), and using (19)--(20), gives
\[
\mathbb E_{R,W}\Delta
 \le\frac{2m\log2}{\beta}+4\beta mN.
\]
As \(\ell\ge K\sqrt N\),
\[
\mathbb E_{R,W}\frac{\Delta}{P}
 \le\frac{4\log2}{\beta\ell}+\frac{8\beta N}{\ell}
 \le\frac{4\log2/c+8c}{K}\le\frac1{1024}.                       \tag{21}
\]
Markov therefore gives
\[
\Pr\{\Delta\le P/256\}\ge3/4.                                  \tag{22}
\]
Call this event \(\mathcal G_{\rm th}\). Its definition uses both full
Gibbs phases for the realized filler.

For sufficiently large \(N\), (13) gives
\(\Pr(\mathcal G_{\rm op})\ge3/4\). Therefore
\(\Pr(\mathcal G_{\rm op}\cap\mathcal G_{\rm th})\ge1/2\).
No independence of these two events is needed. Fix any realized complete
signing in the intersection.

Put \(\epsilon=1/128\). A module is good if
\(d_b^+,d_b^-\le\epsilon\ell^2\). Equation (22) implies that at least
\(s/2\) module indices are good for both phases simultaneously.

## 5. Symmetry makes module trace contributions additive

Let \(P_b\) exchange every twin pair in module \(b\). The transformation
\(F_b:x_b\mapsto-P_bx_b\), leaving all other spins unchanged, preserves
the full Hamiltonian. It preserves internal energy, while each twin
difference \(u_p=x_{p,+}-x_{p,-}\) is fixed, so it preserves (9) and
the unchanged old base.

Let \(\mathcal E_b\) be the pair-even subspace of module \(b\), and
let \(\mathcal R\) be the direct sum of the old-coordinate subspace and
the new pair-odd subspace. On \(\mathcal E_b\), \(F_b\) acts as \(-I\);
on \(\mathcal R\) and all other \(\mathcal E_d\), it acts as \(I\).
Both Gibbs laws are invariant under every \(F_b\). Their means vanish
by global spin reversal, and their actual covariance matrices therefore
have the exact block forms
\[
U=\left(\bigoplus_b U_b\right)\oplus U_{\mathcal R},\qquad
V=\left(\bigoplus_b V_b\right)\oplus V_{\mathcal R}.              \tag{23}
\]
The matrix \(A\) preserves the same orthogonal decomposition:
\[
A=\left(\bigoplus_b C_b\right)\oplus A_{\mathcal R}.
\]
Indeed, the filler annihilates pair-even vectors and internal matrices
preserve the twin parity subspaces. Consequently,
\[
\begin{aligned}
T_A
 &=\sum_b\operatorname{tr}(C_bU_bC_bV_b)
   +\operatorname{tr}(A_{\mathcal R}U_{\mathcal R}
                      A_{\mathcal R}V_{\mathcal R})\\
 &\ge\sum_b\operatorname{tr}(C_bU_bC_bV_b).                     \tag{24}
\end{aligned}
\]
All terms are nonnegative because
\(\operatorname{tr}(CUCV)=\|U^{1/2}CV^{1/2}\|_F^2\).
No independence of full module configurations, and no independence
within the remainder sector, has been asserted.

## 6. A good module supplies a fourth-power contribution

Take independent full-host replicas \(x\sim\mu_+\), \(y\sim\mu_-\),
and fix a good module. Write its magnetizations as
\[
(S_1(x),S_2(x))=(a,a+e),\qquad
(S_1(y),S_2(y))=(-b+f,b).
\]
Equation (17) gives
\[
\begin{gathered}
\mathbb E_+a^2,\ \mathbb E_-b^2\ge(1-\epsilon)\ell^2,\qquad
\mathbb E_+e^2,\ \mathbb E_-f^2\le2\epsilon\ell^2,\qquad
|a|,|b|\le\ell.                                                \tag{25}
\end{gathered}
\]
The bilinear form from the constant community blocks is
\[
W_0=(a,a+e)J(-b+f,b)^T=-2ab-2eb+2af+ef.
\]
Independence of the two replicas and the \(L^2\) triangle inequality
yield
\[
\begin{aligned}
\|W_0\|_2
&\ge2\sqrt{\mathbb E_+a^2\,\mathbb E_-b^2}
 -2\ell\sqrt{\mathbb E_+e^2}
 -2\ell\sqrt{\mathbb E_-f^2}
 -\sqrt{\mathbb E_+e^2\,\mathbb E_-f^2}\\
&\ge(2-4\epsilon-4\sqrt{2\epsilon})\ell^2
 =\frac{47}{32}\ell^2.                                        \tag{26}
\end{aligned}
\]

Let \(x_{\mathcal E},y_{\mathcal E}\) be the orthogonal pair-even
projections in this module. Projection preserves the community sums.
The true internal matrix is the constant community-block matrix minus
the diagonal \(D\) with entries \(+1\) in the first community and
\(-1\) in the second. Thus, with \(C_b\) interpreted as the restricted
internal operator in the inherited Euclidean inner product,
\[
x_{\mathcal E}^TC_by_{\mathcal E}
 =W_0-x_{\mathcal E}^TDy_{\mathcal E},\qquad
|x_{\mathcal E}^TDy_{\mathcal E}|\le2\ell.                       \tag{27}
\]
Since \(\ell\ge8\), (26)--(27) imply
\[
\operatorname{tr}(C_bU_bC_bV_b)
=\mathbb E_{+,-}(x_{\mathcal E}^TC_by_{\mathcal E})^2
\ge\ell^4.                                                     \tag{28}
\]
Combining (24), (28), and the good-module count proves
\[
\boxed{\quad
T_A\ge\frac s2\ell^4
 \ge\frac{K^4}{4}N^{9/4}
\quad}                                                        \tag{29}
\]
for sufficiently large \(N\), using \(s\ge N^{1/4}/2\).
Together with (16), this establishes the universal extension theorem:
every fixed old complete signing \(B\) has an extension with the uniform
energy error (16) and the superquadratic trace (29).

## 7. The two separate near-minimizer corollaries

First choose \(B\) to attain \(m_L\). Restriction gives \(m_L\le m_N\):
for any complete signing, averaging over the outside spins preserves
the induced quadratic form, so its Boolean norm cannot increase.
Equation (16) now gives
\[
0\le\Phi(A)-m_N
 \le m_L+E_N-m_N\le E_N=O_c(N^{11/8}).                          \tag{30}
\]

Alternatively, choose \(B\) to attain \(R_L(\beta_N)\). The inverse
temperature here is exactly \(\beta_N\), not \(c/\sqrt L\).
For any extension \(A\) of any old signing \(B\), every added quadratic
term averages to zero over independent uniform new spins. Jensen
therefore yields, pointwise in the old spins,
\[
\mathbb E_{\rm new}e^{\pm\beta(Q_A-Q_B)}\ge1.
\]
Hence \(Z_A(\pm\beta)\ge Z_B(\pm\beta)\),
\(a_A(\beta)\ge a_B(\beta)\), and
\(R_N(\beta)\ge R_L(\beta)\). Conversely, the uniform error (16)
implies \(a_A(\beta)\le a_B(\beta)+\beta E_N\). Thus
\[
0\le a_A(\beta_N)-R_N(\beta_N)
 \le\beta_NE_N=O_c(N^{7/8}).                                   \tag{31}
\]
No unproved relation between different-order normalized optima is used.
The two choices of \(B\), and the successful fillers, need not coincide.

## 8. Norm and pressure scales, radial divergence, and scope

The auxiliary scale claims can also be obtained without outside results.
Averaging complete edge signs independently, and applying Jensen to
each one-sided log partition, gives
\[
R_N(\beta)\le\binom N2\log\cosh\beta.
\]
At \(\beta=\beta_N\), this is at most \(c^2(N-1)/4\).
For any signing \(A\), the uniform mean of \(Q_A\) is zero, so both
one-sided log partitions are nonnegative. A state attaining \(\Phi(A)\)
then gives
\[
\Phi(A)\le\frac{2a_A(\beta)+N\log2}{\beta}.                      \tag{32}
\]
Applying (32) to a pressure minimizer shows \(m_N=O_c(N^{3/2})\).
Equation (30) therefore gives this norm cap to the first family, and
its pressure is at most \(\beta_N\Phi(A)=O_c(N)\).
For the second family, (31) gives \(a_A(\beta_N)=O_c(N)\), and
(32) gives its norm cap \(O_c(N^{3/2})\).

For every nonconstant complete quadratic Hamiltonian and \(\beta>0\),
\[
0<a_A'(\beta)
 =\tfrac12(\mathbb E_+Q_A-\mathbb E_-Q_A)\le\Phi(A).
\]
Thus \(2a_A'(\beta_N)/\beta_N=O_c(N^2)\) in either family, whereas
(29) is \(\Omega_c(N^{9/4})\). This proves (4), with a ratio lower
bound of order \(N^{1/4}\).

Finally, these examples do not satisfy the actual edge-local
half-product inequalities. In a good module, (25) implies
\[
\mathbb E_+S_1S_2\ge\frac{111}{128}\ell^2,\qquad
\mathbb E_-S_1S_2\le-\frac{111}{128}\ell^2.                      \tag{33}
\]
Summing over its \(\ell^2\) cross-community edges gives an edge with
\(A_{ij}=+1\) and \(U_{ij}-V_{ij}\ge111/64\).
At an edge-local half-product minimum, the exact ratio under flipping
an edge \(e\) would instead require
\[
\begin{aligned}
&(\cosh(2\beta)-A_eU_e\sinh(2\beta))
 (\cosh(2\beta)+A_eV_e\sinh(2\beta))\ge1,\\
&\hspace{8mm}
A_e(U_e-V_e)\le\tanh(2\beta)(1-U_eV_e)\le2\tanh(2\beta).
                                                                    \tag{34}
\end{aligned}
\]
Since \(\beta_N\to0\), (33) violates (34) for all sufficiently large
orders. Exact global half-product minima, actual edge-local minima,
and exact norm-minimizer assertions remain outside the conclusions.

