# Bipartite SDP optimality, spectral mass, and complementarity

2026-09-05. This note proves structural constraints on an actual complete
cross signing B and on its actual optimal SDP diagonal. It distinguishes
the Boolean cut norm beta(B), the vector SDP value tau(B), and original
conditional signing optimality. No small-gap hypothesis is inferred from
conditional optimality, and no sharp Gaussian upper is concluded.

All arguments are analytic. No scalar calculation, signing search,
simulation, or other tool-based mathematical computation is used.

## 1. Actual SDP and its balanced optimal diagonal

Let B be a real n by n sign matrix, n>=1, and define
\[
 \beta(B)=\max_{x,y\in\{-1,1\}^n}|x^TBy|,\qquad
 H_B=\begin{pmatrix}0&B\\B^T&0\end{pmatrix}.
\]
Let
\[
 \tau(B)=\max_{\|u_i\|=\|v_j\|=1}
                   \sum_{ij}B_{ij}\langle u_i,v_j\rangle
 =\frac12\max_{X\succeq0,\ \operatorname{diag}X=1}
                   \operatorname{tr}(H_BX).             \tag{1}
\]
A Gram realization in dimension at most 2n suffices. The finite diagonal
dual gives an attained optimum
\[
 2\tau(B)=\min\{\operatorname{tr}D:D\text{ diagonal},\
                                      D-H_B\succeq0\}. \tag{2}
\]
Both programs are strictly feasible; dual objective sublevels are
compact because the diagonal coordinates are nonnegative.

Fix ANY optimal diagonal
\(D=\operatorname{diag}(D_r,D_c)=\operatorname{diag}(d_1,\ldots,d_n,
e_1,\ldots,e_n)\). Every d_i,e_j is positive: the corresponding
2-by-2 principal minor gives \(d_i e_j\ge B_{ij}^2=1\).
Bipartite sign conjugation also gives \(D+H_B\succeq0\).

For each a>0, the diagonal \(\operatorname{diag}(aD_r,a^{-1}D_c)\)
is feasible, by congruence with \(\operatorname{diag}(\sqrt a I,
a^{-1/2}I)\). Optimality at a=1 therefore implies
\[
                 \sum_i d_i=\sum_j e_j=\tau(B).         \tag{3}
\]
This balance holds for every optimum, not merely after selecting one.

Schur complements in (2), together with every \(B_{ij}^2=1\), give
\[
 d_i\ge\sum_j e_j^{-1},\qquad
 e_j\ge\sum_i d_i^{-1}.                                  \tag{4}
\]
Consequently, putting \(a=\tau(B)\),
\[
 d_i,e_j\ge n^2/a,\qquad
 \sum_i d_i^{-1},\ \sum_j e_j^{-1}\le a/n,\qquad
                         a\ge n^{3/2}.                 \tag{5}
\]
For example, the first bound follows from Cauchy--Schwarz on the
opposite block, and summing it yields \(a^2\ge n^3\).

## 2. A canonical primal and the actual third singular moment

Let \(|B|=(B^TB)^{1/2}\), and let \(\sigma_1,\ldots,\sigma_n\) be
the singular values, including zeros. Define
\[
             Z=\frac1{\sqrt n}\begin{pmatrix}B\\|B|\end{pmatrix},
                    \qquad X_0=ZZ^T.                  \tag{6}
\]
Every row of Z has norm one: both \(BB^T\) and \(|B|^2=B^TB\)
have diagonal n. Thus X_0 is genuinely feasible in (1). Directly,
\[
 \frac12\operatorname{tr}(H_BX_0)
                    ={1\over n}\operatorname{tr}|B|^3.
\]
In particular
\[
 \boxed{\tau(B)\ge {1\over n}\sum_j\sigma_j^3
                         \ge n^{3/2}.}                 \tag{7}
\]
The last step uses \(\sum_j\sigma_j^2=n^2\) and the power-mean
inequality. Define the actual nonnegative canonical-primal gap
\[
              g=\tau(B)-{1\over n}\sum_j\sigma_j^3.     \tag{8}
\]
There is no assertion that g is small for a conditional signing
optimizer. It is a specific, well-defined SDP gap.

## 3. Quantified diagonal and singular-value diffuseness

Write
\[
 \eta={\tau(B)\over n^{3/2}}\ge1,\qquad
 x_i=d_i/\sqrt n,\quad y_j=e_j/\sqrt n,\quad
 r_j=\sigma_j/\sqrt n.
\]
Let E denote the uniform average within the indicated n-coordinate
list. From (3)--(5), both x and y satisfy
\[
 E x=\eta,\quad E(1/x)\le\eta,\quad x_i\ge1/\eta,
\]
and hence
\[
 \boxed{E{(x-1)^2\over x}\le2(\eta-1),\qquad
 E|x-1|,\ E|x^{-1}-1|\le\sqrt{2\eta(\eta-1)}.}           \tag{9}
\]
The same bounds hold for y. The first identity uses
\((x-1)^2/x=x+x^{-1}-2\); the last two follow by Cauchy--Schwarz,
once with the weight x and once with the weight 1/x.

The singular values satisfy
\[
                       E r^2=1,\qquad E r^3\le\eta.
\]
Cauchy--Schwarz in \(E r^2=E(r^{1/2}r^{3/2})\) gives
\(E r\ge1/\eta\). Therefore
\[
 \boxed{E(r-1)^2\le {2(\eta-1)\over\eta},\qquad
 E[r^2\,1_{\{r>R\}}]\le{\eta\over R}\quad(R>0).}         \tag{10}
\]
Also
\[
                    {\operatorname{rank}B\over n}
                                     \ge\eta^{-2}.     \tag{11}
\]
Indeed Holder on the nonzero singular values gives
\(1=E r^2\le(E r^3)^{2/3}\Pr(r>0)^{1/3}\).

Thus an SDP value close to its minimum forces both singular-value
flatness in normalized squared distance and optimal-diagonal flatness
in the weighted and inverse senses (9). The latter includes the
uniform inverse bound \(\|D^{-1}\|_{\rm op}\le\eta/\sqrt n\).
A bounded eta controls total high-singular-value squared mass by (10);
it does not forbid a few exceptional singular values.

## 4. Complementarity residuals with exact constants

Put \(Q=D-H_B\succeq0\). Equations (3), (6), and (8) give
\[
                   \operatorname{tr}(Z^TQZ)=2g.        \tag{12}
\]
Since \(D+H_B\succeq0\), one also has \(0\preceq Q\preceq2D\).
Equivalently \(0\preceq D^{-1/2}QD^{-1/2}\preceq2I\), and squaring
that PSD contraction proves \(QD^{-1}Q\preceq2Q\). Consequently
\[
 \boxed{\|D^{-1/2}QZ\|_F^2\le4g.}                       \tag{13}
\]
No maximum-diagonal or operator cap is used in (13).

For a second quantitative statement let
\(d_{\max}=\max(\max_i d_i,\max_j e_j)\) and \(L_B=\|B\|_{\rm op}\).
Define
\[
 R_1=D_rB-B|B|,\qquad R_2=D_c|B|-|B|^2.
\]
The identity \(QZ=n^{-1/2}(R_1,R_2)^T\) and
\(\|Q\|_{\rm op}\le d_{\max}+L_B\) give
\[
                  \|R_1\|_F^2+\|R_2\|_F^2
                     \le2n(d_{\max}+L_B)g.              \tag{14}
\]
Choose an orthogonal polar factor U with \(B=U|B|\); a completed
singular-value decomposition supplies such U even when B is singular.
Then
\[
 U R_2^T=B D_c-B|B|,\qquad
                  D_rB-BD_c=R_1-U R_2^T.
\]
Since every \(B_{ij}^2=1\), (3) implies
\[
 \begin{aligned}
 \|D_rB-BD_c\|_F^2
 &=\sum_{ij}(d_i-e_j)^2\\
 &=n\left[\sum_i(d_i-\tau/n)^2+\sum_j(e_j-\tau/n)^2\right].
 \end{aligned}
\]
Combining with (14) proves
\[
 \boxed{\sum_i(d_i-\tau/n)^2+\sum_j(e_j-\tau/n)^2
                         \le4(d_{\max}+L_B)g.}          \tag{15}
\]
Unlike (13), this unweighted consequence explicitly retains its
maximum-diagonal/operator factor. It is not advertised as cap-free.

## 5. Exact zero-gap classification, including deficient rank

The following two conditions are equivalent:

- \(g=0\), so the actual canonical primal (6) is SDP-optimal;
- all nonzero singular values of B equal one scalar ell.

When these conditions hold, every optimal diagonal in (2) equals
\(\ell I_{2n}\). A uniform optimal diagonal alone is not claimed
to imply zero canonical gap.

Here is the rank-deficient necessity in full. From g=0 and (12),
positivity implies QZ=0. Thus
\[
               D_rB=B|B|,\qquad D_c|B|=|B|^2.           \tag{16}
\]
Taking the transpose of the second equation gives
\(|B|D_c=|B|^2\). Therefore D_c commutes with |B| and preserves
its range and kernel. On the range of |B|, its restriction agrees
with |B|, since |B| is invertible there. B annihilates the kernel.
It follows that \(BD_c=B|B|\), without inverting |B| on its kernel.
Equivalently this follows immediately from the polar identity
\(U R_2^T=BD_c-B|B|\) with R_2=0.

Now \(D_rB=BD_c\). Every entry B_ij is nonzero, so
\(d_i B_{ij}=B_{ij}e_j\) forces all d_i and e_j to have the same
value ell. Equation (16) becomes \(|B|^2=\ell|B|\), proving that
its positive eigenvalues are all ell. This applies to every optimal D.

Conversely, if the positive singular values equal ell, then
\(D=\ell I_{2n}\) is feasible. If r=rank B,
\[
 r\ell^2=n^2,\qquad
 {1\over n}\sum_j\sigma_j^3={r\ell^3\over n}=n\ell.
\]
The feasible diagonal gives \(\tau\le n\ell\), while (7) gives
the reverse inequality. Thus g=0 and the preceding argument
also proves uniqueness of the optimal diagonal.

In this class
\[
 \ell={n\over\sqrt r},\qquad
 \eta=\sqrt{n/r}.
\]
In particular \(\tau=n^{3/2}\) is equivalent to full rank with
\(B^TB=nI\): the complete sign matrix is a real Hadamard matrix.
Zero canonical gap is more general; rank-deficient scaled partial
isometries are not silently excluded.

## 6. Boolean versus SDP values and the spectral rounding constraint

The quantities beta and tau are distinct. Certainly beta<=tau.
For completeness the elementary tensor-rounding argument gives
\[
 \tau\le K_G\beta,\qquad
                   K_G={\pi\over2\log(1+\sqrt2)}.        \tag{17}
\]
One direct proof lifts unit vectors into odd tensor powers with
coefficients \(\sqrt{c^{2k+1}/(2k+1)!}\), using alternating signs
in the second family, where \(c=\log(1+\sqrt2)\) and sinh c=1.
The lifted vectors are unit and their cross inner products are
\(\sin(c\langle u_i,v_j\rangle)\). Gaussian sign rounding gives
expected product exactly \(2c\langle u_i,v_j\rangle/\pi\).
The Gram matrix has a finite-dimensional realization, so no
infinite-dimensional Gaussian is invoked. Taking expectations
proves (17).

There is also the more spectral Boolean rounding bound. Put
\(\kappa=2/\pi\), distinct from the tensor constant c above. The
unit row families \(B/\sqrt n\) and \(|B|/\sqrt n\) have cross
correlation matrix \(W=B|B|/n\), whose entries lie in [-1,1].
Gaussian sign rounding and the positive-coefficient arcsine series give
\[
 |\arcsin w-w|\le(\pi/2-1)w^2,\qquad |w|\le1.
\]
The linear term is \(\operatorname{tr}|B|^3/n\), and
\(\sum W_{ij}^2=\operatorname{tr}|B|^4/n^2\). Hence
\[
 \boxed{\beta(B)\ge{\kappa\over n}\sum_j\sigma_j^3
                -{1-\kappa\over n^2}\sum_j\sigma_j^4.}  \tag{18}
\]
The errors retain their signs and the full nonnegative moment term.

If beta(B)<=C n^(3/2), the elementary L_B<=n and
\(\sum\sigma^4\le L_B\sum\sigma^3\) first imply
\[
 \sum\sigma^3\le{\beta n\over2\kappa-1}=O_C(n^{5/2}).
\]
It follows that \(L_B=O_C(n^{5/6})=o(n)\); substituting this
back into (18) improves the coefficient to
\[
 \beta(B)\ge(\kappa-o_C(1)){1\over n}\sum_j\sigma_j^3.     \tag{19}
\]
Thus the cube bound supplies a uniform third-moment/tail budget,
while the actual SDP gap and complementarity supply the additional
diagonal structure in (13)--(16). These are different pieces of data.

For an actual conditional cross minimizer, \(\beta(B)\le F_A^*\).
Any available O(n^(3/2)) conditional cap can therefore be inserted
into (17)--(19). This does NOT identify tau with beta, force
\(\eta-1=o(1)\), or force \(g=o(n^{3/2})\). In particular no
near-minimal SDP hypothesis is imported from original norm optimality.

The proved relations constrain actual candidate diagonals and spectral
masses for subsequent resolvent analysis. Evaluating the resulting
Gaussian upper sharply on all attainable optimizer shells remains open.


