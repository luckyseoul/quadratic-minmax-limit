# Exact signed-automorphism identity for the sampled cross-block winner

## Scope and finite provenance

The bounded order-six finite-step experiment stored in
`original_mo_finite_cross_gpu_result.json` has SHA-256
`15e5c0ee119fcb90dec32b584217c6c8998cde67e5e43e347138926bf5d455a8`.
It sampled 8,192 Gaussian-sign and 8,192 independent cross blocks, with
two coherent references, on 20 prescribed temperature/step profiles.
The Gaussian sample's reported winning index is 1067 on every profile.

This note identifies that explicit winner algebraically. It proves
equality with the coherent block \(A-I\) for all coupling parameters,
not only numerical agreement on the finite grid. It does not assert
that the sample minimum is the exhaustive cross-block minimum, classify
global host minima, or draw an asymptotic conclusion.

## 1. Explicit host, winner, and signed permutation

The fixed host and sampled winning cross block are
\[
A=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&1&-1&-1\\
1&1&0&-1&1&-1\\
1&1&-1&0&-1&1\\
1&-1&1&-1&0&1\\
1&-1&-1&1&1&0
\end{pmatrix},\qquad
B=\begin{pmatrix}
1&-1&1&1&-1&-1\\
1&1&-1&-1&-1&1\\
-1&1&-1&1&1&-1\\
-1&1&1&1&1&1\\
-1&1&1&-1&-1&1\\
1&1&-1&-1&1&-1
\end{pmatrix}.
\]
Indices in the following formulas are one-based. Define the signed
permutation matrix by its columns:
\[
P=(e_4,e_1,e_6,e_2,-e_3,-e_5).
\]
Thus \(Pe_j=\varepsilon_j e_{\pi(j)}\), where
\[
\pi=(4,1,6,2,3,5),\qquad
\varepsilon=(1,1,1,1,-1,-1),
\]
and, explicitly,
\[
Pz=(z_2,z_4,-z_5,z_1,-z_6,z_3).
\]
It is orthogonal and bijects the Boolean cube. Let \(C=A-I\).
Comparing the six displayed columns gives
\[
(B_1,B_2,B_3,B_4,B_5,B_6)
=(C_4,C_1,C_6,C_2,-C_3,-C_5),
\qquad\boxed{B=CP=(A-I)P.}                                    \tag{1}
\]

For completeness, direct substitution into
\((P^TAP)_{jk}=\varepsilon_j\varepsilon_kA_{\pi(j),\pi(k)}\)
gives the following upper-triangular entries, grouped by row:
\[
(1,1,1,1,1;\ 1,1,-1,-1;\ -1,1,-1;\ -1,1;\ 1).
\]
These are exactly those of \(A\); the diagonal remains zero. Hence
\[
\boxed{P^TAP=A.}                                               \tag{2}
\]
The two permutation cycles are \((1\ 4\ 2)\) and \((3\ 6\ 5)\),
each with sign product one, so also \(P^3=I\).

## 2. Exact finite-pressure and Gibbs-moment equality

Write \(Q_A(x)=x^TAx/2\). Equation (2) implies
\(Q_A(Py)=Q_A(y)\). Under the cube bijection \(y'=Py\), equation
(1) gives
\[
x^TBy=x^TCy',\qquad Q_A(x)-Q_A(y)=Q_A(x)-Q_A(y').              \tag{3}
\]
Therefore, for every real \(\eta,\gamma\),
\[
\begin{split}
&\log\mathbb E_{x,y}\cosh\!\left[
 \eta(Q_A(x)-Q_A(y))+\gamma x^TBy\right]\\
&\hspace{12mm}=
\log\mathbb E_{x,y}\cosh\!\left[
 \eta(Q_A(x)-Q_A(y))+\gamma x^T(A-I)y\right].
\end{split}                                                    \tag{4}
\]
The expectations here are uniform. In particular, (4) holds at every
balanced choice \(\eta=(c/\sqrt6)\sqrt{1-t/2}\),
\(\gamma=(c/\sqrt6)\sqrt{t/2}\), for all \(c\ge0\) and
\(0\le t\le1\).

The same substitution preserves every actual Gibbs law
\(\mu_{A,s}(y)\propto e^{sQ_A(y)}\), for every real \(s\).
Thus, for independently sampled actual host spins at any two
temperatures, the entire joint distribution of the two internal
energies and the cross bilinear quantity agrees for \(B\) and \(C\).
In particular, if \(U,V\) are the actual opposite-temperature
covariance matrices, then
\[
\overline q(B):=\tfrac12\left[
 \operatorname{tr}(B^TUBV)+\operatorname{tr}(B^TVBU)\right]
=\overline q(A-I).                                             \tag{5}
\]
This follows either from (3) under the two phase orderings or directly
from their invariance under \(P\).

Consequently, the reported Gaussian winner is a known coherent block
in different signed coordinates. Its agreement with the coherent
reference is exact, and is not a new noncoherent selection mechanism.
The factorization specifically uses \(A-I\); no separate equivalence
with \(A+I\) is needed here. The previously proved conference
second-moment minimum is not rederived or promoted by this observation.
