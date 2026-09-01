# Tetrahedral tripling frame

**Status:** Proposition 6.7 is a proved exact reduction. The multiplier-three
ray, and therefore the original limit problem, remain **OPEN**.

## Exact reduction

For an order-\(n\) signing \(A\), three skew signings \(P,Q,T\), and diagonal
sign matrices \(D_{12},D_{13},D_{23}\), Proposition 6.7 constructs

\[
{\cal S}=\begin{pmatrix}
A&P+D_{12}&Q+D_{13}\\
-P+D_{12}&A&T+D_{23}\\
-Q+D_{13}&-T+D_{23}&-A
\end{pmatrix}.
\]

All four projective constant states of each three-vertex cloud induce exactly
\(A\). Conversely, a \(3\times3\) sign block constant on those four states
must have skew off-diagonal part and diagonal consisting of two copies of the
endpoint sign and one negative copy. Thus the construction removes the
tripling endpoint-selection problem rather than assuming it away.

Writing

\[
I=Q_A(x)+Q_A(y)-Q_A(z),\quad
b=x^TPy,\quad c=x^TQz,\quad d=y^TTz,
\]

the cross-layer maximum is exactly

\[
K_3=\max_{x,y,z}\max\{|I+d|+|b+c|,\ |I-d|+|b-c|\}.
\]

The only omitted terms are the three internal edges in each cloud:

\[
\Delta=x^TD_{12}y+x^TD_{13}z+y^TD_{23}z,\qquad |\Delta|\le3n.
\]

Consequently \(|\Phi({\cal S})-K_3|\le3n\). The sharp remaining target is

\[
K_3\le3^{3/2}m_n+n^{3/2}\Omega(n)
\]

with a vanishing dyadic-Dini tail. This would give
\(H(3n)\le3H(n)+O(n\Omega(n)+\sqrt n)\).

## Single-\(R\) specialization and proved shields

For \(P=R,Q=-R,T=R\),

\[
K_3=\max_{x,y,z}\bigl(
|Q_A(x)+Q_A(y)-Q_A(z)|+
|x^TRy+y^TRz+z^TRx|\bigr).
\]

The cyclic term has the cancellation identity

\[
C_R=(x-y)^TR(y-z).
\]

If \(\Lambda=\|R\|_{\rm op}\), this proves the distance-product shield

\[
|C_R|\le4\Lambda\min_{\rm cyc}
\sqrt{d_H(x,y)d_H(y,z)}.
\]

The exact tetrahedral partition in (6.35), together with
\(w=u_1+e^{2\pi i/3}u_2+e^{4\pi i/3}u_3\), proves the distinguished-endpoint
shield

\[
|C_R|\le\frac4{\sqrt3}\Lambda
\bigl(|\operatorname{supp}u_1|+
|\operatorname{supp}u_2|+|\operatorname{supp}u_3|\bigr).
\]

Either shield closes a triple when its right side is at most
\(3\sqrt3m_n-|I|\). No claim is made for the complement. That complement is
the open tetrahedral diamond. The three-independent-signing frame is exact,
but no additional coupled shield for it is claimed here.

## Ruled-out immediate shortcuts

The multiplier-two diamond does not algebraically imply the tetrahedral
diamond, even if all three skew blocks satisfy the pairwise bound. At order
two, take

\[
A=\begin{pmatrix}0&-1\\-1&0\end{pmatrix},\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Here \(M=1\) and the pair score is \(2<2\sqrt2M\). Nevertheless every one
of the eight choices \(P,Q,T\in\{R,-R\}\) has \(K_3=7>3\sqrt3M\). Thus
closing Proposition 6.5 three times cannot by itself close Proposition 6.7;
a genuinely three-state cancellation or budget is required.

Independent random skew blocks followed by a statewise union bound are also
insufficient. For a fixed triple, each cross-sign combination has the exact
law

\[
2\sum_{j=1}^{N}\epsilon_j,\qquad
N=\sum_{uv\in\{xy,xz,yz\}}d_H(u,v)(n-d_H(u,v)).
\]

The uniform tetrahedral Hamming profile simultaneously has
\(N=(3/4+o(1))n^2\) and \(\exp((3\log2+o(1))n)\) states. More strongly,
every signing has at least \((1/2+o(1))2^n\) states with \(|Q_A(x)|\le n\),
by the exact second moment. Typical triples from that actual low-energy set
have the uniform profile. Exact moderate-deviation tails therefore make the
sum of the individual failure probabilities at
\(L=\ell n^{3/2}\) at least

\[
\exp\left((3\log2-\ell^2/6-o(1))n\right).
\]

A literal first-moment certificate needs
\(\ell\ge\sqrt{18\log2}=3.532230\ldots\), whereas the desired threshold is
at most \((3\sqrt3/2+o(1))n^{3/2}=2.598076\ldots n^{3/2}\). This rules out
only the per-state union bound, not correlated random blocks analyzed as one
process. A viable probabilistic attack must use those dependencies, an
\(A\)-dependent construction, or a genuine theorem on energy level sets.

## Verification

The deterministic test **tests/test_prop67_tetrahedral_tripling.py** exhausts
all sign blocks in the endpoint converse and a fixed order-three frame. It
checks the four equal endpoints, the four-sign collapse, the exact
single-skew diamond, the internal-edge formula, the \(3n\) comparison, the
cyclic cancellation, both identities in the tetrahedral partition, and the
pairwise-shortcut counterexample above.
