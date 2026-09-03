# Top-degree QR pair and antipodal-fibre audit

**Date:** 2026-09-02
**Status:** proved common-graph inequality and proved method barrier; residual
(ii) remains open

Work at the lower endpoint of branch B,

\[
 p=4r+1,\qquad q={p+1\over2}=2r+1,\qquad
 T=\sum_U\epsilon_U P_U=p+4,
\]

with the common pencil point translated to the origin. There are \(q\)
hard and \(q\) opposite directions. The hard parallel counts are \(r+2\),
except for three counts \(r+3\), and every opposite count is \(r\).

For an edge \(e=\{a-d,a+d\}\), put

\[
 \tau_e=\chi(d)=\epsilon_{\langle d\rangle}.
\]

In direction \(L\), the normalized off-diagonal coefficient contributed by
this edge is

\[
 W_L(s,t)=\epsilon_L\tau_e,
 \qquad s=L(a+d),\quad t=L(a-d),
\]

provided \(L(d)\ne0\). Define

\[
 R_L=\sum_{s<t,\ s+t\ne0}W_L(s,t),\qquad
 M_L=\sum_{s<t,\ st\ne0}W_L(s,t)\eta(st).
\tag{1}
\]

The outer \(\epsilon_L\) in the projective moment equation is essential.

## The top-degree identity and its exact quotient

Put \(m_0=(p-1)/2\), \(X=(s+t)^2\), and \(Y=st\). In
\(\mathbf F_p[s,t]\),

\[
\begin{aligned}
 F(s,t)
  &:=(s+t)^{p-1}-(st)^{m_0}\\
  &=\sum_{k=0}^{m_0-1}4^k
    (s-t)^2(st)^k(s+t)^{p-3-2k}.
\end{aligned}
\tag{2}
\]

Indeed, the right side telescopes to \(X^{m_0}-4^{m_0}Y^{m_0}\),
and \(4^{m_0}=2^{p-1}=1\). On an off-diagonal fibre pair, \(F\) is the
indicator of \(s+t\ne0\) minus \(\eta(st)\). Proposition 15.759 therefore
gives

\[
 Q-Q_M\equiv0\pmod p,
 \qquad
 Q=\sum_L\epsilon_LR_L,\quad Q_M=\sum_L\epsilon_LM_L.
\tag{3}
\]

There is a useful exact description of the integer quotient

\[
 \mathcal D={Q-Q_M\over p}.
\tag{4}
\]

For a fixed edge, if \(a,d\) are independent, the ratios
\(L(a)/L(d)\) run once through \(\mathbf F_p\) over the \(p\) transverse
directions. Hence

\[
 \sum_x\bigl({\bf1}_{x\ne0}-\eta(x^2-1)\bigr)=p.
\]

If \(a=0\), all \(p\) ratios are zero and the sum is \(-p\). If
\(a=cd\ne0\), all ratios equal \(c\), and the sum is
\(p(1-\eta(c^2-1))\). Consequently the contribution of \(e\) to
\(\mathcal D\) is exactly

\[
 \boxed{
 \tau_e\begin{cases}
  1,&a,d\text{ independent},\\
  -1,&a=0,\\
  1-\eta(c^2-1),&a=cd\ne0.
 \end{cases}}
\tag{5}
\]

The last value lies in \(\{0,1,2\}\). The factor \(\tau_e\) in (5) must
not be dropped: the source coordinate in Proposition 15.759 is
\(\tau_e\mathbf1_{e\in H}\), not the unsigned edge indicator. Thus
\(\mathcal D\) is not a nonnegative count.

## Exact QR/NQR cell normalization

For the two profiles centred at zero, write

\[
 b_L=B_L(\mathrm{QR})+B_L(\mathrm{NQR}),\qquad
 c_L=C_L(\mathrm{QR})+C_L(\mathrm{NQR}).
\]

There are \(q-2\) retained hard pencil fibres, whose exceptional label is
\(j=0\), and two outlier hard fibres, whose exceptional label is nonzero.
Since both QR profiles contain zero, the literal \(x_j\) has paired value
two in a retained direction and paired value one in an outlier direction.
Using

\[
 A_L=x_j+2B_L\quad\text{(hard)},\qquad
 A_L=2C_L\quad\text{(opposite)},
\]

in \(\epsilon_LS_H=3+2A_L=P_L+\sum W_Lz_sz_t\) gives

\[
 M_L=\begin{cases}
  5-P_L+2b_L,&L\text{ retained hard},\\
  4-P_L+2b_L,&L\text{ outlier hard},\\
  3-P_L+2c_L,&L\text{ opposite}.
 \end{cases}
\tag{6}
\]

Since the hard-minus-opposite parallel total is \(T\), (6) yields the
particularly simple identity

\[
 \boxed{Q_M=-5+2(B-C)},\qquad
 B=\sum_{L\ {\rm hard}}b_L,\quad
 C=\sum_{L\ {\rm opposite}}c_L.
\tag{7}
\]

This repairs the erroneous formula obtained by treating all hard
exceptional labels as nonzero.

For completeness, the \(2p\) additive translates of QR and NQR form a

\[
 2\text{-}\left(p,{p+1\over2},{p+1\over2}\right)
\]

design: every point occurs \(p+1\) times and every pair occurs
\((p+1)/2\) times. Thus a nonnegative quadratic of scaled mass
\(4p\mathbb E B=M\) satisfies

\[
 \sum_{k\in\mathbf F_p}
   \bigl(B(\mathrm{QR}+k)+B(\mathrm{NQR}+k)\bigr)={M\over2}.
\tag{8}
\]

At the lower endpoint, summing (8) over all directions gives

\[
 \sum_{L,k}b_L(k)={t(p+1)\over2},\qquad
 \sum_{L,k}c_L(k)={q(r-1)(p-3)\over2},
\]

and hence

\[
 \sum_{L,k}\bigl(b_L(k)-c_L(k)\bigr)=-q^2.
\tag{9}
\]

Equivalently, if \(n(x)\) is the number of exceptional hard fibres through
a common affine centre \(x\), then

\[
 Q_M(x)=-2r-4+n(x)+2(B(x)-C(x)),
 \qquad \sum_xQ_M(x)=-p^2T.
\tag{10}
\]

Here \(n(0)=q-2\), recovering (7). Equations (8)--(10) are aggregate
identities; they do not bound the value at the pencil centre because the
nonnegative masses in (8) need not be equidistributed.

## A positive lower bound for the antipodal-fibre term

For a radial direction \(U\), let \(E_U\) be the number of selected edges
whose two endpoints lie on the line \(U\) through the origin, and let
\(A_U\) be the number of selected antipodal edges \(\{z,-z\}\subset U\).
An edge is counted by \(R_L\) in \(p-1\) directions when \(a,d\) are
independent, in \(p\) directions when \(a\parallel d\ne0\), and in no
direction when \(a=0\). Therefore

\[
 \boxed{
 Q=(p-1)T+\sum_U\epsilon_UE_U-p\sum_U\epsilon_UA_U.}
\tag{11}
\]

Let \(A_h\) be the number of hard antipodal edges. On either removed hard
ray there are at most \(P_U\le r+3\) of them. On a retained hard ray, its
boundary points are the intersections with the two outlier affine lines.
They form an antipodal pair only when one is the unique point of
\(\ell_1\cap(-\ell_2)\). Hence at most one retained hard ray contains a
boundary-boundary antipodal pair.

The remaining hard antipodal edges touch \(V\setminus D\). Such edges form
a matching. If \(k_1\) of them have one endpoint outside \(D\) and \(k_2\)
have two, their \(k_1+2k_2\) distinct outside endpoints all have positive
even degree. Since

\[
 \sum_{v\notin D}\deg_H(v)
 \le\sum_v(\deg_H(v)-\mathbf1_D(v))=2s,
\]

we have \(k_1+k_2\le s\). It follows that

\[
 A_h\le2r+7+s.                                      \tag{12}
\]

In a hard direction, \(E_U-pA_U\ge-(p-1)A_U\). In an opposite direction,
\(-E_U+pA_U\ge-P_U\). Since the opposite parallel total is
\(q r=2r^2+r\), (11)--(12) give

\[
 Q\ge(p-1)(p+4)-(p-1)(2r+7+s)-(2r^2+r)
   =r(6r-9-4s).                                    \tag{13}
\]

Thus

\[
 \boxed{
 Q\ge r(6r-37)>0\quad(s=7),\qquad
 Q\ge r(6r-45)>0\quad(s=9)}                       \tag{14}
\]

for every live \(r\ge13\). This is a genuine common-graph inequality,
not a separate-direction capacity estimate.

## Why positivity does not close the QR pair

The quotient (5) can change sign while preserving exactly the data used in
(14). Choose an opposite radial line \(U\). Away from the at most four
points toggled by the two outlier lines and their negatives, choose
\(\alpha,\beta\in\mathbf F_p^*\) with \(\eta(\alpha\beta)=-1\). All four
points \(\pm\alpha u,\pm\beta u\) lie in \(D\). Compare the two matchings

\[
 \{\{\alpha u,-\alpha u\},\{\beta u,-\beta u\}\}
 \quad\text{and}\quad
 \{\{\alpha u,\beta u\},\{-\alpha u,-\beta u\}\}.
\tag{15}
\]

They have the same four-vertex boundary, two edges in the same opposite
direction, the same signed total, and the same boundary phase. Their
symmetric difference is a four-cycle, so this is an exact
boundary-, edge-count-, and parallel-quota-preserving switch.

For the antipodal matching in (15), (5) gives

\[
 (\mathcal D,Q,Q_M)=(2,0,-2p).
\]

For the crossed matching,

\[
 \eta\!\left(\left({\alpha+\beta\over\alpha-\beta}\right)^2-1\right)
 =\eta(\alpha\beta)=-1,
\]

so (5) gives

\[
 (\mathcal D,Q,Q_M)=(-4,-2p,2p).                  \tag{16}
\]

The switch changes \(\mathcal D\) by \(-6\), \(Q\) by \(-2p\), and
\(Q_M\) by \(4p\), while preserving every coarse near-pencil datum above.
In the cell identity (7), this is a change of \(2p\) in \(B-C\).

Combining (3), (7), and (14) therefore does **not** contradict
\(B\ge0\):

\[
 Q=-5+2(B-C)+p\mathcal D,
\tag{17}
\]

and both \(C\) and the signed quotient \(\mathcal D\) occur with
uncontrolled signs. The top-degree QR/NQR relation and the
antipodal-fibre count are a proved method barrier, not a closure of
residual (ii). A successful continuation must constrain the signed
four-cycle freedom in (15), or use a different common-graph condition that
does so.
