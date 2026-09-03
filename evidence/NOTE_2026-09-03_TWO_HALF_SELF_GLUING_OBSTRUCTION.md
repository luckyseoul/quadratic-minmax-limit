# Two-half phase self-gluing: an exact obstruction

**Status:** proved obstruction to one specific recursive construction. This
does not prove the MathOverflow limit or rule out other gluings.

Let

\[
 \mathcal Q=\{(\pm1\pm i)/\sqrt2\},\qquad
 \nu_4(H)=\max_{w\in\mathcal Q^n}|w^*Hw|
\]

for a zero-diagonal Hermitian matrix $H$.  For $H=A+iR$, where $A$
is a symmetric signing and $R$ is a skew signing, the fourth-phase
identity in Proposition 6.5 is

\[
 \nu_4(H)=B(A,R)
 =\max_{x,y}\bigl(|Q_A(x)+Q_A(y)|+|x^TRy|\bigr).       \tag{1}
\]

## Exact rank-one amplification

Put

\[
 K=\begin{pmatrix}1&i\\-i&1\end{pmatrix}.
\]

Then every zero-diagonal Hermitian $H$ satisfies

\[
 \boxed{\nu_4(K\otimes H)=4\nu_4(H).}                  \tag{2}
\]

Indeed, for $u,v\in\mathcal Q^n$,

\[
 (u,v)^*(K\otimes H)(u,v)=(u+iv)^*H(u+iv).
\]

Coordinatewise, $u_j+iv_j$ lies in

\[
 Z=\mathcal Q+i\mathcal Q
   =\{a+ib:a,b\in\{-\sqrt2,0,\sqrt2\}\}.
\]

The convex hull of $Z$ is the square with vertices $2\mathcal Q$.
Because $H_{jj}=0$, after the other coordinates are fixed,
$z^*Hz$ is affine in the two real coordinates of $z_j$; its absolute
value is convex.  Maximizing one coordinate at a time therefore pushes
every coordinate to a vertex of that square.  Hence

\[
 \max_{z\in Z^n}|z^*Hz|
 \le \max_{z\in(2\mathcal Q)^n}|z^*Hz|=4\nu_4(H).
\]

Every vertex $2\mathcal Q$ belongs to $Z$, proving the reverse
inequality and (2).

## The natural second lift misses the required factor

For arbitrary diagonal sign matrices $D,E$, define

\[
 S_D=\begin{pmatrix}A&D-R\\D+R&A\end{pmatrix},\qquad
 U_E=\begin{pmatrix}R&A+E\\-A-E&R\end{pmatrix}.
\]

These are respectively a symmetric complete signing and a skew complete
signing of order $2n$, and

\[
 S_D+iU_E
 =K\otimes(A+iR)+
   \begin{pmatrix}0&D+iE\\D-iE&0\end{pmatrix}.          \tag{3}
\]

For a fourth-phase vector $(u,v)$, the correction in (3) is

\[
 2\operatorname{Re}\sum_j\overline{u_j}(d_j+ie_j)v_j.
\]

Here $\overline{u_j}v_j\in\{\pm1,\pm i\}$, so every summand has absolute
value exactly $2$.  The sharp perturbation bound is therefore

\[
 \boxed{|B(S_D,U_E)-4B(A,R)|\le2n.}                    \tag{4}
\]

The analogous real paired-edge calculation gives

\[
 \boxed{|\Phi(S_D)-B(A,R)|\le n.}                     \tag{5}
\]

Combining (4)--(5),

\[
 B(S_D,U_E)-2\sqrt2\,\Phi(S_D)
 \ge (4-2\sqrt2)B(A,R)-(2+2\sqrt2)n.                 \tag{6}
\]

Since $B(A,R)\ge2\Phi(A)$ and
$\Phi(A)\ge n\sqrt{n-1}/\pi$, the right side is at least

\[
 {2(4-2\sqrt2)\over\pi}n\sqrt{n-1}-(2+2\sqrt2)n.    \tag{7}
\]

Thus the defect is positive of order $n^{3/2}$.  No normalized vanishing
error, and in particular no dyadic-Dini error, absorbs it.  Reusing the real
and imaginary halves through this phase-compatible rank-one self-gluing
cannot bootstrap the multiplier-two bound.  The failure is the exact
$4$ versus $2\sqrt2$ amplification, not paired-edge signs or error
bookkeeping.

## What remains under a general bisection

For a vertex split $L\sqcup R$, write

\[
 A=\begin{pmatrix}A_L&C\\C^T&A_R\end{pmatrix},\qquad
 S=\begin{pmatrix}S_L&T\\-T^T&S_R\end{pmatrix}.
\]

If $U=P\sqcup Q$, the outward half-cut consists of the two diagonal
half-cuts together with the cross edges

\[
 \mathcal F_T(P,Q)=
 \{ij:T_{ij}=1, i\in P, j\notin Q\}
 \cup
 \{ij:T_{ij}=-1, i\notin P, j\in Q\}.              \tag{8}
\]

Consequently the exact bisection recurrence is

\[
 \Lambda_A(S)=\max_{P,Q}
 \Phi\!\begin{pmatrix}
 A_L^{F_{S_L}(P)}&C^{\mathcal F_T(P,Q)}\\
 (C^{\mathcal F_T(P,Q)})^T&A_R^{F_{S_R}(Q)}
 \end{pmatrix}.                                      \tag{9}
\]

Equation (9) is a decomposition, not yet a recurrence estimate: it leaves a
uniform family of coupled cross-rectangle flips.  A viable recursive proof
needs a new theorem coupling those rectangular terms to the two diagonal
energies.  Bounding them independently discards the endpoint payment from
the two-half diamond, while the phase-compatible coupling above is excluded
by (2)--(7).

**Duplication guard.** Do not retry $K\otimes(A+iR)$, the matrices
$S_D,U_E$ above, or an equivalent rank-one phase-compatible reuse of the
same two halves as a multiplier-two bootstrap.  A changed premise must alter
the cross-block coupling, not just the paired-edge diagonals.
