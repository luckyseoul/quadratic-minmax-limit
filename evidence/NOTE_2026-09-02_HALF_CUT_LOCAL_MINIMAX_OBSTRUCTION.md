# Directed half cuts: exact one-edge optimality and the fractional obstruction

**Status:** proved an exact nearest-neighbor optimality certificate and an
exact plateau obstruction.  The certificate contains no numerical upper
bound and can become tautological, so it does not establish the
`sqrt(2) Phi(A)` bound.  The ordinary convex relaxation has optimum exactly
`Phi(A)` and therefore deletes the integral orientation problem.  This note
does not prove the multiplier-two diamond and uses no finite census.

Let `A` be a symmetric signing, let `S` be a tournament, and let

\[
 H_U(S)=A^{F_S(U)},\qquad
 \Lambda_A(S)=\max_{U\subseteq[n]}\Phi(H_U(S)),                 \tag{1}
\]

where `F_S(U)` consists of the cut edges directed by `S` from `U` to its
complement.  The directed half-cut identity gives

\[
 \Lambda_A(S)={1\over2}B(A,A\circ S).                          \tag{2}
\]

Thus the desired integral statement is

\[
 \min_S\Lambda_A(S)
 \le\sqrt2\,M+o_{\rm Dini}(n^{3/2}),
 \qquad M=\Phi(A)=m_n.                                         \tag{3}
\]

## 1. Exact effect of reversing one arc

Fix an unordered edge `e={p,q}` and write `S^e` for the tournament obtained
by reversing that arc.  Directly from the definition of a directed half
cut,

\[
 H_U(S^e)=
 \begin{cases}
 H_U(S),&e\notin\delta(U),\\
 H_U(S)^e,&e\in\delta(U),
 \end{cases}                                                    \tag{4}
\]

where the superscript `e` on a symmetric signing means reversal of its
`e` coefficient.  Define

\[
\begin{aligned}
 P_e(S)&=\max_{U:e\notin\delta(U)}\Phi(H_U(S)),\\
 Q_e(S)&=\max_{U:e\in\delta(U)}\Phi(H_U(S)),\\
 \widehat Q_e(S)&=\max_{U:e\in\delta(U)}\Phi(H_U(S)^e).
\end{aligned}                                                    \tag{5}
\]

Then (4) proves the exact identities

\[
 \Lambda_A(S)=\max\{P_e(S),Q_e(S)\},\qquad
 \Lambda_A(S^e)=\max\{P_e(S),\widehat Q_e(S)\}.                \tag{6}
\]

Consequently `S` is a one-arc local minimum if and only if, for every edge
`e`,

\[
 \boxed{\quad
 P_e(S)=\Lambda_A(S)\quad\hbox{or}\quad
 \widehat Q_e(S)\ge\Lambda_A(S).
 \quad}                                                         \tag{7}
\]

This is the full consequence of comparison with all single-arc reversals;
there is no differentiability assumption or omitted tie case.

The second alternative has a useful witness form.  Put
`L=Lambda_A(S)>2`.  If `P_e(S)<L`, local minimality supplies a cut `U`
separating `e` and a Boolean state `y` such that, with `H=H_U(S)` and
`h=H_pq`,

\[
 |Q_{H^e}(y)|=|Q_H(y)-2h y_p y_q|\ge L.                         \tag{8}
\]

Since `|Q_H(y)|<=L`, (8) forces

\[
 |Q_H(y)|\ge L-2,
 \qquad Q_H(y)h y_p y_q<0.                                    \tag{9}
\]

Thus a genuinely constraining edge needs a near-maximizer of a crossing
half-cut neighbor on which that edge is anti-aligned with the total energy.
This is an additive-two, near-active witness condition; it is not an upper
bound on `L/M`.

## 2. The plateau obstruction

Let

\[
 \mathcal U_*(S)=
 \{U:\Phi(H_U(S))=\Lambda_A(S)\}                               \tag{10}
\]

be the active cut family.  If an active cut does not separate `e`, its term
in (1) is unchanged when `e` is reversed.  Hence

\[
 e\notin\bigcap_{U\in\mathcal U_*(S)}\delta(U)
 \quad\Longrightarrow\quad P_e(S)=\Lambda_A(S).                \tag{11}
\]

In particular,

\[
 \boxed{\quad
 \bigcap_{U\in\mathcal U_*(S)}\delta(U)=\varnothing
 \quad\Longrightarrow\quad
 S\text{ is automatically a one-arc local minimum.}
 \quad}                                                         \tag{12}
\]

Only edges simultaneously crossing every active cut can impose the witness
condition (9).  Therefore single-edge descent can become completely
tautological merely because several tied constraints have different cut
supports.  It need not express balance of the directed part of any active
cut.

For a globally optimal coefficient signing `A`, minimality adds only

\[
 \Phi(H_U(S))\ge M                                               \tag{13}
\]

for every `U,S`, since every `H_U(S)` is another coefficient signing.  This
lower bound neither rules out the plateau (12) nor turns (7)--(9) into the
upper bound (3).  A proof based on local orientation descent would need an
additional theorem controlling the geometry and common crossing set of the
active family.

## 3. The convex minimax relaxation is exactly trivial

This obstruction is not repaired by relaxing the skew signs.  Let `R` range
over real skew matrices with `|R_ij|<=1`, and use the exact two-half form

\[
 \Lambda_A^{\rm frac}(R)
 =\max_{U,y}\left(
 |I_A(U,y)|+
 \left|\sum_{i\in U,j\notin U}R_{ij}y_i y_j\right|
 \right),                                                       \tag{14}
\]

where `I_A(U,y)` is the sum of the two within-part `A` energies.  If `x` is
obtained from `y` by switching `U`, then

\[
 I_A(U,y)={Q_A(x)+Q_A(y)\over2},
\]

so `|I_A(U,y)|<=M`.  Equality is attained after taking `U` empty and `y` a
maximizer of `A`.  It follows that

\[
 \boxed{\quad
 \min_{|R_{ij}|\le1}\Lambda_A^{\rm frac}(R)
 =\Lambda_A^{\rm frac}(0)=M.
 \quad}                                                         \tag{15}
\]

The lower bound in (15) again comes from the empty cut and holds for every
fractional `R`; the upper bound uses `R=0`.  Equivalently, the linear-program
dual can be supported entirely on an endpoint maximizer and contains no
orientation information.

Thus convex minimax calculus reaches the center `R=0`, while (3) asks for a
vertex of the sign cube.  Passing from (15) to (3) is an integral rounding
theorem with factor `sqrt(2)` and Dini-scale error; it is not a consequence
of convexity.  The independent random-rounding attempt and its central-cut
entropy failure are already recorded separately, so (15) does not reopen
that route.

## 4. Verdict

The strongest exact one-edge conclusion is the disjunction (7), sharpened
to the near-maximizer condition (9) only on the common crossing set of all
active cuts.  Outside that set local optimality is a plateau tautology.  The
fractional optimum (15) is exactly `M`, but it supplies no integral
tournament.  Neither fact proves
`min_S Lambda_A(S)<=sqrt(2)M`; a successful variational proof must add a
global multi-edge rounding or active-family structure theorem not contained
in coefficient minimality or one-edge local descent.

The multiplier-two ray and the original MO limit remain open.
