# Threshold valley and active-witness geometry

Date: 2026-09-04.

Status: proved conditional reductions and a sufficient-method barrier. This
note does **not** prove the all-size minimal-four-gap shell bridge, close
residual (ii), or change a closure predicate. It makes no assertion that all
first-shell Boolean points have Boolean or affine eigenanchors.

## 1. Hypotheses, normalization, and existing inputs

Let $p\geq5$ be an odd prime, $n=p^2+1$, and $C$ the symmetric Paley
conference matrix. For a symmetric zero-diagonal real matrix $B$, set

\[
q_B(x)=\sum_{i<j}B_{ij}x_ix_j,
\qquad \Phi(B)=\max_{x\in\{\pm1\}^n}|q_B(x)|,
\qquad \Phi_0=\Phi(C)=pn/2.
\]

Write $A=C\oplus H=C-2(C\circ1_H)$, $h=|H|$, and suppose

\[
\Phi(A)=M=\Phi_0-4.
\tag{1}
\]

When invoking **threshold inclusion-minimality**, we mean the explicit
strong condition

\[
\Phi(C\oplus F)>M\quad\text{for every proper }F\subsetneq H.
\tag{2}
\]

The one-deletion identities used in Proposition 15.764 do not by themselves
assert (2). They suffice for the singleton conclusions below, but not for
the full all-subsets hierarchy.

A signed state is $z=(\varepsilon,[x])$, where

\[
\varepsilon\in\{\pm1\},\quad [x]=\{x,-x\},\quad
a_e(z)=\varepsilon C_{ij}x_ix_j\quad(e=\{i,j\}).
\]

Set

\[
\Delta_z=M-\varepsilon q_A(x)\in2\mathbb Z_{\geq0},
\qquad
\delta_z=\Phi_0-\varepsilon q_C(x)
          =p\|P_{-\varepsilon}x\|^2,
\quad P_{\varepsilon}=(I+\varepsilon C/p)/2.
\tag{3}
\]

The integrality of $\Delta_z/2$ follows because every Boolean quadratic
score of a signing has the same parity as $\binom n2$, which is odd here.
An **active** state has $\Delta_z=0$. Direct subtraction gives

\[
\sum_{e\in H}a_e(z)=2+\Delta_z/2-\delta_z/2.
\tag{4}
\]

The following are existing results, not new claims of this note.

- Propositions 15.629--15.630, used in Proposition 15.755, give
  $\delta_z=0$ on the signed Boolean eigenshell and
  $\delta_z\geq2p$ everywhere else. Equality $\delta_z=2p$ has
  the integral-anchor representation used in Section 7 below.
- Proposition 15.764 gives the parity bridge: odd $h$ has a residual-(ii)
  deletion at shell level two exactly when some signed $H$-shell row has
  score three. Even $h$ has the Type-I deletion at shell level one
  exactly when some such row has score two. Failure therefore gives signed
  shell floors five and four, respectively.
- Section 3, equation (6), of
  [the 2026-08-29 global-minimality note](NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md)
  already gives the all-subsets restoration-witness hierarchy in its
  closest-global-minimizer setting. Section 2 below applies the same
  algebra to the threshold-minimal hypothesis (2); it is not a new
  hierarchy.

## 2. Restoration witnesses and the exact rounding target

For $D\subseteq H$, restore its edges toward $C$:

\[
A_D=A+2(C\circ1_D)=C\oplus(H\setminus D).
\]

For every signed state,

\[
\varepsilon q_{A_D}(x)=M-\Delta_z+2\sum_{e\in D}a_e(z).
\tag{5}
\]

Consequently (2), together with score parity, is equivalent to the existing
all-subsets witness condition

\[
\forall\,\varnothing\ne D\subseteq H\quad\exists z:\qquad
\sum_{e\in D}a_e(z)\geq1+\Delta_z/2.
\tag{6}
\]

For a singleton, (6) forces $\Delta_z=0$ and $a_e(z)=1$. Thus the
positive $H$-edge sets of the active signed antipodal states cover $H$.
For a pair, its witness has $\Delta_z\in\{0,2\}$ and both edge signs
positive. These are also part of the existing hierarchy.

The exact missing integral-rounding statement is the existence of a
**nonempty** $D\subseteq H$ satisfying

\[
2\sum_{e\in D}a_e(z)\leq\Delta_z
\quad\text{for every signed Boolean state }z.
\tag{7}
\]

By (5), (7) is precisely $\Phi(A_D)\leq M$, contradicting (2).
Replacing its state-dependent right side by a single worst-case norm bound
loses essential information, as Section 4 quantifies.

## 3. A uniform fractional valley using every signed state

Suppose every signed Boolean $C$-eigenvector has $H$-score at least

\[
T_H^{\varepsilon}(y)=\sum_{e\in H}a_e(\varepsilon,y)\geq r,
\qquad r>2.
\tag{8}
\]

This is a bound over both phases, not only a selected frame. For

\[
A_\lambda=(1-\lambda)A+\lambda C,
\qquad
\lambda=\frac{r-2}{p+r-2},
\qquad
\gamma=2(p-2)\lambda>0,
\tag{9}
\]

one has the exact uniform estimate

\[
\boxed{\Phi(A_\lambda)\leq M-\gamma.}
\tag{10}
\]

Indeed, every signed state obeys the identity

\[
\varepsilon q_{A_\lambda}(x)
=M+4\lambda-(1-\lambda)\Delta_z-\lambda\delta_z.
\tag{11}
\]

On the eigenshell, $\delta_z=0$ and (8) gives
$\Delta_z\geq2r-4$. Off the eigenshell,
$\Delta_z\geq0$ and $\delta_z\geq2p$. The two resulting upper
bounds coincide at (9), since

\[
(1-\lambda)(2r-4)=2p\lambda.
\]

The useful specializations are

| Hypothesis | $r$ | $\lambda$ | $\gamma$ |
| --- | ---: | ---: | ---: |
| Odd $h$, using only (1) and parity | $3$ | $1/(p+1)$ | $2(p-2)/(p+1)$ |
| Even $h$, no signed level-two row | $4$ | $2/(p+2)$ | $4(p-2)/(p+2)$ |
| Odd $h$, no signed level-three row | $5$ | $3/(p+3)$ | $6(p-2)/(p+3)$ |

This is a fractional statement; $A_\lambda$ is generally not a signing.
In particular it does not supply (7).

## 4. Parseval blocks the norm-only sufficient rounding method

The triangle inequality would imply (7) if some nonempty restoration had

\[
\Phi(A_D-A_\lambda)\leq\gamma.
\tag{12}
\]

But the degree-two Boolean characters are orthonormal. For $d=|D|$,

\[
\begin{aligned}
\Phi(A_D-A_\lambda)^2
&\geq\mathbb E_x q_{A_D-A_\lambda}(x)^2\\
&=4\bigl((1-\lambda)^2d+\lambda^2(h-d)\bigr)\\
&=4\bigl(\lambda^2h+(1-2\lambda)d\bigr).
\end{aligned}
\tag{13}
\]

For $0<\lambda\leq1/2$, which includes all three rows of the table for
$p\geq5$, every nonempty $D$ therefore has

\[
\Phi(A_D-A_\lambda)\geq
2\sqrt{\lambda^2h+1-2\lambda}.
\tag{14}
\]

In particular (12) is impossible whenever

\[
h>(p-2)^2-\frac{1-2\lambda}{\lambda^2}.
\tag{15}
\]

The simpler sufficient barrier $h>(p-2)^2$ follows as a weaker
corollary. For the universal odd choice $r=3$, the right side of (15)
is $5-4p<0$: **this norm-only sufficient method fails for every nonempty
restoration at every size** with that choice of fractional valley.

This is not a counterexample to (7), nor to another fractional choice or a
rounding proof that exploits $\Delta_z$ separately for each state.
It only excludes establishing (7) through (10) and the norm tolerance (12).

## 5. Active parity, cover counts, and inter-state cuts

For an active state, (4) and $\sum_Ha_e\equiv h\pmod2$ give

\[
\delta_z/2\equiv h\pmod2.
\tag{16}
\]

Thus even $h$ excludes every active first-shell state
$\delta_z=2p$, including the bad-anchor branch. Odd $h$ excludes
active eigenshell states. For even $h$, excluding shell level two also
excludes active eigenshell states.

An active state has exactly

\[
|\{e\in H:a_e(z)=1\}|=\frac{h+2}{2}-\frac{\delta_z}{4}
\tag{17}
\]

positive $H$-edges. In a branch with no active eigenshell state,
$\delta_z\geq2p$; no two active states can cover $H$, because their
positive counts sum to at most $h+2-p<h$. More generally a covering
family $z_1,\ldots,z_k$ must satisfy

\[
\sum_{j=1}^k\delta_{z_j}\leq2(k-2)h+4k.
\tag{18}
\]

These constraints do not bound $h$ above.

There is also genuine cross-state cut structure. For active
$z=(\varepsilon,[x])$, $w=(\eta,[y])$, put
$\sigma_i=x_iy_i$. Then

\[
a_{ij}(z)a_{ij}(w)=\varepsilon\eta\sigma_i\sigma_j.
\tag{19}
\]

Their edge-sign disagreement is therefore a cut when the phases agree,
and the complement of a cut when the phases differ. With
$T=\{i:x_i\ne y_i\}$, the signed $A$-cut in the $x$-gauge equals
zero for the same phase and $M$ for opposite phases. As
$n$ is even and $M$ is odd, $|T|$ is even in the first case and
odd in the second. No extra hereditary restriction for an individual
one-edge deletion is claimed; Proposition 15.755 already identifies the
redundancy of those individual cut restrictions.

## 6. Odd common-anchor bridge

Fix a Boolean eigenanchor $y$ with $Cy=\varepsilon p y$. Write

\[
S=T_H^{\varepsilon}(y),\qquad
\ell_i=\varepsilon y_i(Ay)_i,\qquad
b_{ij}=\varepsilon A_{ij}y_iy_j.
\]

Then $\varepsilon q_A(y)=\Phi_0-2S$. If two distinct single-coordinate
flips $y^i,y^j$ are active in phase $\varepsilon$, their equations give
$\ell_i=\ell_j=2-S$, and the two-flip bound gives

\[
\varepsilon q_A(y^{ij})=\Phi_0+2S-8+4b_{ij}\leq M,
\qquad S+2b_{ij}\leq2.
\tag{20}
\]

For odd $h$, $S\geq3$ is odd, so necessarily $S=3$ and
$b_{ij}=-1$. Thus an odd no-level-three configuration has at most one
active single-flip neighbor per Boolean eigenanchor, in each fixed phase.

For even $h$, the premise already fails: each $\ell_i$ is odd, whereas
$2-S$ is even. A formal $S=4$, negative-clique equality calculation
must not be presented as a realizable active-single-flip configuration.

## 7. Uniform first-shell distance, including nonaffine bad anchors

Let $\delta_{\varepsilon}(x)=2p$. The minimum-shell theorem gives

\[
v=x-2a e_i\in\ker_{\mathbb Z}(C-\varepsilon pI),\quad a\in\{\pm1\},
\qquad g=ax_i\in\{\pm1\}.
\tag{21}
\]

For $g=1$, $v=x^i$ is Boolean. For $g=-1$, $v_i=3x_i$; nothing
here requires this integral eigenanchor to be affine. Set

\[
a_{uv}=\varepsilon C_{uv}x_ux_v,
\qquad c_j=\varepsilon x_j(Cx)_j.
\]

Multiplying (21) by $C$ gives the exact local-field formulas

\[
c_i=p-2pg,\qquad c_j=p+2g a_{ij}\quad(j\ne i).
\tag{22}
\]

For example, if $x$ is also active and
$d_H^x(j)=\sum_{\{j,k\}\in H}a_{jk}$, its odd signed $A$-local fields
are at least one. Hence

\[
d_H^x(i)\leq\frac{p-2pg-1}{2},\qquad
d_H^x(j)\leq\frac{p-1}{2}+g a_{ij}\quad(j\ne i).
\tag{23}
\]

The exceptional bound is $-(p+1)/2$ in the good case and
$(3p-1)/2$ in the bad case.

Now let $y$ be a distinct same-phase first-shell point, let
$T=\{j:x_j\ne y_j\}$, and put $d=|T|$. Equal $C$-scores imply

\[
W_C^x(T)=\sum_{j\in T}c_j-2\sum_{\{u,v\}\subseteq T}a_{uv}=0.
\tag{24}
\]

If $i\notin T$, (22) and the number of internal pairs give

\[
W_C^x(T)\geq d(p-d-1).
\tag{25}
\]

If $i\in T$, put $U=T\setminus\{i\}$. In the good case,

\[
W_C^x(T)=p(d-2)-2\sum_{\{u,v\}\subseteq U}a_{uv}
\geq(d-2)(p-d+1).
\tag{26}
\]

Here $d=1$ gives exactly $-p$, not zero; $d=2$ gives zero and means
that $x$ and $y$ are the two single-coordinate flips of the same Boolean
anchor $v$. For $3\leq d\leq p$, (26) is strictly positive. In the bad
case,

\[
\begin{aligned}
W_C^x(T)
&=p(d+2)-4\sum_{j\in U}a_{ij}
       -2\sum_{\{u,v\}\subseteq U}a_{uv}\\
&\geq(d+2)(p-d+1),
\end{aligned}
\tag{27}
\]

which is strictly positive for $1\leq d\leq p$, with exact value
$3p$ when $d=1$. Equations (25)--(27) prove: two distinct same-phase
first-shell points either have distance two and a common Boolean anchor,
or have distance at least $p-1$. The $d=1$ exception is explicitly
excluded in all cases.

For signed antipodal classes use
$d_{\pm}([x],[y])=\min(d_H(x,y),n-d_H(x,y))$, taking distinct classes.
The preceding argument applies to whichever representative realizes this
minimum. Section 6 eliminates its distance-two option for odd $h$ with
no level-three row. Therefore its same-phase active first-shell classes obey

\[
\boxed{d_{\pm}([x],[y])\geq p-1.}
\tag{28}
\]

There is a separate exact opposite-phase constraint. If
$v=x-2a e_i$ and $w=y-2b e_j$ are the respective integral anchors,
their eigenspaces are orthogonal, whence

\[
x\mathbin\cdot y=2a y_i+2b x_j-4ab\,1_{i=j}.
\tag{29}
\]

For $i\ne j$, the dot product belongs to $\{-4,0,4\}$. For $i=j$,
writing $g=ax_i$, $h'=by_i$, and $t=x_iy_i$ gives
$x\cdot y=2t(g+h'-2gh')\in\{0,\pm4,\pm8\}$. Thus opposite-phase
first-shell representatives have distance between $n/2-4$ and
$n/2+4$; their antipodal distance is at least $n/2-4$.

## 8. What remains open

These deductions constrain the actual Paley cut-code states, not only a
formal maximum of affine rows. They neither force a common eigenshell
anchor across deletion witnesses nor construct a nonempty restoration
satisfying (7). Higher-defect active states remain allowed. No code-packing
bound or completeness of the good-anchor branch is assumed. The unresolved
all-size task is still to prove the shell bridge or otherwise exclude the
threshold-minimal configurations using their joint state-dependent slacks.

## 9. Exact finite illustration: pairs do not replace all subsets

The separate [bounded probe and exact repair record](threshold_valley_probe.json)
contains a complete-signing example at order six. This is **not** a Paley
example of order $p^2+1$ and is not a counterexample to the restoration
target (7). In the script's lexicographic unordered-edge mask convention,
take $A=2393$, $C=7810$, and $H=A\mathbin{\mathrm{xor}}C$, with $|H|=10$.
Exact evaluation of all 32 antipodal representatives, retaining both
energy signs, gives

\[
\Phi(A)=5,\qquad \Phi(C)=9,\qquad
\Phi\bigl(\tfrac34 A+\tfrac14 C\bigr)=4.
\]

Every one-edge and every two-edge restoration has norm greater than five.
Nevertheless the three-edge matching
\[
D=\{\{0,4\},\{1,5\},\{2,3\}\}
\]
restores to mask $2641$, with norm five. The complete 1023-subset check
finds eleven norm-five restorations, exactly two of minimum size three.
This particular $D$ is neither a cut nor the complement of a cut, so the
repair is not a labelled switching/global-sign equivalence. Both $A$
and the repaired signing have exact square $5I$; $C$ is not asserted
to be a conference reference. The one-pair integer/Fraction verifier was
replayed on soulkiller, NUKA, and Jellyfin independently of the GPU score
table.

Thus even the actual complete-graph Boolean norm, strict fractional
descent, and failure of all one- and two-edge restorations do not imply
full threshold inclusion-minimality. This finite example supplies the
repair rather than contradicting its existence. The original probe's
128-reference null result for the full condition remains inconclusive;
its signed-$C$-maximizer filter $M-4$ is stronger than the universal
odd-floor implication $M-2$. No broader probe is an inference from this
example.
