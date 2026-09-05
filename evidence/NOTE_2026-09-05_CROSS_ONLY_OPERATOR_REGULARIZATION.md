# Cross-only operator regularization with fixed internal blocks

2026-09-05. This is a same-order theorem for the conditional objective:
A and -A are kept EXACTLY fixed. The loss includes an explicit internal
boundary penalty. It is not the previously proved whole-source trimming
theorem, which may change the internal blocks.

The bounded-source corollary gives a vanishing cap error when
||A||/sqrt(n) is bounded. It does not prove the analogous assertion,
uniformly in n, for every arbitrary exact original minimizer A.

## 1. Objective and statement

Let n>=2, let A be a complete symmetric zero-diagonal order-n signing,
and let B be an n by n sign matrix. For rectangular real matrices put
\[
 \beta(T)=\max_{x,y\in\{-1,1\}^n}|x^TTy|,
\]
with the dimensions of x,y adjusted for rectangular submatrices. Set
\[
 F_A(B)=\max_{x,y}|Q_A(x)-Q_A(y)+x^TBy|
       =\max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TBy|\bigr).     \tag{1}
\]
The equality follows by replacing x by -x. In particular
\(\beta(B)\le F_A(B)\). Write
\[
 \kappa=\log(1+\sqrt2),\qquad \Lambda=\pi/\kappa.
\]

**Theorem.** For every K>0 there are row and column sets S_1,S_2,
with \(T_i=[n]\setminus S_i\) and \(s=|S_1|+|S_2|\), and a complete
cross signing B' such that
\[
 s\le{\Lambda\beta(B)\over K\sqrt n},\qquad
 B'_{T_1,T_2}=B_{T_1,T_2},\qquad
 \|B'\|_{\rm op}\le(K+8)\sqrt n,                         \tag{2}
\]
and
\[
 \boxed{F_A(B')\le F_A(B)
      +\beta(A_{S_1,T_1})+\beta(A_{S_2,T_2})+2n\sqrt s.} \tag{3}
\]
Empty rectangular cuts have beta zero. If s=0, take B'=B.

Thus, if \(F_A(B)\le Cn^{3/2}\) and
\(\|A\|_{\rm op}\le K_A\sqrt n\), then
\[
 \boxed{F_A(B')\le F_A(B)
      +(2+\sqrt2 K_A)\sqrt{\Lambda C/K}\ n^{3/2}.}        \tag{4}
\]
Both internal blocks remain A,-A throughout (2)--(4).

## 2. The precise existing SDP prerequisite and its cross specialization

The only imported result is the finite tensor-rounding inequality (6)
in Section 2 of
NOTE_2026-09-05_NORM_CAP_FIELD_RESPONSE.md:
\[
 \left|\sum_{ij}T_{ij}\langle u_i,v_j\rangle\right|
       \le{\pi\over2\kappa}\,\beta(T)                    \tag{5}
\]
for every real array T and unit vectors u_i,v_j. Its proof in that
note constructs tensor lifts with inner product
\(\sin(\kappa\langle u_i,v_j\rangle)\), applies finite-dimensional
Gaussian sign rounding, and obtains the exact factor \(2\kappa/\pi\).
In particular (5) is not restricted to complete symmetric signings.

For clarity, here is the full diagonal-majorizer specialization needed
for the present restricted scope. Put
\[
                  J_B=\begin{pmatrix}0&B\\B^T&0\end{pmatrix}.
\]
For the semidefinite program
\[
 \max\{\operatorname{tr}(J_BX):X\succeq0,\ \operatorname{diag}X=1\},
\]
a Gram representation and (5) give the bound
\(2(\pi/(2\kappa))\beta(B)=\Lambda\beta(B)\).
The diagonal dual therefore supplies a nonnegative diagonal D with
\[
            D-J_B\succeq0,\qquad
            \operatorname{tr}D\le\Lambda\beta(B).        \tag{6}
\]
The primal X=I and a sufficiently large scalar diagonal dual are
strictly feasible. Strong duality applies, and the dual optimum is
attained because its nonnegative diagonal sublevels are compact.

Conjugating (6) by \(\operatorname{diag}(I_n,-I_n)\) leaves D fixed
and changes J_B to -J_B. Thus the SAME D also satisfies
\(D+J_B\succeq0\); there is no need to add a second majorizer.

Select S_1 from the first block and S_2 from the second block by
the threshold \(D_{ii}>K\sqrt n\). The trace bound proves the first
inequality in (2). Principal restriction to T_1 together with T_2
gives
\[
       \|B_{T_1,T_2}\|_{\rm op}\le K\sqrt n.              \tag{7}
\]
This is the cross-block form of the already established SDP trim,
not an invocation of a complete-signing theorem on J_B's zero blocks.

## 3. One cross-only filler satisfies both estimates

Let B_0 agree with B on \(T_1\times T_2\) and be zero elsewhere.
On the remaining
\[
       e=n^2-|T_1||T_2|=ns-|S_1||S_2|\le ns
\]
cross entries, choose independent fair signs, with zero entries on
the retained rectangle. Call the resulting real matrix G. The
candidate is \(B'=B_0+G\).

For real unit u,v, the squared sum of the independent coefficients
of \(u^TGv\) is at most
\(\sum_{ij}u_i^2v_j^2=1\). Hence
\[
                  \Pr\{|u^TGv|>t\}\le2e^{-t^2/2}.
\]
Take 1/4-nets of both real unit spheres, each of size at most 9^n.
The bilinear net maximum controls the operator norm by a factor two.
At \(t=4\sqrt n\),
\[
 \Pr\{\|G\|_{\rm op}>8\sqrt n\}
                \le2e^{-(8-2\log9)n}.                  \tag{8}
\]

For each Boolean pair x,y, \(x^TGy\) is a sum of e independent
fair signs. A union bound over the \(2^{2n}\) pairs at
\(t=2\sqrt{ne}\) gives
\[
 \Pr\{\beta(G)>2\sqrt{ne}\}
                \le2e^{-(2-2\log2)n}.                  \tag{9}
\]
The bounds (8) and (9) sum to less than one for n>=2.
At n=2 their sum is \(2\cdot9^4e^{-16}+32e^{-4}<1\),
and both decrease thereafter. Thus one filler simultaneously has
\[
           \|G\|_{\rm op}\le8\sqrt n,\qquad
           \beta(G)\le2\sqrt{ne}\le2n\sqrt s.            \tag{10}
\]
The two success events need not be independent. This is an existence
proof; no sampling experiment or unchanged numerical rerun is needed.
Equations (7) and (10) prove the operator bound in (2).

## 4. Why deleting the cross strips has an INTERNAL boundary cost

For S subset [n], write A^(S) for A with just the edges between
S and its complement set to zero. For fixed x,y, multiply all
coordinates of x in S_1 by a single independent fair sign epsilon,
and all coordinates of y in S_2 by another independent fair sign eta.
Averaging the ORIGINAL energy gives exactly
\[
 Q_{A^{(S_1)}}(x)-Q_{A^{(S_2)}}(y)+x^TB_0y.
\]
The same-set internal edges survive these flips. The internal cuts
and every removed cross entry average to zero. Convexity of absolute
value therefore gives
\[
 \max_{x,y}|Q_{A^{(S_1)}}(x)-Q_{A^{(S_2)}}(y)+x^TB_0y|
                                                   \le F_A(B).
\]
Returning to the original fixed internal blocks costs at most
\[
                   \beta(A_{S_1,T_1})+\beta(A_{S_2,T_2}).
\]
Finally adding G costs at most beta(G), proving (3).

This is the step that is NOT supplied by whole-source regularization.
It would be incorrect to claim \(F_A(B_0)\le F_A(B)\) by averaging
while silently retaining the original internal cuts.

For \(L_A=\|A\|_{\rm op}\),
\[
 \beta(A_{S_i,T_i})\le L_A\sqrt{|S_i|(n-|S_i|)}.
\]
Consequently their sum is at most
\(L_A\sqrt{2ns}\le\sqrt2 K_A n\sqrt s\).
Substituting the size bound in (2) proves (4).

## 5. Conditional minima, cap selection, and Gaussian-floor slack

Let \(F_A^*=\min_{B\in\{-1,1\}^{n\times n}}F_A(B)\).
The elementary independent-cross-sign estimate proves
\[
             F_A^*\le2\Phi(A)+2\sqrt{\log2}\ n^{3/2}.     \tag{11}
\]
Indeed the rank-one sign coefficient class has at most \(2^{2n-1}\)
members, already includes both signs, and the expected maximum of
its n-subgaussian linear forms is at most
\(n\sqrt{2(2n-1)\log2}\).

Suppose \(\Phi(A)\le C_A n^{3/2}\) and
\(\|A\|_{\rm op}\le K_A\sqrt n\), and set
\(C_B=2C_A+2\sqrt{\log2}\).
For epsilon_B>0, applying (4) to an actual conditional minimizer gives, for
\[
                 K\ge
       {\Lambda C_B(2+\sqrt2 K_A)^2\over\varepsilon_B^2},
\]
a complete B' with
\[
 \|B'\|_{\rm op}\le(K+8)\sqrt n,\qquad
 F_A^*\le F_A(B')\le F_A^*+\varepsilon_B n^{3/2}.          \tag{12}
\]
The cap is selected independently of n when C_A,K_A,epsilon_B are fixed.

This B' need not itself be an exact conditional optimizer. Therefore
an intrinsic sign-to-Gaussian floor based at B' has the EXTRA slack
\(\varepsilon_B n^{3/2}\) if written with \(F_A(B')\) on its left.
The valid reasoning is: every rounded competitor has norm at least
\(F_A^*\), and (12) compares \(F_A(B')\) with that minimum. No exact
optimality property is silently transferred to the regularized B'.

## 6. What prior source regularization actually permits

Section 1 of
NOTE_2026-09-05_SAME_ORDER_SPECTRAL_REGULARIZATION.md applies to an
exact original order-n minimizer, with a uniform norm cap C_0.
For \(0<\varepsilon_A\le1\), it supplies an ACTUAL order-n signing
A' such that
\[
 \Phi(A')\le m_n+\varepsilon_A n^{3/2},\qquad
 \|A'\|_{\rm op}\le K_A\sqrt n,\qquad
 K_A=8+{4\Gamma C_0\over\varepsilon_A^2},\quad
 \Gamma=4\pi/\kappa.                                   \tag{13}
\]
One may take the elementary universal C_0=2 used in that theorem.
Thus (12) can be applied with this A' kept fixed,
\(C_A=C_0+\varepsilon_A\), and the printed K_A. This selects a
near-original source and then a near-conditional cross block, with
all losses and cap dependencies retained:
\[
      K_A=O(\varepsilon_A^{-2}),\qquad
      K=O(\varepsilon_A^{-4}\varepsilon_B^{-2}).
\]
In particular equal epsilon choices give a cross cap O(epsilon^(-6)).

No relation \(F_{A'}^*\approx F_A^*\) for a prescribed exact source
A is asserted. A' is only near-optimal for the ORIGINAL norm and
need not preserve any theorem hypothesis requiring exact source
optimality. Any sharp comparison subsequently used must apply to
this near-original bounded-operator class, with its stated slacks.

## 7. The unresolved arbitrary-exact-source scope and cap competition

For an arbitrary exact original minimizer the presently proved
\(L_A^2\le8\Phi(A)=O(n^{3/2})\) gives only \(K_A=O(n^{1/4})\).
With that bound inserted into (4), forcing the loss to be little-o
requires K to grow faster than n^(1/2), at which point
\((K+8)\sqrt n\) is weaker than the trivial cross-operator bound n.
Therefore (4) alone does NOT solve uniform cross-only regularization
for every exact original minimizer.

The exact boundary-sensitive formula (3) remains valid without any
operator assumption on A. A genuinely stronger selected-boundary
estimate could improve this conclusion; none is proved here.

Finally, the regularization slack decays like
\((K_A+1)K^{-1/2}\). It cannot automatically be absorbed into a
proposed sharp Gaussian improvement that decays like K^(-1),
K^(-2), or faster. The cap dependencies must be checked in the
evaluated upper; existence of very large fixed caps is not itself
a sharp comparison or a proof of original convergence.

