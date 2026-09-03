# Residual-(ii) closure DAG after grouped uncertainty

**Date:** 2026-09-03

**Scope:** logical propagation only.  This note does not change a canonical
status document or a live acceptance predicate.

## 1. What the grouped theorem proves

Let

\[
  \Delta=(\mathbf F_p^2\setminus\{0\})/\{\pm1\},\qquad
  h={p-1\over2},\qquad d=p+1,
\]

and let `D=(C,Phi)` be the halved symmetric binary map of (E.23).  The new
square-remainder argument proves, for every odd prime and every nonzero point
word `f`,

\[
 \operatorname {wt}(f)+
 \#\{A:(M^{\mathsf T}f)|_{B_A}\ne0\}\ge p+1.       \tag{1}
\]

The already-proved active-direction symmetric-difference argument (E.29a)
therefore gives

\[
                  d_{\rm row}(D)\ge ph.             \tag{2}
\]

The fixed-transverse word `X_(L,beta)` of (E.24) is a nonzero row word of
weight `ph`.  Hence the exact conclusion is

\[
                  \boxed{d_{\rm row}(D)=ph}.         \tag{3}
\]

The subsequent row-code gap theorem sharpens this: every word of weight
`ph` is a fixed-transverse rectangle and no weight lies strictly between
`ph` and `|Delta|`.

If `D_U` retains precisely the columns outside a puncture `U`, the dual
criterion now gives the unconditional implication

\[
                     |U|<ph\quad\Longrightarrow\quad D_U\text{ is onto}.
                                                               \tag{4}
\]

Thus (1)--(3) settle the minimum-distance question and (4) settles every
punctured **linear parity** problem below that distance.  They do not put a
solution in a prescribed Hamming slice.

## 2. Exact balanced branch-C threshold

For `p=4r+3`,

\[
 h=2r+1,\qquad |\Delta|=(p+1)h,
 \qquad ph=|\Delta|-h.
\]

In the all-active balanced branch-C ledger, put

\[
 t_{\max}=4r^2-2r-5,
 \qquad N={p+1\over2}(p-1)=|\Delta|.
\]

Equation (E.26) says that any support capable of passing the necessary
Hamming equation has

\[
 \kappa\ge t_{\max}-t+1,
 \qquad
 |U|=N-2\kappa
      \le N-2(t_{\max}-t+1).                       \tag{5}
\]

Since `h=2r+1` is odd,

\[
\begin{aligned}
 N-2(t_{\max}-t+1)<ph=N-h
 &\iff 2(t_{\max}-t+1)>h\\
 &\iff t_{\max}-t+1\ge r+1\\
 &\iff \boxed{t\le t_{\max}-r=4r^2-3r-5}.         \tag{6}
\end{aligned}
\]

Consequently the distance bound alone makes the all-active puncture linearly
harmless throughout (6).  At `p=31` (`r=7`), this is `t<=170`.
The later gap theorem removes this threshold: any physically extendable
all-active support has `|U|<=|H|<|Delta|`, and every possible word below
`|Delta|` is a minimum rectangle which the Mobius union cannot contain.
If only `q<=h` centres are nonzero, `|U|<=q(p-1)<ph` directly.
Thus the actual structured punctured map is onto throughout the balanced
zero-odd branch-C regime.

The hypotheses in that sentence matter.  Formula (5) is the all-active
support ledger.  It must not be silently promoted to zero hard centres,
unbalanced allocations, nonzero global forms, branch B, or non-isolated
later shells.

## 3. The remaining symmetric Boolean fibre

Surjectivity in (4) says only that the parity equation has some solution.
The live problem still requires a solution in one exact Boolean slice:

\[
 \sum_{O\notin U} b_O\widehat B_O=\widehat T_U\pmod2,
 \qquad b_O\in\{0,1\},                              \tag{7}
\]

\[
 2\sum_{O\notin U}b_O
   =|H|-|U|-|a(T_U)|,                                \tag{8}
\]

together with the directionwise integral counts

\[
 n_L={P_L-u_L-f_L\over2},                            \tag{9}
\]

their admissible bounds, and all exact signed directional targets.  Equation
(4) does not imply (8) or (9), does not prove nonnegative `0/1` realizability,
and does not construct or exclude one common simple graph.

Therefore the strongest honest status after (1) is:

* grouped uncertainty for all supports and all odd primes: **proved**;
* `d_row(D)=ph`: **proved**;
* `D_U` onto whenever `|U|<ph`: **proved**;
* minimum-word classification and empty interval `ph<wt<|Delta|`: **proved**;
* actual structured branch-C punctured map over `F_2`: **proved onto**;
* prescribed-weight and direction-sliced Boolean fibre (7)--(9): **open**;
* residual (ii), E(1), and `L=1/2`: **open**.

## 4. Global residual predicate and stale wiring

The grouped theorem belongs only to the balanced symmetric-box branch.  It
does not prove the global no-descent statement for all primes, all admissible
sizes, unbalanced allocations, nonzero global forms, branch B, zero-centre
cases, or non-isolated shells.

Moreover, `residual_ii_k_ge_4p_ND_closed()` currently conjoins historical
route predicates that are known false independently of the new result:
the dual-bad mass predicate at `k=4p`, the retracted one-sided level-four
dichotomy, the Lipschitz wrapper that imports that dichotomy, and the open
rigidity/slice close.  Merely flipping the umbrella multi-level flag would
not propagate.  A future genuine global theorem must replace or bypass that
obsolete conjunction; the retracted predicates must not be relabelled as
proved.

The correct eventual DAG is

\[
 \text{global residual-(ii) no-descent theorem}
 \Longrightarrow
 \texttt{residual\_ii\_k\_ge\_4p\_ND\_closed}
 \Longrightarrow
 \begin{cases}
  \texttt{residual\_ii\_full\_closed},\\
  \texttt{e1\_open\_residuals}=\varnothing,\\
  \texttt{four\_e1\_units\_closed}=\mathrm{True},
 \end{cases}
 \Longrightarrow E(1)\Longrightarrow L=\tfrac12.             \tag{10}
\]

The grouped theorem proves a strict local ancestor of one branch of the
first node in (10), not that node itself.
