# Equianharmonic threshold: exact excess classes and the even-syndrome barrier

## Scope

Work over \(\mathbb Q(q)/(q^2+q+1)\), with \(X=2x\), at the first
compact count allowed by the component-packing theorem,

\[
 b=B=2L-1={2r+7\over3},\qquad \Delta=L-1.
\]

This note classifies the component excess at that threshold and derives
exact degree-six/eight formulas for two four-compact trade families.  It
does **not** produce a finite-field zero-syndrome matching, a common global
form, a Boolean lift, or a proof of residual (ii).

## 1. Excess-one classification

For a connected pairing component with \(K\) compact vertices, \(E\)
all-equal vertices, cycle rank \(\mu\), and \(Z\) caps, the previously
proved identity is

\[
 w:=K-2\delta=4-K-4E-4\mu-2Z.                 \tag{1}
\]

Globally \(\sum w=b-2\Delta=1\).  Positive component mass is at most
three by the disjoint-support packing theorem, so the positive/negative
mass pair is exactly one of

\[
 (P,N)=(1,0),(2,1),(3,2).                     \tag{2}
\]

The positive blocks are respectively

\[
\begin{array}{c|c}
P&\text{possible disjoint blocks}\\ \hline
1&\mathcal P\text{ or }\mathcal F_j\\
2&\mathcal H\text{ or }\mathcal F_j+\mathcal P\\
3&\mathcal H+\mathcal P.
\end{array}                                    \tag{3}
\]

The nonpositive tuples needed by (2), written \((K,E,\mu,Z)\), are

\[
\begin{array}{c|l}
w=0 &(0,1,0,0),(2,0,0,1),(4,0,0,0)\\
w=-1&(1,0,0,2),(1,0,1,0),(1,1,0,0),(3,0,0,1),(5,0,0,0)\\
w=-2&(0,1,0,1),(2,0,0,2),(2,0,1,0),(2,1,0,0),
       (4,0,0,1),(6,0,0,0).
\end{array}                                    \tag{4}
\]

This is a finite integer consequence of (1), not a prime census.

## 2. Two genuine zero-excess trade families

Put

\[
 \Phi(z)=qz+1-q,\quad
 A_2=x,\ A_0=\Phi(x),\ A_1=\Phi^2(x),\quad
 B_2=-x,\ B_0=\Phi(-x),\ B_1=\Phi^2(-x).
\]

Write \(K(a,b;c)\) for a compact atom whose positive edge is \(ab\).
The following blocks have four compact atoms and replace the two
all-equal triangles \(A\) and \(B\) exactly in every odd edge channel:

\[
\begin{split}
U(x)=\{&K(B_2,-A_1;-A_0),K(A_1,A_2;-B_1),\\
       &K(A_1,-B_1;-B_0),K(B_0,B_2;-A_1)\},       \tag{5}\\
V(x)=\{&K(A_0,A_2;-B_0),K(-A_1,B_0;-A_0),\\
       &K(-B_0,A_1;-B_1),K(B_2,B_1;-A_1)\}.       \tag{6}
\end{split}
\]

For \(d\in\{6,8\}\), let \(M_{d,j}(a,b)=(a-b)^2(ab)^j
(a+b)^{d-2-2j}\), and let \(D^T_{d,j}(x)\) be the sum over the
compact block \(T\in\{U,V\}\) minus the two all-equal triangles.  Each
\(D^T_{d,j}\) is an explicit polynomial in \(x\) over
\(\mathbb Q(q)/(q^2+q+1)\).  The complete low-to-high coefficient vectors
are emitted by
`trade_deviation_polynomials("U")` and
`trade_deviation_polynomials("V")` in the source verifier; the three
degree-six channels precede the four degree-eight channels.

For \(U\) alone there is the affine identity

\[
\begin{split}
&{84461\over12096}D^U_{6,0}-{187801\over6048}D^U_{6,1}
+{619\over63}D^U_{6,2}-{2509\over5184}D^U_{8,0}\\
&\quad+{95099\over36288}D^U_{8,1}-{3557\over1134}D^U_{8,2}
+D^U_{8,3}={10176\over7}.                       \tag{7}
\end{split}
\]

## 3. Exact mixed-family barrier

Take four \(U\) parameters at \(-1,-2,1,0\) and three \(V\) parameters
at \(-3,-2,-1\).  The Jacobian of the seven syndrome coordinates has

\[
 \det J=
4128623683475967290061619200
=2^{32}3^{26}5^2\,7\,2161.                    \tag{8}
\]

Consequently the mixed seven-channel polynomial map is dominant in
characteristic zero, and the displayed integral specialization remains
nondegenerate in every characteristic outside
\(\{2,3,5,7,2161\}\).  Thus the component-excess invariant by itself
cannot yield a universal affine degree-six/eight obstruction once both
trade families occur.  What remains is genuinely finite-field and global:
one must use rational parameter availability, disjoint \(\Phi\)-cycle
matching, or an additional invariant that excludes a mixed family.

As a formula check, the exact atom witnesses give

\[
\begin{array}{c|c|c}
(p,T,x)&(D_{6,0},D_{6,1},D_{6,2})&(D_{8,0},D_{8,1},D_{8,2},D_{8,3})\\ \hline
(31,U,18)&(3,9,8)&(20,14,21,4)\\
(43,U,38)&(35,23,14)&(36,12,8,25)\\
(43,V,7)&(36,21,18)&(8,14,20,40).
\end{array}                                     \tag{9}
\]

These replays validate (5)--(8); neither witness has simultaneous zero
degree-six/eight syndrome.

## 4. Replay

```bash
PYTHONPATH=src python -m e1_gmin_m4_equianharmonic_threshold_even_barrier
pytest -q tests/test_equianharmonic_threshold_even_barrier.py
```

The returned record deliberately keeps all finite-field zero matching,
global-lift, Boolean-lift, and residual-(ii) flags false.
