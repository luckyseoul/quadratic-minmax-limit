# Source fixed words of localized Mobius halves

Date: 2026-09-03

Status: proved the auxiliary-independent parity of the source-side
\(\Phi\) word of one localized Mobius half and the exact \(\Phi\) word of
the rigid sharp two-cancellation pair. A previous draft incorrectly
identified the full branch-C target fixed word with the hard literal alone.
That identification is retracted: a centrally symmetric compact residual
can have odd fixed-cell coefficients. Consequently no Hamming-parity
obstruction, symmetric completion, or closure of residual (ii) is claimed.

## 1. Auxiliary-independent parity of one half

Fix one hard functional \(L\), a nonzero center \(j\), and any auxiliary
\(M\) independent of \(L\). In the standard \(L,M\) coordinates, the
selected half of the localized trade is

\[
 u_t=\left(j,{jt\over t+1}\right),\qquad
 v_t=(jt,jt),\qquad t\in\mathbf F_p\setminus\{-1\}.       \tag{1}
\]

The midpoint and half-difference of an edge are parallel exactly when its
two endpoints are collinear with the origin. Here

\[
                         \det(u_t,v_t)={j^2t\over t+1}.   \tag{2}
\]

Thus exactly \(t=0\) has zero fixed word. Every other parameter contributes
one nonzero paired-affine-line block \(\Phi(O)\), whose weight is \(p\).
Hamming-weight parity is a linear functional over \(\mathbf F_2\), so the
XOR of all the half's fixed words has parity

\[
                    (p-2)p\equiv1\pmod2.                 \tag{3}
\]

This statement uses neither the auxiliary nor the orientation of the half.
Paley source signs also disappear modulo two.

Equation (3) is only a source-side fact. It does not determine the parity of
the forced fixed vector

\[
                       a(T_U)=a_Y+\sum_{O\in U}\Phi(O).
\]

The missing target contribution is isolated in Section 3.

## 2. Exact rigid-pair Phi word

Take distinct hard functionals \(L_1,L_2\) with nonzero centers \(j_1,j_2\)
and normalize

\[
                       x={L_1(v)\over j_1},\qquad
                       y={L_2(v)\over j_2}.
\]

The two-cancellation locus theorem forces

\[
 q=r={1\over2},\qquad A=B={3\over2}.
\]

For the first selected half, parameter \(t\ne-1\) has endpoints

\[
 u_t=\left(1,{1\over2}-{3\over2(t+1)}\right),\qquad
 v_t=\left(t,{t\over2}\right).
\]

At \(t=0\) the edge meets the origin and its fixed word is zero. For
\(t\ne0,-1\), formula (13) of the fixed-edge elimination note is the paired
affine line through \(v_t\) and \(-u_t\). A normalized equation for it is

\[
 \ell_t(x,y)=
 { \{3-(t+1)^2\}x+2(t+1)^2y\over3t}=\pm1.                \tag{4}
\]

Put \(D=2y-x\). For a chosen right-hand sign
\(\varepsilon\in\{\pm1\}\), equation (4) is

\[
 D t^2+(2D-3\varepsilon)t+2(x+y)=0.                       \tag{5}
\]

The parity of the root count of a nondegenerate quadratic is one exactly
when its discriminant vanishes. Here the discriminant is

\[
                         9-12D(x+\varepsilon).
\]

The two spurious \(t=0\) roots occur together and cancel modulo two. The
excluded value \(t=-1\) occurs once exactly when \(x^2=1\). The linear case
\(D=0\) contributes twice and also cancels. Thus the exact one-half word is

\[
 \boxed{
 F(x,y)=\mathbf1_{x^2=1}
 +\sum_{\varepsilon=\pm1}
  \mathbf1_{(2y-x)(x+\varepsilon)=3/4}\pmod2.}            \tag{6}
\]

The second half has word \(F(y,x)\). Therefore

\[
                         \Phi_{\rm pair}=F(x,y)+F(y,x).   \tag{7}
\]

The two conics in (6) have \(p-1\) points each and are disjoint. Combining
them with the two lines \(x=\pm1\) gives \(4p-6\) nonzero points, hence one
half has fixed-word weight

\[
                              |F|=2p-3.                   \tag{8}
\]

For branch primes \(p\ge31\), direct intersection counting gives

\[
 |F(x,y)\cap F(y,x)|_{\rm points}
 =20+2\eta(3)+4\eta(6)+4\eta(-2).                         \tag{9}
\]

Indeed, away from \(x,y=\pm1\), the four conic-pair systems contribute
\(8+2\eta(3)\). The four corners contribute four. The cases with exactly
one boundary coordinate contribute
\(4\{2+\eta(6)+\eta(-2)\}\). These cases are disjoint for \(p\ge31\).
Consequently

\[
 \boxed{
 |\Phi_{\rm pair}|
 =4p-26-2\eta(3)-4\eta(6)-4\eta(-2).}                    \tag{10}
\]

At \(p=31\), this weight is 108. Formula (10) is invariant under the choice
of the two distinct directions and their nonzero centers, because normalized
coordinates act by a permutation of the antipodal fixed-edge classes.

## 3. Exact target-side blocker

Modulo two, the hard literal does contribute the unit-star chain \(S_j\).
For nonzero \(j\), its fixed cells alone would map under the explicit
fixed-edge inverse to

\[
                         \{[v]:L(v)^2=j^2\},
\]

an odd \(p\)-element word. It is not, however, the full hard target.

The actual hard coefficient chain in the antisymmetric note is

\[
                      C_L=-S_{j_L}+K_L^{\rm compact}.
\]

The theorem that \(K_L^{\rm compact}\) is centrally symmetric only kills
its antisymmetric part. It does not make its fixed-cell coefficients even.
For example, the compact atom \(K(v,-v;0)\) is centrally symmetric and has
an odd coefficient on the fixed antipodal cell \(\{v,-v\}\) modulo two.
Equivalently, the factor \(2B_L\) in the value formula does not imply that
the coefficient chain of \(4B_L\) is even.

Thus the correct target decomposition has an additional, presently
uncontrolled word

\[
       a_Y=a_{\rm literal}+a_{\rm compact},\qquad
       a(T_U)=a_{\rm literal}+a_{\rm compact}+\Phi(U).    \tag{11}
\]

The source parity (3) may cancel the literal parity, but
\(a_{\rm compact}\) can change both parity and weight. No theorem currently
sets it to zero or bounds its coset weight. This is the exact flaw in the
discarded Hamming-parity argument.

Accordingly, this note does not exclude the one-trade-per-hard-direction
Mobius ansatz. The next target-side question would have to control the
fixed-cell parity word of the actual compact residual, not merely its
centrality.

## 4. Phi weights depend on the rigid pairing

Even the source-side \(\Phi\) weight is not determined solely by \(p\), the
number of pairs, or the unpaired hard data.

At \(p=31\), take centers one and the same four directions

\[
 X=(1,0),\quad Y=(0,1),\quad Z=(1,1),\quad W=(1,2).
\]

All three rigid pairings give ternary four-trade sums, exactly four
cancellation units, and 112 used inversion orbits. Their source
\(\Phi\)-weights are nevertheless

\[
\begin{array}{c|c}
\text{pairing}&|\Phi(U)|\\ \hline
XY\mid ZW&174\\
XZ\mid YW&176\\
XW\mid YZ&172.
\end{array}
\]

This is an exact same-data counterexample to a pairing-independent
\(\Phi\)-weight. It says nothing about the full target coset (11).

## Reproduction

The implementation records the symbolic formulas. The one \(p=31\)
four-direction replay is an exact fail-when-wrong counterexample, not a
prime or configuration census.

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_rigid_pair_fixed_word.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_rigid_pair_fixed_word.py
