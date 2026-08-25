# The complete scaled-norm 3p-6 shell is a quartic saddle

Date: 2026-08-25. Proposition 15.640. Retain

\[
n=p^2+1,\qquad d={n\over2},\qquad
P={1\over2}(I+C/p),\qquad L^*=P\mathbb Z^n,
\]

and let \(p\ge11\) be an odd prime. Proposition 15.639 proves that the
complete shell

\[
\mathcal X=\{x\in L^*:2p\lVert x\rVert^2=3p-6\}
\]

is the disjoint union of negative signed triples \(\mathcal T\) and
point--square-circle vectors \(\mathcal O\). This note computes its complete
degree-four harmonic operator.

Let \(W\) be symmetric and admissible:

\[
PWP=W,\qquad \operatorname{diag}W=0,
\qquad F=\lVert W\rVert_F^2.
\]

As usual,

\[
H_W(x)=(x^TWx)^2
-{4\lVert x\rVert^2\over d+4}x^TW^2x
+{2F\over(d+2)(d+4)}\lVert x\rVert^4.             \tag{1}
\]

The result is that \(\sum_{x\in\mathcal X}H_W(x)\) has three eigenspaces,
the same square-circle tensor channels as Proposition 15.634. Its
eigenvalues are

\[
\begin{aligned}
\lambda_0&=-{p^4+2p^3-69p^2+136p+26\over p^2+5},\\
\lambda_-&={p^4-14p^3+89p^2-196p+24\over p^2+5},\\
\lambda_+&={p^4-10p^3+69p^2-176p-76\over p^2+5}.
                                                               \tag{2}
\end{aligned}
\]

For every \(p\ge11\),

\[
\boxed{\lambda_0<0<\lambda_-,\lambda_+.}           \tag{3}
\]

Thus this finite shell is not a spherical 4-design. It is a genuine
quartic saddle.

## A tight frame of square circles through one point

For each unoriented square circle \(S\), choose either signed complement
\(w_S\); the tensor \(w_Sw_S^T\) is independent of that choice. Proposition
15.633 gives

\[
\sum_S w_Sw_S^T=p^2(p-1)P.                        \tag{4}
\]

The refinement needed here is, for every coordinate \(i\),

\[
\boxed{
\sum_{S\ni i}w_Sw_S^T
=p^2\left(P-2(Pe_i)(Pe_i)^T\right).}              \tag{5}
\]

To prove (5), use a signed Paley automorphism to move \(i\) to infinity.
The square circles through infinity are the affine lines in the
\(R=(p+1)/2\) square directions. For one fixed direction there are \(p\)
parallel lines. Their signed complements have Gram matrix

\[
p(pI-J):                                             \tag{6}
\]

the diagonal inner product is \(p(p-1)\), and two distinct parallel lines
give \(-p\). Two lines in distinct square directions meet in two points, so
their complement words are orthogonal. Therefore the nonzero frame
eigenvalue in each parallel class is \(p^2\), with multiplicity \(p-1\),
and distinct classes span orthogonal spaces.

All these words lie in

\[
U_\infty=\{u\in\operatorname{im}P:u_\infty=0\}.
\]

Their total span has dimension

\[
R(p-1)={p^2-1\over2}=d-1=\dim U_\infty.
\]

Hence their frame operator is \(p^2\) times the orthogonal projector onto
\(U_\infty\), namely \(P-2(Pe_\infty)(Pe_\infty)^T\). Signed Paley
transport proves (5) at every point.

## The negative-triangle quartic sum

Put

\[
B=C\circ W,
\qquad S_{ijk}=B_{ij}+B_{ik}+B_{jk},
\qquad \tau_{ijk}=C_{ij}C_{ik}C_{jk}.
\]

Since \(CW=pW\) and \(\operatorname{diag}W=0\),

\[
B\mathbf1=\operatorname{diag}(CW)=0.              \tag{7}
\]

For each triangle with \(\tau_{ijk}=-1\), there is one antipodal pair of
unit signings whose three signed conference edges are negative. For either
signing,

\[
x^TWx=-2S_{ijk}.
\]

Thus

\[
\sum_{x\in\mathcal T}(x^TWx)^2
=4\sum_{i<j<k}(1-\tau_{ijk})S_{ijk}^2.             \tag{8}
\]

The unweighted sum is

\[
\sum_{i<j<k}S_{ijk}^2={n-4\over2}F.               \tag{9}
\]

Indeed, every edge square occurs \(n-2\) times, while (7) makes the sum of
the wedge cross terms equal \(-F\). In the \(\tau\)-weighted sum, all edge
squares vanish because the off-diagonal entries of \(C^2\) vanish. The
wedge terms give

\[
\sum_{i<j<k}\tau_{ijk}S_{ijk}^2
=\operatorname{tr}(CW^2)=pF.                     \tag{10}
\]

Substituting (9)--(10) into (8) yields

\[
\boxed{
\sum_{x\in\mathcal T}(x^TWx)^2
=2(p-3)(p+1)F.}                                  \tag{11}
\]

So the entire negative-triangle family is quartically scalar.

## The point--circle quartic sum

Fix one orientation \(w=w_S\) for each unoriented square circle. For every
point \(i\notin S\), the uniquely oriented representative with value
\(-1\) at \(i\), together with its negative, can be written as

\[
\mathord\pm x_{S,i},\qquad
x_{S,i}=Pe_i-{w_iw\over p}.                       \tag{12}
\]

Let

\[
q_S=w^TWw,qquad a_i=(Ww)_i.
\]

The zero diagonal of \(W\) gives

\[
x_{S,i}^TWx_{S,i}={q_S\over p^2}-{2w_i a_i\over p}.
\]

There are \(p(p-1)\) points off \(S\), and
\(\sum_iw_ia_i=q_S\). Summing both signs in (12) therefore gives

\[
2\sum_{i\notin S}(x_{S,i}^TWx_{S,i})^2
=2{p-5\over p^3}q_S^2
+{8\over p^2}\sum_{i\notin S}(Ww_S)_i^2.         \tag{13}
\]

It remains to sum the last term. Set \(v_i=We_i\). From (4)--(5),
\(Pv_i=v_i\), and

\[
(Pe_i)^Tv_i=e_i^TW e_i=0.
\]

Hence

\[
\begin{aligned}
\sum_{S,i\notin S}(Ww_S)_i^2
&=\sum_i v_i^T\left(\sum_{S\not\ni i}w_Sw_S^T\right)v_i\\
&=p^2(p-2)\sum_i\lVert v_i\rVert^2\\
&=p^2(p-2)F.                                      \tag{14}
\end{aligned}
\]

Equations (13)--(14) prove

\[
\boxed{
\sum_{x\in\mathcal O}(x^TWx)^2
=8(p-2)F+{2(p-5)\over p^3}\sum_S(w_S^TWw_S)^2.} \tag{15}
\]

Combining (11) and (15), the complete quartic evaluation operator is

\[
\boxed{
\sum_{x\in\mathcal X}(x^TWx)^2
=2(p^2+2p-11)F
+{2(p-5)\over p^3}\sum_S(w_S^TWw_S)^2.}          \tag{16}
\]

## Harmonic correction and spectrum

Every vector in \(\mathcal X\) has squared norm

\[
q={3(p-2)\over2p},
\]

and Proposition 15.639 gives

\[
N=|\mathcal X|={p^2(p-1)(p+7)(p^2+1)\over6}.
\]

The shell is invariant under signed \(\operatorname{PSL}(2,p^2)\).
Irreducibility on \(\operatorname{im}P\), followed by the trace, gives

\[
\sum_{x\in\mathcal X}xx^T={Nq\over d}P.           \tag{17}
\]

Using (17) in the two radial terms of (1), their combined scalar subtraction
is

\[
{2Nq^2\over d(d+2)}F
={3(p-1)(p+7)(p-2)^2\over p^2+5}F.                \tag{18}
\]

Therefore the scalar offset in (16) becomes

\[
a=2(p^2+2p-11)
-{3(p-1)(p+7)(p-2)^2\over p^2+5}.                 \tag{19}
\]

Proposition 15.634 gives the spectrum of the square-circle evaluation
operator \(W\mapsto\sum_S(w_S^TWw_S)^2\):

\[
0,\qquad p^3(p-1),\qquad p^3(p+1),                \tag{20}
\]

with multiplicities

\[
{n(p-1)(p-3)\over8},\qquad
{n(p-1)\over4},\qquad
{n(p-3)\over4}.                                  \tag{21}
\]

Substituting (20) into (16), then (18), gives exactly (2).

To certify the signs uniformly, put \(x=p-11\ge0\). The three relevant
positive numerators are

\[
\begin{aligned}
-\lambda_0(p^2+5)
 &=x^4+46x^3+723x^2+4668x+10476,\\
\lambda_-(p^2+5)
 &=x^4+30x^3+353x^2+2004x+4644,\\
\lambda_+(p^2+5)
 &=x^4+34x^3+465x^2+3036x+7668.
\end{aligned}                                    \tag{22}
\]

Every coefficient in (22) is positive, proving (3).

At \(p=11\), the spectrum and multiplicities are

\[
\boxed{
\left(-{582\over7}\right)^{1220},\qquad
\left({258\over7}\right)^{305},\qquad
\left({426\over7}\right)^{244}.}                \tag{23}
\]

This is the precise version of the “sphere with a dent” picture: the shell
is radial and quadratically isotropic, but its fourth moment has a large
negative generic channel and two positive square-circle channels.

## Norm-parity-twisted shadow and boundary

The scaled norm \(3p-6\) is odd. Proposition 15.631's Poisson phase is
therefore \(-1\), while homogeneity gives \(H_W(x/2)=H_W(x)/16\). The
transformed dual-shadow eigenvalues are

\[
-{\lambda_0\over16}>0,qquad
-{\lambda_-\over16}<0,qquad
-{\lambda_+\over16}<0.                           \tag{24}
\]

Thus the new shell supplies cancellation in exactly the two circle-image
channels and reinforcement on the circle kernel. This is exact information
for the multi-scale R1 attack, not a tail bound. For \(p\ge17\), possible
intervening even shells below \(3p-6\) are still unclassified, and all later
shells remain uncontrolled. R1, global QVAR, and the limit remain open.

Independent finite audits construct every square circle and verify (5)
entrywise at \(p=3,5,7,11\). Direct summation of (11) and (15) for random
admissible \(W\), over every classified family member, agrees at
\(p=5,7,11\). These computations check the algebra but are not inputs to the
uniform proof.

Targeted searches for the three rational functions in (2), the exact
\(p=11\) spectrum, the through-point frame identity (5), and the combined
Paley/Miquelian harmonic-shell description found no matching theorem.
Individual OEIS searches for the larger raw numerator values
`3456108`, `2461860`, and `2743668`, as well as a multi-value search, found
no entries. General harmonic-ETF and Miquelian-plane sources are adjacent
but do not state (2) or (5). This is a search record, not an unqualified
priority claim.

Evidence:

- src/e1_gmin_m4_prop15640.py;
- evidence/e1_gmin_m4_prop15640.json;
- tests/test_prop15640.py.
