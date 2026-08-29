# Global minimality, local-stability no-go, and the surviving global routes

**Date:** 2026-08-29
**Status:** proved algebra and a general counter-mechanism; kills the local/product
route but does not prove \(k_\star=o(n^2)\) or settle \(L\)

Let \(C\) be a fixed signing on \(K_n\), put \(N=\binom n2\), and for a
signing \(A\) define \(W=A\circ C\). Distance in this note includes both
Seidel switching and global sign:

\[
d_\pm(A,C)=\min_{\varepsilon\in\{\pm1\},\ s\in\{\pm1\}^n}
d_H(A,\varepsilon D_sCD_s).
\]

The global sign is necessary: \(-A\) has the same value of \(\Phi\) as
\(A\), while switching alone can put those representatives quadratically
far apart.

## 1. Exact product transform

For every \(A,C\),

\[
\boxed{d_\pm(A,C)=\frac{N-\Phi(A\circ C)}2.}
\tag{1}
\]

Indeed, the correlation of \(A\) with the switching \(D_sCD_s\) is
\(Q_W(s)\). Minimizing Hamming distance while allowing global sign is
therefore the same as maximizing \(|Q_W(s)|\).

Now take \(C\) to be the Paley conference signing of order
\(n=p^2+1\). The general Paley frame identity proved in Proposition 15.189 is

\[
\mathbb E_{y\in\mathrm{Max}_+}y_i y_j=C_{ij}/p
\qquad(i\ne j).
\]

Coordinatewise multiplication gives, for every \(x\in\{\pm1\}^n\),

\[
\boxed{
\mathbb E_{y\in\mathrm{Max}_+}Q_A(xy)=\frac{Q_W(x)}p.}
\tag{2}
\]

Consequently

\[
\Phi(A)\ge\frac{\Phi(W)}p
=\Phi(C)-\frac{2d_\pm(A,C)}p,
\]

which is the pointwise form of Max-Lipschitz. There is also the exact Walsh
convolution

\[
\boxed{Q_W(x)=\mathbb E_z[Q_A(z)Q_C(xz)].}
\tag{3}
\]

Only equal degree-two Walsh characters survive the expectation.

The complete second-order content of (2) is flat:

\[
\mathbb E_x\left(\mathbb E_yQ_A(xy)\right)^2
=\frac{N}{p^2}=\frac n2,
\]

and orthogonal projection onto the conditional mean gives

\[
\mathbb E_{x,y}\left(Q_A(xy)-\frac{Q_W(x)}p\right)^2
=N\left(1-\frac1{p^2}\right).
\tag{4}
\]

Both quantities are independent of \(A\). Thus the product transform, the
Paley two-frame, and second moments alone cannot recover Paley distance.

## 2. Far edge-local minima at the correct scale

For every fixed signing \(C\), there exists an edge-flip local minimum
\(A_\ast\) for \(\Phi\) such that

\[
\Phi(A_\ast)=O(n^{3/2}),\qquad
\Phi(A_\ast\circ C)=O(n^{3/2}),
\]

but

\[
\boxed{d_\pm(A_\ast,C)=\frac N2-O(n^{3/2})=\Theta(n^2).}
\tag{5}
\]

Set

\[
T_n=\sqrt{2N(n+3)\log2}.
\]

Choose the edges of \(A\) independently and uniformly. Both \(A\) and
\(A\circ C\) are marginally uniform signings. Hoeffding's inequality and a
union bound over the two signings and all Boolean \(x\) give positive
probability that

\[
\Phi(A),\Phi(A\circ C)\le T_n.
\]

By (1), initially \(d_\pm(A,C)\ge(N-T_n)/2\). Greedily flip an edge whenever
it strictly decreases \(\Phi\). All values of \(\Phi\) have the parity of
\(N\), so every step decreases it by at least \(2\), and at most \(T_n/2\)
steps occur. Each edge flip changes distance from a fixed orbit by at most
one. The terminal edge-local minimum therefore satisfies

\[
d_\pm(A_\ast,C)\ge\frac{N-T_n}{2}-\frac{T_n}{2}
=\frac N2-T_n.
\]

The same flips increase \(\Phi(A\circ C)\) by at most \(2\) apiece, so its
terminal value is at most \(2T_n\). After a best gauge,
\(W=A_\ast\circ C\) also satisfies every switching-minimal cut inequality

\[
0\le\sum_{e\in\delta(S)}W_e\le\Phi(W).
\]

Therefore edge-local minimality, switching/cut stability, Max-Lipschitz, and
the correct \(n^{3/2}\) scale do **not** imply \(d_\pm=o(n^2)\). Any surviving
stability proof must use exact global/cardinality minimality or a genuinely
higher-order invariant.

## 3. What closest global minimality adds

Let \(A\) be a global \(\Phi\)-minimizer chosen to minimize
\(d_\pm(A,C)=k\), and gauge it as

\[
A=C\oplus F,\qquad |F|=k.
\]

Write \(M=m_n\). For **every** nonempty \(H\subseteq F\), there is a signed
cut \(z_e=\sigma x_i x_j\) such that, with

\[
d_z=M-\langle A,z\rangle\in2\mathbb Z_{\ge0},
\]

one has

\[
\boxed{
\sum_{e\in H}C_ez_e\ge1+\frac{d_z}{2}.}
\tag{6}
\]

To prove this, flip \(H\) from \(A\) toward \(C\). The resulting signing is
strictly closer to the signed switching orbit, so it cannot also be a global
minimizer. Its \(\Phi\)-value is at least \(M+2\). Expanding a signed cut
that realizes that value gives (6).

In particular:

- for each \(e\in F\), an exact maximizer of \(A\) has \(C_ez_e=1\);
- for each pair \(e,f\in F\), a cut at level \(M\) or \(M-2\) has
  \(C_ez_e=C_fz_f=1\);
- if \(g=\Phi(C)-M\), an exact maximizer covers at most

  \[
  \frac{k}{2}+\frac g4\le\frac{p+1}{2p}k,
  \]

  where the last inequality uses Max-Lipschitz, \(g\le2k/p\).

This all-subsets/top-two-level witness hierarchy is the surviving Hamming
stability object. It is stronger than inclusion-minimality and edge-local
optimality, but no argument here proves that it forces \(k=o(n^2)\).

## 4. Permanent-gap counterexample shape

For an arbitrary signing \(A=C\oplus F\), put \(H=\Phi(A)\),
\(g=\Phi(C)-H\), and

\[
S_F(x)=\sum_{e=\{i,j\}\in F}C_ex_ix_j.
\]

The condition \(\Phi(A)\le H\) is equivalent to the global slab

\[
\frac{Q_C(x)-H}{2}\le S_F(x)\le\frac{Q_C(x)+H}{2}
\qquad\text{for every Boolean }x.
\tag{7}
\]

On \(\mathrm{Max}_+\) and \(\mathrm{Max}_-\), respectively, this gives

\[
g/2\le S_F\le\Phi(C)-g/2,
\qquad
-\Phi(C)+g/2\le S_F\le-g/2.
\]

Averaging the Paley frames yields \(k\ge pg/2\). Thus a fixed relative
improvement

\[
\Phi(A)\le(1/2-\varepsilon)n^{3/2}
\]

requires

\[
k\ge(\varepsilon/2+o(1))n^2.
\]

Hamming stability is therefore sufficient for \(L=1/2\), not equivalent to
it. A permanent-gap counterexample must be a dense signing satisfying (7) on
**every** Boolean direction, including all non-Max spike directions.

There is also a conference-class warning. Spectral defect zero identifies
the union of all conference switching classes, not the chosen Paley class.
When \(n=p^2+1\), a conference class has \(\rho=1\) exactly when it can be
switched to constant row sum \(p\) or \(-p\). Craigen's regularity lemma
proves only the forward implication that regularity forces square \(n-1\);
the converse is not known. The classified square orders \(2,10,26\) are
regularizable, and order \(50\) is the first unresolved case. Thus a spectral
rigidity argument landing at an arbitrary square-order conference class
still needs a Boolean-radius or regularizability theorem.

Mathon's construction supplies symmetric conference matrices of order
\(5r^2+1\) for prime powers \(r\equiv3\pmod4\). Their eigenvalues
\(\pm r\sqrt5\) are irrational, so they have no Boolean eigenvector and hence
\(\rho<1\) at each order. Irrationality gives no uniform gap. More sharply,
Momihara--Suda's maximum-excess bound gives, for odd
\(k\le r\sqrt5<k+2\),
\[
\rho(C_r)\le
\frac{k^2+2k+5r^2}{2(k+1)r\sqrt5},
\qquad
1-U_r=\frac{1-(k+1-r\sqrt5)^2}{2(k+1)r\sqrt5}.
\]
This upper bound itself approaches one, sometimes at the
\(\Theta(r^{-3})=\Theta(n^{-3/2})\) scale along the Pell/Fibonacci
approximants. Thus integrality, parity, irrationality, nonregularizability,
and the general excess theorem cannot produce a fixed gap.

Mathon's six-dimensional block quotient instead gives the lower bound
\[
\rho(C_r)\ge
\frac{r(8r+2)}{(5r^2+1)\sqrt5}
\longrightarrow\frac8{5\sqrt5}.
\]
The missing quantitative theorem is uniform anti-flatness of both full
conference eigenspaces, equivalently switched-row variance
\(\Omega(r^2)\) for every Boolean switching. If one could prove
\(\rho\le1-\eta\) uniformly on the ratio-dense prime subfamily, that would
be a direct permanent-gap construction and would refute \(L=1/2\).

Primary sources: Rudolf Mathon, [*Symmetric Conference Matrices of Order
\(pq^2+1\)*](https://doi.org/10.4153/CJM-1978-029-1), Canadian Journal of
Mathematics 30 (1978), 321--331; and Momihara--Suda, Proposition 1.1,
arXiv:1611.01305.

## 5. Exact cut-code/free-energy bridge

In binary edge coordinates define the augmented cut code

\[
\mathcal D_n=\{\delta u+s\mathbf1_E:u\in\mathbb F_2^n,\ s\in\mathbb F_2\}.
\]

Its dual consists exactly of the even-size Eulerian subgraphs of \(K_n\):

\[
\mathcal D_n^\perp
=\{H\subseteq E(K_n):\deg_H(v)\equiv0\pmod2,\ |H|\equiv0\pmod2\}.
\]

For a signing \(a\), define the signed even-Eulerian polynomial

\[
P_a(t)=\sum_{H\in\mathcal D_n^\perp}
\left(\prod_{e\in H}a_e\right)t^{|H|}.
\]

Expanding the exponential and averaging over the Boolean cube gives the exact
high-temperature identity

\[
\boxed{
\mathbb E_x\cosh(\beta Q_a(x))
=(\cosh\beta)^N P_a(\tanh\beta).}
\tag{8}
\]

Since \(\mathbb E\cosh(\beta Q_a)\le e^{\beta\Phi(a)}\), (8) shows that the
following initially proposed estimate would settle the lower bound:

\[
\boxed{
\inf_a\log P_a\!\left(\tanh\frac2{\sqrt n}\right)\ge-o(n).}
\tag{9}
\]

It is, however, false. For \(0<\theta<1\), concavity and
\(\cosh(u)^\theta\le e^{\theta u}+e^{-\theta u}\) give
\[
\mathbb E_a\!\left[\mathbb E_x\cosh(\beta Q_a(x))\right]^\theta
\le2^{n(1-\theta)+1}\cosh(\theta\beta)^N.
\]
Hence some deterministic signing satisfies
\[
\log P_a(\tanh\beta)
\le\frac{[n(1-\theta)+1]\log2}{\theta}
+\frac N\theta\log\cosh(\theta\beta)-N\log\cosh\beta.
\tag{10}
\]
For \(\beta=c/\sqrt n\), optimizing at
\(\theta=2\sqrt{\log2}/c\) yields
\[
\inf_a\log P_a(\tanh(c/\sqrt n))
\le-\left(\frac c2-\sqrt{\log2}\right)^2n+o(n).
\tag{11}
\]
At \(c=2\), (11) contradicts (9) by a linear margin
\((1-\sqrt{\log2})^2n\).

The corrected fixed-\(c\) sufficient target is
\[
\boxed{
\inf_a\log P_a(\tanh(c/\sqrt n))
\ge\left(\frac c2-\frac{c^2}{4}\right)n-o(n).}
\tag{12}
\]
Indeed,

\[
\Phi(a)\ge
\frac{N\log\cosh(c/\sqrt n)+
\log P_a(\tanh(c/\sqrt n))}{c/\sqrt n}
=\left(\frac12-o(1)\right)n^{3/2}.
\]

The fractional-moment construction rules out (12) for
\[
c<\frac{\log2}{\sqrt{\log2}-1/2}=2.0843108\ldots.
\]
The clean surviving target is \(c=3\):
\[
\inf_a\log P_a(\tanh(3/\sqrt n))\ge-\frac34n-o(n).
\]

Classical covering-radius information does not approach this scale. The dual
already contains a broad range of Eulerian weights, so Delsarte external
distance loses the \(n^{3/2}\) term. Nor can the ordinary fixed-\(L^q\)
moment-norm lower bound see it:
for each fixed \(q\), degree-two hypercontractivity gives
\(\|Q_a\|_q=O_q(n)\), whereas the target is order \(n^{3/2}\). A useful
cut-code attack therefore has to preserve the multiplicative phases among
many Eulerian shells, at a temperature above the fractional-moment barrier;
termwise absolute values destroy the critical cancellation.

## Consequence for the attack

The old local/product route is closed as a path. The replacement choices are:

1. exploit (6) using the exact global/cardinality-minimal quantifier;
2. prove the corrected signed-Eulerian free-energy bound (12), for example
   at \(c=3\), possibly through a high-level character/SOS hierarchy;
3. prove an eventual Paley character-transport estimate strong enough for the
   actual \(o(p^3)\) deficit target; or
4. on the adversarial side, establish a uniform Boolean-radius gap for a
   ratio-dense non-Paley conference family.

The exact all-prime gap-2 Max+/Max- architecture remains a separate, stronger
sufficient route. Nothing in this note flips an E(1) predicate or settles
\(L\).
