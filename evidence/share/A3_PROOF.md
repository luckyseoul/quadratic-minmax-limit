# A3: every triple of good lines is a \(k=3\) Max+

Stand-alone existence and 2-plane writeup for Lemma D
(`evidence/share/denseness_path_package.md`). No census is used as a
proof. The count \(M_3=C(m,3)\,p^2(p-1)\) is recorded only as a check.
Live construction: `src/e1_gmin_m4_prop15276.py`.

**Priors used, not reproved.** Lemma B / 15.269 B / 15.272
`theorem_Vplus_omega_support`: \(v\in V_+\) iff \(\operatorname{supp}\hat z\subseteq\{0\}\cup\Omega\)
and \(\hat z(0)=p\,v_\infty\). Fourier transform is an isomorphism
\(V_+\cong\mathbb C^{\{0\}\cup\Omega}\). Paley \(C\) is the conference
matrix on \(\{\infty\}\cup\mathbb F_q\), \(q=p^2\), \(p\ge5\) prime.
Good \(F_p\)-forms are those whose dual punctured line lies in \(\Omega\).
There are \(m=(p+1)/2\) of them. 15.271 A5: edge-product phase vectors
on a triangle are \((c_1,c_2)\), \((c_1-c_3,-c_3)\), \((-c_3,c_2-c_3)\);
any two equal forces some \(c_i=0\).

Throughout, representatives in \(\mathbb F_p\) are \(\{0,1,\ldots,p-1\}\).

---

## 1. Occupancy constraint from three-line support

Let \(L_0,L_1,L_2\) be three distinct good \(F_p\)-forms
\(\mathbb F_q\to\mathbb F_p\). Any two are linearly independent, so
\((L_0,L_1):\mathbb F_q\cong\mathbb F_p^2\) and
\(L_2=\alpha L_0+\beta L_1\) with \(\beta\neq0\) and \(L_2\) not a
scalar multiple of \(L_0\) or \(L_1\).

If a boolean \(z:\mathbb F_q\to\{\pm1\}\) has Fourier support on the
three dual lines of \(L_0,L_1,L_2\), then in coordinates
\((x,y)=(L_0(u),L_1(u))\)
\[
z(x,y)=A(x)+B(y)+C(\alpha x+\beta y)
\]
for three functions \(\mathbb F_p\to\mathbb C\). The fibre occupancy
\[
N_i(t)=\#\{u:L_i(u)=t,\,z(u)=+1\}\in\{0,\ldots,p\}
\]
satisfies, writing \(G_i(t)=2N_i(t)-p=\sum_{L_i=t}z\),
\[
G_0(x)+G_1(y)+G_2(\alpha x+\beta y)
=p\bigl(z(x,y)+2(\mu_A+\mu_B+\mu_C)\bigr)
\]
with \(\mu\) the means. The global mean of \(z\) is \(1/p\) as soon as
\(\hat z(0)=p\), so
\[
N_0(x)+N_1(y)+N_2(\alpha x+\beta y)\in\{p+1,\,2p+1\}
\]
and \(\sum_t N_i(t)=p(p+1)/2\). Conversely, any triple of occupancies
obeying the two-valued sum defines a boolean
\[
z(u)=\begin{cases}
+1 & \text{if }\sum_i N_i(L_i(u))=2p+1,\\
-1 & \text{if }\sum_i N_i(L_i(u))=p+1,
\end{cases}
\]
which is a sum of three 1-variable functions (Section 4) and therefore
has three-line Fourier support.

---

## 2. Steps, sumsets, and 2-term APs

Write \(N,M,K\) for the three occupancies, and
\[
\Delta N(x)=N(x+1)-N(x)\in\mathbb Z\cap[-p,p],
\qquad
\Delta_\alpha K(w)=K(w+\alpha)-K(w).
\]
The two-valued constraint implies
\[
\Delta N(x)+\Delta_\alpha K(\alpha x+\beta y)\in U:=\{0,\pm p\}
\]
for every \(y\), hence for every \(w\in\mathbb F_p\). Let \(A=\operatorname{im}\Delta N\)
and \(B=\operatorname{im}\Delta_\alpha K\). Then \(A+B\subseteq U\), so
\(|A+B|\le3\).

**Case I.** \(A\subseteq U\) and \(B\subseteq U\). Each step of \(N\) is
\(0\) or \(\pm p\). The sum of all \(p\) steps is \(0\), and
\(N\) takes values in \(\{0,\ldots,p\}\), so \(N\) is either constant
or \(\{0,p\}\)-valued. Constant occupancy means \(A\) is constant, so
\(z\) is a 2-line function: \(A(x)+B(y)\in\{\pm1\}\) for all \(x,y\)
forces one summand constant, hence one-line (a cylinder). A
\(\{0,p\}\)-valued occupancy is likewise a cylinder (\(k=1\)). Neither
is a genuine 3-line Max+. Case I is impossible for \(k=3\).

**Case II.** Some \(a\in A\setminus U\). Then
\(B\subseteq(U-a)\cap[-p,p]\). For \(a\in\{1,\ldots,p-1\}\) the only
legal steps in that set are \(\{-a,\,p-a\}\); for
\(a\in\{-(p-1),\ldots,-1\}\) they are \(\{-a,\,-p-a\}\). In either
subcase \(B\) is a 2-term arithmetic progression of common difference
\(\pm p\). The same argument with the roles of \(N\) and \(K\) reversed
gives that \(A\) is a 2-term AP of the same difference. Write
\(A=\{\lambda,\,\lambda-p\}\) with \(\lambda\in\mathbb F_p^*\). Then
\[
N(x+1)\equiv N(x)+\lambda\pmod p,
\]
so \(N(x)\equiv N(0)+\lambda x\pmod p\).

---

## 3. The only legal lift is the sawtooth

Among integer-valued lifts of an affine residue class
\(N(x)\equiv\mu x+s\pmod p\) with values in \(\{0,\ldots,p\}\) and
\(\sum N=p(p+1)/2\), the unique solution is
\[
N(x)=1+(\lambda x+s)\bmod p,\qquad\lambda\in\mathbb F_p^*,\ s\in\mathbb F_p.
\]
(The lift \(0+(\lambda x+s)\bmod p\) has sum \(p(p-1)/2\), short by
\(p\). The lift that uses \(p\) in place of \(0\) and \(0\) in place of
\(p\) is the same function as \(1+(\,\cdot\,)\bmod p\).)

Check: as \(x\) runs through \(\mathbb F_p\), \((\lambda x+s)\bmod p\)
is a permutation of \(\{0,\ldots,p-1\}\), so
\[
\sum_x N(x)=p+\sum_{k=0}^{p-1}k=p+\frac{p(p-1)}{2}=\frac{p(p+1)}{2}.
\]
The same \(\lambda\) (up to the Plücker scale of the form) appears for
\(M\) and \(K\), because \(A\) and \(B\) share the common difference
\(\pm p\) and the linear relation \(L_2=\alpha L_0+\beta L_1\) couples
the slopes.

Let \(\alpha=(\alpha_0,\alpha_1,\alpha_2)\in\mathbb F_p^3\) be a
Plücker kernel vector:
\(\alpha_0 L_0+\alpha_1 L_1+\alpha_2 L_2\equiv0\), all
\(\alpha_i\neq0\) (any two forms are independent). The occupancies are
\[
N_i(t)=1+(\lambda\alpha_i\,t+s_i)\bmod p.
\]

---

## 4. Phase lock \(s_0+s_1+s_2\equiv-2\pmod p\)

At a point \(u\in\mathbb F_q\) write
\(r_i=(\lambda\alpha_i L_i(u)+s_i)\bmod p\in\{0,\ldots,p-1\}\). Then
\[
\sum_i N_i(L_i(u))=3+\sum_i r_i.
\]
The linear combination \(\sum_i\lambda\alpha_i L_i(u)=0\), so
\(\sum_i r_i\equiv\sum_i s_i\pmod p\). The occupancy sum lies in
\(\{p+1,2p+1\}\equiv1\pmod p\), hence
\[
3+\sum_i s_i\equiv1\pmod p\qquad\Rightarrow\qquad
s_0+s_1+s_2\equiv-2\pmod p.
\]
(The same identity is 15.271 A2.) Thus \(s_0,s_1\) are free
(\(p^2\) choices) and \(s_2\equiv-2-s_0-s_1\) is determined. The
common fibre scale \(\lambda\) runs through \(\mathbb F_p^*\)
(\(p-1\) choices).

Because each \(r_i\in\{0,\ldots,p-1\}\), one has
\(\sum r_i\in\{0,\ldots,3p-3\}\) and \(\sum r_i\equiv-2\pmod p\), so
\(\sum r_i\in\{p-2,2p-2\}\) and the occupancy sum is automatically
one of \(\{p+1,2p+1\}\). The boolean \(z\) of Section 1 is therefore
defined at every \(u\).

---

## 5. Majority identity and three-line Fourier support

The threshold that defines \(z\) is the majority of the three
sawtooth occupancies: the midpoint of \(\{p+1,2p+1\}\) is
\((3p+2)/2\), and
\[
z=\frac{2(N_0+N_1+N_2)-(3p+2)}{p}\in\{\pm1\}.
\]
Equivalently, with \(r_i=N_i-1\),
\[
z=\frac2p\bigl(r_0\circ L_0+r_1\circ L_1+r_2\circ L_2\bigr)-3+\frac4p,
\]
a constant plus a sum of three 1-variable functions. Hence \(\hat z\)
is supported on the three dual lines of \(L_0,L_1,L_2\). Those lines
lie in \(\{0\}\cup\Omega\) because the forms are good.

The \(\pm1\) majority algebra (15.271 A) is the identity
\(\operatorname{maj}(s_1,s_2,s_3)=(s_1+s_2+s_3-s_1 s_2 s_3)/2\)
for \(s_i=\pm1\). The occupancy threshold is that same boolean
function of the three sawtooths.

---

## 6. \(\hat z(0)=p\), hence \(Cy=py\)

Each \(r_i\circ L_i\) takes every residue \(p\) times (fibres of size
\(p\), and \(r_i\) is a permutation of \(\{0,\ldots,p-1\}\) because
\(\lambda\alpha_i\neq0\)). Thus
\(\sum_u r_i(L_i(u))=p\cdot p(p-1)/2=p^2(p-1)/2\), and
\begin{align*}
\hat z(0)=\sum_u z(u)
&=\frac2p\cdot 3\cdot\frac{p^2(p-1)}{2}-3p^2+4p
=3p(p-1)-3p^2+4p=p.
\end{align*}
Lemma B therefore gives \(Cy=py\) for \(y=(1,z)\) on
\(\{\infty\}\cup\mathbb F_q\). (Individual frequencies on the three
lines need not all appear; the support statement is an inclusion.)

---

## 7. Every unordered triple occurs

The construction starts from an arbitrary triple of distinct good
forms. Every unordered triple therefore occurs as a \(k=3\) Max+.
The locked family attached to one triple has \(p^2(p-1)\) members
(\(s_0,s_1\) free, \(\lambda\in\mathbb F_p^*\)). Over all triples
\[
M_3=C(m,3)\,p^2(p-1).
\]
This equals \(100,1176,24200\) at \(p=5,7,11\), matching the complete
maj-3 enumeration. That numerical match is a check, not the existence
argument above.

---

## 8. Sawtooth DFT (Fejer kernel)

Let \(\varphi(t)=(\mu t+s)\bmod p\) with \(\mu\in\mathbb F_p^*\) and
\(G(t)=2N(t)-p=2\varphi(t)-(p-2)\). The one-variable DFT
\(\widehat G(k)=\sum_{t\in\mathbb F_p}G(t)\,\omega^{kt}\),
\(\omega=\exp(2\pi i/p)\), evaluates as
\[
\widehat G(0)=p,
\qquad
\widehat G(k)
=F_{\mu,s}(k)
:=\frac{2p\,\omega^{-k\mu^{-1}s}}{\omega^{k\mu^{-1}}-1}
\quad(k\neq0).
\]
(The arithmetico-geometric sum
\(\sum_{u=0}^{p-1}u\,\omega^{cu}=p/(\omega^c-1)\) for \(c\neq0\),
followed by the substitution \(u=\mu t+s\).) The denominator is
nonzero: \(\omega^c=1\) iff \(p\mid c\). Thus the sawtooth DFT never
vanishes off frequency \(0\) (Fejer).

If \(\xi\) lies on the dual of \(L_i\), so
\(\operatorname{Tr}(\xi u)=\kappa L_i(u)\) for a unique
\(\kappa\in\mathbb F_p\), the fibre-sum theorem (15.271 A3) gives
\(\hat z(\xi)=\widehat G_i(\kappa)\). In particular \(\kappa=0\) iff
\(\xi=0\), which is not an edge frequency.

---

## 9. Edge amplitudes, as a function of \((c_1,c_2,c_3,s_0,s_1,s_2,\lambda)\)

A locked triple \(T\) has three edges \(E(T)\). For a frequency triple
\((c_1,c_2,c_3)\in(\mathbb F_p^*)^3\) the A5 edge frequencies are
\((c_1,c_2)\), \((c_1-c_3,-c_3)\), \((-c_3,c_2-c_3)\). The edge-product
of two sawtooth DFTs (times the dual character that identifies the
line with \(\mathbb F_p\)) is
\begin{align*}
A_{01}
&=F_{\lambda,s_0}(c_1)\,F_{\lambda,s_1}(c_2),\\
A_{02}
&=F_{\lambda,s_0}(c_1-c_3)\,F_{\lambda,s_2}(-c_3),\\
A_{12}
&=F_{\lambda,s_1}(-c_3)\,F_{\lambda,s_2}(c_2-c_3).
\end{align*}
Explicitly,
\[
F_{\lambda,s}(c)
=\frac{2p\,\omega^{-c\lambda^{-1}s}}{\omega^{c\lambda^{-1}}-1}
\qquad(c\neq0).
\]
The 3-vector of amplitudes is
\(\mathbf A(c,s,\lambda)=(A_{01},A_{02},A_{12})\).

**No \(c_i=0\).** An edge frequency \(c_i=0\) is the DC term on that
line, i.e. \(\xi=0\). For a bad \(\mu\) the unique pair on an edge
has both endpoints nonzero in \(\Omega\) (if \(\xi=0\) then
\(\mu-\xi=\mu\in\Omega\), contradicting badness). Fejer then says
every coordinate of \(\mathbf A\) is nonzero. Convolution on the
three-line support puts \(\mathbf A\) in the plane
\(\{x+y+z=0\}\subset\mathbb C^{E(T)}\). If two coordinates vanished
the third would vanish by that plane relation, contradicting Fejer.

**Not \(\mathbb C\)-parallel.** (not $\mathbb C$-parallel.) 15.271 A5: the three characters on
\(\mathbb F_p^2\) with phase vectors \((c_1,c_2)\),
\((c_1-c_3,-c_3)\), \((-c_3,c_2-c_3)\) are pairwise distinct (any
equality forces some \(c_i=0\)). Three distinct characters times
nonzero amplitudes are not \(\mathbb C\)-parallel as vectors in
\(\mathbb C^{E(T)}\) unless two amplitudes vanish, which they do not.

Independently of characters, the projective class of \(\mathbf A\)
depends on \((c_1,c_2,c_3)\). With \(s=0\), \(\lambda=1\),
\[
\frac{A_{01}}{A_{02}}
=\frac{F(c_1)F(c_2)}{F(c_1-c_3)F(-c_3)}
=\frac{(\omega^{c_1-c_3}-1)(\omega^{-c_3}-1)}{(\omega^{c_1}-1)(\omega^{c_2}-1)}.
\]
At \((c_1,c_2,c_3)=(1,2,3)\) versus \((1,2,4)\) (all six A5
frequencies nonzero for every prime \(p\ge5\)) the ratios are equal
iff \(F(-2)=F(-4)\) iff \(\omega^{-2}=\omega^{-4}\) iff \(p\mid 2\),
which is false. Hence two amplitude 3-vectors are not
\(\mathbb C\)-parallel, and they span the full 2-plane
\(\{x+y+z=0\}\).

(The same non-constancy is visible in the 2-homothety
\(F(2k)/F(k)\propto 1/(\omega^k+1)\), which depends on \(k\) because
\(\omega+1\neq\omega^2+1\). For odd \(p\), \(\omega^k=-1\) never
occurs, so no extra vanishing.)

A single locked phase orbit at fixed \((c_1,c_2,c_3)\) (or fixed bad
\(\mu\)) only multiplies \(\mathbf A\) by one character of the
\((s_0,s_1)\)-torus; that orbit is a ray in the plane. The second
direction is the dependence of the Fejer product on the frequency
triple. Rank-\(2\) at \(p=5,7,11\) is a check of this identity, not
the proof.

---

## 10. Line graph (already in Lemma D)

The line graph of \(K_m\) is connected for \(m\ge3\). The 2-plane of
any triangle containing a given pair of adjacent edges contains
\(\chi_e-\chi_f\). Therefore all differences \(\chi_e-\chi_{e_0}\)
lie in the span, which is \(1^\perp\) on \(E(K_m)\). That is the
bad-\(\mu\) isotypic (Lemma C).

---

## 11. What is not claimed

- Aut-Schur / Jacquet / PSL-span of the \(k=3\) Singer subspace.
- The cotangent pairing \(1^\top K^{-1}v\).
- That \(k=3\) alone fills \(\mathcal W_{++}^0\) (false at \(p=5\):
  rank \(61/65\)).
- Anything about bi-tight, residual (ii), or multi-level Type I.
