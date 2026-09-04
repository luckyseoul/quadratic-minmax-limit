# Denseness-path package (intended \(L=\tfrac12\); the limit is OPEN)

**Audience:** independent check of the denseness-path argument. Use **this file only**.  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-09-04. The residual-(i) two-level hinge is 15.272, not Aut-Schur and not Gsum. \(L=\tfrac12\) is **not proved**.

## Verdict

**15.272** writes a Max+-free spanning argument that \(k=1\cup k=3\) fills \(\mathcal W_{++}^0\) (Lemmas B–G), hence \(G_+\succ0\), hence dual-eq is empty on \(\mathrm{sc}\) (H–I). Aut-Schur is **false**. Gsum is unused.

Independent review found two remaining load-bearing predicates **outside**
the 15.272 span: residual (ii) and the minimal-four-gap implication bridge
exposed by Proposition 15.764. The former Lemma D writeup holes were closed by
15.276 and `A3_PROOF.md`; the former spectral/bi-tight hole is closed directly
by the 15.720 degree congruence, and the Type-I multi-level hole is closed by
15.750. The live remainder is listed under **Caveats**. The live and
legacy-named expanded acceptance gates now test the same predicates.

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** (`solution.md` §4–5) |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** (`evidence/PROOF_rho_eq_1.md`) |
| Denseness (Prop 6.1) | **Proved** (below) |
| Required bi-tight levels 2 and 3, all \(p\ge5\) | **Proved** (15.720; degree congruence + 15.272/15.207 kernel) |
| Residual (ii), affine + even \(k\le4p-2\) | **Proved** (15.179+236+237) |
| Residual (ii), even \(k\ge4p\) | **Open** — 15.734--15.749 close the first three all-prime shells and several later rows; see Caveat 2 for the exact remainder |
| Minimal-four-gap implication bridge | **Open** — 15.764 proves odd `|H|<=5p`; failure ranges start at even `|H|>=4p+2` and odd `|H|>=5p+2` |
| Residual (i) Type I, two-level Max− | **Proved** (15.272 → 15.207 → 15.249 → 15.216) |
| Residual (i) Type I, multi-level Max− | **Proved** (15.750) |
| Aut-Schur / Jacquet | **False** (\(p=5\) \(k=3\) rank \(61/65\); unused) |
| Gsum disj LB | **False** / unused |
| Pairing \(1^\top K^{-1}v\) | **Open** / unused |
| Lemma D existence / two-plane | **Proved** (15.276; `A3_PROOF.md`) |
| **E(1) / \(L=\tfrac12\)** | **Not settled** — see Caveats |

---

## Setup

\[
m_n=\min_{a_{ij}=\pm1}\max_{x=\pm1}\Bigl|\sum_{i<j}a_{ij}x_ix_j\Bigr|,\qquad
\alpha_n=m_n/n^{3/2},\qquad
L\stackrel{?}{=}\lim_n\alpha_n.
\]

Paley conference \(C\) of order \(n=p^2+1\) on \(\{\infty\}\cup\mathbb F_q\), \(q=p^2\), \(p\ge5\) prime. \(C_{\infty x}=1\), \(C_{xy}=\chi(x-y)\). \(\mathrm{Max}_+=\{y\in\{\pm1\}^n:Cy=py\}\). \(\Omega=\{\xi\in\mathbb F_q^\times:\hat\chi(\xi)=p\}\). \(m=(p+1)/2\). \(V_+=\{v:Cv=pv\}\). \(\mathcal W_{++}^0\) is the space of zero-diag \(++\) forms, dimension \(n(n-6)/8\).

E(1) on this family means \(m_n\ge\Phi(C)-2\) with \(\Phi(C)=\tfrac12 n\sqrt{n-1}\). Then \(\alpha_n\to\tfrac12\) along \(n=p^2+1\).

---

## Named priors (one paragraph each)

**Sandwich limsup.** A conference matrix of order \(n\) (Paley, or \(q\equiv1\pmod4\)) has \(\|C\|_{\mathrm{op}}=\sqrt{n-1}\), hence \(\Phi(C)=\tfrac12 n\sqrt{n-1}\) and \(\alpha_n\le\tfrac12\sqrt{1-1/n}\). Infinitely many such \(n\) with consecutive ratio \(\to1\), so \(\limsup\alpha_n\le\tfrac12\). (Lower sandwich \(\liminf\ge1/\pi\) is the dual-Gaussian arcsine bound, used only as a floor, not for the value \(1/2\).)

**\(\rho=1\).** On Paley \(n=p^2+1\), the halfspace boolean vector \(x_\infty=1\), \(x_u=\sigma(L(u))\) for an \(F_p\)-form \(L\) and \(S\subset F_p\) of size \(m\) satisfies \(Cx=px\) (`evidence/PROOF_rho_eq_1.md`: fibre character sums \(p-1\) on \(\ker L\) and \(-1\) off). Thus \(\rho(C)=1\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\).

**Bi-tight correction.** The 15.167 majorization arithmetic is conditional, but its final use of 15.55 is invalid: `ker(G-(n/2)P1)=span{1}+ker G`, not `span{1}`. Proposition 15.720 supplies the valid replacement. A centered level-\(s\) bi-tight indicator lies in `scheme+cross`; its degrees satisfy \(d_i+d_j\equiv2ps\pmod{(p^2-1)/2}\). The resulting common degree residue contradicts the handshake identity for required levels \(s=2,3\) at every prime \(p\ge5\). It also excludes bi-tight level 4, but not one-sided tight level 4. See Lemma K below.

**Residual (ii), proved range.** Affine two-level branch empty (15.179: forces \(k=3p-1\), impossible for \(k\ge3p\)). Even \(k\le4p-2\) Max− dichotomy (15.236). Dual-bad pair-span \(\{S=-4,f_e=-1\}\) cannot be a star, pair-slice, or triangle (15.237). Even \(k\ge4p\) is **not** in those ND lemmas (Caveat 2). Exhaustiveness that freeness-fail forces \(S\in\{2,4\}\) is not claimed.

---

## Lemma A (denseness).

If \(n_{k+1}/n_k\to1\), then \(\liminf_n\alpha_n=\liminf_k\alpha_{n_k}\) and likewise for \(\limsup\).

*Proof.* \(m_n\) is nondecreasing. For \(n_k\le N\le n_{k+1}\),
\[
\alpha_{n_k}(n_k/N)^{3/2}\le\alpha_N\le\alpha_{n_{k+1}}(n_{k+1}/N)^{3/2}.
\]
Both ratios tend to \(1\). \(\square\)

Along \(n_k=p_k^2+1\) one has \(n_{k+1}/n_k\to1\). Combined with \(\limsup\alpha_n\le\tfrac12\) and \(\alpha_{n_k}\to\tfrac12\) from E(1)+\(\rho=1\), one gets \(L=\tfrac12\).

---

## Lemma B (\(V_+\) Fourier).

Every \(v\in V_+\) has \(\operatorname{supp}\hat z\subseteq\{0\}\cup\Omega\) and \(\hat z(0)=p\,v_\infty\). The Fourier transform is an isomorphism \(V_+\cong\mathbb C^{\{0\}\cup\Omega}\).

*Proof.* Write \(v=(v_\infty,z)\) on \(\{\infty\}\cup\mathbb F_q\). \(Cv=pv\) is \(v_\infty+(\chi*z)=pz\) on \(\mathbb F_q\) and \(\sum z=p v_\infty\). FT at \(\xi\neq0\): \(\hat z(\xi)(\hat\chi(\xi)-p)=0\), so \(\operatorname{supp}\hat z\cap\mathbb F_q^\times\subseteq\Omega\). At \(0\): \(\hat z(0)=p v_\infty\). Individual vectors may omit frequencies. Both sides have dimension \((q+1)/2\), and FT on \(\mathbb F_q\) is injective, so the image is the full coordinate space. \(\square\)

(The \(\xi=0\) equation is **not** \((\hat\chi(0)-p)\hat z(0)=0\). \(\hat\chi(0)=0\) does not force \(\hat z(0)=0\).)

---

## Lemma C (isotypic = pair-hyperplane).

A \(++\) form is a symmetric kernel \(\Gamma\) on \((\{0\}\cup\Omega)^2\). Translation multiplies \(\widehat B(\xi,\eta)\) by \(\psi((\xi+\eta)t)\), so the \(\mu\)-isotypic is \(\Gamma\) on
\[
S_\mu=\{\xi\in\Omega\cup\{0\}:\mu-\xi\in\Omega\cup\{0\}\}.
\]
Symmetry ⇒ functions on unique unordered pairs. Zero diagonal:
\[
B_{xx}=\sum_\mu\psi(-\mu x)\sum_{\xi+\eta=\mu}\Gamma(\xi,\eta),
\]
hence \(B_{xx}\equiv0\) iff each pair-sum vanishes (convolution hyperplane). Unique-pair counts: \((p^2+7)/8\) (good \(\mu\)) and \(C(m,2)=(p^2-1)/8\) (bad). After the convolution cut,
\[
\dim(\text{good }\mu)=C(m,2),\qquad\dim(\text{bad }\mu)=C(m,2)-1,
\]
and \(\tfrac12(p^2-1)(2C(m,2)-1)=\dim F^\perp\). This is an isomorphism, not a dimension leap.

---

## Lemma D (bad \(\mu\)).

Write \(m=(p+1)/2\). The set \(\Omega\) is a union of \(m\) punctured \(F_p\)-lines (the *good* lines). A frequency \(\mu\neq0\) is **bad** when \(F_p\cdot\mu\cap\Omega=\emptyset\). Then \(S_\mu\) contains no same-line pair, so unique unordered pairs \(\{\xi,\mu-\xi\}\) are exactly the edges of the complete graph \(K_m\) on the \(m\) good lines: there are \(C(m,2)=(p^2-1)/8\) of them. (Check: \(|S_\mu|=(p^2-1)/4\) ordered, hence \(C(m,2)\) unique.)

**Every triple of good lines occurs.** Let \(L_1,L_2,L_3\) be any three distinct good \(F_p\)-forms. Any two are independent, so after \((L_1,L_2):\mathbb F_q\cong\mathbb F_p^2\) one has \(L_3=\alpha L_1+\beta L_2\) with \(\beta\neq0\). Occupancies \(N,M,K:\mathbb F_p\to\{0,\ldots,p\}\) of the three fibres satisfy \(N(x)+M(y)+K(\alpha x+\beta y)\in\{p+1,2p+1\}\) and \(\sum N=p(p+1)/2\). Steps of \(N\) lie in \(\{0,\pm p\}\); a constant occupancy is a \(2\)-line function, impossible for Max+. Integer sumsets with \(|A+B|\le3\) force both step-sets to be \(2\)-term APs of common difference \(\pm p\). The only legal reconstruction with the occupancy sum is the *sawtooth*
\[
N(x)=1+(\lambda x+s)\bmod p,\qquad\lambda\in\mathbb F_p^*,\ s\in\mathbb F_p
\]
(and likewise \(M,K\)). Boolean \(z=\mathrm{maj}\) of the three sawtooths has \(\hat z\) supported on the three dual \(\Omega\)-lines and \(\hat z(0)=p\), hence \(Cy=py\) by Lemma B. The Plücker relation \(\alpha_1 L_1+\alpha_2 L_2+\alpha_3 L_3\equiv0\) plus \(\sum N\equiv1\pmod p\) forces the phase lock \(s_0+s_1+s_2\equiv-2\pmod p\). Thus \(s_0,s_1\) are free (\(p^2\) choices) and \(s_2\) is determined; the common fibre scale \(\lambda\) runs through \(\mathbb F_p^*\) (\(p-1\) choices). Count:
\[
M_3=C(m,3)\,p^2(p-1).
\]
This is \(100,1176,24200\) at \(p=5,7,11\), matching the complete maj-3 enum (a check, not the existence proof). The construction uses an arbitrary triple, so every unordered triple of good lines occurs. Full occupancy-sum writeup: `evidence/share/A3_PROOF.md`. Live \(Cy=py\): `src/e1_gmin_m4_prop15276.py`.

**One triple spans the \(2\)-plane.** A locked triple contributes products on its three edges \(E(T)\), summing to \(0\) (convolution). The edge-product of two sawtooth DFTs is the Fejer kernel
\[
F_{\lambda,s}(c)=\frac{2p\,\omega^{-c\lambda^{-1}s}}{\omega^{c\lambda^{-1}}-1}\quad(c\neq0),
\]
and the amplitude 3-vector is
\(A_{01}=F_{\lambda,s_0}(c_1)F_{\lambda,s_1}(c_2)\),
\(A_{02}=F_{\lambda,s_0}(c_1-c_3)F_{\lambda,s_2}(-c_3)\),
\(A_{12}=F_{\lambda,s_1}(-c_3)F_{\lambda,s_2}(c_2-c_3)\).
No \(c_i=0\) (Fejer: \(\omega^k=1\) iff \(p\mid k\)). Shifts with the lock give these Fejer-nonzero coordinates times three characters on \(\mathbb F_p^2\) with phase vectors \((c_1,c_2)\), \((c_1-c_3,-c_3)\), \((-c_3,c_2-c_3)\). Equality of any two forces some \(c_i=0\), contradicting Fejer. Hence the three characters are distinct, the three edge-vectors are not \(\mathbb C\)-parallel (two amplitudes vanishing would force the third to vanish by \(x+y+z=0\)), and they span the full \(2\)-plane \(\{x+y+z=0\}\) on \(E(T)\). The ratio \(A_{01}/A_{02}\) depends on \((c_1,c_2,c_3)\) (\(F(-2)\neq F(-4)\)). Rank-\(2\) at \(p=5,7,11\) is a check.

**Line graph.** The line graph of \(K_m\) is connected for \(m\ge3\). The \(2\)-plane of any triangle containing a given pair of adjacent edges contains \(\chi_e-\chi_f\). Therefore all differences \(\chi_e-\chi_{e_0}\) lie in the span, which is \(1^\perp\) on \(E(K_m)\). That is the bad-\(\mu\) isotypic (Lemma C).

---

## Lemma E (Johnson same-line).

**Cylinders.** If \(z=f\circ L\) for a good form \(L\) and \(f=2\,1_S-1\) with \(|S|=m\), then \(\sum f=2m-p=1\), so \(\hat z(0)=p\sum f_{\mathrm{fibre}}=p\), and \(\hat z\) is supported on the dual line \(\subset\{0\}\cup\Omega\). Lemma B gives \(Cy=py\). Conversely a Max+ with one-line support is of this form, and every \(m\)-subset occurs. Count \(N_1=m\,C(p,m)\).

**Products.** On a good \(\mu\)-line identify the line with \(\mathbb F_p\) so \(\mu\) corresponds to \(\alpha\neq0\). Write \(\hat f(k)=\sum_x f(x)\omega^{kx}\). Then
\[
P_S(k)=\hat f(k)\hat f(\alpha-k),\qquad
\sum_k P_S(k)=\sum_x f(x)^2\,p\,\omega^{\alpha x}=p\sum_x\omega^{\alpha x}=0
\]
(\(\alpha\neq0\), \(f^2=1\)). Also \(P_S(k)=P_S(\alpha-k)\).

**Annihilator \(\Rightarrow\) quadratic form.** If \(\sum_k c_k P_S(k)=0\) for every \(m\)-subset \(S\) and \(c_k=c_{\alpha-k}\),
\begin{align*}
\sum_k c_k\hat f(k)\hat f(\alpha-k)
&=\sum_{x,y}f(x)f(y)\sum_k c_k\omega^{kx}\omega^{(\alpha-k)y}\\
&=\sum_{x,y}f(x)f(y)\,\omega^{\alpha y}\,\hat c(x-y)
=g^\top Bg,
\end{align*}
with \(B_{xy}=\omega^{\alpha y}\hat c(x-y)\) and \(g=f\).

**Evenness \(\Rightarrow\) symmetry.** \(\hat c(d)=\sum_k c_k\omega^{kd}\). The substitution \(k\mapsto\alpha-k\) and \(c_k=c_{\alpha-k}\) give \(\hat c(d)=\omega^{\alpha d}\hat c(-d)\). Then \(B_{yx}=\omega^{\alpha x}\hat c(y-x)\). With \(d=x-y\), \(y=x-d\),
\[
B_{xy}=\omega^{\alpha(x-d)}\hat c(d)=\omega^{\alpha x}\omega^{-\alpha d}\hat c(d)=\omega^{\alpha x}\hat c(-d)=B_{yx}.
\]

**1-swap \(\Rightarrow\) \(\beta+u_i+u_j\).** Let \(M=B\) (already symmetric). For \(a\in S\), \(b\notin S\), \(S'=(S\setminus\{a\})\cup\{b\}\), the identity \(g^\top Mg=0=(g')^\top Mg'\) and \(g'=g+2(e_b-e_a)\) imply that all \((m-1)\)-subset sums of \(\delta_i=M_{bi}-M_{ai}\) on \(X=\mathbb F_p\setminus\{a,b\}\) are constant. Here \(|X|=p-2\) and \(k=m-1=(p-1)/2\), so \(1<k<n\). The incidence \(W\) of \(k\)-subsets vs ground set has
\[
(WW^\top)_{ii}=C(n-1,k-1),\quad(WW^\top)_{ij}=C(n-2,k-2)\ (i\neq j),
\]
hence \(WW^\top|_{1^\perp}\) has eigenvalue \(C(n-1,k-1)-C(n-2,k-2)=C(n-2,k-1)>0\). So \(\delta\) is constant on \(X\): \(M_{bi}-M_{ai}\) is independent of \(i\neq a,b\). That is \(M_{ij}=\beta+u_i+u_j\) off-diagonal.

**\(\sigma=1\) kills \(u\).** For \(g\in\{\pm1\}^p\) with \(\sum g=\sigma\),
\[
g^\top(u1^\top+1u^\top)g=2\sigma(u\cdot g).
\]
Here \(\sigma=2m-p=1\neq0\), so \(u\cdot g\) constant on all Johnson \(g\) ⇒ \(u\) constant (same subset-sum lemma). Absorb \(u\) into \(\beta\). Diagonal is free because \(g_x^2=1\). Thus \(B=D+\beta(J-I)\).

**Fourier annihilator.** \(B_{x,x+d}=\omega^{\alpha(x+d)}\hat c(-d)\). Off-diagonals constant and \(\alpha\neq0\) force \(\hat c(-d)=0\) for \(d\neq0\). So \(c\) is constant — the convolution-hyperplane normal — and \(\{P_S\}\) spans the \((m-1)\)-dimensional same-line hyperplane.

---

## Lemma F (good \(\mu\)).

A good \(\mu\) lies on one good line \(L_0\). Unique pairs: \(m\) same-line pairs \(\{k,\alpha-k\}\) on \(L_0\), plus \(C(m-1,2)\) mixed pairs among the other good lines. Total \(m+C(m-1,2)=C(m,2)+1\). Convolution cuts one dimension, so the isotypic has dim \(C(m,2)\).

Lemma E contributes the same-line hyperplane (dim \(m-1\), mixed coordinates \(0\)).

**Complementary mixed.** The mixed pairs are \(E(K_{m-1})\) on the other good lines. For \(p\ge7\), \(m-1\ge3\). Relative to those \(m-1\) lines, \(\mu\) lies *off* them, so it is bad for that \(K_{m-1}\). Lemma D therefore spans \(1^\perp\) on the mixed edges (dim \(C(m-1,2)-1\), same-line coordinates \(0\)).

These two subspaces intersect at \(0\) (disjoint supports). Sum of dimensions: \((m-1)+(C(m-1,2)-1)=C(m,2)-1\), one short of the isotypic.

**Through-\(L_0\).** Take a \(k=3\) triple that includes \(L_0\) and two other good lines. Convolution on unique pairs of this good \(\mu\) reads \(2\sum_{\mathrm{same}}\Gamma+\sum_{\mathrm{mixed}}\Gamma=0\). The mixed coordinates of this triple are a Fejer product on an edge of \(K_{m-1}\) through a third line; A5+Fejer give a nonzero mixed amplitude \(\beta\). Then \(\sum_{\mathrm{same}}w f=-2\beta\neq0\): the same-line component is *not* in the hyperplane \(\sum P=0\) that Lemma E spans. Hence the vector is outside \(k{=}1\oplus\) mixed and supplies the last direction.

**\(p=5\).** Here \(m=3\), so \(m-1=2\) and \(K_2\) has no triangle: the complementary-mixed step is vacuous. Only \(k=1\) and \(k=3\) exist. The Veronese \(\{yy^\top-S:y\in\mathrm{Max}_+\}\) of that finite set has rank \(65=\dim\mathcal W_{++}^0=26\cdot20/8\) (full SVD; \(k=3\) alone is \(61/65\)). A complete rank computation at one prime is a proof at that prime.

---

## Lemma G (\(G_+\succ0\)).

\(\dim F=(p^2-5)/4\) (even mean-zero functions on \(\Omega\)). \(\dim\mathcal W_{++}^0=n(n-6)/8\). \(\dim F^\perp=(p^2-5)(p^2-1)/8\).

**Singer PD on \(F\), every prime \(p\ge7\).** The \(k=3\) locked-triple circulant on the \(m\) good lines:
- \(\mu=0\) block eigenvalues \(3(m-1)(m-2)/2\) and \((m-2)(m-3)/2\), both \(>0\) iff \(m\ge4\) i.e. \(p\ge7\). At \(p=5\), \(m=3\), the second eigenvalue is \(0\).
- Even-class Fejer stretch on each good line is a circulant with off-diagonal row-sum \((p^2-1)/8-|\hat s(1)|^2\) and \(|\hat s(1)|^2=1/(4\sin^2(\pi/(2p)))\). Diagonal dominance is \(\sin(\pi/(2p))<2/\sqrt{p^2-1}\), from \(\sin x<x\) and \(\pi\sqrt{p^2-1}<4p\) (\(p\ge3\)). So even stretches are invertible.
- Off \(\mu=0\): locked pair-sum \(g=\tfrac12(U+T)+\varepsilon_\infty\) with \(|T|\le2\sqrt{p}+1\) (Weil / Perel'muter on a squarefree quadratic of degree \(2\)) and \(|g|\le\sqrt{p}+2\). Gershgorin PD iff \(\sqrt{p}+2<(p-3)/4\), i.e. \(p-4\sqrt{p}-11>0\), which holds for every prime \(p\ge37\).
- Finitely many primes \(7\le p\le31\): direct DFT of the Singer Gram, no vanishing coefficient.

Thus \(G_{k=3}\succ0\) on \(F\) for all primes \(p\ge7\).

Lemmas D–F give \(F^\perp\)-injectivity of \(k=1\cup k=3\). Together with Singer PD on \(F\) (and the \(p=5\) rank-\(65\) fill), the Veronese spans \(\mathcal W_{++}^0\). By 15.212, \(G_+\succ0\) on \(\mathcal W_{++}^0\) iff that Veronese has rank \(\dim\mathcal W_{++}^0\).

Aut-Schur / Jacquet is **false** and unused: at \(p=5\) the \(k=3\)-only Veronese has rank \(61<65\).

---

## Lemma H (\(\ker=\mathrm{sc}\)).

15.207: \(\ker(\mathrm{Gsum})=\mathrm{scheme}\oplus\mathrm{cross}\) iff \(G_+\succ0\) on \(\mathcal W_{++}^0\). Proof: \(\ker Q_+=\) scheme-image iff the only zero-diag \(B=P_+BP_+\) with \(y^\top By=0\) on Max+ is \(B=0\), which is \(G_+\succ0\). Lemma G supplies that.

---

## Lemma I (\(\mathrm{cost}_D<2-\alpha\)).

After Comm + Comm(diag) repair (15.249, any conference): \(W_e=(p^4-4p^2+1)/(2p^2(p^2-2))\), \(\mathrm{sum}_{ne}^{(1)}=p^2(p^2+1)(p^2-2)/(p^4-4p^2+1)-1\). Paley Weil \(|Q_{ij}|\le 2p\) on far \(4\)-point sums gives \(t\le 2(2p-1)/\mathrm{den}\). Dual cost is \(\mathrm{cost}_D=\alpha\cdot\mathrm{sum}_{ne}\) with \(\alpha=6/(p(p^2+1))>0\). The Weil repair satisfies \(\mathrm{sum}_{ne}\le S_1+tN\), hence \(\mathrm{cost}_D\le\alpha(S_1+tN)\). The majorant identity
\[
(2-\alpha)-\alpha(S_1+tN)
=\frac{2(p^4-3p^3-10p^2+9p+1)}{p^4-4p^2+1}
\]
has positive denominator for \(p\ge3\). The numerator polynomial \(f(x)=x^4-3x^3-10x^2+9x+1\) has \(f(5)=46>0\) and \(f'(x)=4x^3-9x^2-20x+9>0\) for \(x\ge5\) (equivalently \(f'(x)\ge x(4x^2-9x-20)\) and \(4x^2-9x-20\ge35\)). Hence \(\mathrm{cost}_D<2-\alpha\) for every prime \(p\ge5\).

---

## Lemma J (what the 15.272 hinge actually closes).

Lemma H + I imply that the free-\(e\) maximum on `sc` is
\(<2-\alpha\), so dual equality is empty for every prime \(p\ge5\)
(15.216 via 15.249). This closes the **two-level** Type-I/residual-(i)
slice. It does not itself close the multi-level Type-I bad case; Proposition
15.750 later closes that case.

Propositions 15.179 and 15.236--15.237 initially prove the affine branch and
even \(k\le4p-2\). Propositions 15.734--15.749 later close the first three
all-prime higher shells and several later rows, but the exact non-Walsh
remainder in Caveat 2 remains open. The required bi-tight levels are
closed by Lemma K, so the spectral floor, global mixed-\(k\) QVAR, and
principal R1 are not needed by this implication chain. This package still
does **not** prove E(1) or \(L=\tfrac12\) until residual (ii) and the
minimal-four-gap implication bridge both close.

---

## Lemma K (required bi-tight levels).

Let `H` be bi-tight of level `s`, `|H|=sp`, with degree sequence `d_i`, and
let `kappa=1_H-2s/(np) 1`. Tightness on Max+ and Max− gives
`kappa in ker(Gsum)`. Lemmas G--H identify this kernel with
`scheme+cross`. Writing `B=C odot 1_H`, projecting `B` onto matrices
commuting with `C`, and comparing an off-diagonal entry gives

\[
d_i+d_j\equiv2ps\pmod{(p^2-1)/2}\qquad(i\ne j).
\]

Hence all degrees have one common residue modulo `M=(p^2-1)/2`. For `s=2`,
`M>2p` for every `p>=5`, so all degrees are equal, but the handshake identity
would give the noninteger `4p/(p^2+1)`. For `s=3`, the same argument works for
`p>=7`; at `p=5`, `M=12`, `n=26`, and total degree `30`, so a common residue
`r` must be `0` or `1`, and neither `30-26r` is divisible by `12`. Thus the
level-2 and level-3 bi-tight alternatives are empty for every prime `p>=5`.
The same residue calculation excludes bi-tight level 4, but this does not
assert generic Max+- or Max−-tight-cover emptiness. In particular it does not
close the one-sided level-4 branch in residual (ii), nor all bi-tight levels.

---

## Not used (do not revive)

- Aut-Schur / Jacquet / PSL-span of \(k=3\) \(F\).
- Cotangent pairing \(1^\top K^{-1}v\) (same-line rank \(m\); unused).
- Gsum disjoint lower bound (still False).
- Envelope / reflection / \(K_4\le\mathrm{Wick}_{hi}\) / \(\lvert\mu\rvert\le2/n\) (counterexamples).
- Path-C / \(16N\) (optional, independent).

Historical `CLOSED` claims in Props. 15.167--15.171, the older 15.272
writeup, and derived evidence use obsolete scope or the false floor premise.
The current claim is \(L\) **OPEN**.

---

## Caveats

A hostile review for public scrutiny found the following. **None of these is Aut-Schur or Gsum.**

1. **Spectral floor (not an acceptance caveat).** \(\lambda_{\min}(\Phi)\ge6\) remains open for all \(p\), and the old Kneser-mask proof is false. Proposition 15.720 bypasses this entirely for E(1), so QVAR/R1 work cannot be counted as closing a remaining gate.

2. **Residual (ii), exact current remainder.** Proposition 15.236/237 covers even \(k\in[3p+1,4p-2]\), and affine two-level is dead for all \(k\ge3p\) (15.179). At **\(k\ge4p\)**, Propositions 15.734--15.737 close the first three shells, Proposition 15.751 closes \(k=4p+6\) for every \(p\ge13\), Proposition 15.752 closes \(k=4p+8\) for every \(p\ge23\) plus its stated contiguous higher band, and Proposition 15.753 closes the p17/p19 fifth-shell endpoints. Propositions 15.744--15.749 close \(u=0,3,4\) at \(p=13,k=60\), and Proposition 15.754 closes the remaining \(u=6\) by an exact finite aggregate/common-form certificate. Thus the fifth shell is closed for every \(p\ge13\). The exact open remainder is: critical \(p=5,7\); \(p=11,k\ge50\); later \(p=13\) layers; layers beyond 15.752's band; and the separate positive \(p=7,z=7\) branch. The global residual predicate remains false.

   **Separate implication bridge.** Proposition 15.764 proves that an
   all-deletions minimal four-gap set enters the official residual-(ii) unit
   when `|H|` is odd and at most `5p`. It also identifies the exact remaining
   failure ranges: even `|H|>=4p+2` and odd `|H|>=5p+2`. This bridge is a
   separate acceptance predicate; closing residual (ii) alone cannot close E(1).

3. **Type I dual-eq is the two-level Max− law \(S\in\{-1,-3\}\).** The 15.169 bad case \(f_e\equiv-1\) on \(\{S=-1\}\) only gives \(\Phi(H)\ge\Phi-4\) if Max− is multi-level. 15.275 writes the mass \(2a+c(3+\mu_c)=2/p\), the pairing min \(E[Sf_e]\ge3/p-2\), \(E[R^2]=E[S^2]-5+4E[Sf]\), integrality \(n_{-1}=M+n_c+t\), the unique 2-orbit Aut\(_e\) collapse \(\mu_{\mathrm{far}}=-2(2p-3)/(p(p^2-1))<0\), and the 3-weight Max− identity \(F_-|_{f=+1}=-(p+1)/(p-1)+\mu_{\mathrm{far}}p(p+1)\). Paley Aut\(_e\) has **two** star orbits (\(\sigma_\square=(p-1)(1+f_e)\), \(\sigma_\boxtimes=(p+1)(1-f_e)\)) and several far orbits. The 3-weight slices \(\mu_\square=\mu_\boxtimes\) and \(\mu_\boxtimes=0\) are empty (negative weight). The slice \(\mu_{\mathrm{far}}=0\) stays a \([0,1]\) solution of \(F_+\bar x=3-2f_e\), but it (and the whole \(\mu_{\mathrm{far}}\ge0\) 3-weight family) **cannot realise the bad case**: \(F_-|_{f=+1}\ge-(p+1)/(p-1)>-2\), while a gap-2 undercutter with \(f_e\equiv-1\) on \(\{S=-1\}\) needs \(S\le-3\) on \(\{f=+1\}\). Star-supported 0-1 Type I graphs Aut\(_e\)-average to that point (\(n_\square=(p+1)/2\), \(n_\boxtimes=5(p-1)/2\)). Dual-eq empty kills only the two-level / pairing-min slice. **Historical leftover:** split far Aut\(_e\) classes (unequal 4-set interpolants) did not reduce to the collapsed-far bound. Proposition 15.750 supersedes that leftover and proves `type_I_multilevel_bad_case_ND_closed=True` for every prime \(p\ge5\).

4. **Lemma D existence / 2-plane — closed.** Written in `evidence/share/A3_PROOF.md` and checked live in 15.276: occupancy sumset \(\to\) sawtooth \(N(x)=1+(\lambda x+s)\bmod p\), majority \(z\), three-line support, \(\hat z(0)=p\), \(Cy=py\), phase lock \(s_0+s_1+s_2\equiv-2\). Amplitudes are the Fejer products \(F_{\lambda,s}(c)=2p\,\omega^{-c\lambda^{-1}s}/(\omega^{c\lambda^{-1}}-1)\), nonzero off \(0\), and the 3-vector is not \(\mathbb C\)-parallel as \((c_1,c_2,c_3)\) varies. \(M_3\) matching enum and rank-\(2\) at \(p=5,7,11\) remain checks, not the proof. This item is no longer a blocker.

5. **Lemma E Johnson** (same-line hyperplane) was expanded independently and has **no algebraic GAP** relative to 15.269 B + the \(WW^\top\) identity. See `evidence/share/lemma_E_johnson.md`.

Until Caveat (2), including its separate implication bridge, is closed,
\(L=\tfrac12\) is not established. Caveat (3)
is discharged by 15.750, and Caveat (4) by `A3_PROOF.md` + 15.276; finite
rank-2 remains an independent check of the latter proof.

---

## AI-test questions (use this file only)

**Q1 (residual i).**  
Is residual (i) (Type I freeness-fail dual-eq empty for all primes \(p\ge5\)) essentially proved by this package?

**Q2 (\(L\)).**  
Is \(\lim\alpha_n=\tfrac12\) proved by this package? (The honest answer is no until Caveat 2 closes.)
