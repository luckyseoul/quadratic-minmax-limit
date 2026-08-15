# Denseness-path package: \(L=\lim\alpha_n=\tfrac12\)

**Audience:** independent AI / human check (Paata AI-test). Use **this file only**.  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-15. Close is 15.272, not Aut-Schur and not Gsum.

## Verdict

**15.272** writes a Max+-free spanning argument that \(k=1\cup k=3\) fills \(\mathcal W_{++}^0\) (Lemmas B–G), hence \(G_+\succ0\), hence dual-eq is empty on \(\mathrm{sc}\) (H–I). Aut-Schur is **false**. Gsum is unused.

**Do not send this to Paata as a finished \(L=1/2\) proof.** Independent hostile review (2026-08-15) found three load-bearing holes **outside** the 15.272 span, plus two writeup holes **inside** Lemma D. They are listed under **Public caveats**. Live code still gates E(1) as True; that is a wiring fact, not a substitute for those holes.

| Piece | Status |
|-------|--------|
| Sandwich \(1/\pi\le\liminf\le\limsup\le\tfrac12\) | **Proved** (`solution.md` §4–5) |
| \(\rho=1\) on Paley \(n=p^2+1\) | **Proved** (`evidence/PROOF_rho_eq_1.md`) |
| Denseness (Prop 6.1) | **Proved** (below) |
| Bi-tight empty, all \(p\ge5\) (15.167) | **Proved** |
| Residual (ii) ND | **Proved** (15.179+236+237) |
| Residual (i) Type I | **Proved** (15.272 → 15.207 → 15.249 → 15.216) |
| Aut-Schur / Jacquet | **False** (\(p=5\) \(k=3\) rank \(61/65\); unused) |
| Gsum disj LB | **False** / unused |
| Pairing \(1^\top K^{-1}v\) | **Open** / unused |
| **E(1) / \(L=\tfrac12\)** | **Not prize-ready** — see Public caveats |

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

**Bi-tight, all \(p\ge5\).** \(\mathrm{mult}(\lambda_{\max}(\Phi))\ge d-1\) (PSL min irrep). \(\lambda_{\min}(\Phi)\ge6\) because \(Q_4\) is a Gram Rayleigh of Max+ edge signs, hence \(Q_4\ge0\Rightarrow\mathrm{ray}\ge0\Rightarrow\lambda_{\min}\ge6\). Majorization then gives \(\lambda_{\max}\le L_*=(p^4+24p^2-1)/(2(p^2-1))\), and \(2d-L_*=(p^4-24p^2-1)/(2(p^2-1))>0\) for \(p\ge5\) (\(f(x)=x^2-24x-1\), \(f(25)=24\), \(f'>0\)). So \(\lambda_{\mathrm{cycle}}<d\) and bi-tight covers are empty.

**Residual (ii).** Affine two-level branch empty (15.179: forces \(k=3p-1\), impossible for \(k\ge3p\)). Even \(k\le4p-2\) Max− dichotomy (15.236). Dual-bad pair-span \(\{S=-4,f_e=-1\}\) cannot be a star, pair-slice, or triangle (15.237). Exhaustiveness that freeness-fail forces \(S\in\{2,4\}\) is **not** claimed and is not required.

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
This is \(100,1176,24200\) at \(p=5,7,11\), matching the complete maj-3 enum (a check, not the existence proof). The construction uses an arbitrary triple, so every unordered triple of good lines occurs.

**One triple spans the \(2\)-plane.** A locked triple contributes products on its three edges \(E(T)\), summing to \(0\) (convolution). Shifts with the lock give three Fejer-nonzero coordinates times three characters on \(\mathbb F_p^2\) with phase vectors \((c_1,c_2)\), \((c_1-c_3,-c_3)\), \((-c_3,c_2-c_3)\). Equality of any two forces some \(c_i=0\), contradicting Fejer (sawtooth DFT never vanishes off \(0\)). Hence the three characters are distinct, the three edge-vectors are not \(\mathbb C\)-parallel, and they span the full \(2\)-plane \(\{x+y+z=0\}\) on \(E(T)\).

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

## Lemma J (Type I / E(1) / \(L\)).

Lemma H + I ⇒ free-\(e\) max on \(\mathrm{sc}\) is \(<2-\alpha\) ⇒ dual-eq empty for every prime \(p\ge5\) (15.216 via 15.249). That closes Type I (15.170 dual-eq path; Gsum unused). Residual (ii) is ND-closed for every \(p\ge5\). Bi-tight is empty for every \(p\ge5\) (15.167: \(\mathrm{mult}(\lambda_{\max})\ge d-1\), \(\lambda_{\min}(\Phi)\ge6\) from \(Q_4\ge0\), and \(L_*<2d\)). Hence E(1) on the whole Paley family \(n=p^2+1\), hence \(\alpha_{p_k^2+1}\to\tfrac12\), hence \(L=\tfrac12\) by Lemma A. The live `e1` gate samples bi-tight at \(p=5\); the algebra is the same \(L_*<2d\) for all \(p\ge5\).

---

## Not used (do not revive)

- Aut-Schur / Jacquet / PSL-span of \(k=3\) \(F\).
- Cotangent pairing \(1^\top K^{-1}v\) (same-line rank \(m\); unused).
- Gsum disjoint lower bound (still False).
- Envelope / reflection / \(K_4\le\mathrm{Wick}_{hi}\) / \(\lvert\mu\rvert\le2/n\) (counterexamples).
- Path-C / \(16N\) (optional, independent).

Historical remarks “\(L\) OPEN” in Props 15.20–15.171 refer to those older routes and are not the current claim.

---

## Public caveats (do not skip)

A hostile review for public scrutiny found the following. **None of these is Aut-Schur or Gsum.**

1. **\(\lambda_{\min}(\Phi)\ge6\) is not proved for all \(p\).** 15.167 majorization needs this floor. The citation “\(Q_4=Be^\top G_{u,\mathrm{disj}}Be\) is a Gram” is **false**: \(G_{u,\mathrm{disj}}\) is a Hadamard mask of a Gram and has negative eigenvalues (\(\approx-30\) at \(p=3\), \(\approx-635\) at \(p=5\)). Floor certified only at \(p=5,7\) from exact \(\Phi\) spectra. Without a uniform floor, bi-tight need not be empty for \(p\ge11\), so E(1) is not proved on the whole Paley family.

2. **Residual (ii) ND is written for even \(k\in[3p+1,4p-2]\)** (15.236/237). Affine two-level is dead for all \(k\ge3p\) (15.179). Two-level non-affine at **\(k\ge4p\)** (where \(E_-[S]\le-4\) and \(N_2/N\le0\)) is not in those ND lemmas.

3. **Type I dual-eq is the two-level Max− law \(S\in\{-1,-3\}\).** The 15.169 bad case \(f_e\equiv-1\) on \(\{S=-1\}\) only gives \(\Phi(H)\ge\Phi-4\) if Max− is multi-level. Dual-eq empty does not by itself write ND for that multi-level bad case.

4. **Lemma D existence / 2-plane (writeup).** The sawtooth+lock count \(M_3=C(m,3)p^2(p-1)\) is a schema. A pointwise \(z=\mathrm{maj}\) with three-line Fourier support is not fully written in-repo (`A3_PROOF.md` is missing). Distinct Fejer characters on a triangle cannot all lie in a \(2\)-plane and stay non-parallel unless amplitudes depend on the phase; that dependence is not written. Rank-\(2\) at each triple is certified at \(p=5,7,11\) only.

5. **Lemma E Johnson** (same-line hyperplane) was expanded independently and has **no algebraic GAP** relative to 15.269 B + the \(WW^\top\) identity. See `evidence/share/lemma_E_johnson.md`.

Until (1)–(3) are closed, **do not claim \(L=1/2\) as a prize result.** Until (4) is written, do not treat bad-\(\mu\) spanning as a complete general proof (finite rank-2 + \(p=5\) Veronese 65 remain).

---

## AI-test questions (use this file only)

**Q1 (residual i).**  
Is residual (i) (Type I freeness-fail dual-eq empty for all primes \(p\ge5\)) essentially proved by this package?

**Q2 (\(L\)).**  
Is \(\lim\alpha_n=\tfrac12\) proved by this package?
