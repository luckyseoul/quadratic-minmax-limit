# Denseness-path package: \(L=\lim\alpha_n=\tfrac12\)

**Audience:** independent AI / human check (Paata AI-test). Use **this file only**.  
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit  
**Date:** 2026-08-15. Close is 15.272, not Aut-Schur and not Gsum.

## Verdict

**\(L=\lim\alpha_n=\tfrac12\) is proved** by sandwich + Paley denseness + E(1).  
E(1) is Type I ∧ residual (ii) ∧ bi-tight, each for every prime \(p\ge5\). Residual (ii) is ND-closed. Residual (i) Type I is dual-eq empty via \(G_+\succ0\) on \(\mathcal W_{++}^0\) (15.272), not the 15.216 \(K_4\) path and not Gsum. Bi-tight is empty for all \(p\ge5\) by 15.167 (mult\((\lambda_{\max})\ge d-1\), \(\lambda_{\min}(\Phi)\ge6\), \(L_*<2d\)).

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
| **E(1) / \(L=\tfrac12\)** | **Proved** |

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

For bad \(\mu\), unique pairs are the edges of \(K_m\). Every triple of good lines occurs as a \(k=3\) Max+ (A3 sawtooth + phase lock \(s_0+s_1+s_2\equiv-2\); \(M_3=C(m,3)p^2(p-1)\) matches the complete enum \(100,1176,24200\) at \(p=5,7,11\)). Each triple spans the \(2\)-plane \(\{x+y+z=0\}\) on its three edges (A5: distinct phases; Fejer \(\neq0\)). The line graph of \(K_m\) is connected for \(m\ge3\), so \(\chi_e-\chi_f\) generate \(1^\perp\). That is the bad-\(\mu\) isotypic.

---

## Lemma E (Johnson same-line).

\(k=1\) Max+ are \(z=f\circ L\) for a good form \(L\) and every \(f=2\,1_S-1\), \(|S|=m\) (15.269 B). On a good \(\mu\)-line, \(P_S(k)=\hat f(k)\hat f(\alpha-k)\) is even and lies in \(\sum P=0\). If \(\sum_k c_k P_S(k)=0\) for every such \(S\) with \(c_k=c_{\alpha-k}\), then \(g^\top Bg=0\) on all Johnson \(g\), where
\[
B_{xy}=\omega^{\alpha y}\,\hat c(x-y).
\]
Evenness \(\hat c(d)=\omega^{\alpha d}\hat c(-d)\) makes \(B\) symmetric. A 1-swap plus the subset-sum lemma (\(WW^\top|_{1^\perp}\) has eigenvalue \(C(n-2,k-2)(n-k)/(k-1)>0\) at \(n=p-2\), \(k=(p-1)/2\)) gives off-diagonal form \(\beta+u_i+u_j\). The quadratic of the rank-1 piece is \(2\sigma(u\cdot g)\) with \(\sigma=2m-p=1\neq0\), so \(u\) is constant and \(B=D+\beta(J-I)\). Then \(B_{x,x+d}=\omega^{\alpha(x+d)}\hat c(-d)\) is \(x\)-independent iff \(\hat c=0\) off \(0\) (good \(\mu\) has \(\alpha\neq0\), so \(x\mapsto\omega^{\alpha x}\) is nonconstant). Hence \(\{P_S\}\) spans the \((m-1)\)-dimensional same-line hyperplane.

---

## Lemma F (good \(\mu\)).

Unique pairs \(= m\) same-line \(+ C(m-1,2)\) mixed \(= C(m,2)+1\); convolution ⇒ hyperplane dim \(C(m,2)\). Lemma E contributes dim \(m-1\). Complementary mixed: Lemma D on \(K_{m-1}\) for \(p\ge7\) (\(m-1\ge3\)) contributes \(C(m-1,2)-1\). Intersection \(0\), sum \(C(m,2)-1\). A through-\(L_0\) \(k=3\) triple has mixed Fejer \(\beta\neq0\), so \(\sum_{\mathrm{same}}wf=-2\beta\neq0\), the last direction. At \(p=5\), only \(k=1,3\) exist and the Veronese rank is \(65=\dim\mathcal W_{++}^0\).

---

## Lemma G (\(G_+\succ0\)).

Singer / Fejer / Weil / Gershgorin+DFT: the \(k=3\) Gram is PD on the Aut_∞-circulant space \(F\) for every prime \(p\ge7\) (at \(p=5\) the \(\mu=0\) block is singular). Lemmas D–F give \(F^\perp\)-injectivity of \(k=1\cup k=3\). Hence the Veronese spans \(\mathcal W_{++}^0\), so \(G_+\succ0\) (15.212). Aut-Schur is **not** used: Jacquet does not force the \(k=3\)-only space \(F\) to meet every irrep (\(p=5\) rank \(61/65\)).

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

## AI-test questions (use this file only)

**Q1 (residual i).**  
Is residual (i) (Type I freeness-fail dual-eq empty for all primes \(p\ge5\)) essentially proved by this package?

**Q2 (\(L\)).**  
Is \(\lim\alpha_n=\tfrac12\) proved by this package?
