# Size-\(p\) Max-covers: tight cover, affine-line stars, and spike

**Date:** 2026-07-27  
**Status:** Structural lemmas proved; full census certified at \(p=3\); star census at \(p=5\).  
**Existence of \(\lim\alpha_n\) remains OPEN.**

## Setup

Let \(C\) be Paley of order \(n=p^2+1\) (\(\rho=1\)), \(\Phi=\tfrac12 np\), \(p=\sqrt{n-1}\).
Write \(\mathrm{Max}_{\pm}=\{x\in\{\pm1\}^n:Cx=\pm px\}\) and, for an edge set \(F\),
\[
S_F(x)=\sum_{\{i,j\}\in F}C_{ij}\,x_ix_j.
\]
A **Max-cover** is an \(F\) with \(S_F(y)\ge1\) for every \(y\in\mathrm{Max}_{+}\). Every strict undercutter of \(\Phi(C)\) is a Max-cover (integer drop of \(Q\) by at least \(2\)).

Frame identities (Prop 15.27 + Max\(_{-}\) dual):
\[
\mathbb E_{\mathrm{Max}_{+}}[yy^\top]=I+\frac Cp,\qquad
\mathbb E_{\mathrm{Max}_{-}}[yy^\top]=I-\frac Cp.
\]
Consequently for every edge \(\{i,j\}\),
\[
\mathbb E_{+}[C_{ij}y_iy_j]=\frac1p,\qquad
\mathbb E_{-}[C_{ij}y_iy_j]=-\frac1p,
\]
so \(\mathbb E_{+}[S_F]=|F|/p\) and \(\mathbb E_{-}[S_F]=-|F|/p\).

## Lemma A (tight size-\(p\) cover)

If \(F\) is a Max-cover with \(|F|=p\), then \(S_F\equiv1\) on \(\mathrm{Max}_{+}\).

*Proof.* \(\mathbb E_{+}[S_F]=1\) and \(S_F\ge1\) pointwise force equality. \(\square\)

(The fractional Max-cover LP of Prop 15.27 has value exactly \(p\), so size \(p\) is the absolute minimum possible for any Max-cover.)

## Lemma B (spike from Max\(_{-}\))

If some \(z\in\mathrm{Max}_{-}\) has \(S_F(z)\ge1\), then
\[
\Phi(C\oplus F)\;\ge\;\Phi(C)+2,
\]
because \(Q_{C\oplus F}(z)=Q_C(z)-2S_F(z)=-\Phi-2S_F(z)\) and \(|Q|\ge\Phi+2\).

In particular, if \(S_F(z)=p\) then \(\Phi(C\oplus F)\ge\Phi+2p\).

## Certified census

### \(p=3\) (\(n=10\), \(\Phi=15\)): all size-\(p\) Max-covers

| Quantity | Value |
|----------|------:|
| Size-\(p\) Max-covers | **405** |
| Of which stars | 60 |
| Non-stars | 345 |
| \(\max_{\mathrm{Max}_{-}}S_F\) histogram | \(\{1{:}330,\;3{:}75\}\) |
| \(\Phi(C\oplus F)\) histogram | \(\{17{:}330,\;21{:}75\}\) |
| Undercuts of \(\Phi=15\) | **0** |

Every size-\(p\) Max-cover has \(\max_{\mathrm{Max}_{-}}S_F\ge1\), hence \(\Phi\ge17=\Phi+2\) by Lemma B. Card-min Max-covers never undercut at \(n=10\); undercutters first appear at \(k=5\) (matchings, N10-S).

### \(p=5\) (\(n=26\), \(\Phi=65\)): covering \(p\)-stars

| Quantity | Value |
|----------|------:|
| Covering \(p\)-stars | **390** |
| \(S_F\) on \(\mathrm{Max}_{+}\) | \(\equiv1\) (all 390) |
| \(\max_{\mathrm{Max}_{-}}S_F\) | \(=5=p\) (all 390) |
| Exact \(\Phi\) sample (10 of 390, MITM) | all **75** \(=\Phi+2p\) |
| Max\(_{-}\) frame \(I-C/p\) | certified |

### Affine-line geometry (star at \(\infty\))

Vertices \(\{\infty\}\cup\mathbb F_{p^2}\cong\{\infty\}\cup\mathrm{AG}(2,p)\). A \(p\)-star at \(\infty\) has leaf set \(L\subset\mathbb F_{p^2}\), \(|L|=p\).

| \(p\) | Covering leaves at \(\infty\) | Affine lines in \(\mathrm{AG}(2,p)\) |
|------:|------------------------------:|------------------------------------:|
| 3 | **6** | 12 |
| 5 | **15** | 30 |

Every covering leaf set at \(\infty\) is an affine line; exactly half of all affine lines cover (one character class of directions). Spike vectors \(z\in\mathrm{Max}_{-}\) achieving \(S_F(z)=p\) are constant on \(L\cup\{\infty\}\).

JSON: `evidence/e1_size_p_maxcover.json`.

## Proposition 15.28 (shipped in `solution.md`)

1. **(Proved)** Size-\(p\) Max-cover \(\Rightarrow S_F\equiv1\) on \(\mathrm{Max}_{+}\).
2. **(Proved)** \(\max_{\mathrm{Max}_{-}}S_F\ge1\Rightarrow\Phi(C\oplus F)\ge\Phi+2\).
3. **(Proved for Paley \(\infty\)-stars)** Affine line of nonsquare direction \(\Rightarrow\) tight Max-cover with \(\mathbb E_{-}[S^2]=2p-1\neq1\Rightarrow\max_{\mathrm{Max}_{-}}S\ge1\Rightarrow\Phi\ge\Phi+2\). Square directions are non-covers (\(S\equiv-1\) on \(\mathrm{Max}_{+}\)). Uses \(\chi|_{\mathbb F_p^\times}\equiv1\) and adjacent-edge frame moments only.
4. **(Certified \(p=3,5\))** Every covering \(p\)-star (any centre) has \(\max_{\mathrm{Max}_{-}}S_F=p\), hence \(\Phi=\Phi+2p\).
5. **(Certified \(p=3\))** All 405 size-\(p\) Max-covers have \(\mathbb E_{-}[S^2]\ge7/3>1\) and \(\Phi\ge\Phi+2\).

**Still open for E(1):** generalise (4)–(5) to all \(p\); control larger Max-covers; or exact optimality for \(p\ge5\).

## Relevance

- Card-min Max-covers cannot undercut on the certified range: they spike.
- At \(n=10\), undercut requires \(k=5>p=3\) (matchings).
- At \(n=26\), no matching Max-cover is known (`E1_MATCHING_MAXCOVER.md`); covering stars all spike to 75.
- Does **not** by itself force \(m_n=\Phi(C)\) or \(k_\star=O(n^{3/2})\). **Do not mark Main Theorem settled.**
