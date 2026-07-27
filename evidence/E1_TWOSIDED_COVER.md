# Two-sided Max-covers: necessary structure of undercutters

**Date:** 2026-07-27  
**Status:** Structural lemmas proved; census/SA certified. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Necessary conditions for undercut

Let \(C\) be \(\rho=1\) conference, \(\Phi=\tfrac12 np\), \(p=\sqrt{n-1}\). Let \(F\) be a flip set and \(A=C\oplus F\).

**Lemma U1 (Max\(_{+}\) cover).** If \(\Phi(A)<\Phi(C)\), then \(S_F(y)\ge1\) for every \(y\in\mathrm{Max}_{+}\).

*Proof.* \(Q_A(y)=\Phi-2S_F(y)\). If some \(y\) has \(S_F(y)\le0\), then \(|Q_A(y)|\ge\Phi\), so \(\Phi(A)\ge\Phi\). \(\square\)

**Lemma U2 (Max\(_{-}\) anti-cover).** If \(\Phi(A)<\Phi(C)\), then \(S_F(z)\le-1\) for every \(z\in\mathrm{Max}_{-}\).

*Proof.* \(Q_A(z)=-\Phi-2S_F(z)\). If some \(z\) has \(S_F(z)\ge0\), then \(|Q_A(z)|\ge\Phi\). Since \(S_F\) is integer-valued, \(S_F(z)<0\Rightarrow S_F(z)\le-1\). \(\square\)

A flip set satisfying both is a **two-sided Max-cover**. Every strict undercutter is a two-sided Max-cover (and must additionally control non-eigenvector spikes).

## LP bound

By the Max\(\pm\) frames (Prop 15.27),
\[
\mathbb E_{+}[S_F]=\frac{|F|}p,\qquad\mathbb E_{-}[S_F]=-\frac{|F|}p.
\]
Two-sided \(\Rightarrow |F|/p\ge1\Rightarrow|F|\ge p\). Equality requires \(S_F\equiv1\) on \(\mathrm{Max}_{+}\) and \(S_F\equiv-1\) on \(\mathrm{Max}_{-}\).

## Gap depth vs margin

If \(\Phi(A)\le\Phi-2t\) with \(t\ge1\) integer and the max on \(\mathrm{Max}_{\pm}\) realises \(\Phi(A)\), then necessarily \(S_F\ge t\) on \(\mathrm{Max}_{+}\) and \(S_F\le-t\) on \(\mathrm{Max}_{-}\), hence \(|F|\ge tp\).

In particular gap \(\ge4\) requires \(|F|\ge2p\).

## Certified facts

### \(p=3\) (\(n=10\))

| Object | Count | Notes |
|--------|------:|-------|
| Size-\(p\) Max\(_{+}\) covers | 405 | Prop 15.28; all \(\Phi\ge17\) |
| Size-\(p\) two-sided covers | **0** | all have \(\max_{\mathrm{Max}_{-}}S\ge1\) |
| Perfect matchings | 945 | |
| Two-sided perfect matchings | **144** | exactly the N10-S undercutters; all \(\Phi=13=m_{10}\) |
| Gap \(\Phi-m_{10}\) | 2 | \(=2t\) with \(t=1\), \(|F|=5>p\) |

### \(p=5\) (\(n=26\))

| Object | Result |
|--------|--------|
| Covering \(p\)-stars | 390; all spike \(\max_{\mathrm{Max}_{-}}S=p\); \(\Phi=75\) |
| Size-\(p\) two-sided | 0 among stars (Prop 15.28) |
| Perfect matching Max\(_{+}\) covers | **none found** (5k random + 20×15k SA; this file + `E1_MATCHING_MAXCOVER.md`) |
| Two-sided covers \(k\le20\) | SA ongoing / none in sparse MITM census |
| Undercut of \(\Phi=65\) | **none found** (intensive SA+MITM) |

## Consequence for E(1)

- Undercutters = two-sided Max-covers + spike control.
- Card-min two-sided covers have size \(\ge p\); size \(p\) is impossible for \(\infty\)-stars (Prop 15.28) and for all size-\(p\) sets at \(n=10\).
- At \(n=10\), minimum two-sided cover size is **5** (perfect matchings), gap exactly 2.
- If for all \(p\ge5\) one has either no undercutter or \(k_\star=O(n^{3/2})\), then \(L=\tfrac12\) by Prop 15.20d / Max-Lipschitz.
- **Not proved in general.**

## Not established

- Absence of undercutters for all \(p\ge5\)
- \(k_\star=O(n^{3/2})\) in general  
- \(\lim\alpha_n\) exists
