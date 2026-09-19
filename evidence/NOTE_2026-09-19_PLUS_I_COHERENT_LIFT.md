# Plus-I coherent lift: the n=26 Paley-beater is [[B, B+I],[B+I, −B]]

**Status:** identified construction; exact at n=10 and n=26; spectrum at n=50.
Does not prove \(\gamma_p\to c>0\) or \(\gamma_p\to0\). Original limit OPEN.

## Construction

Let \(B\) be a Seidel matrix of odd order \(m\). The *plus-I coherent lift* is
the order-\(2m\) Seidel matrix
\[
K(B)=\begin{pmatrix} B & B+I \\ B+I & -B \end{pmatrix}.
\]
The off-block diagonal of \(B+I\) is admissible (\(\pm1\)); the big diagonal
remains zero.

For \((x,y)\in\{\pm1\}^m\times\{\pm1\}^m\),
\[
Q_{K}(x,y)=Q_B(x)-Q_B(y)+x^TBy+x\cdot y.
\]

## Exact values on Paley orders \(n=p^2+1=2m\)

| p | m | B | Φ(B) | n=2m | Φ(K(B)) | Paley conference Φ | undercut |
|---|---|---|------|------|---------|-------------------|----------|
| 3 | 5 | 5-cycle (exact \(m_5=4\)) | 4 | 10 | **13** \(=m_{10}\) | 15 | yes, and optimal |
| 5 | 13 | exact \(m_{13}=20\) | 20 | 26 | **61** | 65 | yes; this *is* the campaign record |
| 5 | 9 | exact \(m_9=12\) | 12 | 18 | 41 | 33 | no (worse) |
| 3 | 9 | Paley Seidel of \(\mathbf F_9\) | 12 | 18 | **33** | 33 | tie |
| 7 | 25 | Paley Seidel of \(\mathbf F_{25}\) (record 60) | 60 | 50 | \(\ge169\), SDP \(\le180.3\), spectral \(\le25\sqrt{61}\approx195.3\) | 175 | unresolved |

The n=26 matrix in `ns_port_n26_undercut_20260912/stage01.npz` equals
\(K(B)\) for the order-13 block \(B=A[:13,:13]\), and \(B+I\) is the off-block.
So the “distant almost-conference” is this lift, not a matching of Paley-26.

## Spectrum when \(B\) has eigenvalues \(\pm p\) and \(0\)

Let \(v\neq0\) with \(Bv=\mu v\). The ansatz \((v,\alpha v)\) is an eigenvector of \(K(B)\) iff
\[
\alpha^2(\mu+1)+2\alpha\mu-(\mu+1)=0
\]
(when \(\mu\neq-1\)). For the Paley Seidel of \(\mathbf F_{25}\) one has
\(\mu\in\{5,-5,0\}\) (multiplicities 12, 12, 1). Then
\[
\operatorname{spec}(K)=\bigl\{\pm\sqrt{61}^{\,(12)},\ \pm\sqrt{41}^{\,(12)},\ \pm1\bigr\}.
\]
In particular \(\|K\|_{\mathrm{op}}=\sqrt{61}\) and
\(\Phi(K)\le\tfrac12\cdot50\cdot\sqrt{61}=25\sqrt{61}\). The same block
ansatz recovers the six-eigenvalue pattern of the n=26 lift
(\(\pm5.849^6,\pm4.448^6,\pm1\)).

The CORE Gaussian bound cannot see this undercut: \(K\) is not conference
(\(\|K\|_{\mathrm{op}}>\sqrt{n-1}\)).

## What this does not prove

The lift beats Paley conference at the two exact Paley orders \(n=10,26\).
It is not known to beat \(C_{50}\). A uniform \(c>0\) would require
\(\Phi(K(B_p))/n^{3/2}\le\tfrac12-c\) along \(n=p^2+1=2m\), which is not
proved. Using a smaller \(\Phi(B)\) can *worsen* the lift (n=18: \(m_9=12\)
lifts to 41>33).

Receipts: `evidence/coherent_BI_lift_20260919/best_n25.json`, n=26 npz,
scripts `lift_BI_phi_search.py`, `lift_BI_exact_pairs.py`.
Local-search lower bound 169 on \(K_{50}\) is not an upper bound.
