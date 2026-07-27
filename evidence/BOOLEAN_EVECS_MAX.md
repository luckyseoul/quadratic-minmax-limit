# Maximizers of ρ=1 conference matrices are boolean eigenvectors

**Date:** 2026-07-27  
**Status:** Proved characterization; counts certified. Useful for E(1) covering analysis. Existence of \(\lim\alpha_n\) remains **OPEN**.

## Theorem

Let \(C\) be a symmetric conference matrix of order \(n\) with \(\rho(C)=1\), so \(\|C\|_{\mathrm{op}}=\sqrt{n-1}=:p\) and \(\Phi(C)=\tfrac12 np\). Write \(P_\pm\) for the orthogonal projectors onto the \(\pm p\) eigenspaces (each of dimension \(n/2\)). Then for every \(x\in\{\pm1\}^n\),
\[
x^\top C x
=
p\bigl(2\|P_+x\|_2^2-n\bigr),
\]
and therefore
\[
\bigl|x^\top C x\bigr|=np
\quad\Longleftrightarrow\quad
x\in\operatorname{range}(P_+)\ \text{or}\ x\in\operatorname{range}(P_-),
\]
i.e. \(Cx=\pm px\). Equivalently: the maximizers of \(\lvert Q_C\rvert=\Phi(C)\) are exactly the **boolean eigenvectors** of \(C\) for eigenvalues \(\pm p\).

*Proof.* \(C=pP_+-pP_-\) and \(P_++P_-=I\), so
\[
x^\top Cx=p\bigl(\|P_+x\|^2-\|P_-x\|^2\bigr)=p\bigl(2\|P_+x\|^2-\|x\|^2\bigr).
\]
For \(x\in\{\pm1\}^n\) one has \(\|x\|^2=n\). Thus \(\lvert x^\top Cx\rvert=np\) iff \(\|P_+x\|^2\in\{0,n\}\) iff \(x\) lies in one of the two eigenspaces. On \(\operatorname{range}(P_+)\) one has \(C=pI\), so \(Cx=px\). \(\square\)

## Certified counts (Paley \(n=p^2+1\))

Boolean vectors in the \(+p\) eigenspace (including \(\{x,-x\}\) pairs), by exact free-variable enumeration of \(\ker(C-pI)\):

| \(p\) | \(n\) | \(\#\{x\in\{\pm1\}^n:Cx=px\}\) | \(\#/n\) | \(\#/n^{3/2}\) |
|------:|------:|----------------------------------:|---------:|---------------:|
| 3 | 10 | 12 | 1.20 | 0.38 |
| 5 | 26 | 260 | 10.0 | 1.96 |
| 7 | 50 | 11452 | 229 | 32.4 |

JSON: `evidence/boolean_evec_counts.json`.

With \(x_0=+1\) folding, \(|\mathrm{Max}|=12,260,11452\) (both eigenspaces identified via \(x\leftrightarrow-x\)).

## Relevance to E(1)

Undercutting \(\Phi(C)\) requires \(S_F(y)\ge1\) for every boolean \(+\)-eigenvector \(y\) (and spike control on the rest of the cube). A naive “one private edge per maximizer” bound gives only
\[
k_\star\le|\mathrm{Max}|,
\]
but \(|\mathrm{Max}|/n^{3/2}\) is **increasing** through \(p=7\), so this upper bound is **not** \(o(n^{3/2})\) and does **not** by itself prove E(1).

A matching-sized (\(|F|\le n/2\)) simultaneous cover of Max would suffice for Prop 15.20d; that remains open for general \(p\).

## Not established

- \(k_\star=o(n^{3/2})\)
- closed form for the boolean-evec counts
- \(\lim\alpha_n\) exists
