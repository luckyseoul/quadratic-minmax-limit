# Paley-gap attack, 2026-09-19

**Status:** original limit still OPEN. Neither
\(m_{p^2+1}\ge\Phi(C)-o(n^{3/2})\) nor a ratio-dense definite-\(c\)
family is proved. This note records what the attack *did* close.

## 1. Plus-I is not a definite-\(c\) construction

For odd-order Seidel \(B\) and \(K(B)=\bigl(\begin{smallmatrix}B&B+I\\B+I&-B\end{smallmatrix}\bigr)\),

\[
\|K(B)\|_{\mathrm{op}}\le\sqrt2\,\|B\|_{\mathrm{op}}+1.
\]

Always \(\|B\|_{\mathrm{op}}\ge\sqrt{m-1}\). If \(B\) is near-min-op
(\(\|B\|_{\mathrm{op}}=(1+o(1))\sqrt{m}\)), the bound yields
\(\Phi(K)\le n\|K\|_{\mathrm{op}}/2\) hence
\(\Phi(K)/n^{3/2}\le\tfrac12+o(1)\). If \(\|B\|_{\mathrm{op}}\) is
larger, the same upper bound is worse, so \(K\) is a weaker beater.
The Paley-\(\mathbf F_{25}\) lift saturates the picture:
\(\|K\|_{\mathrm{op}}=\sqrt{61}\le\sqrt2\cdot5+1\).

Certified undercuts of Paley conference remain \(p-1\) at \(n=10,26\)
(and \(\ge6\) locally at \(n=50\)): \(o(n^{3/2})\) after
normalization. This family cannot produce a uniform \(c>0\).

Tests: `tests/test_plus_i_opnorm_bound.py`.

## 2. Best-response local fields of the known beaters

A local maximizer satisfies \(x=\mathrm{sign}(Ax)\) and
\(\Phi=\tfrac12\|Ax\|_1\). Write \(w=(Ax)\circ x\ge0\). Then \(w_i\)
are odd positive integers of the same parity as \(n-1\).

| matrix | \(\Phi\) | \(w\)-multiset |
|---|---|---|
| Paley \(n=10\) | 15 | \(3^{10}\) |
| plus-I \(C_5\) | 13 | \(1^6,5^4\) |
| Paley \(n=14\) | 21 | \(1^7,5^7\) |
| Paley \(n=26\) | 65 | \(5^{26}\) |
| plus-I \(n=26\) record | 61 | \(1^6,5^{16},9^4\) |

Plus-I \(A^2\) off-diagonals stay in \(\{0,\pm4\}\) at \(n=10,26,50\),
with \(2,6,12\) nonzeros per row. Entrywise-bounded \(A^2-(n-1)I\)
is the only observed beater shell; it is not a classification.

Tests: `tests/test_beater_magnitude.py`.

## 3. Schur on matrix Aut does not give the Max+ frame

The permutation automorphism group of a standard Paley conference
representative (AΓL translations, *square* multiplications, Frobenius,
∞ fixed) commutes with \(C\) but has \(\mathrm{Hom}_{\mathrm{Aut}}(V_+,V_+)\)
of dimension **2** at \(p=3,5,7,11\). So Aut is reducible on \(V_+\)
and the frame identity \(\mathbb E[yy^T]=I+C/p\) does not follow from
Schur on this group. The identity remains a finite certificate at
\(p=3,5,7\).

Tests: `tests/test_paley_Vplus_commutant.py`. Receipt:
`evidence/paley_Vplus_irrep_20260919/commutant.json`.

## 4. Gaussian \(u_{ij}\)-correction is too small

CORE §4 already saturates at \(1/\pi\) on conference (\(u_{ij}=0\)).
The next Taylor term in \(\arcsin(u\pm v)\) improves the floor by
\(\Theta(\mathrm{mean}\,u^2)=\Theta(\|D\|_F^2/n^4)\). On the \(n=26\)
beater this is \(<10^{-3}\). It cannot raise \(1/\pi\) to \(1/2\).

## 5. Referee directions (OpenAI `suggest_direction`)

After one 180s timeout, a second call returned three directions, not
theorems. Branch only on BLOCK; this was direction-only.

1. Fourth-order stability dichotomy on \(D=A^2-(n-1)I\), with a
   *non-Gaussian* rounding (Gaussian is the avoid list).
2. Zero-temperature interpolation / Aut-averaged Gibbs (pays \(k_\star\)
   unless the Gibbs leaves Max+).
3. Amplify the \(p-1\) defect on affine lines of \(\mathbf F_{p^2}\)
   (would be a definite-\(c\) construction if coherent; the natural
   amplification of plus-I is ruled out in §1 at the spectral level).

Avoid list matched the dead routes: op-norm only, entrywise Max-Lipschitz
with \(k_\star=\Theta(n^2)\), unamplified plus-I searches.

## What remains

The missing inequality is still the lower bound
\(m_{p^2+1}\ge\Phi(C)-o(n^{3/2})\) for arbitrary Seidel \(A\), or a
construction that is *not* a plus-I lift of a near-min-op block.
The live non-tautological object is a Boolean (or algebraic) vector
in the top space of an arbitrary almost-conference \(A\), not a
Gaussian in that space and not Max+ of Paley when \(k_\star=\Theta(n^2)\).
