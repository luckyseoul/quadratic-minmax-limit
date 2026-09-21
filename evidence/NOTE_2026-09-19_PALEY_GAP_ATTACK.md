# Paley-gap attack, 2026-09-19

**Status:** original limit still OPEN. Neither
\(m_{p^2+1}\ge\Phi(C)-o(n^{3/2})\) nor a ratio-dense definite-\(c\)
family is proved. Scope corrected 2026-09-21: the spectral calculation
below limits one upper certificate; it does not exclude the construction.

## 1. The coarse spectral certificate does not certify a definite gap

For odd-order Seidel \(B\) and \(K(B)=\bigl(\begin{smallmatrix}B&B+I\\B+I&-B\end{smallmatrix}\bigr)\),

\[
\|K(B)\|_{\mathrm{op}}\le\sqrt2\,\|B\|_{\mathrm{op}}+1.
\]

Always \(\|B\|_{\mathrm{op}}\ge\sqrt{m-1}\). If \(B\) is near-min-op
(\(\|B\|_{\mathrm{op}}=(1+o(1))\sqrt{m}\)), the bound yields
\(\Phi(K)\le n\|K\|_{\mathrm{op}}/2\) hence
\(\Phi(K)/n^{3/2}\le\tfrac12+o(1)\). If \(\|B\|_{\mathrm{op}}\) is
larger, this upper certificate is weaker; the actual Boolean norm need
not be larger. An upper bound tending to 1/2 does not exclude an actual
norm tending to a smaller constant.
The Paley-\(\mathbf F_{25}\) lift satisfies this estimate:
\(\|K\|_{\mathrm{op}}=\sqrt{61}\le\sqrt2\cdot5+1\).

The precise limitation can be stated for every block order \(m\ge2\).
The normalized certificate obtained by substituting this operator bound is
\[
U(B)=\frac{\sqrt2\|B\|_{\rm op}+1}{2\sqrt{2m}}
\ge\frac12\sqrt{1-1/m}+\frac1{2\sqrt{2m}}\ge\frac12.
\]
For the last inequality, use \(\sqrt{1-1/m}\ge1-1/m\) and
\(1/\sqrt{2m}\ge1/m\). Thus this particular certificate cannot prove
a sub-half constant. This is not a lower bound on \(\Phi(K)/n^{3/2}\).

The three certified construction undercuts at \(n=10,26,50\) are
\(2,4,6\). They fit \(p-1\), but no all-orders formula is proved by
these three values. Even an all-orders \(O(n)\) construction undercut
would only give a lower bound on the optimal gap \(\gamma_p\), not an
upper bound. No definite-gap exclusion for the whole family follows.

Tests: `tests/test_plus_i_opnorm_bound.py`.

## 2. Best-response local fields of the known beaters

For a state attaining the absolute norm, choose \(s\in\{-1,1\}\) so
that \(Q_{sA}(x)=\Phi(A)\). Single-spin stability gives
\(w=(sAx)\circ x\ge0\) and \(\Phi(A)=\tfrac12\|Ax\|_1\).
At the even orders in the table, the fields are odd, hence strictly
positive. A merely local maximizer need not attain \(\Phi(A)\).

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
   (would be a definite-\(c\) construction if coherent;
   the coarse certificate in §1 cannot certify it; §1 does not rule out
   an actual Boolean improvement).

Avoid list matched the dead routes: op-norm only, entrywise Max-Lipschitz
with \(k_\star=\Theta(n^2)\), unamplified plus-I searches.

## What remains

One sufficient target is the lower bound
\(m_{p^2+1}\ge\Phi(C)-o(n^{3/2})\) for arbitrary Seidel \(A\); that
would prove the particular value 1/2. Convergence to a smaller value, or
another direct cross-order argument, remains possible. No requirement
to avoid every plus-I construction follows from this note. This scope
correction does not authorize restarting the user-closed family scans.
One proposed object on the 1/2 route is a Boolean (or algebraic) vector
in the top space of an arbitrary almost-conference \(A\), not a
Gaussian in that space and not Max+ of Paley when \(k_\star=\Theta(n^2)\).
