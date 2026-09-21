# Central-shell exchange, then a corrected leftover

2026-09-20. Ran the gpt-6-astra direction (confirm central-shell /
superlevel exchange; no Φ-drop from a global minimizer; no raw nauty
census). Original limit OPEN. This note is a diagnostic of that
direction, not a 15.xxx and not a doubling proof.

## 1. The Φ-drop branch is empty

Let \(A\) be a global minimizer, \(\Phi(A)=m_n\). For every
\(A'\in\mathcal S_n\), \(\Phi(A')\ge m_n\). Coordinated Seidel flips
cannot decrease \(\Phi\). The only remaining exchange on the fibre
\(\{A:\Phi(A)=m_n\}\) is \(\Phi\)-preserving descent of a secondary
functional \(W\).

## 2. Finite numerical evidence against the zero-error bound

Write \(W(A)\) for the energy-shell maximum in
`NOTE_2026-09-05_GAUSSIAN_ENERGY_SHELL_UPPER` display (3), at the
spectral-midpoint covariance. On four recorded exact minimizers
(Monte Carlo widths, \(x_0=+1\) representatives):

| matrix | \(n\) | \(\Phi\) | \(W\) | \(2\sqrt2\Phi\) | \(W/n^{3/2}\) |
|---|---|---|---|---|---|
| C5 | 5 | 4 | 16.83 | 11.31 | 1.50 |
| Paley conference | 6 | 5 | 21.89 | 14.14 | 1.49 |
| plus-I of C5 | 10 | 13 | 54.55 | 36.77 | 1.72 |
| \(B_{13}\) | 13 | 20 | 80.71 | 56.57 | 1.72 |

The maximizing pair \((h,k)\) is opposite near-extreme shells, not
the equatorial cell. Receipt:
`evidence/shell_W_minimizer_probe_20260920.json`,
`scripts/shell_W_minimizer_probe.py`.

These fixtures challenge the zero-error bound at the listed orders.
They do not contradict
\(W(A_n)\le2\sqrt2\Phi(A_n)+o(n^{3/2})\), nor a Dini-summable
remainder. Explicitly, for any finite set of exceptions with largest
order \(N_0\), choose C at least their largest positive normalized
excess and set \(r(n)=C\) for \(n\le N_0\), \(r(n)=0\) otherwise.
This is nonnegative and nonincreasing; its dyadic tail is zero once
the starting order exceeds \(N_0\). Such an envelope does not prove
the untested orders, but proves that finite exceptions alone cannot
exclude that asymptotic target. The previous Dini exclusion was invalid.
Moreover, the receipt uses Monte Carlo for the Gaussian quantities;
its precision must not be presented as an exact all-orders certificate.

## 3. Display (3) is not the leak at these orders

The linearized midpoint model of the expected-paired-norm package
has covariance \(C=I+\kappa(A\otimes A)/\mu\) when \(a=b\). Monte
Carlo of \(\mathbb E\mathcal M_A(Z)\) on the same four matrices
agrees with \(W\) to a few percent, and both have positive measured
normalized excess at these four orders. No asymptotic \(\Theta\) claim
follows from that finite table:

| matrix | \(\mathbb E\mathcal M\) | \(\mathbb E\mathcal M-2\sqrt2\Phi\) |
|---|---|---|
| C5 | 17.48 | 6.16 |
| Paley6 | 23.03 | 8.89 |
| plus-I C5 | 55.00 | 18.23 |
| \(B_{13}\) | 80.05 | 23.48 |

Tightening the two-field split (10) does not recover the dyadic
constant at these \(n\): the actual linearized objective is already
too large. The Boolean pairing \(B=\mathrm{sign}(G)\) at threshold
\(h=0\) is the same size (C5: \(\mathbb E\Phi\approx 16.49\)).

The ratio \(n^{16/11}/n^{3/2}=n^{-1/22}\) is not small at these orders
(and the implicit constant is separate), so this table is **not** an asymptotic
counterexample to the reduction \(m_{2n}\le\mathbb E\Phi+\cdots\).
It is numerical evidence against a zero-error \(2\sqrt2\Phi\)
comparison at the tested orders, not an exact or asymptotic certificate.

## 4. Shifted thresholds collapse to the coherent pairing

For \(B_h=\mathrm{sign}(G+hA)\) (diagonal zeroed), at \(n=5\) and
\(n=6\) the Monte Carlo \(\mathbb E\Phi([[A,B_h],[B_h^\top,-A]])\)
decreases on the tested \(|h|\) grid toward the deterministic coherent
pairing \(\Psi(A):=\Phi([[A,A],[A,-A]])\):

- C5: the observed best value is \(\Psi=12\), still \(12>11.31=2\sqrt2\Phi\).
- Paley6: \(\Psi=18>14.14\).

Always \(\Psi(A)\ge 2\Phi(A)\), by the diagonal sections \(y=\pm x\).
Plus-I of C5 is 13, one more than \(\Psi\).

These samples do not certify an infimum over every real threshold h.
The coherent inequality \(\Psi\le 2\sqrt2\Phi\) fails at
the tested orders \(n=5,6\). The excess at
\(n=5\) is \(0.69\); if \(\Psi=2\Phi+O(n)\) then
\(\Psi/n^{3/2}=2\alpha+O(n^{-1/2})\) and
\(2\alpha<2\sqrt2\alpha\) with a uniform margin since the proved
lower bound keeps \(\alpha\) bounded away from zero. That hypothesized
identity would imply the comparison eventually, but is not proved here.

## 5. Corrected leftover (one implication)

For an exact minimizer \(A\), prove
\[
\Psi(A)=\Phi\begin{pmatrix} A & A \\ A & -A \end{pmatrix}
\le 2\sqrt2\,\Phi(A)+n^{3/2} r(n)
\]
with \(r\ge 0\) nonincreasing and \(\sum_j r(2^j)<\infty\), **or**
exhibit a threshold \(h=h_n\) at which the Boolean pairing of the
shifted-threshold theorem beats \(\Psi\) by a definite
\(n^{3/2}\) term.

Do not infer an asymptotic exclusion from the \(n\le13\) tables or
claim a Dini bound from them. Neither the general shell comparison nor
the displayed coherent-pairing comparison is settled here. The user-closed
Hadamard/family scans remain closed; this correction is not a restart.
Limit OPEN.
