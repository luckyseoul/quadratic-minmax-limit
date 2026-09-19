# Paley \(n=p^2+1\) determines the cluster set of \(\alpha_n\)

**Status:** corollary of CORE §5 plus the archived \(\rho=1\) theorem.
Not a proof that \(\alpha_n\) converges. Not E(1). Original limit OPEN.

**Implication.** Existence of \(\lim\alpha_n\) is equivalent to existence
of \(\lim m_{p^2+1}/(p^2+1)^{3/2}\) along odd primes \(p\). The two-ray
criterion is not needed for this equivalence. The o-versus-\(\Theta\)
Paley gap along this single sequence is the whole question.

## 1. Ratio-dense saturating sequence

Odd primes satisfy \(p_{k+1}/p_k\to1\) by the prime-number theorem.
Hence \(n_k=p_k^2+1\) satisfies \(n_{k+1}/n_k\to1\). CORE §5 transfers
liminf and limsup from \((n_k)\) to every \(n\).

The Paley conference \(C\) of order \(n=p^2+1\) over \(\mathbf F_{p^2}\)
has a Boolean eigenvector \(Cx=px\) (halfspace construction, archived
`evidence/PROOF_rho_eq_1.md`, shipped
`paley_conference_prime_power` / `halfspace_boolean_vector`). Thus
\(\Phi(C)=\tfrac12 n\sqrt{n-1}\) exactly, not merely spectrally.

## 2. The gap sequence

\[
\gamma_p=\bigl(\Phi(C)-m_n\bigr)/n^{3/2}
=\tfrac12\sqrt{1-1/n}-\alpha_n.
\]

| p | n | Φ(C) | m_n | α_n | γ_p |
|---|---|------|-----|-----|-----|
| 3 | 10 | 15 | 13 exact | 0.41110 | 0.06325 |
| 5 | 26 | 65 | ≤ 61 record | ≤ 0.46012 | ≥ 0.03017 |

The n=26 value 61 is a verified upper bound, not a reviewed exact minimum.
Paired-field \(B\) only gives \(\gamma_p\le\tfrac12-B+o(1)\approx0.174\).

Relative undercut \(\Phi-m\) over \(\Phi\) is \(2/15=0.133\) at p=3 and at
least \(4/65=0.062\) at p=5. That is compatible with \(\gamma_p\to0\) and
does not exhibit a uniform \(c>0\).

## 3. What would close which side

- A construction, for infinitely many p, with
  \(\Phi(A)\le\Phi(C)-c n^{3/2}\) and \(c>0\) fixed, on this sequence,
  would give \(\limsup\gamma_p\ge c\), hence \(\liminf\alpha_n\le\tfrac12-c\).
  If the same \(c\) is also a limit of \(\gamma_p\), then
  \(\lim\alpha_n=\tfrac12-c\).
- An \(o(n^{3/2})\) bound on \(\Phi(C)-m_n\) along this sequence
  (in particular any \(O(n)\) Lipschitz undercut size) would give
  \(\gamma_p\to0\) and \(\lim\alpha_n=\tfrac12\).
- Oscillation of \(\gamma_p\) with two cluster points would prove
  nonexistence.

Matching flips at n=10 undercut by 2, which is \(O(n^{-1/2})\) after
normalization. They cannot produce a uniform \(c>0\). Gap-2 / residual
cover trees are not used here.

OpenAI referee `math_review` of the CORE §8 derivation: PASS
(do_not_branch). Claude was not used.

## 4. What this does not do

It does not prove \(\gamma_p\to0\), does not produce a definite \(c>0\),
and does not replace the paired-field lower bound. It names the only
sequence whose cluster set can differ from that of \(\alpha_n\).
