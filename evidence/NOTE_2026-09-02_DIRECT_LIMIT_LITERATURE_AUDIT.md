# Direct-limit literature audit after the geometric calculator

**Status:** current primary-source audit; no cited theorem proves convergence
of `m_n/n^(3/2)`.  Two exact encodings below identify hypotheses that would
prove it and explain why the nearest 2026 results do not supply them.

## 1. A sufficient Bernoulli lower-tail theorem

Put `E_n=binom(n,2)`.  For a uniformly random edge signing `J`, define

\[
 G_n(J)=n^{-3/2}\max_x Q_J(x),\qquad
 W_n(J)=\max\{G_n(J),G_n(-J)\}.
\]

The deterministic outer minimum is exactly

\[
 {m_n\over n^{3/2}}=\min_J W_n(J).                 \tag{1}
\]

Suppose one could prove, at speed `E_n`,

\[
 -{1\over E_n}\log\Pr\{W_n\le t\}\longrightarrow I(t)               \tag{2}
\]

with a unique threshold `t_*` across which `I(t)` crosses `log 2`.  Then
`m_n/n^(3/2) -> t_*`.  Indeed, every nonempty event in the uniform signing
space has mass at least `2^(-E_n)`.  Thus `I(t)>log 2` forces eventual
emptiness, whereas `I(t)<log 2` forces a witness.

The required event is joint and lower-tailed:

\[
 G_n(J)\le t\quad\hbox{and}\quad G_n(-J)\le t.       \tag{3}
\]

Chen--Guionnet--Ko--Lacroix-A-Chez-Toine--Mourrat,
[arXiv:2603.06368](https://arxiv.org/abs/2603.06368), instead prove an upper
large-deviation principle for the maximal energy with Gaussian disorder.
They explicitly distinguish the zero-field lower tail, expected at speed
`n^2`, and note that the rigorous lower-tail results available there are for
spherical models.  Their theorem is neither Bernoulli nor the joint event
(3).  Guerra--Toninelli,
[cond-mat/0204280](https://arxiv.org/abs/cond-mat/0204280), proves the usual
thermodynamic limit for typical disorder; averaging over disorder does not
survive the outer minimum in (1).

## 2. Exact matrix-discrepancy encoding, wrong dimension regime

Index rows by projective Boolean states `[x]` and define diagonal matrices

\[
 D_{ij}=\operatorname{diag}(x_i x_j)_{[x]}.
\]

Then, exactly,

\[
 m_n=\min_{\varepsilon\in\{\pm1\}^{E_n}}
 \left\|\sum_{i<j}\varepsilon_{ij}D_{ij}\right\|_{op}.                \tag{4}
\]

Akbas--Sra's new Matrix Spencer theorem,
[arXiv:2608.28816](https://arxiv.org/abs/2608.28816), gives `O(sqrt N)` for
`N` symmetric `N`-by-`N` contractions.  Encoding (4) has `E_n` matrices but
dimension `2^(n-1)`.  Their hereditary theorem in the rectangular regime
has radius

\[
 \kappa\sqrt{E_n}\left(1+\log{2^n\over E_n}\right)^2,
\]

which is much larger than the needed `Theta(n^(3/2))` scale.  Their earlier
algebraic theorem, [arXiv:2606.16005](https://arxiv.org/abs/2606.16005),
requires a generated `C^*`-algebra of dimension `O(E_n)`, while the diagonal
Walsh algebra in (4) has dimension `2^(n-1)`.

## 3. Audit verdict

Ordinary dense graph limits see the `n^2` scale and collapse every optimizer
to the half-density graphon; they do not see the second-order `n^(3/2)`
quantity.  Standard matrix discrepancy sees the correct signing problem but
not its exponential family of Boolean tests.  Existing spin-glass
interpolation sees typical disorder rather than the rare optimal signing.
The dimension-free bounded-degree influence estimates of
Filmus--Hatami--Keller--Lifshitz,
[arXiv:1404.3396](https://arxiv.org/abs/1404.3396), control the sharp constant
in Proposition 5.3 of `solution.md`, but supply no relation between dimensions
`n` and `n+1` and hence no convergence theorem.

The literature therefore leaves two theorem-shaped direct targets:

1. the joint Bernoulli lower-tail threshold (2)--(3); or
2. a dimension-sensitive discrepancy/influence theorem specialized to the
   degree-two Walsh test family.

Neither target is presently proved, so this note is an audited reduction,
not a closure claim.
