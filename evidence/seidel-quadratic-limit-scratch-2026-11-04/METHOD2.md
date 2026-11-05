# Experiment 2 — do existing majorants give a vanishing one-sided modulus?

## Question

Same open limit: does \(\alpha_n=F(n)/n^{3/2}\) converge?

## Motivation

Experiment 1 showed the finite exact excess \(E_n(15)\) fits under
\(\varepsilon_n=n^{-1/2}\), but is not eventually monotone. HANDOFF
2026-11-04 says adjacent-step control is not enough unless it upgrades to
a uniform all-future bound. CORE.md already has:

- \(F(n+1)\le F(n)+n-2\) for \(n\ge5\),
- \(\alpha_n\le 1/2\) in the limsup sense (and every finite upper construction
  used here is at most \(1/2\) after normalization only asymptotically;
  we use the analytic ceiling \(1/2\) as the best uniform majorant on \(\alpha\)).

## Working hypothesis (H2)

Neither available unconditional majorant supplies a vanishing one-sided
all-future modulus along the exact table:

1. **Padding / iterated step.** Extending with the sharp step bound
   \[
   F(n+d)\;\le\;F(n)+\sum_{k=0}^{d-1}(n+k-2)
   = F(n)+d(n-2)+\tfrac{d(d-1)}{2}
   \]
   produces
   \[
   P_n(d)=\frac{F(n)+d(n-2)+d(d-1)/2}{(n+d)^{3/2}}-\alpha_n
   \]
   with \(\sup_{d\le 500} P_n(d)\ge 0.1\) for every \(n\in\{5,\dots,15\}\).
   (In fact \(P_n(d)\to\infty\) as \(d\to\infty\).)

2. **Conference ceiling.** The residual \(\rho_n=1/2-\alpha_n\) on
   \(n=8,\dots,15\) stays \(\ge 0.03\) and is not monotone decreasing to 0.
   So \(\varepsilon_n=\rho_n\) does not visibly tend to 0 on the exact prefix.

3. **Neutral comparison.** The constant-\(\alpha\) neutral increment
   \(r_n=F(n)\bigl((1+1/n)^{3/2}-1\bigr)\) is beaten from above
   (\(\delta_n=F(n+1)-F(n)>r_n\)) for at least three orders in \(5..14\).
   Super-neutral steps keep reopening local excess.

## Method

Exact \(F(1..15)\) from A399208; no local search. Pure arithmetic over the
table, plus a scan of \(P_n(d)\) for \(d=1..500\).

## Expected result

All three H2 clauses hold: padding excess large, conference residual
stubborn, several super-neutral steps. This does not prove nonconvergence;
it delimits why the next ingredient must improve cross-order constructions
below the conference gap or force a rate \(\alpha_n\to L\).

## Limitations

- Uses the analytic ceiling \(1/2\), not a proved \(\alpha_n\le 1/2\) for each
  finite \(n\) (limsup bound). Finite constructive ceilings may be worse.
- Exact prefix only through 15.
- Super-neutral counts are descriptive, not a Dini-series proof.
