# Experiment 3 — exact doubling vs neutral composition

## Question

Same open limit for \(\alpha_n=F(n)/n^{3/2}\).

## Motivation

Experiments 1–2 showed a vanishing one-sided modulus is plausible on the
exact prefix under \(\varepsilon_n=n^{-1/2}\), but recursive \(n-2\) padding
and the conference ceiling cannot produce it. The missing object is a
cross-order construction with normalized error \(\to 0\). The simplest
cross-order probe is **doubling**: compare \(F(2n)\) to the constant-\(\alpha\)
neutral value \(2^{3/2} F(n)\).

## Working hypothesis (H3)

1. Among exact pairs \((n,2n)\) with \(2n\le 15\), at least one is
   **sub-neutral**:
   \[
   F(2n) < 2^{3/2}\,F(n)
   \qquad\bigl(\Leftrightarrow\;\alpha_{2n}<\alpha_n\bigr).
   \]
2. Restricting to \(n\ge 5\), the worst exact doubling excess
   \(\eta_n=\alpha_{2n}-\alpha_n\) is \(\le 0.10\).
3. The open case \(F(16)\in\{28,30\}\) sits on the neutral knife-edge for
   \(n=8\):
   - \(2^{3/2} F(8)=10\cdot 2^{3/2}\approx 28.284\),
   - so \(F(16)=28\) is sub-neutral (\(\eta_8<0\)),
   - while \(F(16)=30\) is super-neutral (\(\eta_8>0\)).

## Method

Exact A399208 values; arithmetic only. Also scan all pairs
\((n,m)\) with \(8\le n<m\le 15\) and \(m/n\in[1.8,2.2]\) as near-doublings.

## Expected result

H3 holds. In particular, an exact sub-neutral doubling already exists
(\(7\to 14\): \(F(14)=21 < 9\cdot 2^{3/2}\approx 25.46\)), and deciding
\(F(16)\) is deciding whether the first unavailable doubling stays
sub-neutral.

## Analysis / implication

Sub-neutral doublings are finite witnesses that cross-order improvement
below constant-\(\alpha\) is possible. They do not prove an all-future
modulus, but they identify \(F(16)\) as the immediate computational
bottleneck for the doubling form of the one-sided criterion.

## Limitations

- Only \(n\le 15\) exact; \(F(16)\) still open.
- Doubling is one scale ratio; a full all-future bound needs all \(m\ge n\).
- Sub-neutrality at sparse \(n\) can coexist with nonconvergence.
