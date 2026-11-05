# Seidel quadratic limit

## Major result

**Theorem:** given \(F(8)=10\) and \(F(16)\in\{28,30\}\), the CORE §7 doubling
inequality \(H(16)\le 2H(8)\) holds if and only if \(F(16)=28\).
Proof: \(28^2 < 8\cdot10^2 < 30^2\) (i.e. \(784<800<900\)).

See [RESULT.md](RESULT.md).

## Experiments

Finite diagnostics for whether Ivanisvili’s sequence

\[
\alpha_n = n^{-3/2}
\min_{a_{ij}=\pm1}
\max_{x\in\{\pm1\}^n}
\Bigl|\sum_{i<j}a_{ij}x_ix_j\Bigr|
\]

can obey a one-sided all-future modulus (sufficient for convergence; see
[HANDOFF 2026-11-04](https://github.com/luckyseoul/quadratic-minmax-limit)).

The limit itself remains open.

## Reproduce

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File experiment/Run-Experiment1.ps1 -MaxExact 6
powershell -NoProfile -ExecutionPolicy Bypass -File experiment/Run-Experiment2.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File experiment/Run-Experiment3.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File experiment/Run-Experiment4.ps1
```

Outputs: `results/experiment{1,2,3}.*`.

Method cards: [METHOD.md](METHOD.md), [METHOD2.md](METHOD2.md), [METHOD3.md](METHOD3.md](METHOD3.md), [METHOD4.md](METHOD4.md).

## Upstream research

- https://github.com/luckyseoul/quadratic-minmax-limit
- https://github.com/Robby955/mo-413935-research
- https://oeis.org/A399208
- https://mathoverflow.net/questions/413935

