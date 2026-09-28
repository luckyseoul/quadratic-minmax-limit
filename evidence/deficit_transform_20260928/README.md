# Deficit-transform verification and next discriminator

## Completed verification

`regression.json` records one external CPU execution of the 13 small exact
tests in `tests/test_deficit_transform_certificate.py`. All passed, exit 0,
at source commit `fd05d10b334e45c1d3b291cb8d89c7d3aa257eea`, with one worker.
The test runner reported 0.067 seconds. The staged input hashes matched the
reviewed files before execution.

The receipt SHA-256 is
`ef767d48405195f2c42480f4cfe91c3870f5851ce9be5369f554cc6315b12ab9`.
This verifies implementation fixtures. It is not a research-source run or
an all-orders certificate. No major milestone or new backup is claimed.

## Prepared research calculation: not executed

Question: on the recorded order-15 signing of norm 27, does the **full**
deficit-aware section construction recover an incident row with extended
norm at most 30? This includes positive deficits, both phases, retained
section slabs, and reverse sign recovery. The earlier zero-deficit
same-phase depth computation did not answer that question and will not be
rerun.

Inputs:

- `scripts/deficit_transform_certificate.py`, SHA-256
  `5e864e105f59e33e6cb1e803371de0763237de8b2cacebf1c0a35bf0de8735a4`.
- `evidence/K15_exact_minimizer_20260911.json`, SHA-256
  `4de6bae5ce2c3fd144748407b98c3b26bd1bfa397f557130fa7568aa59ad9f37`.

From a checkout containing these inputs, the single external command is:

```sh
timeout 120s python3 -B scripts/deficit_transform_certificate.py --matrix evidence/K15_exact_minimizer_20260911.json --enumerate-stable --target-phi 30 --max-states 16384 --max-pairs 200000 --max-facets 20000 --save-skeleton evidence/deficit_transform_20260928/k15-stable.json --output evidence/deficit_transform_20260928/k15-result.json
```

The program is serial. It enumerates the 16,384 projective states once,
saves the stable skeleton, and never overwrites an existing output. Any
later target or ordering change must reuse the skeleton. A resource stop
or timeout is inconclusive and does not authorize an automatic larger run.
Even successful recovery would be a finite method discriminator, not a
new value of `m_16` or a convergence proof.

Preparation checks found no matching worker on the four compute hosts,
and no existing full-deficit result in local evidence or matching top-level
remote `/tmp` staging names. This is a scoped search, not an audit of every
remote directory. The inputs were copied to a fresh external staging
directory and their hashes were verified there. The research
command has **not** been launched; the handoff's user-run instruction is
pending explicit direction. The convergence target remains an all-orders
extension estimate with summable excess.
