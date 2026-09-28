# Deficit-transform verification and first full-deficit execution

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

## Completed bounded research calculation

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

The single external command ran from a fresh staging directory containing
these inputs with their relative paths preserved:

```sh
timeout 120s python3 -B scripts/deficit_transform_certificate.py --matrix evidence/K15_exact_minimizer_20260911.json --enumerate-stable --target-phi 30 --max-states 16384 --max-pairs 200000 --max-facets 20000 --save-skeleton k15-stable.json --output k15-result.json
```

It finished in 0.906 seconds, exit 2, with `RESOURCE_LIMIT`: more than
20,000 distinct facets appeared in the first section. The trace is empty
because no full transform step completed. No row was recovered. This is
an inconclusive resource stop, not a failed extension inequality. Neither
the facet limit nor the time limit was increased.

The run did enumerate all 16,384 projective states and save the complete
340-class stable skeleton. `execution.json` records the command, budgets,
input hashes, output hashes, empty stdout/stderr, and exit status. The
retained outputs have SHA-256 hashes:

- `execution.json`:
  `d40795271d5483be489307c9a6540f8127bcc59bd1bf992c372444988585c19b`.
- `k15-stable.json`:
  `fd997da804aaf70855cbca14e31610f872216e83e49cb32f58aa7a8da55e0ff6`.
- `k15-result.json`:
  `647763aa5579fb3d6f4aed21967c169b74c46541f418d56353e553564f0c26d1`.

Reuse the skeleton for any justified later calculation. Cached-mode output
retains its conditional-completeness label; this execution receipt and its
matching cache hash supply the earlier complete-enumeration provenance.
Do not repeat the source scan to obtain another completeness receipt.

Preparation checks found no matching worker on the four compute hosts,
and no existing full-deficit result in local evidence or matching top-level
remote `/tmp` staging names. This is a scoped search, not an audit of every
remote directory. The inputs were copied to a fresh external staging
directory and their hashes were verified there. The one-worker run ended
normally with its explicit resource status; no background worker remains.

## Cache-only core profile

`core-profile.json` was computed from the saved skeleton, without another
state scan. It counts classes by `(27-energy, phases)`. For each cutoff
`d=0,2,4`, its core is exactly the listed stable states of deficit at most
`d`; coordinates are grouped by the core columns, with the first signature
entry made positive. The displayed signature classes and odd-class count
are direct data summaries, not a search over smaller cores.

The all-zero-deficit core has 66 states. Its signature classes are the pair
`{0,14}` and 13 singleton coordinates, leaving `u=13`; this specific choice
does not meet the balanced-core requirement `u<7/2`. Including the deficit-2
or deficit-4 states leaves that partition unchanged. This neither excludes
other cores nor proves failure of the local-lemma construction.

The all-orders extension estimate with summable excess remains unproved.
No larger run, new signing census, alternate-engine replay, or new bound on
`m_16` is claimed or authorized by this record.
