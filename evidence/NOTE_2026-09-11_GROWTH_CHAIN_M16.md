# Growth-chain order-16 bound: m_16 <= 30 (independently verified)

**Status:** draft session note (2026-09-11). Sourced from the archived growth
checkpoint (`evidence/original_mo_growth_checkpoint_20260906/checkpoint.tar.gz`,
dir `original-mo-two-vertex-cpu.0eqNgg/`); independently re-verified this
session on two engines (NumPy and CuPy). Matches the checkpoint README line
"their best endpoint-optimal intermediate norms are 34 and 30". No all-orders
claim.

## Object

`K16 = [[S15, b], [b^T, 0]]`, where `S15` is the order-15 source with
`Phi(S15) = 27` and `b` is the tie-optimized two-vertex extension vector
(`nuka15.tie.json`), with the order-17 endpoint (through vector `c`, `d = 1`)
at its optimum `34`.

## Verified facts (this session)

- `K16` is a complete signing (`zero diag`, `+-1` off-diag, symmetric);
  `Phi(K16) = 30` by exact enumeration over `2^15` projective states,
  reproduced identically on CPU NumPy and on CuPy (GPU) — equals the
  recorded value.
- Reconstruction check: `[[S15, b], [b^T, 0]]` equals the archived matrix
  `nuka15.tie.first.txt` entrywise (`True`).
- Endpoint check: `[[S15, b, c], [b^T, 0, 1], [c^T, 1, 0]]` equals
  `nuka15.tie.second.txt` entrywise (`True`); `Phi = 34`.
- Regenerated chain (`nuka16.json` + `.first/.second.txt`): order-16 source
  `Phi = 30`; first addition `Phi = 36`; second `Phi = 37`; both rebuilt
  entrywise from the recorded `(b, c, d)`.

## Consequence

`m_16 <= 30` (the broad-campaign winner table recorded `32`). With
`m_16 >= m_15 = 27` by monotonicity: `m_16 in [27, 30]`.

Note for the multiplier-2 program: `2 sqrt(2) * m_8 = 28.28`, so the
order-8 doubling inequality `m_16 <= 2 sqrt(2) m_8` is not yet established
(`30 > 28.28`); the bracket `[27, 30]` straddles it.

## Repro

- Verifier: `~/scratch/verify_growth_chain.py` (session scratch); output
  `~/scratch/growth_chain_verify.json`.
- Witness landed: `evidence/m16_phi30_witness_20260911.json` (`A`, `Phi`,
  provenance).
- Archive: `checkpoint.tar.gz` (hash in `archive_receipt.json` in the same
  evidence directory).
