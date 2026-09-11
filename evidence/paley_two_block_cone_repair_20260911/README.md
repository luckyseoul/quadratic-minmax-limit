# Exact two-block cone-repair diagnostic

This is a finite test of the cross-block construction

\[
K(B)=\begin{pmatrix}C&B\\B^T&-C\end{pmatrix},
\]

where `C` is the order-14 (`q=13`) Paley conference signing.  It uses the
imported Navier correction pattern only to rank integral flips of the 196
entries of `B`.  A proposed block is scored by an exhaustive integer
enumeration of all `2^27` projective Boolean states of `K(B)`.

The initial one-witness experiment at deterministic seed 731 used 88 OpenMP
workers. At norm 94, none of its 16 ranked single flips or 120 ranked pairs
strictly lowered the complete Boolean norm. This exposed a method defect:
one tied maximizer is not an adequate stress set.

The corrected version retains up to 512 exact maximizers in the band
`|E|>Phi-2`. Soulkiller then ran seeds 731 through 752 concurrently, with
four scorer workers per seed (88 workers total). Every candidate still had a
complete `2^27`-state score. Nineteen of 22 starts strictly improved; the
best two-stage result, for seeds 741 and 742, is `Phi=90` from `Phi=92`.
NUKA independently compiled the scorer and replayed seed 741's retained
block with 16 workers, obtaining norm 90:

```json
{"n":28,"source_phi":21,"phi":90,"energy":90,"xmask":8160,"ymask":12451}
```

| item | SHA-256 |
| --- | --- |
| Soulkiller result | `ef85e6dcd636bd586f736b98ad6355c811a14b461f980bbc9887644c4d90383d` |
| exact scorer | `b1b6c95f582e86f807b181aa9b56599f8a56cb331f1f6be8edc28a86659d8b3c` |
| multi-active cone repair driver | `04bba6cc9b6c7b92372b87ffedf7e921fcecec21d53b9a0f037163722c6aa7bc` |
| HIP parity scorer | `15e81fdeabfd4bab1de1b2620b242845ff1cd1aa80cf792dacc90af346dcc763` |

The 22 raw batch files are in `multi_active_batch/`; seed 741 has SHA-256
`1af0205684c9b3f02815a8a4bf963c0f2e0e32328c7b1ff293716312f935d2b8`.

## NUKA ROCm capability check

NUKA's upgraded ROCm 10.0 installation exposes an AMD Radeon RX 9070 XT
(`gfx1201`) through `rocminfo`, and `hipcc` successfully compiled the
integer HIP scorer `scripts/paley_two_block_score.hip.cpp`.  On seed 741,
the HIP and independent 16-thread CPU scorers both returned exact norm 90.
For this order-28 single-score kernel, HIP measured 0.14 seconds and CPU
measured 0.11 seconds.  Consequently NUKA remains the CPU replay host for
this small exact workload; the HIP backend is retained for a larger kernel
where launch and reduction overhead no longer dominates. `rocm-smi` reports
an unavailable WSL driver, so dispatch capability is determined from
`rocminfo` and an actual HIP parity run, not that utility.

## Interpretation boundary

The best retained block provides only the finite statement `m_28<=90`.
The batch does not bound the optimum over cross blocks, rule out another
construction, or have an all-orders consequence. It does show that the
imported local correction mechanism has been applied to the actual two-block
lift with exact, rather than sampled, acceptance.
