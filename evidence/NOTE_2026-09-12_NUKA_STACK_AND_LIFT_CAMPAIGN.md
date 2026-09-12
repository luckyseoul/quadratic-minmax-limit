# Nuka ROCm stack audit, cupy reduction pathology, and the diagonal-lift campaign

**Status:** draft session note (2026-09-12). Measurements and finite results.
Independent review pending. No asymptotic claim.

## 1. ROCm stack audit (user prompt: "make sure you're using ROCm 10.0, not 7.2.x")

- **The cupy path on nuka is ROCm 7.2.4-locked.** The `amd_cupy-13.5.1` wheel
  (installed from `~/.cache/rocm-wheels/`) hard-links `libamdhip64.so.7`,
  `librocblas.so.5`, `libhipblaslt.so.1` from `/opt/rocm/lib`; forcing
  `LD_LIBRARY_PATH=/opt/rocm/core-10.0/lib` does **not** change the loaded
  runtime (`runtimeGetVersion` remains 71526333). Evidence: `ldd` records.
- **ROCm 10.0 is installed** (`/opt/rocm/core-10.0`, version 10.0.0) with full
  toolchain (`hipcc`, `amdclang++`) and math libraries (hipBLAS/hipBLASLt,
  MIOpen with gfx1201 kernels). AMD publishes **no cupy wheel for ROCm 10**
  (repo.radeon.com/rocm/manylinux tops out at `rocm-rel-7.2.4`; PyPI
  `amd-cupy` is a stub). Conclusion: **ROCm 10 is reachable via hipcc/hipBLAS
  C++ only**, unless cupy is built from source.
- hipBLAS under ROCm 10 on the wide-GEMM hot shape `(256x18)@(18x2^18)`:
  **0.766 ms/call = 1577 GMAC/s** (vs ~0.4 ms/call via cuBLAS on the V100).
  The GEMM was never the bottleneck.
- **Reduction pathology (the actual blocker):** full-array cupy reductions on
  the 7.2 stack are 100–450x slow on gfx1201: `cp.max` 75–280 ms, `cp.sum`
  120 ms, `cp.argmax` 284 ms over 67M f32 (V100 reference: ~1 ms).
  **Fix:** two-stage row reduction
  `av.reshape(-1, 4096).max(axis=1).max()` = **0.62 ms (~450x)**. After the
  fix the order-36 lift pass runs **42 s -> 5 s** on nuka (V100: 3 s).

## 2. Diagonal-lift D-optimization campaign (both GPUs)

Object `K(D) = [[A, A+D],[A+D, -A]]`, D diagonal ±1; exact Phi over
`2^(2n-1)` states via the associated wide-GEMM pass (`M2 = (A+D) @ Yt`
precomputed by associativity; one fast GEMM per chunk; anchors validated).

```
order  source             anchor(D=I)  best found  record  status
28     C14 (m14=21)       70           70          70      D=I stands (no gain)
30     K15 (m15=27)       83           81          75      record stands
32     bank n16 (32)      96           94          80      record stands (80 = different 16-source)
36     C18 (33)           108          108         108     D=I stands (~350 iterations)
38     K19 witness (39)   121          121         109     anchor = archived 121 exactly; C38=109 stands
```

The order-38 anchor (121) reproduces the archived "coherent order-38
diagonal-family minimum 121" precisely — a cross-validation of both. No run
improved its anchor; the diagonal-lift records stand as family optima so far.
No new finite records this campaign.

## 3. Ops summary

- Use the wide-GEMM formulation + two-stage reductions for this kernel family
  on both nodes.
- nuka: for new GPU kernels prefer hipcc/hipBLAS against
  `/opt/rocm/core-10.0`; the cupy path remains ROCm 7.2 (two-stage reduction
  workaround required); best single-core CPU on the mesh as before.

Evidence: `evidence/nuka_ops_and_lift_campaign_20260912.json`; tool
`scripts/lift_opt.py`.
