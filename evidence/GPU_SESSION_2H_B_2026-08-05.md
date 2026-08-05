# GPU session B (V100, ~2h window) — 2026-08-05

**Backend:** CuPy / Tesla V100-SXM2-16GB, free ~15 GiB. Real device wall on Max+ Gram/m4 (not theater).
**L status:** OPEN. `residual_closed_general=false`. No soft-close.

## Shipped: Prop 15.165

| Item | Result |
|------|--------|
| Exact global Es4 p=3 | 1680 = Es4_* (Φ flat at 16; 16N equality) |
| Exact global Es4 p=5 | 120400/13 ≤ Es4_*=492752/53 (ratio 0.99616) |
| Exact global Es4 p=7 | 5218435600/167281 ≤ Es4_*=8116400/251 (ratio 0.96473) |
| p=7 W_CENSUS trap | single-root W gives 12835984/409 **≠** global (not 1-homogeneous) |
| Closed Es4_* | 4(3p⁸−6p⁶−148p⁴−10p²+129)/(p⁴−8p²−1) |
| Closed η_* | (p²−1)n(19p⁴−152p²−3)/(6p²(p⁴−8p²−1)); η_*/n→19/6 |
| GoG ↔ Φ | spec(GoG)={2nN}∪{N·spec(Φ)}∪{0}^{N−1−m} |
| m4 C-eigen | p·m4(a,b,c,d)=∑_j C_aj m4(j,b,c,d) (cert p=5,7) |
| H-sat p=5 | λ_max=176/13=6+2H |
| E[s]=0, E[s²]=2n | 2-design + central symmetry |

## GPU jobs (device wall)

| Job | Wall |
|-----|------|
| Angle multiset / Es4 Gram p=5,7 | ~3s (batched matmul) |
| Full m4 on all 4-sets p=5 (14950), p=7 (230300) | ~1.5s |
| m4 by C-signature stratification | ~1.5s |
| Φ spectrum rebuild p=3,5 | <1s |

## Dead this session (for general Es4≤Es4_*)

- Design LP / crude |s|≤n−4 angle LP (far above Es4_*)
- C-edge signature alone does **not** determine m4 (multiple m4 per sig)
- Single-root W_k at p≥7 (not 1-homogeneous)

## Still OPEN

`E[s⁴] ≤ Es4_*(p)` for all primes p≥5 (or ∑η²≤η_*). Then 16N (15.164) → bi-tight → Main; L closed only with full chain.

Preferred: Weil/character-sum on η after Wick, or Aut-irrep closed Φ spectrum.

## Artifacts

- `src/e1_gmin_m4_prop15165.py`, `tests/test_prop15165.py`
- `evidence/e1_gmin_m4_prop15165.json`
- `evidence/gpu_es4_angles_2h.json`
- `evidence/gpu_m4_orbits_2h.json`
- `evidence/gpu_m4_by_Csig_2h.json`
