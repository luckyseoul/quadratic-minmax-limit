# GPU session (V100, ~2h window) — 2026-08-05

**Backend:** CuPy 14.1 / Tesla V100-SXM2-16GB, free ~15 GiB at start.  
**L status:** OPEN. `residual_closed_general=false`. No soft-close.

## Compute used (real device wall, not theater)

| Job | Result |
|-----|--------|
| `e1_residual_h_close.certify_live_H` p=5,7 | ray≤H, 16N, G≽I on **CuPy GPU** |
| Full Φ|Z eigensystem + maximizer ID | p=5 wall ~2s; p=7 ~31s |
| Haar sampling unit sphere in Z | p=5: 1e5 samples; p=7: 4e4 (batched) |
| Max+ Gram / angle multiset | p=5,7 on GPU matmul |
| p=11 affine Max+ sample (2772) | Φ on sample **not** a certificate (λmax≃240; incomplete Max+) |

## Certified (unchanged)

- p=5: ray = H = 49/13, λ_max = 176/13 = thr, mult = d = 13, G eigs {1,2,4}
- p=7: ray = 933/409 < H = 81/25, λ_max = 4320/409 < thr, mult = d = 25, G_min > 4

## Structure findings (useful; not a residual close)

1. **CBC / CB+BC trivial on Z**  
   For B ∈ Z ⊂ End(V₊): CBC = p² B, CB+BC = 2p B. Do not distinguish Φ levels.

2. **Edge zero-diag rank-1 family sits at μ̄**  
   All sampled edge zd projectors have Rayleigh = μ_bar (p=5: 9.6; p=7: 8.727…). Not maximizers.

3. **Maximizers**  
   Full rank on V₊, pairwise non-commuting, F-orthonormal d-plane, **orthogonal to C** in Frobenius (not in conference BM algebra). Ambient B annihilates V₋ (d zeros) as expected.

4. **Haar on Z is far from the maximizer**  
   p=5: Haar max ray ≃ 2.78 ≪ exact 3.77 = H; mean ≃ 1.80.  
   p=7: Haar max ≃ 1.53 ≪ exact 2.28 < H=3.24; mean ≃ 1.36.  
   Extremal ray lives in a thin d-dimensional locus.

5. **Max+ angle multiset (GPU Gram)**  
   - p=5: 9 off-diag angles {±2,±6,±10,±14,−26}; **distance-regular** from every point; valencies  
     `1,13,20,36,60 | 60,36,20,13,1` (self / antipode); v(±14)=d=13.  
   - p=7: 17 angles; **not** distance-regular (per-point type split); antipodal pair always present.  
   - Both: s ≡ 2 (mod 4), closed under negation; ∑Max+ = 0; 2-design moment ∑_z s² = N n²/d holds at p=5.

6. **||κ_B||² not constant on unit Z**  
   Varies on ON basis (~0.5%–few %); maximizers sit near the upper end. CS-with-constant-κ dead.

7. **Aut(C) on maximizer space (p=5 sample)**  
   AGL elements that preserve C act orthogonally (det=1) on the d-plane; consistent with PSL d-irrep embedding story (15.159.D), not a full character table.

8. **p=11**  
   Full free-coord enum impossible (nul=61). Affine halfspaces give 2772 Max+ vectors; sample-Φ is not comparable to full Max+ Φ.

## Algebra (CPU, accompanying)

- `thr − (6+2H) = 2(p−5)(3p+7)/(p²+1)` (matches 15.160 thr_ray−H up to the ray scaling).
- Majorization with only mult≥d and λ_bulk≥6 is **too weak** for H (needs bulk LB ~8.6 at p=5; actual λ_min=80/13≃6.15).
- H−α = (p⁴+8p³+5p²−2)/(p²(p²+1)).

## Still OPEN (acceptance)

Prove for **all primes p≥5** (Max+-free):

`ray_max ≤ H(p) = 2(p+2)²/n = (p+2)²/d`

(or equivalent δ²≤room_hyp/24 / orth≤room / 16N direct), then wire 16N→bi-tight→E(1)/Main and close L in docs.

Preferred structural targets remain: dual gap G≽I closed form; Weil/μ_G4≤μ_G4_suf; mult≥d embedding + κ²≤96n; Max+-free m₄ on G-orbits.

## Artifacts

- `evidence/e1_residual_h_close.json` (GPU cert p=5,7)
- `evidence/gpu_residual_structure_2h.json`
- `evidence/gpu_haar_ray_samples.json`
- `evidence/gpu_maxplus_angles.json`
- `evidence/gpu_p11_sample_phi.json` (non-certifying sample)
- `/tmp/gpu_residual_2h/Bmax_p{5,7}.npy`, `Phi_evals_p{5,7}.npy`, `haar_rays_p{5,7}.npy`
- `/tmp/maxplus_p5.npy`, `/tmp/maxplus_p7.npy` (full); `/tmp/maxplus_p11_sample.npy` (affine only)
