# Session handoff — 2026-07-28 (evening)

> **Superseded for resume:** use `evidence/SESSION_HANDOFF_2026-07-29.md` (post-shutdown handoff).

**Repo:** `/home/nick/quadratic-minmax-limit`  
**HEAD:** `8f4de72` (Prop 15.49); later `5ace906` handoff commit  
**L = lim α_n:** still **OPEN**

## Settlement chain

\[
m_n\ge\Phi-2\text{ on }n=p^2+1
\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12
\]
via denseness (Prop 6.2). Needs: (i) bi-tight empty all \(p\ge5\); (ii) deep non-tight ND / no undercut.

## Best residual edge (Prop 15.49)

**Uniform LB candidate (proves bi-tight empty if established):**
\[
g_{\min}\;\ge\;L(p)\;:=\;-\frac{p-2}{2p^2}
\]
- Algebra: \(L(p)>\) bi-tight threshold for all odd \(p>2\).
- Certified: holds at \(p=5,7\); fails at \(p=3\) (correct).
- Evidence: `e1_gmin_uniform_lb.json`, `src/e1_gmin_uniform_lb.py`.

**CR classification (certified p=3,5,7):**
- \(g_{\min}=-\alpha_\star\) on constant-\(m_4\) classes with \(|\kappa|=1\).
- Values: \(-\frac13,-\frac3{65},-\frac{109}{2863}\).
- Evidence: `e1_gmin_cr_classify.json`, `src/e1_gmin_cr_classify.py`.

**G spectrum:** rank \(\binom{d}{2}-d+1\); simple eigenvalue \(n/2\); full spectrum known at p=5,7.

**Matching margin (p=5):** \(\mathbf1_M^\top G\mathbf1_M\ge9.96>4\) for all sampled matchings of size \(2p\).

## Dead ends (do not reopen)

4-point LP; Chebyshev; Wick-as-LB; \(-3/\Phi\) general; bare C-types; affine halfspace orbit; pure deg pigeon; min-norm V+ alone; star tautology; Ising maxent; incomplete Max+ samples for gmin.

## Next steps (priority)

1. **Prove \(g_{\min}\ge-(p-2)/(2p^2)\) for all primes \(p\ge5\)**  
   Routes: character sum on min CR class; association-scheme eigenmatrix; joint residual CLT with error bound; SOS on 4-point marginal with evec extension.
2. **Deep non-tight:** ND or always \(\Phi\ge\Phi\) for \(p\ge5\) (p=5 certified spike).
3. **Only then:** Main Theorem \(L=\tfrac12\); full pytest; HANDOFF settled.

## Commands

```bash
cd /home/nick/quadratic-minmax-limit
# caches: /tmp/maxplus_p5.npy, /tmp/e1_p7/maxplus.npy
python3 src/e1_gmin_cr_classify.py
python3 src/e1_gmin_uniform_lb.py
OMP_NUM_THREADS=1 python3 -m pytest tests/test_minmax.py tests/test_gmin_residual.py -q
```

## Do not

- Soft-close Main Theorem
- Claim \(m_n\ge\Phi-2\) from 15.40 alone (F13)
- Re-run dead loops above
