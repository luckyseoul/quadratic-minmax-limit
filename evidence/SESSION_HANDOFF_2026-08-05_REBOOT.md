# Session handoff — 2026-08-05 (reboot)

**Why this exists:** User rebooting host to restore Tesla V100. Next agent continues residual/16N work without replaying the full transcript.

**Workspace:** `/home/nick/quadratic-minmax-limit/`  
**Master research handoff (do not re-read whole chat):** `HANDOFF.md` (status through Prop **15.166**)  
**UTC stamp:** 2026-08-05T07:16Z

---

## 1. Goal (unchanged)

Close lim α_n for MO 413935 via Path C:

**residual → 16N → bi-tight → E(1)/Main**

Acceptance: residual for all primes p≥5; bi-tight via 16N (or g_min>T); E(1)/Main; **L closed** in `solution.md` / `HANDOFF.md` **only if full chain holds**.

**Forbidden:** soft-close (F3); class_key thrash (F19); GPU theater (F20).

**Totality status:** **NOT MET.** `residual_closed_general=false` everywhere. **L OPEN.**

---

## 2. What shipped this session (Props 15.161–15.166)

| Prop | Module | Content | residual_closed_general |
|------|--------|---------|-------------------------|
| 15.164 | `src/e1_gmin_m4_prop15164.py` | mult≥d−1 + `E[s⁴]≤Es4_*` ⇒ 16N | **false** |
| 15.165 | `src/e1_gmin_m4_prop15165.py` | closed Es4_*/η_*; exact Es4 p=3,5,7; GoG↔Φ; m4 C-eigen; p=7 **not** 1-homog | **false** |
| 15.166 | `src/e1_gmin_m4_prop15166.py` | Max+ unit 2-design; Wick C-eigen; **16N ⇔ λ_max(Q₂)≤4N/(d(d−1))** | **false** |

Tests (green last run): `tests/test_prop15164.py` … `test_prop15166.py` (xdist 86 workers).

Evidence JSON: `evidence/e1_gmin_m4_prop1516{4,5,6}.json`  
Dead-end log: `evidence/RESIDUAL_ATTACK_2026-08-05_GPU_B.md`  
GPU session notes: `evidence/GPU_SESSION_2H_B_2026-08-05.md`

**Single open inequality (equivalent forms):**
- `E[s⁴] ≤ Es4_*(p)` with mult≥d−1 already proved, **or**
- `λ_max(Q₂) ≤ 4N/(d(d−1))`, **or**
- `∑η² ≤ η_*`, **or** `δ² ≤ room_hyp/24` (older δ form).

Es4_exact: p=3 equality with star; p=5 `120400/13`; p=7 Gram/Φ (not W_CENSUS).  
Ratios: Es4/Es4_* ≈ 0.996 (p=5), 0.965 (p=7); Q2/thr ≈ 0.846 / 0.660.

---

## 3. Dead attacks (do not retry)

Spherical Delsarte LP; two-level angle majorization; C-sig/CR alone; BM(C) (Bmax ⊥ I,C,J,C²); Gershgorin; Gegenbauer dual a_{k≥4}≤0; CS on ∑ρ κ_B; Q2 majorization with λ_min only; entrywise |δ|≤3/p²; P₊(Wick)=true m4; Wick-orth⇒δ² mismatch; continuous Veronese without finite-design defect control; class_key as m4-equitable at p=7.

---

## 4. Preferred next attacks (still open)

1. **Weil / Jacobi-sum** evaluation of Aut-orbit m4 on PG(1,p²) + N ≥ |PSL|/|Stab| LB  
2. **Aut₀-isotype** V₊↪V_max with closed Rayleigh of the 1-dim line  
3. **SOS / Putinar** on Q₄ ≤ 10N‖B‖² uniform in p  
4. **Equidiagonal rank-d projector** λ₂ bound with theorem N ≥ c d² for Max+  
5. Parallel bi-tight: prove `|m4| ≤ (p−2)/(2p²)` on |κ|=1 (g_min path; census only p=5,7)  
6. Or **abandon residual** and close E(1)/deep ND on n=p²+1, then wire Main

---

## 5. GPU / compute at leave-off

| Item | State |
|------|--------|
| Hardware | Tesla V100-SXM2-16GB @ `0000:03:00.0` (PCI visible) |
| Driver | NVRM 580.173.02; `/dev/nvidia*` present |
| Usability | **BROKEN:** `nvidia-smi` → “Unable to determine the device handle”; CUDA `cudaErrorNoDevice` |
| Cause | Stale/hung GPU after prior jobs (`realesrgan-ncnn` held device; freed later). **User rebooting to restore.** |
| CPU | 88 threads, ~55 GiB avail, full_workers=86 — fine for non-GPU work |
| Budget script | `~/.grok/skills/use-available-compute/scripts/compute-budget.sh` |

**After reboot — first actions:**
```bash
nvidia-smi -L
python3 -c "import cupy; print(cupy.cuda.runtime.getDeviceCount())"
~/.grok/skills/use-available-compute/scripts/compute-budget.sh
```
If GPU is healthy, rebuild Max+ caches (see §6), then resume residual attack with **real** GPU wall (no theater).

---

## 6. Caches lost on reboot (/tmp)

These will **disappear** on reboot; regenerate before GPU census:

| Path | Content |
|------|---------|
| `/tmp/maxplus_p3.npy` | Max+ rows p=3 |
| `/tmp/maxplus_p5.npy` | Max+ p=5 |
| `/tmp/maxplus_p7.npy` | Max+ p=7 |
| `/tmp/maxplus_p11_sample.npy` | sample only |
| `/tmp/e1_p7/maxplus.npy` | alternate p=7 path used by some older scripts |

Look for builders under `src/` / `scripts/` that dump Max+ (prior sessions used Paley Max+ enumeration into those paths). Do **not** point package caches / HOME at scratch dirs.

Scratch from goal harness (if still present):  
`/tmp/grok-goal-02c2ba1753a2/implementer/` — logs: `residual_close.log`, `bitight_chain.log`, `main_theorem.log` (all report OPEN).

---

## 7. Honest flags (must stay false until real proof)

```text
residual_closed_general = false   # e1_residual_h_close, prop15164–166
L = OPEN                            # HANDOFF.md, solution.md
```

Never set true without general-p math + tests.

---

## 8. Next concrete steps (post-reboot)

1. Confirm V100 healthy; re-budget compute.  
2. Rebuild Max+ caches p=5,7 (and p=3 if needed).  
3. **Ship one residual-close path** for all primes p≥5 as pure-math functions + pytest:  
   Weil/Aut m4 **or** equidiag λ₂+N≥c d² **or** SOS Q₄  
   — **or** explicitly abandon residual and close E(1)/deep ND.  
4. Only then: wire 16N → bi-tight → Main; update `solution.md` + `HANDOFF.md` so L is not OPEN.  
5. Verify with goal-verifier: residual/bi-tight/Main logs must show **close**, not OPEN.

Do **not** open new prop numbers that only restate census/dead LPs.

---

## 9. Suggested skills

- `use-available-compute` — mandatory before heavy/GPU jobs  
- `graph-engineered-completion` — residual attack graph; avoid F19 thrash  
- `goal-verifier` — when claiming chain close  
- `handoff` — if context fills again  
- optional: `check-work` before claiming shipped props green  

---

## 10. Pointers (read these, not the chat)

| Path | Why |
|------|-----|
| `HANDOFF.md` | Full prop history through 15.166; L OPEN |
| `evidence/RESIDUAL_ATTACK_2026-08-05_GPU_B.md` | Dead ends R1+R2 |
| `evidence/E1_FAILURE_GRAPH.md` | F1–F20 |
| `src/e1_gmin_m4_prop15164.py` … `15166.py` | Current algebra API |
| `src/e1_residual_h_close.py` | residual_closed_general gate |
| `src/e1_gmin_m4_bound.py` / `e1_gmin_m4_proof.py` | |m4|≤L on |κ|=1 bi-tight alternate |
| MO 413935 | External problem statement |

---

**End of reboot handoff.** Resume at §8 after GPU is back.
