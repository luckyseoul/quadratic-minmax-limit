# Agent bug report: F17 single-core thrash (recurrence)

Filed: 2026-07-29T03:04:05.833052Z
Session feedback: `/home/nick/.grok/sessions/%2Fhome%2Fnick/019f9af7-3128-71c1-984e-2a7102bec72d/feedback.jsonl`

BUG / agent compute policy (F17 recurrence): Grok agent repeatedly pegs a single core on multi-minute research jobs despite prior user rebukes and an explicit project rule (E1_FAILURE_GRAPH F17, workers.py require_workers).

Concrete failures this session:
1. e1_gmin_m4_pseudo.py ran at ~97% CPU, NLWP=1 for many minutes (serial p=7 path on main after a small ProcessPool for p=3,5 only).
2. Earlier in the same goal: same F17 pattern promised fixed twice; user had to re-intervene ("pegging a single core again?", "I shouldnt need to keep telling you").
3. require_workers() only gates worker *count availability* — it does not prevent the agent from writing scripts that (a) spawn a pool then (b) run the heavy leftover work serially on the parent, or (c) pure-Python O(binom(n,4)) loops on one core while 87 cores idle.

Expected:
- Any job expected >~10s CPU on this 88-core host must use real fan-out (ProcessPool/xdist W≈nproc-2, OMP=1 per worker).
- Serial fallback on main after parallel phase is BANNED unless the log explicitly marks the algorithm inherently serial.
- Agent should self-check with `ps` (pcpu + nlwp) after launching heavy work; kill+rewrite if NLWP=1 and pcpu≈100 while nproc>>4.

Impact: user trust destroyed; wasted wall time and compute; blocks progress on MathOverflow 413935 residual proof (P0).

Workspace: /home/nick/quadratic-minmax-limit
Session: 019f9af7-3128-71c1-984e-2a7102bec72d
Related: use-available-compute skill, F17, scripts/pytest_full.sh, src/workers.py

## Recurrence 2026-07-29 ~03:44 UTC

**Job:** `src/e1_gmin_m4_midproof.py` → `_job_hstar_budget(p=7)`  
**Symptom:** PID at **99.9% CPU, nlwp=1**; load ~1.3 on 88 cores.  
**Cause:** "multi-worker" only ran **2** top-level jobs (p=5 and p=7), each building sparse `T` with a **serial pure-Python loop over all 4-sets** + serial `eigsh`/`lsqr`. When p=5 finished, one orphan process pegged one core.  
**Not fixed by** calling `require_workers()` at main — W was unused for the heavy path.  
**Fix required:** shard COO build of T across W workers; never leave a multi-minute pure-Python `for S in combinations` on one process without concurrent independent jobs filling the machine.

