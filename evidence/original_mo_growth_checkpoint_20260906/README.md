# Low-quota growth-search checkpoint

2026-09-06 09:17 UTC. The user reported only 3% quota remaining.
New hunts and the proposed delayed-repair path search were stopped before
launch. The research subagents were interrupted. This is a preservation
checkpoint, not a solved theorem or a new large-drive backup.

The original question remains OPEN: whether m_n/n^(3/2) converges.
No new uniform all-orders upper/lower comparison was obtained from these
hunts. The existing reviewed bounds in CORE.md are unchanged. Growth,
repair, conference optimality, and a particular limit value are optional
routes, not necessary conditions for a solution.

## What was actually obtained

- A complete two-vertex extension identity and O(3^n) integer max-plus
  optimizer, with independently checked CPU/CUDA/OpenCL implementations.
  The optional full endpoint/intermediate tie optimization can cost O(4^n).
- A fixed order-11 greedy trap: 17 -> 18 cannot do better than 24 on the
  second addition, while 17 -> 20 -> 22 attains the complete-family optimum.
- Two norm-27 order-15 sources have the same optimal order-17 endpoint 34,
  but their best endpoint-optimal intermediate norms are 34 and 30.
- ONE retained order-64 signing has exactly Phi=248, normalized 31/64.
  This is not a proof that m_64=248 or an asymptotic construction theorem.
- All 32 tested grown order-19 sources have exact two-add endpoint minimum
  48. The twelve norm-39 sources are one explicitly verified equivalence
  class, not twelve independent construction types.
- A reverse source from a known C26-derived norm-44 order-21 target has
  norm 41 and exact two-add optimum 44. Its path 41 -> 44 -> 44 decreases
  at both steps after normalization. Restriction is not a forward bridge.
- At target 44, retained support capacities are 6+5<19 for the norm-39
  control, versus a compatible 9+10 split with seven complementary pairs
  for the reverse control. These are exhaustive table consequences, not
  analytic capacity bounds for arbitrary sources.
- A proved rank/inertia lemma gives a structural barrier between C18 and
  C26-derived targets: d+k>=5, with d original vertices absent and k final
  edge mismatches among surviving originals. This is not a path-loss bound.
- The actual norm-17 warm A11 does not embed in the specified reverse A19
  or original A21 even up to switching, permutation, and global negation.

## Last completed result, after the longer outline

One exact signed-injection optimization found that the minimum number of
old-edge changes needed to embed this A11 in this A21 is TWO. All 42
root/global-sign cases completed: 578,965 search nodes in 0.146 seconds;
4,120 independent brute comparisons passed, including 1,536 positive-distance
instances. One optimal witness changes zero-based edges (3,8) and (5,10).

Receipt inside the archive:
`original_mo_path_v100_20260906_Q86jOz/signed_injection_repair_results/result.json`
SHA256 `f8f1ba92c01be0ba2e951f601733325e08e22d29bc57b45a6053ecde0235c544`.
It includes the actual source, target, signed injection, and full aligned
order-21 target. No intermediate norm or optimal repair schedule is claimed.

The proposed next finite check was a fixed-witness state graph (S,F), with
ten target vertices and two required old-edge repairs: 4,096 states covering
all addition/repair orders, scored by exact rational Phi^2/n^3 minimax.
IT WAS NOT LAUNCHED. It would optimize one finite path, not solve the
original convergence problem. Do not automatically restart it or another
maze hunt without respecting the user's remaining-quota constraint.

The main mathematical need remains a genuine argument valid at arbitrarily
large orders. None of the finite failures above establishes impossibility
of the original problem, and none of the successes establishes proximity
to a solution. Do not repackage a Dini/amplification sufficient condition
or the existence of isolated good targets as the missing proof.

## Preservation and verification

`checkpoint.tar.gz` preserves the exact order-64 certificate bundle,
two-extension proof/CPU results, source pool and four backend receipts,
V100 scratch/results, reverse/embedding/repair evidence, and the bounded
Gram/support diagnostics. Original temporary directory names are retained
as relative archive roots. The archive was compared against those frozen
source files; see `archive_receipt.json` for its hash and size.

The archive's `original-mo-fixed-target-path.2amt4q/` contains the detailed
reviewed proof notes and an older milestone outline. This README supersedes
that outline about the final two-edge repair result and quota stop. Scratch
files in the archive are preserved, not blanket-certified as theorems.
The fixed-target DP proof is a method only; no DP run was made for this pair.

Some wrapper scripts retain original absolute paths and provenance hashes.
To replay, extract to a fresh directory and supply the archived inputs
explicitly; do not overwrite the original /tmp directories. The independent
field verifier source is included in the grown19-pool tree. Some historical
raw-journal provenance dependencies are not bundled; the checked receipts
retain their source records and hashes. This is not a full archive of every
historical raw search log or every temporary file.

Large raw scratch remains separately in `/tmp/original-mo-pathwalkers.0Suy9n`
(about 4.7 GiB), `/tmp/mo-nuka-pathwalk-20260906.AvNLid` (about 577 MiB),
and `/tmp/original-mo-resident-batch.Nuy2te` (about 80 MiB). Nothing was
deleted. The earlier major campaign backup recorded in HANDOFF.md remains
unchanged; no fresh major backup is claimed here.

At creation, branch main remained at
`902539250598fd21e78a1d9a1f6dbf0233ed0f16`. Seven previously staged metadata
changes were preserved. The checkpoint was initially left uncommitted;
the user subsequently requested that it and the pending metadata be committed
and pushed. The Git commit containing this directory identifies that
publication. The archive receipt's commit/push flags describe its creation
time, not the later publication. No new full regression run was performed.
