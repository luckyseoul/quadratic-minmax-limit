# Research rules

## The actual question

Determine whether `alpha_n = m_n / n^(3/2)` converges, where
`m_n = min_A max_{x in {+1,-1}^n} |sum_{i<j} A_ij x_i x_j|`
and `A` ranges over complete symmetric zero-diagonal signings.

No candidate value or proof architecture is mandatory. In particular,
residual (ii), Paley optimality, a uniform gap of two, and amplification
at multipliers two and three are not necessary conditions established for
this question. The original MO limit is OPEN because there is no reviewed
complete proof, not because a particular route has an open lemma.

## Start here

Read `CORE.md`, `STATUS.md`, `HANDOFF.md`, and `ARTIFACTS.md`.
The old 2026-09-04/05 research instructions are archived, not binding.
Consult `solution.md`, the proposition audit, and the relevant source or
proof note when reusing a result. Their local hypotheses and retractions
matter; historical declarations of a preferred route do not set the goal.

Before an attack, identify one implication it would establish. Search the
existing evidence and code for that exact object and invariant. Reuse
stored certificates. A new backend, larger finite census, longer timeout,
or renamed obstruction is not a new argument.

## Proof and status discipline

- Separate proved theorems, exhaustive finite certificates, open reductions,
  counterexamples, and retracted claims.
- A sufficient condition is not a necessary condition. A conditional
  equivalence inside one route is not a global acceptance criterion.
- Do not assume `L=1/2`, optimizer classification, asymptotic Paley
  optimality, or local-to-global entry bridges.
- Global claims are recorded in `src/original_mo_status.py` with an explicit
  reviewed proof reference. Legacy Paley flags cannot close or veto them.
- Read the complete proof and implementation before running its checker.
  A finite check cannot prove an all-orders assertion.
- Update current status and provenance after a substantive result. Preserve
  old work on artifact branches; never silently erase a counterexample.

## Execution discipline

Preserve dirty or reviewed work. Keep edits scoped and archive before a
structural reset. Do not commit another agent's unfinished changes.
Offload proof computations and test runs from the controller. Use a fresh,
hashed staging directory, explicit test files, and a host-sized worker
budget. Ordinary CUDA belongs on soulkiller/V100; use Orin only for a
capability the V100 lacks. Recheck host availability before dispatch.
Inspect actual failing logs before changing code or test expectations.
Do not rerun unchanged mathematics merely to obtain a cleaner receipt.

## Mandatory backups and milestone publication

Nick explicitly requires large-drive backups across Codex sessions and a
commit and normal push after every reviewed research milestone (2026-09-05).

- Before risky changes and at every reviewed milestone, create a fresh,
  verified backup on the actual large drive at `/mnt/storage`. During a
  long run, also back up material uncommitted progress at least hourly.
- Read and use `scripts/milestone_backup.sh` (also installed as
  `/home/nick/.local/bin/codex-milestone-backup`). It preserves the entire
  working tree, all-ref history bundle, dirty patches and checksums.
  Verify the mount, free space, archive comparison and bundle integrity.
  Request scoped escalation for a sandbox read-only view; never remount
  automatically, silently substitute the OS disk, or overwrite an old backup.
- Quiesce source edits while taking the snapshot. Record its exact path,
  checksum and covered HEAD. Failed/partial snapshots are not completed
  backups and must not be silently removed.
- After verification, commit only reviewed in-scope work, update status
  and provenance, then push to the verified intended remote/branch. Preserve
  unfinished agent work and unrelated dirty files. Never force-push.
  Verify the remote SHA and record the checkpoint and next unresolved target.
- A backup, commit, push, scanner pass or finite certificate does not
  change an open theorem into a proved result. Retain every known gap.
