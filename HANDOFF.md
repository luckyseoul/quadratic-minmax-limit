# Handoff: min-max ±1 quadratic form

**Date:** 2026-08-20 (updated; prior 2026-08-18 content preserved below)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit (`main`)
**Statement:** [MathOverflow 413935](https://mathoverflow.net/questions/413935).
\(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\lvert\sum_{i<j}a_{ij}x_ix_j\rvert\).

## 2026-08-20 post-d470078: exact PSL multiplicity reduction

Work on `codex/leftover-moment-attack` gives the multiplicity-free
decomposition

`Z = W_e direct-sum ((p^2-9)/8 distinct degree-(p^2+1) principal series)`,

where `dim W_e=(p^2+1)/2`. Thus every Phi eigenvalue has multiplicity at
least `n=p^2+1` except one exceptional scalar. That scalar is exactly

`lambda_exc = 32 E|Z_psi|^2/[q(q-1)]`, `psi^2=chi`,

so its floor is the single quartic variance inequality
`E|Z_psi|^2 >= 3q(q-1)/16`. The remaining principal floor route is the
existing delta room `||delta||^2 <= n(n+10)^2/[6(n-6)^2]`. Neither inequality
is proved generally; no flag is flipped. See
`evidence/NOTE_2026-08-20_psl_and_stratum_floor_reduction.md` and Prop 15.589.

Two routes are now explicitly dead: restricted Phi does not have floor 6 on
every profile stratum, and at p=7 a full PSL orbit of size 1,176 has
`Z_psi=0`, so the quartic bound is not pointwise/orbitwise.

**Settled.** Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) (`solution.md`). Paley \(\rho=1\) on \(n=p^2+1\) (`evidence/PROOF_rho_eq_1.md`).

**Not settled.** \(L=\lim\alpha_n=1/2\). That needs E(1): Paley \(C\) is a \(\Phi\)-minimizer for every prime \(p\ge5\). Soft-close from sandwich plus denseness, without E(1), is not a proof.

Live `e1_closed_general` is True only by the old wiring (affine residual (ii) plus two-level Type I). That is not E(1) for all \(p\ge5\). Aut-Schur, Gsum disj LB, and the cotangent pairing stay unused and False.

## 2026-08-20 session — read this first

**No flag flipped. `Max+` at `p=11` is now fully enumerated**, breaking the wall
`fable.md` names as the likely common blocker to leftovers 1 and 3 ("Max+ is
enumerable only for `p<=7`"). Full derivations, exact data, and eight corrections
to claims made mid-session: `evidence/TECHNICAL_NOTES_2026-08-20_maxplus_p11.md`.
Read `evidence/METHOD_why_500_props_never_moved_a_flag.md` before writing a new
numbered proposition — it names the loop that produced most of the ~500 that
never moved a flag, and gives four rules to avoid repeating it.

**Leftover 1** (row below still says "Open, ensemble unnamed" — still true, but
narrower now): reduces to `var(spec Phi) <= 32(n+10)^2/(n-6)^3` together with
`mult(lambda_min) >= n`. Both hold at `p=7` (+32.7% slack) and `p=11` (+91.0%
slack); `p=5` stays a finite check. Proved: no bound using only `dim Z, tr(Phi),
tr(Phi^2)` can ever close this — the best such bound is `10-16x` too weak, and the
gap widens with `p`. The variance is computable from leftover 3's own four-point
tensor: `tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2`, verified exactly at
`p=5,7,11`. **Single remaining target for this route: bound `||M||_F^2` at
general `p`.** Detail: `evidence/NOTE_leftover1_variance_multiplicity_route.md`.

**Leftover 3**: `mu/L` slack grows `0.769 -> 0.746 -> 0.307` (i.e. **growing**
headroom) from `p=5,7,11`. No structural change, just a third confirming point
with more room, not less.

**k=4 stratum of Max+ terminates at p=19** (confirmed both by GPU and by
independent uncapped CPU DFS on all three `p=19` subsets, with a `p=17` positive
control matching exactly — see technical notes §6). Not itself a leftover, but was
this session's original target before the enumeration wall turned out to be the
more useful thing to attack.

**Correction to the "Bottom multiplicity is exactly n" claim two sections below**:
false at `p=11`, where it is `244 = 2n`. The weaker `mult >= n` form survives and
is what the leftover-1 route above uses.

**Repo housekeeping**: GPU code that made the p=11 run possible
(`gpu_inner.py`'s auto-splitting flip resolution, memory-pool capping) had been
sitting **uncommitted on a tmpfs RAM disk** for a full session — a reboot would
have destroyed it. It is committed now. Branch `prop15586-maxplus-gram-reduction`
(`1fa0301`) is an ancestor of `main`; nothing from it was lost.

**Conflict worth knowing about**: the "Do not commit... 15.496, or 15.530" line
further down was written by a prior session and left unexplained beyond being
grouped with other dead ends. Both are now tracked on `main` per an explicit
live instruction in the 2026-08-20 session, with 15.530's two failing tests
quarantined as `xfail(strict=True)` rather than left red. If you are the one who
wrote that original instruction and had a reason beyond "these are dead ends",
that reason is not recorded anywhere — leaving this note so it is not silently
overridden twice.

## Open for E(1)

| Item | Flag | Status |
|---|---|---|
| \(\lambda_{\min}(\Phi\rvert_F)\ge6\) \(\Leftrightarrow\langle\delta,\psi\rangle\le2\) | `phi_F_ge_6_proved_general=False` | Open. Ensemble \(Q_\tau\) (or \(D=\lvert H_+\rvert/(2p)\)) unnamed in \(p\). |
| Residual (ii), even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | Open. Affine and even \(k\le4p-2\) are closed (15.179, 15.236, 15.237). |
| Type I, Max− not two-level \(\{-1,-3\}\) | `type_I_multilevel_bad_case_ND_closed=False` | Open. Remainder is \(A_{\mathrm{full}}\). |
| Lemma D | True | Closed. Do not unflip. |

**Next identity.** Name ensemble \(Q(r)=E[u(\xi)u(r\xi)]\) on the 15.290 Paley×norm types by a Gauss/Jacobi (or equivalent) formula in \(p\), with a fail-when-wrong test; or prove \(\langle\delta,\psi\rangle\le2\) for every even \(\psi\notin\{1,\chi\}\) without that name. Then import `phi_F_ge_6` only if the pairing bound follows. Do not import from a \(p=5\) evaluation or a two-point interpolant.

## Floor (leftover 1)

Wick: \(Q(\pm1)=8q^2\), off-diagonal \(4q^2\). \(\delta=4-Q/q^2\). Floor \(S_\square\ge6q^2\Leftrightarrow\langle\delta,\psi\rangle\le2\).

Live ensemble (not a general proof): \(Q_{++}/q^2=48/13\) at \(p=5\), \(1544/409\) at \(p=7\). \(Q=8A/D\) with \(D=\lvert H_+\rvert/(2p)=13,409\) (\(2^2+3^2\), \(3^2+20^2\); not a polynomial in \(p\)).

Named pieces:

- \(S(\lambda)=\mathrm{Kl}(1,\lambda^2/4)\) (15.550).
- \(F=-2(3p^2+2)\), \(Q_{3,02}=-4N(2p^2+1)/p\) (15.564).
- \(n_{1d}=m\binom{p}{m}\), \(n_{k=3}=\binom{m}{3}(p-1)q\), \(A_{1d}=-4p^3/(p-2)\).
- \(\mu_{k=3}=96p^4 P(r)/(p^2-1)\) (15.574). \(\mu_{1d}=2p^4(p^2-3p-2)/(p-2)\) (15.575).
- Exclusive 1D / \(k=3\) / full mix reconstructs live \(Q\) on every Paley×norm type at \(p=5,7\) (15.573). It is not a \(p\)-identity: the 1D+\(k=3\)-only mix is \(4.68>4\) at \(p=7\); \(\mu_{k=3}/q^2>4\) at \(p\ge11\).
- At \(p=5\), \(n_{\mathrm{full}}=0\), so \(Q_{++}/q^2=48/13<26/7\) (15.581; 15.507 \(J_{N^*}=2\)).
- 1D 4-point vanishes for \(r\notin\mathbb F_p\); \(p=5\) \(Q_{N^*}/q^2=32/13\) (15.582).
- Pointwise \(Q_y^{++}\le4q^2\) is false: about \(23\%\) of Max+ at \(p=5,7\) have \(Q_y/q^2>4\) (max \(5.33\), \(16\)). Any identity that uses only \(z_i^2=1\) and \(Cy=py\) cannot force the ensemble bound.
- Paley×norm types split into many \(\langle\mathrm{Frob},\mathrm{inv}\rangle\) orbits at \(p\ge11\) (++sub has \((p-3)/2\) orbits of size \(2\)). Two-type constancy is certified only for \(p\le7\).
- \(\mathrm{Gal}(\mathbb F_q/\mathbb F_p)\) acts on \(H_+\) with orbits of size \(1\) or \(2\). Orbit masses are \(1/\lvert H_+\rvert\) or \(2/\lvert H_+\rvert\), i.e. they name \(D\).
- \(n_{\mathrm{full}}=\lvert H_+\rvert-n_{1d}-n_{k=3}\) is \(0\) at \(p=5\) and \(90q\) at \(p=7\). \(\mu_{\mathrm{full}}\) is not a single formula in \(p\) (15.578).

## Residual (ii)

Official class is leftover Max− together with \(s_+\ge2\). leftover-only (\(s_+=0\)) at \(p=5\), \(k=20\) exists and is not this class.

- leftover+\(s_+\) empty for all \(n_F\) at \(p=5\), \(k=20\) (15.528).
- 15.585: leftover+\(s_+\) at \(k=4p\) forces \(\min_+=2\); \(\{2,4,6\}\) cannot have \(1_{S=2}\) a plus pair-slice.
- No identity that leftover+\(s_+\) is empty at every even \(k\ge4p\).

## Type I

Two-level Max− is closed (15.272). Multi-level is open. Dead as a multi-level kill: Aut\(_e\) (15.559), Max± of \(C\) (15.565), Type+ 1D Johnson (15.577), Galois support plus \(F\) (15.580), \(\lvert\mu\rvert\le\lvert L\rvert\) on \(\lvert\kappa\rvert=1\) (unsigned \(\lvert\nu_{\mathrm{part}}\rvert\) exceeds \(\lvert L\rvert\)). Remainder is \(A_{\mathrm{full}}\).

## Do not reopen

Occupancy / Aut-involution pairing of \(T_{\mathrm{ns}}\) / \(\bar n_0\) interpolants / half-net census as a \(p\)-law / Aut\(_e\) as a name of \(A_{\mathrm{full}}\) / \((p-5)/15\) / \(10p-46\) / \(16(p-4)/D\) / Paley type as a \(Q\)-constant (false at \(p=7\)) / exclusive mix as a general \(Q_\tau\) / pointwise Wick or Boolean collision as a proof of \(Q_{++}\le4q^2\) / Gsum as a Gram / Aut-Schur.

15.495 catalogs, 15.496, and 15.530 **are now committed** (2026-08-20, see note
at top of file) — this line originally said not to. Left visible rather than
deleted so the reversal is traceable.

## Files

| File | Role |
|---|---|
| `STATUS.md` | Claim table |
| `GOAL.md` | Acceptance for E(1) / \(L=1/2\) |
| `solution.md` | Sandwich; Main Theorem (limit) stays OPEN |
| `evidence/share/denseness_path_package.md` | Stand-alone path; § Caveats |
| `evidence/SESSION_HANDOFF_2026-08-18_leftovers.md` | Named identities 15.550–15.585 |
| `evidence/TECHNICAL_NOTES_2026-08-20_maxplus_p11.md` | p=11 enumeration, derivations, corrections (2026-08-20) |
| `evidence/NOTE_leftover1_variance_multiplicity_route.md` | Leftover 1 reduction, single remaining target |
| `evidence/METHOD_why_500_props_never_moved_a_flag.md` | Read before writing a new numbered proposition |
| `evidence/maxplus_p11/` | Scripts + logs for the p=11 spectrum/moment computations |
| `LONG_HORIZON_GOAL.md` | Terminal states |

Large `.npy` arrays (Max+ at p=11, 4.5 GB each) are **not in git** — they live at
`/mnt/storage/e1work/maxplus_p11/` on soulkiller, verified by md5 against the
original computation. Scripts there have hardcoded `/tmp/e1work` paths; repoint
before rerunning.

**HEAD:** `d877359` (tracks stranded props + method note; ancestor chain includes 15.585 → 15.588 → p=11 enumeration).
