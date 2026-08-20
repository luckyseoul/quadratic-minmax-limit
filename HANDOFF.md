# Handoff: min-max ±1 quadratic form

**Date:** 2026-08-20 (current branch; use `git log -1` for the exact checkpoint)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit
**Current branch:** `codex/leftover-moment-attack` (use `git log -1` for the
current hash).  The latest continuation adds Prop 15.589 Theorems J--K, the
exact p=11 active-subset and top-profile-degree counterexamples described
below.  All 94 focused Prop
15.588/15.589 tests pass.
**Statement:** [MathOverflow 413935](https://mathoverflow.net/questions/413935).
\(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\lvert\sum_{i<j}a_{ij}x_ix_j\rvert\).

## Binding status after e23edef: exact PSL and high-stratum reduction

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

The exceptional target has since narrowed again.  For `p=3 mod 4`, `Z_psi`
is exactly a signed sum of nonnegative directional profile energies whose
pointwise total is `p(p^2-1)/4`.  Combining that identity with the affine
profile classification, and using the Euler-product lower bound for
`L(2,chi_p)` when `p=1 mod 4`, proves `(QVAR)` on every `k=1` and `k=3`
stratum for every prime.  Thus the exceptional scalar remains open only on
the union `k>=4`; the principal delta-variance target is unchanged.  See the
updated reduction note and Prop 15.589.

The same degree theorem now gives exact arithmetic on every genuine profile:
`a_L in 2p Z`.  With `b_L=a_L/(2p)` and `T=(p^2-1)/8`, one has
`sum b_L=T` and `Z_psi/(2p)=T (mod 2)`.  Thus QVAR is the integer
anti-concentration target `E|sum psi(L)b_L|^2 >= 3T/8`.  Parity alone is far
too weak, but this normalization is binding and explains the p=11 histogram.

There is also an exact lattice-coset reformulation.  Max+ is the first shell
of the odd coset `y0+2 ker_Z(C-pI)`, not the ordinary lattice's first shell
(the latter has explicit norm-`p+1` Baer-line vectors).  The radius-sphere
benchmark for the exceptional quartic moment exceeds QVAR by
`q(q-1)(q-11)/(16(q+5))`.  Thus the live exceptional target is equivalently a
lower bound on one degree-4 odd-coset harmonic coefficient; proving that
coefficient nonnegative would suffice, but is still open.

The full p=11 directional covariance also kills a tempting shortcut: although
the quartic direction is top for the complete mixture, it is bottom among the
nonzero modes on k=4 and not top on k=6.  Stratum invariance alone cannot prove
QVAR by a top-eigenmode argument.

The latest attack kills a broader profile-only shortcut for every
`p=3 mod 4`, `p>=7`.  There are artificial full-support energy ensembles with
the exact conserved total, cyclic directional symmetry, equal means, integer
energies, separately admissible line profiles in every direction, and the
actual divisibility `a_L in 2p Z`, but quartic variance zero or `4p^2`,
below QVAR.  The fake profiles can also satisfy the individual polynomial
degree bound from Prop 15.588.
Therefore those facts cannot prove the bound even in combination.  Any
surviving profile argument must use the cross-direction coefficient kernels
and simultaneous Boolean ridge reconstruction, or an equivalent coupling
among directions (Prop 15.589 I).

The p=11 k=4 pure-parabola census further shows that the live inequality is
not fixed-active-subsetwise.  Each of nine genuine balanced four-direction
families has normalized moment `E B^2=5<45/8`; six unbalanced families have
`E B^2=63` and rescue the count-weighted aggregate to `39/2`.  Therefore the
proof must mix projective direction configurations even before mixing profile
strata (Prop 15.589 J).

The full-support coefficient descent kills another tempting induction.  At
`p=7` the top-degree-zero class is empty and every nonzero class clears QVAR.
At `p=11`, however, all `2,090,880` top-degree-zero vectors have actual profile
degree exactly three and moment `E B^2=137/36<45/8`; each of their twelve
projective leading-coefficient classes also fails.  The ten degree-four
nonzero classes each have the identical moment `111483/14039>45/8` and rescue
the exact mixture to `114771/14903`.  Thus QVAR cannot be proved separately by
actual profile degree or leading-coefficient class; adjacent degree families
must be mixed in their exact ensemble proportions (Prop 15.589 K).

There is also one positive high-prime closure.  A centered-Fourier bound for
nonconstant quadratic line profiles, plus six exact two-character-class
checks at `p=41,43,47,53,59,61`, proves that every active quadratic profile
uses more than one quarter of the conserved profile energy for every
`p>=41`.  Four active profiles are therefore impossible: the `k=4` stratum is
empty for all primes `p>=41` (Prop 15.589 L).  This does not remove `k>=5`.

The same centered-Fourier argument plus Weil's additive-character estimate
works at arbitrary profile degree: a `k>=4` stratum is empty whenever
`p>4k^2`.  Hence every surviving stratum is `k=1`, `k=3`, or
`k>=sqrt(p)/2`; the first two already satisfy QVAR.  The unresolved
exceptional scalar is therefore asymptotically a high-activity problem
(Prop 15.589 M).

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

**Leftover 1, current form.**  The older variance-plus-multiplicity route below
has been refined by Prop 15.589.  Multiplicity is now proved exactly: every
principal block has degree `n`, and the only smaller block is one exceptional
degree-`n/2` scalar.  The live route therefore has two separate targets:

1. exceptional QVAR on `k>=4`,
   `E|Z_psi|^2 >= 3q(q-1)/16` (the `k=1,3` strata are proved); and
2. the principal sufficient room
   `||delta||^2 <= n(n+10)^2/[6(n-6)^2]`.

The identity `tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2` remains a live
way to attack target 2, but bounding `||M||_F^2` alone is no longer the single
floor target because the exceptional block must be handled by QVAR.  The
finite `p=5` floor is already checked.  See the PSL reduction note and the
older `NOTE_leftover1_variance_multiplicity_route.md` for the variance history.

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
| \(\lambda_{\min}(\Phi)\ge6\) on \(Z\) | `phi_F_ge_6_proved_general=False` | Open. Exceptional QVAR remains only on `k>=4`; principal blocks retain the delta-variance target. |
| Residual (ii), even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | Open. Affine and even \(k\le4p-2\) are closed (15.179, 15.236, 15.237). |
| Type I, Max− not two-level \(\{-1,-3\}\) | `type_I_multilevel_bad_case_ND_closed=False` | Open. Remainder is \(A_{\mathrm{full}}\). |
| Lemma D | True | Closed. Do not unflip. |

**Next attack.**  For the exceptional block, prove QVAR directly on the union
`k>=4`, or prove the equivalent odd-coset degree-four harmonic excess is at
least `-q(q-1)(q-11)/(16(q+5))`.  Do not use a pointwise/orbitwise floor,
restricted-stratum PSD, ordinary minimum-shell design, or “quartic is top on
every stratum”: each is now disproved.  Positivity, conserved total, cyclic
symmetry, full support, and coarse divisibility are also insufficient; a
profile proof must now exploit the cross-direction coefficient kernels and
simultaneous Boolean realizability, with active-direction configurations mixed
before taking the second moment.  In parallel, an upper bound on
`||M||_F^2` may close the principal delta room.  Import `phi_F_ge_6` only when
both block types are controlled generally, never from finite-p data.

## Floor (leftover 1)

**Current block decomposition.** `Z=W_e direct-sum principal series`; QVAR is
the exact exceptional condition and the delta room is the current sufficient
principal condition.  The formulas below are the older equivalent Fourier
description and remain useful, but “name all of Q(r)” is no longer the binding
next step.

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
| `evidence/NOTE_2026-08-20_psl_and_stratum_floor_reduction.md` | Binding PSL decomposition, QVAR, low-stratum theorem, odd-coset route, and killed routes |
| `evidence/quartic_profile_attack.py` | Direct quartic/profile-energy diagnostic |
| `evidence/maxplus_p11/directional_energy_covariance_p11.{py,json}` | Full p=11 directional covariance and top-mode counterexample |
| `evidence/maxplus_p11/k4_active_subset_quartic_p11.{py,json}` | Exact p=11 k=4 active-subset split: balanced families fail QVAR, aggregate clears |
| `evidence/NOTE_leftover1_variance_multiplicity_route.md` | Historical principal variance route and `||M||_F^2` target |
| `evidence/METHOD_why_500_props_never_moved_a_flag.md` | Read before writing a new numbered proposition |
| `evidence/maxplus_p11/` | Scripts + logs for the p=11 spectrum/moment computations |
| `LONG_HORIZON_GOAL.md` | Terminal states |

Large `.npy` arrays (Max+ at p=11, 4.5 GB each) are **not in git** — they live at
`/mnt/storage/e1work/maxplus_p11/` on soulkiller, verified by md5 against the
original computation. Scripts there have hardcoded `/tmp/e1work` paths; repoint
before rerunning.

**Checkpoint entering the latest continuation:** `54fd110` (exact `2p`
arithmetic, normalized integer QVAR, and the general coarse-profile
countermechanism).  The latest commit(s) after that checkpoint add the p=11
fixed-active-subset counterexample in Prop 15.589 J; use `git log -1` for HEAD.
