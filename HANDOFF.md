# Handoff: min-max ±1 quadratic form

**Date:** 2026-08-22 (15.597 Theorem A* proved; **no leftover flag flipped**)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit
**HEAD:** `e4e359b` on `main`. Working brain is ALWAYS main.

**No leftover flag flipped.** Leftover 1/2/3 False. L OPEN. Aut-Schur /
Gsum / pairing False. `e1_closed_general` True only by the old incomplete
wiring. p=13 orbits / mesh k=6 are not a close.

## 15.597 Theorem A* (new, proved — not a census)

`src/e1_gmin_m4_prop15597.py`. Distinct from 15.108's old "Theorem A*"
(16N). For every prime \(p\ge5\), on \(Z\):

\[
\Phi_{\mathrm{part}}=\bar\lambda\,I,\qquad
\bar\lambda=8(n-2)/(n-6).
\]

Closed-form contractions (lemmas L1–L3 on \(W\in Z\)): \(\sum\kappa t=(n+1)\|W\|^2/4\),
\(\sum\phi t=-n\|W\|^2/4\), \(\sum\star t=-p\|W\|^2\). Direct quadratic-form
checks at \(p=5,7,11,13\) match. The particular solution is spectrally
**invisible**: all deviation of \(\Phi\) is \(\Phi_\delta\).

**Sharpened leftover 1** (equivalence, inequality still open):

\[
\lambda_{\min}(\Phi)\ge6
\iff
\Phi_\delta\succeq -\frac{2n+20}{n-6}\,I
\quad\text{on }Z.
\]

The \(n/12\) R1 bound is the multiplicity-floor op-norm form of the same
statement. Global QVAR is the exceptional isotype of \(\Phi_\delta\);
principal room is the rest.

**Corollaries (proved, buy nothing toward 6):** \(\mathrm{tr}(\Phi_\delta)=0\),
so \(\lambda_{\min}(\Phi_\delta)\le0\). \(\Phi\) is Gram, so \(\Phi\succeq0\)
and \(0\le\lambda_{\min}(\Phi)\le\bar\lambda\). Target 6 sits strictly inside
that window. The remaining content of leftover 1 is pushing the proven
lower bound from 0 up to 6; nothing short of a genuine bound on \(\delta\)
crosses it. Do not re-derive Gram \(\ge0\) or tracelessness as a floor
argument.

Do **not** add another equivalent identity with `inequality_proved=False`.
A* removes \(m_4^{\mathrm{part}}\) from the spectral problem; it does not
bound \(\delta\in\ker(4pI-T)\).

**Killed (2026-08-22):** Aut-invariant 4-point master equation
\(Tm=4pm-4\kappa/p\) plus \(|m_4|\le1\) cannot prove QVAR.
`src/e1_gmin_qvar_box_master.py`: permutation Aut quotient of \(T\) is
exact (\(T\kappa=-6\star\)); at \(p=5\), \(\dim E_{4p}^{\mathrm{Aut}}=2\)
and \(\min\langle m,\kappa_{A_\psi}\rangle=-285/4<0\), while true Max+
pairing is \(+14.13\). Same kill at \(p=7\) (ker 7, LP min \(\approx-2708\)).
Need a constraint outside linear 4-point theory (15.589 I).
Checked and redundant/useless for recovering the sign: \((1^\top y)^4=(1+p)^4\)
is already in the master affine space; \(P_{E_{4p}}\kappa_{A_\psi}\) is
indefinite (both signs at \(p=5\)); \(\ell^1\) CS overshoots (\(-172\) vs need
\(-6.75\)). Next constraint is simultaneous Boolean ridge / 6-point coupling,
not another 4-point rewrite.
Linear 6-point contractions + box still miss: p=5 min \(-101/4\); p=7
15.590 joint deg-6 ker 4, box min \(-10633/8\). SOS-4 along the deg-4
kernel is feasible at pairing \(-45/4\). True pairings positive. Need
strictly stronger than SOS-4 / linear 6-point (Boolean support or SOS-6).

## What shipped (wiring, not a close)

Leftover-1 conjunct is now **GLOBAL QVAR**, not per-stratum k≥7:

| File | Role |
|---|---|
| `src/e1_gmin_global_qvar.py` | Mixed-k floor identities A–I, P. `global_qvar_proved_general()` is A.proved and (P or G or H or I `inequality_proved`). All inequalities False. No handwritten True. |
| `src/e1_gmin_leftover1_qvar_principal.py` | `leftover1 = global_qvar AND principal_delta_room`. Does **not** AND `qvar_k_ge_7`. |
| `src/e1_gmin_qvar_k_ge_7.py` | Per-stratum k≥7 is **false** at (41,7) E=0 Cy=py and (13,7) pointwise. Flag False. Not a leftover-1 import. |
| `src/e1_gmin_r1_principal_pge11.py` | R1 L² recorded, **not** proved. `r1_l2_bound_for_p_ge_11` False. Interpolant 4/(p−3)² killed. |
| `src/e1_gmin_m4_prop15278.py` | `phi_F_ge_6` imports leftover1 AND. |

Gating tests: `tests/test_global_qvar.py`, `tests/test_leftover1_qvar_principal.py`,
`tests/test_qvar_k_ge_7.py`, `tests/test_r1_principal_pge11.py`.

## Live obstruction

Prove \(\Phi_\delta\succeq -(2n+20)/(n-6)\,I\) on \(Z\) (15.597). Split as:

- exceptional isotype = GLOBAL QVAR: \(\hat F(\psi)\ge 0\) / Gauss 4-distinct
  pairing of \(m_4\) on the full Max+ mixture (all \(k\), \(\lambda=0\) unsplit);
- rest = R1 / \(\|\delta\|^2\le n/12\).

Import `inequality_proved` only if that sign is actually proved.

Do **not** add another equivalent identity with `inequality_proved=False`.
Do not set leftover 1 True until GLOBAL QVAR **and** R1. Do not set
`qvar_k_ge_7` True. Leftover 2/3 / L stay their own units.

Killed as proofs: per-stratum k≥7, CS on \(\langle\delta,\kappa_A\rangle\),
two-level occupancy as a p-law (fails p=11 k=4/5), B-weighted 15.588
tautology, energy-only (15.589 I), Delsarte, cyclotomy \(E[E_0 E_r]\),
2-point fit-as-proof, pointwise SOS (Z=0 attained).

Scratch Aut/Torb probes on nuka (`/tmp/qvar-nuka/`,
`/tmp/grok-goal-f38dc225339a/implementer/probe_nuka_*.py`) are **not
shipped**. Torb vs \(T|_V\) still undiagnosed; do not import as a theorem.

## Mesh (out of QVAR scope)

PR #4 merged `mesh/k6-p13-enum` **into main** (`99cbf09`), so
`scripts/maxplus_profile_enum/` now lives on main. The live enum campaign
is still not a leftover close. cpu44 was hard-closed (2072 `orb*.npy`
kept; stop flag still set). Do **not** restart cpu44 unless the user
names cpu44 in the same turn. Worktree
`/home/nick/quadratic-minmax-limit-k6-mesh` may lag main.

## Compute

nuka (5700X3D) for serial / vcache. lucky is DNS only. Never 86 workers
on Orin. Soft-stop mesh unless the user hard-closes a named node.

## Suggested skills

`agent-cost-optimization`, `graph-engineered-completion`,
`use-available-compute`, `claude-referee` then `openai-referee`,
`verification-before-completion`, `handoff`, `scientific-critique`,
`grill-me`, `self-refine-loop`, `research`, `arxiv`.

---

**Date:** 2026-08-21 (R1 L² unit recorded, **not** proved; **no flag flipped**)

`tr(Phi^2) = 4||M||_F^2 - 3n^2 + 2n^2(n-1)/p^2` is identity (I) in
`TECHNICAL_NOTES.md` §4 / `METHOD.md`: leftovers 1 and 3 are moments of
one four-point tensor. It is proved (index split + E[(y·z)^2]=2n). It
does **not** close leftover 1: CS on dim Z, tr Φ, tr(Φ²) alone cannot
get λ_min≥6, and the exceptional block still needs GLOBAL QVAR (mixed-k;
per-stratum k≥7 is false).

Leftover 3 next sufficient target, not imported: p=5 is already a finite
from-C theorem (`type_I_p5_through_e_3AB_positive`). For every prime
p≥7, `|μ|≤2/n` on |κ|=1 is strictly stronger than `|μ|≤L=(p-2)/(2p^2)`
and would close Type I (`2/n < L` iff p≥7). Census slack at p=7 is
tight: 109/2863 vs 2/50 (~5%). Do **not** use `|μ|≤|f4|` (false at p=7,
15.191 I). `|μ|≤maj` remains false at p=7 and p=11. Triangle
`|R̄₄|+2|φ|` is too weak for `|μ|≤L` even at p=5. L2 conversion of
3A+B is rejected. p=5 μ=f4 on each (κ,φ) class; p=7,11 split inside
(κ,φ).

Leftover 2: leftover+splus empty all nF at p=5 k=20 (15.528);
leftover-only is not residual (ii). Walsh cannot flip leftover 2
(interior 4-level only). Uniform Paley E_-[S²]<20+12/p is false.

R1 (`‖P_{E_{4p}} m₄⁺‖² ≤ n/12`) is the binding leftover-1/3 face
(15.595). `src/e1_gmin_r1_principal_pge11.py` records exact measured
‖δ‖² vs R1 (exceeds at p=5,7; census-holds at p=11). p=5 measured
equals κ_hyp_δ (1536/65). The interpolant κ_hyp_δ·4/(p−3)² is **killed
as a retained δ-bound** (equality law false at p=7,11; Aut-dim
ν_G-ratio dies at p=7; no operator identity). `r1_l2_bound_for_p_ge_11`
is False. `principal_delta_room_moment_proved` imports that unit only.

`src/e1_gmin_leftover1_qvar_principal.py` is the Max+-free leftover-1
import. `phi_F_ge_6_proved_general` is `leftover1_qvar_and_principal_proved()`,
which is True only if **GLOBAL** QVAR (mixed-k) **and** the principal
δ-moment both hold for every prime p≥5. Both estimates stay False.
Identities A–D (QVAR iff, V_sph>threshold, D+ room formula,
⟨m4,κ_B⟩≥0 iff floor 6) have fail-eqs and tests. p=13 orbits are not
imported. Aut-Schur / Gsum / pairing False. L OPEN. Leftovers 2 and 3
False. Live `e1` is still the old AND.

---

**Date:** 2026-08-20 (current branch; use `git log -1` for the exact checkpoint)
**Repo:** https://github.com/luckyseoul/quadratic-minmax-limit
**Current branch:** `codex/leftover-moment-attack` (use `git log -1` for the
current hash).  The latest continuation adds Prop 15.589 Theorems J--U,
including the all-prime `k=4,5,6` QVAR closures described below.  All 124
focused Prop 15.588/15.589 and k=7 probe tests pass.
**Statement:** [MathOverflow 413935](https://mathoverflow.net/questions/413935).
\(\alpha_n=n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\lvert\sum_{i<j}a_{ij}x_ix_j\rvert\).

## Binding status: exact PSL and high-stratum reduction

Work on `codex/leftover-moment-attack` gives the multiplicity-free
decomposition

`Z = W_e direct-sum ((p^2-9)/8 distinct degree-(p^2+1) principal series)`,

where `dim W_e=(p^2+1)/2`. Thus every Phi eigenvalue has multiplicity at
least `n=p^2+1` except one exceptional scalar. That scalar is exactly

`lambda_exc = 32 E|Z_psi|^2/[q(q-1)]`, `psi^2=chi`,

so its floor is the single quartic variance inequality
`E|Z_psi|^2 >= 3q(q-1)/16`. The remaining principal floor route is the
exception-removed delta room
`||delta||^2 <= n(n+10)^2/[6(n-14)(n-6)]`.  This is sharper than the old
room by the factor `(n-6)/(n-14)` and is equivalent to
`E[(y dot z)^4] <= 4n(3n^2-37n+2)/(n-14)`.  It applies after QVAR proves
`lambda_exc>=6`. Neither open inequality is proved generally; no flag is
flipped. See
`evidence/NOTE_2026-08-20_psl_and_stratum_floor_reduction.md` and Prop 15.589.

Two routes are now explicitly dead: restricted Phi does not have floor 6 on
every profile stratum, and at p=7 a full PSL orbit of size 1,176 has
`Z_psi=0`, so the quartic bound is not pointwise/orbitwise.

The pointwise repair also fails inside the first live stratum itself.  An
exact coupled coefficient/profile/Boolean CP-SAT model at `p=13,k=7` found a
translation-gauged Max+ representative with `Z_psi=-28-42i`, hence
`|Z_psi|^2=2548 < 10647/2`.  The stored support is independently checked
against `Cy=13y`, all seven active square directions, and the direct quartic
kernel.  Thus QVAR on `k>=7` must use the ensemble mixture even at its first
unresolved prime.  See `evidence/k7_p13_cpsat_{probe.py,witness.json}`.

Orbit averaging makes real progress on that counterexample.  Its signed-PSL
lift orbit is free and has size `4,826,640`, split evenly between the two
global signs.  Jellyfin's A380 evaluated all `2,413,320` epsilon-plus vectors:

```
E_orbit |Z_psi|^2 = 806468/85 > 10647/2,
lambda_exc(orbit) = 19088/1785 > 6.
```

The orbit contains `14,196`, `28,392`, and `2,370,732` vectors at activities
`k=5,6,7`.  Crucially, its `k=7` slice alone clears QVAR with exact moment
`1606124/167`.  Independent profile reconstruction finds 12,852 depressed
nonzero-quintic representatives, exactly 1,071 in each of the twelve scalar
classes.  This does **not** close `p=13,k=7`: after forbidding all 1,071
scalar-7 representatives, independent eight-worker CP-SAT runs on Nuka and
Jellyfin found a second-orbit seed with `Z_psi=-132-198i` and
`|Z_psi|^2=56628`.  That seed also has a free 4,826,640-vector signed orbit.
Its full epsilon-plus orbit has the same exact mean `806468/85`; its
2,384,928-vector `k=7` slice has moment `198692/21`, again above QVAR.  The
next finite task is therefore iterative orbit decomposition, not pointwise
minimization.  Both packed orbits and the second census are in the platter
backup named below.

The exceptional target has since narrowed again.  For `p=3 mod 4`, `Z_psi`
is exactly a signed sum of nonnegative directional profile energies whose
pointwise total is `p(p^2-1)/4`.  Combining that identity with the affine
profile classification, and using the Euler-product lower bound for
`L(2,chi_p)` when `p=1 mod 4`, proves `(QVAR)` on every `k=1` and `k=3`
stratum for every prime.  Theorems L--O close `k=4`, and Theorem P makes
`k=5` empty for every `p>=41`.  Exact finite sieves close every remaining
`k=5` case as well.  Exact quartic-profile energies make `k=6` empty for every
`p>=47`; coefficient sieves close every remaining finite case, while its
`p=11` case already clears by complete census.  Thus the exceptional scalar
is closed through `k=6` for every prime and remains only on `k>=7` from
`p=13`; the principal delta-variance target is sharpened as above.  See the
updated reduction note and Prop 15.589.

The high-activity coefficient attack now has a general normal form.  Whenever
the top degree `k-2` is nonzero, translation uniquely kills the full
two-dimensional degree-`k-3` level.  Every lower degree `d` has
`k-d-1` free coefficients, with any `d+1` direction coefficients serving as
invertible pivots.  For `k=7` this gives depressed quintics
`a*s^5+c*s^3+d*s^2+e*s+f`.  Exact probes at `p=13,17,19` find minimum profile
energies `1,3,4`; energy alone does not eliminate the stratum, while all
`1,36,120` direction subsets have the predicted kernel ladder.  The next
exceptional computation should therefore be a recursive coupled coefficient
sieve, not a Cartesian product (`evidence/k7_quintic_profile_probe.py`).

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
empty for all primes `p>=41` (Prop 15.589 L).

The same centered-Fourier argument plus Weil's additive-character estimate
works at arbitrary profile degree: a `k>=4` stratum is empty whenever
`p>4k^2`.  Hence every surviving stratum is `k=1`, `k=3`, or
`k>=sqrt(p)/2`; the first two already satisfy QVAR.  The unresolved
exceptional scalar is therefore asymptotically a high-activity problem
(Prop 15.589 M).

For `p=19,23,31`, the exact energy partitions and degree-2/degree-1
coefficient kernels leave zero constant-compatible candidates on all
`210,495,1820` direction subsets, respectively.  Together with the `p>=41`
barrier, this proves that when `p=3 mod 4`, `k=4` exists only at `p=7,11`;
its QVAR moment clears the target at both.  Thus the `k=4` contribution is
fully closed in this congruence class, and for `p>=19` QVAR starts at `k=5`
(Prop 15.589 N).

In the complementary congruence class, the same sieve proves `k=4` empty at
`p=29,37`; `p>=41` was already eliminated.  It regenerates the complete
`p=13,17` families and gives exact quartic moments `8788` and `314432/3`, both
above QVAR.  Since `p=5` has only three square directions, QVAR is now proved
on `k=4` for every prime.  The exceptional target starts at `k=5`
(Prop 15.589 O).

For `k=5`, a zero cubic-kernel scalar is impossible by the degree-at-most-two
energy bounds.  With nonzero scalar, all five profiles are cubic; translating
the input reduces exact energy enumeration to depressed cubics
`a s^3+c s+d`.  Their minima eliminate every prime `41<=p<101`, including the
sole numerical exception `p=43`, where all 28 relevant types have `b=45` and
sum to `225`, not the required `T=231`.  The general activity barrier handles
`p>=101`.  Hence `k=5` is empty for every `p>=41` (Prop 15.589 P;
`evidence/k5_cubic_energy_barrier.{py,json}`).

The finite cubic coefficient sieves close four more primes.  At `p=29`, all
736,828,092 low-energy type tuples give zero coefficient candidates; at
`p=37`, all 9,348 admissible leading patterns fail.  The `p=31` stratum is
nonempty but has only 8,000 translation representatives and exact moment
`E B^2=16704/5>45`.  The existing complete `p=11` census gives
`E B^2=163/9>45/8`.  Consequently `k=5` remains open only at
`p=13,17,19,23` (Prop 15.589 Q;
`evidence/k5_p{29,31,37}_coefficient_sieve.{py,json}`).

The same complete sieve closes those last four cases.  At `p=13,17`, direct
Gaussian-integer quartic evaluation gives moments `297468/31` and
`1650768/29`, above thresholds `10647/2` and `15606`.  At `p=19,23`, the
signed-energy moments are `29417/65>135/8` and `8908/19>99/4`.  Consequently
QVAR is proved on `k=5` for every prime (Prop 15.589 R;
`evidence/k5_p{13,17,19,23}_coefficient_sieve.{py,json}`).

For `k=6`, translating the input depresses every genuine quartic profile to
`a s^4+c s^2+d s+e`.  Exact minima at every prime `47<=p<=139` use more than
one sixth of the conserved energy; the general `p>4k^2` theorem handles
`p>=149`.  A vanishing quartic scalar is separately impossible by the cubic
and quadratic energy bounds.  Hence `k=6` is empty for every `p>=47`
(Prop 15.589 S; `evidence/k6_quartic_energy_probe.py` and
`evidence/k6_quartic_energy_probe_{low,high}.json`).

Exact coefficient-kernel sieves close three further `k=6` primes.  At
`p=37`, 8,189,942,400 raw type tuples leave no candidate on any of 27,132
direction subsets.  The same conclusion holds on all 54,264 subsets at
`p=41` and all 74,613 subsets at `p=43`.  Consequently the finite `k=6`
problem is reduced to `p=13,17,19,23,29,31` (Prop 15.589 T;
`evidence/k6_p{37,41,43}_coefficient_sieve.{py,json}`).

The accelerated exact elimination then closes all six residual cases at
once.  The `p=13,17,19` aggregate quartic moments are respectively
`8896212/955`, `149941632/2879`, and `10591740/103`, all above QVAR.  The
`p=23,29,31` coefficient sieves are empty after scanning 71.207, 20.937, and
2.971 trillion raw tuples.  Together with the complete `p=11` census and the
high-prime energy barrier, QVAR now holds on `k=6` for every prime
(Prop 15.589 U; `evidence/k6_p{13,17,19,23,29,31}_coefficient_sieve.json`).

**Settled.** Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le1/2\) (`solution.md`). Paley \(\rho=1\) on \(n=p^2+1\) (`evidence/PROOF_rho_eq_1.md`).

**Not settled.** \(L=\lim\alpha_n=1/2\). That needs E(1): Paley \(C\) is a \(\Phi\)-minimizer for every prime \(p\ge5\). Soft-close from sandwich plus denseness, without E(1), is not a proof.

Live `e1_closed_general` is True only by the old wiring (affine residual (ii) plus two-level Type I). That is not E(1) for all \(p\ge5\). Aut-Schur, Gsum disj LB, and the cotangent pairing stay unused and False.

## 2026-08-21 — remaining estimates isolated; no flag flipped

The p=13 k=7 signed-PSL census is not a leftover close. Three Max+-free
estimates remain; none is proved. See
`evidence/NOTE_2026-08-21_remaining_general_p_estimates.md`.

| Leftover | Predicate | Remaining estimate |
|---|---|---|
| 1 | `phi_F_ge_6_proved_general` | QVAR on k≥7 (all p≥13) **and** principal `\|\|δ\|\|^2 ≤ n(n+10)^2/[6(n-14)(n-6)]`. Crude `E[s^4]≤2n^3` is too weak. |
| 2 | `residual_ii_k_eq_4p_empty` / `multilevel_ND_k_ge_4p_proved` | leftover+splus at k=4p: Paley majorant `E_-[S^2]<20+12/p`, or Walsh pair-slice span (interior 4-level only). Leftover-only with min_+<2 exists and is not residual (ii). |
| 3 | `type_I_multilevel_bad_case_ND_closed` | `\|μ\|≤(p-2)/(2p^2)` on \|κ\|=1, equivalently `\|R̄₄\| ≤ \|L\|(p^4-1)+4(p-2)`. `\|μ\|≤\|T\|` does not close; `\|μ\|≤maj` is false at p=7. |

Live dump (twice): all three leftovers False; Gsum False; pairing False.
`e1_closed_general` is still True only by the old incomplete wiring.
Aut-Schur unused. L=1/2 OPEN. Do not import a leftover from p=13 orbits.

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

1. exceptional QVAR on `k>=7` from `p=13`,
   `E|Z_psi|^2 >= 3q(q-1)/16`; and
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
| \(\lambda_{\min}(\Phi)\ge6\) on \(Z\) | `phi_F_ge_6_proved_general=False` | Open. Exceptional QVAR is closed through `k=6` and remains on `k>=7` from `p=13`; principal blocks retain the delta-variance target. |
| Residual (ii), even \(k\ge4p\) | `residual_ii_k_eq_4p_empty=False` | Open. Affine and even \(k\le4p-2\) are closed (15.179, 15.236, 15.237). |
| Type I, Max− not two-level \(\{-1,-3\}\) | `type_I_multilevel_bad_case_ND_closed=False` | Open. Remainder is \(A_{\mathrm{full}}\). |
| Lemma D | True | Closed. Do not unflip. |

**Next attack.**  Complete the `p=13,k=7` orbit decomposition: accumulate the
depressed scalar-7 representatives from each new free orbit, ask CP-SAT for a
solution outside their union, and compute each orbit's exact `k=7` moment on
Jellyfin until infeasibility certifies exhaustion.  The first orbit contributes
1,071 scalar-7 representatives and its `k=7` slice clears QVAR; a second orbit
seed is already stored.  In parallel, prove QVAR generally on `k>=7`, or prove
the equivalent odd-coset degree-four harmonic excess is at
least `-q(q-1)(q-11)/(16(q+5))`.  Do not use a pointwise/orbitwise floor,
restricted-stratum PSD, ordinary minimum-shell design, or “quartic is top on
every stratum”: each is now disproved.  Positivity, conserved total, cyclic
symmetry, full support, and coarse divisibility are also insufficient; a
profile proof must now exploit the cross-direction coefficient kernels and
simultaneous Boolean realizability, with active-direction configurations mixed
before taking the second moment.  Quantization of exceptional projection norms
does not help by itself: the full p=11 census has 37 shells, no zero shell, and
minimum `4304/15 < 366=3n`, so even the nonzero-shell pointwise repair is
false.  The exact `p=13,k=7` witness `|Z_psi|^2=2548<10647/2` kills that
repair directly in the first unresolved stratum.  The coupled CP-SAT model is
now available for an orbit/weighted complete `k=7` second-moment computation;
more pointwise minimization cannot prove QVAR.  In parallel, an upper bound on
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
| `evidence/k4_p3mod4_coefficient_sieve.{py,json}` | Exact p=19,23,31 coefficient-level emptiness certificate |
| `evidence/k4_p1mod4_closure.{py,json}` | Exact p=13,17 quartic moments and p=29,37 emptiness; completes all-prime k=4 QVAR |
| `evidence/k6_coefficient_sieve_fast.py` | Exact Numba quadratic-three/linear-two elimination with orbit sharding |
| `evidence/merge_k6_coefficient_shards.py` | Validates and merges complete k=6 shard certificates |
| `evidence/k6_p{13,17,19,23,29,31}_coefficient_sieve.json` | Residual finite k=6 moments/emptiness; completes all-prime k=6 QVAR |
| `evidence/k7_quintic_profile_probe.py` | Exact depressed-quintic lift minima and universal seven-direction kernel audit |
| `evidence/k7_p13_cpsat_probe.py`, `k7_p13_cpsat_witness.json` | Exact coupled p=13,k=7 model and independently checked pointwise-QVAR counterexample |
| `evidence/k7_p13_signed_psl_orbit.py`, `k7_p13_signed_psl_orbit.json` | Packed-bit signed-PSL traversal: first witness has a free 4,826,640-vector signed orbit |
| `evidence/k7_p13_orbit_quartic_xpu.py`, `k7_p13_orbit_quartic_xpu.json` | Exact A380 quartic/activity census; first orbit and its k=7 slice both clear QVAR |
| `evidence/k7_p13_extract_orbit_representatives.py`, `k7_p13_orbit_completeness.py`, `k7_p13_second_orbit_seed.json` | Extracts 1,071 representatives per scalar and searches for the next orbit; a second seed was found and independently checked |
| `evidence/exceptional_projection_shell_probe.{py,json}` | Full p=5,7 exceptional shells and an explicit p=11 nonzero shell below `3n`; kills the nonzero-shell pointwise repair |
| `evidence/maxplus_p11/exceptional_projection_shells_p11_xpu.{py,json}` | Full 37.46M-row p=11 exceptional shell census: 37 shells, no zero, minimum `4304/15`, exact mean |
| `evidence/NOTE_leftover1_variance_multiplicity_route.md` | Historical principal variance route and `||M||_F^2` target |
| `evidence/METHOD_why_500_props_never_moved_a_flag.md` | Read before writing a new numbered proposition |
| `evidence/maxplus_p11/` | Scripts + logs for the p=11 spectrum/moment computations |
| `LONG_HORIZON_GOAL.md` | Terminal states |

Large `.npy` arrays (Max+ at p=11, 4.5 GB each) are **not in git** — they live at
`/mnt/storage/e1work/maxplus_p11/` on soulkiller, verified by md5 against the
original computation. Scripts there have hardcoded `/tmp/e1work` paths; repoint
before rerunning.

The p=13 orbit attack is centrally backed up at
`/mnt/storage/e1work/maxplus_p13/orbit_attack_2026-08-20/`.  It contains the
111 MB packed first and second orbits, the 1,071 packed scalar-7
representatives, both independent second-orbit CP-SAT certificates, both
orbit metadata records, both A380 censuses, and a SHA-256 manifest.  The
first-orbit packed hash is
`7223169420a18477dbdf95f6c3685186fbc6a7a1916ac875d761b22800c01eb2`.
The second is
`a3ce4e19e68770b41951b4ba28153fd5ed23884d0bcd912eeb43c421fa0e31c3`.

Jellyfin (`192.168.1.191`) now has a validated Intel Arc A380 environment at
`/home/nick/.venvs/mo-intel`: `torch 2.13.0+xpu`, `dpctl 0.22.1`, and
`pyopencl 2025.1`.  Torch XPU, Level Zero, and OpenCL all see the A380.  During
the shell census the soulkiller repo and p=11 data were mounted read-only under
`/home/nick/mnt/soulkiller-{repo,e1work}` via SSHFS.

**Checkpoint entering the latest continuation:** `54fd110` (exact `2p`
arithmetic, normalized integer QVAR, and the general coarse-profile
countermechanism).  The latest commit(s) after that checkpoint add the p=11
fixed-active-subset counterexample in Prop 15.589 J; use `git log -1` for HEAD.
