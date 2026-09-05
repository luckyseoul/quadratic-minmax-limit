# Global minimal-four-gap bridge: Type-I scope correction

Date: 2026-09-04.

Status: **scope audit and narrow retraction**, not a new numbered
proposition, a new support band, or a global closure. The proved
15.750 restricted Type-I theorem, all 15.774 shell-capacity and odd
residual results, and all 15.775 results remain intact.

## 1. What global cardinality minimum actually supplies

Let `C` be the Paley conference signing of order `n=p^2+1`, prime
`p>=5`, and `Phi_0=p*n/2`. If any signing has norm at most
`Phi_0-4`, choose `H` of minimum cardinality among **all** edge sets
with that property, and write `A=C triangle H`.

Then `Phi(A)=Phi_0-4`: a norm at most `Phi_0-6` would permit restoring
any one edge while staying below the four-gap threshold. Every
deletion has norm exactly `Phi_0-2`, by cardinality minimum, score
parity, and the one-edge Lipschitz bound. The minimum is global, not
merely an inclusion minimum. It additionally rules out any improving
deletion/insertion exchange that preserves the threshold.

For an odd `H`, a signed shell score three gives a deletion `G` with
even size, signed minimum two, deep two-sidedness, and the distinguished
edge frozen positive on every critical row. The phase normalization and
size lower bound of 15.764 give the complete official residual-(ii)
entry. This arrow is unaffected by the correction.

For an even `H`, a signed shell score two gives instead an odd `G`
with signed minimum one and the same critical-edge alignment. It gives
`|G|>=2p+1`, not an exact size. This is **general level-one shell
alignment**, not entry into the restricted proved Type-I unit.

## 2. The literal hypotheses of 15.750

The header and proof of
[`src/e1_gmin_m4_prop15750.py`](../src/e1_gmin_m4_prop15750.py)
assume, after phase normalization,

```text
|G|=3p-2,   e notin G,
S_G=3-2 f_e                         on every Max+ row,
S_G<=-1 and S_G<=-3 f_e             on every Max- row.
```

In particular the multiset `W=G+2e` has exactly `3p` edges counted
with multiplicity and is constant of level three on Max+. Both facts
are used in its isolated-chart proof. The all-prime closure Boolean
certifies this box, not every odd cardinality with shell minimum one.

The even-H ledger of 15.764 proves neither `|G|=3p-2` nor the full
Max+ affine identity. Calling that ledger "Type I" cannot supply
the missing hypotheses. Its corrected metadata therefore distinguishes
`shell_level_entry_proved=True` from `official_entry_proved=False`.

The counting argument does not repair the gap. Write
`a=N_1/N`, `theta=(p+1)/(2p)`, where `N_1` counts the Max+ rows with
`S_G=1`. Critical-edge freezing gives `a<=theta`; the mean bound gives

```text
|G|/p >= a+3(1-a) >= 3-2theta = 2-1/p.
```

This is a lower bound `|G|>=2p-1`, not an upper or equality bound.
Even adding `|G|=3p-2` and `a=theta` fixes only the mean of `S_G`
on the complementary rows to five. It does not force their scores
all to be five. For example, at the scalar level, scores one on
the `f_e=+1` fraction `theta`, and scores three and seven on equal
halves of the `f_e=-1` fraction, satisfy that mean, counting equality,
and `S_G+f_e>=2`, but not `S_G=3-2f_e`. This is an algebraic
illustration of the missing implication, **not a Paley graph example**.
The existing `{1,5}` support hypothesis cannot be silently dropped.

## 3. Exactly what is retracted and retained

The old 15.774 consequence claiming every even minimal four-gap `H`
has size at least `6p-10` for `p=29,31` or at least `6p+6` for
`p>=37` used the restricted Type-I Boolean as if it closed general
level-one alignment. That unconditional consequence is retracted.

The valid conclusion is conditional on **no signed level-two row**:
then both signed shell floors are four, and the unchanged r4 capacity
theorem gives exactly those size lower bounds. With contact allowed,
the basic frame and bi-tight-level-two exclusion still give
`|H|>=2p+2`; the stronger unconditional claim is not proved.

The source retains the old unconditional field as `None` with status
`RETRACTED_SCOPE_MISMATCH`, exposes the no-level-two bound under its
own explicit key, and separately retains the basic frame bound.
Literal cardinality and affine-identity guards prevent a true restricted
15.750 Boolean from being promoted to a general even-H conclusion.

No part of the local small-mass theorem, the r3/r4/r5 shell-capacity
proof, the odd minimal-H bound, the two residual layers of 15.774,
or the eventual layer and cubic band of 15.775 uses the retracted
general Type-I implication. Those theorems are unchanged.

## 4. The corrected single global acceptance gate

The global minimal-four-gap bridge still requires an argument valid for
all possible support sizes. For a selected minimum-cardinality witness,
an odd shell contact at three supplies the official residual entry.
For even support, merely forcing contact at two does not finish the
bridge: one must also exclude the resulting general odd-cardinality
level-one dangerous-edge class, or force the exact hypotheses of a
proved unit. This missing implication is explicitly part of
`minimal_gap4_shell_bridge_closed_general()`, which remains false.

There is a useful check against an overly weak acceptance condition.
For any hypothetical exact four-gap signing, take an active signed
Boolean state and switch it, with an overall sign change if necessary,
onto one fixed positive Boolean C-eigenvector `y`. The resulting signing
`A'` has the same norm and `q_A'(y)=Phi_0-4`. Its relative edge set
`H'` therefore satisfies `T_H'(y)=2` and has even cardinality. This
normalization need not preserve minimum distance to `C`. Thus arbitrary
even contact is readily obtained from any exact four-gap signing; it
is not by itself an independently closed class.

The live proof target must exploit global minimum cardinality, common
Boolean-state structure, or another genuinely stronger invariant. No
new census or next-layer calculation is licensed by this audit.

Source and regression scope:
[`e1_gmin_m4_prop15764.py`](../src/e1_gmin_m4_prop15764.py),
[`e1_gmin_m4_prop15774.py`](../src/e1_gmin_m4_prop15774.py),
[`test_prop15764.py`](../tests/test_prop15764.py), and
[`test_prop15774.py`](../tests/test_prop15774.py).
The changed JSON receipts were regenerated on soulkiller. The offloaded
scope-correction replay passed 310 focused tests with zero failures,
errors, or skips; it took 80.506 seconds. No controller proof computation
or test run was used. Final documentation checks, commands, and hashes
are in `global_bridge_scope_mesh_replay.json` and
`global_bridge_scope_regression.json`; older byte-specific receipts are
retained as historical records, not silently reasserted for changed files.
