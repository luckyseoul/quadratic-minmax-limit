# Session handoff: residual (ii) audit fix (2026-08-06)

## Cold verdict accepted

Auditor: main theorem still open; residual (ii) affine two-level branch closed;
**full residual (ii) not audit-complete** without exhaustiveness.

## What we did

1. **Audited exhaustiveness** from score inequalities:
   - Freeness-fail ⇒ \(f_e\equiv+1\) on \(U_2\) (definition) — proved.
   - Under freeness-fail + \(s_+=2\): Max+ weak-ND witness fails (\(S_H\ge3\)) — proved.
   - Freeness-fail **does not** force support \(S\subseteq\{2,4\}\) (multi-level mass examples exist).
   - Freeness-fail + two-level **does not** force \(f_e\equiv-1\) on \(U_4\) (non-affine escapes freeze).

2. **Shipped Prop 15.193** (`src/e1_gmin_m4_prop15193.py`):
   - `residual_ii_affine_branch_closed()` = True (via 15.179)
   - `residual_ii_exhaustiveness_proved()` = False
   - `residual_ii_full_closed()` = False
   - Open subcases (ii-a) multi-level, (ii-b) non-affine two-level

3. **Retracted full residual-(ii) CLOSED**:
   - `deep_s2_freeness_fail_k_ge_3p_ND_closed()` now wires to `residual_ii_full_closed()` → **False**
   - Affine pieces still checkable via `residual_ii_affine_branch_pieces_ok()` / 15.179

4. **Stale wording fixed**:
   - 15.168 OPEN block (was flipped: residual i CLOSED / ii open sketch)
   - 15.169.I, 15.171 header/predicates, 15.179.G
   - Hinge statuses 15.180–192 `residual_ii_closed: False`
   - STATUS.md, HANDOFF.md, denseness_path_package.md

## Predicates (honest)

| Predicate | Value |
|-----------|-------|
| `residual_ii_dual_twolevel_affine_closed` | True |
| `residual_ii_exhaustiveness_proved` | False |
| `deep_s2_freeness_fail_k_ge_3p_ND_closed` | False |
| `type_I_k_3p_minus_2_closed_general` | False |
| `gsum_disj_lb_proved_general` | False |
| `e1_closed_general` | False |
| `L` | OPEN |

## Recommended next math

**Exhaustiveness lemma** (desired): under residual-(ii) hypotheses, failure of weak ND implies \(S\in\{2,4\}\) and \(f_e=3-S\) on every Max+ vector.  
If proved line-by-line from score inequalities, flip `residual_ii_exhaustiveness_proved` and residual (ii) full closes via 15.179.  
If false, the multi-level / non-affine counter-configuration is the remaining residual-(ii) attack surface.

Residual (i) hinge unchanged: \(|\mu|\le2/n\) or ker-box empty for all \(p\ge5\).

## Tests

`pytest tests/test_prop15193.py tests/test_prop15171.py tests/test_prop15179.py tests/test_prop15169.py tests/test_prop15170.py tests/test_prop15180.py tests/test_gsum_hinge_honesty.py` green.
