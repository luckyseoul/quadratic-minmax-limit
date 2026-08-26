# Long-horizon goal (binding)

**Done means exactly one thing:** MathOverflow [413935](https://mathoverflow.net/questions/413935) is settled.

That is: prove that

\[
L=\lim_{n\to\infty}\alpha_n
\]

exists and determine its value, **or** prove that the limit does not exist.

A session, plan, prop, test suite, or handoff update is **not** completion of this goal.

## The only acceptable terminal states

| Terminal state | What must be true |
|----------------|-------------------|
| **\(L=1/2\) CLOSED** | E(1) on the Paley family \(n=p^2+1\) is proved Max+-free for all primes \(p\ge5\), denseness (Prop 6.1–6.2) is applied, and `e1_closed_general` is True via **real imports** from the hinge modules (not a handwritten `return True`). `solution.md` Main Theorem states \(L=1/2\). |
| **\(L=c\) CLOSED for a specific \(c\neq 1/2\)** | A proved value with the same wiring standard. |
| **Non-existence CLOSED** | Two dense subsequences with unequal proved \(\lim\alpha\). Denseness (Prop 6.2) is mandatory. |

Anything else — including “honest OPEN,” “structure shipped,” “census holds at \(p=5,7\),” “AI-test later,” “good increment” — is **not done**.

## What is already proved (do not re-derive)

- Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le 1/2\).
- Denseness: the limit exists globally iff it exists along Paley orders \(n_k=p_k^2+1\).
- \(\rho=1\) on that Paley family, so \(\Phi(C_n)=\frac12 n\sqrt{n-1}\) and \(\alpha_{n_k}\to 1/2\) **if** E(1) holds there.
- 15.167 majorization algebra: bi-tight empty **if** \(\lambda_{\min}(\Phi)\ge6\). The floor is **not** proved (\(G_{u,\mathrm{disj}}\) is not a Gram).
- Exact Paley optimality is **false** (\(m_{10}=13<15\)). E(1) is asymptotic.

Therefore the remaining theorem is **E(1)** on \(n=p^2+1\). The live
acceptance gate is the four-unit dictionary returned by
`src/e1_main_chain_status.py`; as of 2026-08-26 it reads:

1. **Spectral floor:** OPEN. The current route requires both global mixed-
   \(k\) QVAR and principal R1.
2. **Residual (ii), even \(k\ge4p\):** OPEN. The Walsh slice, Eulerian
   boundary, every two-point boundary, and every four-point boundary are
   closed; every six-point boundary is also closed for `p>=11`, and both
   `p=7` infinity-plus-five signs are closed by 15.658--15.659. The
   six-finite `p=7` branch, all `p=5` size-six cases, and boundaries of size
   at least eight remain.
3. **Type I, multi-level Max−:** OPEN. The general \(3A+B>0\) estimate is
   not proved.
4. **Lemma D:** TRUE.

The older statement that residual (ii) was closed by 15.179/15.236/15.237
covered an earlier, narrower split and is not the live multi-level predicate.
Likewise, “residual (i) is the only leftover” is retired shorthand. The
authoritative functions are `residual_ii_k_ge_4p_ND_closed()`,
`type_I_multilevel_bad_case_ND_closed()`, and
`phi_F_ge_6_proved_general()`.

Path-C / 16N / 15.193 exhaustiveness is **not** required.

Closing any one open unit is **necessary progress**. It is **not** goal
completion. Goal completion is a terminal state in the table above.

## Forbidden translations (these are the small-failure mode)

Do **not** replace this goal with any of:

- Ship a Max+-free identity / dead-path / layer rewrite whose `residual_i` / `e1` / \(L\) flags stay False, then treat the session as done.
- “Continue residual-(i)” as a plan whose acceptance criteria are “a new increment exists.”
- Recertifying that unsigned \(\sum|\mathrm{per}|\) exceeds \(B\) (15.231 already did this).
- Updating STATUS / HANDOFF / session notes as the deliverable.
- Census \(p\le 7\), SA, or Lipschitz scale-counting as E(1).
- Soft-close: sandwich + denseness + \(\rho=1\) \(\Rightarrow L=1/2\) without E(1).
- Flipping predicates without a general Max+-free hinge imported for real.

If a proof attempt fails, **stop**. Record the failed mechanism in one paragraph. Do not wrap the failure as Prop 15.xxx with `proved: True` on an identity that does not flip a leftover.

## Session rule

A session may end only if one of these happened:

- one of the three open E(1) units **actually closed** (predicate True via a
  real import), or
- a live route was **killed as a path** by a general counter-mechanism (not a
  small-\(p\) census), and the writeup names the replacement route, or
- A terminal state in the table was reached, or
- The human redirected the goal.

“We shipped structure and left \(L\) OPEN on purpose” is **not** an allowed end condition for this goal.

## After a terminal state

The writeup must stand on its own: independent cold reads of `evidence/share/denseness_path_package.md` should call the argument essentially correct. Channel is X + GitHub, not a MathOverflow answer.
