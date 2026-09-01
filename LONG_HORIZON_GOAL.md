# Long-horizon goal (binding)

**Done means exactly one thing:** MathOverflow [413935](https://mathoverflow.net/questions/413935) is settled.

That is: prove that

\[
L=\lim_{n\to\infty}\alpha_n
\]

exists (whether or not the proof identifies its value), **or** prove that the
limit does not exist.  The original MathOverflow question asks for existence;
requiring the value as well is a strictly stronger, optional objective.

A session, plan, prop, test suite, or handoff update is **not** completion of this goal.

## The only acceptable terminal states

| Terminal state | What must be true |
|----------------|-------------------|
| **Existence CLOSED; value unidentified** | A proved convergence theorem applies to the actual sequence.  For example, both Dini-summable amplification rays in Proposition 6.3 are proved, so \(H(n)=m_n^{2/3}\) satisfies the two-ray criterion.  A reduction with an open hypothesis is not enough. |
| **\(L=1/2\) CLOSED** | E(1) on the Paley family \(n=p^2+1\) is proved Max+-free for all primes \(p\ge5\), denseness (Prop 6.1–6.2) is applied, and `e1_closed_general` is True via **real imports** from the hinge modules (not a handwritten `return True`). `solution.md` Main Theorem states \(L=1/2\). |
| **\(L=c\) CLOSED for a specific \(c\neq 1/2\)** | A proved value with the same wiring standard. |
| **Non-existence CLOSED** | A rigorous proof that \(\alpha_n\) does not converge.  Two ratio-dense subsequences with unequal proved limits are one sufficient mechanism, not an extra requirement on every proof. |

Anything else — including “honest OPEN,” “structure shipped,” “census holds at \(p=5,7\),” “AI-test later,” “good increment” — is **not done**.

## What is already proved (do not re-derive)

- Sandwich \(1/\pi\le\liminf\alpha_n\le\limsup\alpha_n\le 1/2\).
- Denseness: the limit exists globally iff it exists along Paley orders \(n_k=p_k^2+1\).
- Two-ray convergence reduction (Prop 6.3): Dini-summable normalized defects
  at multipliers 2 and 3 suffice. In particular,
  errors `O(n/(log n)^(1+epsilon))` in `H`, or power-saving errors, force
  convergence. Propositions 6.4--6.5 identify the exact four-state
  equal-endpoint diamond for the all-Hadamard multiplier-2 construction;
  the hereditary endpoint conditions are automatic. Proposition 6.6 proves
  the diamond outside the explicit Hamming-central/joint-energy residue
  (6.20), but does not close the ray.
- \(\rho=1\) on that Paley family, so \(\Phi(C_n)=\frac12 n\sqrt{n-1}\) and \(\alpha_{n_k}\to 1/2\) **if** E(1) holds there.
- 15.167's majorization algebra is conditionally valid, but its final bi-tight implication is **retracted** because `ker(G-(n/2)P1)` contains `ker G`. Proposition 15.720 instead excludes the required bi-tight levels 2 and 3 for every prime `p>=5` by a degree congruence. Its bi-tight level-4 corollary does not exclude one-sided tight level 4.
- Propositions 15.726--15.727 historically narrowed the first general
  residual-(ii) endpoint. Proposition 15.733 closes the former
  `p=31,R=10` case, and Proposition 15.734 supersedes that endpoint route by
  closing the whole `k=4p` shell for every boundary and every `p>=13`.
  Propositions 15.735--15.737 extend the first-three-shell close to every
  `p>=11`.
- Exact Paley optimality is **false** (\(m_{10}=13<15\)). E(1) is asymptotic.

There are now two honest routes.  The direct route to the original question
is the two-ray amplification theorem of Proposition 6.3 (or any other proof
of convergence/non-existence).  The value-specific route to \(L=1/2\) is
**E(1)** on \(n=p^2+1\).  On that optional Paley route, the live gate is the
corrected dictionary returned by
`src/e1_main_chain_status.py`; as of 2026-09-01 through Proposition 15.750 it
reads:

1. **Required bi-tight levels 2 and 3:** TRUE by 15.720. The spectral floor,
   global mixed-\(k\) QVAR, and principal R1 are no longer acceptance gates.
2. **Residual (ii), even \(k\ge4p\):** OPEN. Propositions 15.734--15.737
   close `k in {4p,4p+2,4p+4}` for every boundary and every prime `p>=11`.
   Propositions 15.738--15.742 close `p=13,k=58`, and Proposition 15.743
   closes `p=17,k=74`. At `p=13,k=60`, the complete residue sieve leaves
   `u in {0,3,4,6}`; Propositions 15.744--15.745 close `u=3,0`, and
   Propositions 15.746--15.749 close `u=4`, leaving exactly `u=6`. The live
   remainder includes critical `p=5,7`, `p=11,k>=50`, that sole
   `p=13,k=60` residue and later p13 layers,
   every `p>=17,t>=4` layer (beginning with `p=17,k>=76`), generic branch B
   at `p>=29,t=3`, and the positive `p=7,z=7` subbranch. Historical
   endpoint/profile artifacts are not live gates; only the lemmas and
   certificates explicitly retained by the dedup audit remain valid evidence.
   Read `evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md` before reopening one.
3. **Type I, multi-level Max−:** TRUE by Proposition 15.750. Isolated-chart
   square rigidity and parity halving close `p>=11`; tracked exact integer
   Farkas identities close `p=5,7`. The old `3A+B` route remains incomplete
   but is no longer a global gate.
4. **Lemma D:** TRUE.

The older statement that residual (ii) was closed by 15.179/15.236/15.237
covered an earlier, narrower split and is not the live multi-level predicate.
Likewise, “residual (i) is the only leftover” is retired shorthand. The
authoritative open function is `residual_ii_k_ge_4p_ND_closed()`;
`type_I_multilevel_bad_case_ND_closed()` is now True. The valid bi-tight predicate is
`required_bitight_levels_empty_all_primes()`.

Path-C / 16N / 15.193 exhaustiveness is **not** required.

Closing the sole open residual-(ii) unit completes the remaining E1 ledger,
but it is not necessary for a different proof of convergence and is not goal
completion. The direct convergence front is to close the multiplier-two
diamond on the exact residue (6.20), then prove multiplier three or the
`1:2` split in Proposition 6.3. Endpoint selection, independent skew-norm
budgets, and an all-pairs composition theorem are not the live target. The
disk surrogate is strictly stronger, not equivalent: its asymptotic form
would prove a new `1/sqrt(2pi)` lower bound, and only its zero-error form is
currently disproved.

## Forbidden translations (these are the small-failure mode)

Do **not** replace this goal with any of:

- Ship a Max+-free identity / dead-path / layer rewrite whose `residual_i` / `e1` / \(L\) flags stay False, then treat the session as done.
- “Continue residual-(i)” as a plan whose acceptance criteria are “a new increment exists.”
- Recertifying that unsigned \(\sum|\mathrm{per}|\) exceeds \(B\) (15.231 already did this).
- Updating STATUS / HANDOFF / session notes as the deliverable.
- Census \(p\le 7\), SA, or Lipschitz scale-counting as E(1).
- Treating the all-prime, gap-2 Paley architecture as the only way to answer
  the original existence question.
- Soft-close: sandwich + denseness + \(\rho=1\) \(\Rightarrow L=1/2\) without E(1).
- Flipping predicates without a general Max+-free hinge imported for real.

If a proof attempt fails, **stop**. Record the failed mechanism in one paragraph. Do not wrap the failure as Prop 15.xxx with `proved: True` on an identity that does not flip a leftover.

## Session rule

A session may end only if one of these happened:

- one of the two fixed amplification rays in Proposition 6.3 is actually
  proved with a Dini-summable error, and the other ray remains the named direct
  convergence gate, or
- the sole remaining open E(1) unit **actually closed** (predicate True via a
  real import), or
- a live route was **killed as a path** by a general counter-mechanism (not a
  small-\(p\) census), and the writeup names the replacement route, or
- A terminal state in the table was reached, or
- The human redirected the goal.

“We shipped structure and left \(L\) OPEN on purpose” is **not** an allowed end condition for this goal.

## After a terminal state

The writeup must stand on its own.  If closure uses the Paley route, an
independent cold read of `evidence/share/denseness_path_package.md` must call
that argument essentially correct; a different convergence proof needs its
own self-contained package. Channel is X + GitHub, not a MathOverflow answer.
