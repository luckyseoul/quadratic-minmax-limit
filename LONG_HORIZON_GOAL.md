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
- Bi-tight empty for all primes \(p\ge5\) (15.167).
- Exact Paley optimality is **false** (\(m_{10}=13<15\)). E(1) is asymptotic.

Therefore the remaining theorem is **E(1)** on \(n=p^2+1\). Residual **(ii)** is **CLOSED** (affine 15.179 + (ii-b) 15.236 + (ii-a) 15.237). \(\nu=0\) on every \(|\kappa|=1\) four-set is **proved** (15.268), so \(m_4^+=m_4^-=\mu=\mu_{\mathrm{part}}+2\delta_+\) there. The only leftover:

1. Residual **(i)** — Type I freeness-fail \(k=3p-2\): control the **even** \(\delta\in E_{\pm4p}\) so \(|\mu_4|\le 1/(2p)\) on every \(|\kappa|=1\) four-set (or envelope / reflection / \(\|m_4\|_2^2\le n(n-2)/4\) / \(K_4\le\mathrm{Wick}_{hi}\) / \(\ker=\mathrm{sc}\) + free-\(e\), or Gsum\(\ge-1/p\), dual-eq empty). Particular majorant is already under the threshold. Local 4–5 point Grams only give \(|\mu|\le1-2/p\) (dead).

Path-C / 16N / 15.193 exhaustiveness is **not** required.

Closing residual (i) is **necessary progress**. It is **not** goal completion. Goal completion is a terminal state in the table above.

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

- Residual (i) **actually closed** (predicate True via real import), or
- Residual (i) was **killed as a path** by a general counter-mechanism (not a small-\(p\) census), and the writeup names the replacement leftover, or
- A terminal state in the table was reached, or
- The human redirected the goal.

“We shipped structure and left \(L\) OPEN on purpose” is **not** an allowed end condition for this goal.

## After a terminal state

The writeup must stand on its own: independent cold reads of `evidence/share/denseness_path_package.md` should call the argument essentially correct. Channel is X + GitHub, not a MathOverflow answer.
