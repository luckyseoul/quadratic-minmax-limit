# Matching spike criterion: when perfect matchings cannot undercut

**Date:** 2026-07-27  
**Status:** Implication **proved**; n=10 complete classification certified; p=5 extensive census certified.  
**Existence of \(\lim\alpha_n\) remains OPEN** (needs criterion for all \(M\) when \(p\ge5\), plus control of non-matching undercutters / \(k_\star\)).

## Setup

Paley conference \(C\) of order \(n=p^2+1\), \(\Phi=\tfrac12 np\). Perfect matching \(M\), \(A=C\oplus M\), \(S_M(x)=\sum_{e\in M}C_ex_ix_j\). Note \(|M|=n/2\) is odd, so \(S_M\) is always odd (Prop 15.29).

## Proposition 15.30a (spike criterion ⇒ no undercut) — **proved**

If there exists \(x\in\{\pm1\}^n\) with either
\[
S_M(x)=-p\quad\text{and}\quad Q_C(x)\ge\Phi-2p,
\]
or
\[
S_M(x)=+p\quad\text{and}\quad Q_C(x)\le-(\Phi-2p),
\]
then \(\Phi(A)\ge\Phi(C)\).

*Proof.* In the first case \(Q_A(x)=Q_C(x)-2S_M(x)=Q_C(x)+2p\ge\Phi-2p+2p=\Phi\). The second case is symmetric. \(\square\)

## Single-bit drop (second level of \(Q_C\))

For any \(y\in\mathrm{Max}_{+}\) and coordinate \(i\), the one-bit flip \(y^{\oplus i}\) satisfies
\[
Q_C(y^{\oplus i})=\Phi-2p
\]
(formula \(Q(y^{\oplus F})=\Phi-2p|F|+4\sum_{i<j\in F}C_{ij}y_iy_j\) at \(|F|=1\)). Thus the second maximiser level is **at least** \(\Phi-2p\). Certified exact at \(p=3\): second max of \(|Q_C|\) on non-max vectors is exactly \(9=\Phi-2p\) (no values in \(\{11,13\}\)). At \(p=5\), sampling finds second max \(55=\Phi-2p\) (no values in \((55,65)\)).

## Multi-bit clique formula

For \(y\in\mathrm{Max}_{+}\) and \(F\subset[n]\),
\[
Q_C(y^{\oplus F})=\Phi-2p\,|F|+4\sum_{i<j\in F}C_{ij}y_iy_j.
\]
If \(F\) is a clique in the \(y\)-switched Paley graph (\(C_{ij}y_iy_j=+1\) all pairs) then
\[
Q_C(y^{\oplus F})=\Phi-2|F|(p-|F|+1).
\]
For \(|F|=p\) this equals \(\Phi-2p\). Combined with a matching-sign pattern giving \(S_M(y^{\oplus F})=-p\), the spike criterion fires.

Affine lines of **square** direction in \(\mathrm{AG}(2,p)\) are \(p\)-cliques of the Paley graph on \(\mathbb F_{p^2}\) (differences \((t-s)d\) have \(\chi\equiv\chi(d)=+1\)). This supplies the cliques for a constructive programme (not yet closed for every \(M\)).

## Certified classification

### \(p=3\) (\(n=10\)) — complete over all 945 perfect matchings

| \(\Phi(C\oplus M)\) | Count | Spike criterion |
|--------------------:|------:|:----------------|
| 13 (undercut) | 144 | **False** (max \(Q_C\) on \(S_M=-3\) is 7 \(<9\)) |
| 17 | 405 | True |
| 21 | 360 | True |
| 25 | 36 | True |

Spike criterion holds **if and only if** \(M\) does not undercut. Scripts: inline in tests + this note.

### \(p=5\) (\(n=26\)) — census

| Family | Spike criterion | Notes |
|--------|:----------------:|-------|
| 15 random PMs | **all True** | full enum of \(S_M=\pm p\) levels (\(\binom{13}{4}\cdot2^{13}\) vectors each); max \(Q_C\) on \(S=-p\) always \(\ge55\) |
| 3 SA Max-covers (two-sided) | **True** | exact MITM \(\Phi=65=\Phi(C)\) |

On \(S_M=-p\), identity \(\mathbb E[Q_C]=-p\) (matching contribution; cross terms average 0). Tail \(\max Q_C\ge\Phi-2p\) is the nontrivial content.

## What would finish matching non-undercut for \(p\ge5\)

Prove: for every perfect matching \(M\) and every odd prime \(p\ge5\),
\[
\max\bigl\{Q_C(x):S_M(x)=-p\bigr\}\;\ge\;\Phi-2p.
\]
Then by 15.30a, \(\Phi(C\oplus M)\ge\Phi(C)\) for all \(M\). Combined with \(p=3\) gap \(2=O(1)\) this yields E(1) **along the matching dichotomy**. Full E(1) still needs non-matching undercutters controlled (or shown absent for \(p\ge5\)).

## Not established

- Spike criterion for every \(M\) at general \(p\ge5\) (only certified samples)
- Matching dichotomy \(m_n=\min(\Phi,\min_M\Phi(C\oplus M))\)
- \(\lim\alpha_n\) exists

**Do not mark Main Theorem settled.**
