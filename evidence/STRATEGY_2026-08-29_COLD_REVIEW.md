# Cold review: the shape of the problem and the shortest live attack

**Date:** 2026-08-29
**Verdict:** \(L\) is still open, but the search space is materially narrower
and the next multi-agent pass has three bounded targets

## 1. The original problem, stripped back down

For a signing \(a\) of the edges of \(K_n\),

\[
Q_a(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
m_n=\min_a\|Q_a\|_{L^\infty(\{\pm1\}^n)}.
\]

Equivalently, \(m_n/2\) is the deficit of the covering radius of the
antipodal cut code. The target \(m_n\sim\frac12n^{3/2}\) is therefore a
global cube-covering statement, not primarily a finite-field classification
problem.

The exact amount needed for \(L=1/2\) is only

\[
\Phi(C_p)-m_{p^2+1}=o(p^3)
\tag{1}
\]

on any ratio-dense Paley tail. Proposition 15.20e makes (1) equivalent to
the global conclusion. The current four-unit, every-prime, exact gap-2 gate
is a much stronger sufficient architecture. It remains useful, but it is not
the theorem itself.

## 2. The geometric shape

There are two dual pictures.

### Cube versus conference eigenspaces

A conference matrix has the smallest possible operator norm. Paley at
order \(p^2+1\) also has Boolean \(\pm p\)-eigenvectors, so its cube maximum
hits the spherical ceiling exactly. A proof of \(L=1/2\) must show that no
other signing can reduce the cube maximum by a fixed proportion while paying
for its larger spectral norm.

A permanent-gap counterexample has a very rigid shape. Relative to Paley it
must differ on \(\Theta(n^2)\) edges and satisfy the global slab

\[
\frac{Q_C(x)-H}{2}\le S_F(x)\le\frac{Q_C(x)+H}{2}
\quad\text{for every Boolean }x.
\tag{2}
\]

Sparse flips, bounded-degree deformations, and local minima cannot do this.

### Signed even-Eulerian gas

For the augmented cut code, let

\[
P_a(t)=\sum_{H\text{ even Eulerian}}
\left(\prod_{e\in H}a_e\right)t^{|H|}.
\]

Then exactly

\[
\mathbb E_x\cosh(\beta Q_a(x))
=(\cosh\beta)^{\binom n2}P_a(\tanh\beta).
\]

The initially proposed critical-temperature estimate

\[
\inf_a\log P_a(\tanh(2/\sqrt n))\ge-o(n)
\tag{3}
\]

would prove \(m_n\ge(1/2-o(1))n^{3/2}\), but it is false. For every
\(0<\theta<1\), a fractional-moment argument produces a deterministic
signing with
\[
\log P_a(\tanh(c/\sqrt n))
\le-\left(\frac c2-\sqrt{\log2}\right)^2n+o(n).
\]
At \(c=2\) this is a negative linear upper bound, contradicting (3).

The corrected fixed-\(c\) sufficient target is
\[
\inf_a\log P_a(\tanh(c/\sqrt n))
\ge\left(\frac c2-\frac{c^2}{4}\right)n-o(n).
\tag{3'}
\]
The fractional-moment obstruction rules this out for
\(c<\log2/(\sqrt{\log2}-1/2)=2.0843108\ldots\). The clean surviving choice
is \(c=3\):
\[
\inf_a\log P_a(\tanh(3/\sqrt n))\ge-\frac34n-o(n).
\]

So the likely final solution is not a circle or a finite census. Its shape is
either a global signed free-energy inequality such as (3'), or a rigidity
theorem saying that the only signings able to suppress that free energy are
asymptotically conference-like in the relevant Boolean sense.

## 3. What the cold audit killed

1. **Local/product stability.** For every fixed \(C\), random signings
   followed by greedy edge descent produce edge-local minima with both cube
   norms \(O(n^{3/2})\) but signed Paley distance
   \(\binom n2/2-O(n^{3/2})\). Product/frame second moments are exactly
   independent of the signing. Local optimality cannot imply
   \(k_\star=o(n^2)\).
2. **Spectral defect identifies Paley.** False. Defect zero identifies all
   conference classes. At square order, \(\rho=1\) is equivalent to the
   switching class having a constant-row-sum representative; square
   \(n-1\) is necessary but is not known to be sufficient. Order 50 is the
   first unresolved case.
3. **Generic character cancellation.** A sharp Weil bound still misses the
   normalized Paley dilation target by a factor \(\Theta(p)\). The full
   Max+ multiplicity average has to remain intact.
4. **Quartic parity majorants improve the bulk floor.** The truncated-moment
   criterion gives \(M_4=1\) uniformly for \(p\ge43\) in the bulk, so the
   constant majorant is already optimal. This route cannot raise the existing
   \(2p\) floor there.

The precise proofs and replacement targets are in
`NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md`.

## 4. New positive structure

### Closest-global hierarchy

If \(A=C\oplus F\) is a global minimizer chosen closest to the signed Paley
orbit, then every nonempty \(H\subseteq F\) has a signed-cut witness \(z\)
with

\[
\sum_{e\in H}C_ez_e
\ge1+\frac{m_n-\langle A,z\rangle}{2}.
\tag{4}
\]

This all-subsets condition is the first stability invariant that genuinely
uses global/cardinality minimality. It survives the local no-go theorem.

### One character-energy lemma now controls two fronts

Let \(q=p^2\), \(H=(\mathbb F_q^\times)^2\),
\(K=H/\{\pm1\}\), and let \(L(a)\) be the full-Max+ dilation correlation.
After auditing the split-torus class sizes, the centered character energy is
exactly

\[
\frac{24\|\delta\|^2}{n}
=\frac1{q-1}\sum_{t\in H}|\gamma(t)|^2,
\]

with \(\gamma(1)=\gamma(-1)=0\). Here
\(\bar L=q(q-1)(q+5)/16\), and the unresolved estimate is

\[
\sum_{a\in K\setminus\{1\}}
\left[
\frac{16}{q}(L(a)-\bar L)+\frac{8(q-1)}{q-5}
\right]^2\le q-1.                                  \tag{5}
\]

If \(S_K\) denotes the left side, then exactly
\[
S_K=12(q-1)\frac{\|\delta\|^2}{n}.
\]
Thus (5) is the strong R1 target \(\|\delta\|^2\le n/12\) itself in
dilation coordinates, not a new downstream lemma. Representation theory,
positivity, equivariance, trace, fixed norms, and additive-autocorrelation
PSD do not imply it: explicit abstract spectra and PSD autocorrelations
violate the bound. A proof must use the Boolean rank-one identity
\(B_y+2P_+=yy^\top\) and the exact uniform Max+ orbit mixture. Pointwise
\(|\gamma|\le2\) is also false at \(p=5,7\).

### The residual front just shortened

Proposition 15.688 proves the sharp all-prime lift bound

\[
4p\mathbb E B\ge p-3
\]

for every nonzero nonnegative integral quadratic on the middle slice. It
deletes every positive-residue row at the live \(p=19,s=16\) endpoint. The
deficit-minimizing pair is

\[
5[b=0]+5[b=16],\qquad9[b=2]+[b=16],qquad
\text{putative pair slack }34.                     \tag{6}
\]

The pair (6) is impossible modulo four, but it is not the complete row.
Exact completion leaves 143 phase-labelled residue-zero profiles.
Proposition 15.689 excludes all 129 profiles of slack at most twelve,
leaving fourteen with histogram
\(\{16:7,20:4,24:1,28:1,32:1\}\).

## 5. The bounded Ultra attack

The next pass should use three proof teams and one red team, with no branch
allowed to widen its target.

1. **Character team:** prove or decisively obstruct the eventual dilation
   energy inequality (5). This is the highest leverage target because it
   joins principal R1 and Type I.
2. **Residual front:** attack only the fourteen \(p=19\) profiles left by
   15.689, then stop. If they close, move to the existing \(p=17\) row; do
   not generate later shells.
3. **Global team:** attack either the all-subsets hierarchy (4) or the
   signed-Eulerian free-energy inequality (3), with the acceptance condition
   being the actual asymptotic deficit (1), not exact gap two.
4. **Red team:** test counterexample geometry only: Mathon's ratio-dense
   conference family for a uniform \(\rho<1\) gap, and square-order
   nonregularizability beginning at order 50. A fixed gap wins; isolated
   nonsaturation does not.

That is the current shape of the problem and the shortest honest attack map.
It is ready for the Ultra multi-agent pass.
