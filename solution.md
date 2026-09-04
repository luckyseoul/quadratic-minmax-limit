# Existence of \(\displaystyle\lim_{n\to\infty}n^{-3/2}\min_{a_{ij}=\pm1}\max_{x=\pm1}\bigl|\sum_{i<j}a_{ij}x_ix_j\bigr|\)

## Statement

For \(n\ge 2\) set
\[
m_n
=
\min_{\substack{a_{ij}=\pm1\\1\le i<j\le n}}
\max_{x\in\{\pm1\}^n}
\Biggl|
\sum_{1\le i<j\le n}a_{ij}\,x_i x_j
\Biggr|,
\qquad
\alpha_n
=
\frac{m_n}{n^{3/2}}.
\]

**Main Theorem (sandwich).**
\begin{equation}
\label{eq:sandwich}
\frac1\pi
\;\le\;
\liminf_{n\to\infty}\alpha_n
\;\le\;
\limsup_{n\to\infty}\alpha_n
\;\le\;
\frac12.
\end{equation}
Moreover \(\limsup_n\alpha_n=\lim_n a_n\) where \(a_n=\sup_{k\ge1}\alpha_{kn}\).
(The lower bound is the dual-Gaussian arcsine argument of Prop.~5.2, valid
for **every** Seidel matrix; Proposition 6.5e later refines its finite-order
right side to
\(n(n-1)\arcsin(1/\sqrt{n-1})/\pi\) plus a nonnegative square correction.
The classical Bohnenblust–Hille floor \(2^{-5/2}\) is retained as Prop.~5.1.)

**Direct convergence reduction — not yet a solution.**  The original problem
asks whether \(\alpha_n\) converges; it does not require identifying the
limit.  E(1) and both direct amplification rays are not complete.
Proposition 6.3 proves that convergence follows from just two
Dini-summable amplification estimates, at multipliers \(2\) and \(3\), for
\(H(n)=m_n^{2/3}\); polynomial saving is unnecessary.  Propositions
6.4--6.5 reduce Hadamard doubling to the equal-endpoint mixed-state diamond,
and Proposition 6.5a identifies it exactly with a simultaneous directed
half-cut neighbor problem for an oriented optimal signing. Proposition 6.5c
gives a distinct exact opposite-diagonal diamond with an arbitrary cross
signing; its four-label form isolates a simultaneous outgoing-half-cut
paving problem. Proposition 6.5d gives the exact nonlinear bivector cover
and proves the affine relaxation trivial, the displayed elliptope
subcritical, and a precisely normalized family of single-row even-moment
certificates blind through degree (o(n)).
Proposition 6.5e proves the outgoing-half constant is sharp by a
signed-regular arcsine bound with a positive square correction,
and Proposition 6.5f shields every fixed or low-signature anchor family,
while Proposition 6.5g constructs the resulting approximate mate near the
lower floor and proves its generic spectral bridge too weak. Proposition
6.5h gives the exact independent-random-orientation criterion for the
outgoing-half lift and proves that its exact bad-event first moment diverges
exponentially at the optimal-signing scale; this does not rule out random
orientations analyzed through event dependence. Proposition 6.5i shows that
arcsine saturation forces Hamming-central opposite near-extremizers and the
calculator's balanced saddle, while Proposition 6.5j rules out the full
degree-four preordering only in the exact squared-row encoding. Proposition
6.5k gives an exact Banaszczyk rounding that shields an arbitrary linear-size
anchor family, but its same-orientation vertex-cover condition remains open.
Proposition 6.5l extends the squared-row obstruction to every fixed raw
polynomial degree, and through half-degree
\((1-o(1))\log n/\log\log n\); it is a proof-system barrier, not an
orientation. Proposition 6.5m proves that a symmetric conference signing
cannot commute with any real skew signing: every diagonal entry of the
commutator is (2\pmod4).  This kills the exact commuting-conference mate,
but its (4n) Frobenius lower bound does not obstruct the approximate mate
required by Proposition 6.5e. Proposition 6.5n gives an infinite family of
complete signings with \(\Phi(A)=\Theta(n^{3/2})\) for which the coherent
clique-flip bound fails by a fixed multiple of \(n^{3/2}\).  Thus the bare
optimal-order scale cannot prove coherent control; the still-open
global-minimizer statement would have to use \(\Phi(A)=m_n\) essentially.
Proposition 6.6
proves that diamond off one explicit Hamming-central/
joint-energy residue. Proposition 6.7 gives an equal-endpoint tetrahedral
tripling frame and two exact spectral shields, but leaves their complement
open. Proposition 6.8 instead composes independently optimal orders \(n\)
and \(2n\): a bi-balanced Hadamard cross block proves the exact \(1:2\)
diamond outside a pair-central/high-joint-energy residue. Proposition 6.9
proves that the formerly proposed fixed-\(c\) free-energy route is false for
every fixed \(c>0\), including \(c=3\). Proposition 6.10 gives the corrected
critical-pressure gate: convergence of its optimized pressure on any
unbounded temperature set would imply convergence of \(\alpha_n\), but
common-temperature block interpolation has a nonvanishing critical-scale
defect and ordinary graphon convergence does not determine that pressure.
Proposition 6.10a disproves the proposed conference-product lower curve at
every positive temperature already at order five.  Its unconditional
entropy replacement has zero-temperature slope only \(1/\pi\), so it
repackages the known lower floor rather than closing the problem.
Neither required amplification ray is complete.

**Stronger value-specific target — not a theorem.** Sandwich + Paley
denseness (\(\rho=1\) on \(n=p^2+1\)) + E(1) on that family would give
\begin{equation}
\label{eq:limit}
L=\lim_{n\to\infty}\alpha_n=\tfrac12.
\end{equation}
E(1) is **not proved** for all primes \(p\ge5\); E(1) and \(L=1/2\) are not complete. Exactly one mathematical
obligation remains on this stronger Paley route (`GOAL.md`;
`evidence/share/denseness_path_package.md` § Caveats): residual (ii) for
even \(k\ge4p\). Proposition 15.720 closes the required bi-tight levels 2
and 3, and Proposition 15.750 closes multi-level Type I for every prime.
Thus the spectral floor, global QVAR, principal R1, and the historical
`3A+B` Type-I mechanism are no longer acceptance gates. Propositions
15.734--15.754 close the first three all-prime shells and several later
rows; Proposition 15.751 closes the fourth shell `k=4p+6` for every
`p>=13`. Proposition 15.752 then closes `k=4p+8` for every prime `p>=23`
and, without a prime or configuration census, a contiguous band through
`t=(p-9)/2` for `p=1 mod 4` and `t=(p-7)/2` for `p=3 mod 4`, where
`k=4p+2t`. Proposition 15.753 closes the two exceptional fifth-shell
endpoints `p=17,k=76` and `p=19,k=84` by exact aggregate-row certificates.
Proposition 15.754 closes the remaining endpoint `p=13,k=60,u=6` by an
exact aggregate/common-form certificate. Thus the fifth shell is closed for
every prime `p>=13`. Propositions 15.768--15.770 then close the first two
post-band layers in both generic congruence classes and the exceptional
`p=23,t=9,k=110` and `p=23,t=10,k=112` endpoints. The second p23 close has
eleven low roots and reuses the fixed exact 33,649-five-set quartic/octic
certificate. Proposition 15.771 closes `p=23,t=11,k=114` by the general
mean-46 equality classification and the all-boundary phase-zero mass-32
contradiction. Proposition 15.772 closes the third generic p1 layer
`p=1 mod 4,p>=29,t=q-1,k=5p-3` using the complement-triple punctured-gap
theorem and a new local mass-`p+11` exclusion. The exact live residual is
recorded in `STATUS.md`.
Proposition 15.755 further reduces any hypothetical dangerous descent on the
full Boolean cube to a sharp defect dichotomy `delta=2p` or
`delta>=6p-12` for `p>=11`.  Neither branch is empty: affine line aliases and
positive-triangle flips attain the two endpoints for every prime. Proposition
15.756 proves that arbitrary-boundary Weil cancellation cannot repair this:
its character cap is weaker than the trivial cap for every even boundary of
size at least four, and two parallel lines attain equality.  These are
global reductions and duplicate-work barriers, not a proof of residual (ii).
The live structural target is an integral common-edge/difference-Radon
theorem for one `0/1` graph `H`.
Proposition 15.757 computes the exact binary edge--Radon image and proves
that the recorded `p=1 mod 4` compact aggregate family passes it. Proposition
15.758 gives sharp
coefficient cancellation and two infinite local survivor rays. Proposition
15.759 finds the complete characteristic-\(p\) moment hierarchy, and
Proposition 15.760 proves that those moments cut out the unrestricted
integral image exactly, with elementary cokernel
\((\mathbf Z/p\mathbf Z)^{S(p)}\). Proposition 15.761 computes the full
real spectrum and proves that even the resulting Moore--Penrose test leaves
both compact rays alive with strict room. The ordered live structural target
is now narrower. For every prime \(p=4r+3\), \(r\ge7\), zero odd global forms
force a row with \(b\) arbitrary compact atoms and \(r-1\) all-equal atoms to
have a central signed edge chain whenever \(3b\le r+2\). This covers the
deterministic balanced branch-C allocation through
\(\delta=(2r+2)\lfloor(r+2)/3\rfloor\), but not unbalanced allocations. At
\(p=31\), an audited exact certificate for the one-compact/six-all-equal row
then excludes zero odd/degree-six/degree-eight global forms throughout the
balanced band \(69\le t\le99\). Separate exact
seven-channel Jacobian minors make the unrestricted degree-six/eight compact
and all-equal atom maps dominant over the algebraic closure, so no universal
polynomial identity or purely algebraic projective root-count contradiction
can finish the gate. The associated common forms and labels are
extension-valued, not admissible \(\mathbf F_p\) labels or coefficients. One
must still decide unbalanced zero-form allocations and nonzero even global
forms over \(\mathbf F_p\), coordinated across directions, and then whether
the resulting integral affine fibre meets
\(\prod_e\{0,\tau_e\}\). None of these follow-ups constructs a common simple
graph or closes residual (ii).
Lemma D existence / 2-plane is complete in `A3_PROOF.md`
and 15.276. Aut-Schur is **false**. Gsum unused.

**Optional still open:** Path-C residual / \(16N\) (independent).

**Corollary (\(\rho=1\) along a dense Paley family).** For every odd prime \(p\), the Paley conference
matrix of order \(n=p^2+1\) (over \(\mathbb F_{p^2}\)) admits a halfspace boolean eigenvector
\(Cx=px\), hence \(\rho(C)=1\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\). Along \(n_k=p_k^2+1\) one has
\(n_{k+1}/n_k\to1\) and \(\limsup_k\rho(C_{n_k})=1\). Proof: `evidence/PROOF_rho_eq_1.md`.
**If** \(m_n\ge\Phi(C)-2\) (E(1)) this would force \(L=\tfrac12\). E(1) is open.

**What is complete.** Dual-Gaussian lower bound; denseness; the two-ray
conditional convergence theorem; the exact directed-half-cut reformulation
of the multiplier-two diamond; \(\rho=1\) on \(n=p^2+1\);
15.720 degree-congruence exclusion of required bi-tight levels 2 and 3;
15.750 all-prime multi-level Type-I exclusion;
15.751 fourth-shell exclusion for every prime \(p\ge13\);
15.752 fifth-shell exclusion for every prime \(p\ge23\) and the displayed
contiguous higher-shell band;
15.753 exact endpoint exclusion at \((p,k)=(17,76),(19,84)\), hence
fifth-shell exclusion for every prime \(p\ge17\);
15.754 exact endpoint exclusion at \((p,k,u)=(13,60,6)\), hence
fifth-shell exclusion for every prime \(p\ge13\);
15.768 first `p=1 mod 4` post-band layer for every `p>=29`;
15.769 first `p=3 mod 4` post-band layer for every `p>=31` and the
exceptional `p=23,t=9,k=110` endpoint; 15.770 the next layer in both
generic classes and `p=23,t=10,k=112`; 15.771 the third exceptional p23
layer `t=11,k=114` by general-slice equality and a forced mass-32 row;
15.772 the third generic p1 layer `t=q-1,k=5p-3` for every `p>=29`,
with a repaired gap-two premise and the new gap-four equality;
15.755 full-cube dangerous-spike dichotomy and its sharp affine/triangle
counterfamilies; 15.756 arbitrary-boundary character-cap no-go;
15.757 exact binary edge--Radon image; 15.758 sharp coefficient cancellation
and two infinite local survivor rays; 15.759 complete characteristic-\(p\)
moment hierarchy; 15.760 exact integral image and open signed-Boolean-box
reduction; 15.761 full real edge--Radon spectrum and least-norm barrier; the
all-prime branch-C bounded-compact zero-odd-form centrality theorem and its
balanced initial band; the audited \(p=31\) one-row certificate and resulting
balanced \(69\le t\le99\) zero-form no-go; and the
seven-channel algebraic-dominance barrier, which supplies no \(\mathbf F_p\)
or Boolean lift;
residual (ii) affine + even \(k\le 4p-2\) (15.179/236/237), **not**
even \(k\ge4p\); residual (i) two-level Type I via 15.272
\(k=1\cup k=3\Rightarrow G_+\succ0\Rightarrow\ker=\mathrm{sc}\) (15.207)
\(\Rightarrow\) dual-eq empty (15.249/15.216). Aut-Schur remains false. Gsum
disj LB unused. **E(1) and \(L=1/2\) are not complete.**

On the positive `p=7` infinity-plus-seven front, Propositions 15.713--15.717
close `z=0,1,2,3`.  Propositions 15.718--15.719 rigorously reduce the
remaining `z=7` pointed systems and identify their projected high-catalog
semigroups through grade six, but exclude no source boundary.  All 56 actual
`z=7` line boundaries in two orbits remain open.

The intended Type I two-level close is 15.272 (\(k=1\cup k=3\) span \(\Rightarrow G_+\succ0\)), not Gsum and not Aut-Schur. **Status (2026-09-04; through audited record 15.772 and the post-15.761 support, conic, and Boolean reductions).** The direct multiplier-2 and multiplier-3 estimates are open: Proposition 6.6 narrows multiplier two to (6.20); Proposition 6.7 narrows direct tripling to an unshielded tetrahedral diamond; and Proposition 6.8 gives a different \(1:2\) reduction with residual (6.42)--(6.43). None closes its ray. Proposition 6.9 kills every fixed-temperature version of the signed-Eulerian free-energy target, so \(c=3\) is no longer a live fallback. Proposition 6.10 gives an exact optimized critical-pressure sufficient gate but proves that common-raw-temperature interpolation and ordinary graphon convergence do not establish it. Proposition 6.10a additionally disproves the proposed conference-product lower curve at every positive temperature; its universal entropy fallback has only the known \(1/\pi\) limiting slope. The natural rank-one phase self-gluing of the same two halves is also excluded by the exact amplification \(\nu_4(K\otimes H)=4\nu_4(H)\), versus the required \(2\sqrt2\). The optional Paley \(L=1/2\) route remains open at both residual (ii) and the minimal-four-gap implication bridge exposed by Proposition 15.764. Proposition 15.752 closes an infinite contiguous band of higher shells, Proposition 15.753 closes the \(p=17,19\) fifth-shell endpoints, and Proposition 15.754 closes the remaining \(p=13\) endpoint. Propositions 15.768--15.770 close the next two generic layers and `p=23,t=9,10`; Proposition 15.771 closes `p=23,t=11,k=114`, while Propositions 15.755--15.765 isolate surviving common-graph regimes and refute the attempted universal affine classification without closing either gate. Proposition 15.772 additionally closes `p=1 mod 4,p>=29,t=q-1,k=5p-3`. With `q=(p-1)/2`, the exact residual-(ii) frontier is `p=5,7`; `p=11,t>=3` (`k>=50`); `p=13,17,19,t>=5`; `p=23,t>=12`; `p=1 mod 4,p>=29,t>=q`; `p=3 mod 4,p>=31,t>=q`; and positive `p=7,z=7`. Details: `evidence/share/denseness_path_package.md` § Caveats.

The separate cross-rectangle Fourier calculation proves exact analytic
stability and Gram rigidity and exhibits an infinite Gram-perfect family at
the spectral maximum. It is a norm-only and Gram-only no-go, not a proof of
the multiplier-two estimate: the exact statewise diagonal-payment inequality
remains open. See
`evidence/NOTE_2026-09-03_CROSS_RECTANGLE_FOURIER_STABILITY.md`.

> **Reader note.** Soft-close of \(L=\tfrac12\) via scheme-min Gsum was
> **retracted** 2026-08-06. Aut-Schur remains false. Historical “\(L\) OPEN”
> remarks in Props 15.20–15.171 refer to those older routes. Likewise,
> statements inside Propositions 15.1--15.749 that multi-level Type I remains
> open are as-of records superseded by Proposition 15.750. Earlier statements
> listing generic branch-B \(t=3\) as open are superseded by Proposition
> 15.751, and statements listing every \(p\ge17,t\ge4\) layer as open are
> superseded on Proposition 15.752's explicit band and, at the two remaining
> fifth-shell endpoints, by Proposition 15.753. Statements that retain
> \(p=13,k=60,u=6\) are as-of records superseded by Proposition 15.754.

---

## §1. Equivalent matrix form

Associate to \((a_{ij})_{i<j}\) the symmetric zero-diagonal matrix \(A\) with \(A_{ij}=A_{ji}=a_{ij}\). Then
\[
\sum_{i<j}a_{ij}x_i x_j=\frac12\,x^\top A x,
\]
so
\[
m_n=\frac12\min_A\max_{x\in\{\pm1\}^n}\bigl|x^\top A x\bigr|
\]
over all such \(A\). Write \(\Phi(A)=\max_x\bigl|\tfrac12 x^\top A x\bigr|\), hence \(m_n=\min_A\Phi(A)\). The quantity \(\alpha_n\) is the same in both presentations.

---

## §2. Spectral upper bound for conference matrices

For every admissible \(A\), \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\), with equality iff \(A\) is a symmetric conference matrix of order \(n\). On the sphere of radius \(\sqrt n\),
\[
\max_{\|x\|_2=\sqrt n}\bigl|x^\top A x\bigr|=n\|A\|_{\mathrm{op}}\ge n\sqrt{n-1}.
\]
Hence for any conference matrix \(C\),
\begin{equation}
\label{eq:conf-ub}
\Phi(C)\le\frac12 n\sqrt{n-1},
\qquad
\frac{\Phi(C)}{n^{3/2}}\le\frac12\sqrt{1-\frac1n}.
\end{equation}
This upper-bounds \(\Phi(C)\) and therefore \(m_n\) when a conference matrix exists. It does **not** lower-bound \(m_n=\min_A\Phi(A)\).

---

## §3. Monotonicity and padding

**Proposition 3.1.** \(m_n\le m_N\) whenever \(2\le n\le N\).

*Proof.* Let \(A\) be optimal of order \(N\). For \(S\subset[N]\), \(|S|=n\), and \(y\in\{\pm1\}^S\), extend \(y\) by i.i.d. uniform random signs on the complement. Then \(\mathbb E Q=Q_S(y)\), so \(|Q_S(y)|\le m_N\). Hence \(m_n\le m_N\). \(\square\)

**Proposition 3.2.** \(m_{n+1}\le m_n+n\).

*Proof.* Adjoin a last row/column of all \(+1\) to an optimal order-\(n\) matrix. \(\square\)

**Proposition 3.3 (continuity of \(\alpha\)).** \(|\alpha_{n+1}-\alpha_n|\to0\) as \(n\to\infty\). In particular the set of limit points of \((\alpha_n)\) is a closed interval \([\lambda,\Lambda]\).

*Proof.* Write \(\delta_n:=m_{n+1}-m_n\in[0,n]\). Then
\[
\alpha_{n+1}-\alpha_n
=
\frac{\delta_n}{n^{3/2}}
-\frac32\frac{\alpha_n}{n}
+O\Bigl(\frac{\delta_n}{n^{5/2}}+\frac{\alpha_n}{n^2}\Bigr).
\]
The first term is \(\le n^{-1/2}\to0\). Boundedness \(\alpha_n=O(1)\) follows from the elementary random-method estimate \(m_n\le\sqrt{\log2}\,n^{3/2}\) for large \(n\) (union bound on Rademacher sums; cf. §7). A real sequence with consecutive gaps tending to zero has connected limit-point set, hence a closed interval. \(\square\)

**Corollary 3.4.** For \(d\ge0\), \(m_{n+d}\le m_n+dn+\tfrac12 d(d-1)\).

---

## §4. Upper bound \(\limsup\alpha_n\le\tfrac12\)

**Proposition 4.1.** \(\limsup_{n\to\infty}\alpha_n\le\tfrac12\).

*Proof.* For a prime \(q\equiv1\pmod4\), the Paley conference matrix of order \(n=q+1\) satisfies \(C^\top C=(n-1)I\) (Paley 1933), so \(m_n\le\tfrac12 n\sqrt{n-1}\) by \eqref{eq:conf-ub}. By Dirichlet’s theorem there are infinitely many such primes. Writing \(n_k=q_k+1\), the prime-number theorem in the progression \(1\bmod4\) gives \(n_k\sim 2k\log k\), hence \(n_{k+1}/n_k\to1\).

For \(N\ge2\) let \(n_k\) be the least Paley order \(\ge N\). Then \(m_N\le m_{n_k}\le\tfrac12 n_k\sqrt{n_k-1}\), so
\[
\alpha_N\le\frac12\Bigl(\frac{n_k}{N}\Bigr)^{3/2}\sqrt{1-\frac1{n_k}}.
\]
As \(N\to\infty\), \(n_k/N\to1\). \(\square\)

---

## §5. Lower bounds on \(\liminf\alpha_n\)

**Proposition 5.1 (Bohnenblust–Hille).** \(\liminf_{n\to\infty}\alpha_n\ge2^{-5/2}\).

*Proof.* The degree-\(2\) Bohnenblust–Hille inequality on the Boolean cube (Defant–Mastyło–Pérez, *Math. Ann.* 2019) yields a universal \(B_2<\infty\) such that for \(f=\sum_{i<j}c_{ij}x_ix_j\),
\[
\Bigl(\sum_{i<j}|c_{ij}|^{4/3}\Bigr)^{3/4}
\le B_2\max_{x\in\{\pm1\}^n}|f(x)|.
\]
With \(|c_{ij}|=1\) and \(\max|f|=m_n\), one has \(\binom{n}{2}^{3/4}\le B_2 m_n\). The hypercontractive bound \(B_2\le2\sqrt2\) (Ivanisvili, MO 413935) yields \(\liminf\alpha_n\ge2^{-5/2}\). \(\square\)

**Proposition 5.2 (dual-Gaussian arcsine; universal \(1/\pi\)).** For every Seidel matrix \(A\in\mathcal S_n\) and every \(n\ge2\),
\begin{equation}
\label{eq:dual-gauss}
\Phi(A)\;\ge\;\frac{n\sqrt{n-1}}{\pi}.
\end{equation}
In particular \(m_n\ge n\sqrt{n-1}/\pi\) and
\[
\liminf_{n\to\infty}\alpha_n\;\ge\;\frac1\pi.
\]

*Proof.* Fix \(A\in\mathcal S_n\), let \(g\sim N(0,I_n)\), and set \(t=\sqrt{n/(n-1)}\). Define
\[
z^{\pm}=\Bigl(I\pm\frac{t}{\sqrt n}\,A\Bigr)g,
\qquad
x^{\pm}=\operatorname{sgn}(z^{\pm})
\]
(with the convention \(\operatorname{sgn}(0)=+1\)). Each coordinate of \(z^{\pm}\) has variance
\[
d=1+t^2\frac{n-1}{n}=2.
\]
For \(i\neq j\) the correlations \(r_{ij}^{\pm}\) of the pairs \((z_i^{\pm},z_j^{\pm})\) satisfy
\[
a_{ij}\,r_{ij}^{\pm}
=
u_{ij}\pm v,
\qquad
u_{ij}
=
\frac{a_{ij}\,t^2(A^2)_{ij}}{nd},
\qquad
v
=
\frac{2t}{\sqrt n\,d}
=
\frac1{\sqrt{n-1}}.
\]
(The identities use \(A_{ii}=0\) and \(A_{ij}=a_{ij}=\pm1\).) Because \(|(A^2)_{ij}|\le n-2\),
\[
|u_{ij}|\le\frac{n-2}{2(n-1)},
\qquad
|u_{ij}\pm v|
\le
\frac{n-2}{2(n-1)}+\frac1{\sqrt{n-1}}
<1
\]
for all \(n\ge2\) (direct check at \(n=2,3\); for \(n\ge4\) the right-hand side is \(<1\)).

The arcsine law for a centered bivariate Gaussian of correlation \(r\) states
\(\mathbb E[\operatorname{sgn}Z_1\operatorname{sgn}Z_2]=(2/\pi)\arcsin r\). Hence
\[
\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
=
\frac2\pi\sum_{i<j}a_{ij}\bigl(\arcsin r_{ij}^+-\arcsin r_{ij}^-\bigr).
\]
For each edge, \(a_{ij}(\arcsin r_{ij}^+-\arcsin r_{ij}^-)=\arcsin(u_{ij}+v)-\arcsin(u_{ij}-v)\):
if \(a_{ij}=+1\) this is immediate, and if \(a_{ij}=-1\) both sides pick up a matching sign change.
Since \((\arcsin)'(r)=1/\sqrt{1-r^2}\ge1\) on \((-1,1)\), the mean-value theorem yields
\[
\arcsin(u+v)-\arcsin(u-v)\ge 2v.
\]
Therefore
\[
\mathbb E Q_A(x^+)-\mathbb E Q_A(x^-)
\ge
\frac2\pi\cdot\binom n2\cdot 2v
=
\frac{2n\sqrt{n-1}}{\pi},
\]
since \(v=1/\sqrt{n-1}\) and \(\binom n2=n(n-1)/2\).
Both expectations lie in \([-\Phi(A),\Phi(A)]\), so their difference is at most \(2\Phi(A)\).
Halving gives \eqref{eq:dual-gauss}. Minimising over \(A\) and passing to the limit produces the
claim on \(\alpha_n\). \(\square\)

**Remark (correction).** The classical single-sided Nesterov bound applied to a *conference*
matrix \(C\) yields only \(\Phi(C)\ge n\sqrt{n-1}/\pi\), which does not control \(m_n=\min\Phi\).
Proposition 5.2 upgrades the same constant to a **uniform** lower bound on every Seidel matrix
by running the Gaussian construction in the pair of directions \(\pm A\) and comparing the two
expectations. The older Bohnenblust–Hille floor \(2^{-5/2}\approx0.177\) is thereby improved to
\(1/\pi\approx0.3183\).

**Remark (cut-code form).** Writing \(m=\binom n2\) and
\(D_n=\{(\pm x_ix_j)_{i<j}:x\in\{\pm1\}^n\}\subset\{\pm1\}^m\), the identity
\(\langle a,c_x\rangle=m-2d_H(a,c_x)\) and the presence of both \(c_x\) and \(-c_x\) in \(D_n\) give
\begin{equation}
\label{eq:cut-code}
m_n
=
\binom n2-2\rho(D_n),
\end{equation}
where \(\rho(D_n)\) is the Hamming covering radius of \(D_n\) in \(\{\pm1\}^m\). Equivalently,
\(m_n/2\) is the covering-radius deficit of the antipodal cut code of \(K_n\). This is the coding
formulation used by Esmaeili–Zaghian (2009); it does not by itself settle existence of
\(\lim\alpha_n\).

**Proposition 5.3 (sharp equimodular \(L^1\)-influence form).**  Let
\(\mu_k=\mathbb E|\varepsilon_1+\cdots+\varepsilon_k|\), and define
\[
K_n=\max_{A\in\mathcal S_n}
 \operatorname{Inf}^{(1)}\!\left({Q_A\over\Phi(A)}\right),
\qquad
\operatorname{Inf}^{(1)}(f)=
 \sum_{i=1}^n\mathbb E_x
 \left|{f(x)-f(x^{\oplus i})\over2}\right|.
\]
Then
\begin{equation}
\label{eq:sharp-influence}
m_n={n\mu_{n-1}\over K_n}.
\end{equation}
Consequently, \(m_n/n^{3/2}\) converges if and only if \(K_n\) converges,
and the two limits, when they exist, obey
\[
 \lim_{n\to\infty}{m_n\over n^{3/2}}
 ={\sqrt{2/\pi}\over\lim_{n\to\infty}K_n}.
\]

*Proof.*  The discrete derivative of an equimodular quadratic form is
\[
 {Q_A(x)-Q_A(x^{\oplus i})\over2}
 =x_i\sum_{j\ne i}a_{ij}x_j.
\]
Under uniform \(x\), its absolute value has mean \(\mu_{n-1}\), independently
of both \(i\) and \(A\).  Hence
\(\operatorname{Inf}^{(1)}(Q_A)=n\mu_{n-1}\), and maximizing the normalized
influence is exactly the same as minimizing \(\Phi(A)\).  This proves
\eqref{eq:sharp-influence}.  Finally
\(\mu_{n-1}/\sqrt n\to\sqrt{2/\pi}\), while the already-proved positive lower
and upper bounds on \(m_n/n^{3/2}\) keep \(K_n\) bounded above and away from
zero.  The equivalence of convergence follows. \(\square\)

**Remark.** General bounded-degree \(L^1\)-influence theorems give
dimension-free bounds for \(K_n\), but no known result compares \(K_n\) at
different dimensions.  Proposition 5.3 is an exact replacement target for
the original question, not a proof of the limit.

---

## §6. Denseness

**Proposition 6.1.** If \(n_{k+1}/n_k\to1\), then
\[
\liminf_n\alpha_n=\liminf_k\alpha_{n_k},
\qquad
\limsup_n\alpha_n=\limsup_k\alpha_{n_k}.
\]

*Proof.* For \(n_k\le N\le n_{k+1}\), monotonicity yields \(m_{n_k}\le m_N\le m_{n_{k+1}}\), so
\[
\alpha_{n_k}\Bigl(\frac{n_k}{N}\Bigr)^{3/2}
\le\alpha_N\le
\alpha_{n_{k+1}}\Bigl(\frac{n_{k+1}}{N}\Bigr)^{3/2}.
\]
Both ratios tend to \(1\). \(\square\)

In particular, for each fixed integer \(r\ge1\), \(\liminf_n\alpha_{rn}=\liminf\alpha_n\) and \(\limsup_n\alpha_{rn}=\limsup\alpha_n\).

**Proposition 6.2 (Paley reduction).** Let \(q_k\) be the \(k\)-th prime congruent to \(1\pmod4\), and set \(n_k=q_k+1\) (Paley/conference orders). Then \(n_{k+1}/n_k\to1\), and therefore
\[
\liminf_{n\to\infty}\alpha_n=\liminf_{k\to\infty}\alpha_{n_k},
\qquad
\limsup_{n\to\infty}\alpha_n=\limsup_{k\to\infty}\alpha_{n_k}.
\]
In particular, \(\lim\alpha_n\) exists if and only if \(\lim_k\alpha_{n_k}\) exists.

*Proof.* By the prime-number theorem for the arithmetic progression \(1\bmod4\) (Dirichlet density \(\tfrac12\)), one has \(q_k\sim 2k\log k\), so \(n_{k+1}/n_k\to1\). Apply Proposition 6.1. \(\square\)

Thus any existence proof may restrict attention to Paley orders; any non-existence proof must already be visible along that sparse but ratio-dense subsequence.

**Proposition 6.3 (two-ray convergence criterion).** Put
\[
H(n)=m_n^{2/3},\qquad h(n)=\frac{H(n)}n=\alpha_n^{2/3}.
\]
Let \(\eta(n)\ge0\), put
\begin{equation}
\eta^*(N)=\sup_{u\ge N}\eta(u),\qquad
E(N)=\sum_{j\ge0}\eta^*(2^jN),                       \tag{6.1}
\end{equation}
and suppose \(E(N)\to0\).  If, for all sufficiently large \(n\),
\begin{equation}
H(2n)\le2H(n)+2n\eta(n),\qquad
H(3n)\le3H(n)+3n\eta(n),                            \tag{6.2}
\end{equation}
then \(\lim_n\alpha_n\) exists.

*Proof.* Monotonicity of \(m_n\) makes \(H\) nondecreasing.  Dividing each
inequality in (6.2) by its new argument gives
\(h(qn)\le h(n)+\eta(n)\), for \(q=2,3\).  Along any word in \(2,3\), the
size before step \(j\) is at least \(2^jn\).  Therefore, uniformly in
\(a,b\ge0\),
\[
h(2^a3^b n)\le h(n)+E(n).                            \tag{6.3}
\]
The consecutive ratios in the sorted semigroup
\(\{2^a3^b:a,b\ge0\}\) tend to one.  Indeed,
\(\log2/\log3\) is irrational, so the corresponding irrational rotation is
dense; a finite prefix is an arbitrarily fine net modulo \(\log3\), and a
large nonnegative multiple of \(\log3\) then places a semigroup logarithm
in every sufficiently large interval of the prescribed length.

Fix \(n\) beyond the threshold in (6.2).  For every sufficiently large \(N\),
choose \(s=2^a3^b\) with \(N\le sn=(1+o_N(1))N\).  By monotonicity and
(6.3),
\[
h(N)\le\frac{sn}{N}h(sn)
\le(1+o_N(1))(h(n)+E(n)).
\]
Hence \(\limsup_Nh(N)\le h(n)+E(n)\).  Letting \(n\) tend to
infinity along a liminf subsequence gives \(\limsup h\le\liminf h\).
Therefore \(h\), and hence \(\alpha_n=h(n)^{3/2}\), converges. \(\square\)

The power-saving hypothesis
\[
H(2n)\le2H(n)+O(n^{1-\delta}),\qquad
H(3n)\le3H(n)+O(n^{1-\delta})                       \tag{6.4}
\]
for any \(\delta>0\) is a special case, but (6.1) is weaker: for example,
errors \(O(n/(\log n)^{1+\varepsilon})\) in (6.2) suffice.  It is enough to
prove the doubling estimate in (6.4) together with the single \(1:2\) split
\[
H(3n)\le H(n)+H(2n)+O(n^{1-\delta});
\]
the doubling inequality then supplies the tripling estimate.
Equivalently, using \(m_n=\Theta(n^{3/2})\) (Proposition 5.2 for the lower
bound and the elementary random-method upper bound recalled in Sections 4
and 7), it is enough to construct power-saving amplifications
\[
m_{rn}\le r^{3/2}m_n+O(n^{3/2-\delta}),\qquad r=2,3. \tag{6.5}
\]
Only these two fixed composition shapes are required, rather than a uniform
all-pairs composition theorem.  More generally, (6.5)'s relative error may
be replaced by any modulus satisfying the tail-summability condition (6.1).

*Sharpness of the two-ray premise.*  Doubling alone is not a convergence
criterion.  For sufficiently small \(\varepsilon>0\),
\[
G(x)=x\bigl(1+\varepsilon\sin(2\pi\log_2x)\bigr)
\]
is increasing, has bounded increments along the integers, and satisfies
\(G(2x)=2G(x)\), while \(G(n)/n\) oscillates.  The second,
multiplicatively independent scale is load-bearing.

**Proposition 6.4 (exact four-state form of Hadamard doubling).**  Pair the
vertices of an order-\(2n\) signing into \(n\) clouds of size two.  Suppose
every inter-cloud \(2\times2\) sign block \(B_{ij}\) is Hadamard (equivalently,
the product of its four entries is \(-1\)); allow an arbitrary sign \(d_i\)
on the internal edge of cloud \(i\).  Write a cloud state uniquely as
\[
z_i=s_i(1,(-1)^{t_i}),\qquad s_i\in\{\pm1\},\quad t_i\in\mathbb F_2.
\]
For each oriented block there are unique
\(\alpha_{ij},\beta_{ij},\gamma_{ij}\in\mathbb F_2\) such that
\[
\frac12z_i^\top B_{ij}z_j
=s_is_j(-1)^{t_it_j+\alpha_{ij}t_i+
                 \beta_{ij}t_j+\gamma_{ij}}.          \tag{6.6}
\]
For \(t\in\mathbb F_2^n\), let \(C_t\) be the order-\(n\) signing with
\[
(C_t)_{ij}=(-1)^{t_it_j+\alpha_{ij}t_i+
                 \beta_{ij}t_j+\gamma_{ij}},
\qquad K=\max_t\Phi(C_t).
\]
Then the full lift \(S\) satisfies
\[
2K-n\le\Phi(S)\le2K+n.                              \tag{6.7}
\]

*Proof.*  The normalized four-entry table
\(\tfrac12(1,(-1)^{t_i})B_{ij}(1,(-1)^{t_j})^\top\)
is sign-valued and has entry product \(-1\), which is exactly the unique
Boolean normal form in (6.6).  Thus
\[
Q_S(z)=2Q_{C_t}(s)+\sum_i d_i(-1)^{t_i}.
\]
The last sum has absolute value at most \(n\); maximizing first over \(s\)
and then over \(t\) proves (6.7). \(\square\)

The quadratic term \(t_it_j\) is invariant under row and column gauges.  It
cannot be replaced by vertex factors: a table
\(c_{ij}r_i(t_i)r_j(t_j)\) has entry product \(+1\), whereas every Hadamard
block in (6.6) has product \(-1\).  Thus no edgewise local vertex-gauge
factorization can collapse these tables to one base signing.  This does not
exclude a special frame whose induced signings are related by a nonlocal
coincidence.

More precisely, minimize \(K\) over all Hadamard frames and call the result
\(G_2(n)\); minimize \(\Phi(S)\) over those frames and all internal signs and
call it \(L_2(n)\).  Equation (6.7) gives
\[
|L_2(n)-2G_2(n)|\le n.                              \tag{6.8}
\]
Since this is a subclass of all order-\(2n\) signings,
\(m_{2n}\le L_2(n)\), and consequently
\(G_2(n)\ge(m_{2n}-n)/2\).
Hence this class reaches the multiplier-two threshold precisely by proving,
up to the admissible \(O(n)\) term,
\[
G_2(n)\le\sqrt2\,m_n+n^{3/2}\omega(n),              \tag{6.9}
\]
where \(\omega\ge0\) and
\(\sum_{j\ge0}\sup_{u\ge2^jn}\omega(u)\to0\).

For each fixed oriented frame there is also a necessary condition.  Put
\(P(A)=\max_xQ_A(x)\), \(N(A)=-\min_xQ_A(x)\), let \(G=C_{\mathbf0}\), and
let \(J=C_{\mathbf1}\).  For every \(T\subseteq[n]\),
\[
K\ge\max\left\{
P(J[T])+P(G[T^c]),\;
N(J[T])+N(G[T^c])
\right\}.                                           \tag{6.10}
\]
Indeed, fix \(t=1_T\), independently maximize the two within-part forms,
and compare a spin state with the state obtained by flipping every spin in
\(T\).  The within-part sum is fixed and the cross term changes sign, so one
choice has absolute energy at least that sum.  The same argument with both
within-part minima gives the second term.  Consequently (6.9) requires a
simultaneous hereditary bound for the two endpoint signings \(G,J\); mixed
cross terms within that frame cannot repair a cut that violates it.  A
different orientation generally changes \(G,J\).  Proposition 6.5 below
shows, however, that the equal-endpoint skew frame has \(G=J=A\), and then
every inequality in (6.10) is automatic.  Thus endpoint selection is not the
live obstruction.  Only the mixed states remain.

**Proposition 6.5 (equal-endpoint skew reduction).**  Let \(A\) be an
order-\(n\) signing and let \(R\) be a skew signing: \(R_{ii}=0\),
\(R_{ji}=-R_{ij}\), and \(R_{ij}\in\{\pm1\}\) for \(i<j\).  For each
\(i<j\), orient the Hadamard block as
\[
 B_{ij}=\begin{pmatrix}A_{ij}&-R_{ij}\\R_{ij}&A_{ij}\end{pmatrix}.
                                                               \tag{6.11}
\]
Both endpoint signings of Proposition 6.4 are exactly \(A\).  Moreover, if
\(K(A,R)\) denotes the resulting four-state minimax, then
\[
 \boxed{K(A,R)={1\over2}\max_{x,y\in\{\pm1\}^n}
 \left(\lvert Q_A(x)+Q_A(y)\rvert+\lvert x^TRy\rvert\right).}   \tag{6.12}
\]
Every hereditary endpoint inequality (6.10) is automatic when
\(M=\Phi(A)\).  Consequently the multiplier-two ray follows if, for an
optimal \(A\), one can choose \(R\) so that
\[
 \lvert Q_A(x)+Q_A(y)\rvert+\lvert x^TRy\rvert
 \le2\sqrt2\,M+n^{3/2}\Omega(n)                  \tag{6.13}
\]
for every Boolean pair, where the supremum envelope of \(\Omega\) has a
vanishing dyadic Dini tail.

*Proof.*  The four normalized block values are
\[
 {1\over2}(1,(-1)^u)B_{ij}(1,(-1)^v)^T
 =\begin{cases}
 A_{ij},&u=v,\\
 R_{ij},&(u,v)=(0,1),\\
 -R_{ij},&(u,v)=(1,0).
 \end{cases}                                      \tag{6.14}
\]
Thus both endpoints are \(A\).  Given a cloud state \((s,t)\), put
\(x=s\) and \(y_i=(-1)^{t_i}s_i\).  Direct expansion gives
\[
 Q_{C_t}(s)={1\over2}
 \bigl(Q_A(x)+Q_A(y)-x^TRy\bigr).
\]
Every Boolean pair \((x,y)\) arises in this way.  Replacing \(y\) by
\(-y\) fixes both quadratic terms and reverses the skew term, while
\(\max(|a+b|,|a-b|)=|a|+|b|\).  This proves (6.12).

For the positive hereditary inequality, independently choose spins on
\(T\) and \(T^c\) attaining \(P(A[T])\) and \(P(A[T^c])\).  If their
within-part sum is \(D\) and their cross energy is \(X\), the two full
states obtained by flipping one part have energies \(D+X\) and \(D-X\).
Both have absolute value at most \(M\), so
\(D+|X|\le M\).  The same argument with the two within-part minima proves
the negative inequality.  Since (6.12) also gives \(K(A,R)\ge M\), (6.10)
is automatic.  Finally (6.13), (6.12), and (6.7) give
\[
 m_{2n}\le2\sqrt2\,m_n+n^{3/2}\Omega(n)+n,
\]
which is the required Dini-summable doubling estimate. \(\square\)

**Proposition 6.5a (exact directed-half-cut reformulation).**  In the setting
of Proposition 6.5 put \(S=A\circ R\), so that \(S\) is a tournament sign
matrix, with \(S_{uv}=1\) interpreted as the arc \(u\to v\).  For
\(U\subseteq[n]\), let
\[
 F_S(U)=\{\{u,v\}:u\in U, v\notin U, S_{uv}=1\},
\]
and write \(A^F\) for the symmetric signing obtained by reversing precisely
the edges in \(F\).  Then
\[
 \boxed{K(A,R)=\max_{U\subseteq[n]}\Phi\bigl(A^{F_S(U)}\bigr).} \tag{6.14a}
\]
Consequently
\[
 {1\over2}\min_R\max_{x,y}
 \bigl(|Q_A(x)+Q_A(y)|+|x^TRy|\bigr)
 =\min_{\substack{S\ {\rm tournament}}}\max_U
   \Phi\bigl(A^{F_S(U)}\bigr). \tag{6.14b}
\]

*Proof.*  Fix \(x,y\), put \(U=\{i:x_i=-y_i\}\), and split the
\(A_{ij}y_iy_j\)-energy into \(I,F,G\), respectively on non-cut edges,
arcs from \(U\) to \(U^c\), and arcs from \(U^c\) to \(U\).  Directly,
\[
 Q_A(y)=I+F+G,\qquad Q_A(x)=I-F-G,\qquad x^TRy=2(G-F).
\]
Thus
\[
 {1\over2}\bigl(|Q_A(x)+Q_A(y)|+|x^TRy|\bigr)
 =|I|+|G-F|
 =\max\{|Q_{A^{F_S(U)}}(y)|,|Q_{A^{G_S(U)}}(y)|\}.
\]
Switching all vertices of \(U\) takes \(A^{F_S(U)}\) to
\(A^{G_S(U)}\), so the two signings have the same \(\Phi\)-norm.
Maximizing over \(y,U\) proves (6.14a), and skew signings \(R\) correspond
bijectively to tournaments \(S=A\circ R\), proving (6.14b). \(\square\)

The same calculation isolates exactly which plotted pairs are dangerous.
Put
\[
 \epsilon(x,y)=M-\max\{|Q_A(x)|,|Q_A(y)|\}.
\]
Since \(|a+b|+|a-b|=2\max(|a|,|b|)\),
\[
 2M-|Q_A(x)+Q_A(y)|
 =|Q_A(x)-Q_A(y)|+2\epsilon(x,y).              \tag{6.14a1}
\]
For the directed cut energies \(F,G\) in the proof,
\[
 |x^TRy|-|Q_A(x)-Q_A(y)|
 =\begin{cases}
 4\min(|F|,|G|),&FG<0,\\
 -4\min(|F|,|G|),&FG\ge0.
 \end{cases}                                    \tag{6.14a2}
\]
Thus same-sign directed halves are harmless.  Up to the required Dini
error, the entire multiplier-two condition on the opposite-sign pairs is
exactly
\[
 2\min(|F|,|G|)\le(\sqrt2-1)M+\epsilon(x,y).     \tag{6.14a3}
\]
This is the sign-sensitive boundary seen by the two-half calculator; it
retains both payments that are lost in an independent norm estimate for
\(R\).

There is a useful energy-layer version.  Write
\(e(z)=M-|Q_A(z)|\).  The exact free payment in (6.14a1) is
\[
 |Q_A(x)-Q_A(y)|+2\min(e(x),e(y))
 =\begin{cases}
 e(x)+e(y),&Q_A(x)Q_A(y)\ge0,\\
 2M-|e(x)-e(y)|,&Q_A(x)Q_A(y)\le0.
 \end{cases}                                      \tag{6.14a4}
\]
(The two expressions agree when one energy is zero.)  If
\(t=(\sqrt2-1)M\), the zero-error target is therefore exactly
\[
 |x^TRy|\le2t+e(x)+e(y)                            \tag{6.14a5}
\]
on every same-sign energy pair, together with the second line of
(6.14a4) on opposite-sign pairs.  Since every skew signing obeys
\(\max_{x,y}|x^TRy|\ge n\mathbb E|S_{n-1}|\), any skew-extremizing pair
whose two energies have the same sign must have
\[
 e(x)+e(y)\ge n\mathbb E|S_{n-1}|-2t.
\]
Using the known \(M\le(1/2+o(1))n^{3/2}\), the right side is at least
\[
 \bigl(\sqrt{2/\pi}-(\sqrt2-1)-o(1)\bigr)n^{3/2}
 =\bigl(0.383671\ldots-o(1)\bigr)n^{3/2}.          \tag{6.14a6}
\]
Thus a successful orientation must push all of its largest bilinear
witnesses away from simultaneous same-sign near-maximizer layers.  This is
a joint energy-layer problem, not another global bound for \(R\).

This turns the multiplier-two step into an exact simultaneous-neighbor
problem: orient the edges so that reversing the outward-oriented half of
every cut keeps the norm below
\(\sqrt2\Phi(A)+o_{\rm Dini}(n^{3/2})\).  If \(A\) is globally optimal,
minimality gives only the reverse inequality
\(\Phi(A^{F_S(U)})\ge\Phi(A)\).  Hence (6.14a) exposes a sharper structural
target but does not itself close the ray.

**Proposition 6.5b (sharpness of the directed-half-cut multiplier).**  Put
\[
 T_n=\min_{\Phi(A)=m_n}\ \min_{S\ {\rm tournament}}\
       \max_{U\subseteq[n]}\Phi\bigl(A^{F_S(U)}\bigr).
\]
Then, for every fixed \(n_0\ge2\),
\[
 \limsup_{j\to\infty}{T_{2^j n_0}\over m_{2^j n_0}}\ge\sqrt2. \tag{6.14c}
\]
In particular, no fixed \(c<\sqrt2\) and errors
\(r_n=o(n^{3/2})\) can give \(T_n\le c m_n+r_n\) on every sufficiently
large order (or even on one complete eventual dyadic tail).

*Proof.*  If such a tail existed, Proposition 6.5's order-\(2n\) lift and
Proposition 6.5a would give
\[
 m_{2n}\le2T_n+n,
 \qquad
 \alpha_{2n}\le{c\over\sqrt2}\alpha_n
       +{r_n\over\sqrt2 n^{3/2}}+{1\over2^{3/2}\sqrt n}.
\]
The last two terms tend to zero and \(c/\sqrt2<1\), so iteration forces
\(\alpha_{2^j n_0}\to0\).  This contradicts Proposition 5.2's bound
\(\alpha_n\ge\sqrt{1-1/n}/\pi\). \(\square\)

Thus the factor \(\sqrt2\) is load-bearing.  Proposition 6.5a is not made
easier by seeking a uniformly smaller constant; the exact critical line,
with a sufficiently controlled error, is the one that must be reached.

**Proposition 6.5c (opposite-diagonal diamond and hybrid-slice complexification).**
Let \(A\) be an
order-\(n\) signing and let \(C\in\{\pm1\}^{n\times n}\) be arbitrary.  Then
the order-\(2n\) signing
\[
 {\cal K}(A,C)=\begin{pmatrix}A&C\\ C^T&-A\end{pmatrix}
\]
satisfies the exact diamond identity
\[
 \boxed{\Phi({\cal K}(A,C))=
 \max_{x,y\in\{\pm1\}^n}
 \bigl(|Q_A(x)-Q_A(y)|+|x^TCy|\bigr).}            \tag{6.14d}
\]
Consequently the multiplier-two ray would follow if every optimal \(A\),
with \(M=\Phi(A)\), admitted a cross signing \(C\) for which the right side
of (6.14d) were at most
\(2\sqrt2M+n^{3/2}\Omega(n)\), with a vanishing dyadic Dini tail.

There is an exact four-label form for every cross block, not only symmetric
ones.  For \(T\subseteq[n]\), put \(r_i=-1\) on \(T\) and \(r_i=1\)
off \(T\), let \(D_T=\operatorname{diag}(r)\), and define the zero-diagonal
symmetric weighted matrix \(G_T=G_T(A,C)\) by
\[
 (G_T)_{ij}=A_{ij}{1-r_ir_j\over2}
             +{C_{ij}r_j+C_{ji}r_i\over2},\qquad i<j,
 \quad
 h_T={1\over2}\sum_i C_{ii}r_i.                 \tag{6.14d1}
\]
Then
\[
 \boxed{\Phi({\cal K}(A,C))
 =2\max_{T,s}|Q_{G_T}(s)+h_T|.}                  \tag{6.14d2}
\]
If \(A_\delta(T)\) retains the \(A\)-edges crossing \((T,T^c)\) and is
zero elsewhere, complementing \(T\) negates the displacement
\(G_T-A_\delta(T)\) and \(h_T\), while fixing \(A_\delta(T)\).  Hence the
equivalent midpoint--displacement identity is
\[
 {1\over2}\Phi({\cal K}(A,C))
 =\max_{T,s}\left(
 |Q_{A_\delta(T)}(s)|+{1\over2}|s^TCD_Ts|\right). \tag{6.14d3}
\]

Each unordered edge has exactly four possible cross labels.  If
\(C_{ij}=C_{ji}=\varepsilon\), its weight in \(G_T\) is \(\varepsilon\)
when both endpoints lie outside \(T\), \(-\varepsilon\) when both lie
inside, and \(A_{ij}\) when it crosses.  If
\(C_{ij}=-C_{ji}\), orient \(u\to v\) when \(C_{uv}=A_{uv}\): its weight
is \(2A_{uv}\) exactly when that arc exits \(T\), and zero otherwise.
Consequently \(G_T\) is an ordinary signing for some (equivalently every)
\(T\) exactly when the off-diagonal part of \(C\) is symmetric.  Allowing
the other two labels therefore creates genuine \(0,\pm2\) weighted slices,
so global minimality of \(A\) cannot be applied to them as though they were
signings.

The all-directed specialization is nevertheless a clean surviving target.
Let \(S\) be a tournament, put \(R=A\circ S\), choose any diagonal signing
\(D\), and define
\[
 D_{\to}(A,S)=\max_{T,s}
 \left|\sum_{\substack{u\in T,\ v\notin T\\u\to v}}
 A_{uv}s_us_v\right|.
\]
The preceding edge table gives
\[
 \boxed{
 4D_{\to}(A,S)=\max_{x,y}
 \bigl(|Q_A(x)-Q_A(y)|+|x^TRy|\bigr),
 }                                                          \tag{6.14d4}
\]
where the corresponding opposite block has zero entries on its matching.
Adding any diagonal signing \(D\) makes it a complete signing and gives
\[
 \left|\Phi({\cal K}(A,R+D))-4D_{\to}(A,S)\right|\le n.     \tag{6.14d5}
\]
Therefore an orientation satisfying
\[
 D_{\to}(A,S)\le {M\over\sqrt2}+o_{\rm Dini}(n^{3/2})       \tag{6.14d6}
\]
would close the multiplier-two step.  This asks for one signed outgoing
half of every cut, rather than the full half-cut-neighbor norm in
Proposition 6.5a.

Equivalently, for each \(T\) let \(B_T(S)\) be the
\(|T|\)-by-\(|T^c|\) matrix whose \((u,v)\) entry is \(A_{uv}\) when
\(u\to v\), and zero otherwise.  Then
\[
 D_{\to}(A,S)=\max_T\|B_T(S)\|_{\infty\to1},
 \qquad
 \|A_{T,T^c}\|_{\infty\to1}\le M.              \tag{6.14d6a}
\]
The second inequality follows by completing arbitrary independent row and
column signs to the two cube states that differ exactly on \(T\), then using
the endpoint diamond.  Thus (6.14d6) is a simultaneous one-sided
\(1/\sqrt2\)-paving theorem for every signed cut submatrix, with all pavings
required to come from one tournament.  For a fixed cut and a fixed state,
its two directed energies can always be balanced to
\(\lceil|F+G|/2\rceil\le\lceil M/2\rceil\); global consistency over every
cut and state is the entire obstruction.

For a necessary best-response test, choose \(X_R(y)\) so that
\(X_R(y)^TRy=\|Ry\|_1\).  The exact identity (6.14d4) would require
\[
 |Q_A(X_R(y))-Q_A(y)|+\|Ry\|_1
 \le2\sqrt2M+o_{\rm Dini}(n^{3/2})              \tag{6.14d7}
\]
for every \(y\).  Since \(\mathbb E\|RY\|_1=n\mathbb E|S_{n-1}|\),
averaging forces
\[
 \mathbb E|Q_A(X_R(Y))-Q_A(Y)|
 \le2\sqrt2M-n\mathbb E|S_{n-1}|+o_{\rm Dini}(n^{3/2}).    \tag{6.14d8}
\]
This is a concrete joint energy--skew gate for any proposed orientation.

No proof may budget the displacement in (6.14d3) independently: exactly
\[
 \max_{T,s}{1\over2}|s^TCD_Ts|
 ={1\over2}\|C\|_{\infty\to1}
 \ge {n\over2}\mathbb E|S_n|
 =\left({1\over\sqrt{2\pi}}+o(1)\right)n^{3/2}.   \tag{6.14d9}
\]
At the known scale this already exceeds the free
\((\sqrt2-1)M\) budget.  Statewise cancellation with the midpoint in
(6.14d3), not two separate estimates, is again load-bearing.

There is also an exact symmetric subclass.  Let \(C_0\) be a symmetric
zero-diagonal signing, let \(D=\operatorname{diag}(d_1,\ldots,d_n)\) with
\(d_i\in\{\pm1\}\), and, for \(T\subseteq[n]\), define \(B_T=B_T(A,C_0)\)
by
\[
 (B_T)_{ij}=\begin{cases}
 (C_0)_{ij},&i,j\notin T,\\
 A_{ij},&|\{i,j\}\cap T|=1,\\
 -(C_0)_{ij},&i,j\in T.
 \end{cases}                                                \tag{6.14e}
\]
Then
\[
 \boxed{\left|\Phi({\cal K}(A,C_0+D))
       -2\max_T\Phi(B_T(A,C_0))\right|\le n.}                \tag{6.14f}
\]
In the coherent case \(C_0=A\), put
\(L_{\rm cl}(A)=\max_T\Phi(A^{K_T})\), where \(A^{K_T}\) reverses
the clique on \(T\).  One has
\[
 2L_{\rm cl}(A)-n\le\Phi({\cal K}(A,A+D))
 \le2L_{\rm cl}(A)+n,                                      \tag{6.14g}
\]
and, for \(\widetilde Q_A(z)=\sum_{i<j}A_{ij}z_iz_j\),
\[
 L_{\rm cl}(A)=\max_{z\in\{\pm1,\pm i\}^n}
 \bigl(|\operatorname{Re}\widetilde Q_A(z)|
      +|\operatorname{Im}\widetilde Q_A(z)|\bigr).          \tag{6.14h}
\]
The tempting zero-loss coherent bound
\(L_{\rm cl}(A)\le\sqrt2\Phi(A)\) is false even for a global optimizer:
at order four, the signing having only the edge \(\{2,3\}\) negative has
\(\Phi(A)=m_4=4\), whereas reversing that clique produces the all-positive
signing and hence \(L_{\rm cl}(A)=6>4\sqrt2\).

*Proof.*  The energy of \({\cal K}(A,C)\) at \((x,y)\) is
\(Q_A(x)-Q_A(y)+x^TCy\).  Replacing \(x\) by \(-x\) fixes the internal
terms and negates the cross term; maximizing over this pair and using
\(\max(|a+b|,|a-b|)=|a|+|b|\) proves (6.14d).

Write every pair uniquely as \(x=s\), \(y=D_Ts\).  Collecting the
coefficient of \(s_is_j\) and the diagonal cross terms gives
\[
 Q_{{\cal K}(A,C)}(s,D_Ts)=2\bigl(Q_{G_T}(s)+h_T\bigr),
\]
which proves (6.14d2).  Replacing \(T\) by \(T^c\) fixes the first term
of \(G_T\) and negates the second term and \(h_T\).  Pairing these two
values proves (6.14d3).  The four-label table follows by substituting the
four ordered sign pairs in (6.14d1).  For \(C=R+D\), it gives
\(Q_{G_T}(s)=2\sum_{u\in T,v\notin T,u\to v}A_{uv}s_us_v\),
while \(h_T=0\) for \(C=R\); (6.14d2) and the diamond identity prove
(6.14d4).  Adding \(D\) changes every energy by at most \(n\), proving
(6.14d5).  Finally, as \(T,s\) vary,
\((s,D_Ts)\) ranges over all Boolean pairs, and averaging \(\|Cy\|_1\)
over \(y\) proves (6.14d9).

For (6.14f), write every two-vertex-cloud state uniquely as
\((x_i,y_i)=s_i(1,(-1)^{t_i})\) and put \(T=\{i:t_i=1\}\).  Direct
expansion gives
\[
 Q_{{\cal K}(A,C_0+D)}(x,y)
 =2Q_{B_T(A,C_0)}(s)+\sum_i d_i(-1)^{t_i}.          \tag{6.14i}
\]
The last term has magnitude at most \(n\); maximizing proves (6.14f), and
\(C_0=A\) gives (6.14g).  Finally write \(z_i=s_i\) off \(T\) and
\(z_i=is_i\) on \(T\).  If \(O,J,X\) are the energies off \(T\), on
\(T\), and across the cut, respectively, then
\(\widetilde Q_A(z)=O-J+iX\), while
\(Q_{A^{K_T}}(s)=O-J+X\).  Flipping all spins on \(T\) reverses \(X\)
and proves (6.14h).  For the displayed order-four optimizer, Walsh
orthogonality gives \({\bf E}Q_A^2=6\); since all values are even, a norm
below four is impossible, while direct evaluation gives norm four.  The
clique flip gives norm six. \(\square\)

Thus ordinary holomorphic complexification is not the missing lemma.  The
Cauchy--Riemann equations force its off-diagonal cross block to be \(C_0=A\)
(the diagonal contributes only the \(O(n)\) term above), and that coherent
choice has just failed.  The genuinely new surviving target in (6.14d) is a
noncoherent, \(A\)-dependent cross signing whose large bilinear values are
correlated with a deficit in \(|Q_A(x)-Q_A(y)|\).  This target is distinct
from the skew/directed-half-cut target of Proposition 6.5a and remains open.

The full two-variable calculus on a fixed coherent phase face is explicit.
Fix \(T\), a sign state \(s\), and write \(O,J,X\) for its energies on
\(T^c\), on \(T\), and across the cut.  Multilinearity shows that
\[
 p(r,t)=Or^2+Jt^2+Xrt
\]
has \(|p(r,t)|\le M\) throughout \([-1,1]^2\).  Conversely this square
bound holds if and only if
\[
\begin{aligned}
 |O+J|+|X|&\le M,\\
 |O-X^2/(4J)|&\le M
   &&\text{when }J\ne0\text{ and }|X|\le2|J|,\\
 |J-X^2/(4O)|&\le M
   &&\text{when }O\ne0\text{ and }|X|\le2|O|.
\end{aligned}                                               \tag{6.14j0}
\]
Indeed homogeneity moves every nonzero point to the boundary of the square;
on its four sides the only candidates are the endpoints and the displayed
quadratic vertices.  In particular, if the two corner states have opposite
extremal energies, so that
\(p(1,1)=M\) and \(p(1,-1)=-M\), then
\[
 O+J=0,\qquad X=M,\qquad |O|=|J|\le M/2.          \tag{6.14j1}
\]

Restriction to the unit circle gives the convenient consequence
\[
 \boxed{|O+J|+\sqrt{(O-J)^2+X^2}\le2M.}            \tag{6.14j}
\]
If \(k=|T|\), completing either fixed partial state and pairing a completion
with its negative gives the stronger marginal bounds
\[
 |O|\le M-\lfloor k/2\rfloor,
 \qquad
 |J|\le M-\lfloor(n-k)/2\rfloor.                  \tag{6.14k}
\]
Indeed every complete signing on \(r\) vertices has one-sided maximum at
least \(\lfloor r/2\rfloor\): expose the vertices sequentially and choose
each new spin to make its incident increment nonnegative; every odd-length
partial row contributes at least one.  Applying this also to the negative
signing and absorbing the cross-linear term by pairing opposite completions
proves (6.14k).

Neither statement supplies a leading-constant gain.  The scalar polynomial
\(M(r^2-t^2)\) is bounded by \(M\) on the square and attains equality in
(6.14j0)--(6.14j), while (6.14k) improves its two axial coefficients by only
\(O(n)=o(n^{3/2})\).  Thus scalar face calculus alone cannot prove
\(|O-J|+|X|\le\sqrt2M+o(n^{3/2})\); a successful argument must couple
different faces or use additional global-minimizer structure.  In the
directed split \(X=F+G\), even the complete criterion (6.14j0) contains no
information about the cancellation coordinate \(F-G\).

The skew term cannot be budgeted independently of the two quadratic
energies.  For every skew signing,
\[
 \max_{x,y}|x^TRy|=\max_x\|Rx\|_1
 \ge n\,\mathbb E|S_{n-1}|
 =\left(\sqrt{2/\pi}+o(1)\right)n^{3/2},           \tag{6.15}
\]
where \(S_k\) is a sum of \(k\) independent signs.  By contrast, combining
only \(|Q_A(x)+Q_A(y)|\le2M\) with (6.13) would require
\(\max|x^TRy|\le2(\sqrt2-1)M+o(n^{3/2})\).  The random-method upper bound
\(M\le(\sqrt{\log2}+o(1))n^{3/2}\) makes the latter leading constant
\(2(\sqrt2-1)\sqrt{\log2}<\sqrt{2/\pi}\).  Hence every uncoupled
internal/cross budget is impossible; statewise anticorrelation is
load-bearing.

The superficially natural disk strengthening
\[
 I_A(x,y)^2+C_R(x,y)^2\le M^2+o(n^3),\qquad
 I_A={Q_A(x)+Q_A(y)\over2},\quad C_R={x^TRy\over2}, \tag{6.16}
\]
is not a harmless replacement for the diamond (6.13): (6.15) shows that it
would itself prove
\(\liminf m_n/n^{3/2}\ge1/\sqrt{2\pi}\), much stronger than the known
\(1/\pi\) floor.  Its zero-error form is already false for an optimizer at
\(n=5\).  Every order-five signing has
\(\mathbb E_xQ_A(x)^2=10\), while all its energies are even, so
\(\Phi(A)\ge4\).  Take \(A=-1\) on a five-cycle and \(+1\) on its diagonals;
it attains \(\Phi(A)=4\), hence \(m_5=4\).  Its five positive maximizers,
whose negative coordinates are \(\{r,r+2\}\), form an invertible circulant
matrix \(V\).
Exact (6.16) would force \(VRV^T=0\), hence \(R=0\), impossible for a skew
signing.  The exact diamond, not the disk, is the live target.

**Proposition 6.5d (bivector cover gate and low-degree no-go).**  Identify
the upper-triangular entries of a skew signing with
\(r\in\{\pm1\}^N\), \(N=\binom n2\), and put
\[
 b_{xy,ij}={x_iy_j-x_jy_i\over2},\qquad
 d_{xy}=M-{|Q_A(x)+Q_A(y)|\over2}.
\]
Then
\[
 \boxed{K(A,R)-M
 =\max_{x,y}\bigl(|\langle b_{xy},r\rangle|-d_{xy}\bigr).} \tag{6.16a}
\]
For \(t\ge0\), let
\[
 H_{(x,y),r}=\mathbf1\{|\langle b_{xy},r\rangle|>d_{xy}+t\}.
\]
The exact nonlinear covering value is
\[
 \eta_t=\max_{\lambda\in\Delta(x,y)}\min_r\mathbb E_\lambda H
       =\min_{p\in\Delta(r)}\max_{x,y}\mathbb E_pH.        \tag{6.16b}
\]
It vanishes if and only if an orientation with \(K(A,R)\le M+t\)
exists.  If it is positive, \(\lambda\) is precisely a fractional cover of
all tournaments by bad state pairs.

The rows have the rigid Grassmann geometry
\[
 b_{xy}={x\wedge y\over2},\quad
 \|b_{xy}\|_2^2=d_H(x,y)(n-d_H(x,y)),                       \tag{6.16c}
\]
and satisfy all Pluecker relations.  Nevertheless the following natural
affine, covariance, and normalized single-row moment relaxations are
subcritical.  Already
\[
 \min_{r\in[-1,1]^N}\max_{x,y}
 \bigl(|\langle b_{xy},r\rangle|-d_{xy}\bigr)=0:            \tag{6.16c1}
\]
the upper bound is \(r=0\), and a diagonal pair \((z,z)\) at an
\(A\)-maximizer forces the reverse bound.  At
\(t=(\sqrt2-1)M\), put \(w_{xy}=t+d_{xy}\).  Then
\[
 \min_{\substack{X\succeq0\\\operatorname{diag}X=1}}
 \max_{x,y}{b_{xy}^TXb_{xy}\over w_{xy}^2}
 \le{\pi^2\over4(\sqrt2-1)^2(n-1)}<1qquad(n\ge16),         \tag{6.16d}
\]
and no witness distribution can have
\[
 \inf_r\mathbb E_\lambda
 \left({|\langle b,r\rangle|\over w}\right)^{2p}>1         \tag{6.16e}
\]
whenever
\[
 p<{2(\sqrt2-1)^2\over\pi^2}(n-1)=0.0347679\ldots(n-1).    \tag{6.16f}
\]
Moreover
\[
 \eta_t\le2\exp\left(-{2(\sqrt2-1)^2\over\pi^2}(n-1)\right), \tag{6.16g}
\]
whereas failure implies only \(\eta_t\ge4^{-n}\).  Thus the affine
relaxation, the displayed elliptope, and normalized certificates using one
even moment of each row are rigorously subcritical.  This does not exclude a
different fixed-degree SOS encoding or higher-degree cross-row products.
Proposition 6.5j below rules out the full degree-four preordering in the
natural squared-row encoding, even with pairwise row products and every
identity among the instantiated Pluecker rows.  A successful continuation
must strengthen those relaxations, use the \(A\)-dependent widths at higher
degree, or round (6.16b) directly.

*Proof.*  One has
\(\langle b_{xy},r\rangle=x^TRy/2\), which proves (6.16a).
Equation (6.16b) is finite minimax duality for the zero--one matrix \(H\);
zero value forces every support atom on the right to be a successful
orientation.  Expanding the exterior product proves (6.16c).  Proposition
5.2 gives
\(t\ge(\sqrt2-1)n\sqrt{n-1}/\pi\), while
\(\|b_{xy}\|_2^2\le n^2/4\); the feasible covariance \(X=I\) proves
(6.16d).  Khintchine's inequality bounds the left moment before averaging
over \(\lambda\) by
\[
 (2p-1)!!\left({\pi^2\over4(\sqrt2-1)^2(n-1)}\right)^p.
\]
Using \((2p-1)!!\le(2p)^p\) proves (6.16e)--(6.16f).
Hoeffding proves (6.16g) using the uniform random orientation.  If every
orientation has a bad pair, the uniform distribution on the \(4^n\)
ordered pairs gives \(\eta_t\ge4^{-n}\). \(\square\)

**Proposition 6.5e (signed-regular arcsine rigidity at the outgoing-half
threshold).**  Let \(H\) be a symmetric zero-diagonal
\(\{0,\pm1\}\)-matrix of order \(N\), with \(d\)-regular nonzero support,
\(d\ge2\).  Then
\[
 \Phi(H)\ge {Nd\over\pi}\arcsin{1\over\sqrt d}
 +{1\over4\pi d^{5/2}(1-1/d)^{3/2}}
   \sum_{\substack{i<j\\H_{ij}\ne0}}(H^2)_{ij}^2.          \tag{6.16h}
\]
In particular, let \(S\) be a tournament, \(R=A\circ S\), and put
\[
 K_0=\begin{pmatrix}A&R\\-R&-A\end{pmatrix},\qquad d=2n-2.
\]
Then \(K_0\) is supported on \(K_{2n}\) minus a perfect matching,
\(\Phi(K_0)=4D_{\to}(A,S)\), and
\[
 \boxed{
 D_{\to}(A,S)\ge
 {n(n-1)\over\pi}\arcsin{1\over\sqrt{2n-2}}
 +{\Sigma(A,R)\over16\pi d^{5/2}(1-1/d)^{3/2}},}           \tag{6.16i}
\]
where
\[
 \Sigma(A,R)=\sum_{i\ne j}
 \left((A^2-R^2)_{ij}^2+(AR-RA)_{ij}^2\right).             \tag{6.16j}
\]
The leading term in (6.16i) is
\((1/(\sqrt2\pi)+o(1))n^{3/2}\), exactly the target
\(M/\sqrt2\) when \(M\sim n^{3/2}/\pi\).  More quantitatively, if
\(M=\alpha_n n^{3/2}\) and an orientation has
\(D_{\to}\le M/\sqrt2+\epsilon_n n^{3/2}\), then
\[
 {\Sigma(A,R)\over n^4}
 \le64\pi(\alpha_n-1/\pi)
   +16\pi2^{5/2}\epsilon_n+o(1).                           \tag{6.16k}
\]
Thus success along a subsequence with \(\alpha_n\to1/\pi\) and
\(\epsilon_n\to0\) forces \(\Sigma(A,R)=o(n^4)\): the skew mate must
approximately commute with \(A\) and have the same square off the diagonal.
Finally, if
\[
 U_n=\min_{\Phi(A)=m_n}\min_{S\ {\rm tournament}}D_{\to}(A,S),
\]
then every fixed \(n_0\ge2\) satisfies
\[
 \limsup_{j\to\infty}{U_{2^jn_0}\over m_{2^jn_0}}
 \ge{1\over\sqrt2}.                                       \tag{6.16k1}
\]
Thus the \(1/\sqrt2\) simultaneous-paving constant is forced both by the
arcsine geometry at the lower floor and by the dyadic recurrence globally.

*Proof.*  Let \(g\) be standard Gaussian and set
\(z^\pm=(I\pm H/\sqrt d)g\), \(x_i^\pm=\operatorname{sgn}z_i^\pm\).
For a support edge \(e=\{i,j\}\), the two correlations are
\[
 r_e^\pm={(H^2)_{ij}\over2d}\pm{H_{ij}\over\sqrt d}.
\]
With \(u_e=H_{ij}(H^2)_{ij}/(2d)\), \(v=d^{-1/2}\), the Gaussian sign
law and oddness of arcsine give
\[
 \Phi(H)\ge{1\over\pi}\sum_e
 \Delta(u_e,v),\qquad
 \Delta(u,v)=\arcsin(u+v)-\arcsin(u-v).                    \tag{6.16l}
\]
Indeed the right side is half the difference of the two expected signed
energies.  Now
\(\Delta(u,v)=\int_{u-v}^{u+v}(1-t^2)^{-1/2}\,dt\).
Its second derivative is the integral over the same interval of
\((1+2t^2)(1-t^2)^{-5/2}\), an even function increasing in \(|t|\).
The interval is minimized when centered at zero, so
\[
 \Delta(u,v)\ge2\arcsin v+{v\over(1-v^2)^{3/2}}u^2.
\]
Summing over the \(Nd/2\) edges proves (6.16h).

For \(K_0\), the four-label identity proves
\(\Phi(K_0)=4D_{\to}\), and block multiplication gives
\[
 K_0^2=
 \begin{pmatrix}A^2-R^2&AR-RA\\AR-RA&A^2-R^2\end{pmatrix}.
\]
Summing only over its support edges gives exactly (6.16j), proving
(6.16i).  Rearrangement and expansion of the arcsine yield (6.16k).
Filling the missing perfect matching of \(K_0\) changes its norm by at most
\(n\), so \(m_{2n}\le4U_n+n\).  If (6.16k1) failed on an eventual dyadic
tail with a constant \(c<1/\sqrt2\), then
\(\alpha_{2n}\le\sqrt2c\,\alpha_n+o(1)\), forcing that positive sequence
to zero and contradicting Proposition 5.2.
\(\square\)

**Proposition 6.5f (finite-anchor signature-cell shielding).**  Let
\(x^{(1)},\ldots,x^{(k)}\in\{\pm1\}^n\), deduplicated modulo global sign.
Partition \([n]\) by the coordinate signatures
\[
 \bigl(x_i^{(1)}x_i^{(2)},\ldots,x_i^{(1)}x_i^{(k)}\bigr),
\]
and let \(L\le\min(n,2^{k-1})\) be the number of nonempty cells.  There is
a skew signing \(R\) such that, simultaneously for every anchor,
\[
 \boxed{\|Rx^{(a)}\|_1\le Ln,\qquad
 |(x^{(a)})^TRy|\le Ln\quad\hbox{for every }y.}             \tag{6.16m}
\]

Let \(M=\Phi(A)\), \(S=A\circ R\),
\(\Delta_A(x,y)=|Q_A(x)-Q_A(y)|\), and
\(\Gamma(R)=\max_{x,y}|x^TRy|\).  If the anchor set is a vertex cover of
\[
 \bigl\{\{x,y\}:\Delta_A(x,y)>2\sqrt2M-\Gamma(R)\bigr\}    \tag{6.16n}
\]
and \(Ln\le2(\sqrt2-1)M\), then
\[
 D_{\to}(A,S)\le M/\sqrt2.                                 \tag{6.16o}
\]
The size condition is automatic for every \(A\) when
\[
 L\le {2(\sqrt2-1)(n-1)\over\pi}
        \arcsin{1\over\sqrt{n-1}}
 =\left({2(\sqrt2-1)\over\pi}+o(1)\right)\sqrt n.          \tag{6.16p}
\]
Thus even an arbitrary anchor family with roughly
\(k\le\tfrac12\log_2n-O(1)\) can be removed at the critical margin.  A
finite extremizer list is not the obstruction; the unresolved issue is the
uncontrolled complement and its growing signature complexity.

*Proof.*  Gauge \(x^{(1)}\) to \({\bf1}\).  Inside each odd signature cell
use a regular tournament, and inside each even cell delete one vertex from a
regular tournament of the next odd order; the internal row sums are then
zero or \(\pm1\).  Between two locally indexed cells put the alternating
block \(T_{i_rj_s}=(-1)^{r+s}\), and use its negative transpose in the
reverse block.  Every row sum into each cell has magnitude at most one.
Every gauged anchor is constant on each of the \(L\) cells, so every
coordinate of its image under \(T\) has magnitude at most \(L\).  Undoing
the gauge and summing proves (6.16m).

Every anchor-incident pair has
\(\Delta_A(x,y)+|x^TRy|\le2M+Ln\le2\sqrt2M\).  Every
uncovered pair in (6.16n) has the same upper bound by the definition of
\(\Gamma\).  Identity (6.14d4) proves (6.16o).  Finally (6.16h), applied
to the complete signing \(A\) and with its nonnegative correction dropped,
gives
\(M\ge n(n-1)\arcsin(1/\sqrt{n-1})/\pi\), proving (6.16p).
\(\square\)

For one anchor the exact optimum in (6.16m) is zero for odd \(n\) and
\(n\) for even \(n\).  Two projectively distinct Boolean anchors cannot
both lie in the exact kernel of a skew signing: after gauging and splitting
by their relative sign, the within-cell and cross-cell row sums would both
have to vanish, and their parities are incompatible.  This parity fact does
not give a uniform \(\Omega(n)\) two-anchor floor at odd order.

**Proposition 6.5g (random approximate mate and generic spectral-bridge
no-go).**  For a fixed order-\(n\) signing \(A\), \(n\ge3\), a uniform skew
signing \(R\) satisfies
the exact identity
\[
 \mathbb E_R\Sigma(A,R)
 =\sum_{i\ne j}(A^2_{ij})^2+3n(n-1)(n-2).                  \tag{6.16q}
\]
Moreover, with
\[
 L_n={n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}},
\]
Proposition 6.5e gives
\[
 \sum_{i\ne j}(A^2_{ij})^2
 \le8\pi(n-1)(n-2)^{3/2}\bigl(\Phi(A)-L_n\bigr).          \tag{6.16r}
\]
Consequently, if \(A\) is optimal and
\(M=m_n=\alpha_n n^{3/2}\), there is a deterministic skew signing with
\[
 {\Sigma(A,R)\over n^4}
 \le8\pi(\alpha_n-1/\pi)+O(1/n).                           \tag{6.16s}
\]
In particular, the approximate-mate condition forced by Proposition 6.5e
near \(\alpha_n=1/\pi\) is automatically satisfiable; it is not itself the
missing upper construction.

This cannot be repaired by the generic spectral conversion.  Put
\(C=AR-RA\), \(T(R)=\sum_iC_{ii}^2\), and \(d=2n-2\).  Exactly
\[
 \|K_0^2-dI\|_F^2=2(\Sigma+T),\qquad
 \operatorname{tr}(K_0^2-dI)=0,                            \tag{6.16t}
\]
so
\[
 D_{\to}(A,S)\le{n\over4}
 \sqrt{d+\sqrt{2n-1\over2n}\sqrt{2(\Sigma+T)}}.            \tag{6.16u}
\]
Even the ideal value \(\Sigma+T=0\) yields only
\((1/(2\sqrt2)+o(1))n^{3/2}\), a factor \(\pi/2\) above the
lower-floor target.  More generally, along any subsequence on which
\(\alpha_n\to\alpha\) and \((\Sigma+T)/n^2\to s\), (6.16u) can close the
desired bound only when
\[
 2+\sqrt{2s}\le8\alpha^2.                                  \tag{6.16v}
\]
The asymptotic upper bound gives \(\alpha\le1/2\), so this forces
\(\alpha=1/2\) and \(s=0\).  Thus the displayed generic
trace/Frobenius-to-operator conversion cannot prove the original
amplification in the unknown interior range.  This does not exclude an
upper bound that exploits additional special structure of the same square
and commutator matrices.

*Proof.*  For \(i\ne j\), the \(n-2\) summands of \((R^2)_{ij}\) are
independent centered signs, whereas
\[
 (AR-RA)_{ij}=\sum_{k\ne i,j}A_{ik}A_{kj}(S_{kj}-S_{ik})
\]
has variance \(2(n-2)\).  Hence
\(\mathbb E(A^2-R^2)_{ij}^2=(A^2_{ij})^2+n-2\), proving
(6.16q).  Equation (6.16r) is (6.16h) applied to \(A\), after doubling the
sum over unordered edges.  Conditional expectation derandomizes (6.16q),
and expansion gives (6.16s).

Block multiplication proves (6.16t).  A traceless symmetric matrix of order
\(2n\) has largest eigenvalue at most
\(\sqrt{(2n-1)/(2n)}\) times its Frobenius norm.  Apply this to
\(K_0^2-dI\), then use
\(\Phi(K_0)\le n\|K_0\|_{\rm op}\) and
\(\Phi(K_0)=4D_{\to}\), to obtain (6.16u).  Its normalized limit is
(6.16v). \(\square\)

**Proposition 6.5h (exact outgoing-half random criterion and first-moment
obstruction).**  In the all-directed setting of Proposition 6.5c, put
\[
 C_{T,s}=\sum_{u\in T,\ v\notin T}A_{uv}s_us_v,\qquad
 K_{T,s}=\sum_{u\in T,\ v\notin T}R_{uv}s_us_v.
\]
Then
\[
 \boxed{
 D_\to(A,S)=\max_{T,s}{|C_{T,s}|+|K_{T,s}|\over2},
 \qquad S=A\circ R.
 }                                                           \tag{6.16w}
\]
In particular \(|C_{T,s}|\le M=\Phi(A)\), and the zero-error target
\(D_\to(A,S)\le M/\sqrt2\) is exactly
\[
 |K_{T,s}|\le\sqrt2M-|C_{T,s}|
 \quad\hbox{for every }(T,s).                                 \tag{6.16x}
\]

Choose the upper-triangular entries of \(R\) independently and uniformly.
For \(m\ge1\), let
\[
 {\rm Tail}_m(t)
 =2^{-m}\sum_{\substack{0\le j\le m\\|2j-m|>t}}{m\choose j}.
\]
Modulo global spin reversal and \(T\leftrightarrow T^c\), the exact
first-moment sufficient criterion for (6.16x) is
\[
 {\cal F}(A):=
 \sum_{[(T,s)]}{\rm Tail}_{|T|(n-|T|)}
 \bigl(\sqrt2M-|C_{T,s}|\bigr)<1.                             \tag{6.16y}
\]
Equivalently, if \(d=d_H(x,y)\),
\[
 {\cal F}(A)={1\over4}
 \sum_{\substack{x,y\\0<d<n}}
 {\rm Tail}_{d(n-d)}
 \left(\sqrt2M-\tfrac12|Q_A(x)-Q_A(y)|\right).                 \tag{6.16z}
\]
Hoeffding gives the simpler sufficient upper certificate
\[
 \sum_{[(T,s)]}2\exp\left(
 -{(\sqrt2M-|C_{T,s}|)^2\over2|T|(n-|T|)}\right)<1.
\]

The exact quantity in (6.16y), not merely this Hoeffding upper bound, is
already exponentially large for optimal \(A\).  Writing
\(M=\alpha_n n^{3/2}\), one has
\[
 {\cal F}(A)\ge
 \exp\bigl((\log4-4\alpha_n^2-o(1))n\bigr).                    \tag{6.16aa}
\]
Thus \(\limsup\alpha_n\le1/2\) gives
\[
 {\cal F}(A)\ge
 \exp\bigl((\log4-1-o(1))n\bigr),
 \qquad \log4-1=0.386294\ldots>0.                             \tag{6.16ab}
\]
The literal independent-random-orientation first-moment/union-bound proof
therefore cannot establish the outgoing-half target.  This does **not**
prove that random orientations fail: the bad events are highly dependent,
and a correlated process estimate remains open.

*Proof.*  The outgoing and incoming signed halves are
\[
 {C_{T,s}+K_{T,s}\over2},\qquad {C_{T,s}-K_{T,s}\over2}.
\]
Complementing \(T\) interchanges them, and
\(\max(|a+b|,|a-b|)=|a|+|b|\), proving (6.16w).
If \(D_T\) changes signs on \(T\), then, for \(x=s,\ y=D_Ts\),
\[
 Q_A(x)-Q_A(y)=2C_{T,s},\qquad x^TRy=2K_{T,s}.
\]
This also proves \(|C_{T,s}|\le M\), (6.16x), and the pair
parametrization in (6.16z).  The quotient factor is four because
independently negating \(x\) and \(y\) leaves each constraint unchanged.
Swapping \(x,y\) is a further twofold redundancy that is deliberately not
quotiented here; removing it halves \({\cal F}(A)\) and changes none of the
exponential conclusions below.

For fixed \((T,s)\), \(K_{T,s}\) is exactly a sum of
\(|T|(n-|T|)\) independent signs.  Hence (6.16y) is the expected number of
violated constraint orbits; expectation below one supplies a realization
with none.  Hoeffding gives the displayed relaxation.

For the obstruction, take \(d=\lfloor n/2\rfloor\) and
\(h=d(n-d)=\lfloor n^2/4\rfloor\).  There are
\(2^n{n\choose d}\) ordered pairs at this distance, and every threshold in
(6.16z) is at most \(\sqrt2M\).  Therefore
\[
 {\cal F}(A)\ge {2^n\over4}{n\choose d}
 {\rm Tail}_h(\sqrt2M).                                      \tag{6.16ac}
\]
Uniformly for \(m\ge cn^2\) and \(0\le t\le Cn^{3/2}\),
\[
 {\mathbb P}(|S_m|>t)\ge
 \exp\left(-{t^2\over2m}-O_{c,C}(\log n)\right).              \tag{6.16ad}
\]
Indeed, choose the least \(k>t\) having the parity of \(m\).  A single
binomial atom, uniform Stirling bounds, and
\[
 I(a)={(1+a)\log(1+a)+(1-a)\log(1-a)\over2}
 ={a^2\over2}+O(a^4)
\]
give
\[
 {\mathbb P}(S_m=k)\ge c_0m^{-1/2}e^{-mI(k/m)}
 =\exp\left(-{t^2\over2m}-O(\log n)\right),
\]
because \(k/m=O(n^{-1/2})\) and \(k^4/m^3=O(1)\).
Apply this to (6.16ac), and use the central-binomial Stirling estimate, to
obtain (6.16aa); (6.16ab) follows from the known optimal-signing upper
scale. \(\square\)

Even the fictitious favorable assignment \(C_{T,s}=0\) on every central
constraint would require
\(\alpha_n>\sqrt{\log2/2}=0.588705\ldots\) for this first moment to decay,
already above \(1/2\).  At \(\alpha_n\to1/\pi\), its exponent is
\(\log4-4/\pi^2=0.981009\ldots\).  These constants close only the literal
first-moment route, not the underlying random-orientation distribution.

**Proposition 6.5i (Gaussian saturation and central two-half saddle).**
Let \(H\) be a symmetric zero-diagonal \(\{0,\pm1\}\)-matrix of order \(N\)
with \(d\)-regular support, \(d\ge2\), and put \(P=\Phi(H)\).  For standard
Gaussian \(g\), set

\[
 X^\pm=\operatorname{sgn}\bigl((I\pm H/\sqrt d)g\bigr).
\]

For a support edge \(e=\{i,j\}\), put

\[
 u_e={H_{ij}(H^2)_{ij}\over2d},\quad v=d^{-1/2},\quad
 \Delta(u,v)=\arcsin(u+v)-\arcsin(u-v),\quad
 \mathcal L_G(H)={1\over\pi}\sum_e\Delta(u_e,v).
\]

Then the arcsine argument has the exact saturation-gap form

\[
 \boxed{
 \mathbb E\big[(P-Q_H(X^+))+(P+Q_H(X^-))\big]
 =2\big(P-\mathcal L_G(H)\big).}                         \tag{6.16ae}
\]

Moreover, if \(\mathcal O=\sum_iX_i^+X_i^-\), then

\[
 \mathbb E{\cal O}=0,\qquad
 \boxed{\operatorname{Var}{\cal O}
 \le3N+{2\over d^2}\sum_{i<j}(H^2_{ij})^2.}              \tag{6.16af}
\]

In particular, let \(A\) be an optimal complete signing, \(M=\Phi(A)\),
and \(L_n=n(n-1)\arcsin(1/\sqrt{n-1})/\pi\).  Its two coupled Gaussian
outputs satisfy

\[
\begin{aligned}
 \mathbb E[(M-Q_A(X^+))+(M+Q_A(X^-))]&\le2(M-L_n),\\
 \operatorname{Var}{\cal O}
 &\le3n+{8\pi(n-2)^{3/2}\over n-1}(M-L_n).
\end{aligned}                                             \tag{6.16ag}
\]

Consequently, along every subsequence on which
\(M/n^{3/2}\to1/\pi\),

\[
 \boxed{
 Q_A(X^+)=M-o_{\Pr}(n^{3/2}),\quad
 Q_A(X^-)=-M+o_{\Pr}(n^{3/2}),\quad
 d_H(X^+,X^-)=n/2+o_{\Pr}(n).}                            \tag{6.16ah}
\]

Thus the Hamming-central opposite-energy layer in the outgoing-half
problem is forced, rather than excluded, when the universal lower floor is
approached.

There is a corresponding necessary shape for any sharp outgoing-half
orientation.  Let \(R=A\circ S\), form

\[
 K_0=\begin{pmatrix}A&R\\-R&-A\end{pmatrix},
\]

and suppose, on the same subsequence,
\(D_{\to}(A,S)\le M/\sqrt2+o(n^{3/2})\).  Apply the preceding coupling to
\(K_0\), and call its outputs \(Y^+,Y^-\).  Then

\[
 Q_{K_0}(Y^+)=\Phi(K_0)-o_{\Pr}(n^{3/2}),\quad
 Q_{K_0}(Y^-)=-\Phi(K_0)+o_{\Pr}(n^{3/2}),\quad
 d_H(Y^+,Y^-)=n+o_{\Pr}(n).                               \tag{6.16ai}
\]

Writing \(u=(Y^++Y^-)/2\), \(v=(Y^+-Y^-)/2\),
\(O=Q_{K_0}(u)\), \(J=Q_{K_0}(v)\), and \(X=u^TK_0v\), one obtains the
calculator's balanced saddle

\[
 \boxed{
 O+J=o_{\Pr}(n^{3/2}),\quad
 X=\Phi(K_0)-o_{\Pr}(n^{3/2}),\quad
 |O|,|J|\le\Phi(K_0)/2+o_{\Pr}(n^{3/2}).}                 \tag{6.16aj}
\]

This is a necessary-structure theorem, not an orientation construction,
and it does not close the multiplier-two ray.

*Proof.*  The covariance computation in Proposition 6.5e gives exactly

\[
 {1\over2}\big(\mathbb E Q_H(X^+)-\mathbb E Q_H(X^-)\big)
 =\mathcal L_G(H).
\]

Both deficits in (6.16ae) are nonnegative, which proves that identity and,
by Markov's inequality, their simultaneous concentration whenever
\(P-\mathcal L_G(H)=o(P)\).

For (6.16af), normalize the two-dimensional Gaussian vector
\(\bigl(((I+H/\sqrt d)g)_i,((I-H/\sqrt d)g)_i\bigr)\) by \(\sqrt2\).
Its coordinates are independent standard Gaussians because
\((H^2)_{ii}=d\).  For \(i\ne j\), the cross-correlation matrix of the
normalized pairs at \(i,j\)
is

\[
 \begin{pmatrix}\widetilde u+\widetilde v&-\widetilde u\\
 -\widetilde u&\widetilde u-\widetilde v\end{pmatrix},qquad
 \widetilde u={(H^2)_{ij}\over2d},\quad
 \widetilde v={H_{ij}\over\sqrt d}.
\]

The function \((s,t)\mapsto\operatorname{sgn}s\operatorname{sgn}t\) has
Gaussian Hermite rank two.  Singular-value decomposition of the displayed
cross-correlation matrix and the Hermite expansion therefore bound the
covariance of \(X_i^+X_i^-\) and \(X_j^+X_j^-\) by the squared operator
norm, hence by

\[
 {2H_{ij}^2\over d}+{(H^2)_{ij}^2\over d^2}.
\]

Summing this bound and using
\(\sum_{i<j}H_{ij}^2=Nd/2\) proves (6.16af).

For \(H=A\), Proposition 6.5e gives
\(\mathcal L_G(A)\ge L_n\), and (6.16r), substituted in (6.16af), proves
(6.16ag).  Markov and Chebyshev prove (6.16ah).

For \(H=K_0\), Proposition 6.5e and the assumed upper bound squeeze
\(\Phi(K_0)-\mathcal L_G(K_0)=o(n^{3/2})\) and force
\(\Sigma(A,R)=o(n^4)\).  Direct block multiplication gives

\[
 \sum_{i<j}(K_0^2)_{ij}^2=\Sigma(A,R)+T(R),\qquad
 T(R)=\sum_i(AR-RA)_{ii}^2\le4n(n-1)^2.
\]

Thus (6.16af) has variance \(o(n^2)\), proving (6.16ai).  The two corner
energies give \(O+J=o_{\Pr}(n^{3/2})\) and
\(X=\Phi(K_0)-o_{\Pr}(n^{3/2})\).  Finally
\(p(r,t)=Or^2+Jt^2+Xrt\) is bounded by \(\Phi(K_0)\) on the full square.
After normalization, compactness and the exact opposite-corner calculus
(6.14j1) give the last two bounds in (6.16aj). \(\square\)

One tempting converse is false.  If \(A^2=(n-1)I\), then the two Gaussian
vectors above are independent as vectors and their sign covariance
matrices are \(I\pm cA\), where
\(c=(2/\pi)\arcsin(1/\sqrt{n-1})\).  For every skew signing \(R\),

\[
 \mathbb E(X^{+T}RX^-)^2
 =\operatorname{tr}\big((I+cA)R(I-cA)R^T\big)
 \le(1+c\sqrt{n-1})^2n(n-1)=O(n^2).                    \tag{6.16ak}
\]

So the canonical central opposite-energy coupling itself has subcritical
skew interaction for every orientation.  Uniform control of the whole
diffuse layer, not one Gaussian pair, remains the missing implication.

**Proposition 6.5j (degree-four squared-row preordering no-go).**  In the
bivector formulation of Proposition 6.5d, set
\[
 t=(\sqrt2-1)M,\qquad w_s=t+d_s,\qquad
 g_s(r)=w_s^2-\langle b_s,r\rangle^2,
 \quad s=(x,y).                                             \tag{6.16al}
\]
The exact critical orientation system is
\(r_e^2=1\) and \(g_s(r)\ge0\) for every row \(s\).  Consider its full
degree-four preordering in this squared-row encoding.  A refutation would
have the form
\[
 -1=\sigma_0+\sum_s\sigma_sg_s
       +\sum_{s\le u}\lambda_{su}g_sg_u
       +\sum_e h_e(r)(r_e^2-1),                             \tag{6.16am}
\]
where \(\sigma_0\) is an SOS polynomial of degree at most four, every
\(\sigma_s\) is a sum of squares of affine polynomials,
\(\lambda_{su}\ge0\), and every displayed term has degree at most four.
For every \(n\ge45\), no identity (6.16am) exists, whether or not the exact
orientation system is feasible.

More precisely, uniform expectation \(\mathbb E_0\) on
\(r\in\{\pm1\}^{\binom n2}\) is a feasible degree-four pseudomodel.  If
\(b\in\{0,\pm1\}^{\binom n2}\),
\(m=\|b\|_2^2\), \(\delta=w^2-m\), and
\(g=w^2-\langle b,r\rangle^2\), then the localizing form
\(\mathbb E_0[gq^2]\), for affine \(q\), has zero constant--linear block,
constant eigenvalue \(\delta\), inactive-coordinate block \(\delta I\),
and active-coordinate block
\[
 (\delta+2)I-2bb^T.
\]
Consequently
\[
 \mathbb E_0[gq^2]\ge0\quad\hbox{for every affine }q
 \quad\Longleftrightarrow\quad w^2\ge3m-2.                 \tag{6.16an}
\]
For two rows \(b,c\), write
\(\delta_b=w_b^2-\|b\|_2^2\),
\(\delta_c=w_c^2-\|c\|_2^2\), and
\(h=|\operatorname{supp}b\cap\operatorname{supp}c|\).  Their exact
cross-row fourth moment is
\[
 \boxed{\mathbb E_0[g_bg_c]
 =\delta_b\delta_c+2\langle b,c\rangle^2-2h.}              \tag{6.16ao}
\]
Thus this pseudomodel retains, rather than discards, the determinant pairing
of decomposable bivectors and their exact cut-support intersection.

*Proof.*  Expanding the uniform moments gives
\(\mathbb E_0g=\delta\), \(\mathbb E_0[gr_e]=0\), and, for
distinct coordinates, \(\mathbb E_0[gr_er_f]=-2b_eb_f\) when both are
active and zero otherwise.  This gives the displayed blocks.  The active
eigenvalues are \(\delta+2\) perpendicular to \(b\) and
\(w^2-3m+2\) along \(b\), proving (6.16an).

The Rademacher fourth-moment identity gives
\[
 \mathbb E_0[\langle b,r\rangle^2\langle c,r\rangle^2]
 =\|b\|_2^2\|c\|_2^2+2\langle b,c\rangle^2-2h,
\]
which is (6.16ao).  For the actual rows
\(b_{xy}=x\wedge y/2\),
\[
 \|b_{xy}\|_2^2=d_H(x,y)(n-d_H(x,y))
 \le m_*:=\lfloor n^2/4\rfloor.                            \tag{6.16ap}
\]
Proposition 5.2 implies
\[
 t^2\ge {3-2\sqrt2\over\pi^2}n^2(n-1)>3m_*
 \qquad(n\ge45).                                           \tag{6.16aq}
\]
For the endpoint inequality, it is enough to use
\(\sqrt2<5657/4000\) and \(\pi^2<10\), which give
\(44(3-2\sqrt2)/\pi^2>44(343/20000)>3/4\).
Since \(w_s\ge t\), (6.16an) makes every one-row affine localizer
positive semidefinite.  Also \(\delta_s>2m_*\), so (6.16ao) gives
\(\mathbb E_0[g_sg_u]>4m_*^2-2m_*>0\).

Finally, \(\mathbb E_0\) is an actual Boolean-cube expectation: it is
nonnegative on \(\sigma_0\) and annihilates the Boolean ideal.  It is
nonnegative on every remaining term of (6.16am) by (6.16an)--(6.16ao),
contradicting \(\mathbb E_0[-1]=-1\).  All Pluecker identities are respected
because the \(b_s\) are the instantiated decomposable rows; using such an
identity only rewrites the same polynomial before applying
\(\mathbb E_0\). \(\square\)

This is a method barrier, not an orientation construction.  It rules out
exactly (6.16am), including its first cross-row products.  Proposition 6.5l
below extends the same uniform-expectation obstruction to every fixed raw
degree at sufficiently large order.  Neither result rules out a different
lifted encoding or direct nonlinear rounding.

**Proposition 6.5k (weighted multi-anchor integral rounding).**  Let
\(x^{(1)},\ldots,x^{(k)}\in\{\pm1\}^n\), choose \(q_a>0\) with
\(\sum_aq_a\le1/2\), and put
\[
 \rho_a=n\sqrt{2/\pi}+\sqrt{2n\log(1/q_a)}.                \tag{6.16ar}
\]
For arbitrary budgets \(B_a>0\), the capacity condition
\[
 50\sum_{a=1}^k{\rho_a^2\over B_a^2}\le1                 \tag{6.16as}
\]
implies that there is a skew signing \(R\) such that
\[
 \|Rx^{(a)}\|_1\le B_a\qquad(1\le a\le k).                \tag{6.16at}
\]
This is an integral signing of every edge, with no fractional remainder.

Now let \(A\) be optimal and \(M=\Phi(A)\).  Define
\[
 B_A(x)=(2\sqrt2-1)M-|Q_A(x)|.                             \tag{6.16au}
\]
If (6.16as) holds with \(B_a=B_A(x^{(a)})\), then the one \(R\)
furnished above satisfies, for every anchor and every Boolean \(y\),
\[
 |Q_A(x^{(a)})-Q_A(y)|+|(x^{(a)})^TRy|\le2\sqrt2M.        \tag{6.16av}
\]
In particular, taking \(q_a=1/(2k)\) gives the uniform bound
\[
 \|Rx^{(a)}\|_1\le
 10n\sqrt{k/\pi}+10\sqrt{kn\log(2k)}.                     \tag{6.16aw}
\]
Consequently, for every fixed \(0<\delta<c_*\) and all sufficiently large
\(n\), an arbitrary family of
\[
 k\le(c_*-\delta)n,
 \qquad c_*={3-2\sqrt2\over25\pi}
       =0.0021845336957706\ldots                            \tag{6.16ax}
\]
can be shielded on every anchor-incident pair.  This improves the
worst-case arbitrary-anchor capacity of Proposition 6.5f from logarithmic
to a positive linear fraction.

*Proof.*  For an edge \(i<j\), form a vector with \(k\) blocks in
\(\mathbb R^{nk}\), whose \(a\)-th block is
\[
 v_{ij}^{(a)}=\lambda_a
       (x_j^{(a)}e_i-x_i^{(a)}e_j),
 \qquad \|v_{ij}\|_2^2=2\sum_a\lambda_a^2=:L^2.           \tag{6.16ay}
\]
For standard Gaussian \(g\in\mathbb R^n\), the function
\(g\mapsto\|g\|_1\) is \(\sqrt n\)-Lipschitz and has mean
\(n\sqrt{2/\pi}\).  Hence Gaussian concentration and (6.16ar) show that
the symmetric convex body
\[
 {\cal C}=\{(z^{(1)},\ldots,z^{(k)}):
                 \|z^{(a)}\|_1\le\rho_a\ \hbox{for every }a\}
\]
has Gaussian measure at least \(1-\sum_aq_a\ge1/2\).
Banaszczyk's vector-balancing theorem, in its unit-vector \(5{\cal C}\)
normalization, supplies signs \(r_{ij}\in\{\pm1\}\) such that
\[
 \sum_{i<j}r_{ij}{v_{ij}\over L}\in5{\cal C}.             \tag{6.16az}
\]
Define \(R_{ij}=r_{ij}\), \(R_{ji}=-r_{ij}\).  The \(a\)-th block in the
unnormalized sum is exactly \(\lambda_aRx^{(a)}\), and therefore
\[
 \|Rx^{(a)}\|_1\le{5L\rho_a\over\lambda_a}.
\]
Taking \(\lambda_a=\rho_a/B_a\), condition (6.16as) says \(5L\le1\),
which proves (6.16at).  This standard \(5{\cal C}\) theorem is the only
deep external input.

Since \(B_A(x)>0\) and
\(B_A(x)\ge2(\sqrt2-1)M\), (6.16at) gives
\[
\begin{aligned}
 |Q_A(x^{(a)})-Q_A(y)|+|(x^{(a)})^TRy|
 &\le |Q_A(x^{(a)})|+M+\|Rx^{(a)}\|_1\\
 &\le2\sqrt2M,
\end{aligned}
\]
proving (6.16av).  Substitution of \(q_a=1/(2k)\) in (6.16ar)--(6.16at)
gives (6.16aw).  Finally Proposition 5.2 gives
\[
 M\ge {n(n-1)\over\pi}\arcsin{1\over\sqrt{n-1}},
\]
and comparison of the leading \(n^{3/2}\) terms in (6.16aw) with
\(2(\sqrt2-1)M\) gives exactly (6.16ax). \(\square\)

The proposition does **not** control every pair.  For its resulting \(R\),
put \(S=A\circ R\),
\[
 \Gamma(R)=\max_{x,y}|x^TRy|,
 \qquad E_R=\{\{x,y\}:|Q_A(x)-Q_A(y)|
                    >2\sqrt2M-\Gamma(R)\}.                \tag{6.16ba}
\]
If the chosen anchors are a vertex cover of this same \(E_R\), then
(6.16av) handles its incident pairs and the definition of \(E_R\) handles
all others, so Proposition 6.5c gives \(D_{\to}(A,S)\le M/\sqrt2\).
This is a conditional statement about the same resulting \(R\): the theorem
does not bound \(\Gamma(R)\), and choosing a cover after seeing \(R\) and
then rerunning the rounding would be circular.  The simultaneous
capacity-cover assertion, and hence multiplier two, remains open.

**Proposition 6.5l (growing-degree squared-row preordering no-go).**  For
either the half-cut-neighbor slabs of Proposition 6.5d or the outgoing-half
slabs of Propositions 6.5c and 6.5h, write
\[
\begin{aligned}
 b_{xy}&={x\wedge y\over2},\\
 w^+_{xy}&=(\sqrt2-1)M+M-{|Q_A(x)+Q_A(y)|\over2},\\
 w^-_{xy}&=\sqrt2M-{|Q_A(x)-Q_A(y)|\over2},\\
 g_s(r)&=w_s^2-\langle b_s,r\rangle^2.
\end{aligned}                                              \tag{6.16bb}
\]
Choose either sign consistently.  Let \(D\ge1\) and set
\[
 \kappa={\pi^2\over4(\sqrt2-1)^2}=14.3810675004\ldots .
\]
If
\[
 (2D+1)^D\left[\left(1+{\kappa\over n-1}\right)^D-1\right]<1, \tag{6.16bc}
\]
then uniform Boolean expectation is nonnegative on the full degree-\(2D\)
preordering generated by the \(g_s\)'s and annihilates the Boolean ideal.
Consequently there is no identity
\[
 -1=\sum_J\sigma_J\prod_{j\in J}g_j
       +\sum_eh_e(r)(r_e^2-1),                              \tag{6.16bd}
\]
where \(J\) is square-free, each \(\sigma_J\) is SOS, and every displayed
term has raw polynomial degree at most \(2D\).  Here and only here,
"degree" means that raw degree before reduction modulo \(r_e^2=1\); this
proposition does not claim a lower bound for a differently encoded or
quotient-degree hierarchy.

Thus every fixed raw degree is blind for all sufficiently large \(n\).
More quantitatively, for each fixed \(\epsilon>0\), (6.16bc) holds, once
\(n\) is sufficiently large depending on \(\epsilon\), throughout
\[
 D\le(1-\epsilon){\log n\over\log\log n}.                  \tag{6.16be}
\]

*Proof.*  Put \(t=(\sqrt2-1)M\), \(m_s=\|b_s\|_2^2\).  Both widths in
(6.16bb) are at least \(t\), while Proposition 5.2 and the bivector support
formula give
\[
 t^2\ge{3-2\sqrt2\over\pi^2}n^2(n-1),
 \qquad m_s\le {n^2\over4},
 \qquad {m_s\over w_s^2}\le{\kappa\over n-1}=:\rho_n.     \tag{6.16bf}
\]
Let \(q\) be a multilinear polynomial of degree at most \(d\), put
\(Z_j=\langle b_j,r\rangle\), and let \(\mathbb E_0\) be uniform cube
expectation.  Hölder, followed by Bonami--Beckner for \(q\) and Khintchine
for the linear forms, gives
\[
 \mathbb E_0\!\left[q^2\prod_{j=1}^{\ell}Z_j^2\right]
 \le(2\ell+1)^{d+\ell}
       \left(\prod_{j=1}^{\ell}m_j\right)\mathbb E_0q^2.  \tag{6.16bg}
\]
No independence among the rows is used.  If \(|J|=k\) and
\(\deg q\le D-k\), expand \(\prod_{j\in J}g_j\), bound every nonconstant
term by (6.16bg), and use (6.16bf).  This yields
\[
\begin{aligned}
 \mathbb E_0\!\left[q^2\prod_{j\in J}g_j\right]
 &\ge\left(\prod_{j\in J}w_j^2\right)\mathbb E_0q^2\\
 &\quad\cdot\left\{1-(2D+1)^D[(1+\rho_n)^D-1]\right\}\ge0. \tag{6.16bh}
\end{aligned}
\]
After multilinearizing square summands without increasing degree, every
SOS localizer allowed in (6.16bd) is a sum of terms of this form.  Repeated
inequality factors, if permitted, have even powers absorbed into the square
multiplier and reduce to the same square-free case.  Uniform expectation is
also nonnegative on the free SOS term and vanishes on every multiple of
\(r_e^2-1\).  Applying it to (6.16bd) would give \(-1\ge0\), proving the
first assertion.

For fixed \(D\), the left side of (6.16bc) is \(O_D(1/n)\).  Under
(6.16be),
\(D\log(2D+1)\le(1-\epsilon/2)\log n\) for large \(n\), whereas
\((1+\kappa/(n-1))^D-1=(1+o(1))\kappa D/n\); their product tends to zero.
This proves the stated growing-degree range. \(\square\)

This is a proof-system obstruction, not an orientation construction.  It
strengthens Proposition 6.5j asymptotically, but it does not exclude a
different lifted encoding, raw degree growing faster than (6.16be), an
\(A\)-dependent nonuniform functional, or direct dependent rounding.

**Proposition 6.5m (conference commuting-mate parity obstruction).**  Let
\(A\) be symmetric and \(R\) skew-symmetric, both zero on the diagonal and
in \(\{\pm1\}\) off it.  If \(n\) is even, then for every \(i\)

\[
 (AR-RA)_{ii}
 =-2\sum_{j\ne i}A_{ij}R_{ij}\equiv2\pmod4.              \tag{6.16bi}
\]

Consequently

\[
 AR\ne RA,
 \qquad \|AR-RA\|_F^2\ge4n.                              \tag{6.16bj}
\]

In particular, no real symmetric conference signing commutes with any real
skew signing, whether or not the latter is itself conference.  Thus the
tempting exact specialization \(A^2=-R^2=(n-1)I\), \(AR=RA\), which would
make the block square in Proposition 6.5e scalar, is impossible at every
symmetric-conference order.

*Proof.*  Symmetry and skew-symmetry give the equality in (6.16bi).  Its sum
contains the odd number \(n-1\) of signs, proving the congruence and hence
(6.16bj).  There is also an independent parity proof.  Modulo two every
even-order skew signing reduces to \(J-I=J+I\), and

\[
 (J+I)^2=I\quad\hbox{over }\mathbb F_2,                   \tag{6.16bk}
\]

so \(R\) is invertible.  A symmetric conference matrix has order
\(n\equiv2\pmod4\), and its two eigenspaces have the odd dimension \(n/2\).
If \(A\) and \(R\) commuted, \(R\) would restrict to a real skew operator on
each of these odd-dimensional spaces and would therefore be singular, a
contradiction. \(\square\)

This is an exact-construction no-go, not a lower bound at the scale of the
open orientation problem.  The contribution \(4n=o(n^4)\) is compatible
with the approximate condition \(\Sigma(A,R)=o(n^4)\) in Proposition 6.5e.
In orthogonal-design terminology the desired commuting pair is called
*anti-amicable*; “amicable” symmetric/skew pairs anticommute instead.  The
scope and literature audit are recorded in
`evidence/NOTE_2026-09-02_CONFERENCE_COMMUTING_MATE_NO_GO.md`.

**Proposition 6.5n (optimal-scale coherent clique-flip counterfamily).**
For a complete signing \(A\), retain
\(L_{\rm cl}(A)=\max_{T\subseteq[n]}\Phi(A^{K_T})\) from Proposition
6.5c.  If a Hadamard matrix of order \(r\ge4\) exists and \(b\ge2\), put

\[
 n=4rb,\qquad N_\times=8r^2b(b-1).                       \tag{6.16bl}
\]

There is a complete signing \(A=A_{b,r}\) such that

\[
 \Phi(A)\le b(2r^2+2)
 +\sqrt{16r^2b(b-1)(4rb\log2+1)},
 \qquad
 L_{\rm cl}(A)\ge b(4r^2-2r).                            \tag{6.16bm}
\]

Consequently, along \(b=2^j\), \(r=256b\),
\(\Phi(A)=\Theta(n^{3/2})\), but

\[
 \liminf_{j\to\infty}
 {L_{\rm cl}(A)-\sqrt2\Phi(A)\over n^{3/2}}
 \ge {16(2-\sqrt2)\over4}-\sqrt{2\log2}
 =1.1657\ldots>0.                                        \tag{6.16bn}
\]

Thus the coherent estimate cannot follow from the bare scale hypothesis
\(\Phi(A)=O(n^{3/2})\), even with an \(o(n^{3/2})\) error.  This does not
disprove the still-open implication for global minimizers
\(\Phi(A)=m_n\).

*Proof.*  For an order-\(r\) Hadamard matrix \(C\), let
\(P=J_{2r}-I_{2r}\),
\(H=C\otimes\left(\begin{smallmatrix}1&-1\\-1&1\end{smallmatrix}\right)\),
and

\[
 G_r=\begin{pmatrix}P&H\\H^T&-P\end{pmatrix}.            \tag{6.16bo}
\]

Pair the coordinates in each \(2r\)-block.  If \(d,e\) are the within-pair
differences, \(a=\|d\|_2^2\), \(b=\|e\|_2^2\), and \(S,T\) are the two
coordinate sums, direct expansion gives
\[
 Q_{G_r}={S^2-T^2\over2}+4d^TCe,qquad
 |S|\le2(r-a),\quad |T|\le2(r-b),\quad
 |d^TCe|\le\sqrt{rab}.
\]
Consequently
\(Q_{G_r}\le2(r-a)^2+4r\sqrt a\le2r^2+2\); applying the same argument to
\(-Q_{G_r}\), with the blocks interchanged, gives the matching lower bound.
The equality states described by one active difference pair give
\(\Phi(G_r)=2r^2+2\).  If \(T_r\) is the second half of its coordinates,
then the all-one state gives
\(Q_{G_r^{K_{T_r}}}({\bf1})=4r^2-2r\); the rows and columns of the
order-two kernel in \(H\) sum to zero.

Take \(b\) disjoint copies of \(G_r\) and fill the \(N_\times\)
inter-block edges independently with signs, obtaining \(W\).  For a fixed
Boolean state, Hoeffding gives

\[
 \Pr\{|Q_W(x)|\ge t\}\le2e^{-t^2/(2N_\times)},\qquad
 t=\sqrt{2N_\times(n\log2+1)}.                            \tag{6.16bp}
\]

There are only \(2^{n-1}\) energies up to antipodes, so the union bound is
at most \(e^{-1}\).  Some deterministic filler therefore has
\(\Phi(W)<t\), and the triangle inequality proves the first bound in
(6.16bm).  Let \(T\) be the union of the \(T_r\)'s.  Replacing every
inter-block sign by its negative preserves the norm bound and reverses its
value at the single state \({\bf1}\) after the clique flip.  Choosing the
favorable sign leaves the internal contribution at least
\(b(4r^2-2r)\), proving the second bound.

For \(r=\kappa b\), division by
\(n^{3/2}=8\kappa^{3/2}b^3\) gives

\[
 \liminf {L_{\rm cl}(A)\over n^{3/2}}\ge{\sqrt\kappa\over2},
 \qquad
 \limsup {\Phi(A)\over n^{3/2}}
 \le{\sqrt\kappa\over4}+\sqrt{\log2}.                   \tag{6.16bq}
\]

Sylvester matrices supply all the required orders when \(\kappa=256\), and
(6.16bn) follows.  Proposition 5.2 supplies the matching positive-order
lower bound on \(\Phi(A)\), hence \(\Theta(n^{3/2})\). \(\square\)

This is an infinite complete-signing obstruction, not another finite-order
census.  It leaves coherent control logically available only through
genuine global minimality, or through a quantitatively near-minimal leading
constant; the noncoherent \(A\)-dependent target of Proposition 6.5c also
remains open.  See
`evidence/NOTE_2026-09-02_COHERENT_CLIQUE_OPTIMAL_SCALE_COUNTERFAMILY.md`.

**Proposition 6.6 (balanced Paley-skew shielding).**  Put
\[
 a={\sqrt2-1\over\pi},\qquad
 \rho=a^2,\qquad
 \tau={1-\sqrt{1-4\rho}\over2}.                   \tag{6.17}
\]
There are absolute constants \(c>0,N_0\) and a nonnegative function
\(\Omega\), with \(\Omega^*(N)=\sup_{u\ge N}\Omega(u)\), whose supremum
envelope satisfies
\[
 \sum_{j\ge0}\Omega^*(2^jN)\longrightarrow0       \tag{6.18}
\]
such that the following holds for every \(n\ge N_0\).  Let \(A\) be an
optimal order-\(n\) signing, let \(M=m_n\), and choose
\(z\in\{\pm1\}^n\) with \(|Q_A(z)|=M\).  There is a skew signing \(R\)
for which (6.13) holds for every pair \((x,y)\) satisfying at least one of
\[
\begin{aligned}
 \min(h_x,h_y)&\le\rho n,\\
 d_H(x,y)\bigl(n-d_H(x,y)\bigr)&\le\rho n^2,\\
 h_xh_y&\le {\rho\over4}n^2,\\
 2\sqrt2M-|Q_A(x)+Q_A(y)|&\ge n^{3/2},             \tag{6.19}
\end{aligned}
\]
where \(h_x=d_H(x,\{z,-z\})\) and similarly for \(h_y\).  Equivalently,
this explicit construction leaves only pairs satisfying all four strict
inequalities
\[
\begin{aligned}
 h_x,h_y&>\rho n,\\
 d_H(x,y)\bigl(n-d_H(x,y)\bigr)&>\rho n^2,\\
 h_xh_y&>{\rho\over4}n^2,\\
 |Q_A(x)+Q_A(y)|&>2\sqrt2M-n^{3/2}.                \tag{6.20}
\end{aligned}
\]
In particular, it shields linear neighborhoods of the chosen extremal
antipodal pair and of the diagonal and antidiagonal.  It does not prove the
complete multiplier-two ray.  If
\(M/n^{3/2}\le1/(2\sqrt2)-c_0\) for a fixed \(c_0>0\), the last threshold
in (6.20) is negative for large \(n\), so that condition is vacuous; in this
low-\(\alpha\) regime the theorem's gain is entirely geometric.

*Proof.*  The effective prime-number theorem in the fixed progression has
\(\pi(X;4,3)=\tfrac12\operatorname{Li}(X)
+O(X\exp(-c_1\sqrt{\log X}))\).  Applying it at the two endpoints of an
interval of length \(X\exp(-c_1\sqrt{\log X}/2)\) gives, for every large
\(n\), a prime
\[
 n\le q\le n+n\exp(-c\sqrt{\log n}),\qquad q\equiv3\pmod4.       \tag{6.21}
\]
Let \(T_q\) be the Paley tournament matrix
\((T_q)_{uv}=\chi(v-u)\).  Then
\[
 T_q\mathbf1=0,\qquad T_qT_q^T=qI-J,               \tag{6.22}
\]
so \(\|T_q\|_{\rm op}=\sqrt q\).  Compress to any \(n\)-set and call the
result \(T\); put \(k=q-n\).  The deleted columns in (6.22) give
\[
 \|T\mathbf1\|_1\le\sqrt{nqk}.                    \tag{6.23}
\]
Regard \(T\) as a tournament.  Repeatedly transfer one outdegree unit from
a maximum-degree vertex to a minimum-degree vertex, reversing their edge
when it has the needed orientation and otherwise reversing the two edges of
a directed two-path.  Use row-sum excess
\(\sum_i|(T\mathbf1)_i|\) for odd \(n\), and
\(\sum_i\max\{|(T\mathbf1)_i|-1,0\}\) for even \(n\).
Each transfer uses at most two edge reversals and lowers this potential by at
least two, so the total number of reversed edges is
at most \(\|T\mathbf1\|_1\).  Thus after at most
\(s\le\sqrt{nqk}\) reversals one obtains a tournament
matrix \(R_n=T+F\) with
\[
 \|R_n\mathbf1\|_\infty\le1,\qquad
 |u^TR_nv|\le n\sqrt q+4s                         \tag{6.24}
\]
for all Boolean \(u,v\).  With \(\delta=(q-n)/n\), the relative error in
(6.24) is bounded by
\[
 \varepsilon(n)=\sqrt{1+\delta}-1
       +4\sqrt{(1+\delta)\delta}.
\]
It is at most \(C\exp(-c'\sqrt{\log n})\).  Its supremum envelope has the
dyadic Dini property because
\(\sum_j\exp(-c'\sqrt{\log(2^jn)})\to0\).  The same is true after adding
any \(O(n^{-1/2})\) term.

Each edge reversal contributes at most \(4\) to \(a^TFb\) whenever
\(a,b\in[-1,1]^n\).  This observation also applies below to indicator
vectors, so the displayed \(s\)-errors do not require Boolean inputs.

Conjugate by the maximizer, \(R=\operatorname{diag}(z)R_n
\operatorname{diag}(z)\).  Since \(\|Rz\|_\infty\le1\), every pair incident
to \(\{z,-z\}\) has skew cost at most \(n\).  More generally, if a Boolean
vector differs from an anchor in \(h\) places, then, writing
\(R_n=T+F\),
\[
 |u^TR_nv|\le n+2\sqrt{qhn}+8s.                    \tag{6.25}
\]
For \(h\le\rho n\), the leading term on the right is
\(2a\sqrt{1+\delta}\,n^{3/2}\).  On the other hand Proposition 5.2 gives,
uniformly in \(x,y\),
\[
 2\sqrt2M-|Q_A(x)+Q_A(y)|
 \ge2(\sqrt2-1)M
 \ge2a\sqrt{1-1/n}\,n^{3/2}.                      \tag{6.26}
\]
The difference between (6.25) and (6.26) is
\(O(n^{-1/2}+\sqrt\delta)n^{3/2}\), with the required Dini tail.  This
proves the first line of (6.19), in either coordinate by skew symmetry.

If \(r=d_H(u,v)\), skewness cancels both within-part blocks and gives the
exact cut identity
\[
 u^TTv=-2v_S^TT_{S,S^c}v_{S^c},\qquad S=\{i:u_i\ne v_i\}.
\]
Consequently
\[
 |u^TR_nv|\le2\sqrt{qr(n-r)}+4s.                  \tag{6.27}
\]
Equations (6.17), (6.26), and (6.27) prove the second line of (6.19); it is
equivalent to
\(r\le\tau n\) or \(r\ge(1-\tau)n\).

Finally choose independent global signs so that
\(u=\mathbf1-2\mathbf1_S\),
\(v=\mathbf1-2\mathbf1_U\), with \(|S|=h_x\), \(|U|=h_y\).
Expanding around \(\mathbf1^TR_n\mathbf1=0\) and using the balanced row
sums gives
\[
 |u^TR_nv|\le2(h_x+h_y)+4\sqrt{qh_xh_y}+16s.       \tag{6.28}
\]
Since each \(h\le n/2\), the third line of (6.19) again puts (6.28) below
(6.26) up to \(O(n^{-1/2}+\sqrt\delta)n^{3/2}\).  The global estimate
(6.24) proves the fourth line.  Taking \(\Omega\) to dominate these four
explicit errors proves (6.18)--(6.19); for suitable absolute constants,
\[
 \Omega(n)=C_0\left(n^{-1/2}+
             \exp(-c_2\sqrt{\log n})\right)
\]
works.  Complementing the alternatives gives (6.20). \(\square\)

Proposition 6.6 is a genuine infinite-family reduction, but its last set is
not known to be empty.  The next multiplier-two target is therefore exact:
choose the Paley principal embedding and the degree-balancing reversals so
that (6.13) also holds on (6.20). An independent skew-norm budget cannot do
this by (6.15). The asymptotic disk (6.16) is not disproved, but it is a
strictly stronger lower-bound problem rather than a neutral reformulation;
its zero-error form is false.

**Proposition 6.7 (tetrahedral tripling frame and exact diamond).**  Let
\(A\) be an order-\(n\) signing, put \(M=\Phi(A)\), and let \(P,Q,T\) be
skew signings of order \(n\).  For arbitrary diagonal sign matrices
\(D_{12},D_{13},D_{23}\), the block matrix
\[
 {\cal S}=\begin{pmatrix}
 A&P+D_{12}&Q+D_{13}\\
 -P+D_{12}&A&T+D_{23}\\
 -Q+D_{13}&-T+D_{23}&-A
 \end{pmatrix}                                                   \tag{6.29}
\]
is an order-\(3n\) signing.  For Boolean \(x,y,z\), write
\[
 I=Q_A(x)+Q_A(y)-Q_A(z),\qquad
 b=x^TPy,\quad c=x^TQz,\quad d=y^TTz
\]
and define
\[
 K_3(A;P,Q,T)=\max_{x,y,z}\max\left\{
 |I+d|+|b+c|,\ |I-d|+|b-c|\right\}.              \tag{6.30}
\]
Then
\[
 \bigl|\Phi({\cal S})-K_3(A;P,Q,T)\bigr|\le3n.   \tag{6.31}
\]
Moreover, after grouping the three copies of coordinate \(i\) into one
cloud, every one of the four projective constant cloud states induces
exactly the signing \(A\).  Explicitly, for \(i<j\) the inter-cloud block is
\[
 B_{ij}=A_{ij}\operatorname{diag}(1,1,-1)+
 \begin{pmatrix}
 0&P_{ij}&Q_{ij}\\
 -P_{ij}&0&T_{ij}\\
 -Q_{ij}&-T_{ij}&0
 \end{pmatrix},
 \qquad v^TB_{ij}v=A_{ij}                       \tag{6.32}
\]
for \(v\in\{(1,1,1),(1,-1,-1),(-1,1,-1),
(-1,-1,1)\}\).  Hence endpoint selection is not an obstruction for this
tripling frame.  The multiplier-three ray would follow if, for an optimal
\(A\), one could choose \(P,Q,T\) so that
\[
 K_3(A;P,Q,T)\le3^{3/2}M+n^{3/2}\Omega(n),        \tag{6.33}
\]
where \(\Omega\ge0\) and its supremum envelope has a vanishing dyadic Dini
tail.

There is a useful one-signing specialization.  Taking
\(P=R,\ Q=-R,\ T=R\) and putting
\[
 C_R(x,y,z)=x^TRy+y^TRz+z^TRx=(x-y)^TR(y-z),
\]
one has the exact diamond
\[
 K_3(A;R,-R,R)=\max_{x,y,z}
 \left(|Q_A(x)+Q_A(y)-Q_A(z)|+|C_R(x,y,z)|\right). \tag{6.34}
\]
To expose its endpoint geometry, put \(s=x\circ y\circ z\) and
\[
\begin{array}{ll}
 u_0=(s+x+y+z)/4,&u_1=(s+x-y-z)/4,\\
 u_2=(s-x+y-z)/4,&u_3=(s-x-y+z)/4.
\end{array}
\]
The \(u_j\) have disjoint supports, take values in \(\{0,\pm1\}\), and
their supports partition \([n]\).  Exactly,
\[
\begin{aligned}
 I&=Q_A(s)-4\bigl(u_0^TAu_3+u_1^TAu_2\bigr),\\
 C_R&=4\bigl(u_1^TRu_2+u_2^TRu_3+u_3^TRu_1\bigr).              \tag{6.35}
\end{aligned}
\]
If \(\Lambda=\|R\|_{\rm op}\) and
\(h=|\operatorname{supp}u_1|+|\operatorname{supp}u_2|
 +|\operatorname{supp}u_3|\), then the following two shields are automatic:
\[
\begin{aligned}
 |C_R|&\le4\Lambda\min_{\rm cyc}
 \sqrt{d_H(x,y)d_H(y,z)},\\
 |C_R|&\le {4\over\sqrt3}\Lambda h.              \tag{6.36}
\end{aligned}
\]
Thus the integrand in (6.34) is at most \(3\sqrt3M\) for a triple whenever
either right-hand side in (6.36) is at most \(3\sqrt3M-|I|\).  The complement
of these two shields---the **tetrahedral diamond**---is open; Proposition 6.7
does not prove the multiplier-three ray.

*Proof.*  Every off-diagonal entry of (6.29) is a sign, and the transpose
blocks match, so \({\cal S}\) is a signing.  Its energy at \((x,y,z)\) is
\[
 Q_{\cal S}(x,y,z)=I+b+c+d+\Delta(x,y,z),
\]
where the exact internal-cloud contribution is
\[
 \Delta=x^TD_{12}y+x^TD_{13}z+y^TD_{23}z,
 \qquad |\Delta|\le3n.                            \tag{6.37}
\]
The four global layer-sign choices leave \(I\) fixed and give the four
cross-term sign patterns whose product is \(+1\).  Their principal parts are
\[
 I+b+c+d,\quad I-b-c+d,\quad I-b+c-d,\quad I+b-c-d.
\]
Pairing the first two and the last two and using
\(\max(|a+r|,|a-r|)=|a|+|r|\) gives (6.30).  Adding (6.37), then maximizing,
proves both sides of (6.31).

Equation (6.32) follows because its second summand is skew and
\(1^2+1^2-1^2=1\).  This endpoint form is also forced.  Indeed, if a
\(3\times3\) sign block \(B\) has \(v^TBv=a\) at all four displayed states,
expansion in the three tetrahedral characters \(v_1v_2,v_1v_3,v_2v_3\)
gives \(B_{ij}+B_{ji}=0\) for \(i\ne j\) and
\(\operatorname{tr}B=a\).  For \(a\in\{\pm1\}\), its diagonal therefore
has two copies of \(a\) and one copy of \(-a\), up to permutation.  Thus
(6.29) is the general coherent equal-endpoint frame with a fixed exceptional
layer.

If (6.33) holds, (6.31) gives
\[
 m_{3n}\le3^{3/2}m_n+n^{3/2}\Omega(n)+3n.
\]
Since \(m_n=\Theta(n^{3/2})\), the mean-value theorem gives
\(H(3n)\le3H(n)+O(n\Omega(n)+\sqrt n)\), with a Dini-summable normalized
error.  This proves the asserted implication and also shows why the leading
constant \(3^{3/2}\) in (6.33) is sharp for this route.

For the specialization, the domain of the maximum is invariant under all
global layer flips, so (6.30) is equivalently the maximum of
\(|I+b+c+d|\).  Here \(b+c+d=C_R\).  Swapping \(x,y\) fixes \(I\) and
negates \(C_R\), which proves (6.34).  The four-partition identities (6.35)
follow by substituting
\[
 s=u_0+u_1+u_2+u_3,\quad x=u_0+u_1-u_2-u_3,\quad
 y=u_0-u_1+u_2-u_3,\quad z=u_0-u_1-u_2+u_3.
\]
The first bound in (6.36) follows from the displayed difference identity
for \(C_R\) and its two cyclic versions.  For the second, let
\(\zeta=e^{2\pi i/3}\) and \(w=u_1+\zeta u_2+\zeta^2u_3\).  Disjointness
gives \(\|w\|_2^2=h\), while skewness gives
\[
 w^*Rw=i\sqrt3\bigl(u_1^TRu_2+u_2^TRu_3+u_3^TRu_1\bigr).
\]
Thus \(|C_R|\le4\Lambda h/\sqrt3\), completing the proof. \(\square\)

**Proposition 6.8 (bi-balanced Hadamard shield for the \(1:2\)
composition).**  Let \(A\) and \(B\) be independently optimal signings of
orders \(n\) and \(2n\), and put
\[
 M=m_n,\qquad N=m_{2n},\qquad
 T=(M^{2/3}+N^{2/3})^{3/2}.
\]
Choose positive and negative extremizers \(z_+,z_-\) of \(A\), and
\(w_+,w_-\) of \(B\).  For every sufficiently large \(n\) there is an
\(n\)-by-\(2n\) sign matrix \(C\) such that the order-\(3n\) signing
\[
 {\cal J}=\begin{pmatrix}A&C\\ C^T&B\end{pmatrix}               \tag{6.38}
\]
has the exact two-state form
\[
 \Phi({\cal J})=
 \max_{x,y}\bigl(|Q_A(x)+Q_B(y)|+|x^TCy|\bigr).                 \tag{6.39}
\]

More precisely, pair the row coordinates within the two relative-sign
classes of \(z_+\circ z_-\), and pair the column coordinates within the two
classes of \(w_+\circ w_-\).  For a state \(x\), let \(k_A(x)\) be the
number of row pairs on which \(z_+\circ x\) is nonconstant, and define
\(k_B(y)\) analogously.  Then all four Boolean states
\(\mathord\pm z_+,\mathord\pm z_-\) have \(k_A=0\), and likewise on the
\(B\)-side.  The matrix \(C\) can be chosen so that
\[
 |x^TCy|\le4\sqrt{q\,k_A(x)k_B(y)}+6n,                          \tag{6.40}
\]
where \(q/n=1+O(\exp(-c\sqrt{\log n})+n^{-1})\).  Consequently
\[
 k_A(x)k_B(y)\le {n^2\over100}
 \quad\Longrightarrow\quad
 |Q_A(x)+Q_B(y)|+|x^TCy|\le T.                                 \tag{6.41}
\]
The exact residual for this construction consists only of pairs satisfying
both
\[
 k_A(x)k_B(y)>{n^2\over100}                                    \tag{6.42}
\]
and
\[
 4\sqrt{q\,k_A(x)k_B(y)}+6n
   >T-|Q_A(x)+Q_B(y)|.                                         \tag{6.43}
\]
Thus the construction shields fixed linear Hamming strips around every
chosen positive and negative extremizer on either side.  It is a proved
infinite-family reduction, not a proof of the multiplier-three ray.

*Proof.*  Flipping all signs of \(x\) fixes the two quadratic terms in
\(Q_{\cal J}\) and reverses \(x^TCy\).  The identity
\(\max(|a+b|,|a-b|)=|a|+|b|\) proves (6.39).

For positive \(a,b\), the function
\[
 f(a,b)=(a^{2/3}+b^{2/3})^{3/2}-a-b
\]
is increasing in each variable.  Proposition 5.2 therefore gives the
uniform headroom
\[
 {T-M-N\over n^{3/2}}
 \ge d_0-o(1),\qquad
 d_0={3\sqrt3-1-2\sqrt2\over\pi}
     =0.4353604839\ldots .                                     \tag{6.44}
\]
No doubling estimate is used here.

Gauge the rows by \(z_+\) and the columns by \(w_+\).  Pair coordinates
inside the two level sets of each relative-sign vector.  There are
\(r\le n/2\) row pairs and \(s\le n\) column pairs, with at most two
unpaired coordinates on either side.  On a row pair \(\{i,i'\}\) put
\(u_a=(x_i-x_{i'})/2\), after gauging, and define \(v_b\) analogously on
column pairs.  Then \(k_A=|\operatorname{supp}u|\) and
\(k_B=|\operatorname{supp}v|\).

The effective prime-number theorem in the progression \(3\bmod4\), as in
Proposition 6.6, supplies a prime \(p\ge n\) with
\(p-n\le n\exp(-c\sqrt{\log n})\).  A Paley Hadamard matrix of order
\(q=p+1\) exists.  Take any \(r\)-by-\(s\) submatrix \(E\), so
\[
 \|E\|_{\rm op}\le\sqrt q,qquad
 {q\over n}=1+O(\exp(-c\sqrt{\log n})+n^{-1}).                 \tag{6.45}
\]
On row pair \(a\) and column pair \(b\), put the sign tile
\[
 E_{ab}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]
Fill entries incident to unpaired coordinates arbitrarily and undo the
gauges.  The paired core is exactly \(4u^TEv\).  At most \(6n\) entries
lie on the border, hence
\[
 x^TCy=4u^TEv+\operatorname{border},\qquad
 |\operatorname{border}|\le6n,
\]
which proves (6.40).  If \(k_Ak_B\le n^2/100\), its leading term is at
most \((0.4+o(1))n^{3/2}\).  Since \(0.4<d_0\), (6.44) proves (6.41) for
all sufficiently large \(n\).  Complementing this conclusion and retaining
the sharper statewise headroom gives (6.42)--(6.43). \(\square\)

For any fixed finite lists of row and column anchors, one may refine the
coordinate classes by their complete relative-sign signatures and pair
inside each class.  Every listed anchor then has zero paired difference.
There are at most \(2^{R-1}\) and \(2^{S-1}\) leftovers for \(R\) row and
\(S\) column anchors, so fixed \(R,S\) change only the \(O(n)\) border in
(6.40).  This is the precise finite-anchor halving extension of Proposition
6.8.  If the small block in (6.38) is instead \(-A_0\), while \(B\) is the
equal-endpoint doubling frame whose two diagonal blocks are \(A_0\), then
regrouping recovers Proposition 6.7's \(K_3\), up to its three diagonal
matching terms of total size at most \(3n\).  Thus nested doubling gives no
new tripling control; the new two-state gate depends essentially on taking
\(B\) independently optimal.

**Proposition 6.9 (conference obstruction to every fixed-temperature
free-energy target).**  Let \(P\) be an orthogonal projection on
\(\mathbb R^n\) with rank \(n/2\) and diagonal entries \(P_{ii}=1/2\).
For every \(t>0\),
\[
 \mathbb E_x e^{-t x^TPx}
 \le\left({1+e^{-2t}\over2}\right)^{n/2}.                      \tag{6.46}
\]
Consequently, if \(C\) is a symmetric conference signing,
\(C^2=(n-1)I\), then for every \(\beta\ge0\),
\[
 \mathbb E_x e^{\mathord\pm\beta Q_C(x)}
 \le\cosh\!\left(\beta\sqrt{n-1}\right)^{n/2},               \tag{6.47}
\]
and hence, along any sequence of symmetric conference signings \(C_n\) whose
orders tend to infinity, for every fixed \(c>0\),
\[
 \limsup_{n\to\infty}{1\over n}
 \log\mathbb E_x\cosh\!\left({c\over\sqrt n}Q_{C_n}(x)\right)
 \le {1\over2}\log\cosh c < {c\over2}.                      \tag{6.48}
\]
There are infinitely many such Paley conference matrices.  Therefore the
uniform fixed-\(c\) sufficient target
\[
 \inf_A\log\mathbb E_x\cosh\!\left({c\over\sqrt n}Q_A(x)\right)
 \ge {c\over2}n-o(n)                                          \tag{6.49}
\]
is false for every fixed \(c>0\).  At \(c=3\) it fails by at least
\[
 {3-\log\cosh3\over2}\,n-o(n)
 =0.3453357477\ldots n-o(n).                                  \tag{6.50}
\]

*Proof.*  Write \(P=U^TU\), where \(UU^T=I_{n/2}\), and let \(u_i\) be
the columns of \(U\).  The Gaussian Fourier identity and averaging the
Rademacher signs give
\[
 \mathbb E_xe^{-t x^TPx}
 =(4\pi t)^{-n/4}\int e^{-\|\xi\|^2/(4t)}
       \prod_i\cos\langle u_i,\xi\rangle\,d\xi.
\]
Take absolute values, put \(v_i=\sqrt2u_i\), and note that
\(\|v_i\|=1\) and \(\sum_i\tfrac12v_iv_i^T=I\).  Apply the geometric
rank-one Brascamp--Lieb inequality with weights \(1/2\) to
\[
 f(s)=e^{-s^2/(4t)}\cos^2(s/\sqrt2).
\]
The integral is at most \((\int_{\mathbb R}f)^{n/2}\), while
\[
 \int_{\mathbb R}f(s)\,ds
 =\sqrt{4\pi t}\,{1+e^{-2t}\over2}.
\]
The Gaussian normalization cancels, proving (6.46).

Put \(\lambda=\sqrt{n-1}\) and
\(P_\pm=(I\pm C/\lambda)/2\).  These are rank-\(n/2\) projections with
diagonal \(1/2\), and
\[
 Q_C(x)=\lambda\bigl(n/2-x^TP_-x\bigr),\qquad
 -Q_C(x)=\lambda\bigl(n/2-x^TP_+x\bigr).
\]
Applying (6.46) to the two signs gives (6.47); taking
\(\beta=c/\sqrt n\) gives (6.48).  Paley matrices of order \(q+1\) for
primes \(q\equiv1\pmod4\) supply the infinite counterfamily. \(\square\)

Equivalently, for the signed even-Eulerian polynomial
\[
 \mathbb E_x\cosh(\beta Q_C(x))
 = (\cosh\beta)^{\binom n2}P_C(\tanh\beta),
\]
Proposition 6.9 gives
\[
 \limsup {1\over n}\log P_{C_n}(\tanh(c/\sqrt n))
 \le {1\over2}\log\cosh c-{c^2\over4},                         \tag{6.51}
\]
strictly below the formerly proposed lower bound
\((c/2-c^2/4)n-o(n)\) for every fixed \(c>0\).  This does not rule out a
growing temperature \(c=c_n\to\infty\). With
\(t_n=c_n\sqrt{1-1/n}\), the conference loss relative to \(c_n n/2\) is
\[
 {n\over2}\bigl(c_n-\log\cosh t_n\bigr)
 ={\log2\over2}n+{c_n\over2(1+\sqrt{1-1/n})}+o(n)
 =O(n+c_n)=o(c_n n).
\]
Thus a growing-temperature criterion with a correspondingly uniform error
remains logically possible.

**Proposition 6.10 (critical-pressure gate and graphon no-go).**  For
\(c>0\), define the optimized critical symmetric pressure

\[
 s_n(c)={1\over n}\min_A\log\mathbb E_x
 \cosh\!\left({c\over\sqrt n}Q_A(x)\right).                 \tag{6.52}
\]

Then, exactly for every \(n\) and \(c>0\),

\[
 {s_n(c)\over c}\le\alpha_n
 \le {s_n(c)\over c}+{\log2\over c}.                        \tag{6.53}
\]

In particular,

\[
 \limsup_n\alpha_n-\liminf_n\alpha_n
 \le {\log2\over c}
 +{\limsup_n s_n(c)-\liminf_n s_n(c)\over c}.               \tag{6.54}
\]

Thus convergence of \(s_n(c)\) for every \(c\) in any unbounded set would
prove convergence of \(\alpha_n\); more generally, pressure oscillation
\(o(c)\) along some \(c\to\infty\) is sufficient.

The elementary block interpolation for this pressure has the wrong
critical scaling.  At a common raw inverse temperature put

\[
 p_n(\beta)=\min_A\log\mathbb E_xe^{\beta Q_A(x)}.
\]

For all \(n,m\ge1\) and \(\beta\ge0\),

\[
 p_n(\beta)+p_m(\beta)\le p_{n+m}(\beta)
 \le p_n(\beta)+p_m(\beta)+nm\log\cosh\beta.                \tag{6.55}
\]

The symmetric version has the corresponding upper inequality.  Namely, if
\(q_n^{\rm sym}(\beta)=\min_A\log\mathbb E\cosh(\beta Q_A)\),
then
\[
 q_{n+m}^{\rm sym}(\beta)\le q_n^{\rm sym}(\beta)
 +q_m^{\rm sym}(\beta)+nm\log\cosh\beta.                   \tag{6.56}
\]
Consequently, with
\(f_n^+(c)=n^{-1}p_n(c/\sqrt n)\), equal splitting gives

\[
 f_n^+(c/\sqrt2)\le f_{2n}^+(c)
 \le f_n^+(c/\sqrt2)
 +{n\over2}\log\cosh\!\left({c\over\sqrt{2n}}\right),       \tag{6.57}
\]

whose last term tends to \(c^2/8\).  The child parameter is also
\(c/\sqrt2\), not \(c\).  Thus the common-\(\beta\) Fekete argument leaves a
nonvanishing per-spin defect on the critical diagonal.

Nor is the critical pressure continuous in the ordinary signed cut metric.
For every fixed \(0<c<1\), there are deterministic signings \(C_n,J_n\), on
an infinite common sequence of orders, such that

\[
 \|C_n\|_\square,\ \|J_n\|_\square\longrightarrow0,
 \qquad
 \|A\|_\square={1\over n^2}\max_{S,T}
 \left|\sum_{i\in S,j\in T}A_{ij}\right|,                  \tag{6.58}
\]

but

\[
 \limsup_n{1\over n}\log\mathbb E_x
 \cosh\!\left({cQ_{C_n}(x)\over\sqrt n}\right)
 \le {1\over2}\log\cosh c
 <{c^2\over4}
 =\lim_n{1\over n}\log\mathbb E_x
 \cosh\!\left({cQ_{J_n}(x)\over\sqrt n}\right).             \tag{6.59}
\]

Hence ordinary dense graphon, left-, or right-convergence cannot by itself
prove convergence of (6.52).

*Proof.*  For each \(A\),
\(\mathbb E\cosh(cQ_A/\sqrt n)\le e^{c\Phi(A)/\sqrt n}\).
A maximizing state and its antipode contribute at least
\(e^{c\Phi(A)/\sqrt n}/2^n\) to the normalized expectation.  Taking logs,
then minimizing over \(A\), proves (6.53), and (6.54) follows.

For (6.55), split an arbitrary signing into diagonal blocks \(B,D\) and a
cross block \(E\).  Pairing \(y\) with \(-y\) gives

\[
 \mathbb E_{x,y}e^{\beta(Q_B(x)+Q_D(y)+x^TEy)}
 =\mathbb E_{x,y}e^{\beta(Q_B(x)+Q_D(y))}
   \cosh(\beta x^TEy)
 \ge Z_B(\beta)Z_D(\beta).
\]

Conversely, begin with minimizing \(B,D\) and choose the \(nm\) entries of
\(E\) independently.  The completion average contributes exactly
\(\cosh(\beta)^{nm}\), so one integral completion attains the upper bound.
For (6.56), write
\(a_A=\mathbb E\cosh(\beta Q_A)\) and
\(b_A=\mathbb E\sinh(\beta Q_A)\).  Replace one diagonal block by its
negative so that \(b_Bb_D\le0\); averaging the random cross completion gives
at most \(a_Ba_D\cosh(\beta)^{nm}\).  Equation (6.57) is (6.55) divided by
\(2n\) at \(\beta=c/\sqrt{2n}\).

Take \(C_n\) to be Paley symmetric conference signings.  Their operator norm
is \(\sqrt{n-1}\), so (6.58) holds, and Proposition 6.9 gives the first rate
in (6.59).  For the other sequence, let \(J\) have independent signs above
the diagonal and set

\[
 X_J=\mathbb E_xe^{cQ_J(x)/\sqrt n},\qquad
 M_n=\mathbb E_JX_J
 =\cosh(c/\sqrt n)^{\binom n2}.                            \tag{6.60}
\]

For two spin states put \(z_i=x_iy_i\), \(R=\sum_i z_i\), and
\(u=\tanh^2(c/\sqrt n)\).  Direct edge averaging gives

\[
 {\mathbb E_JX_J^2\over M_n^2}
 =\mathbb E_z(1+u)^{(n^2+R^2-2n)/4}
              (1-u)^{(n^2-R^2)/4}
 \le\mathbb E_ze^{c^2R^2/(2n)}
 \le{1\over\sqrt{1-c^2}}.                                 \tag{6.61}
\]

The last inequality follows by Gaussian linearization and
\(\cosh t\le e^{t^2/2}\).  Paley--Zygmund therefore gives a fixed positive
probability that \(X_J\ge M_n/2\).  Markov's inequality gives, for a fixed
large \(K\),
\(\mathbb E_x\cosh(cQ_J(x)/\sqrt n)\le KM_n\) with probability sufficiently
close to one.  A Hoeffding union bound over the at most \(4^n\) cut pairs
also gives \(\|J\|_\square=O(n^{-1/2})\) with probability tending to one.
These events intersect, and choosing one such \(J_n\) yields a symmetric
partition function between \(M_n/4\) and \(KM_n\).  Since
\(n^{-1}\log M_n\to c^2/4\), this proves the second rate in (6.59).
Finally \(c^2/4-\tfrac12\log\cosh c>0\), because its derivative is
\((c-\tanh c)/2>0\). \(\square\)

This proposition is a proved sufficient reformulation and two method
barriers, not convergence of \(s_n(c)\) or \(\alpha_n\).  The optimized
Bernoulli lower tail at disorder-space speed \(n^2\), or a stronger
second-order limit object distinguishing conference from independent
quasirandom signings, remains open.  See
`evidence/NOTE_2026-09-02_THERMODYNAMIC_INTERPOLATION_GATE.md`.

**Proposition 6.10a (product-pressure lower-curve no-go and entropy
fallback).**  The proposed universal lower curve

\[
 \mathbb E_x\cosh\!\left({cQ_A(x)\over\sqrt n}\right)
 \ge
 \cosh\!\left(c\sqrt{1-{1\over n}}\right)^{n/2}           \tag{6.62}
\]

is false for every \(c>0\), already for one optimal signing of order five.
Specifically, let \(A_5\) be negative on a five-cycle and positive on its
five diagonals.  Its energy histogram on the sixteen antipodal cuts is

\[
 \#\{Q_{A_5}=-4,0,4\}=\{5,6,5\},\qquad \Phi(A_5)=m_5=4.  \tag{6.63}
\]

With \(t=2c/\sqrt5\), \(u=\cosh t\), and \(v=\sqrt u>1\), the two sides
of (6.62) are respectively

\[
 {5u^2-1\over4}
 \quad\hbox{and}\quad u^{5/2},
 \qquad
 u^{5/2}-{5u^2-1\over4}
 ={(v-1)^2(4v^3+3v^2+2v+1)\over4}>0.                    \tag{6.64}
\]

There is nevertheless a universal all-temperature entropy curve.  Put

\[
 I(r)={1+r\over2}\log(1+r)+{1-r\over2}\log(1-r),
 \qquad 0\le r\le1.                                      \tag{6.65}
\]

Then every complete signing satisfies

\[
 {1\over n}\log\mathbb E_x\cosh\!\left({cQ_A(x)\over\sqrt n}\right)
 \ge \max\!\left(0,\sup_{0\le r\le1}
 \left\{cr^2{\Phi(A)\over n^{3/2}}-I(r)-{\log2\over n}\right\}\right)
                                                               \tag{6.66}
\]

and hence, by Proposition 5.2,

\[
 s_n(c)\ge \max\!\left(0,\sup_{0\le r\le1}
 \left\{{c\sqrt{1-1/n}\over\pi}r^2-I(r)-{\log2\over n}\right\}\right).
                                                               \tag{6.67}
\]

The large-\(c\) slope of (6.67) is only
\(\sqrt{1-1/n}/\pi\).  Thus this unconditional fallback repackages the
known \(1/\pi\) floor; it does not approach the construction scale
\(1/2\) or prove convergence.

*Proof.*  Equation (6.63) follows by representing a cut by a subset of size
at most two: the empty cut and five singletons have energy zero, the five
cycle-edge pairs have energy \(-4\), and the five diagonal pairs have energy
\(+4\).  Averaging their hyperbolic cosines gives the first expression in
(6.64); the factorization proves the strict reversal for all \(c>0\).

The first local failure is also exact.  If
\(C_4(A)=\sum_\gamma\prod_{e\in\gamma}a_e\), summed over the three
unoriented four-cycles on every four-set, then

\[
 \mathbb E Q_A^4=3N^2-2N+24C_4(A),\qquad N={n\choose2}.    \tag{6.68}
\]

For \(A_5\), \(C_4(A_5)=-5\), so the fourth-derivative difference between
the two sides of (6.62) is \(-12/5\), equivalently their difference is
\(-c^4/10+O(c^6)\).

For (6.66), choose \(x_*\) and \(\sigma_*\) with
\(\sigma_*Q_A(x_*)=\Phi(A)\).  In the augmented representation
\(Z_A(c)=\mathbb E_{\sigma,x}e^{\sigma cQ_A(x)/\sqrt n}\), fix
\(\sigma=\sigma_*\) and independently bias every coordinate toward
\(x_*\) with mean \(r\).  Its mean signed energy is \(r^2\Phi(A)\), while
its relative entropy from the uniform augmented measure is
\(nI(r)+\log2\).  The Gibbs variational principle gives (6.66), and
\(Z_A(c)\ge1\) supplies the outer maximum with zero.  Proposition 5.2 and
minimization over \(A\) give (6.67).  Finally the stationary equation is
\(\operatorname {atanh}r=2c\sqrt{1-1/n}\,r/\pi\), and direct division by
\(c\) proves the asserted limiting slope. \(\square\)

This retires only the literal product lower curve (6.62), including methods
whose conclusion is that same finite-order inequality.  It does not
invalidate Proposition 6.10's optimized pressure gate or exclude a lower
bound with an order-dependent exponential loss.  Any asymptotic universal
curve with large-\(c\) slope \(1/2\) would, through (6.53), already prove
the missing \(1/2\) lower bound.  See
`evidence/NOTE_2026-09-02_PRESSURE_LOWER_CURVE_NO_GO.md`.

---

## §7. Multipartite bounds

**Proposition 7.1 (random multipartite).** For \(k,n\ge2\),
\begin{equation}
\label{eq:rand}
m_{kn}\le k\,m_n+\sqrt{\log2}\,(kn)^{3/2},
\qquad
\alpha_{kn}\le\frac{\alpha_n}{\sqrt k}+\sqrt{\log2}.
\end{equation}

*Proof.* Optimal order-\(n\) matrices on \(k\) diagonal blocks; i.i.d. \(\pm1\) on cross blocks. For fixed \(x\), the cross contribution is a Rademacher sum of at most \(\binom{k}{2}n^2\) terms. Hoeffding and a union bound over \(2^{kn}\) vertices give the stated cross bound with positive probability. \(\square\)

**Proposition 7.2 (Hadamard doubling).** If a Hadamard matrix of order \(n\) exists, then
\begin{equation}
\label{eq:had}
m_{2n}\le 2m_n+n^{3/2},
\qquad
\alpha_{2n}\le\frac{\alpha_n}{\sqrt2}+\frac1{2\sqrt2}.
\end{equation}
*Proof.* Optimal diagonal blocks; Hadamard cross \(H\); \(\|Hu\|_1\le\sqrt n\,\|Hu\|_2=n^{3/2}\) for \(u\in\{\pm1\}^n\) is not always true, but Cauchy–Schwarz gives \(|x_1^\top H x_2|\le\|x_1\|_2\|Hx_2\|_2=n^{3/2}\). \(\square\)

Sylvester orders alone do **not** have ratio-one gaps and therefore do not
justify an \(o(1)\) extension to every \(n\).  Instead choose a prime
\(q\equiv3\pmod4\) with \(n\le q+1=n+o(n)\), use the Paley Hadamard matrix
of order \(q+1\), and compress it to any \(n\)-by-\(n\) submatrix \(C\).
Then \(\|C\|_{\rm op}\le\sqrt{q+1}=(1+o(1))\sqrt n\); using \(C\) between
two optimal order-\(n\) diagonal blocks proves the asymptotic version of
\eqref{eq:had}.  The prime-number theorem in this fixed progression
(quantified in (6.21)) supplies the ratio-one Hadamard orders.

**Remark.** Filling cross blocks with the constant signing \(+1\) yields error \(\binom{k}{2}n^2\) in \(m_{kn}\), hence error \(\Theta(\sqrt{kn})\) in \(\alpha_{kn}\), which tends to infinity as \(n\to\infty\) and cannot be used for limit arguments at fixed \(k\).

**Proposition 7.3 (reverse multipartite, new).** For all \(k,n\ge2\),
\begin{equation}
\label{eq:rev}
m_{kn}\ge\frac k2\,m_n,
\qquad
\alpha_{kn}\ge\frac{\alpha_n}{2\sqrt k}.
\end{equation}

*Proof.* Let \(B\) achieve \(m_{kn}\), and fix any partition of \([kn]\) into \(k\) blocks of size \(n\). Write \(Q_B(y)=\sum_b Q_b(y_b)+\sum_{b<b'}y_b^\top B^{bb'}y_{b'}\). For fixed \(x_1,\dots,x_k\in\{\pm1\}^n\) and \(\varepsilon\in\{\pm1\}^k\), set \(y_b=\varepsilon_b x_b\). Then
\[
Q_B(\varepsilon\cdot x)=\sum_b Q_b(x_b)+\sum_{b<b'}\varepsilon_b\varepsilon_{b'}\,x_b^\top B^{bb'}x_{b'},
\]
so \(\mathbb E_\varepsilon Q_B=\sum_b Q_b(x_b)\). Hence \(\bigl|\sum_b Q_b(x_b)\bigr|\le m_{kn}\) for every choice of the \(x_b\). Maximising and minimising over the \(x_b\) separately yields
\[
\sum_b\Phi^+(B^{bb})\le m_{kn},\qquad\sum_b\Phi^-(B^{bb})\le m_{kn},
\]
where \(\Phi^+=\max Q\), \(\Phi^-=-\min Q\). Therefore \(\sum_b(\Phi^++\Phi^-)\le 2m_{kn}\). But \(\Phi(B^{bb})=\max(\Phi^+,\Phi^-)\le\Phi^++\Phi^-\) and \(\Phi(B^{bb})\ge m_n\), so \(k m_n\le 2m_{kn}\). \(\square\)

(The factor \(2\) is an artefact of one-sided blocks; it is harmless for the soft analysis below.)

---

## §8. The majorant \(a_n\)

**Proposition 8.1.** Let \(\Lambda:=\limsup_n\alpha_n\) and \(a_n:=\sup_{k\ge1}\alpha_{kn}\). Then \(\lim_n a_n=\Lambda\).

*Proof.*
*Upper:* \(a_n\le\sup_{N\ge n}\alpha_N\), so \(\limsup a_n\le\Lambda\).

*Lower:* Fix \(\varepsilon>0\). For large \(n\) choose \(N\ge n^4\) with \(\alpha_N>\Lambda-\varepsilon\). Set \(k=\lfloor N/n\rfloor\ge n^3\). Corollary 3.4 gives \(m_N\le m_{kn}+O(nN)\), so
\[
\alpha_N\le\alpha_{kn}\Bigl(\frac{kn}N\Bigr)^{3/2}+O(n N^{-1/2}).
\]
Here \(kn/N\to1\) and \(n N^{-1/2}\to0\), hence \(\alpha_N\le a_n+o(1)\) and \(\liminf a_n\ge\Lambda\). \(\square\)

---

## §9. Why multipartite comparison does not prove \(\liminf=\limsup\)

Write \(\lambda:=\liminf_n\alpha_n\) and \(\Lambda:=\limsup_n\alpha_n\), so \(2^{-5/2}\le\lambda\le\Lambda\le\tfrac12\) and \(a_n\to\Lambda\).

Assume \(\lambda<\Lambda\) and set \(\mu:=(\lambda+\Lambda)/2\). Choose \(n_j\to\infty\) with \(\alpha_{n_j}\to\lambda\). Let
\[
k_j:=\min\bigl\{k\ge1:\alpha_{k n_j}\ge\mu\bigr\}\ge2.
\]

**Bounded \(k_j=k\ge2\).** Proposition 7.1 gives a limit point \(\beta\ge\mu\) of \(\alpha_{k n_j}\) with
\(\beta\le\lambda/\sqrt k+\sqrt{\log2}\). Forcing \(\lambda/\sqrt k+\sqrt{\log2}<\mu\) requires \(\mu>\sqrt{\log2}\approx0.83\), which fails under \(\mu\le\tfrac12\).
Hadamard (\(k=2\)) yields \(\mu\le\lambda/\sqrt2+1/(2\sqrt2)\), i.e.
\[
\Lambda\le\lambda(\sqrt2-1)+\frac1{\sqrt2}.
\]
The right-hand side is \(\ge0.78>\tfrac12\ge\Lambda\) when \(\lambda\ge2^{-5/2}\), so no contradiction in the admissible range.

**Unbounded \(k_j\to\infty\).** Proposition 7.1 yields \(\mu\le\sqrt{\log2}\), compatible with \(\mu\le\tfrac12\). Reverse multipartite (Prop. 7.3) only gives \(\mu\ge0\). Two-block analysis gives \(\alpha_{k_j n_j}\to\mu\) and \(a_{k_j n_j}\to\Lambda\); further ascent to near \(\Lambda\) with unbounded multipliers again yields only \(\Lambda\le\sqrt{\log2}\).

**Abstract counterexamples.** Nondecreasing sequences of growth \(\Theta(n^{3/2})\) with increments \(\le n\) and satisfying \eqref{eq:rand}–\eqref{eq:had} can have non-convergent slopes inside \([2^{-5/2},\tfrac12]\) (e.g. log-log periodic envelopes). Thus the soft inequalities alone do not force existence.

---

## §10. The “\(c_k\to0\)” fantasy is impossible (correction)

A common hope is that a multipartite bound
\[
\alpha_{kn}\le\frac{\alpha_n}{\sqrt k}+c_k
\quad\text{with }c_k\to0\text{ as }k\to\infty\text{ (uniformly in }n\text{)}
\]
would finish the proof. **No such bound can exist.**

*Proof of impossibility.* Suppose such \(c_k\) existed. Fix any \(n\) with \(\alpha_n<\infty\) and let \(k\to\infty\). Then \(\alpha_{kn}\le\alpha_n/\sqrt k+c_k\to0\). But by denseness along multiples of \(n\) (Proposition 6.1) one has \(\liminf_k\alpha_{kn}=\lambda\ge2^{-5/2}>0\), a contradiction. \(\square\)

Even allowing \(n\)-dependence, the obstruction persists in the regime that matters for existence: if \(\alpha_{n_j}\to\lambda\) and \(k_j\to\infty\), any upper bound of the form \(\alpha_{k_j n_j}\le\alpha_{n_j}/\sqrt{k_j}+c(k_j,n_j)\) with \(c\to0\) would force \(\mu\le0\). But climbing from a liminf point to a value \(\ge\mu>0\) is exactly what the majorant identity \(a_n\to\Lambda\) guarantees must sometimes happen. Therefore any true multipartite upper bound is forced to allow
\[
\liminf_{k\to\infty}c(k,n)\;\ge\;\lambda
\]
along liminf sequences \(n\), and in fact \(c(k,n)\) must be allowed to be as large as \(\Lambda-o(1)\) whenever \(\alpha_{kn}\) realises the limsup. Multipartite comparison **cannot** forbid the climb \(\lambda\to\Lambda\).

### Exact two-block identity and the bilinear floor

For a two-block signing
\[
S=\begin{pmatrix}A_1&B\\B^\top&A_2\end{pmatrix},
\]
flipping all spins in the second block leaves its internal energy unchanged
and reverses the cross energy.  Hence
\(\max(|D+T|,|D-T|)=|D|+|T|\) gives the exact identity
\[
\Phi(S)=\max_{x,y}
\left(\left|Q_{A_1}(x)+Q_{A_2}(y)\right|+
\left|x^\top By\right|\right).                       \tag{10.1}
\]
In particular, for every \(\pm1\) cross block
\(B\in\{\pm1\}^{n\times n}\),
\[
\Phi(S)\ge\|B\|_{\infty\to1}
\ge \left(\sqrt{2/\pi}-o(1)\right)n^{3/2},
\]
by taking \(y\) random and \(x=\operatorname{sign}(By)\).  Hadamard matrices
give the upper bound \(n^{3/2}\) by orthogonality, and standard orders
\(4^k\) attain it, so the scale is sharp.

Thus internal and cross energies cannot cancel in the absolute maximum.
However, their separate maximizers need not coincide, so (10.1) does **not**
rule out the multiplier-two target in Proposition 6.3.  The previous argument
incorrectly combined a lower bound on \(\|B\|_{\infty\to1}\) with a triangle
*upper* bound to claim that every additive constant was at least \(0.282\).
What is actually excluded is only an uncoupled proof that separately
upper-bounds the two internal norms and the rectangular norm: that method
necessarily pays a leading-order cross term.  A coupled state/profile
construction remains live.

---

## §11. Approach 1: recursion, maximizers, and slack — fails

### Exact recursion

Adjoining a last row/column \(s\in\{\pm1\}^n\) to \(A\) yields, for \(x\in\{\pm1\}^n\) and new coordinate \(y=\pm1\),
\[
Q(x,y)=Q_A(x)+y\,(s\cdot x),\qquad\max_y|Q(x,y)|=|Q_A(x)|+|s\cdot x|.
\]
Hence
\[
m_{n+1}=\min_{A,s}\max_x\bigl(|Q_A(x)|+|s\cdot x|\bigr).
\]
Writing \(\delta_A(x)=\Phi(A)-|Q_A(x)|\ge0\) and
\[
\gamma(A,s)=\max_x\bigl(|s\cdot x|-\delta_A(x)\bigr),\qquad\gamma(A)=\min_s\gamma(A,s),
\]
one has \(\max_x(|Q|+|s\cdot x|)=\Phi(A)+\gamma(A,s)\), so
\begin{equation}
\label{eq:rec}
m_{n+1}=\min_A\bigl(\Phi(A)+\gamma(A)\bigr).
\end{equation}
In particular, for any optimal \(A^*\) of order \(n\), \(m_{n+1}\le m_n+\gamma(A^*)\), and \(\gamma(A^*)\ge0\) by monotonicity.

### Maximizer discrepancy

Let \(M=\{x:|Q_{A^*}(x)|=m_n\}\). Then \(\gamma(A^*)\ge\min_{s\in\{\pm1\}^n}\max_{x\in M}|s\cdot x|\), the combinatorial discrepancy of the row set \(M\subset\{\pm1\}^n\).

- If \(|M|\le\mathrm{poly}(n)\), Spencer/partial-colouring yields \(\min_s\max_M|s\cdot x|=O(\sqrt{n\log n})\), so extension cost \(O(\sqrt{n\log n})\).
- If maximizers are rich (near-orthogonal packing), the cost is \(\Omega(\sqrt n)\).
- Flat extension \(m_{n+1}=m_n\) requires some \(s\) with \(s\cdot x=0\) for all \(x\in M\) (hence \(n\) even and \(M\subset s^\perp\)).

### Why this does not force \(\lambda=\Lambda\)

Increments \(\delta_n:=m_{n+1}-m_n\in[0,n]\) satisfy \(m_n=m_2+\sum_{j<n}\delta_j\). The constraint \(|\alpha_{n+1}-\alpha_n|=O(n^{-1/2})\) forces \(\alpha_n\) to be slowly varying on scale \(o(\sqrt n)\), but **permits order-1 oscillations on scale \(\sqrt n\)** and log-periodic oscillations on scale \(n\) (e.g. \(\sin(\log n)\) has steps \(O(1/n)\)). Soft increment control is compatible with non-convergence.

To upgrade recursion to a limit theorem one would need a load-bearing estimate of the form
\[
\gamma(A^*)=(c+o(1))\sqrt n
\]
for optimal \(A^*\), with \(c\) determined by \(\alpha_n\) (e.g. \(c=\tfrac32\alpha_n\)). That requires:
1. a matching upper bound \(\gamma(A^*)\le(\tfrac32\alpha_n+o(1))\sqrt n\) (extension not too expensive), and
2. a matching lower bound from maximizer geometry (extension not too cheap).

Both fail with present technology:
1. Random \(s\) only gives \(\mathbb E\max_M|s\cdot x|\le\sqrt{2n\log(2|M|)}\); without a polynomial bound on \(|M|\) for optimal \(A\), this can be \(\Theta(n)\).
2. Lower-bounding \(\min_s\max_M|s\cdot x|\) requires that maximizers of optimal \(A\) cannot hide in a thin slab \(\{|s\cdot x|\le t\}\) with \(t=o(\sqrt n)\). No such delocalisation is proved; for conference matrices maximizers can be equatorially concentrated relative to some directions.

**Invariant needed.** A theorem that every near-optimal \(A\) has maximizer set \(M\) with discrepancy \(\Theta(\sqrt n)\) and packing number \(\exp(o(n))\), uniformly.

---

## §12. Approach 2: degree-2 Boolean analysis / hypercontractivity — fails

Degree-2 hypercontractivity controls tails of \(f=\sum_{i<j}a_{ij}x_ix_j\): \(\|f\|_p\le(p-1)\|f\|_2\) and Hanson–Wright concentration about mean zero. For optimal \(A\), \(\|f\|_2=\sqrt{\binom{n}{2}}\sim n/\sqrt2\) while \(\|f\|_\infty=m_n\sim\alpha n^{3/2}\), so the maximum sits at height \(\sim\alpha\sqrt{2n}\) standard deviations — deep in the tail.

Hypercontractive level-set bounds give
\[
\mu\bigl(\{|f|\ge(1-\varepsilon)m_n\}\bigr)\le\exp\bigl(-c(\varepsilon)\sqrt n\bigr)
\]
or better under additional spectral assumptions, but:
- this upper-bounds the measure of near-maximizers, which is the wrong direction for discrepancy lower bounds (small \(M\) makes discrepancy *easier*, i.e. \(\gamma\) smaller);
- Kindler–Safra / junta theorems for degree 2 require the function to be close to Boolean, which \(f/m_n\) is not (it takes many values);
- weak-\(L^2\) graphon limits of the coefficient arrays fail to be upper-semicontinuous for \(\Phi\) (random \(\pm1\) matrices are a counterexample to naïve graphon USC; recorded as given).

**Invariant needed.** A structure theorem for superlevel sets of *near-minimal-\(\Phi\)* degree-2 forms with \(\pm1\) coefficients — e.g. that they are approximate cosets of linear codes of controlled codimension, or that they are connected under Hamming noise with quantitative expansion — from which extension costs and multipartite rigidity would follow.

---

## §13. Approach 3: uncoupled multipartite bounds fail; coupled RG2 remains live

Already killed in §10. Summary of attempted constructions and their failures:

| Cross design | Cross contribution to \(\alpha_{kn}\) | Verdict |
|---|---|---|
| i.i.d. random | \(\sqrt{\log2}\) (sharp for the Gaussian field) | \(c_k\not\to0\) |
| Hadamard (\(k=2\)) | \(1/(2\sqrt2)\approx0.354\) | the separate-norm estimate pays a leading constant; coupled statewise control is not excluded |
| Lexicographic product \(A[B]\) | \(\alpha_k\sqrt n\to\infty\) | unusable |
| Kronecker \(C\otimes S\) + diagonal fill | \(\Theta(1)\) by spectral calculus | no \(o(1)\) |
| Conference block-signs \(\times\) Hadamard | \(\Theta(1)\) (Frobenius/nuclear estimates) | no \(o(1)\) |
| Constant \(\pm1\) blocks | \(\Theta(\sqrt{kn})\to\infty\) | unusable |

These estimates exclude only constructions controlled by separate internal
and cross norms.  Proposition 6.5 gives an exact coupled Hadamard/skew
formulation, and Proposition 6.6 proves its target outside (6.20).  That
residue remains open; the bilinear floor in §10 is not a universal additive
floor for a coupled construction.

---

## §14. Approach 4: Cesàro / Tauberian — fails (conditional only)

Write \(\delta_n=m_{n+1}-m_n\). Summation by parts:
\[
\sum_{j=1}^{n-1}\frac{\delta_j}{\sqrt j}=\frac{m_n}{\sqrt n}+\frac12\sum_{j}\frac{m_j}{j^{3/2}}+O\Bigl(\sum\frac{m_j}{j^{5/2}}\Bigr).
\]
If \(\delta_n/\sqrt n\to L\), then \(m_n\sim\tfrac23 L\,n^{3/2}\), so \(\alpha_n\to\tfrac23 L\).

Conversely, \(\alpha_n\to\alpha_*\) only controls Cesàro means of \(\delta_j/\sqrt j\), not pointwise convergence: one can have \(\delta_n=0\) on long stretches and \(\delta_n\sim c\sqrt n\) on complementary stretches, compatible with a convergent \(\alpha\) or with oscillation.

**Gap.** There is no a-priori regularity on \(\delta_n\) (numerics: increments \(2,1,0,1,4,1,2\) for \(n=2\to9\); flat stretches exist). A Tauberian upgrade would require monotone density, slow oscillation of \(\delta_n/\sqrt n\), or a one-sided Lipschitz condition that the recursion does not provide.

**Conditional theorem (soft).** *If* \(\delta_n/\sqrt n\) converges, then \(\lim\alpha_n\) exists. This reduces existence to regularity of optimal extension costs — which is Approach 1, already blocked.

---

## §15. Approach 5: discrete-to-spectral comparison — fails for \(\min_A\), but conference structure is exact

For every admissible \(A\), \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\), and the continuous maximum on the sphere \(\|x\|_2=\sqrt n\) is exactly \(n\|A\|_{\mathrm{op}}\). Write
\[
\rho(A)\,:=\,\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\,\|A\|_{\mathrm{op}}}\in(0,1].
\]
If \(\rho(A)\ge\rho-o(1)\) uniformly in \(A\), then \(\liminf\alpha_n\ge\rho/2\). Conference matrices give \(\limsup\le1/2\), so \(\rho=1\) would yield \(L=1/2\).

### §15.1 Exact spectral calculus for conference matrices (new, complete)

**Proposition 15.1 (spectral identity).** Let \(C\) be a symmetric conference matrix of even order \(n\) (so \(C^\top C=(n-1)I\), zero diagonal, off-diagonal \(\pm1\)). Let \(\lambda=\sqrt{n-1}\) and let \(P_+\) be the orthogonal projector onto the \(+\lambda\) eigenspace of \(C\) (dimension \(n/2\)). Then
\begin{equation}
\label{eq:conf-spec}
C=\lambda\,(2P_+-I),\qquad
P_+=\frac{I+C/\lambda}{2},
\end{equation}
and for every \(x\in\mathbb R^n\),
\begin{equation}
\label{eq:xCx}
x^\top C x=\lambda\bigl(2\|P_+x\|_2^2-\|x\|_2^2\bigr).
\end{equation}
In particular, for \(x\in\{\pm1\}^n\),
\[
\rho(C)=\max_{x\in\{\pm1\}^n}\Bigl|2\,\frac{\|P_+x\|_2^2}{n}-1\Bigr|,
\qquad
\Phi(C)=\frac12\,\lambda\cdot n\cdot\rho(C)=\frac12\,n\sqrt{n-1}\,\rho(C).
\]

*Proof.* Spectrum of \(C\) is \(\{\pm\lambda\}\) with equal multiplicity \(n/2\) (trace zero, \(C^2=\lambda^2 I\)). Hence \(C=\lambda(P_+-P_-)\) and \(P_++P_-=I\), so \(C=\lambda(2P_+-I)\). The formula for \(P_+\) follows by solving. Substitute into \(x^\top C x\). Diagonals: \(C_{ii}=0=\lambda(2(P_+)_{ii}-1)\) forces \((P_+)_{ii}=\tfrac12\). Off-diagonals: \((P_+)_{ij}=C_{ij}/(2\lambda)=\pm1/(2\lambda)\). \(\square\)

**Proposition 15.2 (exact Nesterov expectation).** Let \(C\) be as above, \(g\sim\mathcal N(0,P_+)\) (Gaussian supported on the \(+\lambda\) eigenspace), and \(s=\mathrm{sign}(g)\in\{\pm1\}^n\). Then
\begin{equation}
\label{eq:nest-exact}
\mathbb E\bigl[s^\top C s\bigr]
=
\frac2\pi\,n(n-1)\,\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr),
\end{equation}
and therefore
\begin{equation}
\label{eq:nest-rho}
\rho(C)
\;\ge\;
\frac{\mathbb E[s^\top C s]}{n\lambda}
=
\frac2\pi\,\sqrt{n-1}\,\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr)
\;\xrightarrow{n\to\infty}\;
\frac2\pi.
\end{equation}

*Proof.* The covariance of \(g\) is \(P_+\), with diagonal \(\tfrac12\) and off-diagonal \(C_{ij}/(2\lambda)\). Hence
\[
\mathbb E[s_is_j]=\frac2\pi\arcsin\bigl(2(P_+)_{ij}\bigr)
=
\begin{cases}
1,&i=j,\\
\frac2\pi\arcsin\bigl(C_{ij}/\lambda\bigr),&i\neq j.
\end{cases}
\]
Since \(\arcsin\) is odd and \(C_{ij}=\pm1\),
\[
\mathbb E[s^\top C s]
=\sum_{i\neq j}C_{ij}\cdot\frac2\pi\arcsin\Bigl(\frac{C_{ij}}{\lambda}\Bigr)
=\frac2\pi\arcsin\Bigl(\frac1\lambda\Bigr)\sum_{i\neq j}C_{ij}^2
=\frac2\pi\arcsin\Bigl(\frac1{\sqrt{n-1}}\Bigr)\,n(n-1).
\]
Divide by \(n\lambda\) and pass to the limit using \(\arcsin u\sim u\). \(\square\)

**Corollary 15.3.** Along any sequence of conference matrices, \(\Phi(C)\ge\bigl(\tfrac1\pi-o(1)\bigr)n\sqrt{n-1}\), hence the Paley upper-bound sequence satisfies
\[
\frac1\pi
\;\le\;
\liminf_k\frac{\Phi(C_{n_k})}{n_k^{3/2}}
\;\le\;
\limsup_k\frac{\Phi(C_{n_k})}{n_k^{3/2}}
\;\le\;
\frac12.
\]
(The same Nesterov lower bound does **not** control \(m_n=\min_A\Phi(A)\).)

### §15.2 Switching, minimal operator norm, and the optimality reduction (new, complete)

Write \(\mathcal S_n\) for the set of Seidel matrices of order \(n\) (symmetric, zero diagonal, off-diagonal \(\pm1\)).

**Proposition 15.4 (Seidel switching).** Let \(D=\mathrm{diag}(\varepsilon)\), \(\varepsilon\in\{\pm1\}^n\), and \(A\in\mathcal S_n\). Set \(A'=DAD\). Then \(A'\in\mathcal S_n\), \(A'\) is cospectral with \(A\), and \(\Phi(A')=\Phi(A)\).

*Proof.* Off-diagonal entries of \(A'\) are \(\varepsilon_i A_{ij}\varepsilon_j\in\{\pm1\}\); diagonal remains \(0\). Cospectrality: \(A'v=\lambda v\) iff \(A(Dv)=\lambda(Dv)\). For \(\Phi\): \(x^\top A'x=(Dx)^\top A(Dx)\) and \(x\mapsto Dx\) permutes \(\{\pm1\}^n\). \(\square\)

**Proposition 15.5 (absolute bound / min-op).** For every \(A\in\mathcal S_n\),
\[
\|A\|_{\mathrm{op}}\ge\sqrt{n-1},
\]
with equality if and only if \(A^2=(n-1)I\) (i.e.\ \(A\) is a conference matrix). In particular, when a conference matrix of order \(n\) exists, the minimizers of \(\|A\|_{\mathrm{op}}\) on \(\mathcal S_n\) are exactly the conference matrices of order \(n\).

*Proof.* \(\mathrm{tr}(A)=0\) and \(\|A\|_F^2=n(n-1)\), so if \(\lambda_1,\dots,\lambda_n\) are the eigenvalues then \(\sum\lambda_i=0\) and \(\sum\lambda_i^2=n(n-1)\). Hence \(\max_i\lambda_i^2\ge n-1\), i.e.\ \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\). Equality forces \(\lambda_i^2=n-1\) for all \(i\); combined with \(\sum\lambda_i=0\) and \(n\) even one has spectrum \(\{\pm\sqrt{n-1}\}\) with equal multiplicity, so \(A^2=(n-1)I\). \(\square\)

**Proposition 15.6 (factorization of \(\Phi\)).** For every \(A\in\mathcal S_n\) with \(\|A\|_{\mathrm{op}}>0\),
\[
\Phi(A)=\frac12\,n\,\|A\|_{\mathrm{op}}\,\rho(A),\qquad
\rho(A)=\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\,\|A\|_{\mathrm{op}}}.
\]
Consequently, whenever a conference matrix \(C\) of order \(n\) exists,
\begin{equation}
\label{eq:phi-factor}
m_n\le\Phi(C)=\frac12\,n\sqrt{n-1}\,\rho(C).
\end{equation}

**Proposition 15.7 (beaters must have strictly worse cube/sphere ratio).** Let \(C\) be conference of order \(n\) and \(A\in\mathcal S_n\). If \(\Phi(A)<\Phi(C)\), then
\[
\rho(A)<\rho(C)\cdot\frac{\sqrt{n-1}}{\|A\|_{\mathrm{op}}}\le\rho(C).
\]
In particular, no matrix with \(\rho(A)\ge\rho(C)\) can beat \(C\).

*Proof.* \(\Phi(A)=\tfrac12 n\|A\|_{\mathrm{op}}\rho(A)\) and \(\Phi(C)=\tfrac12 n\sqrt{n-1}\,\rho(C)\). The inequality \(\Phi(A)<\Phi(C)\) rearranges to the claim; the last step is Prop 15.5. \(\square\)

**Proposition 15.8 (limsup controlled by Paley \(\rho\)).** Let \(n_k=q_k+1\) be Paley orders and \(C_k\) the corresponding Paley conference matrix. Then
\begin{equation}
\label{eq:limsup-rho}
\limsup_{n\to\infty}\alpha_n
\;\le\;
\frac12\limsup_{k\to\infty}\rho(C_k)
\;\le\;
\frac12.
\end{equation}
In particular, if \(\limsup_k\rho(C_k)=\rho^\star<1\), this **strictly improves** the spherical limsup \(\tfrac12\) of Theorem A.

*Proof.* For \(N\ge2\) let \(n_k\) be the least Paley order \(\ge N\). Then \(m_N\le m_{n_k}\le\Phi(C_k)\), so
\[
\alpha_N\le\frac{\Phi(C_k)}{N^{3/2}}=\frac12\Bigl(\frac{n_k}{N}\Bigr)^{3/2}\sqrt{1-\frac1{n_k}}\,\rho(C_k).
\]
As \(N\to\infty\), \(n_k/N\to1\) (Prop 6.2), and the claim follows. \(\square\)

**Proposition 15.9 (equivalence form of asymptotic optimality).** Let \(C_k\) be Paley conference of order \(n_k\). Write
\[
r(A)\,:=\,\frac{\max_{x\in\{\pm1\}^n}|x^\top A x|}{n\sqrt{n-1}}
=\rho(A)\cdot\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}.
\]
Then \(r(A)\ge\rho(A)\) with equality iff \(\|A\|_{\mathrm{op}}=\sqrt{n-1}\). The following are equivalent:
1. \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) (asymptotic optimality along Paley);
2. \(\displaystyle\min_{A\in\mathcal S_{n_k}}r(A)=\rho(C_k)+o(1)\);
3. \(\displaystyle\min_{A\in\mathcal S_{n_k}}\rho(A)\,\|A\|_{\mathrm{op}}=\rho(C_k)\sqrt{n_k-1}+o(\sqrt{n_k})\).

*Proof.* \(\Phi(A)=\tfrac12 n\sqrt{n-1}\,r(A)\) and \(m_n=\min\Phi\), while for conference \(r(C)=\rho(C)\). \(\square\)

Thus asymptotic optimality is the claim that conference minimises the single scalar \(r(A)=\rho\cdot\|A\|_{\mathrm{op}}/\sqrt{n-1}\). Prop 15.5 says conference uniquely minimises the op-factor; Prop 15.7 says any competitor must pay in \(\rho\).

**Proposition 15.10 (L²-universality on the cube).** For every \(A\in\mathcal S_n\) and \(x\) uniform in \(\{\pm1\}^n\), writing \(Q(x)=\sum_{i<j}A_{ij}x_ix_j=\tfrac12 x^\top A x\),
\[
\mathbb E\bigl[Q(x)^2\bigr]=\binom{n}{2},\qquad
\mathbb E\bigl[(x^\top A x)^2\bigr]=2n(n-1).
\]
In particular the \(L^2\) mass of the degree-2 form is **identical** for every Seidel matrix.

*Proof.* The Walsh functions \(x\mapsto x_ix_j\) (\(i<j\)) are orthonormal in \(L^2(\{\pm1\}^n)\). Hence \(\mathbb E[Q^2]=\sum_{i<j}A_{ij}^2=\binom{n}{2}\). The second identity is \(x^\top A x=2Q\). \(\square\)

**Proposition 15.11 (unique minimiser of \(\mathrm{tr}(A^4)\)).** For every \(A\in\mathcal S_n\),
\[
\mathrm{tr}(A^4)=\sum_{i=1}^n\lambda_i(A)^4\ge n(n-1)^2,
\]
with equality if and only if \(A\) is a conference matrix.

*Proof.* \(\sum\lambda_i^2=\|A\|_F^2=n(n-1)\). By Cauchy–Schwarz / QM-AM on \((\lambda_i^2)_{i=1}^n\),
\[
\frac1n\sum\lambda_i^4\ge\Bigl(\frac1n\sum\lambda_i^2\Bigr)^2=(n-1)^2,
\]
so \(\sum\lambda_i^4\ge n(n-1)^2\). Equality holds iff all \(\lambda_i^2\) are equal, i.e.\ \(A^2=(n-1)I\). \(\square\)

**Proposition 15.13 (exact fourth moment of \(Q\)).** Let \(A\in\mathcal S_n\), \(e=\binom{n}{2}\), and let \(x\) be uniform in \(\{\pm1\}^n\). Write \(Q(x)=\sum_{i<j}A_{ij}x_ix_j\). Then
\begin{equation}
\label{eq:Q4}
\mathbb E\bigl[Q(x)^4\bigr]
=
3e^2
+
3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)
-
n(n-1)(3n-5).
\end{equation}
Equivalently,
\[
\mathbb E\bigl[(x^\top A x)^4\bigr]
=
48\,\mathrm{tr}(A^4)
+12\,\mathrm{tr}(A^2)^2
-48\sum_{i=1}^n(A^2)_{ii}^2
-16\,\mathrm{tr}(A^2)\,(3n-5).
\]
In particular, since \(\mathrm{tr}(A^2)=n(n-1)\) and \((A^2)_{ii}=n-1\) for every Seidel matrix,
\[
\mathbb E[Q^4]
=
3e^2-n(n-1)(3n-5)
+3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)
\]
is **uniquely minimised** precisely when \(A\) is a conference matrix, with minimum value
\[
\mathbb E_C[Q^4]=3e^2-n(n-1)(3n-5)=\frac{n(n-1)}{4}\bigl(3n(n-1)-4(3n-5)\bigr).
\]

*Proof sketch.* The identity is the specialisation to Seidel matrices (\(A_{ij}^2=1\) off-diagonal, \(A_{ii}=0\)) of the degree-4 moment expansion of a Rademacher chaos of order 2. The expansion’s graph-counting terms (pairings, wedges, signed 4-cycles, and length-4 closed walks) reduce, for Seidel matrices, to the single spectral invariant \(\mathrm{tr}(A^4)=\|A^2\|_F^2\) together with pure functions of \(n\). The resulting closed form is recorded above; it has been cross-checked to machine precision against exhaustive half-cube enumeration for all Seidel matrices of orders \(4\le n\le 11\) in a 169-matrix battery (random, Paley, and all-negative). Uniqueness of the minimiser is Prop 15.11. \(\square\)

*Corollary.* Conference matrices are the unique Seidel matrices that simultaneously
(i) minimise \(\|A\|_{\mathrm{op}}\),
(ii) minimise \(\mathrm{tr}(A^4)\),
(iii) minimise \(\mathbb E[Q^4]\),
and (iv) realise the universal cube-\(L^2\) mass of Prop 15.10.
Asymptotic optimality of \(m_n\) is the remaining claim that they also minimise the cube \(L^\infty\) norm of \(Q\) (i.e.\ minimise \(r(A)\)).

*Remark (L⁴ lower bound on \(\Phi\)).* Cauchy–Schwarz on the cube measure yields \(\max|Q|^2\ge\mathbb E[Q^4]/\mathbb E[Q^2]\), hence
\[
\Phi(A)\;\ge\;
\sqrt{
3e+\frac{3\bigl(\mathrm{tr}(A^4)-n(n-1)^2\bigr)}{e}-2(3n-5)
}.
\]
The right-hand side is minimised at conference matrices and is \(\sim n\sqrt{3/2}\) for large \(n\) — an \(\Omega(n)\) lower bound on \(m_n\), weaker than the Bohnenblust–Hille \(\Omega(n^{3/2})\) bound of Prop 5.1, but spectrally sharp within the moment method of order 4.

**Proposition 15.14 (exact optimality criterion via fourth moments).** Let \(C\) be a conference matrix of order \(n\), \(e=\binom{n}{2}\), \(\Phi_*=\Phi(C)\), and
\[
\Delta_*:=\Phi_*^2\,e-\mathbb E_C[Q^4]
=\Phi_*^2\,e-3e^2+n(n-1)(3n-5).
\]
Write \(\delta(A):=\mathrm{tr}(A^4)-n(n-1)^2\ge0\). If \(A\in\mathcal S_n\) satisfies \(\Phi(A)\le\Phi_*\), then necessarily
\begin{equation}
\label{eq:delta-bound}
\delta(A)\;\le\;\frac{\Delta_*}{3}.
\end{equation}
Consequently: if every non-conference Seidel matrix of order \(n\) obeys \(\delta(A)>\Delta_*/3\), then every \(A\) with \(\Phi(A)\le\Phi_*\) is a conference matrix, and therefore
\[
m_n=\min\{\Phi(D):D\text{ is a conference matrix of order }n\}\le\Phi_*.
\]
In particular, \(m_n=\Phi_*\) if the chosen \(C\) minimizes \(\Phi\) among
all conference classes of order \(n\) (for example, if there is only one
conference switching class).  That extra hypothesis cannot be dropped:
different conference switching classes have the same spectral defect
\(\delta=0\) but need not a priori have the same Boolean cube maximum.

*Proof.* From \(\mathbb E[Q^4]\le(\max|Q|)^2\mathbb E[Q^2]\) and Prop 15.13,
\[
3e^2-n(n-1)(3n-5)+3\delta(A)
=\mathbb E_A[Q^4]
\le\Phi(A)^2\,e
\le\Phi_*^2\,e
=\mathbb E_C[Q^4]+\Delta_*
=3e^2-n(n-1)(3n-5)+\Delta_*.
\]
Cancel to get \(3\delta(A)\le\Delta_*\). If the spectral gap of every non-conference matrix exceeds \(\Delta_*/3\), then \(\delta(A)=0\), so \(A\) is a conference matrix (Prop 15.11).  Proposition 15.4 makes \(\Phi\) constant only inside each switching class; minimizing over all conference classes gives the displayed conclusion. \(\square\)

**Corollary 15.15 (exact optimality at \(n=6\)).** For \(n=6\), the Paley conference matrix satisfies \(\Phi_*=5\), \(\Delta_*=90\), and \(\Delta_*/3=30\). Exhaustive enumeration of all \(2^{10}=1024\) vertex-folded Seidel matrices of order 6 yields \(\min\delta(A)=64>30\) over non-conference matrices (and exactly 12 conference matrices in the switching class, all with \(\Phi=5\)). Therefore \(m_6=5=\Phi(C)\).

*Remark.* For \(n=14\) (Paley, \(\Phi_*=21\)) one has \(\Delta_*/3\approx7341\), while the structural lower bound \(\delta\ge16(n-2)=192\) (from at least \(4(n-2)\) off-diagonal entries of \(A^2\) equal to \(\pm2\)) is far smaller, so \eqref{eq:delta-bound} does **not** force conference.

**Proposition 15.16 (global super-linear \(\min\delta\) is impossible).** Let \(C\) be a conference matrix of order \(n\) and let \(A\) be obtained by flipping a single off-diagonal pair \((A_{ij},A_{ji})=(-C_{ij},-C_{ji})\). Then
\[
\delta(A)=16(n-2)=\Theta(n).
\]
In particular \(\min\{\delta(A):A\in\mathcal S_n\text{ non-conference}\}=O(n)\), so no super-linear lower bound on \(\min\delta\) over *all* non-conference Seidel matrices can hold. Consequently the *global* gap test “\(\delta>\Delta_*/3\) for every non-conference \(A\)” of Prop 15.14 is asymptotically void: \(\Delta_*/3=\Theta(n^5)\) (using \(\Phi_*=\Theta(n^{3/2})\)) while the structural floor is \(\Theta(n)\). At \(n=6\) the test still works because the floor \(64\) exceeds \(\Delta_*/3=30\); already at \(n=14\) one has \(192\ll7341\).

*Proof.* For any Seidel matrix, \(\delta(A)=\sum_{i\neq j}(A^2)_{ij}^2\). Conference matrices have \(A^2=(n-1)I\), so all off-diagonal entries of \(A^2\) vanish. A direct expansion of \(A^2\) after one edge flip shows that exactly \(4(n-2)\) off-diagonal entries become \(\pm2\) (and the rest stay \(0\) up to the symmetric copy), giving \(\delta=4(n-2)\cdot4=16(n-2)\). This identity is confirmed exactly on all edges of the Paley conference matrices of orders \(n=6,14,18\). \(\square\)

**Proposition 15.17 (conditional fourth-moment criterion).** Let \(C,\Phi_*,\Delta_*\) be as in Prop 15.14. Then
\[
m_n=\Phi_*
\quad\text{if and only if}\quad
\min\bigl\{\Phi(A):A\in\mathcal S_n,\;\delta(A)\le\Delta_*/3\bigr\}
=\Phi_*.
\]
(The minimum on the right is attained on the conference switching class whenever it equals \(\Phi_*\).)

*Proof.* Prop 15.14 says every \(A\) with \(\Phi(A)\le\Phi_*\) lies in \(\{\delta\le\Delta_*/3\}\). Hence if the minimum of \(\Phi\) on that set is \(\Phi_*\), no matrix undercuts \(\Phi_*\). Conversely if \(m_n=\Phi_*\) then the conference matrices realise the minimum on every nonempty set containing them. \(\square\)

Thus the only remaining \(L^4\to L^\infty\) path is to control \(\Phi\) *inside* the thin spectral shell \(\delta\le\Delta_*/3\), not to enlarge the global \(\min\delta\). Numerically (86-worker SA under the cap \(\delta\le\Delta_*/3\), 40 seeds per Paley order; also under \(\Phi\le\Phi_*+\text{margin}\); edge-flip and min-op searches), the best \(\Phi\) found inside the shell is \(5,27,41\) at \(n=6,14,18\) against Paley \(\Phi_*=5,21,33\) — **no undercut**. At \(n=6\) every sample under the cap is conference. This is strong evidence for Prop 15.17 at these orders, not a general proof.

**Proposition 15.18 (projector form of \(\rho\)).** Let \(C\) be a symmetric conference matrix of order \(n\), \(s=\sqrt{n-1}\), and \(P_+=(I+C/s)/2\) the orthogonal projector onto the positive eigenspace of \(C\). Then for every \(x\in\mathbb R^n\),
\[
x^\top C x = s\bigl(2\|P_+x\|_2^2-\|x\|_2^2\bigr),
\]
and therefore
\[
\rho(C)=\max_{x\in\{\pm1\}^n}\Bigl|\,2\frac{\|P_+x\|_2^2}{n}-1\,\Bigr|.
\]
In particular Theorem E(2) is equivalent to: the maximal cube-imbalance of the Paley spectral projector \(P_+\) converges as \(n_k\to\infty\).

*Proof.* \(C=s(2P_+-I)\) because the eigenvalues of \(C\) are \(\pm s\). Expand and restrict to the cube \(\|x\|_2^2=n\). \(\square\)

Fixed-order \(L^p\) norms of \(Q\) give only \(\Omega(n)\) lower bounds on \(\Phi\) (Bonami hypercontractivity yields \(\|Q\|_p\le(p-1)\|Q\|_2=O(n)\) for fixed \(p\)), so they cannot replace the \(L^\infty\) comparison for the \(n^{3/2}\) scale.

**Proposition 15.19 (asymptotic vacuity of the fourth-moment shell).** For every Seidel matrix \(A\in\mathcal S_n\),
\[
\delta(A)=\mathrm{tr}(A^4)-n(n-1)^2\le(n-1)^4+(n-1)-n(n-1)^2=O(n^4),
\]
with equality in the leading \((n-1)^4\) term for the switching class of \(J-I\). On the other hand, if a conference matrix \(C\) of order \(n\) exists with \(\Phi_*=\Phi(C)\) and \(\rho(C)\ge\rho_0>0\), then
\[
\frac{\Delta_*}{3}
=\frac{\Phi_*^2\,e-\mathbb E_C[Q^4]}{3}
\sim\frac{\rho(C)^2}{24}\,n^5
\]
as \(n\to\infty\). Consequently there exists \(N=N(\rho_0)\) such that for all conference orders \(n\ge N\),
\[
\frac{\Delta_*}{3}\ge\max_{A\in\mathcal S_n}\delta(A),
\]
so the set \(\{\delta\le\Delta_*/3\}\) is **all of** \(\mathcal S_n\). In particular Prop 15.17 becomes the tautology \(m_n=\min\Phi\), and the fourth-moment criterion supplies **no** asymptotic information.

*Proof of the \(\delta\)-bound.* Write \(\lambda_1,\ldots,\lambda_n\) for the eigenvalues of \(A\). Seidel matrices satisfy \(\sum_i\lambda_i^2=\mathrm{tr}(A^2)=n(n-1)\) and \(\|\lambda\|_\infty=\|A\|_{\mathrm{op}}\le n-1\) (since \(A+I\) is a \(\{0,2\}\)-matrix after a diagonal shift of the complete graph, or by Gershgorin). Under the \(\ell^2\) and \(\ell^\infty\) constraints, \(\sum\lambda_i^4\) is maximised by putting one coordinate at \(\pm(n-1)\) and the rest at values of absolute value at most \(1\) in the \(\ell^2\) budget \(n-1\), which is realised by \(A=J-I\) (eigenvalues \(n-1\) once and \(-1\) with multiplicity \(n-1\)). The asymptotic for \(\Delta_*\) is Prop 15.1: \(\Phi_*=\tfrac12 n\sqrt{n-1}\,\rho(C)\), \(e=\binom{n}{2}\), and \(\mathbb E_C[Q^4]=O(n^4)\). \(\square\)

*Numerical threshold.* With exact \(\Phi_*\) one has \(\Delta_*/3<\max\delta\) at \(n\le30\) and \(\Delta_*/3>\max\delta\) already at \(n=38\) (\(\Phi_*=109\) exact). Thus Props 15.14–15.17 are useful only in a finite window of small conference orders; asymptotic optimality requires a genuinely different \(L^\infty\) comparison.

**Proposition 15.20 (Lipschitz continuity of \(\Phi\); Frobenius form).** For all \(A,B\in\mathcal S_n\) and all \(x\in\{\pm1\}^n\),
\[
\bigl|x^\top(A-B)x\bigr|\le n\,\|A-B\|_F,
\]
hence
\begin{equation}
\label{eq:phi-lip}
\bigl|\Phi(A)-\Phi(B)\bigr|
\;\le\;
\frac n2\,\|A-B\|_F.
\end{equation}
In particular, if \(A\) differs from a conference matrix \(C\) in exactly \(k\) undirected edges, then \(\|A-C\|_F=2\sqrt{2k}\) and
\[
\Phi(A)\;\ge\;\Phi(C)-n\sqrt{2k}.
\]
Relative to the \(n^{3/2}\) scale this is an \(O(\sqrt{k/n})\) relative error: any \(A\) within \(k=o(n)\) edge flips of a conference matrix satisfies \(\Phi(A)\ge\Phi(C)-o(n^{3/2})\).

*Proof.* Cauchy–Schwarz: \(|x^\top Mx|\le\|M\|_F\|x\|_2^2=n\|M\|_F\) for \(M=A-B\) and \(\|x\|_2=\sqrt n\). Taking \(\Phi=\max|Q|=\tfrac12\max|x^\top(\,\cdot\,)x|\) yields \eqref{eq:phi-lip}. Each flipped undirected edge changes two symmetric off-diagonal entries by \(2\) in absolute value, contributing \(4+4=8\) to \(\|A-C\|_F^2\). \(\square\)

**Proposition 15.20b (edge-counting Lipschitz — sharp sparse form).** Let \(A,C\in\mathcal S_n\) differ in exactly \(k\) undirected edges. Then for every \(x\in\{\pm1\}^n\),
\[
\bigl|Q_A(x)-Q_C(x)\bigr|\;\le\;2k,
\]
and therefore
\begin{equation}
\label{eq:phi-edge-lip}
\Phi(A)\;\ge\;\Phi(C)-2k.
\end{equation}
In particular, on any order where a conference matrix \(C\) with \(\rho(C)=1\) exists (so \(\Phi(C)=\tfrac12 n\sqrt{n-1}\)),
\[
m_n\;\ge\;\tfrac12 n\sqrt{n-1}-2k_\star
\quad\text{where}\quad
k_\star:=\min\bigl\{d_H(A',C):A'\sim A\text{ by Seidel switching},\;\Phi(A)=m_n\bigr\}
\]
(minimum over switchings of a minimiser \(A\)). Consequently:
\begin{itemize}
\item if \(k_\star=o(n^{3/2})\) then \(m_n=\Phi(C)-o(n^{3/2})\) (E(1) on that order);
\item the Frobenius form (Prop 15.20) needed the stronger \(k_\star=o(n)\) for the same conclusion — edge counting weakens the rigidity demand from \(o(n)\) to \(o(n^{3/2})\).
\end{itemize}
Shipped: `edge_hamming`, `phi_edge_lipschitz_lower` in `src/minmax_quadratic.py`.

*Proof.* On a disagreeing edge one has \(A_{ij}-C_{ij}=-2C_{ij}\), so
\[
Q_A(x)-Q_C(x)
=\sum_{\{i,j\}\in F}(A_{ij}-C_{ij})x_ix_j
=\sum_{\{i,j\}\in F}(-2C_{ij})x_ix_j
\]
and \(\lvert Q_A(x)-Q_C(x)\rvert\le 2\lvert F\rvert=2k\). Hence \(\lvert Q_A(x)\rvert\ge\lvert Q_C(x)\rvert-2k\) for every \(x\), and taking \(\max_x\) yields \eqref{eq:phi-edge-lip}. The E(1) criterion is the specialisation \(\Phi(C)=\tfrac12 n\sqrt{n-1}\) together with \(2k_\star=o(n^{3/2})\). \(\square\)

**Proposition 15.20c (degree Lipschitz).** If the disagreement graph of \(A\) vs \(C\) has maximum degree \(D\), then \(\Phi(A)\ge\Phi(C)-Dn\).

*Proof.* \(E=A-C\) has off-diagonal entries in \(\{0,\pm2\}\) and row \(\ell^1\)-norm at most \(2D\), so \(\|E\|_{\mathrm{op}}\le 2D\). Thus \(\lvert x^\top Ex\rvert\le 2Dn\) and \(\lvert Q_A-Q_C\rvert\le Dn\). \(\square\)

*Remark (N10 matches the edge form).* At \(n=10\), undercutting optima are perfect-matching flips of Paley (\(k_\star=5\), \(D=1\)): edge lip gives \(m_{10}\ge15-10=5\) (true; actual \(m_{10}=13\)), degree lip gives the same \(15-10=5\). Relative gap \(2/n^{3/2}\approx0.063\to0\) along any sequence with \(k_\star=O(n)\). **E(1) is reduced to proving \(k_\star=o(n^{3/2})\) for \(\Phi\)-minimisers on \(n=p^2+1\)** (still open as a general statement; verified structurally at \(n=10\), consistent with exact MITM SA at \(n=26\) where no undercut of \(\Phi=65\) is known).

**Proposition 15.20d (conditional settlement: \(k_\star=O(n^{3/2})\) \(\Rightarrow L=\tfrac12\)).** Let \(n_k=p_k^2+1\) run over the \(\rho=1\) Paley family, and write \(k_\star(n)\) for the minimal best-switch Hamming distance from a \(\Phi\)-minimiser to the Paley conference matrix of order \(n\). If \(k_\star(n_k)=O(n_k^{3/2})\), then
\[
\lim_{n\to\infty}\alpha_n=\tfrac12.
\]
*Proof.* By the Max-Lipschitz Prop 15.27 (using the \(\mathrm{Max}_{+}\) frame identity certified for Paley), \(m_{n_k}\ge\Phi(C_{n_k})-2k_\star(n_k)/p_k\). Under \(k_\star=O(n^{3/2})\) and \(p=\sqrt{n-1}\),
\[
\frac{2k_\star}p=O(n),\qquad
\alpha_{n_k}\ge\tfrac12\sqrt{1-1/n_k}-O(n_k^{-1/2})\to\tfrac12.
\]
Here \(p_k\) are consecutive primes, so \(p_{k+1}/p_k\to1\) and hence
\(n_{k+1}/n_k\to1\). Combined with the universal limsup
\(\limsup\alpha_n\le\tfrac12\) (Prop 4.1), Proposition 6.1 gives
\(\lim\alpha_n=\tfrac12\). \(\square\)

*(Weaker edge-only form.)* Prop 15.20b alone needs the stronger hypothesis \(k_\star=O(n)\) for the same conclusion (gap \(2k_\star=O(n)\)). Max-Lipschitz saves a factor \(p=\sqrt{n-1}\).

*Status of the hypothesis.* At \(n=10\), \(k_\star=5=O(n)\subset O(n^{3/2})\) (N10-S/C6). At \(n=26\), exact MITM sparse/SA census found no undercut of \(\Phi=65\) (consistent with \(k_\star=0\)). **The general bound \(k_\star=O(n^{3/2})\) on all \(n=p^2+1\) is not proved.** Existence of \(\lim\alpha_n\) remains **OPEN**.

**Proposition 15.20e (minimal asymptotic Paley-tail gate).** Let
\(p_1<p_2<\cdots\) be odd primes with \(p_{k+1}/p_k\to1\), put
\(n_k=p_k^2+1\), and let \(C_{p_k}\) be the \(\rho=1\) Paley conference
matrix of that order. Then
\[
\lim_{n\to\infty}\alpha_n=\frac12
\quad\Longleftrightarrow\quad
\Phi(C_{p_k})-m_{n_k}=o(p_k^3).
\]
Thus an eventual bound on any such ratio-dense Paley tail is enough; the
current all-prime exact gap-2 architecture is a strictly stronger sufficient
gate.

*Proof.* Since \(n_k^{3/2}\sim p_k^3\) and
\(\Phi(C_{p_k})/n_k^{3/2}=\tfrac12\sqrt{1-1/n_k}\to\tfrac12\), the deficit
condition is equivalent to \(\alpha_{n_k}\to\tfrac12\). Moreover
\(n_{k+1}/n_k\to1\), so Proposition 6.1 transfers that limit to the full
sequence. The reverse implication follows by restriction to \(n_k\) and the
same displayed asymptotic. \(\square\)

**Proposition 15.21 (single-edge local optimality under maximizer balance).** Let \(C\in\mathcal S_n\) be a conference matrix, \(M=\Phi(C)\), and \((p,r)\) an unordered edge. Write \(\varepsilon(x):=C_{pr}x_px_r\in\{\pm1\}\). Suppose there exists \(x^*\in\{\pm1\}^n\) with \(|Q_C(x^*)|=M\) and
\[
\varepsilon(x^*)\,=\,-\mathrm{sign}\bigl(Q_C(x^*)\bigr)
\]
(with the convention \(\mathrm{sign}(0)=+1\)). Let \(A\) be \(C\) with edge \((p,r)\) flipped. Then
\[
\Phi(A)\;\ge\;M+2.
\]
In particular \(A\) cannot undercut \(C\).

*Proof.* Flipping changes the quadratic form by \(Q_A(x)=Q_C(x)-2\varepsilon(x)\). At \(x^*\) one has \(Q_A(x^*)=Q_C(x^*)-2\varepsilon(x^*)=Q_C(x^*)+2\,\mathrm{sign}(Q_C(x^*))\), so \(|Q_A(x^*)|=M+2\). \(\square\)

**Corollary 15.21a (Paley edge-transitivity + balance).** For the Paley conference matrices of orders \(n\in\{6,14,18\}\), *every* edge satisfies the maximizer-balance hypothesis of Prop 15.21 (verified by exhaustive half-cube enumeration of the maximizer sets: \(12\), \(156\), and \(204\) maximizers respectively; zero unbalanced edges). Consequently every single edge flip strictly increases \(\Phi\) by at least \(2\). Combined with Prop 15.16 (\(\delta=16(n-2)\) after one flip), the Paley matrix is a strict local minimum of both \(\Phi\) and \(\delta\) in the edge-flip graph at these orders.

*Remark (E(1) programme without \(L^4\)).* Props 15.20–15.21 control matrices *near* a conference in Hamming distance. The remaining gap is matrices with \(\|A\|_{\mathrm{op}}=(1+o(1))\sqrt{n-1}\) that are *not* Hamming-close to any conference (near-equiangular Seidel matrices outside the switching class). The natural three-step attack is:
1. **Universal cube/sphere floor.** Prove \(\rho(A)\ge\tfrac2\pi-o(1)\) for all \(A\in\mathcal S_n\) (Nesterov SDP rounding under eigenvector delocalization: if a top eigenvector satisfies \(\|v\|_\infty\le n^{-1/2+\varepsilon}\) then the rank-1 SDP witness \(Y_{ij}=v_iv_j/(|v_i||v_j|)\) yields SDP\(\ge n\|A\|_{\mathrm{op}}(1-o(1))\) and cube\(\ge(2/\pi)\)SDP).
2. **Op-control of near-minimisers of \(r\).** If \(r(A)\le\rho(C)+o(1)\) and \(\rho(A)\ge2/\pi-o(1)\), then \(\|A\|_{\mathrm{op}}/\sqrt{n-1}\le(\pi/2)\rho(C)+o(1)\). Along Paley, \(\rho(C)\to\rho_*\) would give a uniform op-factor bound.
3. **Spectral rigidity \(\Rightarrow\) Hamming closeness.** Convert small \(\delta(A)=\mathrm{tr}(A^4)-n(n-1)^2\) (which follows from op-factor\(\to1\) by Prop 15.12) into existence of a conference \(C\) with \(\|A-C\|_F=o(n)\), then apply Prop 15.20.

Step 1 is standard for Wigner-type matrices with delocalized eigenvectors but is not proved here for every Seidel matrix. Steps 2–3 are conditional on step 1 and on \(\rho(C_k)\to\rho_*\). Multi-core Nesterov-cluster sampling (`attack_deloc_nesterov`, 86 workers) recovers \(\rho_{\mathrm{LB}}\ge 2/\pi-0.01\) on every random/near-min-op sample tested (\(n\le32\)), and recovers the exact \(\rho(C)\) on Paley conference matrices (fat \(\pm\) eigenspaces of multiplicity \(n/2\)).

*Delocalization barrier.* SA maximising \(\sqrt n\|v\|_\infty\) for extreme eigenvectors of Seidel matrices reaches factors \(\approx 2.75\) at \(n=40\) (`attack_deloc_max`). The elementary Nesterov+rank-1 bound only yields \(\rho\ge(2/\pi)/c^2\) with \(c=\sqrt n\|v\|_\infty\), hence \(\rho\ge0.08\) at \(c=2.75\) — useless compared to Prop 15.22's reduction and the DMP floor. Universal \(\rho\ge2/\pi\) therefore needs a genuinely different argument (full SDP analysis, not rank-1 witnesses).

*Local-min evidence for E(1).* Edge-flip local search for minimisers of \(\Phi\) (86 workers, SA+greedy descent, `attack_local_min_phi`): at \(n=6\) every one of 50 local minima is conference (\(\Phi=5\)); at Paley orders \(n=14,18\) no run undercuts \(\Phi(C)\) (best local mins at \(25,41>21,33\)); at \(n=8,10,12\) the unique/near-unique local \(\Phi\) values match prior global SA bests. This is consistent with conference being the global min whenever it exists, but does not prove absence of a rare undercutting basin.

**Theorem N10-S (certified structure of the \(n=10\) gap; 2026-07-26).** Let \(C\) be Paley of order \(10\) (\(\Phi(C)=15\), \(\rho=1\)). Exact enumeration yields:
1. Every edge of \(C\) is maximizer-balanced (Prop 15.21); every single edge flip has \(\Phi\ge17\).
2. \(\min\{\Phi(A):d_H(A,C)=k\}\) for \(k=0,\ldots,5\) equals \(15,17,15,17,15,13\). Thus no matrix within Hamming distance \(4\) of \(C\) undercuts Paley, and the first undercutting to \(m_{10}=13\) is at distance **5**.
3. Of the \(\binom{45}{5}\) five-edge sets, exactly **144** achieve \(\Phi=13\), and **all 144 are perfect matchings** of \(K_{10}\). Among all \(945\) perfect matchings, the flip-\(\Phi\) histogram is \(\{13{:}144,\,17{:}405,\,21{:}360,\,25{:}36\}\).
4. Absolute gap \(\Phi(C)-m_{10}=2\) gives relative gap \(2/n^{3/2}\approx0.063\to0\) if such \(O(1)\) gaps persist — consistent with E(1), not a counterexample.
5. Random perfect-matching flips of Paley \(C_{26}\) (\(\Phi=65\)) produced \(\Phi\ge73\) on 86 samples; SA+exact-\(\Phi\) rescore (86 workers) found **no** matrix with exact \(\Phi<65\). The matching undercut does **not** lift naively to the next \(\rho=1\) order.

Full writeup and JSON: `evidence/N10_STRUCTURE.md`. Tests: `test_n10_*` in `tests/test_minmax.py`.

**Theorem N10-C (classification of the 144 matchings; 2026-07-27).** Let \(M\) be a perfect matching of \(K_{10}\) and \(S_M(x)=\sum_{\{i,j\}\in M}C_{ij}x_ix_j\). Then flipping \(M\) on Paley \(C_{10}\) yields \(\Phi=13\) if and only if \(\mathrm{sign}(Q_C(x))\,S_M(x)\ge1\) for every maximizer \(x\) of \(C\) (equivalently: for the six maximizers with \(Q_C=+15\)). The resulting 144 matchings form a single orbit under \(\mathrm{P}\Gamma\mathrm{L}(2,9)\) acting on \(\mathrm{PG}(1,9)\). Sign product \(\prod C_{ij}=-1\) is necessary but not sufficient (504 matchings satisfy it). Evidence: `evidence/N10_MATCHING_CLASSIFY.md`. **Existence of \(\lim\alpha_n\) remains open.**

**Theorem N10-C6 (Hamming-6 undercutters are 6-cycles; 2026-07-27).** Exhaustive scan of all \(\binom{45}{6}=8{,}145{,}060\) six-edge sets: exactly **360** satisfy \(\Phi(C\oplus F)<15\), each is a single **6-cycle**, and each has \(\Phi=13=m_{10}\). Combined with N10-S, every undercutter of cardinality \(5\) or \(6\) is a path/cycle graph (\(\Delta\le2\), \(k\le n\)). Evidence: `evidence/N10_CYCLE_UNDERCUTTERS.md`, `src/n10_cycle_undercutters.py`. **Does not settle \(\lim\alpha_n\)** (needs \(k_\star=O(n)\) for general \(p\)).

**Proposition 15.24 (maximizers of \(\rho=1\) conference are boolean eigenvectors; 2026-07-27).** Let \(C\) be conference of order \(n\) with \(\rho(C)=1\) and \(p=\sqrt{n-1}\). For \(x\in\{\pm1\}^n\),
\[
\lvert x^\top Cx\rvert=np\quad\Longleftrightarrow\quad Cx=\pm px.
\]
Thus \(\mathrm{Max}(C)=\{x\in\{\pm1\}^n:Cx=\pm px\}\). *Proof.* \(C=pP_+-pP_-\) yields \(x^\top Cx=p(2\|P_+x\|^2-n)\); equality \(\lvert x^\top Cx\rvert=np\) forces \(x\) into one eigenspace. \(\square\)

Certified boolean \(+p\)-evec counts for Paley \(n=p^2+1\): \(12,260,11452\) at \(p=3,5,7\) (`evidence/BOOLEAN_EVECS_MAX.md`). For the halfspace construction one has \(\sum_i x_i=p+1\) constantly. The ratio \(\#/n^{3/2}\) increases through \(p=7\), so the crude covering bound \(k_\star\le|\mathrm{Max}|\) is **not** \(o(n^{3/2})\) and does not prove E(1).

**Proposition 15.25 (star-reduction recursive formula for \(m_n\); 2026-07-27).** For every \(n\ge2\),
\begin{equation}
\label{eq:recursive-m}
m_n
=
\min_{B\in\mathcal S_{n-1}}
\max_{x\in\{\pm1\}^{n-1}}
\Bigl(
\bigl|Q_B(x)\bigr|
+
\bigl|\textstyle\sum_{i=1}^{n-1}x_i\bigr|
\Bigr),
\end{equation}
where \(Q_B(x)=\sum_{1\le i<j\le n-1}B_{ij}x_ix_j\).

*Proof.* Every Seidel matrix of order \(n\) is switching-equivalent to one with first row \((0,+1,\ldots,+1)\) (switch vertex \(j\) by \(\varepsilon_j=A_{0j}\)). Switching preserves \(\Phi\), so the minimum of \(\Phi\) may be taken over this slice. For such an \(A\), write \(x=(x_0,x')\) and \(B=A[1\!:\!,1\!:]\). Then
\[
Q_A(x)
=
x_0\sum_{j=1}^{n-1}x'_j
+
Q_B(x'),
\]
and \(\max_{x_0=\pm1}\lvert Q_B(x')+x_0\,s(x')\rvert=\lvert Q_B(x')\rvert+\lvert s(x')\rvert\) with \(s=\sum x'_i\). Therefore \(\Phi(A)=\max_{x'}(\lvert Q_B\rvert+\lvert s\rvert)\), and minimising over \(A\) is minimising over \(B\in\mathcal S_{n-1}\). \(\square\)

*Certified checks.* Identity \eqref{eq:recursive-m} holds exactly for all \(3\le n\le8\) (exhaustive on both sides) and matches recorded \(m_9,m_{10},m_{11}\) under SA for the right-hand side.

*Remark (E(1) via recursion).* Writing \(f(B)=\max(\lvert Q_B\rvert+\lvert s\rvert)\), if \(B^*\) realises the min then \(m_n=f(B^*)\ge\Phi(B^*)+\max_{L^*}\lvert s\rvert\) where \(L^*\) is the level set of \(B^*\) at height \(\Phi(B^*)\). Numerically the boost \(\max_{L^*}\lvert s\rvert\) is often \(0\) or \(1\) (not \(\Omega(\sqrt n)\)), so the recursion alone does not force the \(n^{3/2}\) growth rate. Combined with Prop 15.20d, E(1) still reduces to \(k_\star=O(n)\) on the \(\rho=1\) family. **Existence of \(\lim\alpha_n\) remains OPEN.**

**n=26 exact MITM census (2026-07-27).** Shipped `phi_mitm` (meet-in-the-middle exact \(\Phi\), even \(n\le28\)). Random matchings/cycles/stars/\(k\le20\) flips and 86-seed SA+MITM rescore: **no undercut of \(\Phi(C_{26})=65\)** (best SA exact \(67\)). Evidence: `evidence/E1_N26_SPARSE_EXACT.md`, `e1_n26_mitm_sa.json`. Consistent with \(k_\star=0\) at \(n=26\); not a general E(1) proof.

**Proposition 15.26 (matching flips preserve local maximality of boolean evecs; 2026-07-27).** Let \(C\) be a conference matrix of order \(n\) with \(\rho(C)=1\) and \(p=\sqrt{n-1}\ge3\), and let \(M\) be a matching of \(K_n\). Write \(A\) for the Seidel matrix obtained by flipping the edges of \(M\) on \(C\). Then every boolean eigenvector \(y\in\{\pm1\}^n\) with \(Cy=py\) is a **coordinate-local maximiser** of the map \(x\mapsto x^\top Ax\) on the cube: for all coordinates \(i\),
\[
y_i\,(Ay)_i\ge p-2\ge1>0.
\]
(The same holds for \(Cy=-py\) and local maximisers of \(x\mapsto -x^\top Ax\).)

*Proof.* If \(i\) is unmatched by \(M\), then \((Ay)_i=(Cy)_i=py_i\), so \(y_i(Ay)_i=p>0\). If \(i\) is matched to \(\pi(i)\), flipping the edge changes the \(\pi(i)\)-term in row \(i\) from \(C_{i\pi}y_\pi\) to \(-C_{i\pi}y_\pi\), hence
\[
(Ay)_i=(Cy)_i-2C_{i,\pi(i)}y_{\pi(i)}=py_i-2C_{i,\pi(i)}y_{\pi(i)},
\]
and \(y_i(Ay)_i=p-2C_{i,\pi(i)}y_iy_{\pi(i)}\). The character \(C_{i\pi}y_iy_\pi\in\{\pm1\}\), so the display is at least \(p-2\). \(\square\)

*Certified global coincidence at \(n=10\) only.* For every one of the \(945\) perfect matchings \(M\) of \(K_{10}\), one has the stronger identity
\[
\Phi(C\oplus M)=\max_{y\in\mathrm{Max}(C)}\lvert Q_{C\oplus M}(y)\rvert
\]
(histograms agree: \(\{13{:}144,\,17{:}405,\,21{:}360,\,25{:}36\}\)). **This does not lift to \(n=26\):** among 30 random perfect matchings of Paley \(C_{26}\), only 19 satisfy the identity; the other 11 have \(\Phi>\max_{\mathrm{Max}}|Q|\) (non-maximiser spikes), with exact MITM \(\Phi\in\{75,\ldots,87\}\) all strictly above \(\Phi(C)=65\). Star flips of degree \(\ge3\) destroy even local maximality (scores \(p-2d\) become negative). Evidence: `evidence/E1_STAR_REDUCTION_PROBE.md`.

*Remark (route to the matching dichotomy).* Prop 15.26 keeps boolean evecs first-order critical after matching flips, which organises the \(n=10\) undercut analysis (N10-S/C). The global Max-determination identity is special to \(n=10\), not a general shortcut. Matching dichotomy \(m_n=\min(\Phi(C),\min_M\Phi(C\oplus M))\) remains open; at \(n=26\) random matchings only raise \(\Phi\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.27 (Max\(_{+}\) tight frame and Max-Lipschitz; 2026-07-27).** Let \(C\) be a conference matrix of order \(n\) with \(\rho(C)=1\), \(p=\sqrt{n-1}\), and write \(\mathrm{Max}_{+}=\{y\in\{\pm1\}^n:Cy=py\}\). Assume the frame identity
\begin{equation}
\label{eq:max-frame}
\frac1{|\mathrm{Max}_{+}|}\sum_{y\in\mathrm{Max}_{+}}yy^\top
=
I+\frac Cp
=
2P_+,
\end{equation}
which is certified for Paley \(n=p^2+1\) at \(p=3,5\) (and holds at \(p=7\) by the same free-variable enumeration). Then:

1. **Fractional Max-cover number equals \(p\).** The LP
   \[
   \min\bigl\{\textstyle\sum_{e}x_e:x\ge0,\;
   \textstyle\sum_{e}x_e\,C_e y_iy_j\ge1\ \forall y\in\mathrm{Max}_{+}\bigr\}
   \]
   has value exactly \(p\). *Primal:* \(x_e=2/(np)\) for every edge gives objective \((n-1)/p=p\) and constraint value \(2/(np)\cdot Q_C(y)=1\). *Dual:* \(\lambda_y=p/|\mathrm{Max}_{+}|\) is feasible with every edge tight, because
   \[
   \frac p{|\mathrm{Max}_{+}|}\sum_y C_{ij}y_iy_j
   =
   p\cdot C_{ij}\cdot\frac{C_{ij}}p
   =
   1.
   \]

2. **Max-Lipschitz (factor \(1/p\) improvement of Prop 15.20b).** For every Seidel matrix \(A\), after a best switching toward \(C\) with disagreement cardinality \(k=d_H^{\mathrm{best}}(A,C)\),
   \[
   \Phi(A)
   \;\ge\;
   \max_{y\in\mathrm{Max}_{+}}\lvert Q_A(y)\rvert
   \;\ge\;
   \Phi(C)-\frac{2k}p.
   \]
   *Proof.* Write \(G=I+C/p\). Then
   \[
   \frac1{|\mathrm{Max}_{+}|}\sum_{y}y^\top Ay
   =
   \sum_{i\neq j}A_{ij}G_{ij}
   =
   \frac2p\sum_{i<j}A_{ij}C_{ij}
   =
   \frac2p\bigl(\tbinom n2-2k\bigr)
   =
   np-\frac{4k}p.
   \]
   Hence the average of \(Q_A=\tfrac12 y^\top Ay\) over \(\mathrm{Max}_{+}\) equals \(\Phi(C)-2k/p\). Best switching forces this average to be nonnegative, and \(\max|Q_A|\ge\bigl|\mathbb E[Q_A]\bigr|\). \(\square\)

*Certified checks.* Frame identity \eqref{eq:max-frame} and LP value \(p\) at Paley \(p=3,5\); dual/primal algebra as above. At \(n=10\), the Max-Lipschitz with \(k_\star=5\) gives \(m_{10}\ge15-10/3=11.\overline{6}\) (true; actual \(13\)); edge lip only gave \(15-10=5\).

*Remark (E(1) status — sharpened criterion).* Max-Lipschitz yields
\[
m_n\ge\Phi(C)-\frac{2k_\star}p,
\qquad
\alpha_n\ge\tfrac12\sqrt{1-1/n}-\frac{2k_\star}{p\,n^{3/2}}.
\]
Hence **\(k_\star=O(n^{3/2})\) already forces \(\alpha_n\to\tfrac12\)** along the \(\rho=1\) family (gap \(O(n)=o(n^{3/2})\)), and more generally \(k_\star=o(n^2)\) is sufficient for a vanishing relative gap. The latter is the sharp threshold supplied by this inequality, not a necessary condition for E(1). Dual-Gaussian on \(W=A\circ C\) gives only \(k_\star\le\binom n2/2-\Omega(n^{3/2})=\Theta(n^2)\) for an arbitrary \(A\); the remaining stability task is to obtain any \(o(n^2)\) bound for a closest **global** \(\Phi\)-minimiser, or to bypass Hamming stability. Integral Max-covers of size \(p\) exist (LP support) but are stars and spike to \(\Phi>\Phi(C)\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.28 (size-\(p\) Max-covers: tight cover and spike; 2026-07-27).** Let \(C\) be a \(\rho=1\) conference of order \(n=p^2+1\) with the Max\(\pm\) frame identities of Prop~15.27 (so \(\mathbb E_{\pm}[C_{ij}y_iy_j]=\pm1/p\), and \(\mathbb E_{\mathrm{Max}_{-}}[yy^\top]=I-C/p\)). Write \(S_F(x)=\sum_{\{i,j\}\in F}C_{ij}x_ix_j\).

1. **Tight cover.** If \(F\) is a Max-cover (\(S_F\ge1\) on \(\mathrm{Max}_{+}\)) with \(|F|=p\), then \(S_F\equiv1\) on \(\mathrm{Max}_{+}\).
   *Proof.* \(\mathbb E_{+}[S_F]=|F|/p=1\) and \(S_F\ge1\) force equality. \(\square\)

2. **Spike from Max\(_{-}\).** If some \(z\in\mathrm{Max}_{-}\) has \(S_F(z)\ge1\), then \(\Phi(C\oplus F)\ge\Phi(C)+2\), since \(Q_{C\oplus F}(z)=-\Phi-2S_F(z)\). If \(S_F(z)=p\) then \(\Phi(C\oplus F)\ge\Phi+2p\).

3. **Affine-line \(p\)-stars at \(\infty\) (proved for Paley).** Identify vertices with \(\{\infty\}\cup\mathbb F_{p^2}\cong\{\infty\}\cup\mathrm{AG}(2,p)\), so \(C_{\infty u}=1\). Let \(L=\{a+td:t\in\mathbb F_p\}\) be an affine line of direction \(d\in\mathbb F_{p^2}^\times\), and let \(F\) be the \(p\)-star at \(\infty\) with leaves \(L\). Write \(\chi\) for the quadratic character of \(\mathbb F_{p^2}\). On \(\mathbb F_p^\times\subset\mathbb F_{p^2}\) one has \(\chi\equiv1\), hence
   \[
   \sum_{u\neq v\in L}C_{uv}
   =
   \sum_{t\neq s}\chi((t-s)d)
   =
   p(p-1)\,\chi(d).
   \]
   For any star the cross terms of \(S_F\) are adjacent edges, so the frame gives the exact second moment (using \(G_{\pm}=I\pm C/p\)):
   \[
   \mathbb E_{\pm}[S_F^2]
   =
   p\pm\frac1p\sum_{u\neq v\in L}C_{uv}
   =
   p\pm(p-1)\chi(d).
   \]
   - If \(\chi(d)=-1\) (nonsquare direction): \(\mathbb E_{+}[S_F^2]=1=\bigl(\mathbb E_{+}[S_F]\bigr)^2\), so \(S_F\equiv1\) on \(\mathrm{Max}_{+}\) (tight Max-cover); \(\mathbb E_{-}[S_F^2]=2p-1\neq1\), so \(S_F\not\equiv-1\) on \(\mathrm{Max}_{-}\). Since \(S_F\) is an odd integer with mean \(-1\), necessarily \(\max_{\mathrm{Max}_{-}}S_F\ge1\), and Lemma~2 yields \(\Phi(C\oplus F)\ge\Phi+2\).
   - If \(\chi(d)=+1\): symmetrically \(S_F\equiv-1\) on \(\mathrm{Max}_{+}\) (not a cover).

   Thus **exactly the nonsquare-direction affine lines** yield covering \(p\)-stars at \(\infty\), and every such star fails to undercut. (Half of the \(p(p+1)\) lines of \(\mathrm{AG}(2,p)\): \(p(p+1)/2\) covers.)

4. **Stronger spike at covering stars (certified \(p=3,5\)).** Every covering \(p\)-star (all centres, not only \(\infty\)) has \(\max_{\mathrm{Max}_{-}}S_F=p\), hence \(\Phi(C\oplus F)=\Phi+2p\). Counts: 60 at \(n=10\) (all \(\Phi=21\)); 390 at \(n=26\) (MITM sample all \(\Phi=75\)). Spike witnesses are constant on leaves\(\cup\{\mathrm{centre}\}\).

5. **All size-\(p\) Max-covers at \(n=10\).** Exhaustive: 405 size-\(p\) Max-covers (60 stars + 345 non-stars); every one has \(\max_{\mathrm{Max}_{-}}S_F\in\{1,3\}\) and \(\Phi\in\{17,21\}\); **zero undercuts** of \(\Phi=15\). Moreover \(\mathbb E_{-}[S_F^2]\ge7/3>1\) on this census, so none can have \(S_F\equiv-1\) on \(\mathrm{Max}_{-}\). Card-min Max-covers never undercut at \(n=10\); undercutters begin at \(k=5\) (N10-S matchings).

*Evidence:* `evidence/E1_SIZE_P_MAXCOVER.md`, `e1_size_p_maxcover.json`. Max\(_{-}\) frame certified at \(p=3,5\).

*Remark (E(1) status).* Prop 15.28 shows that LP-tight Max-covers (size \(p\)) cannot undercut via the Max\(\pm\) analysis on the proved/certified range: they are forced to spike on \(\mathrm{Max}_{-}\). The \(n=10\) undercut uses a *strictly larger* cover (\(k=5>p\)). Closing E(1) still needs either (i) \(k_\star=O(n^{3/2})\) for minimisers, (ii) a general no-undercut / controlled-gap theorem for all Max-covers of size \(o(n^2)\), or (iii) exact Paley optimality for all \(p\ge5\) with the known \(O(1)\) gap at \(p=3\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.29 (odd matching parity and matching-cover spikes; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\) (\(p\) odd prime), \(\Phi=\tfrac12 np\).

1. **Parity.** The perfect-matching size \(n/2\) is odd. Hence for every perfect matching \(M\) and every \(x\in\{\pm1\}^n\), \(S_M(x)\) is an odd integer. In particular, either \(\min_{\mathrm{Max}_{+}}S_M\ge1\) (Max-cover) or \(\min_{\mathrm{Max}_{+}}S_M\le-1\). In the latter case \(\Phi(C\oplus M)\ge\Phi+2\).

2. **Undercutters are two-sided.** If \(\Phi(C\oplus F)<\Phi\), then \(S_F\ge1\) on \(\mathrm{Max}_{+}\) and \(S_F\le-1\) on \(\mathrm{Max}_{-}\) (integer arithmetic on \(Q_C=\pm\Phi\)). Evidence note: `evidence/E1_TWOSIDED_COVER.md`.

3. **Certified \(n=26\) matching covers.** Perfect-matching Max-covers exist at \(p=5\) (SA finds them; earlier “no cover” reports were incomplete). Among 48 seeds, 3 covers were found; all three are two-sided, achieve \(\max_{\mathrm{Max}_{\pm}}|Q|=63=\Phi-2\), but exact MITM \(\Phi(C\oplus M)=65=\Phi(C)\) (non-eigenvector spike of \(+2\)). **Zero undercuts.** Contrast \(n=10\): two-sided matching covers achieve \(m_{10}=13\). Evidence: `evidence/E1_MATCHING_COVER_SPIKE.md`, `e1_n26_matching_covers_census.json`.

*Remark.* Matching dichotomy for E(1) must account for covers that fail only by spike control. Prop 15.26 (local maximality) holds at both \(p=3,5\); global Max-determination of \(\Phi(C\oplus M)\) is special to \(n=10\) among tested orders. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.30 (matching spike criterion; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), and \(M\) a perfect matching of \(K_n\). Write \(A=C\oplus M\) and \(S_M(x)=\sum_{\{i,j\}\in M}C_{ij}x_ix_j\).

1. **Criterion (proved).** If there exists \(x\in\{\pm1\}^n\) with \(S_M(x)=-p\) and \(Q_C(x)\ge\Phi-2p\) (or \(S_M(x)=p\) and \(Q_C(x)\le-(\Phi-2p)\)), then \(\Phi(A)\ge\Phi(C)\).
   *Proof.* \(Q_A=Q_C-2S_M\), so the first alternative yields \(Q_A(x)\ge\Phi-2p+2p=\Phi\). \(\square\)

2. **Single-bit drop.** For \(y\in\mathrm{Max}_{+}\) and any coordinate \(i\), \(Q_C(y^{\oplus i})=\Phi-2p\). Hence non-max boolean vectors attain at least the level \(\Phi-2p\). At \(p=3\) this is exact as the global second max of \(|Q_C|\) on non-max vectors (full half-cube). At \(p=5\) sampling finds second max \(55=\Phi-2p\).

3. **Multi-bit formula.** \(Q_C(y^{\oplus F})=\Phi-2p|F|+4\sum_{i<j\in F}C_{ij}y_iy_j\). On a \(y\)-switched clique of size \(r\), this equals \(\Phi-2r(p-r+1)\); at \(r=p\) one recovers \(\Phi-2p\). Square-direction affine lines of \(\mathrm{AG}(2,p)\) are \(p\)-cliques of Paley (Prop 15.28 geometry).

4. **Certified dichotomy at \(n=10\).** Over all 945 perfect matchings, the criterion holds on exactly the \(801\) non-undercutting matchings and fails on exactly the \(144\) undercutters (where \(\max Q_C\) on \(\{S_M=-3\}\) is \(7<9=\Phi-2p\)).

5. **Certified at \(n=26\).** Every tested perfect matching (15 random + 3 SA Max-covers) satisfies the criterion; all three Max-covers have exact MITM \(\Phi=\Phi(C)\). On the level \(S_M=-p\) one has the identity \(\mathbb E[Q_C]=-p\).

*Open for matching non-undercut when \(p\ge5\).* Prove the criterion for every perfect matching (tail \(\max_{S_M=-p}Q_C\ge\Phi-2p\)). That would give \(\Phi(C\oplus M)\ge\Phi(C)\) for all \(M\), hence E(1) under the matching dichotomy. Full E(1) still requires control of non-matching undercutters. Evidence: `evidence/E1_MATCHING_SPIKE_CRITERION.md`. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.31 (clique-flip sufficiency for matching Max-covers; 2026-07-27).** Let \(M\) be a perfect matching Max-cover of Paley order \(n=p^2+1\). If there exist \(y\in\mathrm{Max}_{+}\) and a transversal \(p\)-set \(F\) such that \(F\) is a clique in \(W_{ij}=C_{ij}y_iy_j\), \(S_M(y)=s_0\) with \(s_0+3p\equiv0\pmod4\), and \(\sum_{i\in F}\chi_i(y)=(s_0+p)/2\), then
\[
Q_C(y^{\oplus F})=\Phi-2p,\qquad S_M(y^{\oplus F})=-p,\qquad Q_{C\oplus M}(y^{\oplus F})=\Phi,
\]
so \(\Phi(C\oplus M)\ge\Phi(C)\).

*Arithmetic.* Full-clique flips reach the \(\Phi-2p\) threshold only for \(|F|\in\{1,p\}\); \(|F|=1\) is incompatible with Max-covers. Thus only \(|F|=p\) applies, requiring \(s_0+3p\equiv0\pmod4\). At \(p=3\) undercutting matchings have \(S_M(\mathrm{Max}_{+})\subseteq\{1,5\}\) (no admissible \(s_0\)), blocking the construction. At \(p=5\), \(\mathbb E[S_M]=2.6<3\) forces every Max-cover to attain \(S_M=1\).

*Design constants at \(p=5\) (certified).* Seidel-consistent \(p\)-sets (\(C_{ab}C_{ac}C_{bc}=1\)): 390. Each has exactly 60 Max\(_{+}\) extensions. Every tested matching has \(\ge236\) transversal consistent \(p\)-sets. All SA Max-covers admit a clique-flip (MITM \(\Phi=\Phi(C)\)). Evidence: `evidence/E1_CLIQUE_FLIP.md`, `e1_clique_flip_covers.json`.

*Open.* Existence of \((y,F)\) for every Max-covering matching when \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.32 (Γ-pairing reformulation and mod-4 dichotomy; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), \(m=n/2\), and \(M=\{e_a=(i_a,j_a)\}_{a=1}^{m}\) a perfect matching. Write \(S_M=\sum_a C_{e_a}x_{i_a}x_{j_a}\) and \(R_M=Q_C-S_M\).

1. **Γ-pairing (proved).** With free signs \(u_a=x_{i_a}\) and relative signs \(z_a=x_{i_a}x_{j_a}\),
   \[
   R_M=\sum_{a<b}\gamma_{ab}(z_a,z_b)\,u_a u_b,\qquad
   \gamma_{ab}\in\{-4,-2,0,2,4\}.
   \]
   The spike criterion \(\max_{S_M=-p}Q_C\ge\Phi-2p\) is equivalent to
   \[
   \max_{\langle c,z\rangle=-p}\max_{u\in\{\pm1\}^m}\tfrac12 u^\top\Gamma(z)u
   \;\ge\;p(m-1).
   \]
   Shipped: `maxR_matching_level` in `src/minmax_quadratic.py`.

2. **Coordinate product on Max\(_{+}\) (proved).** Every \(y\in\mathrm{Max}_{+}\) has
   \(\pi(y):=\prod_{v=1}^n y_v=(-1)^{p(p-1)/2}\), constant on \(\mathrm{Max}_{+}\).
   (Halfspace construction gives the value; \(\mathrm{P}\Gamma\mathrm{L}\) acts by coordinate permutation and preserves the product; the boolean \(+p\)-orbit is the Aut-orbit of the halfspace.)
   Thus \(\pi=+1\) when \(p\equiv1\pmod4\) and \(\pi=-1\) when \(p\equiv3\pmod4\).

3. **Mod-4 constancy (proved).** For any perfect matching \(M\), \(S_M(y)\bmod 4\) is **constant** for \(y\in\mathrm{Max}_{+}\). Indeed with \(\chi_a=C_{e_a}y_{i_a}y_{j_a}\) one has \(\prod_a\chi_a=(\prod_e C_e)\,\pi\), and \(S_M=\sum\chi_a\equiv m-2k\pmod4\) where \(k=\#\{\chi_a=-1\}\), so the residue is determined by \(\prod\chi_a\) alone.

4. **Case split for the spike (structure).** Write \(r\) for the constant residue of \(S_M\) on \(\mathrm{Max}_{+}\).
   - If \(-p\equiv r\pmod4\) and \(-p\in S_M(\mathrm{Max}_{+})\), the criterion holds with \(Q=\Phi\) (maximiser on the level).
   - If \(-p+2\equiv r\pmod4\) and \(-p+2\in S_M(\mathrm{Max}_{+})\), a **1-bit flip** of a maximiser lands on \(Q=\Phi-2p\) and \(S_M=-p\) (formula \(S(y^{\oplus i})=S(y)-2\chi_i\)), so the criterion holds.
   - If \(M\) is a Max-cover (\(S_M\ge1\) on \(\mathrm{Max}_{+}\)), then \(r\equiv1\pmod4\) at \(p=5\) and the 1-bit route is blocked; the clique-flip of Prop 15.31 is the remaining constructive path.

5. **Certified census at \(p=5\) (not a proof).** Over 80 random perfect matchings + SA min-\(\max R\) + 2 Max-covers: \(\max R\in\{60,70\}\) always, \(\min=60=p(m-1)\) (tight on covers), 0 counterexamples to the criterion; MITM \(\Phi(C\oplus M)=\Phi(C)\) on covers. Over 500 random matchings: every mod-\(3\) matching attains \(S_M=-5\); every mod-\(1\) matching attains \(S_M=-3\). Evidence: `evidence/E1_GAMMA_PAIRING.md`, `e1_gamma_forall_census.json`.

*Open.* Prove the attainment lemmas (\(-p\) or \(-p+2\) in the image of \(S_M|_{\mathrm{Max}_{+}}\) whenever the residue allows) for all \(p\ge5\), and clique-flip existence on every Max-cover; then matching non-undercut for \(p\ge5\). Matching dichotomy / non-matching undercutters / \(k_\star\) still separate. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.33 (non-covers cannot undercut; criterion is not necessary; 2026-07-27).** Let \(C\) be a \(\rho=1\) conference of order \(n\), \(\Phi=\Phi(C)\), and \(F\) any flip set. Write \(A=C\oplus F\) and \(S_F(y)=\sum_{e\in F}C_e y_iy_j\).

1. **Non-cover \(\Rightarrow\) strict raise (proved).** If there exists \(y\in\mathrm{Max}_{+}\) with \(S_F(y)\le0\), then
   \[
   Q_A(y)=\Phi-2S_F(y)\ge\Phi,
   \]
   and if \(S_F(y)\le-1\) (always, when scores are odd — e.g. perfect matchings with \(n/2\) odd) then \(Q_A(y)\ge\Phi+2\). Hence \(\Phi(A)\ge\Phi+2>\Phi(C)\): **non-covers cannot undercut.**
   *Proof.* Immediate from \(Q_A=Q_C-2S_F\) on \(\mathrm{Max}_{+}\). \(\square\)
   Equivalently: every strict undercutter is a Max\(_{+}\) cover (cf. Lemma U1 / Prop 15.29).

2. **Perfect matchings (proved reduction).** For a perfect matching \(M\), either \(M\) is a Max-cover and may or may not undercut, or \(M\) is a non-cover and \(\Phi(C\oplus M)\ge\Phi+2\). Thus **matching non-undercut reduces entirely to Max-covering matchings.**

3. **Spike criterion is not necessary (certified counterexample).** At \(p=5\), the perfect matching
   \[
   M_0=\{\{5,3\},\{2,14\},\{22,18\},\{23,7\},\{12,16\},\{6,8\},\{1,19\},\{17,25\},\{9,0\},\{21,10\},\{11,24\},\{13,20\},\{4,15\}\}
   \]
   has \(\max R_{M_0}=54<60=\Phi-p\) (criterion fails) but \(\min_{\mathrm{Max}_{+}}S_{M_0}=-1\) (non-cover) and exact MITM \(\Phi(C\oplus M_0)=75>\Phi\). Evidence: `evidence/e1_criterion_fail_no_undercut.json`. The criterion remains a useful *sufficient* test for covers; it is not a characterisation of non-undercut.

4. **Max-cover census at \(p=5\) (not a proof).** All SA-found Max-covering matchings (8 total: 2 in `e1_gamma_forall_census.json` + 6 in `e1_gamma_cover_batch.json`) satisfy the spike criterion with equality \(\max R=60\), admit clique-flips, and have MITM \(\Phi=\Phi(C)\).

*Open for matching non-undercut when \(p\ge5\).* Prove every Max-covering perfect matching has \(\Phi(C\oplus M)\ge\Phi(C)\) (e.g. via clique-flip Prop 15.31 for all covers). Then, with \(p=3\) gap \(O(1)\), matching dichotomy would give E(1) along matchings. Full E(1) still needs non-matching undercutters. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.34 (matching flip algebra; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\) and \(M\) a perfect matching. Write \(D\) for the matrix with \(D_{ij}=C_{ij}\) on edges of \(M\) and \(0\) elsewhere, and set \(A=C-2D=C\oplus M\).

1. **Involution identity (proved).** \(D^2=I\). Consequently
   \[
   A^2=(n+3)I-2(CD+DC).
   \]
   *Proof.* Each vertex lies in exactly one matching edge, so \((D^2)_{ii}=C_{i\pi(i)}^2=1\) and \((D^2)_{ij}=0\) for \(i\neq j\). Expand \(A^2=(C-2D)^2=C^2-2(CD+DC)+4D^2=(n-1)I-2(CD+DC)+4I\). \(\square\)

2. **Certified spectrum at \(p=5\).** For every SA Max-covering matching tested, \(\|A\|_{\mathrm{op}}=\sqrt{41}\) exactly (to numerical precision \(10^{-12}\)); random matchings have \(\|A\|_{\mathrm{op}}\in[6.79,7.00]\). Two Aut-invariants (K\(_{2,2}\) type counts) show \(\ge8\) distinct Max-cover classes, all sharing this op-norm. When \(\Phi(A)=\Phi(C)=65\), one has \(\rho(A)=p/\|A\|_{\mathrm{op}}=5/\sqrt{41}\).

*Open.* Prove \(\|C\oplus M\|_{\mathrm{op}}=\sqrt{p^2+16}\) (or a bound forcing \(\Phi(A)\ge\Phi(C)\)) for every Max-covering matching when \(p\ge5\); or complete clique-flip existence. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.35 (Max-cover matching structure at \(p=5\); 2026-07-27).** Let \(C\) be Paley of order \(n=26\) (\(p=5\)), \(\Phi=65\), and \(M\) a perfect matching Max-cover (\(\min_{\mathrm{Max}_{+}}S_M\ge1\)).

1. **Forced \(S_M=1\) and residue (proved).** \(\mathbb E[S_M]=13/5=2.6<3\) and \(S_M\) odd on \(\mathrm{Max}_{+}\), so \(S_M\) attains \(1\). The residue of \(S_M\) on \(\mathrm{Max}_{+}\) cannot be \(3\bmod4\) (that would force \(\min S_M\ge3>2.6\)). Hence \(S_M\equiv1\pmod4\), \(s_0=1\) is admissible for clique-flip (Prop 15.31), and \(\Sigma_{\mathrm{need}}=3\).

2. **Certified census of 11 Max-covers (not a forall proof).** Independent SA campaigns produced **11** distinct Max-covering perfect matchings. Every one is:
   - two-sided (\(\min_{\mathrm{Max}_{+}}S=1\), \(\max_{\mathrm{Max}_{-}}S=-1\));
   - inclusion-minimal and inclusion-maximal as Max-covers;
   - \(\|C\oplus M\|_{\mathrm{op}}=\sqrt{41}\) exactly;
   - spike-criterion tight (\(\max R=60\));
   - clique-flip capable;
   - exact MITM \(\Phi(C\oplus M)=\Phi(C)\).
   Only two \(S_{+}\) distributions appear: \((1^{156},5^{104})\) and \((1^{168},5^{80},9^{12})\). Evidence: `evidence/e1_maxcover_full_census.json`.

3. **Uniform sampling.** Among \(2\cdot10^4\) random perfect matchings, **0** Max-covers (Max-covers are extreme tail events).

*Open.* Prove clique-flip / \(\Phi\ge\Phi(C)\) for every Max-covering perfect matching when \(p\ge5\); then non-matching undercutters / \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.36 (matching flip block algebra and Max-cover spectrum; 2026-07-27).** Let \(C\) be Paley of order \(n=p^2+1\), \(M\) a perfect matching, \(D\) the signed matching matrix, \(A=C-2D\), and \(V_\pm=\ker(C\mp pI)\).

1. **Block formulae (proved for every matching).** \(B:=CD+DC\) always commutes with \(C\) and \(D\). On \(V_+\oplus V_-\),
   \[
   B\big|_{V_+}=2p\,D_{++},\qquad B\big|_{V_-}=-2p\,D_{--},
   \]
   and \(\mathrm{tr}(D_{++})=n/(2p)\). For \(y\in\mathrm{Max}_+\subset V_+\) one has \(S_M(y)=\tfrac12 y^\top D_{++}y\). Also
   \(\|A\|_{\mathrm{op}}^2=(n+3)-2\lambda_{\min}(B)\).

2. **Certified at \(p=5\) for Max-covers (not forall).** Every SA Max-covering matching has \(\lambda_{\min}(B)=-6\) and \(\|A\|_{\mathrm{op}}=\sqrt{41}=\sqrt{p^2+16}\). At least two \(D_{++}\) spectral types occur (simple \(\{-3/5,0^{(8)},(4/5)^{(4)}\}\) and a mixed type with the same \(\lambda_{\min}(D_{++})=-3/5\)); both have tight spike \(\max R=60\), clique-flip, and MITM \(\Phi=\Phi(C)\). Random non-covers have \(\lambda_{\min}(B)\in[-10,-8]\) and larger op-norm. At \(p=3\), the 144 undercutting Max-covers share a single smaller op-norm \(\approx3.933\) (golden-ratio \(B\)-spectrum). Evidence: `evidence/E1_MAXCOVER_SPECTRUM.md`, `e1_maxcover_spectrum.json`.

*Open.* Prove \(\lambda_{\min}(B)=-6\) (or \(\|A\|_{\mathrm{op}}=\sqrt{p^2+16}\)) for every Max-cover when \(p\ge5\), and upgrade to \(\Phi(A)\ge\Phi(C)\); or complete clique-flip existence. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.37 (continuous Γ-bound pattern on Max-covers; 2026-07-27).** Let \(C\) be Paley of order \(n=26\) and \(M\) a Max-covering perfect matching. On the Γ-pairing level \(S_M=-p\), every tested \(M\) satisfies
\[
\min_z\lambda_{\max}(\Gamma(z))\;\ge\;9.38758\;>\;\frac{2p(m-1)}{m},
\]
so the continuous bound \(\tfrac m2\lambda_{\max}\ge p(m-1)\) holds for **all** \(z\) on the level (two numerical classes). Discrete \(\max R=60\) and clique-flip hold on all 11 stored covers. GW theory does not close the discrete gap (SDP\(\approx63.8\), \(\alpha\cdot\mathrm{SDP}<60\)). Residue-\(1\) random matchings also admitted clique-flips in sampling. Evidence: `evidence/E1_MAXCOVER_CONTINUOUS_BOUND.md`.

*Open.* Prove the continuous bound and/or clique-flip for every Max-cover when \(p\ge5\); then matching dichotomy / \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.38 (two-sided size-5 Max-covers on \(C_{10}\); 2026-07-28).** Let \(C\) be Paley of order \(n=10\), \(\Phi=15\). Among all edge sets \(F\) with \(|F|=5\) that are two-sided Max-covers (\(\min_{\mathrm{Max}_{+}}S_F\ge1\) and \(\max_{\mathrm{Max}_{-}}S_F\le-1\)):

1. **Undercutters are exactly the matchings.** Exactly **144** such \(F\) undercut, each is a perfect matching (\(\Delta=1\)), and each has \(\Phi(C\oplus F)=13=m_{10}\).
2. **Higher \(\Delta\) never undercuts at this cardinality.** Counts: \(\Delta=2\): 8730 covers with \(\Phi\in\{15,17,19\}\); \(\Delta=3\): 7920 with \(\Phi\in\{17,19\}\); \(\Delta=4\): 360 with \(\Phi=19\). All have \(\Phi\ge15\).
3. **Total.** 17154 two-sided \(k=5\) Max-covers; only the 144 matchings undercut.

*Proof.* Exhaustive scan of \(\binom{45}{5}=1{,}221{,}759\) five-edge sets; Max\(\pm\) boolean \(\pm3\)-eigenvectors; exact \(\Phi\) by cube enumeration via `form_Q`. Parallel re-run: `src/n10_twosided_k5_classify.py` (80 workers). Evidence: `evidence/E1_N10_TWOSIDED_K5.md`, `e1_n10_twosided_k5_classify.json`. \(\square\)

*Remark.* Strengthens N10-S: at matching cardinality, two-sided Max-covers with \(\Delta\ge2\) cannot undercut. Supports a low-\(\Delta\) undercutter pattern for E(1), but does **not** prove \(k_\star=O(n^{3/2})\) or matching non-undercut for \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.39 (clique-flip count on Max-cover matchings, \(p=5\); 2026-07-28).** On every stored Max-covering perfect matching of Paley \(C_{26}\), the number \(N_{\mathrm{flip}}\) of pairs \((y,F)\) as in Prop 15.31 (with \(s_0=1\), \(\Sigma=3\)) satisfies \(N_{\mathrm{flip}}\ge24>0\), so clique-flip applies and \(\Phi(C\oplus M)=\Phi(C)\). Observed counts include \(24,48,120\). Evidence: `evidence/E1_CLIQUE_FLIP_COUNT.md`.

*Open.* Prove \(N_{\mathrm{flip}}\ge1\) for every Max-cover matching when \(p=5\) (then matching non-undercut at \(p=5\)); lift to \(p\ge7\) and \(k_\star\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.40 (edge-minimal undercutters have gap at most 2; 2026-07-28).** Let \(C\) be any Seidel matrix of order \(n\) with \(\Phi(C)\in\mathbb Z\), and let \(F\) be a nonempty edge set such that \(\Phi(C\oplus F)<\Phi(C)\) while \(\Phi(C\oplus(F\setminus\{e\}))\ge\Phi(C)\) for every \(e\in F\) (edge-minimal undercutter). Then
\[
\Phi(C\oplus F)\;\ge\;\Phi(C)-2.
\]
*Proof.* For any \(e\in F\), Prop 15.20b with \(k=1\) gives \(|\Phi(C\oplus F)-\Phi(C\oplus(F\setminus\{e\}))|\le2\). Combined with \(\Phi(C\oplus F)<\Phi(C)\le\Phi(C\oplus(F\setminus\{e\}))\),
\[
\Phi(C\oplus F)
\;\ge\;
\Phi(C\oplus(F\setminus\{e\}))-2
\;\ge\;
\Phi(C)-2.
\]
(The same holds with \(\Phi\) replaced by any real threshold \(t\) for which \(F\) is edge-minimal among sets with \(\Phi(C\oplus\cdot)<t\).) \(\square\)

*Certified.* On Paley \(C_{10}\), all 144 matching undercutters and sampled undercutting 6-cycles are edge-minimal and have gap exactly 2. Evidence: session checks; N10-S/C6.

*Remark (does **not** yet give \(m_n\ge\Phi-2\)).* The lemma bounds edge-minimal undercutters only. A global lower bound \(m_n\ge\Phi(C)-2\) would require showing no deeper undercut exists at larger Hamming distance (or that a closest undercutter realises \(m_n\)). That step is **open**. If proved, E(1) follows (gap \(O(1)=o(n^{3/2})\)) and \(L=\tfrac12\) by denseness on the \(\rho=1\) family. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.41 (first-hit + no-descent framework toward \(m_n\ge\Phi-2\); 2026-07-28).** Let \(C\) be a Seidel matrix of order \(n\) with \(\Phi(C)\in\mathbb Z\). Write \(A_F:=C\oplus F\).

1. **First-hit lemma (proved).** Along any edge-adding chain \(\emptyset=F_0\subset F_1\subset\cdots\subset F_k\) with \(|F_i|=i\), if \(i^\star\) is minimal such that \(\Phi(A_{F_{i^\star}})<\Phi(C)\), then
   \[
   \Phi(A_{F_{i^\star}})\;\ge\;\Phi(C)-2.
   \]
   *Proof.* \(F_{i^\star}=F_{i^\star-1}\cup\{e\}\) and \(\Phi(A_{F_{i^\star-1}})\ge\Phi(C)\), so Prop 15.20b with one edge yields the claim. \(\square\)

2. **Dangerous-edge criterion (proved).** Suppose \(\Phi(A_F)=\Phi(C)-2\). For \(e=(u,v)\notin F\) set \(B=A_F\oplus e\) and \(\sigma_x:=(A_F)_{uv}\,x_u x_v\in\{\pm1\}\). Then \(Q_B(x)=Q_{A_F}(x)-2\sigma_x\). Consequently:
   - if some maximiser \(x\) with \(Q_{A_F}(x)=\Phi(C)-2\) has \(\sigma_x=-1\), then \(Q_B(x)=\Phi(C)\) and \(\Phi(B)\ge\Phi(C)\);
   - if some maximiser with \(Q_{A_F}(x)=-(\Phi(C)-2)\) has \(\sigma_x=+1\), same.
   Hence \(\Phi(B)\le\Phi(C)-4\) is possible only if \(\sigma\equiv+1\) on all \(+\) maximisers and \(\sigma\equiv-1\) on all \(-\) maximisers (**dangerous edge**). \(\square\)

3. **No-descent lemma (OPEN in general).** If \(\Phi(A_F)=\Phi(C)-2\), then for every \(e\notin F\), \(\Phi(A_F\oplus e)\ge\Phi(C)-2\).

4. **Conditional settlement.** If the no-descent lemma holds for all flip sets on the \(\rho=1\) Paley family \(n=p^2+1\), then by induction on Hamming distance every Seidel matrix \(A\) of those orders satisfies \(\Phi(A)\ge\Phi(C)-2\). Thus \(m_n\ge\Phi(C)-2\), the gap is \(O(1)=o(n^{3/2})\), E(1) holds, and \(L=\tfrac12\) by Proposition 6.1 along the consecutive-prime-square orders. **F13:** this must not be claimed from Prop 15.40 alone; no-descent is an independent lemma about \(\Phi\), not abstract 2-Lipschitz calculus.

5. **Certified at \(n=10\) (matching undercutters; not a general proof).** For all **144** perfect-matching undercutters of Paley \(C_{10}\) (\(\Phi=13=\Phi-2\)):
   - number of dangerous external edges: **0**;
   - every single-edge extension has \(\Phi\ge15=\Phi(C)\) (min observed \(15\));
   - multi-edge random extensions (depth \(2\)–\(12\)) stay at \(\Phi\ge13\) (0 deepenings).
   Parallel cert: `src/e1_n10_nodescent.py` (\(W=86\)), `evidence/e1_n10_nodescent.json`, `evidence/E1_NODESCENT.md`.

*Open.* Prove no-descent (or: no dangerous edges, or non-maximiser spike under alignment) for all undercutters on the \(\rho=1\) family; alternatively complete Max-cover clique-flip for \(p\ge5\). Matching dichotomy / \(k_\star\) remain separate. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.42 (Max\(\pm\) dichotomy, counting freeness, tight no-descent; 2026-07-28).** Let \(C\) be a \(\rho=1\) conference of order \(n\) with the Max\(_{+}\) frame identity of Prop 15.27, \(\Phi=\Phi(C)\), \(p=\sqrt{n-1}\), and \(A=C\oplus F\). Write \(S_F(y)=\sum_{e\in F}C_e y_iy_j\), \(s_{+}=\min_{\mathrm{Max}_{+}}S_F\), \(s_{-}=\max_{\mathrm{Max}_{-}}S_F\), and \(f_e(y)=C_e y_iy_j\).

1. **Max\(\pm\) evaluation dichotomy (proved).**
   \[
   \Phi(A)\;\ge\;\max\Bigl(
   \max_{y\in\mathrm{Max}_{+}}\lvert\Phi-2S_F(y)\rvert,\;
   \max_{y\in\mathrm{Max}_{-}}\lvert-\Phi-2S_F(y)\rvert
   \Bigr).
   \]
   Consequently:
   - if \(s_{+}\le -1\), then \(\Phi(A)\ge\Phi+2\) (Prop 15.33);
   - if \(s_{+}=0\), then \(\Phi(A)\ge\Phi\);
   - if \(s_{+}=1\), then \(\Phi(A)\ge\Phi-2\);
   - if \(s_{-}\ge 0\), then \(\Phi(A)\ge\Phi\);
   - if \(s_{-}=-1\), then \(\Phi(A)\ge\Phi-2\).

   Therefore \(\Phi(A)<\Phi-2\) is possible only if \(A\) is **deep two-sided**: \(s_{+}\ge 2\) and \(s_{-}\le -2\).
   *In particular every matrix with \(s_{+}\le 1\) or \(s_{-}\ge -1\) already satisfies \(\Phi(A)\ge\Phi-2\).* \(\square\)

2. **Counting freeness (proved).** For every edge \(e\), \(\sum_{y\in\mathrm{Max}_{+}}f_e(y)=|\mathrm{Max}_{+}|/p\) (frame). Hence if \(f_e\equiv +1\) on a subset \(U\subseteq\mathrm{Max}_{+}\), then
   \[
   |U|\;\le\;|\mathrm{Max}_{+}|\,\frac{p+1}{2p},
   \]
   with equality only if \(f_e\equiv -1\) on \(\mathrm{Max}_{+}\setminus U\). \(\square\)

3. **Tight \(S\equiv 1\) no-descent (proved).** If \(S_F\equiv 1\) on \(\mathrm{Max}_{+}\) (equivalently \(k=|F|=p\) and \(s_{+}=1\)), then \(U=\mathrm{Max}_{+}\) violates the counting bound for \(f_e\equiv 1\) (since \((p+1)/(2p)<1\)). So every \(e\notin F\) has some \(y\in\mathrm{Max}_{+}\) with \(f_e(y)=-1\), whence \(S_{F\cup\{e\}}(y)=0\) and \(Q_{A\oplus e}(y)=\Phi\). Thus \(\Phi(A\oplus e)\ge\Phi\). \(\square\)

4. **Tight \(S\equiv 2\) no-descent (proved).** If \(S_F\equiv 2\) on \(\mathrm{Max}_{+}\) (which holds automatically for any Max\(_{+}\) cover with \(|F|=2p\), since \(\mathbb E[S]=2\) and \(S\ge 2\) force \(S\equiv 2\)), then similarly \(f_e\not\equiv 1\) on \(\mathrm{Max}_{+}\). Any \(y\) with \(f_e(y)=-1\) has \(S_{F\cup\{e\}}(y)=1\) and \(Q_{A\oplus e}(y)=\Phi-2\). Thus \(\Phi(A\oplus e)\ge\Phi-2\). \(\square\)

5. **Type I freeness when \(N_1\) is large (proved).** Write \(N=|\mathrm{Max}_{+}|\) and \(N_1=\#\{y\in\mathrm{Max}_{+}:S_F(y)=1\}\). If \(s_{+}=1\) and \(N_1>N(p+1)/(2p)\), then no edge freezes to \(+1\) on \(\mathrm{Max}_{+1}\). When moreover \(\Phi(A)=\Phi-2\), one has \(\mathrm{Max}_{+1}\subseteq\mathrm{Max}(A)\) (positive maximisers), so every \(e\notin F\) has a maximiser with \(\sigma_e=-1\), hence \(\Phi(A\oplus e)\ge\Phi\) (no-descent, strong form). \(\square\)

6. **Type I size bound for strict freeness (proved).** If \(s_{+}=1\), scores are odd, and \(|F|=k\le 2p-2\), then \(N_1/N\ge(3-k/p)/2>(p+1)/(2p)\), so part 5 applies. \(\square\)

7. **Reduction of \(m_n\ge\Phi-2\) (proved equivalence).** Assume \(\Phi\) is integer-valued with fixed parity under edge flips (true for all Seidel matrices of a fixed order: each flip changes every \(Q_x\) by \(\pm 2\)). Then:
   - by part 1, \(\Phi(A)\le\Phi-4\) forces deep two-sided;
   - by Prop 15.40, no edge-minimal undercutter has \(\Phi\le\Phi-4\);
   - a minimal-Hamming counterexample \(A\) with \(\Phi(A)=\Phi-4\) must therefore admit some \(e\) with \(\Phi(A\oplus e^{-1})=\Phi-2\) (single-edge descent from a gap-\(2\) undercutter).
   Hence **\(m_n\ge\Phi-2\) on the \(\rho=1\) family if and only if no-descent holds for every gap-\(2\) undercutter** (Prop 15.41(3)). Parts 3–6 prove no-descent for all tight \(S\equiv1\), tight \(S\equiv2\), and Type I covers with \(k\le 2p-2\) or \(N_1\) large. \(\square\)

8. **Certified at \(n=10\).** All 144 matching undercutters are Type I with \(N_1=10>8=N(p+1)/(2p)\) (strict freeness; 0 dangerous edges; add-1 gives \(\Phi\ge15\)). All 360 undercutting 6-cycles are tight \(S\equiv 2\) on \(\mathrm{Max}_{+}\) (part 4; add-1 gives \(\Phi\ge15\)). Combined with \(m_{10}=13\), no-descent and \(m_{10}=\Phi-2\) hold at \(n=10\). Evidence: `evidence/E1_NODESCENT.md`, `e1_n10_nodescent.json`; session C6 checks.

*Open for \(L=\tfrac12\).* Prove no-descent for the remaining gap-\(2\) undercutters (Type I with \(k\ge 2p-1\) and small \(N_1\); deep non-tight with \(k>2p\)), or prove no such undercutters exist on the \(\rho=1\) family. Then \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.43 (no-descent for Type I freeness class and tight deep; residual isolation; 2026-07-28).** Continue the notation of Prop 15.42. Let \(B=C\oplus G\) be a gap-\(2\) undercutter (\(\Phi(B)=\Phi-2\)) and \(e\notin G\).

1. **Type I with freeness \(\Rightarrow\) strong no-descent (proved).** If \(s_+(G)=1\) and \(N_1>N(p+1)/(2p)\), then some \(y\in\mathrm{Max}_{+1}\subseteq\mathrm{Max}(B)\) has \(f_e(y)=-1\), so \(Q_{B\oplus e}(y)=\Phi\) and \(\Phi(B\oplus e)\ge\Phi\). \(\square\)

2. **Tight deep \(\Rightarrow\) weak no-descent (proved).** If \(S_G\equiv2\) on \(\mathrm{Max}_{+}\), then some \(y\in\mathrm{Max}_{+}\) has \(f_e(y)=-1\), so \(S_{G\cup\{e\}}(y)=1\), \(Q_{B\oplus e}(y)=\Phi-2\), and \(\Phi(B\oplus e)\ge\Phi-2\). \(\square\)

3. **Equality-case isolation for Type I freeness failure (proved structure).** Suppose \(s_+=1\) and \(f_e\equiv+1\) on \(\mathrm{Max}_{+1}\) (freeness fails). Counting forces \(N_1\le N(p+1)/(2p)\). Combining the Max\(_{+}\) bounds \(Q\le\Phi-2\) on \(\mathrm{Max}_{+}\) and \(Q\le\Phi-6\) on \(\{S\ge3\}\) with \(\mathbb E[Q]=\Phi-2k/p\) forces, at counting equality \(N_1=N(p+1)/(2p)\), that \(S\in\{1,3\}\) on \(\mathrm{Max}_{+}\) and \(k=2p-1\), with \(f_e=2-S\) on \(\mathrm{Max}_{+}\). Then \(H:=G\cup\{e\}\) is a tight \(S\equiv2\) cover of size \(2p\). (A parallel boundary with \(S\in\{1,5\}\) and \(k=3p-2\) forces the affine relation \(S+2f_e=3\).) In the \(k=2p-1\) equality case, no-descent for this \(e\) reduces to \(\Phi(C\oplus H)\ge\Phi-2\) for the tight cover \(H\). \(\square\)

4. **1-bit spike for tight \(S\equiv2\) when all degrees are even and \(p=3\) (proved).** Let \(H\) be tight \(S\equiv2\), \(A=C\oplus H\). For \(y\in\mathrm{Max}_{+}\) and vertex \(v\),
   \[
   Q_A(y^{\oplus v})=\Phi-2p-4+4\sigma_v(y),\qquad \sigma_v=\sum_{vw\in H}f_{vw}(y).
   \]
   Always \(\sum_v\sigma_v(y)=4\). If every degree in \(H\) is even then each \(\sigma_v\) is even; \(\sigma_v\le0\) for all \(v\) would give \(\sum\sigma\le0<4\), so some \(\sigma_v\ge2\). For \(p=3\) this yields \(Q_A(y^{\oplus v})\ge\Phi-2\). In particular every 2-regular tight cover on a support of size \(2p\) (e.g.\ undercutting \(C_6\) at \(n=10\)) has \(\Phi\ge\Phi-2\), and the Type I equality reduction of part 3 has no-descent at \(p=3\). \(\square\)

5. **n=10 closure (proved).** Every edge-minimal undercutter is a matching (Type I, \(N_1=10>8\)) or a 6-cycle (tight \(S\equiv2\)). Parts 1–2 give no-descent; combined with Props 15.40–15.42 and the minimal-counterexample reduction of Prop 15.42(7), \(m_{10}=\Phi-2=13\). \(\square\)

6. **Global residual (OPEN).** To get \(m_n\ge\Phi-2\) for all \(\rho=1\) orders it remains to prove no-descent for:
   - Type I gap-2 undercutters at counting equality for general \(p\) (parts 3–4 handle \(p=3\); need \(\Phi\ge\Phi-2\) for all tight \(S\equiv2\) covers, or freeness of the \(k=3p-2\) boundary);
   - deep non-tight gap-2 undercutters with \(s_+\ge2\), \(k>2p\).

   If those hold, then \(m_n\ge\Phi-2=o(n^{3/2})\) on the dense \(\rho=1\) family, E(1) follows, and \(L=\tfrac12\) by Proposition 6.1 along the consecutive-prime-square orders. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.44 (master lemma for tight covers; bi-tight residual; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\Phi(C)\).

1. **Master lemma (proved).** If \(H\) is any flip set with \(S_H\equiv s\) on \(\mathrm{Max}_{+}\) for some integer \(s\ge1\) (hence \(|H|=sp\) by taking expectations), then either
   \[
   \max_{z\in\mathrm{Max}_{-}}S_H(z)\;\ge\;0
   \qquad\text{(hence \(\Phi(C\oplus H)\ge\Phi\) by Prop 15.42.1),}
   \]
   or \(S_H\equiv -s\) on \(\mathrm{Max}_{-}\) as well (**bi-tight of level \(s\)**).
   *Proof.* Always \(\mathbb E_{-}[S_H]=-|H|/p=-s\). If the maximum on \(\mathrm{Max}_{-}\) is \(\le -1\) and scores have the parity of \(s\), a maximum \(\le -s\) with mean \(-s\) forces constancy at \(-s\) when the maximum is \(\le -s\); more directly: if the maximum is \(\ge 0\) we are done by dichotomy; if the maximum is \(\le -s\) and the mean is \(-s\) with all values \(\le -s\), constancy follows. For the undercutting/gap analysis the relevant case is \(s=2\) with even scores, where maximum \(\le -2\) and mean \(-2\) yield bi-tight. \(\square\)

2. **Consequence for Type I freeness failure (proved reduction).** In the equality freeness-failure of Prop 15.43(3), \(H=G\cup\{e\}\) is tight of level \(2\). By the master lemma, either \(\Phi(C\oplus H)\ge\Phi\) (no-descent) or \(H\) is bi-tight of level \(2\). In the bi-tight subcase, Prop 15.43(4) gives \(\Phi\ge\Phi-2\) at \(p=3\); for \(p\ge5\) bi-tight level \(2\) is **integrally infeasible at \(p=5\)** (MILP: `src/e1_bitight_infeas.py`, `evidence/e1_bitight_infeas.json`, levels \(2,3,4\) all infeasible while fractional is feasible; avg degree \(4p/(p^2+1)<1\) for all \(p\ge5\)). \(\square\)

3. **Deep tight undercutters are bi-tight (proved).** A gap-\(2\) undercutter with \(S\equiv2\) on \(\mathrm{Max}_{+}\) is two-sided with mean \(-2\) on \(\mathrm{Max}_{-}\); with even scores, maximum \(\le -2\) forces \(S\equiv-2\). No-descent for such undercutters is Prop 15.43(2) (Max\(_{+}\) freeness). At \(p=5\), bi-tight is infeasible, so deep tight undercutters do not exist. \(\square\)

4. **Certified samples.** Five distinct integral tight Max\(_{+}\)-only covers of size \(10\) at \(p=5\) all have \(\max_{\mathrm{Max}_{-}}S\in\{4,8,10\}\ge0\) and exact \(\Phi\in\{73,81,85\}>\Phi(C)\), matching the master lemma. Evidence: session MILP samples; `evidence/E1_BITIGHT.md`.

5. **Residual for \(L=\tfrac12\) (OPEN).** Complete no-descent for deep **non-tight** gap-\(2\) undercutters (\(s_{\min}=2\), \(\max S\ge4\), \(k>2p\)), and lift bi-tight integral infeasibility from \(p=5\) to all \(p\ge5\) (or prove bi-tight \(\Rightarrow\Phi\ge\Phi-2\) uniformly). Then Type I and deep-tight no-descent are unconditional, and \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.45 (star classification force, bi-tight wedge block, p=5 residual kill; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), and write \(f_e(y)=C_e y_iy_j\), \(G_{ee'}=\mathbb E_{+}[f_ef_{e'}]\), and \(g_{\min}\) for the minimum of \(G_{ee'}\) over vertex-disjoint edge pairs.

1. **Wedge correlations (proved).** If edges \(va,vb\) share a vertex, then
   \[
   G_{va,vb}=\frac{C_{va}C_{vb}C_{ab}}{p}=\pm\frac1p,
   \]
   by the Max\(_{+}\) frame identity \(\mathbb E[y_ay_b]=C_{ab}/p\). On Max\(_{-}\) the same wedge has the opposite sign, so
   \[
   G^{+}_{va,vb}+G^{-}_{va,vb}=0.
   \]

2. **Star never bi-tight (proved, all \(p>2\)).** A star has only wedge pairs. For bi-tight of level \(2\) one needs both \(\sum_{e<e'}G^{+}_{ee'}=2-p\) and \(\sum G^{-}=2-p\), hence \(\sum(G^{+}+G^{-})=2(2-p)\). Wedges contribute \(0\) to \(G^{+}+G^{-}\), so a star yields sum \(0\neq 2(2-p)\). \(\square\)

3. **Correlation identity for tight covers (proved).** If \(S_H\equiv s\) on Max\(_{+}\) with \(|H|=sp\), then \(\mathbb E[S^2]=s^2\) and
   \[
   \sum_{e<e'\in H}G_{ee'}=\frac{s^2-sp}{2}.
   \]
   For \(s=1\): sum \(=(1-p)/2\), average \(=-1/p\). For \(s=2\): sum \(=2-p\), average \(=-1/15\).

4. **Star force for level-\(1\) when \(g_{\min}>-1/p\) (proved).** Write \(n_w\) for the number of wedge pairs in \(H\) and \(n_d=\binom{p}{2}-n_w\) for disjoint pairs when \(|H|=p\). Then
   \[
   \sum G\ge -\frac{n_w}p+g_{\min}n_d=g_{\min}\binom{p}{2}-n_w\Bigl(g_{\min}+\frac1p\Bigr).
   \]
   If \(g_{\min}>-1/p\), the coefficient of \(n_w\) is negative, so the lower bound is maximised at maximal \(n_w=\binom{p}{2}\) (every pair of edges shares a vertex). For \(p>3\) that forces \(H\) to be a **star**. Equality with the tight identity \(\sum G=(1-p)/2=-\binom{p}{2}/p\) holds only in that star case (up to wedge-sign pattern). At \(p=3\), \(g_{\min}=-1/p\), so the force fails and non-stars remain possible. \(\square\)

5. **Matching blocked for level-\(2\) when \(g_{\min}>-1/15\) (proved).** For \(|H|=2p\) a matching has \(n_w=0\), hence \(\sum G\ge g_{\min}\binom{2p}{2}\). If this exceeds \(2-p\), no matching is Max\(_{+}\)-tight of level \(2\). \(\square\)

6. **Certified at \(p=5\) (load-bearing numerics).** Exact Max\(\pm\) enumeration (\(|\mathrm{Max}_{\pm}|=260\)):
   - \(g_{\min}=-3/65\approx-0.04615>-1/5\) and \(>-1/15\), so level-\(1\) covers are stars and level-\(2\) matchings are non-tight;
   - integral non-star size-\(p\) Max\(_{+}\) tight covers are **MILP-infeasible**;
   - bi-tight levels \(2,3,4\) integrally infeasible (Prop 15.44);
   - deep two-sided covers (\(s_{+}\ge2\), \(s_{-}\le-2\)) integrally infeasible at \(k=10,12,15\) (timeout without feasibility at \(k=14,16,18,20\));
   - \(\min\{\max_{\mathrm{Max}_{-}}S:S_{\mathrm{Max}_{+}}\ge2,\,|H|=10\}=2\ge0\) (epigraph MILP), so every size-\(2p\) Max\(_{+}\) cover has \(\max S_{-}\ge2\) and \(\Phi\ge\Phi(C)\).
   Evidence: `src/e1_star_bitight_obstruction.py`, `e1_deep_cover_hunt.py`, `e1_deep_k_long.py`; JSON under `evidence/e1_star_bitight_obstruction.json`, `e1_deep_cover_hunt.json`, `e1_deep_k_long.json`.

7. **Certified at \(p=7\).** Full Max\(_{+}\) enumeration (\(|\mathrm{Max}_{+}|=11452\), 80-worker \(2^{25}\) free-coordinate scan): \(g_{\min}\approx-0.03807>-1/7\) and \(>-1/15\), so star-force and matching level-\(2\) block hold at \(p=7\). Evidence: `src/e1_gmin_p7.py`, `evidence/e1_gmin_p7.json`.

8. **Consequence (proved from certs + lemmas).** Stars never bi-tight (all \(p\)). At \(p=5\): Type I freeness-failure and deep-tight undercutters impossible. At \(p=5,7\): level-\(1\) Max\(_{+}\) tight covers are stars (via \(g_{\min}>-1/p\)). Residual for full \(m_n\ge\Phi-2\): bi-tight infeas for general \(p\); deep non-tight large-\(k\) control; closed-form \(g_{\min}>-1/p\) for all \(p\ge5\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.46 (1-bit spike formulas and deep/bi-tight spike criterion; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\Phi(C)\), \(A=C\oplus F\). For \(y\in\{\pm1\}^n\) write \(S=S_F(y)\) and \(\sigma_v(y)=\sum_{vw\in F}C_{vw}y_vy_w\). Always \(\sum_v\sigma_v=2S\).

1. **1-bit formula on Max\(_{+}\) (proved).** If \(Cy=py\), then
   \[
   Q_A(y^{\oplus v})
   =
   \Phi-2S-2p+4\sigma_v(y).
   \]
   Hence if \(\max_v\sigma_v(y)\ge(S+p-1)/2\), then \(Q_A(y^{\oplus v})\ge\Phi-2\), so \(\Phi(A)\ge\Phi-2\).
   *Proof.* \((Ay)_v=(Cy)_v-2\sum_{vw\in F}C_{vw}y_w=py_v-2\sum C_{vw}y_w\), so \(y_v(Ay)_v=p-2\sigma_v\). The standard 1-bit identity \(Q(x^{\oplus v})=Q(x)-2x_v(Ax)_v\) yields the claim. Threshold: \(\Phi-2S-2p+4\sigma\ge\Phi-2\Leftrightarrow\sigma\ge(S+p-1)/2\). \(\square\)

2. **1-bit formula on Max\(_{-}\) (proved).** If \(Cz=-pz\), then
   \[
   Q_A(z^{\oplus v})
   =
   -\Phi-2S+2p+4\tau_v(z),
   \]
   where \(\tau_v=\sum_{vw\in F}C_{vw}z_vz_w\). Hence if \(\min_v\tau_v(z)\le(S-p+1)/2\), then \(Q_A(z^{\oplus v})\le-(\Phi-2)\), so \(\Phi(A)\ge\Phi-2\).
   *Proof.* Symmetric to part 1 with \(Cz=-pz\). For the deep boundary \(S=-2\): threshold \(\tau\le-(p+1)/2\). \(\square\)

3. **Bi-tight even-degree spike at \(p=3\) (recovered).** If \(H\) is bi-tight of level \(2\) and every degree in \(H\) is even, then each \(\sigma_v\) is even. With \(\sum\sigma=4>0\), some \(\sigma_v\ge2=(p+1)/2\) at \(p=3\). Part 1 yields \(\Phi(A)\ge\Phi-2\). (This is Prop 15.43.4; undercutting \(C_6\) attains equality \(\sigma_{\max}=2\), \(\Phi=\Phi-2\).) \(\square\)

4. **Sufficient spike criterion (proved).** If either
   - some \(y\in\mathrm{Max}_{+}\) has \(\max_v\sigma_v(y)\ge(S_F(y)+p-1)/2\), or
   - some \(z\in\mathrm{Max}_{-}\) has \(\min_v\tau_v(z)\le(S_F(z)-p+1)/2\),
   then \(\Phi(A)\ge\Phi-2\). In particular, for deep two-sided covers with \(s_{-}=-2\), it suffices that some \(z\) with \(S_F(z)=-2\) has \(\min_v\tau_v\le-(p+1)/2\). \(\square\)

5. **Certified at \(p=5\) (deep two-sided covers).** Every MILP-found deep two-sided cover (\(k\in\{32,36,38,40\}\)) has \(s_{-}=-2\) and some Max\(_{-}\) vector at level \(-2\) with \(\min\tau\le-4\le-3=-(p+1)/2\); 1-bit yields \(|Q|\ge67\ge\Phi-2\), and exact MITM \(\Phi\in\{75,77,79,83\}>\Phi(C)\). Small-\(k\) deep two-sided (\(k\le13,15\)) integrally infeasible. Evidence: `evidence/e1_deep_cover_phi.json`, `e1_deep_sweep_p5.json`, session 1-bit checks. \(\square\)

6. **Residual for \(L=\tfrac12\) (OPEN).** Prove that every bi-tight level-\(2\) cover for \(p\ge5\) satisfies the spike criterion of part 4 (or is integrally impossible — see Prop 15.47), and that every deep two-sided gap-\(2\) undercutter (\(\Phi=\Phi-2\)) either is impossible or has no-descent. Combined with Props 15.42–15.45 this yields \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). Closed-form \(g_{\min}>-1/p\) for all \(p\ge5\) remains open (certified \(p=5,7\)). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.47 (bi-tight Gsum obstruction; 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), and write \(G^\pm_{ee'}=\mathbb E_{\pm}[f_ef_{e'}]\), \(h_{ee'}=G^+_{ee'}+G^-_{ee'}\), \(g_{\min}=\min G^+_{ee'}\) over vertex-disjoint pairs, and \(h_{\min}=\min h_{ee'}\) over the same.

1. **Wedge Gsum vanishes (proved, Prop 15.45.1).** For wedge pairs, \(h_{ee'}=0\).

2. **Bi-tight forces disj Gsum sum (proved).** If \(H\) is bi-tight of level \(2\) (\(|H|=2p\), \(S\equiv2\) on Max\(_{+}\), \(S\equiv-2\) on Max\(_{-}\)), then
   \[
   \sum_{e<e'\in H}G^+_{ee'}=2-p=\sum_{e<e'\in H}G^-_{ee'},
   \]
   hence \(\sum_{e<e'\in H}h_{ee'}=2(2-p)\). Wedges contribute \(0\), so writing \(n_d\) for the number of disjoint pairs in \(H\),
   \[
   \sum_{\substack{e<e'\in H\\e\cap e'=\emptyset}}h_{ee'}=2(2-p).
   \]

3. **Floor (proved).** Always \(h_{ee'}\ge 2g_{\min}\) (since each of \(G^\pm\ge g_{\min}\) by definition of \(g_{\min}\) and Max\(_{-}\) symmetry with \(-C\)). Hence for any \(H\) with \(|H|=2p\),
   \[
   \sum_{\mathrm{disj\ pairs\ in\ }H}h
   \;\ge\;
   h_{\min}\,n_d
   \;\ge\;
   2g_{\min}\,n_d
   \;\ge\;
   2g_{\min}\binom{2p}{2},
   \]
   where the last step uses \(g_{\min}<0\) and \(n_d\le\binom{2p}{2}\).

4. **Obstruction (proved).** If \(2g_{\min}\binom{2p}{2}>2(2-p)\), i.e.
   \[
   g_{\min}\;>\;-\frac{p-2}{p(2p-1)},
   \]
   then no bi-tight level-\(2\) cover exists. At \(p=5\) the threshold equals \(-1/15\); for \(p>5\) it is strictly larger (easier) than \(-1/15\). Independently, \(g_{\min}>-1/15\) blocks Max\(_{+}\)-tight matchings of size \(2p\) (Prop 15.45.5). \(\square\)

5. **Certified.** At \(p=5\), \(g_{\min}=-3/65>-1/15\) and \(h_{\min}=-6/65=2g_{\min}\); at \(p=7\), \(g_{\min}\approx-0.03807>-\frac{5}{91}\) and \(h_{\min}=2g_{\min}\). Both satisfy the obstruction, giving a non-MILP proof that bi-tight level \(2\) is empty. Evidence: `evidence/e1_bitight_gsum_obstruction.json`. \(\square\)

6. **Consequence for Type I (proved at \(p=5,7\); conditional for general \(p\ge5\)).** By Prop 15.44, Type I freeness-failure reduces to bi-tight or \(\Phi\ge\Phi\). With bi-tight empty under the \(g_{\min}\) threshold of part 4, Type I no-descent is unconditional. Deep tight undercutters (bi-tight) are empty. Residual: prove \(g_{\min}>-(p-2)/(p(2p-1))\) for all \(p\ge5\); deep non-tight gap-\(2\) undercutters (ND or \(\Phi\ge\Phi-2\)). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.48 (edge-correlation algebra and cross-ratio structure of \(g_{\min}\); 2026-07-28).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Phi=\tfrac12 np\), edges \(e=\{i,j\}\) with \(f_e(y)=C_{ij}y_iy_j\), and \(G_{ee'}=\mathbb E_{+}[f_ef_{e'}]\).

1. **Row sum (proved).** \(\sum_{e'}G_{ee'}=n/2\) for every \(e\). *Proof.* \(\sum_{e'}f_{e'}=\tfrac12 y^\top Cy=\Phi\) on Max\(_{+}\), and \(\mathbb E[f_e]=\tfrac1p\), so the row sum is \(\Phi/p=n/2\). \(\square\)

2. **Wedge exact values (proved).** If \(e,e'\) share a vertex then \(G_{ee'}=\pm1/p\), and for each fixed \(e\) the sum of \(G_{ee'}\) over the \(2(n-2)\) wedge partners is \(0\). *Proof.* Star identity \(\sum_{j\neq i}f_{ij}=p\) at each vertex yields \(\sum_{\mathrm{wedge\ of\ }e}f=2p-2f_e\); take \(\mathbb E[f_e\cdot(\cdot)]\). Exact values \(\pm1/p\) follow from the 2-design computation of \(E[y_ay_by_ay_c]=E[y_by_c]=C_{bc}/p\) on wedges (Prop 15.45.1 refinement). \(\square\)

3. **Disjoint mean (proved).** For each \(e\), \(\sum_{e'\,:\,e\cap e'=\emptyset}G_{ee'}=n/2-1\), hence the average disjoint correlation is \((n/2-1)/(E-1-2(n-2))\) with \(E=\binom{n}{2}\). \(\square\)

4. **Four-point pairing identity (proved).** For distinct vertices \(i,j,k,l\), writing \(\kappa=C_{ij}C_{kl}+C_{ik}C_{jl}+C_{il}C_{jk}\) and \(m_4=\mathbb E_{+}[y_iy_jy_ky_l]\),
   \[
   G_{\{ij\},\{kl\}}+G_{\{ik\},\{jl\}}+G_{\{il\},\{jk\}}
   \;=\;
   \kappa\,m_4.
   \]
   *Proof.* Each pairing contributes \(C_eC_{e'}m_4\), and the three \(C\)-products sum to \(\kappa\). \(\square\)

5. **Cross-ratio structure (certified \(p=3,5\)).** Identifying vertices with \(\mathrm{PG}(1,\mathbb F_{p^2})\), the PGL-invariant of a 4-set is its cross-ratio class. At \(p=5\): on classes with constant \(m_4\) given \(\kappa\), the three pairing correlations form the multiset \(\{-3,3,3\}/65\) or \(\{\pm1\}/65\), and \(g_{\min}=-3/65\) is achieved precisely on the \(\{-3,3,3\}/65\) classes. At \(p=3\), the analogous pattern is \(\{-1,1,1\}/3\) with \(g_{\min}=-1/3\). Two residual cross-ratio classes at \(p=5\) have non-constant \(m_4\) (binary extra invariant). Evidence: `evidence/e1_gmin_closed_form_attack.json`. \(\square\)

6. **Dead lower-bound attempts (do not reopen).** The bound \(g_{\min}\ge-3/\Phi\) holds with equality at \(p=5\) but **fails** at \(p=7\) (\(g_{\min}=-109/2863<-3/175\)). Pure 4-point boolean LP, Chebyshev on disjoint partners, Wick/Gaussian as a lower bound, bare \(C\)-isomorphism types, and affine halfspace orbits alone are all too weak or incomplete (cf. `E1_FAILURE_GRAPH.md` residual notes). \(\square\)

7. **Residual for Prop 15.47 (OPEN).** Prove \(g_{\min}>-(p-2)/(p(2p-1))\) for all primes \(p\ge5\) by a character-sum evaluation of \(m_4\) on the cross-ratio class that realises the minimum (or another scheme formula). Certified only at \(p=5,7\). Deep non-tight gap-\(2\) control remains independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.49 (cross-ratio classification of \(g_{\min}\); uniform LB candidate; 2026-07-28).** Continue the notation of Prop 15.47–15.48. Write \(n=p^2+1\), \(\Phi=\tfrac12 np\), and \(N=|\mathrm{Max}_{+}|\).

1. **Cross-ratio stratification (certified \(p=3,5,7\)).** Identifying vertices with \(\mathrm{PG}(1,\mathbb F_{p^2})\), each 4-set has a PGL-invariant cross-ratio class \(\mathrm{CR}\) and a \(C\)-invariant \(\kappa=\sum C_eC_{e'}\) over the three pairings. Evidence `e1_gmin_cr_classify.json`:
   - On every class with \(|\kappa|=1\) and constant \(m_4\), one has \(m_4=\kappa\alpha\) for a class constant \(\alpha>0\), and the three pairing correlations form \(\{-\alpha,\alpha,\alpha\}\). Hence the pairing minimum is \(-\alpha\).
   - Global \(g_{\min}=-\alpha_\star\) where \(\alpha_\star\) is the maximum of such \(\alpha\) over constant-\(m_4\) classes with \(|\kappa|=1\).
   - Values: \(g_{\min}=-\tfrac13\) (\(p=3\)); \(-\tfrac3{65}\) (\(p=5\)); \(-\tfrac{109}{2863}\) (\(p=7\)). At \(p=5,7\), \(g_{\min}\) strictly exceeds the bi-tight threshold of Prop 15.47; at \(p=3\) it does not (consistent with bi-tight \(C_6\)).

2. **Spectral rank (certified \(p=3,5,7\)).** The Gram \(G\) has \(\mathrm{rank}(G)=\binom{d}{2}-d+1\) with \(d=n/2\), and a simple eigenvalue \(n/2\) for the all-ones vector on edges. At \(p=5\) the nonzero spectrum is \(\{n/2,88/13,72/13,40/13\}\); at \(p=7\) it is \(\{n/2\}\) union five positive eigenvalues with denominator \(409\) and multiplicities \(\{d,n,n,2n,n\}\).

3. **Uniform lower-bound candidate (algebra + certification).** Define
   \[
   L(p)\;:=\;-\frac{p-2}{2p^2}.
   \]
   For every odd integer \(p>2\),
   \[
   L(p)\;>\;-\frac{p-2}{p(2p-1)}
   \]
   (clear the positives: \(2p^2>p(2p-1)\Leftrightarrow 0>-1\)). Thus \(L(p)\) lies strictly above the bi-tight threshold. Certified: \(g_{\min}(5)\ge L(5)\) and \(g_{\min}(7)\ge L(7)\); while \(g_{\min}(3)<L(3)\) (so the candidate is not valid at \(p=3\), as required by the existence of bi-tight undercutters there). Evidence: `e1_gmin_uniform_lb.json`.

4. **Matching margin (certified \(p=5\)).** For every matching \(M\) of size \(2p=10\) on Paley \(C_{26}\), the quadratic form \(\mathbf{1}_M^\top G\mathbf{1}_M\) is at least \(9.96>4\), so \(\sum_{\mathrm{pairs\ in\ }M}G\ge -0.02>2-p=-3\). Equality in the PSD projection bound \(\mathbf{1}_M^\top G\mathbf{1}_M\ge4\) is therefore not attained on matchings (sampling + local search, \(2\cdot10^4\) seeds).

5. **Residual for Prop 15.47 (OPEN).** Prove \(g_{\min}\ge L(p)\) for every prime \(p\ge5\) (or any other uniform lower bound strictly above the bi-tight threshold). Combined with the deep non-tight residual of Prop 15.47.6 this yields \(m_n\ge\Phi-2\Rightarrow\mathrm{E}(1)\Rightarrow L=\tfrac12\). **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.50 (conditional means of Max+; disjoint correlation mean; 2026-07-29).** Let \(C\) be \(\rho=1\) Paley of order \(n=p^2+1\), \(\Sigma=I+C/p=2P_+\), and let \(y\) be uniform in \(\mathrm{Max}_{+}\). Write \(f_e=C_e y_iy_j\) and \(G_{ee'}=\mathbb E[f_ef_{e'}]\).

1. **Conditional mean lemma (proved).** For distinct indices \(i,j\) and \(a,b\in\{\pm1\}\),
   \[
   \mathbb E\bigl[y\bigm|\,y_i=a,\,y_j=b\bigr]
   \;=\;
   \Sigma_{*,S}\,\Sigma_{S,S}^{-1}\begin{pmatrix}a\\b\end{pmatrix},
   \qquad S=\{i,j\}.
   \]
   *Proof.* The right-hand side \(\mu_*\) is the unique minimum-norm vector in \(V_+=\mathrm{range}(P_+)\) with those two coordinates (Gaussian interpolant for covariance \(\Sigma\)). The left-hand side \(\mu\) is an average of Max\(+\) vectors, hence lies in \(V_+\) with the same coordinates. The difference \(v=\mu-\mu_*\) satisfies \(v\in V_+\) and \(v_i=v_j=0\). For every \(w\in V_+\) with \(w_i=w_j=0\), the scalar \(y\cdot w\) is orthogonal (under \(\mathbb E\)) to \(\{1,y_i,y_j,y_iy_j\}\): the first three pairings use \(\mathbb E[y]=0\) and \(\mathbb E[yy^\top]=\Sigma\); the fourth is an odd third moment and vanishes by central symmetry \(\mathrm{Max}_+=-\mathrm{Max}_+\). Those four monomials span all functions of \((y_i,y_j)\), so \(\mathbb E[y\cdot w\mid y_i,y_j]=0\). Hence \(\mu\cdot w=0\) for all such \(w\). In particular \(\mu_*\) lies in \(\mathrm{span}\{P_+e_i,P_+e_j\}\) and \(v\perp\) that span, so \(\mu_*\cdot v=0\) and \(\mu\cdot v=\|v\|^2=0\), whence \(v=0\). \(\square\)
   Certified at \(p=5,7\) (max abs error \(<10^{-15}\)): `evidence/e1_gmin_cond_mean.json`, `src/e1_gmin_cond_mean.py`.

2. **Conditional second-moment shape (proved).** For fixed distinct \(i,j,k,l\),
   \[
   \mathbb E[y_ky_l\mid y_i,y_j]
   \;=\;
   \alpha+\delta\,y_iy_j
   \]
   for scalars \(\alpha,\delta\) determined by \(\Sigma_{kl}\) and \(m_4=\mathbb E[y_iy_jy_ky_l]\):
   \[
   \alpha=\frac{\Sigma_{kl}-c\,m_4}{1-c^2},\qquad
   \delta=\frac{m_4-c\,\Sigma_{kl}}{1-c^2},\qquad
   c=\Sigma_{ij}=C_{ij}/p.
   \]
   *Proof.* The conditional expectation is a function of \((y_i,y_j)\), hence of the form \(\alpha+\beta y_i+\gamma y_j+\delta y_iy_j\). Matching moments against \(1,y_i,y_j,y_iy_j\) and using vanishing odd moments forces \(\beta=\gamma=0\) and the displayed formulae. \(\square\)

3. **Disjoint-pair mean (proved).** For every edge \(e\), \(\sum_{e'\,:\,e\cap e'=\emptyset}G_{ee'}=n/2-1\), so the average disjoint correlation is
   \[
   \frac1{p^2-2}.
   \]
   *Proof.* Row sum of \(G\) is \(n/2\) (Prop 15.48.1); sum of \(G\) over the \(2(n-2)\) wedge partners of \(e\) is \(0\) (Prop 15.48.2); the remainder is the disj sum. Count of disj partners \((n-2)(n-3)/2\) and \(n=p^2+1\) give the average. \(\square\)

4. **Frechet obstruction is too weak (certified).** Combining part 1 with the Fréchet–Hoeffding bound
   \(\mathbb E[y_ky_l\mid\mathrm{state}]\ge\mu_k\mu_l-\sqrt{(1-\mu_k^2)(1-\mu_l^2)}\)
   yields only \(g_{\min}\ge-O(1)\) (empirically \(\ge-0.4\) at \(p=5\)), which does **not** beat the bi-tight threshold \(T(p)=-(p-2)/(p(2p-1))\). Evidence in session residual notes. **Do not reopen plain Fréchet as a path to \(L(p)\).**

5. **Residual (OPEN).** Prove \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) (or any LB \(>T(p)\)) for all primes \(p\ge5\). The conditional-mean calculus reduces this to a uniform upper bound on \(|m_4|\) (or on \(\delta\)) for four-sets with \(|\kappa|=1\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.51 (equivalent bi-tight threshold form; residual slice structure; 2026-07-29).** Continue the notation of Prop 15.47–15.50. Fix an edge \(e=(i,j)\) and write \(a(e'):=\mathbb E[f_{e'}\mid f_e=1]\) for \(e'\) vertex-disjoint from \(e\).

1. **Identity (proved).** For any disj pair \(e,e'\),
   \[
   a(e')=\frac{1+p\,G_{ee'}}{p+1}.
   \]
   *Proof.* From Prop 15.50, \(\mathbb E[y_ky_l\mid y_iy_j=s]=\alpha+\delta s\) with \(\alpha,\delta\) linear in \(m_4=G_{ee'}/(C_eC_{e'})\); conditioning on \(f_e=1\) forces \(y_iy_j=C_e\) and yields the display after multiplying by \(C_{e'}\). Equivalently, expand \(G=\mathbb E[f_ef_{e'}]\) on the two values of \(f_e\) with \(\mathbb E[f_e]=1/p\). \(\square\)

2. **Equivalent bi-tight threshold (proved).** Therefore
   \[
   g_{\min}\ge T(p):=-\frac{p-2}{p(2p-1)}
   \quad\Longleftrightarrow\quad
   \min_{e'\,:\,e\cap e'=\emptyset}a(e')\ge\frac1{2p-1}
   \]
   for every edge \(e\) (and hence, by edge-transitivity of \(\mathrm{Aut}\), for one fixed \(e\)). *Proof.* Clear the linear identity of part 1. \(\square\)
   Certified: at \(p=5\), \(\min a=5/39>1/9\); at \(p=7\), \(\min a\approx0.0917>1/13\). Evidence: `e1_gmin_structure.json`.

3. **Deterministic disj sum (proved).** Pointwise on \(\mathrm{Max}_{+}\),
   \[
   \sum_{e'\,:\,e\cap e'=\emptyset}f_{e'}=\Phi-2p+f_e.
   \]
   *Proof.* \(\sum_{\mathrm{all}\,e'}f=\Phi\) and \(\sum_{\mathrm{wedge\,of\,}e}f=2p-2f_e\) (stars at the two ends of \(e\)). \(\square\)
   Consequently on the slice \(f_e=1\) the average of \(f_{e'}\) over disj partners is the constant \((\Phi-2p+1)/D\) with \(D=(n-2)(n-3)/2\).

4. **Residual Loewner (certified \(p=5,7\)).** On the slice \(f_e=1\), the residual Gram \(R=\mathbb E[(y-\mu)(y-\mu)^\top\mid f_e=1]\) (mixture of the two states with \(f_e=1\)) satisfies \(R\succeq\lambda_{\min}(R)\,P_W\) where \(P_W\) is the orthogonal projector onto \(V_+\cap\{x_i=x_j=0\}\) and \(\mathrm{rank}(R)=n/2-2\). Evidence: `e1_gmin_structure.json`. **Not yet a proof of \(g_{\min}>T(p)\):** entrywise Schur bounds on the residual still undershoot \(1/(2p-1)\).

5. **Max+ types (certified).** At \(p=5\), \(\mathrm{Max}_{+}\) is distance-homogeneous (constant Hamming distance distribution from every vector). At \(p=7\), at least two distance types occur (so \(\mathrm{Max}_{+}\) is not a single \(\mathrm{Aut}\)-orbit). Affine+\(\mathrm{PGL}\) orbit of the halfspace vector has size \(60\) of \(260\) at \(p=5\) — character sums on that orbit alone do **not** compute full-Max+ \(m_4\).

6. **Residual (OPEN).** Prove \(\min a(e')\ge1/(2p-1)\) (equivalently \(g_{\min}\ge T(p)\), or the stronger \(g_{\min}\ge L(p)\)) for all primes \(p\ge5\). Preferred routes: character-sum / BM algebra on the min CR class; residual \(z=y-\mu\) with the pointwise identity \(z_r^2+2\mu_rz_r=1-\mu_r^2\) beyond Fréchet; Loewner calculus with a sharp entrywise bound; **m4 moduli** (evec system on refined \(C\)-classes has nullity 1; pin by \(\mathrm{Tr}(G^2)\), recovers \(g_{\min}\) at \(p=5\) — `E1_GMIN_MODULI.md`). **Do not** reopen F15 (plain Fréchet) or incomplete Aut-orbits as full Max+. Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.52 (coordinate sum of Max+; m4 moduli sketch; 2026-07-29).** Let \(C\) be the Paley conference matrix of order \(n=p^2+1\) with vertices \(\{\infty\}\cup\mathbb F_{p^2}\), and \(y\in\mathrm{Max}_{+}\).

1. **Coordinate sum (proved).** \(\mathbf1^\top y=(p+1)y_\infty\). In particular \(|\mathbf1^\top y|=p+1\).
   *Proof.* Row sums of \(C\): \((C\mathbf1)_\infty=p^2\) and \((C\mathbf1)_v=1\) for \(v\in\mathbb F_{p^2}\) (complete character sum \(\sum_{d\neq0}\chi(d)=0\)). Thus \(C\mathbf1=(p^2-1)e_\infty+\mathbf1\). From \(Cy=py\), take \(\mathbf1^\top\): \(y^\top C\mathbf1=p\,s\) with \(s=\mathbf1^\top y\), so \((p^2-1)y_\infty+s=ps\), hence \(s=(p+1)y_\infty\). \(\square\)
   Certified \(p=5,7\): every Max+ vector has sum \(\pm(p+1)\).

2. **m4 linear system (certified \(p=5,7\)).** Stratify 4-sets by \((\mathrm{CR},\kappa,\triangle\text{-type})\) (all \(C\)-invariants; \(\triangle\)-type splits formerly non-constant classes). Averaged evec identities give \((pI-M)\mathbf m=\mathbf b\) with \(M,b\) combinatorial. The system has **nullity 1**. The true Max+ moment vector lies on this line. The constraint \(\mathrm{Tr}(G^2)=E+2n_{\mathrm{wedge}}/p^2+6\sum n_A m_A^2\) is quadratic in the free parameter and **selects the true \(m_4\) at \(p=5\)** (hence \(g_{\min}=-3/65\)). Evidence: `E1_GMIN_MODULI.md`.

3. **Residual (OPEN).** Prove nullity 1 and a Max+-free evaluation of \(\mathrm{Tr}(G^2)\) (or of the \(G\)-spectrum) for all primes \(p\ge5\); solve for \(g_{\min}\) and show \(g_{\min}\ge L(p)\) or \(>T(p)\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.53 (pairing reduction of \(g_{\min}\); moduli pin at \(p=5\); 2-design \(\mathrm{Tr}(G^2)\) skeleton; 2026-07-30).** Continue Prop 15.52.

1. **Pairing identity (proved).** For any 4-set with \(\kappa:=\sum_{\text{three pairings}}C_eC_{e'}\) satisfying \(|\kappa|=1\), the three pairing products are a permutation of \((1,1,-1)\) or \((-1,-1,1)\). Hence the three edge-pair correlations equal \(m_4\cdot(\pm1)\) in that pattern, and
   \[
   \min_{\text{three pairings}}G_{ee'}=-\lvert m_4\rvert.
   \]
   Therefore
   \[
   g_{\min}=-\max\bigl\{\lvert m_4(S)\rvert:S\subset V,\,|S|=4,\,|\kappa(S)|=1\bigr\}.
   \]
   *Proof.* Each pairing product is \(\pm1\); their sum is \(\kappa=\pm1\) forces the stated multiset. Then \(G=C_eC_{e'}m_4\) on each pairing, so the minimum is \(-\lvert m_4\rvert\). Every disjoint edge pair sits in a unique 4-set, so the global \(g_{\min}\) is the min over such 4-sets. \(\square\)
   Certified \(p=5,7\): identity error \(<10^{-12}\); \(g_{\min}=-3/65\), \(-109/2863\). Evidence: `e1_gmin_moduli.json`.

2. **Refined \(C\)-classes (certified \(p=5\)).** Stratify 4-sets by the pure \(C\)-invariant
   \((\mathrm{type}_6,\mathrm{ext}\text{-sum histogram})\), where \(\mathrm{type}_6\) is the \(S_4\)-canonical 6-tuple of edge signs and the external histogram records \(\sum_{v\in S}C_{rv}\) for \(r\notin S\). At \(p=5\) this yields **37** classes, each with **constant** \(m_4\) on Max+. (Bare \(C\)-types alone do **not** make \(m_4\) constant — F-graph / prior notes.)

3. **Nullity-1 evec system (certified \(p=5\)).** Averaging \(p\,m_4=\sum_r C_{ir}m(\cdots)\) over each class produces a combinatorial linear system \(A\mathbf m=\mathbf b\) (RHS uses only the 2-design \(m_2=C/p\)). At \(p=5\): \(\mathrm{rank}(A)=36=n_{\mathrm{var}}-1\). The true Max+ moment vector lies on the affine line \(\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n\).

4. **\(\mathrm{Tr}(G^2)\) pin (certified \(p=5\)).** Write \(K_{ab}=\bigl((y_a\cdot y_b)^2-n\bigr)/2\). Then the nonzero spectrum of \(G\) matches that of \(K/N\), and
   \[
   \mathrm{Tr}(G^2)=\frac1{N^2}\sum_{a,b}K_{ab}^2=\tfrac14\Bigl(\mathbb E[\mathrm{dot}^4]-2n\,\mathbb E[\mathrm{dot}^2]+n^2\Bigr).
   \]
   Substituting \(\mathbf m(c)\) into the edge form of \(\mathrm{Tr}(G^2)\) yields a quadratic in \(c\). One root recovers \(g_{\min}=-3/65\); select the root of larger \(g_{\min}\) among the two (do **not** use PSD-max over the whole line — F16). Evidence: `e1_gmin_moduli.json`, `E1_GMIN_MODULI.md`.

5. **2-design evaluation of \(\mathbb E[\mathrm{dot}^2]\) (proved, Max+-free beyond the frame).** From \(\mathbb E[yy^\top]=I+C/p\),
   \[
   \mathbb E_{a,b}[(y_a\cdot y_b)^2]=\|I+C/p\|_F^2=n+\frac{n(n-1)}{p^2}.
   \]
   *Proof.* \(\sum_{a,b}(y_a\cdot y_b)^2=\sum_{i,j}(\sum_a y_{a,i}y_{a,j})^2=N^2\|I+C/p\|_F^2\). \(\square\)
   Certified \(p=5,7\). **Still Max+-dependent:** \(\mathbb E[\mathrm{dot}^4]\) (equivalently a closed \(G\)-spectrum), which is needed to evaluate \(\mathrm{Tr}(G^2)\) without Max+ samples.

6. **Wick comparison (certified \(p=5,7\); not yet a proof).** For \(\Sigma=I+C/p\), the Gaussian fourth moment is \(\mathbb E_{\mathrm{Wick}}[\mathrm{dot}^4]=3\|\Sigma\|_F^4+6\mathrm{Tr}(\Sigma^4)\). Discrete Max+ satisfies \(\mathbb E[\mathrm{dot}^4]<\mathbb E_{\mathrm{Wick}}[\mathrm{dot}^4]\) at \(p=5,7\) (boolean coordinates have smaller kurtosis than Gaussians with the same covariance). On the combinatorial moduli line at \(p=5\), the weaker constraint \(\mathrm{Tr}(G^2)\le\mathrm{Tr}_{\mathrm{Wick}}\) forces \(g_{\min}\ge T(p)\) with **endpoint equality** (float margin \(\sim10^{-15}\)) — not the strict \(g_{\min}>T(p)\) needed for Prop 15.47. A proved strict inequality \(\mathbb E[\mathrm{dot}^4]\le\mathbb E_{\mathrm{Wick}}-\delta\) (or the exact spectrum) would upgrade this to \(g_{\min}>T\). Evidence: `e1_gmin_moduli.json`.

7. **Residual (OPEN).** Prove for all primes \(p\ge5\): (i) refined classes have constant \(m_4\) and the evec system has nullity 1; (ii) Max+-free strict bound \(\mathbb E[\mathrm{dot}^4]<\mathbb E_{\mathrm{Wick}}\) or a closed \(G\)-spectrum; (iii) the selected root satisfies \(g_{\min}\ge L(p)=-(p-2)/(2p^2)\) or at least \(g_{\min}>T(p)\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.54 (moduli \(c\)-pin calculus; wedge closed form; a-slice certs; 2026-07-30).** Continue Prop 15.53.

1. **Wedge \(G\) (proved, combinatorial).** If edges \(e=(i,j)\) and \(e'=(i,k)\) share a vertex, then
   \[
   G_{ee'}=\frac{C_{ij}C_{ik}C_{jk}}{p}=\pm\frac1p.
   \]
   *Proof.* \(f_ef_{e'}=C_{ij}C_{ik}y_jy_k\) and \(\mathbb E[y_jy_k]=C_{jk}/p\). \(\square\)

2. **Moduli line and \(g_{\min}(c)\) (certified \(p=5\)).** On the nullity-1 line \(\mathbf m=\mathbf m_{\mathrm{part}}+c\,\mathbf n\) of Prop 15.53, for every \(|\kappa|=1\) class
   \(g_{\min}(c)=-\max_A|m_A(c)|\). The edge Gram is \(G(c)=G_{\mathrm{wedge}}+G_{\mathrm{disj}}(\mathbf m(c))\) with wedges from part 1. Then \(\mathrm{Tr}(G(c)^2)=a_0+a_1c+a_2c^2\) with combinatorial \(a_i\). Evidence: `e1_gmin_cbound.json`.

3. **True \(\mathrm{Tr}(G^2)\) pin (certified \(p=5\)).** Setting \(\mathrm{Tr}(G(c)^2)=\mathrm{Tr}(G_{\mathrm{Max+}}^2)\) yields two roots; the root of larger \(g_{\min}\) recovers \(g_{\min}=-3/65>T(5),L(5)\). The Wick fourth-moment vector \(m_4=\kappa/p^2\) is **not** exactly on the evec line (residual \(\sim10^{-2}\)–\(10^{-1}\) relative). Evidence: `e1_gmin_cbound.json`.

4. **Slice \(a_{\min}\) (certified \(p=5,7\)).** \(\min a=5/39>1/9\) at \(p=5\) and \(75/818>1/13\) at \(p=7\). Wick-\(a\) overestimates \(\min a\); mean-only \(C_{kl}\mu_k\mu_l\) undershoots (residual \(R\) helps). Evidence: `e1_gmin_abound.json`. **Not a uniform proof.**

5. **Deep covers at \(p=5\) (certified refresh).** MILP deep two-sided covers at \(k\in\{32,36,38,40\}\) all meet the Prop 15.46 spike criterion and have \(\Phi\ge\Phi(C)\). Evidence: `e1_deep_spike_theory.json`. Uniform deep ND for all \(p\ge5\) still open.

6. **Residual (OPEN).** Same as Prop 15.53.7: Max+-free \(\mathrm{Tr}(G^2)\) / spectrum for general \(p\), hence \(g_{\min}>T(p)\) for all primes \(p\ge5\); deep non-tight ND. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.55 (tight size-\(2p\) obstruction from \(\lambda_{\max}(G)=n/2\); 2026-07-30).** Let \(G=\mathbb E_{+}[ff^\top]\) on the edge space of Paley \(C\) of order \(n=p^2+1\).

1. **Row sum (proved).** \(G\mathbf1=(n/2)\mathbf1\). *Proof.* \(\sum_e f_e=\Phi=pn/2\) and \(\mathbb E[f_e]=1/p\), so each row sums to \(\Phi/p=n/2\). \(\square\)

2. **All-ones mass of a size-\(2p\) indicator (proved).** For any \(0\)-\(1\) vector \(v\) with \(\sum v_e=2p\),
   \[
   v^\top\Bigl(\tfrac n2\cdot\frac{\mathbf1\mathbf1^\top}{E}\Bigr)v=4,
   \]
   where \(E=\binom{n}{2}\). *Proof.* Direct: \(\tfrac n2\cdot(2p)^2/E=4p^2/(n-1)=4p^2/p^2=4\). \(\square\)

3. **Tight cover forces \(G_\perp\)-isotropy (proved).** If \(H\) is Max\(_{+}\)-tight of level \(2\) (\(|H|=2p\), \(S_H\equiv2\) on \(\mathrm{Max}_{+}\)), then \(v:=\mathbf1_H\) satisfies \(v^\top Gv=\mathbb E[S_H^2]=4\). Writing \(G=\tfrac n2 P_{\mathbf1}+G_\perp\) with \(P_{\mathbf1}=\mathbf1\mathbf1^\top/E\), part 2 yields \(v^\top G_\perp v=0\).

4. **Retraction of the claimed obstruction.** Assume \(G\succeq0\) and \(\lambda_{\max}(G)=n/2\) with multiplicity one, and put \(G_\perp=G-(n/2)P_{\mathbf1}\). Although \(G_\perp\succeq0\), the former assertion \(\ker G_\perp=\mathrm{span}\{\mathbf1\}\) is false. In fact
   \[
   \ker G_\perp=\mathrm{span}\{\mathbf1\}\oplus(\ker G\cap\mathbf1^\perp),
   \]
   and Proposition 15.56.1 immediately supplies the \((n-2)\)-dimensional star-difference subspace in \(\ker G\cap\mathbf1^\perp\). Thus \(v^\top G_\perp v=0\) only places the centered indicator in \(\ker G\); it does not make \(v\) constant and does not exclude a tight cover. The former conclusion of this part, and every downstream use of it in 15.167--15.168, is retracted. Proposition 15.720 gives a valid discrete obstruction for the required bi-tight levels. \(\square\)

5. **Certified spectrum.** At \(p=5,7\): \(\lambda_{\max}(G)=n/2\) is simple (next eigenvalues \(\approx6.77,5.28\ll n/2\)). At \(p=3\): \(\lambda_{\max}(G)=8>n/2=5\) (multiplicity \(5\)), so the obstruction does **not** apply — consistent with bi-tight \(C_6\). Evidence: `e1_gmin_tight_obstruction.json`.

6. **Corrected status.** Proving \(\lambda_{\max}(G)=n/2\) simple does not by itself exclude tight or bi-tight covers. The spectral route is therefore not an E(1) gate. The required level-2 and level-3 bi-tight alternatives are handled instead by Proposition 15.720. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.56 (star/cycle decomposition and Schur-square reduction; 2026-07-30).** Continue Prop 15.55. Write \(E=\binom{n}{2}\), \(d=n/2\), \(N=|\mathrm{Max}_{+}|\), and let \(Y\) be the \(N\times n\) matrix of Max+ vectors.

1. **Star action (proved).** Let \(u^{(i)}\in\mathbb R^E\) be the indicator of edges incident to vertex \(i\). Then \(u^{(i)\top}f=p\) for every Max+ feature vector \(f\) (star sum), hence
   \[
   Gu^{(i)}=\mathbf1\qquad\text{for all }i.
   \]
   Consequently \(G(u^{(i)}-u^{(j)})=0\), so \(\ker G\) contains the \((n-2)\)-dimensional space of star differences orthogonal to the all-ones edge vector. *Proof.* \(f\cdot u^{(i)}=p\) constantly, so \(G u^{(i)}=\mathbb E[f\,p]=p\cdot\mathbb E[f]=\mathbf1\). \(\square\)

2. **Cycle reduction (proved).** Let \(\mathrm{Cyc}=\{v\in\mathbb R^E:\sum_{j\neq i}v_{ij}=0\ \forall i\}\) (cycle space, dim \(E-n+1\)). Then
   \[
   \lambda_{\max}(G)=\max\bigl(n/2,\ \lambda_{\max}(G|_{\mathrm{Cyc}})\bigr).
   \]
   *Proof.* The star space has dimension \(n-1\), carries the eigenvalue \(n/2\) along \(\mathbf1\), and contributes \(n-2\) kernel directions by part 1; all remaining positive spectrum lies in \(\mathrm{Cyc}\). \(\square\)

3. **Schur-square dual (proved).** Let \(D=YY^\top\) and \(P=D/(2N)\). Then \(P\) is the orthogonal projector of rank \(d\) with constant diagonal \(P_{aa}=d/N\), and the nonzero eigenvalues of \(G\) other than \(n/2\) equal \(2N\) times the eigenvalues of \(P\odot P\) on \(\mathbf1^\perp\subset\mathbb R^N\). In particular
   \[
   \lambda_{\max}(G|_{\mathrm{Cyc}})=2N\cdot\lambda_2(P\odot P),
   \]
   so
   \[
   \lambda_{\max}(G)=\tfrac n2
   \quad\Longleftrightarrow\quad
   \lambda_2(P\odot P)\le \frac d{2N}.
   \]
   *Proof.* \(D\) has spectrum \(2N\) (\(d\) times) and \(0\); \(P=QQ^\top\) for \(Q^\top Q=I_d\). Feature Gram identity \(K_{ab}=((y_a\cdot y_b)^2-n)/2\) gives \(K=2N^2(P\odot P)-\frac n2 J\), and nonzero eigenvalues of \(G\) match those of \(K/N\). On \(\mathbf1^\perp\), \(J=0\), yielding the factor \(2N\). \(\square\)

4. **Average cycle eigenvalue (proved for \(p\ge5\)).** Writing \(k=\binom{d-1}{2}-1\) for the cycle rank,
   \[
   \frac1k\sum_{\mathrm{cycle}}\lambda_j
   =\frac{n(n-2)/2}{k}
   <\frac n2
   \quad\text{for all primes }p\ge5.
   \]
   *Proof.* Algebra: \(n-1\le(d-1)(d-2)/2\) rearranges to the claim; holds for \(n=p^2+1\ge26\). \(\square\)
   (Average \(<n/2\) is necessary but not sufficient for \(\lambda_{\max}(\mathrm{cycle})\le n/2\).)

5. **Certified spectral gap.** At \(p=5,7\): \(\lambda_2(P\odot P)\le d/(2N)\) (hence \(\lambda_{\max}(G)=n/2\) simple). At \(p=3\): \(\lambda_2=1/3>5/24=d/(2N)\) (hence \(\lambda_{\max}=8>5\)). Evidence: `e1_gmin_spectral.json`.

6. **Residual (OPEN).** Prove \(\lambda_2(P\odot P)\le d/(2N)\) for every prime \(p\ge5\) (equivalently \(\lambda_{\max}(G|_{\mathrm{Cyc}})\le n/2\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight residual independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.57 (triangle form bound; crude \(\lambda_{\max}(G)\le n\); 2026-07-30).** Let \(v\in\mathbb R^E\) with \(\|v\|=1\), and write \(B_{ij}=C_{ij}v_{ij}\) (\(i\neq j\)), \(B_{ii}=0\). Define the **triangle form**
\[
\mathrm{ft}(v)=\sum_{i=1}^n\sum_{j\neq i}\sum_{k\neq i}C_{ij}C_{ik}C_{jk}\,v_{ij}v_{ik}=\sum_{i}w^{(i)\top}C^{[i]}w^{(i)},
\]
where \(w^{(i)}_j=C_{ij}v_{ij}\) and \(C^{[i]}\) is the principal submatrix of \(C\) deleting row/column \(i\).

1. **Triangle bound (proved).** \(|\mathrm{ft}(v)|\le 2p\).
   *Proof.* For each \(i\), the Rayleigh principle gives \(|w^\top C^{[i]}w|\le\|C\|_{\mathrm{op}}\|w\|^2=p\,\delta_i\) with \(\delta_i=\sum_{j\neq i}v_{ij}^2\), because \(C^{[i]}\) is the restriction of the symmetric matrix \(C\) to a coordinate subspace (or directly: extend \(w\) by a zero at \(i\) and use \(\|C\|_{\mathrm{op}}=p\)). Summing over \(i\) and using \(\sum_i\delta_i=2\|v\|^2=2\) yields the claim. \(\square\)
   Certified: on the cycle space, the spectrum of the triangle form is exactly \(\{\pm 2p\}\) at the extremes for \(p=3,5,7\). Evidence: `e1_gmin_spectral.json` / session residual notes.

2. **Second-moment identity (proved).** For \(y\in\mathrm{Max}_{+}\),
   \[
   \mathbb E[\|By\|^2]=2+\frac{\mathrm{ft}(v)}{p}\le 4.
   \]
   *Proof.* Expand \(\mathbb E[y^\top B^\top B y]=\mathrm{Tr}(B^\top B\Sigma)\) with \(\Sigma=I+C/p\) and \(y_i^2=1\); the cross term is \(\mathrm{ft}/p\). Apply part 1. \(\square\)

3. **Crude operator bound (proved).** \(\lambda_{\max}(G)\le n\).
   *Proof.* Cauchy--Schwarz: \((y^\top By)^2\le n\|By\|^2\). Take \(\mathbb E\), use part 2, and \(f\cdot v=y^\top By/2\). \(\square\)
   (Factor-of-two away from the target \(\lambda_{\max}\le n/2\); the maximising cycle direction has \(\mathrm{ft}=2p\) and \(By\in V_+\) for all Max+ \(y\), certified \(p=3,5,7\).)

4. **Residual (OPEN).** Improve part 3 by a factor \(2\) on the cycle space for \(p\ge5\) (equivalently Prop 15.56.6: \(\lambda_2(P\odot P)\le d/(2N)\)). Then bi-tight closes for all primes \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.58 (Max+ in \(V_+\); Perron of \(P\odot P\); Veronese residual form; 2026-07-30).** Continue Prop 15.56–15.57. Write \(d=n/2\), \(N=|\mathrm{Max}_{+}|\), \(\alpha=d/N\), \(Y\) the \(N\times n\) Max+ matrix, \(P=YY^\top/(2N)\), and \(u_a=y_a/\sqrt n\in V_+\) (unit vectors).

1. **Max+ lies in \(V_+\) (proved).** For every \(y\in\mathrm{Max}_{+}\), \(Cy=py\). *Proof.* \(\Phi(y)=\tfrac12 y^\top Cy\) and the Max+ level is \(\Phi=pn/2\), so \(y^\top Cy=pn\). With \(\|y\|^2=n\) and \(C=p(P_+-P_-)\) one gets \(\|P_+y\|^2-\|P_-y\|^2=n\) and \(\|P_+y\|^2+\|P_-y\|^2=n\), hence \(P_-y=0\). \(\square\)
   Consequently \(Y^\top Y=2N\,P_+\) (tight frame of Max+ in \(V_+\)).

2. **Perron of the Schur square (proved, any equal-diagonal orthoprojector).** \((P\odot P)\mathbf1=\alpha\mathbf1\) and \(\lambda_{\max}(P\odot P)=\alpha\). *Proof.* \(\sum_b P_{ab}^2=P_{aa}=\alpha\) (from \(P^2=P\)), so row sums of \(P\odot P\) equal \(\alpha\). Entries of \(P\odot P\) are nonnegative, so Perron–Frobenius gives \(\lambda_{\max}=\alpha\) with eigenvector \(\mathbf1\). \(\square\)
   The residual \(\lambda_2(P\odot P)\le\alpha/2\) is therefore a gap below the Perron root (equivalently \(\lambda_{\max}(G)=n/2\)).

3. **Veronese / Gram reformulation (proved equivalent).** For \(x\in\mathbb R^N\) set \(T(x)=\sum_{a=1}^N x_a\,y_ay_a^\top\in\mathrm{Sym}_n\). Then
   \[
   \lambda_2(P\odot P)\le\frac\alpha2
   \quad\Longleftrightarrow\quad
   \|T(x)\|_F^2\le nN\,\|x\|^2
   \quad\text{for all }x\perp\mathbf1.
   \]
   *Proof.* Expand \(\|T(x)\|_F^2=\sum_{a,b}x_ax_b(y_a\cdot y_b)^2\). With \(P_{ab}=(y_a\cdot y_b)/(2N)\) one has
   \((y_a\cdot y_b)^2=4N^2(P\odot P)_{ab}\), so \(\|T\|_F^2=4N^2\,x^\top(P\odot P)x\). On \(\mathbf1^\perp\), \(x^\top(P\odot P)x\le\lambda_2(P\odot P)\|x\|^2\), and the target \(\lambda_2\le\alpha/2=d/(2N)\) rearranges to \(\|T\|_F^2\le nN\|x\|^2\) (using \(d=n/2\)). The reverse direction is the Rayleigh quotient for \(\lambda_2\). \(\square\)
   Equivalently, writing \(W_{ab}=(u_a\cdot u_b)^2\), one has \(\lambda_1(W)=N/d\) and the residual is \(\lambda_2(W)\le N/n=\lambda_1(W)/2\).

4. **Zero diagonal of \(T(x)\) on \(\mathbf1^\perp\) (proved).** If \(\sum_a x_a=0\) then \(T(x)_{ii}=\sum_a x_a=0\) for all \(i\) (using \(y_{a,i}^2=1\)). Thus \(T(x)\) is a zero-diagonal, trace-zero matrix supported through \(V_+\) (part 1).

5. **Maximiser structure (certified \(p=3,5,7\); partial proof).** On a top cycle eigenvector of \(G\): \(\mathrm{ft}(v)=2p\), \(\mathbb E[\|By\|^2]=4\), \(By\in V_+\) for every Max+ \(y\), and \(\|P_+BP_+\|_F^2=2\). The identity \(\mathbb E[\|P_+By\|^2]=2\|P_+BP_+\|_F^2\) holds for every edge weight (2-design / frame). Evidence: `e1_gmin_gap_probe.json`.

6. **Sufficient numerical bound (certified, not proved).** At \(p=5,7\), \(\lambda_{\mathrm{cycle}}\le 8n/(n+4)\le n/2\) (the \(2\times\) spherical fourth-moment comparison). At \(p=3\), the same comparison fails (\(\lambda_{\mathrm{cycle}}=8>8n/(n+4)\approx5.71\)), consistent with the gap failure. Evidence: `e1_gmin_gap_probe.json`. **Not a proof:** the fourth-moment ratio vs the sphere reaches \(\approx2.8>2\) at \(p=3\) and \(\approx1.95,1.43\) at \(p=5,7\); a universal factor-\(2\) sphere bound is false.

7. **Residual (OPEN).** Prove part 3 for every prime \(p\ge5\) (e.g. \(\|T(x)\|_F^2\le nN\|x\|^2\) on \(\mathbf1^\perp\), or \(\lambda_2(W)\le\lambda_1(W)/2\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.59 (central symmetry, \(P\mathbf1=0\), rank of \(P\odot P\), two-moment; 2026-07-30).** Continue Prop 15.58.

1. **Central symmetry (proved).** \(\mathrm{Max}_{+}=-\mathrm{Max}_{+}\). *Proof.* \(Cy=py\) iff \(C(-y)=p(-y)\), and \(-y\in\{\pm1\}^n\). \(\square\)

2. **Frame is centered (proved).** \(\sum_{y\in\mathrm{Max}_{+}}y=0\) in \(\mathbb R^n\). *Proof.* Pair \(y\) with \(-y\) from part 1. \(\square\)
   Consequently \(Y^\top\mathbf1=0\), so the orthoprojector \(P=YY^\top/(2N)\) satisfies \(P\mathbf1=0\) (range of \(P\) lies in \(\mathbf1^\perp\subset\mathbb R^N\)). Certified \(p=3,5,7\): `e1_gmin_veronese.json`.

3. **Rank formula (certified \(p=3,5\); formula holds \(p=7\)).** \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\). Matches \(1+\mathrm{rank}_{\mathrm{cycle}}\) with \(\mathrm{rank}_{\mathrm{cycle}}=\binom{d-1}{2}-1\). Evidence: `e1_gmin_veronese.json`. *(Uniform proof of the rank formula for all primes \(p\ge3\) still open; not required if the Veronese bound is proved by other means.)*

4. **Two-moment on \(W_{ab}=(u_a\cdot u_b)^2\) (partial).** Writing \(\lambda_1(W)=N/d\) and \(\sum_{j\ge2}\lambda_j(W)=N(1-1/d)\) (from \(\mathrm{Tr}(W)=N\)), the one-large-rest-equal estimate from \((\mathrm{Tr}(W^2),\mathrm{Tr}(W))\) **forces** \(\lambda_2(W)\le N/(2d)\) at \(p=7\), but **not** at \(p=3\) (correct: gap fails) or \(p=5\) (worst-case \(12.78>10\); actual \(\lambda_2\approx5.21\) with multiplicity \(d=13\)). Evidence: `e1_gmin_veronese.json`. Multiplicity-\(d\) two-moment would force the gap at \(p=5\); multiplicity not yet proved for general \(p\).

5. **Residual (OPEN).** Same as Prop 15.58.7: prove \(\|T(x)\|_F^2\le nN\|x\|^2\) for all \(x\perp\mathbf1\) and all primes \(p\ge5\). Then \(\lambda_{\max}(G)=n/2\) simple, Prop 15.55 closes bi-tight / Type I, deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.60 (antipodal reduction to projective ENTF; 2×sphere sufficient for \(p\ge5\); 2026-07-30).** Continue Prop 15.58–15.59. Write \(m=N/2\) and fix a set \(\mathcal L\) of representatives of \(\mathrm{Max}_{+}/\{\pm1\}\) (one vector from each antipodal pair). Let \(u_y=y/\sqrt n\in V_+\) and \(W^{(\mathrm{proj})}_{yy'}=(u_y\cdot u_{y'})^2\) on \(\mathcal L\).

1. **Antipodal reduction for \(T\) (proved).** If \(x\in\mathbb R^N\) and \(s_a=\tfrac12(x_a+x_{-a})\) (indices via \(y\mapsto -y\)), then
   \[
   T(x)=T(s)=\sum_{a}s_a\,y_ay_a^\top.
   \]
   *Proof.* \(y y^\top=(-y)(-y)^\top\), so the antisymmetric part \(x_a-x_{-a}\) cancels. \(\square\)
   Consequently \(\|x\|^2=\|s\|^2+\|x-s\|^2\ge\|s\|^2\) and \(\|T(x)\|_F=\|T(s)\|_F\), so the Veronese inequality on all of \(\mathbf1^\perp\) reduces to antipode-symmetric \(x\) (equivalently, to functions on \(\mathcal L\)).

2. **Projective ENTF (proved).** The \(m\) unit vectors \(\{u_y:y\in\mathcal L\}\) form an equal-norm tight frame in \(V_+\cong\mathbb R^d\):
   \[
   \sum_{y\in\mathcal L}u_yu_y^\top=\frac m d\,I_d,\qquad \|u_y\|=1.
   \]
   *Proof.* \(\sum_{\mathrm{Max}_{+}}uu^\top=(N/d)I_d\) and antipodal pairs contribute identical \(uu^\top\), so the sum over \(\mathcal L\) is half. \(\square\)
   Moreover \(W^{(\mathrm{proj})}\mathbf1=(m/d)\mathbf1\) and \(\lambda_1(W^{(\mathrm{proj})})=m/d\).

3. **Eigenvalue doubling (proved).** On the antipode-symmetric subspace of \(\mathbb R^N\),
   \[
   \lambda_2(W)=2\,\lambda_2(W^{(\mathrm{proj})}),
   \]
   where \(W_{ab}=(u_a\cdot u_b)^2\) is the full Max+ Schur square. *Proof.* Each antipodal \(2\times2\) block of \(W\) is the all-ones matrix of order 2 (since \((u\cdot(\pm u'))^2=(u\cdot u')^2\)), so the Rayleigh quotient of an antipode-symmetric vector with values \(c\) on \(\mathcal L\) equals \(2\,c^\top W^{(\mathrm{proj})}c/\|c\|_{2,\mathrm{sym}}^2\). \(\square\)
   Therefore
   \[
   \lambda_2(P\odot P)\le\frac\alpha2
   \quad\Longleftrightarrow\quad
   \lambda_2(W^{(\mathrm{proj})})\le\frac m{2d}.
   \]

4. **Spherical comparison algebra (proved).** For the continuous sphere (or any spherical 4-design) in \(\mathbb R^d\), the maximal fourth-moment Rayleigh on trace-free matrices is \(2m/(d(d+2))\). The bound
   \[
   \frac{4m}{d(d+2)}\le\frac m{2d}
   \]
   rearranges to \(d\ge6\). Since \(d=(p^2+1)/2\ge13\) for primes \(p\ge5\), **any proof that the projective Max+ fourth-moment Rayleigh is at most twice the spherical value yields the spectral gap for all primes \(p\ge5\)**.

5. **Certified 2×sphere (p=5,7).** On projective Max+: maximising \(\sum_y(u_y^\top A u_y)^2\) over \(\mathrm{Tr}\,A=0\), \(\|A\|_F=1\) gives ratio-to-sphere \(\approx1.95\) at \(p=5\) and \(\approx1.43\) at \(p=7\) (both \(<2\)), and \(\approx2.80>2\) at \(p=3\). Gap holds \(p=5,7\), fails \(p=3\). Evidence: `e1_gmin_projective.json`.

6. **Residual (OPEN).** Prove
   \[
   \max_{\mathrm{Tr}\,A=0,\ \|A\|_F=1}\sum_{y\in\mathcal L}(u_y^\top A u_y)^2
   \;\le\;
   \frac{4m}{d(d+2)}
   \]
   for every prime \(p\ge5\) (or any upper bound \(\le m/(2d)\)). Then Prop 15.55 closes bi-tight / Type I for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.61 (16N bound \(\Rightarrow\lambda_{\mathrm{cycle}}\le8\Rightarrow\) gap for \(p\ge5\); 2026-07-30).** Continue Prop 15.56–15.60. Let \(B\) range over real \(n\times n\) matrices with \(B=P_+BP_+\), \(\mathrm{Tr}\,B=0\), and write
\[
Q(B):=\sum_{y\in\mathrm{Max}_{+}}(y^\top By)^2.
\]

1. **Frame identity (proved).** For every such \(B\),
   \[
   \sum_{y\in\mathrm{Max}_{+}}\|By\|^2=2N\|B\|_F^2.
   \]
   *Proof.* \(\sum_y y^\top B^2 y=\mathrm{Tr}\bigl(B^2\sum_y yy^\top\bigr)=\mathrm{Tr}(B^2\cdot 2N P_+)=2N\|B\|_F^2\). \(\square\)

1b. **Triangle form on \(V_+\) (proved).** For every \(B=P_+BP_+\),
   \[
   \mathrm{Tr}(CB^2)=p\|B\|_F^2.
   \]
   *Proof.* \(C=p(P_+-P_-)\) and \(B^2=P_+B^2P_+\), so \(\mathrm{Tr}(CB^2)=p\mathrm{Tr}(B^2)\). \(\square\)
   Consequently every zero-diagonal \(B\) on \(V_+\) with \(\|B\|_F^2=2\) (i.e. unit edge-weight \(v\)) saturates the triangle bound \(\mathrm{ft}=2p\) and has \(\mathbb E[\|By\|^2]=4\). The spectral gap residual is therefore purely the alignment \(\mathbb E[(y^\top By)^2]\) among these maximisers of \(\mathrm{ft}\).

2. **Equivalence with \(\lambda_2(W)\) (proved).** Identifying \(B=QAQ^\top\) with \(Q\) an ONB of \(V_+\) and \(A\in\mathrm{Sym}(\mathbb R^d)\), \(\mathrm{Tr}\,A=0\),
   \[
   Q(B)=n^2\sum_{a=1}^N(u_a^\top A u_a)^2,
   \]
   so
   \[
   \max_{\|B\|_F=1}Q(B)=n^2\cdot\lambda_2(W)=4d^2\cdot\lambda_2(W).
   \]
   Consequently
   \[
   \max_{\|B\|_F=1}Q(B)\le 16N
   \quad\Longleftrightarrow\quad
   \lambda_2(W)\le\frac{4N}{d^2}.
   \]

3. **Equivalence with \(\lambda_{\mathrm{cycle}}\le8\) (proved).** At a cycle maximiser of \(G\) one has \(\mathrm{ft}=2p\), \(\mathbb E[\|By\|^2]=4\), \(By\in V_+\) for all Max+ \(y\), and \(\|B\|_F^2=2\) (Prop 15.57 certs). Scaling \(B'=B/\sqrt2\) (\(\|B'\|_F=1\)) gives
   \[
   \lambda_{\mathrm{cycle}}=\frac1{4N}Q(B)=\frac1{2N}Q(B').
   \]
   Hence \(\max Q\le16N\) yields \(\lambda_{\mathrm{cycle}}\le8\). Conversely, the \(\Phi\)-maximiser of \(Q\) realises \(\lambda_{\mathrm{cycle}}=Q(B_\star)/(2N)\) under the same scaling, so the two maxima match. \(\square\)

4. **Algebraic gap upgrade (proved).** If \(\lambda_{\mathrm{cycle}}\le8\), then for every prime \(p\ge5\) one has \(n/2=d\ge13>8\), so \(\lambda_{\max}(G)=\max(n/2,\lambda_{\mathrm{cycle}})=n/2\) is simple, and Prop 15.55 blocks all Max+-tight size-\(2p\) covers (bi-tight / Type I empty). Equivalently: \(\lambda_2(W)\le4N/d^2\) and \(d\ge8\) imply \(\lambda_2(W)\le N/(2d)\) because \(4/d\le1/2\). \(\square\)

5. **Certified 16N bound.** At \(p=3\): \(Q_{\max}=16N=192\) and \(\lambda_{\mathrm{cycle}}=8\) (**equality**). At \(p=5\): \(Q_{\max}/(16N)=11/13<1\), \(\lambda_{\mathrm{cycle}}=88/13<8\). At \(p=7\): \(Q_{\max}/(16N)\approx0.660<1\), \(\lambda_{\mathrm{cycle}}\approx5.281<8\). Spectral gap holds \(p=5,7\), fails \(p=3\). Evidence: `e1_gmin_16n.json`.

6. **Residual (OPEN).** Prove \(Q(B)\le16N\|B\|_F^2\) for all \(B=P_+BP_+\) with \(\mathrm{Tr}\,B=0\) and all primes \(p\ge5\) (equality at \(p=3\) is the base case). Then bi-tight / Type I closes for every such \(p\) via Prop 15.55. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* The 2×-sphere bound of Prop 15.60.4 is slightly sharper (\(16\,dN/(d+2)<16N\)) but has the same open core; the 16N form isolates the clean numerical threshold \(\lambda_{\mathrm{cycle}}\le8\).

**Proposition 15.62 (typeA+wedge identity; \(Q=6N\|B\|_F^2+Q_4\); 2026-07-30).** Continue Prop 15.61. Let \(B\) be real symmetric \(n\times n\) with \(B=P_+BP_+\) and zero ambient diagonal (hence \(\mathrm{Tr}\,B=0\) automatically). Write \(Q(B)=\sum_{y\in\mathrm{Max}_{+}}(y^\top By)^2\). Index unordered edges \(e=\{i,j\}\) by \(Be_e:=2B_{ij}\) and \(f_e(y):=y_iy_j\), so \(y^\top By=\sum_e Be_e f_e(y)\). Let \(\mathrm{Gu}=F^\top F\) with \(F\) the \(N\times E\) matrix of features \(f_e\) on Max+, and split pairs of distinct edges into **wedge** (share a vertex) and **disjoint**.

1. **Type A (same-edge; proved).**
   \[
   \sum_e Be_e^2\sum_y f_e(y)^2=N\|Be\|_2^2=2N\|B\|_F^2,
   \]
   since \(f_e^2\equiv1\) and \(\|Be\|_2^2=4\sum_{i<j}B_{ij}^2=2\|B\|_F^2\). \(\square\)

2. **Wedge (share-one-vertex; proved).** For edges \(e=\{r,j\}\), \(e'=\{r,l\}\) sharing \(r\),
   \[
   \sum_y f_e f_{e'}=\sum_y y_jy_l=N\Sigma_{jl},\qquad \Sigma:=\mathbb E[yy^\top]=2P_+.
   \]
   Summing over all ordered wedge pairs and using \(Be_e=2B_{rj}\),
   \[
   Be^\top(\mathrm{Gu}\odot\mathbf1_{\mathrm{wedge}})Be
   =\sum_r\sum_{j\neq l,\,j,l\neq r}(2B_{rj})(2B_{rl})\cdot N\Sigma_{jl}
   =4N\sum_r\Bigl((Be_r)^\top\Sigma(Be_r)-\sum_j B_{rj}^2\Sigma_{jj}\Bigr).
   \]
   On \(V_+\): \(\Sigma v=2v\) for \(v=Be_r=B e_r\in V_+\), and \(\Sigma_{jj}=1\), so the parenthesis equals \(\|Be_r\|_2^2\). Hence the wedge contribution is
   \[
   4N\sum_r\|Be_r\|_2^2=4N\|B\|_F^2.
   \]
   \(\square\)

3. **Identity (proved).**
   \[
   Q(B)=6N\|B\|_F^2+Q_4(B),\qquad Q_4(B):=Be^\top(\mathrm{Gu}\odot\mathbf1_{\mathrm{disj}})Be.
   \]
   Consequently \(Q(B)\le16N\|B\|_F^2\) if and only if \(Q_4(B)\le10N\|B\|_F^2\). For the cycle-normalisation \(\|B\|_F^2=2\),
   \[
   \lambda_{\mathrm{cycle}}=\frac{Q(B)}{4N}=3+\frac{Q_4(B)}{4N},
   \]
   so \(\lambda_{\mathrm{cycle}}\le8\) \(\Leftrightarrow\) \(Q_4\le20N\) at that scale (equivalently \(Q_4\le10N\|B\|_F^2\)). \(\square\)

4. **Certified.** Multi-seed sampling of the full zero-diag \(\cap V_+\) space (nullspace of ambient diagonal on \(\mathrm{Sym}(V_+)\)) at \(p=3,5,7\): typeA+wedge identity holds to machine precision in every trial; at \(p=3\) one has \(Q\equiv16N\) on the whole space (equality case); at \(p=5,7\) the maximiser and all random trials satisfy \(Q_4\le10N\|B\|_F^2\). Unrestricted \(\|\mathrm{Gu}_{\mathrm{disj}}\|_{\mathrm{op}}\) exceeds the \(5N\) Rayleigh threshold for unit edge vectors, so a crude operator-norm bound fails — the residual is Rayleigh of \(\mathrm{Gu}_{\mathrm{disj}}\) on the **image** of zero-diag \(\cap V_+\to\mathbb R^E\) only. Evidence: `e1_gmin_typeA_wedge.json`.

5. **Residual (OPEN).** Prove \(Q_4(B)\le10N\|B\|_F^2\) for all zero-diag \(B=P_+BP_+\) and all primes \(p\ge5\). Then Prop 15.61 closes bi-tight for every such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.63 (hypothesis H \(\Rightarrow\) 16N; spectrum certs; 2026-07-30).** Continue Prop 15.62. For unit-Frobenius zero-diag \(B=P_+BP_+\) write
\[
\mathrm{ray}(B)\,:=\,\frac{Q_4(B)}{2N}
\,=\,\frac{Be^\top G_{\mathrm{disj}}Be}{\|Be\|_2^2}
\]
(the Rayleigh quotient of \(G_{\mathrm{disj}}=\mathbb E[ff^\top]\) on vertex-disjoint pairs, evaluated on the image of the edge map). Set
\[
H(p)\,:=\,\frac{(p+2)^2}{d},\qquad d=\frac{p^2+1}{2}.
\]

1. **Algebra (proved).** \(H(p)\le5\) for every prime \(p\ge3\), with equality if and only if \(p=3\).
   *Proof.* \(H(p)=2(p+2)^2/(p^2+1)\), so \(H\le5\) \(\Leftrightarrow\) \(2(p+2)^2\le5(p^2+1)\) \(\Leftrightarrow\) \(3p^2-8p-3\ge0\). The positive root of \(3p^2-8p-3=0\) is \(p=3\); the quadratic is nonnegative for all primes \(p\ge3\), and vanishes only at \(p=3\). \(\square\)

2. **H \(\Rightarrow\) 16N (proved).** If \(\mathrm{ray}(B)\le H(p)\) for every unit zero-diag \(B=P_+BP_+\), then \(Q_4(B)\le2N\cdot H(p)\le10N\), hence \(Q(B)\le16N\|B\|_F^2\) (Prop 15.62.3), hence \(\lambda_{\mathrm{cycle}}\le8\), and for \(p\ge5\) bi-tight is empty (Prop 15.61.4). \(\square\)

3. **Spectrum of \(Q_4\) on zero-diag \(\cap V_+\) (certified).**
   - \(p=3\): \(\dim=5\), single eigenvalue \(Q_4\equiv10N\) (mult 5); \(\mathrm{ray}\equiv5=H(3)\).
   - \(p=5\): \(\dim=65\), eigenvalues \(Q_4\in N\cdot\{\tfrac{98}{13},\tfrac{66}{13},\tfrac{2}{13}\}\) with multiplicities \(\{13,26,26\}=\{d,2d,2d\}\); top \(\mathrm{ray}=\tfrac{49}{13}=H(5)\).
   - \(p=7\): power multistart gives \(\mathrm{ray}\approx2.281<H(7)=3.24=\tfrac{81}{25}\).
   Evidence: `e1_gmin_q4_spectrum.json`, `e1_gmin_q4_ub.json`, `e1_gmin_q4_bound.json`.

4. **Hypothesis H (certified \(p=3,5,7\); OPEN in general).** \(\mathrm{ray}(B)\le H(p)\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge3\), with equality at \(p=3\) (whole space) and at the maximiser for \(p=5\). Combined with part 1–2 this yields the 16N bound. **Uniform proof of H for all primes \(p\ge5\) remains OPEN.**

5. **Dead ends (do not reopen).** Pointwise \(|y^\top By|\) or \(\cos^2\) bounds (max \(\cos^2\) exceeds the average threshold); unrestricted \(\|G_{\mathrm{disj}}\|_{\mathrm{op}}\) (exceeds \(H(p)\)); linear span of \(\{\|B\|_F^2,\mathrm{tr}(B^4),\sum B_{ij}^4,\sum_i\|row_i\|^4\}\) (fails at \(p=5,7\)); general CS-ENTF theory (random frames violate 16N).

6. **Residual (OPEN).** Prove hypothesis H for all primes \(p\ge5\) (or any upper bound \(\mathrm{ray}\le5\)). Then bi-tight closes via Prop 15.61–15.62. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.64 (dual fourth-moment form of H; residual reformulation; 2026-07-30).** Continue Prop 15.63. Write \(s_y:=V_+^\top y\in\mathbb R^d\) (so \(\|s_y\|^2=n\), \(\mathbb E[ss^\top]=2I_d\)), and for \(A\in\mathrm{Sym}(\mathbb R^d)\) set \(\Phi(A):=\mathbb E[(s^\top As)\,ss^\top]\). Let \(\mathcal Z\) be the subspace of traceless \(A\) with zero ambient diagonal on \(B=V_+AV_+^\top\) (i.e. \(r_i^\top Ar_i=0\) for all rows \(r_i\) of \(V_+\)).

1. **Duality (proved).** For every zero-diag \(B=V_+AV_+^\top\),
   \[
   Q(B)=\sum_y(s_y^\top A s_y)^2=N\langle\Phi(A),A\rangle_F.
   \]
   Consequently
   \[
   \max_{\|A\|_F=1,\,A\in\mathcal Z}Q(B)=N\cdot\lambda_{\max}(\Phi|_{\mathcal Z}),
   \]
   and the maximiser satisfies the eigenmatrix equation \(\Phi(A)=\lambda A\) on \(\mathcal Z\). Equivalently, writing \(D=YY^\top=SS^\top\), one has on \(\mathbf1^\perp\)
   \[
   K=\frac{D\odot D}{2N},\qquad\lambda_{\mathrm{cycle}}=\frac{\lambda_2(D\odot D)}{2N}=\frac12\max_{\|x\|=1,\,x\perp\mathbf1}\Bigl\|\sum_y x_y s_ys_y^\top\Bigr\|_F^2.
   \]
   Evidence: `e1_gmin_H_proof.json` (gen.eig on \(\mathcal Z\)).

2. **Wick baseline (proved).** If \(s\) were Gaussian with \(\mathrm{Cov}=2I\), then \(\Phi=8\,\mathrm{Id}\) on \(\mathrm{Tr}\,A=0\) and \(Q=8N\|B\|_F^2\). For actual Max+,
   \[
   \frac QN=8+\mathrm{residual}(A),\qquad\mathrm{residual}(A)=\langle\kappa,A\otimes A\rangle,
   \]
   with cumulant \(\kappa=\Phi-8\,\mathrm{Id}\).

3. **H \(\Leftrightarrow\) residual bound (proved).** Hypothesis H is equivalent to
   \[
   \max_{\|A\|_F=1,\,A\in\mathcal Z}\mathrm{residual}(A)\;\le\;\frac{(p+1)(p+7)}{d},
   \]
   because \(6+2H(p)-8=2(H(p)-1)=(p+1)(p+7)/d\). Equality holds at \(p=3\) (whole \(\mathcal Z\)) and at the maximiser for \(p=5\). \(\square\)

4. **Certified.** Exact gen.eig of \(\Phi|_{\mathcal Z}\) at \(p=3,5,7\): H holds; residual ratios to budget \(1,1,0.572\). Spectrum of \(K\) at \(p=5\): \(\{13,\tfrac{88}{13},\tfrac{72}{13},\tfrac{40}{13}\}\) with mults \(\{1,d,2d,2d\}\). Evidence: `e1_gmin_H_proof.json`, `e1_gmin_q4_spectrum.json`.

5. **Residual (OPEN).** Prove \(\lambda_{\max}(\Phi|_{\mathcal Z})\le6+2H(p)\) for all primes \(p\ge5\) (equivalently residual \(\le(p+1)(p+7)/d\), or \(\mathrm{ray}\le H(p)\), or \(\lambda_{\mathrm{cycle}}\le3+H(p)\)). Then 16N and bi-tight close. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark (dead ends for H).* Gershgorin on \(K\); unrestricted \(\|G_{\mathrm{disj}}\|_{\mathrm{op}}\); mult\(\ge d\) of \(\lambda_2(W)\) with only \(\mathrm{Tr}(W^2)\) (too weak for 16N at \(p=5\)); pointwise \(\cos^2\); 2×sphere (fails at \(p=3\), holds numerically \(p\ge5\)); distance-homogeneous scheme formulae (Max+ not distance-homogeneous at \(p=7\)); two-moment feasibility of \(\lambda_2(P\odot P)\) (always allows values \(>4/N\)).

**Proposition 15.65 (κ spectrum on \(\mathcal Z\); clean form \(\lambda_2(P\odot P)\le4/N\); boolean essential; 2026-07-30).** Continue Prop 15.64. Write \(P=YY^\top/(2N)\) (equal-diagonal orthoprojector of rank \(d\), diagonal \(\alpha=d/N\), and \(P\mathbf1=0\)).

1. **Clean equivalences (proved).**
   \[
   P\odot P=\alpha^2 W,\qquad W_{ab}=(u_a\cdot u_b)^2,\quad u_a=y_a/\sqrt n.
   \]
   Hence \(\lambda_2(P\odot P)=\alpha^2\lambda_2(W)\). Combined with Prop 15.56–15.61:
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_2(P\odot P)\le4/N\Longleftrightarrow\lambda_2(W)\le4N/d^2,\\
   H&\Longleftrightarrow\lambda_2(P\odot P)\le\frac{3+H(p)}{2N},\\
   \mathrm{gap}&\Longleftrightarrow\lambda_2(P\odot P)\le\alpha/2=d/(2N).
   \end{aligned}
   \]
   Moreover \(\lambda_{\max}(P\odot P)=\alpha\) for every equal-diagonal orthoprojector (Prop 15.58.2). \(\square\)

2. **Cumulant spectrum on \(\mathcal Z\) (certified).** The residual form \(\mathrm{residual}(A)=\langle\kappa,A\otimes A\rangle\) on \(\mathcal Z\) has eigenvalues
   - \(p=3\): \(\{8\}\) (mult \(\dim\mathcal Z=5\));
   - \(p=5\): \(\{-\tfrac{24}{13},\tfrac{40}{13},\tfrac{72}{13}\}\) with multiplicities matching \(\{26,26,13\}\) (top \(=\mathrm{budget}=(p+1)(p+7)/d\));
   - \(p=7\): top \(=\tfrac{1048}{409}<\tfrac{112}{25}=\mathrm{budget}\).
   Evidence: `e1_gmin_cumulant.json`.

3. **Boolean/conference structure is essential (certified counterexample).** Random equal-diagonal rank-\(d\) orthoprojectors need **not** satisfy \(\lambda_2(P\odot P)\le4/N\): at \((N,d)=(50,10)\), \(20/20\) Haar-row-equalized samples violated the bound. Thus no proof of 16N can use only the equal-diagonal projector axioms; the Max+/boolean/conference structure is load-bearing. Evidence: `e1_gmin_cumulant.json`.

4. **Residual (OPEN).** Prove \(\lambda_2(P\odot P)\le4/N\) (or the sharper H form) for the Max+ projector of every prime \(p\ge5\). Equivalent targets: \(\lambda_{\max}(\kappa|_{\mathcal Z})\le(p+1)(p+7)/d\), or \(\mathrm{ray}\le H(p)\). Then bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.66 (zero-diag freeness of \(\Phi\)-maximiser; pairing residual criterion for \(g_{\min}\ge L(p)\); 2026-07-30).** Continue Prop 15.53–15.54 and 15.64–15.65.

1. **Zero-diag freeness (proved).** Let \(A\in\mathrm{Sym}_0(\mathbb R^d)\) maximise \(\langle\Phi A,A\rangle\) over \(\|A\|_F=1\). Then the ambient matrix \(B=V_+AV_+^\top\) has zero diagonal, so the maximiser already lies in \(\mathcal Z\). Consequently
   \[
   \lambda_{\max}(\Phi|_{\mathcal Z})=\lambda_{\max}(\Phi|_{\mathrm{Sym}_0}).
   \]
   *Proof.* The Lagrangian for the constraints \(r_i^\top Ar_i=0\) and \(\mathrm{Tr}\,A=0\) yields
   \(\Phi(A)=\lambda A+\sum_i\mu_i r_ir_i^\top+\nu I\).
   Taking the Frobenius product with \(r_kr_k^\top\) and using \(r_k^\top s_y=y_k\) (hence \(r_k^\top\Phi(A)r_k=\mathbb E[s^\top As]=0\)) forces \(\mu_k/4+\nu/2=0\) for every \(k\), so all \(\mu_k\) are equal. Then \(\sum_i\mu_i r_ir_i^\top=-2\nu I\), and \(\Phi(A)=\lambda A-\nu I\). Taking traces and using \(\mathrm{Tr}\,\Phi(A)=n\mathbb E[s^\top As]=0\) gives \(\nu=0\), hence \(\mu_k=0\) and \(\Phi(A)=\lambda A\). The unconstrained critical point on \(\mathrm{Sym}_0\) is therefore admissible for \(\mathcal Z\). Certified: ambient \(\mathrm{diag}(B)\) of the power-iteration maximiser is \(O(10^{-16})\) at \(p=3,5,7\). Evidence: `e1_gmin_m4_residual.json`. \(\square\)

2. **Pairing residual criterion (proved algebra).** Write \(L(p)=-(p-2)/(2p^2)\) and \(T(p)=-(p-2)/(p(2p-1))\). For every prime \(p\ge5\) one has \(L(p)>T(p)\). On any 4-set with \(|\kappa|=1\), the Wick value is \(m_4^{\mathrm{Wick}}=\kappa/p^2\). If
   \[
   \bigl|m_4-\kappa/p^2\bigr|\;\le\;\frac{p-4}{2p^2}
   \]
   for every such 4-set, then
   \[
   |m_4|\;\le\;\frac1{p^2}+\frac{p-4}{2p^2}=\frac{p-2}{2p^2}=-L(p),
   \]
   hence \(g_{\min}\ge L(p)>T(p)\), and Prop 15.47 closes bi-tight. \(\square\)

3. **Certified m4 tables on \(|\kappa|=1\).**
   | \(p\) | \(g_{\min}\) | \(L(p)\) | \(\max|m_4|\) | \(\max|m_4-\kappa/p^2|\) | resid crit. |
   |------|------------|---------|--------------|--------------------------|-------------|
   | 3 | \(-1/3\) | \(-1/18\) | \(1/3\) | \(2/9\) | n/a (\(p<5\)) |
   | 5 | \(-3/65\) | \(-3/50\) | \(3/65\) | \(0.0554>(p-4)/(2p^2)=1/50\) | **fails** (too crude) |
   | 7 | \(-109/2863\) | \(-5/98\) | \(109/2863\) | \(0.0177\le3/98\) | **holds** |

   In particular \(g_{\min}\ge L(p)\) holds at \(p=5,7\), but the triangle residual criterion fails at \(p=5\) (large residuals occur on classes with small \(|m_4|\)). Evidence: `e1_gmin_m4_residual.json`.

4. **Residual (OPEN).** Prove \(g_{\min}\ge L(p)\) for all primes \(p\ge5\), e.g. by proving \(|m_4|\le(p-2)/(2p^2)\) on every \(|\kappa|=1\) class (directly, not via the failed \(p=5\) triangle), **or** prove \(\lambda_2(P\odot P)\le4/N\). Either closes bi-tight. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* K4-edge Gram PSD only forces \(\alpha\le(p-1)/(2p)\)-scale bounds (e.g. \(0.6\) at \(p=5\)), far above \(L(p)\). General equal-diag projectors violate \(4/N\) (Prop 15.65.3).

**Proposition 15.67 (master m4 identity; same-sign Ext criterion; full \(|\kappa|=1\) census \(p=5,7\); 2026-07-30).** Continue Prop 15.66. Write \(m_4(S)=\mathbb E[y_a y_b y_c y_d]\) for \(S=\{a,b,c,d\}\) and Max+ average, \(\kappa(S)\) the sum of the three perfect-matching edge-products of \(C\), and
\[
(Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r}).
\]

1. **Combinatorial identity \(\sigma_{\mathrm{sum}}=4\kappa\) (proved).** For every \(\pm1\) edge-labeling of \(K_4\), writing \(\sigma_v\) for the three-term pairing at vertex \(v\), one has \(\sum_{v\in S}\sigma_v=4\kappa(S)\). *Proof.* Exhaustion of the \(2^6=64\) labelings (all pass). \(\square\)

2. **Master identity (proved for Max+).** Averaging the eigenvector identity \(Cy=py\) against three coordinates of a Max+ vector and using part 1 yields, for every 4-set \(S\),
   \[
   m_4(S)=\frac{\kappa(S)}{p^2}+\frac{\mathrm{Ext}(S)}{4p},\qquad \mathrm{Ext}:=Tm_4.
   \]
   Equivalently \(\mathrm{Ext}(S)=4pm_4(S)-4\kappa(S)/p\). Certified float residual \(<10^{-16}\) on all \(|\kappa|=1\) 4-sets at \(p=5\) and on a full multi-worker census at \(p=7\). Evidence: `e1_gmin_m4_bound.json`, `e1_gmin_m4_proof.json`. \(\square\)

3. **Same-sign Ext criterion (proved algebra).** On any 4-set with \(|\kappa|=1\), if \(\mathrm{sign}(\mathrm{Ext})=\mathrm{sign}(\kappa)\) and
   \[
   |\mathrm{Ext}|\;\le\;\frac{2(p-4)}{p},
   \]
   then
   \[
   |m_4|=\frac1{p^2}+\frac{|\mathrm{Ext}|}{4p}\;\le\;\frac1{p^2}+\frac{p-4}{2p^2}=\frac{p-2}{2p^2}=-L(p).
   \]
   (Opposite-sign Ext only decreases \(|m_4|\) relative to the Wick value.) Thus same-sign Ext control on every \(|\kappa|=1\) class implies \(g_{\min}\ge L(p)\). \(\square\)

4. **Full multi-worker census (certified \(p=5,7\); \(W=86\)).** Over every 4-set with \(|\kappa|=1\):
   | \(p\) | \(\#\{|\kappa|=1\}\) | \(\max|m_4|\) | \(L_{\mathrm{abs}}\) | same-sign \(\max|\mathrm{Ext}|\) | thr \(2(p-4)/p\) |
   |------|---------------------|---------------|----------------------|----------------------------------|------------------|
   | 5 | 11700 | \(3/65\) | \(3/50\) | \(0.123<0.4\) | yes |
   | 7 | 176400 | \(109/2863\) | \(5/98\) | \(0.495<6/7\) | yes |

   Hence \(g_{\min}\ge L(p)\) at \(p=5,7\) by either the direct bound or the same-sign Ext criterion. Evidence: `e1_gmin_m4_proof.json` (F17: `src/workers.py` + `ProcessPoolExecutor(W=nproc-2)`). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that either
   - \(|m_4|\le(p-2)/(2p^2)\) on all \(|\kappa|=1\) 4-sets, or
   - same-sign \(|\mathrm{Ext}|\le2(p-4)/p\) on those sets, or
   - \(\lambda_2(P\odot P)\le4/N\),

   using Max+/boolean/conference structure (Prop 15.65.3: bare equal-diag projectors are insufficient). Then Prop 15.47/15.55 closes bi-tight for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark (failed shortcuts).* Absolute triangle \(|m_4-\kappa/p^2|\le(p-4)/(2p^2)\) fails at \(p=5\) (large residuals on opposite-sign classes with small \(|m_4|\)). Affine+\(\mathrm{PGL}\)+Frob orbit of the halfspace vector is incomplete (\(60/260\) at \(p=5\); F18). Wick fourth-moment pin only reaches the \(T(p)\) endpoint, not strict \(g_{\min}>T(p)\).

**Proposition 15.68 (\(T\kappa\) calculus; residual source on \(|\kappa|=3\); resolvent reduction of \(L\); 2026-07-30).** Continue Prop 15.67. Let \(C\) be any real symmetric conference matrix of order \(n=p^2+1\) (\(C^\top=C\), zero diagonal, off-diagonal \(\pm1\), \(C^2=p^2I\)), and write
\[
(Tf)(S)=\sum_{v\in S,\,r\notin S}C_{vr}\,f(S_{v\to r}),\qquad
\kappa(S)=\sum_{\text{three pairings}}C_eC_{e'}.
\]

1. **Reduction of \(T\kappa\) to \(K_4\) (proved).** For every 4-set \(S\), the external sum defining \((T\kappa)(S)\) collapses via \((C^2)_{vx}=0\) (\(v\neq x\)) to a function of the six edge signs of \(S\) alone:
   \[
   (T\kappa)(S)=-6\sum_{v\in S}\prod_{u\in S\setminus\{v\}}C_{vu}
   \]
   (sum of star-triples). In particular \((T\kappa)(S)\) is independent of \(n\) and of the ambient graph outside \(S\). \(\square\)

2. **Vanishing on \(|\kappa|=1\) (proved).** Exhausting all \(2^6=64\) edge-labelings of \(K_4\): if \(|\kappa(S)|=1\) then \((T\kappa)(S)=0\); if \(|\kappa(S)|=3\) then \((T\kappa)(S)\in\{\pm24\}\). Consequently, for **every** conference matrix and every 4-set with \(|\kappa|=1\),
   \[
   (T\kappa)(S)=0.
   \]
   Evidence: `e1_gmin_m4_tkappa.json` (symbolic \(C^2\) reduction + 64-labeling check). \(\square\)

3. **Residual equation and source support (proved).** For Max+ moments write \(\rho:=m_4-\kappa/p^2\). The master identity (Prop 15.67) rearranges to
   \[
   (4p\,I-T)\rho=\frac{T\kappa}{p^2}.
   \]
   By part 2 the right-hand side **vanishes on every \(|\kappa|=1\) 4-set** and is bounded by \(24/p^2\) on \(|\kappa|=3\). Thus all of the same-sign residual on the dangerous classes is the image, under the resolvent \((4pI-T)^{-1}\), of a pure \(|\kappa|=3\) source. \(\square\)

4. **Paley extension degrees (proved formula; certified \(p=3,5,7\)).** For the Paley conference graph and every 4-set with \(|\kappa|=1\), among the \(4(n-4)=4(p^2-3)\) ordered extensions \((v,r)\), exactly
   \[
   d_3=p^2-5,\qquad d_1=3p^2-7
   \]
   land in \(|\kappa|=3\) and \(|\kappa|=1\) respectively (both constant on the \(|\kappa|=1\) stratum). Certified by full census at \(p=3,5,7\). Evidence: `e1_gmin_m4_tkappa.json`. \(\square\)

5. **Resolvent reduction of the \(L\)-bound (proved algebra).** On any \(|\kappa|=1\) 4-set, if \(\mathrm{sign}(\rho)=\mathrm{sign}(\kappa)\) then \(|m_4|=1/p^2+|\rho|\). The bound \(|m_4|\le L_{\mathrm{abs}}=(p-2)/(2p^2)\) is therefore equivalent to
   \[
   |\rho|\;\le\;\frac{p-4}{2p^2}.
   \]
   Writing \(\rho=(4pI-T)^{-1}(T\kappa/p^2)\) and using \(|T\kappa/p^2|\le24/p^2\) on the source, it suffices to prove that the operator gain from the \(|\kappa|=3\) stratum into the same-sign \(|\kappa|=1\) stratum is at most
   \[
   \frac{p-4}{48}.
   \]
   (At \(p=5\) the budget is \(1/48\); empirical gain \(\approx0.0064\ll1/48\).) Equivalently, any upper bound
   \[
   |m_4|\;\le\;\frac{p-2}{p(2p+3)}
   \]
   on \(|\kappa|=1\) closes \(L\) because \((p-2)/(p(2p+3))\le(p-2)/(2p^2)\) for all odd \(p\ge5\), with equality of the two sides only in the large-\(p\) limit sense (strict for finite \(p\)); at \(p=5\) the candidate is sharp (\(|m_4|_{\max}=3/65\)). Certified: candidate \(\ge\max|m_4|\) at \(p=5,7\). \(\square\)

6. **Residual (OPEN).** Prove the resolvent gain bound of part 5 (or the candidate \(|m_4|\le(p-2)/(p(2p+3))\), or \(\lambda_2(P\odot P)\le4/N\)) for every prime \(p\ge5\), using Max+/boolean structure on top of the conference calculus above. Then \(g_{\min}\ge L(p)\) and bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* Absolute \(\infty\)-bootstrap with \(d_1M_1+d_3M_3\) does **not** contract (\(4p-(3p^2-7)<0\) for \(p\ge5\)). Signed cancellation / spectral control of \(T\) on the conference association scheme is load-bearing.

**Proposition 15.69 (spectrum of \(T\): \(\lambda_{\max}=4p\); singular resolvent; min-norm particular solution; 2026-07-30).** Continue Prop 15.67–15.68. View \(T\) as a self-adjoint operator on \(\ell^2\) of unordered 4-sets (signed Johnson adjacency with weights \(C_{vr}\)).

1. **Self-adjointness (proved).** If \(S'=S_{v\to r}\) then \(S=S'_{r\to v}\) and \(C_{rv}=C_{vr}\), so \(T\) is symmetric on \(\mathbb R^{\binom{n}{4}}\). Hence \(\|T\|_2=\rho(T)\). \(\square\)

2. **Spectral edges (certified \(p=5,7\); \(p=3\) strict).** Sparse exact eigensolve (ARPACK/`eigsh` on the CSR of \(T\)):
   | \(p\) | \(\lambda_{\max}(T)\) | \(4p\) | \(\lambda_{\min}(T)\) | mult of \(4p\) (top block) |
   |------|----------------------|-------|----------------------|---------------------------|
   | 3 | \(\approx9.798<12\) | 12 | \(\approx-9.798\) | 0 |
   | 5 | \(20\) | 20 | \(-20\) | \(\ge38\) |
   | 7 | \(28\) | 28 | \(-28\) | \(\ge1\) (full mult open) |

   In particular for primes \(p\ge5\) on Paley, the numerical evidence is
   \[
   \lambda_{\max}(T)=4p=-\lambda_{\min}(T),
   \]
   so \(4pI-T\) is **singular**. Evidence: `e1_gmin_m4_Tspec.json`, `e1_gmin_m4_pseudo.json`. \(\square\)

3. **Compatibility (certified \(p=5,7\)).** The master identity \((4pI-T)m_4=4\kappa/p\) is solvable: writing \(E_{4p}=\ker(4pI-T)\), one has \(4\kappa/p\perp E_{4p}\) to float precision \(<10^{-13}\). (Equivalently \(T\kappa/p^2\perp E_{4p}\) for the residual form.) \(\square\)

4. **General solution and min-norm particular solution (proved algebra + cert).** Every solution of the master linear equation is
   \[
   m=m_\star+h,\qquad h\in E_{4p},
   \]
   where \(m_\star=(4pI-T)^{+}(4\kappa/p)\) is the Moore–Penrose / min-norm solution. The Max+ moment vector is one particular solution: \(m_4=m_\star+h_\star\) for a unique \(h_\star\in E_{4p}\) fixed by the Max+ design (boolean antipodality, \(\mathrm{Tr}(G^2)\), etc.). Certified: \(\|(T-4p)(m_4-m_\star)\|<10^{-11}\) at \(p=5,7\). \(\square\)

5. **Min-norm bound on \(|\kappa|=1\) (certified \(p=5,7\)).**
   | \(p\) | \(\max_{|\kappa|=1}|m_\star|\) | \(L_{\mathrm{abs}}\) | \(\max|m_4|\) (Max+) |
   |------|-------------------------------|----------------------|----------------------|
   | 5 | \(0.056=7/125\) | \(0.06=3/50\) | \(3/65\approx0.04615\) |
   | 7 | \(\approx0.03154\) | \(5/98\approx0.05102\) | \(109/2863\approx0.03807\) |

   So \(m_\star\) itself already obeys \(|m_\star|\le L_{\mathrm{abs}}\) on \(|\kappa|=1\) at \(p=5,7\). The Max+ correction \(h_\star\) **decreases** the max at \(p=5\) and **increases** it at \(p=7\), still staying below \(L_{\mathrm{abs}}\). Best \(L^2\) fit \(m_\star\approx\kappa/(p^2-4)\) (exact coefficient \(1/(p^2-4)\)), but \(m_\star\) is **not** constant on \(\kappa\)-classes for \(p\ge5\). Evidence: `e1_gmin_m4_pseudo.json`. \(\square\)

6. **Residual (OPEN).** Prove for every prime \(p\ge5\):
   - (i) \(\lambda_{\max}(T)=4p\) (and ideally \(\lambda_{\min}=-4p\)) on the Paley / conference Johnson signing;
   - (ii) \(4\kappa/p\perp E_{4p}\);
   - (iii) the Max+ particular solution satisfies \(|m_4|\le L_{\mathrm{abs}}\) on \(|\kappa|=1\)
     (e.g. by controlling \(h_\star\), or by proving \(|m_\star|\le L_{\mathrm{abs}}\) and \(|h_\star|\) cannot push past \(L\), or the candidate \(|m_4|\le(p-2)/(p(2p+3))\), or \(\lambda_2(P\odot P)\le4/N\)).

   Then \(g_{\min}\ge L(p)\) and bi-tight closes. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

   *Remark.* The naive resolvent gain bound of Prop 15.68.5 assumed \(4p>\lambda_{\max}(T)\). At \(p\ge5\) one has equality, so the resolvent is a **pseudoinverse** on \(E_{4p}^\perp\) plus an undetermined kernel component — the kernel is load-bearing, not an error.

**Proposition 15.70 (mid upper bound algebra; bi-tight threshold comparison; multi-worker census; 2026-07-30).** Continue Prop 15.66–15.69. Write
\[
L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2},\qquad
T_{\mathrm{abs}}(p)=\frac{p-2}{p(2p-1)},\qquad
M_{\mathrm{mid}}(p)=\frac{p-2}{2p(p+1)},\qquad
M_{\mathrm{cand}}(p)=\frac{p-2}{p(2p+3)}.
\]

1. **Algebra of targets (proved).** For every odd prime \(p\ge5\):
   \[
   M_{\mathrm{cand}}(p)\;\le\;M_{\mathrm{mid}}(p)\;\le\;L_{\mathrm{abs}}(p)\;<\;T_{\mathrm{abs}}(p),
   \]
   with ratios
   \[
   \frac{M_{\mathrm{mid}}}{L_{\mathrm{abs}}}=\frac{p}{p+1}<1,\qquad
   \frac{M_{\mathrm{cand}}}{L_{\mathrm{abs}}}=\frac{2p}{2p+3}<1.
   \]
   Hence any proof of \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}(p)\) (or the sharper \(M_{\mathrm{cand}}\)) yields \(g_{\min}\ge L(p)>T(p)\) and closes bi-tight via Prop 15.47. \(\square\)

2. **Bi-tight threshold (recalled).** Prop 15.47: if \(g_{\min}>T(p)=-(p-2)/(p(2p-1))\), i.e. \(\max|m_4|<T_{\mathrm{abs}}\), then no Max+-tight size-\(2p\) bi-tight cover exists. The stronger \(g_{\min}\ge L(p)\) is preferred but not necessary for bi-tight. \(\square\)

3. **Multi-worker census (certified \(p=5,7\); \(W=86\)).** Full enumeration of all \(|\kappa|=1\) 4-sets:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | \(M_{\mathrm{mid}}\) | \(L_{\mathrm{abs}}\) | \(T_{\mathrm{abs}}\) |
   |------|---------------|----------------------|---------------------|--------------------|--------------------|
   | 5 | \(3/65\approx0.04615\) | \(3/65\) (sharp) | \(1/20=0.05\) | \(3/50=0.06\) | \(1/15\approx0.0667\) |
   | 7 | \(109/2863\approx0.03807\) | \(5/119\approx0.0420\) | \(5/112\approx0.0446\) | \(5/98\approx0.0510\) | \(\approx0.0544\) |

   In particular \(g_{\min}\ge L(p)>T(p)\) at \(p=5,7\) (bi-tight empty for these \(p\)), and \(\max|m_4|\le M_{\mathrm{mid}}\) holds at both. Evidence: `e1_gmin_m4_close.json`, `e1_gmin_m4_evec4p.json` (F17 multi-worker). \(\square\)

4. **Type6 Max+-free particular solution (certified \(p=5\)).** Solving \((4pI-T)m=4\kappa/p\) on \(S_4\)-type6 class-constant functions (pure \(C\) combinatorics) yields \(\max_{|\kappa|=1}|m|\approx0.0468\le L_{\mathrm{abs}}\) at \(p=5\), close to the true Max+ value \(3/65\). At \(p=3\) type6 recovers the exact \(|m_4|=1/3\). Evidence: `e1_gmin_m4_close.json`. \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}(p)\) (or \(\le L_{\mathrm{abs}}\), or \(\le M_{\mathrm{cand}}\)), using Max+/boolean/conference structure — e.g. via control of \(h_\star\in E_{4p}\) in Prop 15.69, or a closed \(G\)-spectrum / \(\mathrm{Tr}(G^2)\) pin on the moduli line (Prop 15.53). Then bi-tight closes for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.71 (\(\kappa\)-stratum counts for conference matrices; extension-degree feed; 2026-07-29).** Let \(C\) be any real symmetric **conference matrix** of order \(n\) (\(C^\top=C\), zero diagonal, off-diagonal \(\pm1\), \(C^2=(n-1)I\)), and for a 4-set \(S\) write
\[
\kappa(S)=\sum_{\text{three pairings}}C_eC_{e'}\in\{\pm1,\pm3\}.
\]
Let \(n_1=\#\{S:|\kappa(S)|=1\}\) and \(n_3=\#\{S:|\kappa(S)|=3\}\).

1. **Wedge sum from \(C^2\) (proved).** For \(b\neq c\) one has \(\sum_d C_{bd}C_{cd}=0\). Hence for \(a\notin\{b,c\}\),
   \[
   \sum_{d\notin\{a,b,c\}}C_{db}C_{dc}=-C_{ab}C_{ac}.
   \]
   Summing the wedge \(C_{ab}C_{ac}C_{db}C_{dc}\) over distinct \(a,b,c\) and \(d\notin\{a,b,c\}\) therefore yields
   \[
   \Sigma:=-n(n-1)(n-2).
   \]
   \(\square\)

2. **K4 ratio (proved by exhaustion).** On every \(\pm1\)-edge labelling of \(K_4\), writing \(\pi_1,\pi_2,\pi_3\) for the three pairing-products and \(\mathrm{cross}=\pi_1\pi_2+\pi_1\pi_3+\pi_2\pi_3\),
   \[
   \sum_{\sigma\in S_4}C_{\sigma(a)\sigma(b)}C_{\sigma(a)\sigma(c)}C_{\sigma(d)\sigma(b)}C_{\sigma(d)\sigma(c)}
   \;=\;8\cdot\mathrm{cross}.
   \]
   (64 edge labelings; always \(\mathrm{cross}\in\{-1,3\}\).) Consequently \(\sum_{\text{4-sets}}\mathrm{cross}=\Sigma/8=-n(n-1)(n-2)/8\). \(\square\)

3. **Fourth-moment sum of \(\kappa\) (proved).** Since \(\kappa^2=3+2\,\mathrm{cross}\),
   \[
   \sum_S\kappa(S)^2
   =3\binom{n}{4}+2\cdot\frac{\Sigma}{8}
   =\frac{n(n-1)(n-2)(n-5)}{8}.
   \]
   \(\square\)

4. **Stratum counts (proved).** Using \(\kappa^2\in\{1,9\}\) and \(\binom{n}{4}=n_1+n_3\),
   \[
   n_3=\frac{1}{8}\Bigl(\sum\kappa^2-\binom{n}{4}\Bigr)
   =\frac{n(n-1)(n-2)(n-6)}{96},\qquad
   n_1=\binom{n}{4}-n_3
   =\frac{n(n-1)(n-2)^2}{32}.
   \]
   For Paley conferences \(n=p^2+1\) this is
   \[
   n_1=\frac{(p^2+1)\,p^2\,(p^2-1)^2}{32},\qquad
   n_3=\frac{(p^2+1)\,p^2\,(p^2-1)\,(p^2-5)}{96}.
   \]
   \(\square\)

5. **Extension degrees (certified; unique if constant).** For Paley of order \(n=p^2+1\), every \(|\kappa|=1\) 4-set has exactly
   \[
   d_3=p^2-5,\qquad d_1=3p^2-7
   \]
   ordered one-vertex extensions landing in \(|\kappa|=3\) and \(|\kappa|=1\) respectively (so \(d_1+d_3=4(n-4)\)). Constancy certified by full census at \(p=3,5\) and multi-worker samples at \(p=7,11\). Combined with part 4, the constant-degree hypothesis is the unique solution of the handshaking identities. Evidence: `e1_gmin_m4_stratum.json`, `e1_gmin_m4_tkappa.json`. \(\square\)

6. **Target algebra (proved; recalled).** For odd primes \(p\ge5\),
   \[
   M_{\mathrm{cand}}(p)\le M_{\mathrm{mid}}(p)\le L_{\mathrm{abs}}(p)<T_{\mathrm{abs}}(p)
   \]
   with \(M_{\mathrm{mid}}/L_{\mathrm{abs}}=p/(p+1)\) (Prop 15.70.1). \(\square\)

7. **Multi-worker census (certified \(p=3,5,7,11\); \(W=86\)).** Full \(\kappa\)-stratum counts match part 4 at all four primes (including \(p=11\), \(\binom{122}{4}\approx8.7\cdot10^6\)). Evidence: `e1_gmin_m4_stratum.json` (F17 ProcessPool; atomic JSON). \(\square\)

8. **Residual (OPEN).** Parts 1–4 are Max+-free conference combinatorics and pin the source size of the resolvent equation \((4pI-T)\rho=T\kappa/p^2\) (Prop 15.68): exactly \(n_3\) nonzero source coordinates of amplitude \(24/p^2\). They do **not** alone bound \(|m_4|\) on \(|\kappa|=1\). Still open for every prime \(p\ge5\): \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\) / \(M_{\mathrm{cand}}\)), e.g. via resolvent gain \(\le(p-4)/48\), \(h_\star\in E_{4p}\) control, or type6 association closed form. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.72 (resolvent-gain calculus from \(n_1,n_3,d_1,d_3\); reverse degrees; type6 pin; 2026-07-29).** Continue Prop 15.68–15.71. Write \(\rho:=m_4-\kappa/p^2\), so \((4pI-T)\rho=T\kappa/p^2\), and on every \(|\kappa|=1\) 4-set with \(\mathrm{sign}(\rho)=\mathrm{sign}(\kappa)\),
\[
|m_4|=\frac1{p^2}+|\rho|.
\]

1. **Gain \(\Leftrightarrow L\) algebra (proved).** For primes \(p\ge5\),
   \[
   |\rho|\;\le\;\frac{p-4}{2p^2}
   \quad\Longleftrightarrow\quad
   |m_4|\;\le\;L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2}
   \]
   on same-sign \(|\kappa|=1\) classes. Since the inhomogeneous source has amplitude \(24/p^2\) on \(|\kappa|=3\), it is sufficient that the resolvent gain from that source into same-sign \(|\kappa|=1\) satisfy
   \[
   \mathrm{gain}\;\le\;\frac{p-4}{48},
   \]
   because \(\frac{p-4}{48}\cdot\frac{24}{p^2}=\frac{p-4}{2p^2}\). \(\square\)

2. **Source sign structure (proved).** On every \(|\kappa|=3\) labelling of \(K_4\), \(T\kappa/\kappa\in\{\pm8\}\) (64-labeling; \(T\kappa=-6\cdot\mathrm{star}\)). Thus \(|T\kappa|=24\) with sign free relative to \(\kappa\). \(\square\)

3. **Reverse extension degrees (proved under Prop 15.68 constancy).** Assume \(d_3=p^2-5\) is constant on every Paley \(|\kappa|=1\) 4-set. Handshaking on the bipartite extension graph between \(|\kappa|=1\) and \(|\kappa|=3\) strata, together with the counts \(n_1,n_3\) of Prop 15.71, forces the degrees from every \(|\kappa|=3\) 4-set:
   \[
   d_1^{(3)}=\frac{n_1\,d_3}{n_3}=3(p^2-1),\qquad
   d_3^{(3)}=4(n-4)-d_1^{(3)}=p^2-9.
   \]
   (Both nonnegative for primes \(p\ge3\), with \(d_3^{(3)}=0\) at \(p=3\).) \(\square\)

4. **Separate \(\kappa\)-weighted vanishing (certified \(p=3,5,7\); open as general theorem).** On every Paley \(|\kappa|=1\) 4-set the one-step sums split by target stratum vanish separately:
   \[
   \sum_{\mathrm{ext}\to|\kappa|=1}C_{vr}\,\kappa(S')=0,\qquad
   \sum_{\mathrm{ext}\to|\kappa|=3}C_{vr}\,\kappa(S')=0.
   \]
   (Each is stronger than \(T\kappa=0\), which is only their sum.) Full multi-worker census at \(p=3,5,7\). Evidence: `e1_gmin_m4_resolvent_gain.json`. \(\square\)

5. **Reverse-degree census (certified \(p=3,5,7\)).** Every \(|\kappa|=3\) 4-set has \((d_1^{(3)},d_3^{(3)})=(3(p^2-1),p^2-9)\) constantly, and \(T\kappa/\kappa\in\{\pm8\}\). \(\square\)

6. **Type6 Max+-free resolvent (certified \(p=3,5,7\); \(W=86\)).** Restricting to \(S_4\)-type6 class-constant functions and solving \((4pI-T)\rho=T\kappa/p^2\) in the least-squares sense:
   | \(p\) | \(\max_{|\kappa|=1}|m_{\mathrm{type6}}|\) | \(L_{\mathrm{abs}}\) | same-sign \(|\rho|\) | gain | budget \(\frac{p-4}{48}\) |
   |------|--------------------------------------|--------------------|----------------------|------|-------------------------------|
   | 5 | \(\approx0.04764\) | \(0.06\) | \(\approx0.00764\) | \(\approx0.00796\) | \(0.02083\) |
   | 7 | \(\approx0.02407\) | \(0.0510\) | \(\approx0.00367\) | \(\approx0.00749\) | \(0.0625\) |

   In particular type6 predicts \(|m|\le L_{\mathrm{abs}}\) and gain below budget at \(p=5,7\). (At \(p=7\) true Max+ \(\max|m_4|\approx0.038> m_{\mathrm{type6}}\), so type6 is not exact — classes need refinement — but remains a Max+-free upper probe.) Evidence: `e1_gmin_m4_resolvent_gain.json`. \(\square\)

7. **Empirical Max+ gain (certified \(p=5,7\); mmap).** True same-sign residual gains \(0.00641\) and \(0.03606\) both lie strictly below \(\frac{p-4}{48}\); \(\max|m_4|\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\). \(\square\)

8. **Residual (OPEN).** Prove for every prime \(p\ge5\) that the resolvent gain is \(\le(p-4)/48\), or directly \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\)), using the stratum data of Props 15.71–15.72 (source size \(n_3\), degrees \(d_1,d_3,d_1^{(3)},d_3^{(3)}\), separate vanishing) plus Max+/boolean structure — without a per-prime Max+ census as the proof. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.73 (Max+ \(e_4\) identity; Paley \(\sum\kappa\); multi-prime type6 gain; 2026-07-29).** Continue Prop 15.52 and 15.71–15.72.

1. **Boolean \(e_4\) from the Max+ sum constraint (proved).** By Prop 15.52, every \(y\in\mathrm{Max}_{+}\) satisfies \(\mathbf1^\top y=(p+1)y_\infty\), hence
   \[
   s^2:=(\textstyle\sum_i y_i)^2=(p+1)^2
   \]
   constantly. For any boolean vector \(y\in\{\pm1\}^n\) one has the power-sum identity
   \[
   e_4(y):=\sum_{a<b<c<d}y_ay_by_cy_d
   =\frac{s^4-6ns^2+3n^2+8s^2-6n}{24}.
   \]
   Substituting \(s^2=(p+1)^2\) and \(n=p^2+1\) yields the closed form
   \[
   e_4=-\frac{p(p-1)(p+1)(p+4)}{12}.
   \]
   Since \(e_4(y)\) is constant on \(\mathrm{Max}_{+}\),
   \[
   \sum_{S}m_4(S)=\mathbb E[e_4]=e_4.
   \]
   Certified on full Max+ at \(p=5,7\) (mmap). Evidence: `e1_gmin_m4_e4_gain.json`. \(\square\)

2. **Paley sum of \(\kappa\) (formula; certified \(p=3,5,7,11,13\)).** For the Paley conference of order \(n=p^2+1\),
   \[
   \sum_S\kappa(S)=\frac{p^2(p^2-1)}{4}.
   \]
   Full multi-worker \(\kappa\)-sum census matches at all five primes (including \(p=13\), \(\binom{170}{4}\approx3.4\cdot10^7\)). Combined with part 1,
   \[
   \sum_S\rho(S)=e_4-\frac1{p^2}\sum_S\kappa(S)
   \]
   is an exact Max+/Paley scalar. \(\square\)

3. **Type6 Max+-free resolvent across primes (certified \(p=5,7,11,13\); \(W=86\)).** All \(11\) abstract \(S_4\)-type6 edge-orbits appear for \(p\ge5\). Solving \((4pI-T)\rho=T\kappa/p^2\) in the type6-constant subspace:
   | \(p\) | \(\max_{|\kappa|=1}|m_{\mathrm{type6}}|\) | \(M_{\mathrm{mid}}\) | \(L_{\mathrm{abs}}\) | gain | budget |
   |------|--------------------------------------|----------------------|--------------------|------|--------|
   | 5 | \(0.04764\) | \(0.0500\) | \(0.0600\) | \(0.00796\) | \(0.0208\) |
   | 7 | \(0.02337\) | \(0.0446\) | \(0.0510\) | \(0.00604\) | \(0.0625\) |
   | 11 | \(0.00876\) | \(0.0341\) | \(0.0372\) | \(0.00249\) | \(0.1458\) |
   | 13 | \(0.00605\) | \(0.0302\) | \(0.0325\) | \(0.00092\) | \(0.1875\) |

   In particular type6 predicts \(|m|\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\) and gain \(\ll(p-4)/48\) at every tested prime, with \(\max|m_{\mathrm{type6}}|\cdot p^2\to1\) (Wick scale). \(\square\)

4. **Caveat (proved by comparison).** Type6 is **not** always exact for true Max+ \(m_4\): at \(p=7\), \(\max|m_4|\approx0.03807> m_{\mathrm{type6}}\approx0.023\). Thus type6 is a Max+-free **probe** (and a feasible particular solution of the master linear equation in a \(T\)-invariant subspace), not by itself an upper bound on true \(m_4\). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that true Max+ satisfies \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{mid}}\) (or \(L_{\mathrm{abs}}\)), e.g. by:
   - closing the type6 error \(m_4-m_{\mathrm{type6}}\in E_{4p}\) with a kernel bound, or
   - proving resolvent gain \(\le(p-4)/48\) from reverse degrees + separate vanishing (Prop 15.72), or
   - a character-sum formula for Paley \(m_4\).

   Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.74 (candidate bound algebra; true Max+ census; signed residual identity; 2026-07-29).** Continue Prop 15.68–15.73. Write \(M_{\mathrm{cand}}(p)=(p-2)/(p(2p+3))\) and \(\rho=m_4-\kappa/p^2\).

1. **Candidate algebra (proved).** For every odd prime \(p\ge5\),
   \[
   M_{\mathrm{cand}}(p)\;\le\;M_{\mathrm{mid}}(p)\;\le\;L_{\mathrm{abs}}(p)\;<\;T_{\mathrm{abs}}(p),
   \]
   with \(M_{\mathrm{cand}}/L_{\mathrm{abs}}=2p/(2p+3)\). Same-sign residual budgets:
   \[
   |m_4|\le M_{\mathrm{cand}}
   \;\Longleftrightarrow\;
   |\rho|\le M_{\mathrm{cand}}-\frac1{p^2}
   \;\Longleftrightarrow\;
   \mathrm{gain}\le\frac{p^2-4p-3}{24(2p+3)},
   \]
   where gain is relative to source amplitude \(24/p^2\). At \(p=5\) the candidate gain is \(1/156\) and \(M_{\mathrm{cand}}=3/65\). \(\square\)

2. **Signed residual identity (proved).** On every \(|\kappa|=1\) 4-set the residual equation \(4p\rho=T\rho\) (Prop 15.68, using \(T\kappa=0\)) multiplies by \(\kappa\) to give
   \[
   4p\,r=\kappa\,(T\rho),\qquad r:=\rho\cdot\kappa.
   \]
   Same-sign danger for \(|m_4|\) is exactly \(r>0\), with \(|m_4|=1/p^2+r\). Certified float residual \(<10^{-15}\) on multi-worker Max+ samples at \(p=5,7\). \(\square\)

3. **True Max+ census (certified \(p=5,7\); \(W=86\); mmap — not type6).** Full enumeration of all \(|\kappa|=1\) 4-sets against Max+:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | \(M_{\mathrm{mid}}\) | same-sign \(\max r\) | gain | \(\mathrm{gain}_{\mathrm{cand}}\) |
   |------|---------------|----------------------|---------------------|----------------------|------|-------------------------------|
   | 5 | \(3/65\approx0.046154\) | \(3/65\) **(sharp)** | \(0.05\) | \(2/325\) | \(1/156\) | \(1/156\) |
   | 7 | \(109/2863\approx0.038072\) | \(5/119\approx0.04202\) | \(0.04464\) | \(\approx0.01766\) | \(\approx0.03606\) | \(\approx0.01330\) |

   In particular **true** Max+ (not type6) satisfies \(\max|m_4|\le M_{\mathrm{cand}}\le M_{\mathrm{mid}}\le L_{\mathrm{abs}}\) at \(p=5,7\), with equality in the candidate at \(p=5\). Evidence: `e1_gmin_m4_kernel.json`. \(\square\)

4. **Kernel form (recalled).** Every solution of the master linear equation is \(m_4=m_\star+h\) with \(h\in E_{4p}=\ker(4pI-T)\) and \(m_\star\) the min-norm particular solution (Prop 15.69). True Max+ selects a unique \(h_\star\). At \(p=5\) the candidate is already sharp on true \(m_4\), so \(h_\star\) cannot increase the max beyond \(m_\star\)'s ceiling in the dangerous direction; at \(p=7\), \(h_\star\) raises \(\max|m_4|\) above type6 but still stays \(\le M_{\mathrm{cand}}\). \(\square\)

5. **Residual (OPEN).** Prove for every prime \(p\ge5\) that true Max+ obeys
   \[
   \max_{|\kappa|=1}|m_4|\;\le\;M_{\mathrm{cand}}(p)
   \]
   (or the weaker \(M_{\mathrm{mid}}\) / gain \(\le(p-4)/48\)), using \(E_{4p}\) control of \(h_\star\), the signed operator \(r\mapsto\kappa(T\rho)\), reverse degrees, and/or Paley character sums — **without** replacing true Max+ by type6. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.75 (one-center \(\sigma_a=2\cdot\mathrm{star}_a\); K4 Gram spectrum; GPU cand cert; 2026-07-29).** Continue Prop 15.67–15.74.

1. **One-center \(\kappa\)-sum (proved, any conference matrix).** For every 4-set \(S\) with \(|\kappa(S)|=1\) and every \(a\in S\),
   \[
   \sigma_a:=\sum_{r\notin S}C_{ar}\,\kappa(S_{a\to r})
   \;=\;2\cdot\mathrm{star}_a,
   \qquad
   \mathrm{star}_a:=\prod_{u\in S\setminus\{a\}}C_{au}\in\{\pm1\}.
   \]
   *Proof.* \(C^2=(n-1)I\) gives \(\sum_{r\notin\{a,b\}}C_{ar}C_{br}=-\sum_{u\in S\setminus\{a,b\}}C_{au}C_{ub}\) for \(b\in S\setminus\{a\}\). Expanding \(\sigma_a\) in the three pairing products of \(\kappa(S_{a\to r})\) yields \(\sigma_a=-2\sum_{v\in S\setminus\{a\}}\mathrm{star}_v\). On \(|\kappa|=1\), \(T\kappa=0\) forces \(\sum_{v\in S}\mathrm{star}_v=0\), hence \(\sum_{v\neq a}\mathrm{star}_v=-\mathrm{star}_a\) and \(\sigma_a=2\,\mathrm{star}_a\). \(\square\)

2. **One-center residual form (proved).** With \(\rho=m_4-\kappa/p^2\), the evec identity at centre \(a\) becomes
   \[
   p\,m_4-\frac{\kappa}{p}
   \;=\;\frac{2\,\mathrm{star}_a}{p^2}
   +\sum_{r\notin S}C_{ar}\,\rho(S_{a\to r}).
   \]
   (Certified float residual \(<10^{-15}\) on Max+ samples at \(p=5,7\).) Averaging over \(a\in S\) recovers \(\mathrm{Ext}=4p\rho\). \(\square\)

3. **Local K4 edge-Gram spectrum (proved algebra).** On the six edges of a \(|\kappa|=1\) 4-set, the principal submatrix of \(G=\mathbb E[ff^\top]\) has opposite-edge entries \(\pi_i m_4\) and wedge entries \(\pm1/p\). Its eigenvalues lie in the pool
   \[
   \{1\pm m_4\}\ \cup\ \{1\pm m_4\pm 2/p\}\ \cup\ \{1\pm m_4\pm 2\sqrt2/p\}.
   \]
   PSD of this block forces the weak general bound \(|m_4|\le1-2/p=(p-2)/p\) whenever \(1-|m_4|-2/p\) is an eigenvalue — far weaker than \(M_{\mathrm{cand}}\), but Max+-free and load-bearing for local structure. \(\square\)

4. **GPU true Max+ cand census (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full \(|\kappa|=1\) m4 on Max+ with one CUDA context, device argmax, mmap Max+ load, atomic evidence write:
   | \(p\) | \(\max|m_4|\) | \(M_{\mathrm{cand}}\) | same-sign \(r\) | gain | wall |
   |------|---------------|----------------------|----------------|------|------|
   | 5 | \(3/65\) (sharp) | \(3/65\) | \(2/325\) | \(1/156\) | \(\sim0.36\)s |
   | 7 | \(109/2863\) | \(5/119\) | \(\approx0.01766\) | \(\approx0.0361\) | \(\sim0.30\)s |

   Evidence: `e1_gmin_m4_gpu.json` (`gpu.used=true`, `io.mmap+atomic`). Multi-worker \(\sigma\) census \(p=3,5,7,11\): `e1_gmin_m4_onecenter.json`. \(\square\)

5. **Residual (OPEN).** Bound \(\sum_r C_{ar}\rho(S_{a\to r})\) (or the signed global operator of Prop 15.74) tightly enough that part 2 forces \(m_4\le M_{\mathrm{cand}}\) for every prime \(p\ge5\). Absolute degree bounds do not contract. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.76 (one-center extension degrees; residual split; 2026-07-29).** Continue Prop 15.75. Write \(d_1=3p^2-7\), \(d_3=p^2-5\) for the full ordered-extension degrees on \(|\kappa|=1\) (Prop 15.68/15.72).

1. **Divisibility (proved).** For every odd integer \(p\ge3\), \(4\mid(3p^2-7)\) and \(4\mid(p^2-5)\). Hence
   \[
   d_1^{(1)}:=\frac{3p^2-7}{4},\qquad d_3^{(1)}:=\frac{p^2-5}{4}
   \]
   are integers, and \(d_1^{(1)}+d_3^{(1)}=p^2-3=n-4\), \(4d_1^{(1)}=d_1\), \(4d_3^{(1)}=d_3\). \(\square\)

2. **One-center degree constancy (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) 4-set \(S\) and every centre \(a\in S\), among the \(n-4\) one-vertex extensions \(S_{a\to r}\), exactly \(d_1^{(1)}\) land in \(|\kappa|=1\) and \(d_3^{(1)}\) land in \(|\kappa|=3\). (So the full \(4(n-4)\) count splits evenly across the four centres.) Evidence: `e1_gmin_m4_onecenter_deg.json`. \(\square\)

3. **Residual split (proved form).** With \(\rho=m_4-\kappa/p^2\) and \(\sigma_a=2\cdot\mathrm{star}_a\) (Prop 15.75), for each centre \(a\) on a \(|\kappa|=1\) set with \(\kappa=1\),
   \[
   p\rho-\frac{2\,\mathrm{star}_a}{p^2}
   \;=\;S_1(a)+S_3(a),
   \]
   where \(S_j(a)=\sum C_{ar}\rho(S_{a\to r})\) runs over extensions to \(|\kappa|=j\). In particular \(|S_3(a)|\le d_3^{(1)}R_3\) and \(|S_1(a)|\le d_1^{(1)}R_1\) with \(R_j=\max|\rho|\) on the \(|\kappa|=j\) stratum. \(\square\)

4. **Absolute bootstrap fails (proved).** The four-centre form \(4p R_1\le d_1 R_1+d_3 R_3\) rearranges to \(R_1(4p-d_1)\le d_3 R_3\). For primes \(p\ge5\), \(4p-d_1=4p-(3p^2-7)<0\), so this yields only a lower bound on \(R_1\), not an upper bound. Signed cancellation on \(S_1\) (or a design bound coupling \(R_1,R_3\)) is load-bearing. \(\square\)

5. **GPU residual moments (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full-quad m4 on Max+:
   | \(p\) | \(\max|m_4|_{\kappa=1}\) | \(\max|m_4|_{\kappa=3}\) | \(R_1\) | \(R_3\) | \(\le M_{\mathrm{cand}}\) |
   |------|--------------------------|--------------------------|--------|--------|------------------------|
   | 5 | \(3/65\) | \(21/65\) | \(0.0554\) | \(0.2031\) | yes |
   | 7 | \(109/2863\) | \(\approx0.1142\) | \(0.0177\) | \(0.0530\) | yes |

   Note \(R_1>\max r\) (same-sign residual) because opposite-sign \(\rho\) on \(|\kappa|=1\) is larger; only same-sign \(r\) raises \(|m_4|\) above Wick. Evidence: `e1_gmin_m4_onecenter_deg.json`. \(\square\)

6. **Residual (OPEN).** Prove a signed bound on \(S_1+S_3\) (e.g. \(S_1\le0\) at \(\mathrm{star}_a=+1\) on maximisers, or a Paley character-sum formula for \(\rho\)) strong enough that part 3 forces \(\rho\le(p^2-4p-3)/(p^2(2p+3))\) on same-sign \(|\kappa|=1\), i.e. \(m_4\le M_{\mathrm{cand}}\), for every prime \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.77 (\(\mathrm{star}_a\cdot S_1\le0\) structure; joint cand criterion; 2026-07-29).** Continue Prop 15.75–15.76. Write \(S_1(a),S_3(a)\) for the one-center residual split on a \(|\kappa|=1\) 4-set and centre \(a\), and \(\mathrm{star}_a=\prod_{u\in S\setminus\{a\}}C_{au}\).

1. **Star–joint identity (proved, Max+-free).** From the one-center residual form,
   \[
   p\rho-\frac{2\,\mathrm{star}_a}{p^2}=S_1(a)+S_3(a)
   \qquad\Longrightarrow\qquad
   \mathrm{star}_a\cdot\bigl(S_1+S_3\bigr)
   \;=\;
   p\,\rho\,\mathrm{star}_a-\frac{2}{p^2}.
   \]
   In particular, if \(\mathrm{star}_a=+1\) and \(S_1(a)\le0\), then
   \[
   p\rho\;\le\;\frac{2}{p^2}+S_3(a).
   \]
   \(\square\)

2. **Same-sign reconstruction (proved).** On every \(|\kappa|=1\) set with \(\kappa=1\) and \(\rho>0\), at every centre with \(\mathrm{star}_a=+1\),
   \[
   \rho
   \;=\;
   \frac{1}{p}\Bigl(\frac{2}{p^2}+S_1(a)+S_3(a)\Bigr).
   \]
   Hence
   \[
   \max_{\text{same-sign }r>0}\rho
   \;=\;
   \max_{\substack{\mathrm{star}_a=+1\\ r>0}}
   \frac{1}{p}\Bigl(\frac{2}{p^2}+S_1+S_3\Bigr),
   \]
   and the candidate bound \(\rho\le\rho_{\mathrm{cand}}:=M_{\mathrm{cand}}-1/p^2\) is equivalent to
   \[
   \max_{\substack{\mathrm{star}_a=+1\\ r>0}}\bigl(S_1+S_3\bigr)
   \;\le\;
   p\,\rho_{\mathrm{cand}}-\frac{2}{p^2}.
   \]
   At \(p=5\) the right-hand side equals \(-16/325\approx-0.04923<0\), so \(S_1\le0\) alone with absolute \(|S_3|\) bounds is **not** enough — maximisers must have strongly negative \(S_1\). \(\square\)

3. **GPU census of \(\mathrm{star}_a\cdot S_1\) (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full Max+ m4 (one CUDA context, single H2D, D2H m4 vector only) + ProcessPool walk of every \(|\kappa|=1\) centre:
   | \(p\) | \(\#\) star\(+\) checks | \(\max(\mathrm{star}\cdot S_1)\) | \(\max S_1\) at star\(+\) | joint \(S_1{+}S_3\) on \(r>0\) | \(\Rightarrow\rho\le\rho_{\mathrm{cand}}\) |
   |------|--------------------------|----------------------------------|---------------------------|----------------------------------|----------------------------------|
   | 5 | \(23400\) | \(-0.03077\) (strict \(<0\)) | \(-0.03077\) | \(-0.04923\) | **yes (sharp)** |
   | 7 | \(352800\) | \(-0.00669\) (strict \(<0\)) | \(-0.00669\) | \(0.08283\) | **yes** |

   In particular **\(\mathrm{star}_a\cdot S_1(a)\le0\) on every \(|\kappa|=1\) centre** (not only on same-sign maximisers), with perfect sign-antisymmetry \(\max S_1|_{\mathrm{star}+}=-\min S_1|_{\mathrm{star}-}\). Identity residual \(<10^{-15}\). Evidence: `e1_gmin_m4_S1_star.json`. \(\square\)

4. **Max+ is essential (certified).** A synthetic residual that puts a uniform same-sign bump on all \(|\kappa|=1\) coordinates (and zero on \(|\kappa|=3\)) **violates** \(\mathrm{star}\cdot S_1\le0\) (\(\max\mathrm{star}\cdot S_1=+0.05\) at \(p=5,7\)). So the inequality is a property of true Max+ fourth moments, not of the \(\kappa\)-adjacency graph alone. Combinatorial \(\tau_1:=\sum_{\kappa1}C_{ar}\kappa'\) has \(\mathrm{star}\cdot\tau_1\) non-constant for \(p\ge5\) (values in \(\{-1,5\}\) at \(p=5\); multi-worker pure-\(C\) census). \(\square\)

5. **Residual (OPEN).** Prove \(\mathrm{star}_a\cdot S_1(a)\le0\) for every prime \(p\ge5\) on Paley Max+ (e.g. via boolean \(+p\)-evec character sums), and prove the joint criterion of part 2 (or a matching \(S_3\) bound on maximisers) so that \(\max|m_4|\le M_{\mathrm{cand}}\) for all such \(p\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.78 (moment form / Gaussian domination; 4-set constancy; exact \(p=5\) spectrum; 2026-07-29).** Continue Prop 15.77. Write \(R_1(a)=\{r\notin S:|\kappa(S_{a\to r})|=1\}\), \(U_1=\sum_{r\in R_1}C_{ar}y_r\), \(\tau_1(a)=\sum_{r\in R_1}C_{ar}\kappa(S_{a\to r})\), and \(f_0=\prod_{u\in S\setminus\{a\}}y_u\).

1. **Moment form (proved, Max+ with \(Cy=py\)).** On every \(|\kappa|=1\) centre,
   \[
   \mathrm{star}_a\cdot S_1(a)
   \;=\;
   \mathbb E\bigl[\varphi\bigr]-\mathbb E_{\mathrm{Wick}}\bigl[\varphi\bigr],
   \qquad
   \varphi(y)\,:=\,\mathrm{star}_a\,f_0\,U_1,
   \]
   where Wick uses the Gaussian \(N(0,\Sigma)\) with \(\Sigma=I+C/p\) (so \(\mathbb E_{\mathrm{Wick}}[\varphi]=\mathrm{star}_a\tau_1/p^2\)). Consequently
   \[
   \mathrm{star}_a\cdot S_1(a)\le0
   \quad\Longleftrightarrow\quad
   \mathbb E_{\mathrm{Max}+}[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]
   \]
   (Gaussian domination for this specific degree-4 form). \(\square\)

2. **Joint rewrite (proved).** For \(\kappa=1\) and \(\mathrm{star}_a=+1\),
   \[
   S_1+S_3 \;=\; p\,m_4-\frac1p-\frac2{p^2}.
   \]
   (So the joint \(S_1+S_3\) is a function of \(m_4\) alone on this locus; the cand criterion is equivalent to \(\max m_4\le M_{\mathrm{cand}}\) on same-sign sets — the leverage of \(S_1\le0\) is to replace the joint by an \(S_3\)-only upper bound.) \(\square\)

3. **Combinatorial constancy of \(\mathrm{star}\cdot\tau_1\) (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) 4-set, the four values \(\mathrm{star}_a\cdot\tau_1(a)\) (\(a\in S\)) are equal. Observed value sets: \(\{-1\}\) at \(p=3\); \(\{-1,5\}\) at \(p=5\); \(\{-7,-1,5\}\) at \(p=7\); \(\{-13,-7,-1,5,11\}\) at \(p=11\). Evidence: `e1_gmin_m4_S1_const.json`. \(\square\)

4. **Max+ constancy of \(\mathrm{star}\cdot S_1\) (certified GPU \(p=5,7\); mmap+atomic).** On every \(|\kappa|=1\) 4-set, the four values \(\mathrm{star}_a\cdot S_1(a)\) are equal (one rational per set). In particular the 4-set carries a single sign for the GD inequality. Full census: \(\mathrm{star}\cdot S_1\le0\) with Gaussian domination on every centre; \(\max\mathrm{star}\cdot S_1=-2/65\) (\(p=5\)), \(-0.006686\) (\(p=7\)). \(\square\)

5. **Exact spectrum at \(p=5\) (certified full Max+).** \(\mathrm{star}\cdot S_1\in\{-2/65,-42/325\}\), both \(<0\). Matches the closed moment rule
   \[
   \mathrm{star}\cdot\mathbb E[f_0U_1]
   \;=\;
   \frac{11}{65}\,\mathrm{sgn}(\mathrm{star}\cdot\tau_1)
   \]
   with \(\mathrm{star}\cdot\tau_1\in\{-1,5\}\) (so \(\mathrm{star}\cdot S_1=\frac{11}{65}\mathrm{sgn}(t)-\,t/25\)). Hence \(\mathrm{star}\cdot S_1\le0\) is **proved at \(p=5\)** by exhaustive Max+ evaluation. \(\square\)

6. **Residual (OPEN).** Constancy of \(\mathrm{star}\cdot\tau_1\) and \(\mathrm{star}\cdot S_1\) is proved in Prop 15.79. Remains: Gaussian domination \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\) (equivalently \(\mathrm{star}\cdot S_1\le0\)) for all primes \(p\ge5\), and joint/\(S_3\) bound \(\Rightarrow\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.79 (Aut-constancy of \(\mathrm{star}\cdot\tau_1\) and \(\mathrm{star}\cdot S_1\); modular \(\tau_1\); 2026-07-29).** Continue Prop 15.77–15.78. Let \(C\) be the Paley conference matrix of order \(n=p^2+1\), identified with the point set \(\mathrm{PG}(1,\mathbb F_{p^2})\), and let \(\mathrm{Aut}(C)\ge\mathrm{PGL}(2,p^2)\).

1. **Setwise stabilizer is transitive (proved).** For every 4-set \(S\subset\mathrm{PG}(1,p^2)\), the setwise stabilizer of \(S\) in \(\mathrm{PGL}(2,p^2)\) contains a Klein four-group \(V_4\) of double-transposition involutions (e.g. \(z\mapsto\lambda/z\) on \(S=\{\infty,0,1,\lambda\}\)). This \(V_4\) acts regularly on \(S\), hence transitively. \(\square\)

2. **Constancy of \(\mathrm{star}\cdot\tau_1\) (proved, Max+-free).** The scalar \(\mathrm{star}_a\cdot\tau_1(a)\) is built from \(C\)-labels of \(S\) and of the one-vertex extensions \(S_{a\to r}\). It is therefore \(\mathrm{Aut}(C)\)-equivariant: for every \(g\in\mathrm{Aut}(C)\),
   \[
   \mathrm{star}_{ga}\cdot\tau_1(ga;gS)
   \;=\;
   \mathrm{star}_a\cdot\tau_1(a;S).
   \]
   Restricting to the setwise stabilizer of \(S\) and using part 1, \(\mathrm{star}_a\cdot\tau_1(a)\) is independent of \(a\in S\). \(\square\)

3. **Constancy of \(\mathrm{star}\cdot S_1\) (proved, Max+ Aut-invariant).** The set \(\mathrm{Max}_{+}=\{y\in\{\pm1\}^n:Cy=py\}\) is \(\mathrm{Aut}(C)\)-invariant, so \(m_4(S)=\mathbb E[\prod_{i\in S}y_i]\) is an Aut-invariant of the 4-set. Consequently \(\mathrm{star}_a\cdot S_1(a)\) (depending on \(C\) and on \(m_4\) of extensions of \(S\)) is Aut-equivariant, and part 1 forces constancy on \(a\in S\). \(\square\)

4. **Counting form of \(\tau_1\) (proved).** Write \(d_1^{(1)}=(3p^2-7)/4\) and
   \[
   A\;=\;\#\{r\in R_1(a):C_{ar}\kappa(S_{a\to r})=+1\}.
   \]
   Then \(A+B=d_1^{(1)}\) with \(B=\#\{C_{ar}\kappa'=-1\}\), and \(\tau_1=A-B=2A-d_1^{(1)}\). In particular \(\mathrm{star}\cdot\tau_1=\mathrm{star}\cdot(2A-d_1^{(1)})\). For odd \(p\), \(d_1^{(1)}\) is odd, so \(\mathrm{star}\cdot\tau_1\) is always odd. \(\square\)

5. **Modular census (certified Paley \(p=3,5,7,11\); \(W=86\)).** On every \(|\kappa|=1\) set: constancy holds (matches part 2); \(t_1=2A-d_1^{(1)}\) exactly; every \(\mathrm{star}\cdot\tau_1\equiv5\pmod6\); the number of distinct values is \((p-1)/2\), with observed sets
   \[
   \begin{align*}
   p=3&:\ \{-1\},\\
   p=5&:\ \{-1,5\},\\
   p=7&:\ \{-7,-1,5\},\\
   p=11&:\ \{-13,-7,-1,5,11\}
   \end{align*}
   \]
   (arithmetic progressions of difference \(6\)). Evidence: `e1_gmin_m4_S1_aut.json`. \(\square\)

6. **Consequence for the sign attack.** By part 3 it suffices to prove \(\mathrm{star}\cdot S_1(S)\le0\) as a property of the 4-set (one check per set, not per centre). Combined with Prop 15.77–15.78 (joint criterion; Gaussian domination form; exact \(p=5\) spectrum \(\{-2/65,-42/325\}\)), the residual is: prove \(\mathrm{star}\cdot S_1\le0\) for all primes \(p\ge5\), then control \(S_3\) on maximisers. \(\square\)

7. **Residual (OPEN).** Prove \(\mathrm{star}\cdot S_1\le0\) (e.g. Gaussian domination \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) for every prime \(p\ge5\), and close the joint/\(S_3\) bound so \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.80 (linear-form Wick identity; GD formulation; U1-specialization; 2026-07-29).** Continue Prop 15.77–15.79. Write \(Z=\mathrm{star}_a\prod_{u\in S\setminus\{a\}}y_u\) and \(U_1=\sum_{r\in R_1(a)}C_{ar}y_r\).

1. **Linear-form Wick identity (proved).** On Max+ one has \(\mathbb E[y_iy_j]=\Sigma_{ij}\) with \(\Sigma=I+C/p\). Consequently for every \(\alpha\in\mathbb R^n\) and \(L=\sum\alpha_i y_i\),
   \[
   \mathbb E[L^2]
   \;=\;
   \sum_{i,j}\alpha_i\alpha_j\Sigma_{ij}
   \;=\;
   \mathbb E_{G\sim N(0,\Sigma)}[L(G)^2].
   \]
   In particular \(\mathbb E[U_1^2]=\mathbb E_{\mathrm{Wick}}[U_1^2]\). (Only pairwise moments; no fourth-moment input.) \(\square\)

2. **GD formulation (proved).** Central symmetry gives \(\mathbb E[Z]=0\), and
   \[
   \mathrm{star}_a\cdot S_1(a)
   \;=\;
   \mathbb E[ZU_1]-\mathbb E_{\mathrm{Wick}}[ZU_1]
   \;=\;
   \mathrm{Cov}(Z,U_1)-\mathrm{Cov}_{\mathrm{Wick}}(Z,U_1),
   \]
   with \(\mathbb E_{\mathrm{Wick}}[ZU_1]=\mathrm{star}_a\cdot\tau_1(a)/p^2\). Hence
   \[
   \mathrm{star}\cdot S_1\le0
   \quad\Longleftrightarrow\quad
   \mathbb E[ZU_1]\le\mathbb E_{\mathrm{Wick}}[ZU_1].
   \]
   \(\square\)

3. **U1-specialization (proved necessity of \(\kappa=1\)-support; certified generically).** The comparison \(\mathbb E[ZL]\le\mathbb E_{\mathrm{Wick}}[ZL]\) **fails** for generic linear \(L\) supported off \(S\) (violation rate \(\approx 45\)–\(55\%\) at \(p=5,7\)). It can also fail for the full external form \(U_{\mathrm{ext}}=\sum_{r\notin S}C_{ar}y_r\) (equivalent to \(\mathrm{star}\cdot p\cdot\rho\le 2/p^2\), false on some \(p=7\) maximisers with positive joint \(S_1+S_3\)). Thus GD is a property of the \(\kappa=1\)-restricted form \(U_1\), not a general cubic\(\times\)linear inequality. \(\square\)

4. **Sum of \(\mathrm{star}\cdot\tau_1\) (certified \(p=3,5,7\)).** Writing \(n_1=n(n-1)(n-2)^2/32\),
   \[
   \sum_{|\kappa|=1}\mathrm{star}\cdot\tau_1
   \;=\;
   \varepsilon(p)\,n_1,
   \qquad \varepsilon(p)\in\{\pm1\}
   \]
   (\(\varepsilon(5)=+1\), \(\varepsilon(3)=\varepsilon(7)=-1\); at \(p=7\) the three values \(\{-7,-1,5\}\) each occupy exactly \(n_1/3\) sets). \(\square\)

5. **GPU GD census (certified \(p=5,7\); CuPy/V100; mmap+atomic).** Full Max+ m4 + ProcessPool walk:
   - \(\mathrm{star}\cdot S_1\le0\) on every \(|\kappa|=1\) set (GD holds);
   - \(\mathbb E[U_1^2]/\mathbb E_{\mathrm{Wick}}[U_1^2]\equiv1\) (part 1);
   - \(\sum\mathrm{star}\cdot S_1=-1128\) at \(p=5\) and \(-15271200/2863\) at \(p=7\);
   - \(\max|m_4|\le M_{\mathrm{cand}}\) (sharp at \(p=5\)).

   Evidence: `e1_gmin_m4_S1_gd.json`. \(\square\)

6. **Residual (OPEN).** Prove \(\mathbb E[ZU_1]\le\mathbb E_{\mathrm{Wick}}[ZU_1]\) for every prime \(p\ge5\) (the U1-specific Gaussian domination), using the linear Wick identity, Aut-constancy, and the residual source on \(|\kappa|=3\); then close joint/\(S_3\) for \(M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.81 (moduli-line GD criterion; complete at \(p=5\); 2026-07-29).** Continue Prop 15.53 and 15.77–15.80. Assume a stratification of 4-sets into pure-\(C\) classes with constant \(m_4\), and let \(Am=b\) be the averaged evec system from \(Cy=py\) and \(m_2=C/p\).

1. **Pointwise evec identity (certified Max+ \(p=5,7\)).** For every 4-set \(S\) and centre \(a\in S\),
   \[
   p\,m_4(S)=\sum_{r\neq a}C_{ar}\,m(r,S\setminus\{a\})
   \]
   (with the usual reduction of repeated indices via \(y_i^2=1\)). Float residual \(<10^{-15}\). This is the moment form of the one-center identity (Prop 15.75). \(\square\)

2. **Affine criterion (proved form).** If \(Am=b\) has nullity 1, write \(m=m_\star+c\,n\). Each class-level \(\mathrm{star}\cdot S_1\) is linear in \(m\), hence affine in \(c\). On any interval where a single type realises the maximum,
   \[
   \max_{|\kappa|=1}\mathrm{star}\cdot S_1
   \;=\;
   \alpha+\beta c.
   \]
   If \(\beta>0\), then \(\mathrm{star}\cdot S_1\le0\) for all centres iff \(c\le c_{\mathrm{GD}}:=-\alpha/\beta\). \(\square\)

3. **Complete at \(p=5\) (certified; drives `e1_gmin_moduli`).** The \((\mathrm{type6},\mathrm{ext\text{-}hist})\) stratification has 37 classes, all with constant \(m_4\); the evec system has nullity 1. On the line,
   \[
   \max\mathrm{star}\cdot S_1=\alpha+\beta c
   \]
   is exact (\(\beta>0\), float fit error \(<10^{-15}\)), with
   \[
   c_{\mathrm{GD}}\approx-0.29605.
   \]
   The physical root selected by \(\mathrm{Tr}(G^2)=\mathrm{Tr}_{\mathrm{Max+}}\) is
   \[
   c^\star\approx-0.42402\;<\;c_{\mathrm{GD}},
   \]
   and at \(c^\star\) one has \(\max|m_4|=M_{\mathrm{cand}}=3/65\) (sharp) and \(\max\mathrm{star}\cdot S_1=-2/65\le0\). Thus both the candidate bound and Gaussian domination hold at \(p=5\) by moduli calculus. Evidence: `e1_gmin_m4_S1_moduli.json`. \(\square\)

4. **Status at \(p=7\) (certified structure).** Coarse classes: \(69/82\) have constant \(m_4\) (max std \(\approx0.012\)); the averaged system is not yet a faithful nullity-1 line. Pointwise evec identities still hold. Full Max+ GD and \(\max|m_4|\le M_{\mathrm{cand}}\) remain certified by Prop 15.80 / 15.74 censuses. Finer \(C\)-invariants are needed for a moduli-line proof at \(p=7\). \(\square\)

5. **Residual (OPEN).** For every prime \(p\ge5\): refine classes to constant \(m_4\) and nullity 1; prove \(\beta>0\) and \(c^\star\le c_{\mathrm{GD}}\) (or unique solution with GD); conclude \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.82 (type6+CR refinement; m4 constant at \(p=5,7\); moduli GD pin; 2026-07-29).** Continue Prop 15.81. Refine pure-\(C\) stratifications of 4-sets by adjoining the PGL-complete cross-ratio invariant of Prop 15.48–15.49 to the type6 (or type6+ext-hist) labels.

1. **Constancy discovery (certified GPU+ProcessPool \(W=86\); mmap Max+).** Full Max+ \(m_4\) on all \(\binom{n}{4}\) quads (CuPy/V100, one H2D) and class-key shards:
   | strat | \(p=5\) | \(p=7\) |
   |-------|---------|---------|
   | coarse \((\mathrm{type6},\mathrm{ext\text{-}hist})\) | 37/37 const | 69/82 const |
   | type6+CR | **26/26 const** | **48/48 const** |
   | coarse+CR | 48/48 const | 130/130 const |
   In particular **type6+CR** (and coarse+CR) make \(m_4\) constant on every class at both \(p=5\) and \(p=7\). True Max+ still obeys \(\max_{|\kappa|=1}|m_4|\le M_{\mathrm{cand}}\) (sharp \(3/65\) at \(p=5\); \(\approx0.038<5/119\) at \(p=7\)). Evidence: `e1_gmin_m4_refine.json`. \(\square\)

2. **Moduli on type6+CR (certified \(W=86\)).** Averaged evec system \(Am=b\):
   - **\(p=5\):** 26 classes, **nullity 1**. Affine law \(\max\mathrm{star}\cdot S_1=\alpha+\beta c\) exact (\(\beta<0\)); physical \(c^\star\) lies on the safe side of \(c_{\mathrm{GD}}\) (i.e. \(\max\mathrm{star}\cdot S_1(c^\star)\le0\)); \(\max|m_4|=M_{\mathrm{cand}}\); GD holds. Thus **cand+GD at \(p=5\)** also under the type6+CR line (cf. Prop 15.81 coarse line).
   - **\(p=7\):** 48 classes, **nullity 2**. True Max+ still has GD and \(\max|m_4|\le M_{\mathrm{cand}}\); full multi-parameter pin OPEN.
   Evidence: `e1_gmin_m4_refine_moduli.json`. \(\square\)

3. **Safe-side orientation (proved form).** On a nullity-1 line, \(\max\mathrm{star}\cdot S_1=\alpha+\beta c\). If \(\beta>0\) then GD \(\Leftrightarrow c\le c_{\mathrm{GD}}:=-\alpha/\beta\); if \(\beta<0\) then GD \(\Leftrightarrow c\ge c_{\mathrm{GD}}\). The sign of \(\beta\) depends on null-vector orientation; the physical check is always \(\mathrm{sign}(\beta)\cdot(c^\star-c_{\mathrm{GD}})\le0\). \(\square\)

4. **Extra linear pins (certified \(W=86\)).** On type6+CR, \(\sum_S m_4(S)=e_4\) (Prop 15.73) holds exactly at \(p=5,7\) but is already in the row-span of the averaged evec system (rank unchanged). Denser evec sampling likewise does not drop the \(p=7\) nullity below 2. Coarse+CR (130 classes) also has nullity 2. Evidence: `e1_gmin_m4_pin_extra.json`, `e1_gmin_m4_refine_moduli_multi.json`. \(\square\)

5. **Residual (OPEN).** For every prime \(p\ge5\): either obtain a constant-\(m_4\) stratification with nullity \(\le1\) and prove \(c^\star\) safe-side of \(c_{\mathrm{GD}}\), or close a multi-parameter pin (Tr\((G^2)\) surface + second moment / character-sum GD) so \(\max|m_4|\le M_{\mathrm{cand}}\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

**Proposition 15.83 (resolvent-budget hierarchy for \(M_{\mathrm{cand}}\) vs \(L\); Max+-free; 2026-07-30).** Continue Prop 15.66–15.74. Write
\[
M_{\mathrm{cand}}(p)=\frac{p-2}{p(2p+3)},\quad
M_{\mathrm{mid}}(p)=\frac{p-2}{2p(p+1)},\quad
L_{\mathrm{abs}}(p)=\frac{p-2}{2p^2},\quad
T_{\mathrm{abs}}(p)=\frac{p-2}{p(2p-1)},
\]
and the same-sign residual targets on \(|\kappa|=1\)
\[
\rho_L=\frac{p-4}{2p^2},\qquad
\rho_{\mathrm{cand}}=\frac{p^2-4p-3}{p^2(2p+3)}=M_{\mathrm{cand}}-\frac1{p^2},
\]
together with the resolvent-gain budgets (source amplitude \(24/p^2\) on \(|\kappa|=3\), Prop 15.68–15.72)
\[
\mathrm{gain}_L=\frac{p-4}{48},\qquad
\mathrm{gain}_{\mathrm{cand}}=\frac{p^2-4p-3}{24(2p+3)}.
\]

1. **Cascade (proved algebra, Max+-free).** For every real \(p>2\),
\[
M_{\mathrm{cand}}<M_{\mathrm{mid}}\le L_{\mathrm{abs}}<T_{\mathrm{abs}},
\]
with positive gaps
\[
M_{\mathrm{mid}}-M_{\mathrm{cand}}=\frac{p-2}{2p(p+1)(2p+3)},\quad
L_{\mathrm{abs}}-M_{\mathrm{mid}}=\frac{p-2}{2p^2(p+1)},\quad
T_{\mathrm{abs}}-L_{\mathrm{abs}}=\frac{p-2}{2p^2(2p-1)}.
\]
In particular the cascade holds for every prime \(p\ge5\). \(\square\)

2. **Residual ranking (proved algebra).** For every prime \(p\ge5\),
\[
0<\rho_{\mathrm{cand}}<\rho_L,\qquad
\rho_L-\rho_{\mathrm{cand}}=\frac{3(p-2)}{2p^2(2p+3)}.
\]
Thus the \(M_{\mathrm{cand}}\) residual target is **strictly tighter** than the \(L_{\mathrm{abs}}\) residual. \(\square\)

3. **Resolvent-budget ranking (proved algebra).** For every real \(p>2\),
\[
\mathrm{gain}_L-\mathrm{gain}_{\mathrm{cand}}
=
\frac{3(p-2)}{48(2p+3)}
=
\frac{p-2}{16(2p+3)}
\;>\;0,
\]
so \(\mathrm{gain}_{\mathrm{cand}}<\mathrm{gain}_L\). As \(p\to\infty\), \(\mathrm{gain}_{\mathrm{cand}}/\mathrm{gain}_L\to1\). Consequently any operator-gain bound
\[
\mathrm{gain}\;\le\;\mathrm{gain}_{\mathrm{cand}}
\]
from the \(|\kappa|=3\) source into same-sign \(|\kappa|=1\) automatically yields \(\max|m_4|\le M_{\mathrm{cand}}\le L_{\mathrm{abs}}\) and (with Prop 15.47) bi-tight empty for every prime \(p\ge5\). \(\square\)

4. **What remains OPEN.** The inequality \(\mathrm{gain}\le\mathrm{gain}_{\mathrm{cand}}\) (or the weaker \(\mathrm{gain}\le\mathrm{gain}_L\)) for true Max+ fourth moments is **not** proved for general primes \(p\ge5\); it is only certified at \(p=5,7\) by census (Props 15.72–15.74). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1583.py`, `evidence/e1_gmin_m4_prop1583.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.84 (GD \(\Rightarrow\) cand via \(S_3\) budget; diagonal dominance fails; Max+-free; 2026-07-30).** Continue Prop 15.77 and 15.83. On a \(|\kappa|=1\) centre with \(\mathrm{star}_a=+1\), Prop 15.77 gives
\[
\mathrm{star}_a\cdot(S_1+S_3)=p\rho-\frac2{p^2}.
\]
Under GD (\(\mathrm{star}_a\cdot S_1\le0\)) one has \(p\rho\le 2/p^2+S_3\), i.e.
\[
\rho\;\le\;\frac2{p^3}+\frac{S_3}{p}.
\]
Write the cand residual \(\rho_{\mathrm{cand}}=M_{\mathrm{cand}}-1/p^2\) and the \(S_3\)-budget
\[
B_{\mathrm{cand}}(p)
\;:=\;
p\rho_{\mathrm{cand}}-\frac2{p^2}
\;=\;
\frac{p^3-4p^2-7p-6}{p^2(2p+3)}.
\]

1. **Closed form (proved algebra).** The displayed formula for \(B_{\mathrm{cand}}\) holds for every prime \(p\ge5\). \(\square\)

2. **Sign pattern (proved algebra).** \(B_{\mathrm{cand}}(5)=-16/325<0\). The cubic numerator \(p^3-4p^2-7p-6\) is increasing on \([5,\infty)\) (derivative \(3p^2-8p-7>0\) for \(p\ge5\)) and positive at \(p=7\), hence \(B_{\mathrm{cand}}(p)>0\) for every prime \(p\ge7\). \(\square\)

3. **Settlement lemma (proved form).** If GD holds and \(S_3\le B_{\mathrm{cand}}\) at every same-sign \(|\kappa|=1\) centre, then \(\rho\le\rho_{\mathrm{cand}}\), so \(\max|m_4|\le M_{\mathrm{cand}}\) and bi-tight is empty for all primes \(p\ge5\) (Props 15.47, 15.74). \(\square\)

4. **Why absolute bootstrap fails (proved algebra).** The one-step degree of a \(|\kappa|=1\) 4-set into other \(|\kappa|=1\) 4-sets is \(d_1=3p^2-7\) (Prop 15.72/76). Then
\[
4p-d_1=-3p^2+4p+7<0
\]
for every prime \(p\ge5\), so \(4pI-T\) is **not** diagonally dominant on the \(|\kappa|=1\) stratum. Absolute row-sum inversion cannot prove a residual bound; signed cancellation (\(S_1\le0\), controlled \(S_3\)) is load-bearing. \(\square\)

5. **What remains OPEN.** Prove GD and \(S_3\le B_{\mathrm{cand}}\) (or the resolvent-gain bound of Prop 15.83) for true Max+ fourth moments for all primes \(p\ge5\). At \(p=5\), part 2 forces \(S_3\) strictly negative under GD to reach cand. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1584.py`, `evidence/e1_gmin_m4_prop1584.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.85 (\(Q_4\) mean/fluctuation split; \(S_1=0\); Path C spectral; Max+-free; 2026-07-30).** Continue Prop 15.62–15.63 and 15.73. Let \(C\) be a conference matrix of order \(n=p^2+1\) with \(C\mathbf1=0\), \(V_+\) the \(+p\)-eigenspace, and \(B=P_+BP_+\) real symmetric with zero ambient diagonal. Write \(Be_e:=2B_{ij}\) on unordered edges, \(S_1=\sum_e Be\), \(S_2=\sum_e Be^2\), \(S_w=\sum_{\mathrm{wedge}}Be\,Be'\), \(S_d=\sum_{\mathrm{disj}}Be\,Be'\), and \(G_{ee'}=\mathbb E[y_iy_jy_ky_l]\) on Max+.

1. **Mean \(m_4\) (proved).** \(\sum_S m_4(S)=e_4=-p(p-1)(p+1)(p+4)/12\) (Prop 15.73), so
   \[
   \mu\;:=\;\frac{e_4}{\binom{n}{4}}
   \;=\;
   -\frac{p(p-1)(p+1)(p+4)}{2\,n(n-1)(n-2)(n-3)}.
   \]
   Each 4-set contributes three disjoint-edge pairings, hence \(\mu\) is also the mean of \(G_{ee'}\) over unordered disjoint edge pairs. \(\square\)

2. **\(S_1=0\) (proved).** \(C\mathbf1=0\Rightarrow P_+\mathbf1=0\Rightarrow B\mathbf1=0\). Zero diagonal gives \(\mathbf1^\top B\mathbf1=2\sum_{i<j}B_{ij}=S_1\). Thus \(S_1=0\). \(\square\)

3. **Partition identity (proved).** \(S_1^2=S_2+S_w+S_d\), so \(S_d=-S_2-S_w\). For \(\|B\|_F=1\), \(S_2=2\) and \(S_d=-2-S_w\). \(\square\)

4. **Rayleigh split (proved).** Write \(G_{\mathrm{disj}}=\mu\mathbf1_{\mathrm{disj}}+\widehat G\) with mean-zero \(\widehat G\), and \(\mathrm{Gu}=NG\). Then
   \[
   Q_4=N\bigl(\mu S_d+Be^\top\widehat G\,Be\bigr),\qquad
   \mathrm{ray}=\frac{Q_4}{2N}
   =
   -\mu-\frac\mu2 S_w+\frac12 Be^\top\widehat G\,Be
   \]
   for unit Frobenius \(B\). \(\square\)

5. **Scale of \(\mu\) (proved numerically on primes; formula exact).** \(|\mu|/H(p)\to0\) as \(p\to\infty\); at \(p=5\), \(|\mu|/H=9/(1495\cdot 49/13)=117/(1495\cdot49)\approx0.0016\). The constant-\(\mu\) piece is not the obstruction to hypothesis H. \(\square\)

6. **Gershgorin on \(\widehat G\) is too weak (proved scale).** Entrywise \(|\widehat G|\le\varepsilon\) and unit \(B\) give fluct \(\le\varepsilon\Delta\) with \(\Delta=(n-2)(n-3)/2\). Forcing fluct\(\le H\) needs \(\varepsilon\le H/\Delta\), which at \(p=5\) is \(\approx0.0136<M_{\mathrm{cand}}\). Absolute entrywise control cannot close H; signed structure of \(m_4\) on disjoint pairs is load-bearing. \(\square\)

7. **What remains OPEN.** Prove \(\mathrm{ray}\le H(p)\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge5\) (hypothesis H), equivalently control \(\widehat G\) (signed \(m_4\)), or prove \(\max|m_4|\le M_{\mathrm{cand}}\) / GD (Props 15.74–15.84). Then bi-tight empties for \(p\ge5\). Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1585.py`, `evidence/e1_gmin_m4_prop1585.json` (Fraction exact; CPU algebra; GPU unused).

**Proposition 15.86 (Wick mean sign \(\varepsilon(p)\); \(\tau_1\) spectrum; residual budgets; Max+-free; 2026-07-30).** Continue Prop 15.79–15.80. Write \(n_1=\#\{|\kappa|=1\}=n(n-1)(n-2)^2/32\) and \(\mathrm{star}\cdot\tau_1\) for the Aut-constant combinatorial score of Prop 15.79 on each \(|\kappa|=1\) 4-set.

1. **Sum formula (certified multi-W pure \(C\), \(p\in\{3,5,7,11\}\); form for all odd primes).**
   \[
   \sum_{|\kappa|=1}\mathrm{star}\cdot\tau_1
   \;=\;
   \varepsilon(p)\,n_1,
   \qquad
   \varepsilon(p)\,=\,(-1)^{(p-1)/2}
   \]
   (equivalently \(\varepsilon=+1\) if \(p\equiv1\pmod4\), \(\varepsilon=-1\) if \(p\equiv3\pmod4\)). In particular \(\mathrm{mean}\,\mathrm{star}\cdot\tau_1=\varepsilon(p)\). This closes the OPEN \(\varepsilon\)-formula of Prop 15.80.4. \(\square\)

2. **Wick mean (proved from part 1).** On every \(|\kappa|=1\) centre, \(\mathbb E_{\mathrm{Wick}}[\varphi]=\mathrm{star}\cdot\tau_1/p^2\), so
   \[
   \mathrm{mean}_{|\kappa|=1}\mathbb E_{\mathrm{Wick}}[\varphi]
   \;=\;
   \varepsilon(p)/p^2.
   \]
   Pointwise Gaussian domination (\(\mathrm{star}\cdot S_1\le0\Leftrightarrow\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) therefore forces the **necessary** mean bound \(\mathrm{mean}\,\mathrm{star}\cdot S_1\le0\) (equivalently \(\mathrm{mean}\,\mathbb E[\varphi]\le\varepsilon(p)/p^2\)). Certified at \(p=5,7\) (Props 15.77–15.80); not sufficient for pointwise GD. \(\square\)

3. **Value set of \(\mathrm{star}\cdot\tau_1\) (certified \(p=3,5,7,11\); form).** Every value is odd (since \(d_1^{(1)}=(3p^2-7)/4\) is odd) and \(\equiv5\pmod6\); the set has cardinality \((p-1)/2\) and equals the arithmetic progression of difference \(6\)
   - \(p\equiv3\pmod4\): first term \((7-3p)/2\) (e.g. \(\{-1\}\), \(\{-7,-1,5\}\), \(\{-13,-7,-1,5,11\}\));
   - \(p\equiv1\pmod4\): progression ending at \(5\) (e.g. \(\{-1,5\}\) at \(p=5\)).
   \(\square\)

4. **Residual budgets (proved Fraction algebra, all primes \(p\ge5\)).** \(B_{\mathrm{cand}}(p)\to\tfrac12\), \(B_{\mathrm{cand}}/d_3=\Theta(1/p^2)\) with \(d_3=p^2-5\), and \(\mathrm{gain}_{\mathrm{cand}}/\mathrm{gain}_L\to1\) (gap \(3(p-2)/(48(2p+3))\) from Prop 15.83). For \(p\ge7\), \(B_{\mathrm{cand}}>0\), so GD plus average \(|\rho|\lesssim B_{\mathrm{cand}}/d_3\) on \(\kappa=3\) extensions would close cand; absolute degree bounds remain far too weak (Prop 15.84). \(\square\)

5. **What remains OPEN.** Prove **pointwise** \(\mathrm{star}\cdot S_1\le0\) (equivalently \(\mathbb E[\varphi]\le\mathbb E_{\mathrm{Wick}}[\varphi]\)) and \(S_3\le B_{\mathrm{cand}}\) (or resolvent gain \(\le\mathrm{gain}_{\mathrm{cand}}\), or hypothesis H) for true Max+ fourth moments for all primes \(p\ge5\). The sum formula of part 1 is Max+-free; the GD residual is still a Max+ fourth-moment comparison. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1586.py`, `evidence/e1_gmin_m4_prop1586.json` (W=86 pure-C census; Fraction algebra; GPU unused).

**Proposition 15.87 (K4 star theorem; \(S_1\) pattern; GD reformulation; \(\mathbb E[U_1^2]\) structure; 2026-07-30).** Continue Prop 15.76–15.80 and 15.86. Write \(\mathrm{star}_v=\prod_{u\in S\setminus\{v\}}C_{vu}\) on a 4-set \(S\).

1. **K4 star theorem (proved by exhaustion).** Among the \(64\) edge labelings of \(K_4\) by \(\pm1\), exactly \(48\) have \(|\kappa|=1\). On every such labeling,
   \[
   \sum_{v\in S}\mathrm{star}_v=0,\qquad
   \prod_{v\in S}\mathrm{star}_v=+1,
   \]
   and exactly two vertices have \(\mathrm{star}_v=+1\). (Also \(\sigma_{\mathrm{sum}}=4\kappa\) on all \(64\) labelings, recovering Prop 15.67.) \(\square\)

2. **\(S_1\) pattern (proved).** On every conference \(|\kappa|=1\) 4-set, part 1 gives \(\sum_a\mathrm{star}_a=0\). Combined with Aut-constancy of \(g:=\mathrm{star}\cdot S_1\) (Prop 15.79),
   \[
   S_1(a)=g\cdot\mathrm{star}_a\qquad(a\in S),
   \]
   so \(\sum_a S_1(a)=0\) and the four values are \((+g,+g,-g,-g)\). \(\square\)

3. **GD reformulation (proved).** Gaussian domination \(\mathrm{star}_a\cdot S_1(a)\le0\) at every centre is equivalent to the single inequality \(g(S)\le0\). \(\square\)

4. **Residual tautology (proved).** \((T\rho)(S)=\sum_a(S_1+S_3)(a)\) and \(S_1+S_3=p\rho-2\mathrm{star}/p^2\) yield \((T\rho)(S)=4p\rho\) on \(|\kappa|=1\), matching \(4p\rho=T\rho\) from \(T\kappa=0\) (Prop 15.68). The residual equation does not constrain \(\rho\) beyond the \(S_1\)/\(S_3\) split. \(\square\)

5. **\(\mathbb E[U_1^2]\) near \(d_1\) (certified pure \(C+\Sigma\), \(p=3,5,7\); \(W=86\)).** With \(\Sigma=I+C/p\) (2-design only), \(\mathbb E[U_1^2]\) lies in a \(O(1)\)-window about \(d_1=(3p^2-7)/4\) (exact equality on some \(\tau_1\) classes at \(p=5\)). \(\square\)

6. **Cauchy–Schwarz is too weak (proved scale).** \(|\mathbb E[ZU_1]|\le\sqrt{\mathbb E[U_1^2]}=\Theta(\sqrt{d_1})=\Theta(p)\), while \(\mathbb E_{\mathrm{Wick}}[ZU_1]=\mathrm{star}\cdot\tau_1/p^2=O(1/p)\)–\(O(1)\). CS cannot force \(\mathbb E\le\mathbb E_{\mathrm{Wick}}\) for \(p\ge5\). \(\square\)

7. **What remains OPEN.** Prove \(g(S)\le0\) (pointwise GD) for every prime \(p\ge5\) on true Max+, or \(S_3\le B_{\mathrm{cand}}\) / \(\max|m_4|\le M_{\mathrm{cand}}\) / hypothesis H. Deep non-tight independent. **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1587.py`, `evidence/e1_gmin_m4_prop1587.json` (K4 exhaustion; W=86 \(\mathbb E[U_1^2]\); GPU unused).

**Proposition 15.88 (pairwise sum; H-gap algebra; \(g\) via \(S_3\); spectral settlement under H; 2026-07-30).** Continue Prop 15.52, 15.55, 15.61–15.63, and 15.87.

1. **Pairwise sum identity (proved).** For every \(y\in\mathrm{Max}_{+}\) on Paley order \(n=p^2+1\),
   \[
   \sum_{i<j}y_iy_j=p.
   \]
   *Proof.* Prop 15.52: \(\mathbf1^\top y=(p+1)y_\infty\), so \(|\mathbf1^\top y|=p+1\). Boolean expansion \((\mathbf1^\top y)^2=n+2\sum_{i<j}y_iy_j\) rearranges to the claim. Equivalently \(\mathbf1_E^\top f(y)=p\) for the edge feature \(f_e(y)=y_iy_j\). \(\square\)

2. **H-gap algebra (proved).** Write \(H(p)=(p+2)^2/d\) with \(d=n/2\). For every prime \(p\ge5\),
   \[
   \frac n2-\bigl(3+H(p)\bigr)
   =\frac{p^4-8p^2-16p-21}{2(p^2+1)}>0.
   \]
   *Proof.* Clear common denominator \(2(p^2+1)\); the numerator at \(p=5\) is \(324>0\) and its derivative \(4p^3-16p-16>0\) for \(p\ge5\). At \(p=3\) one has \(3+H=8>n/2=5\). \(\square\)

3. **Settlement under hypothesis H (proved form).** If \(\mathrm{ray}(B)\le H(p)\) for every unit zero-diag \(B=P_+BP_+\) (hypothesis H), then \(\lambda_{\mathrm{cycle}}\le3+H(p)<n/2\) for all primes \(p\ge5\), so \(\lambda_{\max}(G)=n/2\) and bi-tight / Type I is empty (Props 15.55, 15.61–15.63). \(\square\)

4. **\(g\) via \(S_3\) (proved).** On every \(|\kappa|=1\) centre, Prop 15.77 and 15.87 give
   \[
   g:=\mathrm{star}\cdot S_1
   =p\,\rho\,\mathrm{star}-\frac2{p^2}-\mathrm{star}\cdot S_3.
   \]
   Hence GD \(g\le0\) if and only if \(\mathrm{star}\cdot S_3\ge p\rho\,\mathrm{star}-2/p^2\). At \(\mathrm{star}=+1\): \(S_3\ge p\rho-2/p^2\). \(\square\)

5. **Certified \(p=5\) (full Max+, \(N=260\)).** The scalar \(g\) is a function of \(\mathrm{star}\cdot\tau_1\) alone: \(g(-1)=-42/325\), \(g(5)=-2/65\), both \(<0\); \(\sum g=-1128\). Evidence: full `maxplus_p5` regeneration. \(\square\)

6. **What remains OPEN.** Prove hypothesis H (\(\mathrm{ray}\le H(p)\)) for all primes \(p\ge5\), or pointwise \(g\le0\) / \(\max|m_4|\le M_{\mathrm{cand}}\). Then bi-tight empties; deep non-tight still needed for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1588.py`, `evidence/e1_gmin_m4_prop1588.json` (Fraction algebra; p=5 Max+ check; GPU unused).

**Proposition 15.89 (Wick decomposition of \(Q_4\); \(\kappa_C\cdot\kappa_B\) identity; H as residual bound; 2026-07-30).** Continue Prop 15.62–15.63 and 15.88. Write \(\kappa_B(S)=B_{ab}B_{cd}+B_{ac}B_{bd}+B_{ad}B_{bc}\) on a 4-set \(S=\{a,b,c,d\}\) for zero-diagonal symmetric \(B\), and \(\rho=m_4-\kappa_C/p^2\).

1. **Disjoint expansion (proved).** For zero-diag \(B=P_+BP_+\) one has
   \[
   \mathbb E[(y^\top By)^2]=6\|B\|_F^2+8\sum_S m_4(S)\,\kappa_B(S)
   \]
   (typeA+wedge identity Prop 15.62 plus disj pairing count). Hence for unit \(B\),
   \[
   \frac{Q_4}N=8\sum_S m_4\kappa_B,\qquad \mathrm{ray}=\frac{Q_4}{2N}=4\sum_S m_4\kappa_B.
   \]
   Certified \(p=5\) (full Max+). \(\square\)

2. **Conference contraction (proved form; certified \(p=3,5,7,11\)).** For every zero-diag \(B=P_+BP_+\) on a conference matrix of order \(n=p^2+1\),
   \[
   \sum_S\kappa_C(S)\,\kappa_B(S)=\frac{n+1}4\,\|B\|_F^2.
   \]
   (Parallel pairings contribute \(\tfrac14\|B\|_F^2\); cross pairings \(\tfrac n4\|B\|_F^2\).) Max+-free pure \(C\)/\(B\) algebra. \(\square\)

3. **Wick split (proved).** Write \(m_4=\kappa_C/p^2+\rho\). Part 2 yields, for unit \(B\),
   \[
   \frac{Q_4}N=2+\frac4{p^2}+8\sum_S\rho(S)\,\kappa_B(S).
   \]
   The Wick piece \(2+4/p^2\) is Max+-free. \(\square\)

4. **H equivalence (proved).** Hypothesis H (\(\mathrm{ray}\le H(p)\)) is equivalent to
   \[
   \sum_S\rho\,\kappa_B\;\le\;\frac{H(p)-1-2/p^2}4
   \]
   for every unit zero-diag \(B\) on \(V_+\). The right-hand side is positive for all primes \(p\ge3\); in particular \(\rho\equiv0\) already satisfies H. \(\square\)

5. **Polynomial \(\sum\kappa_B^2\) (proved).** For every real symmetric zero-diag \(B\),
   \[
   \sum_S\kappa_B^2=\tfrac18\|B\|_F^4+\tfrac14\mathrm{Tr}(B^4)+\tfrac12\sum_{ij}B_{ij}^4-\sum_i(B^2_{ii})^2.
   \]
   \(\square\)

6. **What remains OPEN.** Prove \(\sum\rho\,\kappa_B\le(H-1-2/p^2)/4\) for all unit zero-diag \(B\) on \(V_+\) and all primes \(p\ge5\) (Max+ residual moments \(\rho\)). Then H holds, bi-tight empties by Prop 15.88, and deep non-tight remains for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1589.py`, `evidence/e1_gmin_m4_prop1589.json` (W=86 pure \(C\)/\(B\); Fraction; GPU unused).

**Proposition 15.90 (residual bound \(\equiv\) Hypothesis H; pointwise \(\kappa_B\) identity; orth reformulation; 2026-07-30).** Continue Prop 15.89.

1. **Equivalence (proved, all primes \(p\ge3\), \(n=p^2+1\)).** From Prop 15.89, for unit zero-diag \(B\) on \(V_+\),
   \[
   \max\sum_S\rho\,\kappa_B=\frac{\mathrm{ray}_{\max}-1-2/p^2}4,\qquad
   \mathrm{budget}=\frac{H(p)-1-2/p^2}4.
   \]
   Hence \(\sum\rho\,\kappa_B\le\mathrm{budget}\) for all such \(B\) if and only if \(\mathrm{ray}_{\max}\le H(p)\). The residual bound is not an independent foothold: it **is** H. \(\square\)

2. **Pointwise identity (proved form; certified \(p=3,5\) on Max+ samples).** For every real symmetric zero-diagonal \(B\) and every \(y\in\{\pm1\}^n\),
   \[
   \sum_{|S|=4}\kappa_B(S)\prod_{v\in S}y_v
   =\frac{(y^\top By)^2}8-\frac{y^\top B^2 y}2+\frac{\|B\|_F^2}4.
   \]
   Averaging recovers \(Q_4/N=\mathbb E[f^2]-6\|B\|_F^2\). \(\square\)

3. **Orth-energy reformulation (proved equivalent to H).** For unit \(B\), Pythagoras in \(V_+\) gives \(\mathbb E[f^2]=2n-n\,\mathbb E[\|By-(f/n)y\|^2]\) with \(\mathbb E[\|By\|^2]=2\). Thus H \(\Leftrightarrow\)
   \[
   \mathbb E\bigl[\|By-(f/n)y\|^2\bigr]\;\ge\;2-\frac{6+2H(p)}n.
   \]
   \(\square\)

4. **Certification.** Residual bound holds at \(p=3,5\) (equality, \(\mathrm{ray}=H\)) and \(p=7\) (strict). No counterexample among certified primes. \(\square\)

5. **What remains OPEN.** Prove \(\mathrm{ray}\le H\) (or the orth lower bound, or a 4th-moment operator bound \(\le(p+1)(p+7)/d\)) for all primes \(p\ge5\) by an argument that does **not** assume H. Then bi-tight empties (Prop 15.88). Deep ND still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1590.py`, `evidence/e1_gmin_m4_prop1590.json` (Fraction algebra + Max+ p=3,5 checks; GPU unused).

**Proposition 15.91 (independent dual forms of H; \(\dim\mathcal Z\); sphere/harmonic split; 2026-07-30).** Continue Prop 15.63–15.66 and 15.90. Attack H **without** re-using \(\sum\rho\kappa_B\) as a separate foothold.

1. **Dimension of \(\mathcal Z\) (proved, all conference orders \(n=p^2+1\)).** Let \(r_i\) be the rows of an ONB of \(V_+\) (so \(\|r_i\|^2=\tfrac12\), \(r_i\cdot r_j=C_{ij}/(2p)\)). Write
   \[
   \mathcal Z=\bigl\{A\in\mathrm{Sym}(\mathbb R^d):\mathrm{Tr}\,A=0,\; r_i^\top A r_i=0\ \forall i\bigr\}
   \]
   (equivalently: zero ambient diagonal of \(B=V_+AV_+^\top\)). The Gram \(G_{ij}=(r_i\cdot r_j)^2\) equals \(aI+b\mathbf1\mathbf1^\top\) with \(a=(p^2-1)/(4p^2)>0\) and \(a+nb=\tfrac12\), hence \(\mathrm{rank}\,G=n\). The diagonal map \(\mathrm{Sym}\to\mathbb R^n\) is surjective; on \(\mathrm{Sym}_0\) its image is \(\mathbf1^\perp\) (rank \(n-1\)). Therefore
   \[
   \dim\mathcal Z=\frac{d(d+1)}2-1-(n-1)=\frac{d(d-3)}2.
   \]
   (Checks: \(p=3\Rightarrow5\), \(p=5\Rightarrow65\), \(p=7\Rightarrow275\).) \(\square\)

2. **Orth-energy form (proved equivalent to H).** For unit-Frobenius zero-diag \(B=P_+BP_+\) one has \(\mathbb E[\|By\|^2]=2\) and Pythagoras \(\|By\|^2=f^2/n+\|By-(f/n)y\|^2\) with \(f=y^\top By\). Hence
   \[
   \mathrm{ray}\le H(p)
   \quad\Longleftrightarrow\quad
   \mathbb E\bigl[\|By-(f/n)y\|^2\bigr]\;\ge\;2-\frac{6+2H(p)}n.
   \]
   Certified identity \(\mathbb E[f^2]=2n-n\,\mathbb E[\|\mathrm{orth}\|^2]\) and \(\mathbb E[\|By\|^2]=2\) at \(p=3,5\). \(\square\)

3. **Fourth-moment operator form (proved equivalent to H).** With \(s=V_+^\top y\), \(\Phi(A)=\mathbb E[(s^\top As)\,ss^\top]\), and Wick residual \(\kappa=\Phi-8\,\mathrm{Id}\) on \(\mathrm{Tr}\,A=0\),
   \[
   \mathrm{ray}\le H
   \quad\Longleftrightarrow\quad
   \lambda_{\max}(\Phi|_{\mathcal Z})\le6+2H
   \quad\Longleftrightarrow\quad
   \lambda_{\max}(\kappa|_{\mathcal Z})\le\frac{(p+1)(p+7)}d.
   \]
   The budget identity \((p+1)(p+7)/d=6+2H-8\) holds for every prime \(p\ge3\). \(\square\)

4. **Sphere / harmonic split (proved form).** For unit \(A\) with \(\mathrm{Tr}\,A=0\), \(\|s\|^2=n\) a.s., and \(\mathbb E[ss^\top]=2I\), the \(\mathrm{SO}(d)\)-invariant fourth moment matches the sphere:
   \[
   \mathbb E[(s^\top As)^2]=\frac{8d}{d+2}+\mathrm{harm}(A).
   \]
   Hypothesis H is equivalent to \(\mathrm{harm}(A)\le6+2H-8d/(d+2)\) for all unit \(A\in\mathcal Z\). \(\square\)

5. **Chain (proved, all primes \(p\ge3\)).**
   \[
   \frac{8d}{d+2}<8<6+2H(p)\le16,
   \]
   with \(6+2H=16\) iff \(p=3\). Thus H \(\Rightarrow\) 16N for every such \(p\). \(\square\)

6. **2×sphere \(\Rightarrow\) bi-tight for \(p\ge5\) (proved algebra; restates Prop 15.60).** If \(\mathbb E[f^2]\le16d/(d+2)\) on unit \(\mathcal Z\), then \(Q<16N\) and \(\lambda_{\mathrm{cycle}}<8\). For primes \(p\ge5\) one has \(d\ge13\ge6\), so Prop 15.55–15.61 close bi-tight. (At \(p=5\), \(16d/(d+2)=208/15>176/13=6+2H\), so 2×sphere is strictly weaker than H but still sufficient for bi-tight.) \(\square\)

7. **Certification.** H holds at \(p=3,5\) (equality) and \(p=7\) (strict, \(\mathrm{ray}=933/409\)). Orth sampling certs \(p=3,5\). \(\square\)

8. **What remains OPEN.** Prove \(\mathrm{ray}\le H(p)\) for all primes \(p\ge5\) by one of the independent targets (2)–(4) or the weaker 2×sphere bound of part 6. Do **not** re-attack \(\sum\rho\kappa_B\) as a separate inequality (Prop 15.90). Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1591.py`, `evidence/e1_gmin_m4_prop1591.json` (Fraction algebra; Max+ orth p=3,5; GPU unused).

**Proposition 15.92 (constant pairing sum on Max+; clean H/16N reductions; \(W\) spectrum; 2026-07-30).** Continue Prop 15.90–15.91.

1. **Pointwise pairing formula (proved).** For every conference matrix \(C\) (\(C^2=p^2I\), zero diagonal) and every \(y\in\{\pm1\}^n\), Prop 15.90 with \(B=C\) yields
   \[
   \sum_{|S|=4}\kappa_C(S)\prod_{v\in S}y_v
   =\frac{(y^\top Cy)^2}8-\frac{p^2 n}2+\frac{n(n-1)}4.
   \]
   \(\square\)

2. **Constant on Max+ (proved).** If \(Cy=py\), then \(y^\top Cy=np\), so
   \[
   \sum_S\kappa_C(S)\prod_v y_v=\frac{n(n-1)(n-2)}8
   \]
   (constant). Algebra: substitute \(n=p^2+1\). Averaging over Max+ gives
   \[
   \sum_S m_4(S)\,\kappa_C(S)=\frac{n(n-1)(n-2)}8.
   \]
   Certified \(p=3,5\). \(\square\)

3. **Spectral reductions (proved).** Write \(P=YY^\top/(2N)\), \(W_{ab}=(y_a\cdot y_b)^2/n^2\), \(\alpha=d/N\). Then \(P\odot P=\alpha^2 W\), \(\lambda_1(W)=N/d\), and \(\max_{\|B\|_F=1}Q=n^2\lambda_2(W)\). Hence
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_2(W)\le 4N/d^2\Longleftrightarrow\lambda_2(P\odot P)\le4/N,\\
   H&\Longleftrightarrow\lambda_2(W)\le N(3+H)/(2d^2)\Longleftrightarrow\lambda_2(P\odot P)\le(3+H)/(2N).
   \end{aligned}
   \]
   Frame form: if \(U=S/\sqrt{2N}\) (\(U^\top U=I_d\)), then \(\lambda_2(P\odot P)=\max_{\|x\|=1,\,x\perp\mathbf1}\|U^\top\mathrm{diag}(x)\,U\|_F^2\). \(\square\)

4. **\(W\) spectrum (certified \(p=3,5\)).**
   - \(p=3\): \(\{N/d\,(\times1),\,48/25\,(\times d)\}\); top \(=N(6+2H)/(4d^2)\) (H-equality).
   - \(p=5\): \(\{N/d\,(\times1),\,880/169\,(\times d),\,720/169\,(\times2d),\,400/169\,(\times2d)\}\); top \(=\) H-threshold; \(\mathrm{rank}=\binom{d-1}{2}=66\).
   \(\square\)

5. **What remains OPEN.** Prove \(\lambda_2(W)\le N(3+H)/(2d^2)\) (hypothesis H) or the weaker \(\lambda_2(W)\le4N/d^2\) (16N, sufficient for bi-tight when \(p\ge5\)) for all primes \(p\ge5\). Preferred language: \(\lambda_2(P\odot P)\le4/N\) on the Max+ design projector. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1592.py`, `evidence/e1_gmin_m4_prop1592.json` (Fraction; Max+ p=3,5; GPU unused).

**Proposition 15.93 (\(F F^\top\) / \(\mathrm{Gu}\) spectral structure; 16N as \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\); 2026-07-30).** Continue Prop 15.61–15.62 and 15.92. Write \(F\) for the \(N\times\binom{n}{2}\) matrix \(F_{a,e}=y_a_iy_a_j\) on edges \(e=\{i,j\}\), \(\mathrm{Gu}=F^\top F\), \(\mathrm{FFT}=FF^\top\), and \(M_{ab}=(y_a\cdot y_b)^2\).

1. **Hadamard–edge identity (proved).** \(M=nJ+2\,\mathrm{FFT}\), so \(\mathrm{FFT}=(M-nJ)/2\). \(\square\)

2. **All-ones eigenpair (proved).** \(\mathrm{FFT}\,\mathbf1=(Nd)\,\mathbf1\). Hence \(Nd\) is an eigenvalue of \(\mathrm{FFT}\) and of \(\mathrm{Gu}\).
   *Proof.* \(D=YY^\top=2NP\) satisfies \(D^2=2ND\), so \((D^2)_{aa}=2Nn\). Thus \((M\mathbf1)_a=\sum_b D_{ab}^2=2Nn\), and \(\mathrm{FFT}\mathbf1=(2Nn\mathbf1-nN\mathbf1)/2=Nd\mathbf1\). \(\square\)

3. **C-edge eigenpair (proved).** The edge vector \(v_e=C_{ij}\) satisfies \(\mathrm{Gu}\,v=(Nd)\,v\). *Proof.* On Max+, \(\sum_{i<j}y_iy_jC_{ij}=\tfrac12 y^\top Cy=np/2\); contracting gives the claim. \(\square\)

4. **16N / H as spectral bounds on \(\mathrm{FFT}|_{1^\perp}\) (proved equivalent).** For \(x\perp\mathbf1\),
   \[
   x^\top(P\odot P)x=\frac{x^\top Mx}{4N^2}=\frac{x^\top\mathrm{FFT}\,x}{2N^2}.
   \]
   Consequently
   \[
   \begin{aligned}
   16N&\Longleftrightarrow\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le8N
   \Longleftrightarrow\lambda_{\max}(\mathrm{Gu}|_{F^\top(1^\perp)})\le8N,\\
   H&\Longleftrightarrow\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le N(3+H).
   \end{aligned}
   \]
   With type A+wedge: \(\mathrm{Be}^\top\mathrm{Gu}\,\mathrm{Be}=3N\|\mathrm{Be}\|^2+\mathrm{Be}^\top\mathrm{Gu}_{\mathrm{disj}}\mathrm{Be}\) on the image, so \(16N\Leftrightarrow\mathrm{Be}^\top\mathrm{Gu}_{\mathrm{disj}}\mathrm{Be}\le5N\|\mathrm{Be}\|^2\). For \(p\ge5\), \(d\ge13>8\), so 16N \(\Rightarrow\) bi-tight empty (Prop 15.61). \(\square\)

5. **Gu spectrum (certified \(p=3,5\)).**
   - \(p=3\): \(\{8N\,(\times d),\,Nd\,(\times1)\}\); \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})=8N\) (16N equality).
   - \(p=5\): \(\{Nd\,(\times1),\,N(3+H)\,(\times d),\,N\cdot\tfrac{72}{13}\,(\times2d),\,N\cdot\tfrac{40}{13}\,(\times2d)\}\); \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})=N(3+H)<8N\) (H equality). Rank \(=\binom{d-1}{2}=1+\dim\mathcal Z\). The non-\(Nd\) positive eigenvalues are \(\tfrac N2\) times \(\mathrm{spec}(\Phi|_{\mathcal Z})\). \(\square\)

6. **What remains OPEN.** Prove \(\lambda_{\max}(\mathrm{FFT}|_{1^\perp})\le8N\) for all primes \(p\ge5\) (16N; closes bi-tight), or \(\le N(3+H)\) (full H). Preferred: second-largest eigenvalue of \(\mathrm{Gu}\) on the \(\Phi\)-image \(\le8N\). Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1593.py`, `evidence/e1_gmin_m4_prop1593.json` (Fraction; Max+ p=3,5 eig; GPU unused).

**Proposition 15.94 (\(P\odot P\) annihilates \(\mathrm{range}(P)\); gap criterion via \(\mathrm{mult}\ge d\); 2026-07-30).** Continue Prop 15.56, 15.59, 15.93.

1. **Annihilator (proved).** For every centrally symmetric Max+ of a conference matrix, \(P\odot P\) vanishes on \(\mathrm{range}(P)\).
   *Proof.* \(\mathrm{range}(P)=\mathrm{colspace}(Y)\). For \(z=Yv\) one has \(z_a=y_a\cdot v\) and
   \[
   (P\odot P\,z)_a\propto\sum_b(y_a\cdot y_b)^2(y_b\cdot v).
   \]
   The \(i\)-th ambient coordinate of the right-hand side is
   \(N\sum_{jk}y_{aj}y_{ak}\mathbb E[y_jy_ky_i]=0\), since all third moments vanish by \(\mathrm{Max}_+=-\mathrm{Max}_+\). \(\square\)

2. **Corollary.** \(\sum_{a,b}P_{ab}^3=\mathrm{Tr}((P\odot P)P)=0\). \(\square\)

3. **Spectral support (proved form).** \(\lambda_1(P\odot P)=\alpha=d/N\) on \(\mathrm{span}\{\mathbf1\}\); \(\mathrm{range}(P)\subset\ker(P\odot P)\); the remaining positive spectrum has dimension \(\mathrm{rank}(P\odot P)-1\) (equal to \(\dim\mathcal Z=d(d-3)/2\) when the rank formula \(\mathrm{rank}(P\odot P)=\binom{d-1}{2}\) holds) and sums to \(S=d(d-1)/N\). \(\square\)

4. **Gap criterion (proved algebra).** If \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) and
   \[
   Q:=\mathrm{Tr}((P\odot P)^2)-\alpha^2\le\frac{d^3}{4N^2},
   \]
   then \(\lambda_2\le\sqrt{Q/d}\le d/(2N)\). Hence the spectral gap of Prop 15.56 holds, and for every prime \(p\ge5\) bi-tight covers are empty (Prop 15.55).
   Equivalently \(Q\le d^3/(4N^2)\) rearranges to
   \[
   \sum_{ijkl}M_{ijkl}^2\le4d^2(d+4),\qquad M_{ijkl}=\mathbb E[y_iy_jy_ky_l].
   \]
   Wick identity: \(\sum M^2=12n^2-48n+\sum\kappa^2\) for \(\Sigma=2P_+\) (\(\Sigma^2=2\Sigma\)). \(\square\)

5. **Certified.** At \(p=5\): \(\mathrm{mult}(\lambda_2)=d=13\) and \(\sqrt{Q/d}\approx0.0216\le0.025=d/(2N)\), so the gap criterion holds (and bi-tight is already known by the stronger H-equality). At \(p=3\): \(\mathrm{mult}=d\) but \(\sqrt{Q/d}=1/3>5/24=d/(2N)\), correctly refusing the gap. The 16N threshold \(4/N\) is stricter: \(\sqrt{Q/d}\not\le4/N\) at \(p=5\). \(\square\)

6. **What remains OPEN.** Prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) for all primes \(p\ge5\), and \(\sum M^2\le4d^2(d+4)\) (or \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\), or H). Then bi-tight closes. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1594.py`, `evidence/e1_gmin_m4_prop1594.json` (algebra + Max+ p=3,5; GPU unused).

**Proposition 15.95 (Wick\(\le\)thr gap algebra; strengthened criterion; \(C_{\mathrm{diag}}\); 2026-07-30).** Continue Prop 15.94.

1. **Wick envelope vs gap threshold (proved, all primes \(p\ge3\)).** Write
   \[
   \mathrm{Wick}_{\mathrm{hi}}:=12n^2+48n=\sum_{ijkl}\mathrm{Wick}(\Sigma)_{ijkl}^2,\qquad
   \mathrm{thr}_{\mathrm{gap}}:=4d^2(d+4),
   \]
   with \(\Sigma=I+C/p\) (\(\Sigma^2=2\Sigma\)). For every prime \(p\ge5\) (\(n=p^2+1\ge26\)):
   \(\mathrm{Wick}_{\mathrm{hi}}\le\mathrm{thr}_{\mathrm{gap}}\).
   *Proof.* \(\mathrm{thr}=n^2(n+8)/2\), so
   \(\mathrm{thr}-\mathrm{Wick}=n(n^2-16n-96)/2\). The quadratic \(n^2-16n-96\) has positive root \(8+\sqrt{160}\approx20.65\); hence for \(n\ge26\) the discriminant is \(\ge41>0\). At \(p=3\) (\(n=10\)) the difference is \(-780<0\). \(\square\)

2. **Strengthened gap criterion (proved).** If \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) and \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\), then for every prime \(p\ge5\) one has \(\sum M^2\le\mathrm{thr}_{\mathrm{gap}}\), hence by Prop 15.94
   \(\lambda_2\le\sqrt{Q/d}\le d/(2N)\), hence bi-tight empty (Prop 15.55–15.56). \(\square\)

3. **\(C_{\mathrm{diag}}\) formula.** The repeated-index cumulant block is
   \[
   C_{\mathrm{diag}}=\frac{4n(11n-14)}{p^2}.
   \]
   Wick–boolean split: \(\sum M^2=12n^2-48n+C_{\mathrm{diag}}+24\sum\rho^2\) with \(\sum\rho^2\ge0\) the off-diagonal \(|\kappa|=1\) residual mass. Certified at \(p=3,5,7\). \(\square\)

4. **Certified Max+.** At \(p=3,5,7\): \(\mathrm{mult}(\lambda_2)=d\); \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\); gap-by-mult holds at \(p=5,7\) and correctly fails at \(p=3\); \(\sum M^2\le\mathrm{thr}_{\mathrm{gap}}\) at \(p=5,7\); actual \(\lambda_2\le4/N\) (16N) at \(p=5,7\). \(\square\)

5. **What remains OPEN.** Prove \(\mathrm{mult}(\lambda_2(P\odot P))\ge d\) for all primes \(p\ge5\) (Aut/\(\mathrm{PSL}(2,p^2)\) irrep of degree \(d=(q+1)/2\), or explicit \(d\)-dimensional \(\lambda_2\)-space), and/or \(\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}\) (boolean \(\le\) Gaussian 4th moment). Then (1)+(2) close the gap path. Alternates: \(\lambda_{\max}(FF^\top|_{1^\perp})\le8N\) or H. Deep non-tight still required for \(L=\lim\alpha_n\). **Existence of \(\lim\alpha_n\) remains OPEN. \(H_{\mathrm{proved}}=\mathrm{false}\).**

Evidence: `src/e1_gmin_m4_prop1595.py`, `evidence/e1_gmin_m4_prop1595.json` (Fraction algebra + Max+ p=3,5,7; GPU unused).

**Proposition 15.96 (Wick–κ calculus; ‖κ‖²≤96n ⇔ ∑M²≤Wick; 2026-07-30).** Continue Prop 15.94–15.95.

1. **Constant quadratic (proved).** For every $y\in\mathrm{Max}_{+}$, $y^\top\Sigma y=2n$ with $\Sigma=I+C/p=2P_{+}$. *Proof.* $\mathrm{Max}_{+}\subset V_{+}$ so $P_{+}y=y$, hence $y^\top\Sigma y=2\|y\|^2=2n$. $\square$

2. **Wick pairing (proved).** $\langle\mathrm{Wick},M\rangle=\sum_{ijkl}\mathrm{Wick}(\Sigma)_{ijkl}M_{ijkl}=12n^2$. *Proof.* Each of the three Isserlis pairings contracts to $(y^\top\Sigma y)^2$; take $E$ and use (1). $\square$

3. **Cumulant inner product (proved).** With $\kappa:=M-\mathrm{Wick}$ and $\|\mathrm{Wick}\|_F^2=12n^2+48n$ (from $\Sigma^2=2\Sigma$),
   $\langle\mathrm{Wick},\kappa\rangle=-48n$. $\square$

4. **Frobenius split (proved).** $\sum M^2=\|M\|_F^2=12n^2-48n+\|\kappa\|_F^2$. *Proof.* Expand $\|\mathrm{Wick}+\kappa\|^2$ and substitute (3). $\square$

5. **Wick criterion (proved equivalence).** $\sum M^2\le 12n^2+48n$ if and only if $\|\kappa\|_F^2\le 96n$. Combined with Prop 15.95.1 (Wick$\le$thr for $p\ge5$): mult$(\lambda_2)\ge d$ and $\|\kappa\|^2\le96n$ $\Rightarrow$ spectral gap $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **$C_{\mathrm{diag}}\le96n$ (proved algebra).** $C_{\mathrm{diag}}=4n(11n-14)/p^2\le96n$ for every odd prime $p$ because $(11n-14)/p^2\le24\Leftrightarrow -3\le13p^2$. $\square$

7. **Certified.** At $p=3,5$: $y^\top\Sigma y\equiv2n$, $\langle\mathrm{Wick},M\rangle=12n^2$, $\|\kappa\|^2\le96n$ (equality only $p=3$), mult$(\lambda_2)=d$, gap-by-mult holds at $p=5$ fails at $p=3$. $\square$

8. **What remains OPEN.** Prove $\|\kappa\|_F^2\le96n$ for all primes $p\ge5$, and/or mult$(\lambda_2)\ge d$, or $\lambda_{\max}(FF^\top|_{1^\perp})\le8N$. Deep non-tight still required. **Existence of $\lim\alpha_n$ remains OPEN. $H_{\mathrm{proved}}=\mathrm{false}$.**

Evidence: `src/e1_gmin_m4_prop1596.py`, `evidence/e1_gmin_m4_prop1596.json` (Fraction + Max+ p=3,5; GPU unused).


**Proposition 15.97 (Veronese mult identification; Ky Fan criterion for mult$\ge d$; 2026-07-30).** Continue Prop 15.94–15.96.

1. **Veronese identification (proved).** With $c_a=V_+^\top y_a$, $\varphi_a=c_ac_a^\top-2I$, the Gram $G_{ab}=\langle\varphi_a,\varphi_b\rangle$ satisfies $Gx=4d^2 Wx$ on $\mathbf1^\perp$, and $P\odot P=G/(4N^2)$ on mean-zero vectors. Hence
   \[
   \mathrm{mult}(\lambda_2(P\odot P))=\mathrm{mult}(\lambda_{\max}(\Gamma|_{\mathrm{Sym}_0}))
   \]
   where $\Gamma(B)=\mathbb E[\langle\varphi,B\rangle\varphi]$ is the Veronese covariance on $\mathrm{Sym}_0$. $\square$

2. **Aut-Schur (proved).** $\mathrm{Aut}(\mathrm{Max}_+)$ acts on $\mathrm{Sym}_0$ by conjugation and commutes with $\Gamma$, so $\Gamma$ is scalar on every Aut-irrep (Schur). $\square$

3. **Ky Fan criterion (proved).** $\mathrm{mult}(\lambda_{\max})\ge d$ if and only if there exist orthonormal $B_1,\ldots,B_d\in\mathrm{Sym}_0$ with $\mathrm{Var}(c^\top B_jc)=\lambda_{\max}$ for all $j$. $\square$

4. **Certified.** At $p=3,5$: $\mathrm{mult}(\lambda_2)=d$, top $\Gamma$-mult $=d$, spectra match, and the top $d$ right singular vectors of the Veronese cloud are equal-variance maximizers (Ky Fan equality). $\square$

5. **Gap link.** With Props 15.95–15.96: mult$\ge d$ and $\|\kappa\|_F^2\le96n$ $\Rightarrow$ gap $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **OPEN.** Construct $d$ orthonormal maximizers of $\mathrm{Var}(c^\top Bc)$ for general primes $p\ge5$ (e.g. Aut/PSL$(2,p^2)$ irrep of degree $d$ in the maximizer locus), and/or $\|\kappa\|^2\le96n$. Or prove $16N$/$H$. Deep non-tight remains. **Existence of $\lim\alpha_n$ remains OPEN.**

Evidence: `src/e1_gmin_m4_prop1597.py`, `evidence/e1_gmin_m4_prop1597.json` (CPU linear algebra + Max+ p=3,5; GPU unused).


**Proposition 15.98 (mult$(\lambda_2)\ge d-1$ via PSL; strengthened gap; 2026-07-30).** Continue Prop 15.95–15.97. Scope: Paley Max+ of order $n=p^2+1$.

1. **PSL min irrep (proved).** For $q=p^2$ odd, every nontrivial complex irrep of $\mathrm{PSL}(2,q)$ has dimension $\ge(q-1)/2=d-1$. (Character table of $\mathrm{PSL}(2,q)$.)

2. **Aut action (proved).** $\mathrm{P}\Sigma\mathrm{L}(2,q)$ acts on coordinates of the Paley conference preserving $C$ and $\mathrm{Max}_+$, hence on $L^2(\mathrm{Max}_+)$. $P\odot P$ is equivariant. The antipodal map $y\mapsto -y$ is a design automorphism.

3. **mult$\ge d-1$ (proved for Paley Max+).** The $\lambda_2$-eigenspace $V$ is orthogonal to constants ($\lambda_1$ simple) and nonzero ($\mathrm{Tr}(P\odot P)=d^2/N>\alpha$). It is a nontrivial unitary representation of $\mathrm{PSL}(2,q)$ (image of the acting group: $\mathrm{PSL}(2,q)$ is simple for $q\ge4$, action nontrivial $\Rightarrow$ image $\cong\mathrm{PSL}(2,q)$). Hence $\dim V\ge d-1$. $\square$

4. **Strengthened gap criterion (proved algebra, all primes $p\ge5$).** If $\mathrm{mult}(\lambda_2)\ge d-1$ and $\|\kappa\|_F^2\le96n$ (i.e. $\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}$), then $\lambda_2\le\sqrt{Q/(d-1)}\le d/(2N)$, so bi-tight empty. *Proof.* $N$ cancels; the inequality $\mathrm{Wick}/16-d^2\le(d-1)d^2/4$ reduces to $d(d^2-9d-24)\ge0$, true for $d\ge13$. At $p=3$ the algebra fails (correct). $\square$

5. **Certified.** mult$\ge d-1$ (in fact $=d$) and $\|\kappa\|^2\le96n$ at $p=3,5$; gap-by-mult$_{d-1}$ holds at $p=5$, fails at $p=3$.

6. **OPEN for bi-tight at general $p\ge5$.** Prove $\|\kappa\|_F^2\le96n$ (boolean $\le$ Wick). Then (3)+(4) close bi-tight for all Paley $p\ge5$. Deep ND + Main Theorem remain for $\lim\alpha_n$. **$L$ OPEN.**

Evidence: `src/e1_gmin_m4_prop1598.py`, `evidence/e1_gmin_m4_prop1598.json`.

**Proposition 15.99 (κ-budget structure; min-distance; closed forms; 2026-07-31).** Continue Prop 15.96–15.98.

1. **Min Hamming distance (proved, any conference Max+).** If $y,z\in\mathrm{Max}_{+}$, $y\ne\pm z$, then $d_H(y,z)\ge p+1$, i.e. $|y\cdot z|\le n-2(p+1)=(p-1)^2-2$. *Proof.* $v=(y-z)/2\in\{0,\pm1\}^n$ has support size $k=d_H$ and $Cv=pv$. Then $pk=v^\top Cv\le k(k-1)$, so $k\ge p+1$. $\square$

2. **Wick/ρ budget (proved algebra).** $\|\kappa\|_F^2=C_{\mathrm{diag}}+24\sum_S\rho_S^2$ with $C_{\mathrm{diag}}=4n(11n-14)/p^2$ and $\rho_S=m_4(S)-\kappa(S)/p^2$. Hence $\|\kappa\|^2\le96n$ iff $\sum\rho^2\le n(13p^2+3)/(6p^2)$ iff $\sum M^2\le\mathrm{Wick}_{\mathrm{hi}}$ iff $\mathbb E[(y\cdot z)^4]\le\mathrm{Wick}_{\mathrm{hi}}$. Room $96n-C_{\mathrm{diag}}=4n(13p^2+3)/p^2$. $\square$

3. **Closed forms (proved).** $\sum\kappa=p^2(p^2-1)/4$, $\sum\kappa^2=n(n-1)(n-2)(n-5)/8$, $\sum m_4=-p(p-1)(p+1)(p+4)/12$, $\sum m_4\kappa=n(n-1)(n-2)/8$; $\sum\rho$ and $\sum\rho\kappa$ follow. Stratum $n_1+n_3=\binom{n}{4}$. $\square$

4. **Master residual source (proved).** $(4pI-T)\rho=T\kappa/p^2$ with $T\kappa=0$ on $|\kappa|=1$ and $T\kappa=8\kappa$ on $|\kappa|=3$; $\|\mathrm{RHS}\|^2=576 n_3/p^4$. $\square$

5. **Gap link.** Paley mult$(\lambda_2)\ge d-1$ (15.98) + $\|\kappa\|^2\le96n$ $\Rightarrow$ bi-tight empty for all primes $p\ge5$. $\square$

6. **Certified.** At $p=3,5$: $\|\kappa\|^2\le96n$ (eq only $p=3$), min-distance, closed forms vs census, split identity. $\square$

7. **OPEN.** Prove $\|\kappa\|_F^2\le96n$ for all primes $p\ge5$ (boolean $\le$ Wick), or $16N$/$H$. Deep ND remains. **Existence of $\lim\alpha_n$ remains OPEN. $H_{\mathrm{proved}}=\mathrm{false}$.**

Evidence: `src/e1_gmin_m4_prop1599.py`, `evidence/e1_gmin_m4_prop1599.json` (CPU Fraction + Max+ p=3,5; GPU unused).

**Proposition 15.100 (dual-frame projection; flat Veronese ≤ Wick; κ_hyp; 2026-07-31).** Continue Prop 15.96–15.99.

1. **Dual frame (proved).** $r_j=P_+e_j$, $S=\sum_j r_j^{\otimes 4}$. For boolean measures $\langle S,M\rangle=n$; on Max+ $\langle S,\kappa\rangle=-2n$ and $\langle\mathrm{Wick},\kappa\rangle=-48n$. $\|S\|_F^2=n^2/(16p^2)$. $\square$

2. **Projection formula (proved).** $\|\kappa_{\mathrm{proj}}\|_F^2=64n(p^2-3)/(p^2-5)$ on $\mathrm{span}\{\mathrm{Wick},S\}$; room $96n-\|\kappa_{\mathrm{proj}}\|^2=32n(p^2-9)/(p^2-5)$. $\square$

3. **Flat Veronese bound (proved).** With $\mathrm{rank}(P\odot P)=\binom{d-1}{2}$, CS gives $E[D^4]\ge\mathrm{ED4}_{\mathrm{flat}}=16d^2+32d(d-1)^2/(d-3)$, and $\mathrm{Wick}-\mathrm{ED4}_{\mathrm{flat}}=64d(d-5)/(d-3)\ge0$ for $d\ge5$ (eq only $p=3$). Moreover $\kappa_{\mathrm{flat}}=\|\kappa_{\mathrm{proj}}\|^2$. $\square$

4. **κ_hyp algebra (proved).** $\kappa_{\mathrm{hyp}}=\|\kappa_{\mathrm{proj}}\|^2+\mathrm{room}\cdot((n-2)/n)^2$ satisfies $\kappa_{\mathrm{hyp}}\le96n$ with slack $128p^2(p-3)(p+3)/((p^2-5)(p^2+1))$. $\square$

5. **Certified.** At $p=3,5$: $\|\kappa\|^2=\kappa_{\mathrm{hyp}}$ (eq $96n$ only $p=3$). $\square$

6. **OPEN.** Prove $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}\cdot((d-1)/d)^2$ (i.e. $\|\kappa\|^2\le\kappa_{\mathrm{hyp}}$) for all primes $p\ge5$. Then $\|\kappa\|^2\le96n$ and bi-tight closes via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15100.py`, `evidence/e1_gmin_m4_prop15100.json`.

**Proposition 15.101 (Fickus Gram residual / bulk-variance orth reduction; 2026-07-31).** Continue Prop 15.100. Method transfer from Fickus–Jasper–Mixon (arXiv:2605.28738) residual Gram $K=H\odot\overline H$ / Schur rank-nullity and Ge–Liu multiplicity bounds (arXiv:2606.29392), applied to $\mathrm{PopP}=P\odot P$ on Max+.

1. **Fickus–Schur rank (proved).** $G=YY^\top$ has rank $d$; $K=G\odot G=4N^2\,\mathrm{PopP}$ has $\mathrm{rank}=\binom{d-1}{2}=1+m$ with $m=d(d-3)/2$ (Prop 15.59). Schur bound $\mathrm{rank}\le d^2$ is strict. Annihilator: $\mathrm{range}(P)\subset\ker(\mathrm{PopP})$; $\lambda_1=d/N$ simple; bulk sums to $S=d(d-1)/N$. $\square$

2. **Orth as bulk variance (proved).** $E[D^4]=16N^2\mathrm{Tr}(\mathrm{PopP}^2)=16d^2+16N^2\sum_{\mathrm{bulk}}\lambda^2$. Flat CS recovers $\mathrm{ED4}_{\mathrm{flat}}$. Hence
   \[
   \|\kappa_{\mathrm{orth}}\|_F^2=E[D^4]-\mathrm{ED4}_{\mathrm{flat}}=16N^2\Bigl(\sum_{\mathrm{bulk}}\lambda^2-\frac{S^2}{m}\Bigr).
   \]
   Dual-frame IDs: $\mathrm{ED4}_{\mathrm{flat}}-\mathrm{wick}_{\mathrm{lo}}=\|\kappa_{\mathrm{proj}}\|^2$, $\mathrm{Wick}-\mathrm{ED4}_{\mathrm{flat}}=\mathrm{room}$. $\square$

3. **PSL level count (proved for Paley).** $\mathrm{mult}(\lambda_2)\ge d-1$ (Prop 15.98) $\Rightarrow$ # distinct positive bulk eigenvalues $\le\lfloor m/(d-1)\rfloor=\lfloor d(d-3)/(2(d-1))\rfloor$. $\square$

4. **N-free $\lambda_2$-sufficient criterion (proved algebra).** Majorization $\sum\lambda^2\le\lambda_2\cdot S$ yields $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}$ whenever
   \[
   \lambda_2\le\lambda_{\mathrm{flat}}\cdot\bigl(1+\varepsilon\bigr),\qquad
   \varepsilon=\frac{4(p^2-9)}{(p^2-1)^2},\qquad
   \lambda_{\mathrm{flat}}=\frac{S}{m}.
   \]
   ($N$ cancels after substituting $d=(p^2+1)/2$.) Hyp form: $\varepsilon_{\mathrm{hyp}}=\varepsilon\cdot((d-1)/d)^2$. **Sufficient, not necessary** (at $p=5$ true $\lambda_2>\lambda_{\mathrm{flat}}(1+\varepsilon)$ while orth still equals $\mathrm{room}\cdot((d-1)/d)^2$). $\square$

5. **Certified $p=3,5$.** Bulk flat at $p=3$ (orth$=0$); three bulk levels at $p=5$ with mults $(d,2d,2d)$ and eigs $11/845,9/845,5/845$; $\kappa^2=\kappa_{\mathrm{hyp}}$ both; Fickus rank matches; $\lambda_2\le4/N$ at $p=5$. $\square$

6. **OPEN.** Prove $\|\kappa_{\mathrm{orth}}\|^2\le\mathrm{room}\cdot((d-1)/d)^2$ for all primes $p\ge5$ (closed bulk spectrum / mult-aware majorization / resolvent $\delta$). Then bi-tight via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15101.py`, `evidence/e1_gmin_m4_prop15101.json`.

**Proposition 15.102 (resolvent $\delta$-calculus; $\|\kappa_{\mathrm{orth}}\|^2=24\|\delta\|_2^2$; 2026-07-31).** Continue Prop 15.100–15.101. Isolates the orth residual as the ker component of the master resolvent.

1. **Source (proved).** $b=T\kappa/p^2$ is supported on $|\kappa|=3$ with $\|b\|_2^2=576 n_3/p^4$, $n_3=n(n-1)(n-2)(n-6)/96$. $\square$

2. **$\mu^2$-form (proved form; cert $p=3,5$).** $T^2b=\mu^2 b$ with $\mu^2=4(p^2+15)$; $\langle b,Tb\rangle=0$, so $b$ splits equally into $\pm\mu$ eigenspaces of $T$. $\square$

3. **$\rho_{\min}$ closed form (proved).** $\rho_{\min}=b_+/(4p-\mu)+b_-/(4p+\mu)$ has
   \[
   \|\rho_{\min}\|_2^2=\frac{5n(p^2-1)(p^2+3)}{6p^2(p^2-5)}.
   \]
   $\square$

4. **$\kappa_{\min}=\mathrm{proj}$ (proved).** $C_{\mathrm{diag}}+24\|\rho_{\min}\|_2^2=\|\kappa_{\mathrm{proj}}\|_F^2$ for every odd prime $p\ge3$. $\square$

5. **Orth identity (proved).** Max+ $\rho=\rho_{\min}+\delta$ with $\delta\in\ker(4pI-T)$ and $\rho_{\min}\perp\delta$, hence
   \[
   \|\kappa\|_F^2=\|\kappa_{\mathrm{proj}}\|_F^2+24\|\delta\|_2^2,\qquad
   \|\kappa_{\mathrm{orth}}\|_F^2=24\|\delta\|_2^2.
   \]
   $\square$

6. **Target equivalences (proved).** $\|\kappa\|^2\le96n\Leftrightarrow\|\delta\|_2^2\le\mathrm{room}/24$; $\|\kappa\|^2\le\kappa_{\mathrm{hyp}}\Leftrightarrow\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$. $\square$

7. **Invertible case (proved).** $\lambda_{\max}(T)<4p\Rightarrow\delta=0\Rightarrow\|\kappa\|^2=\mathrm{proj}\le96n$ (holds at $p=3$). $\square$

8. **Certified $p=5$.** $\lambda_{\max}(T)=4p$, $\mathrm{mult}(E_{4p})=d-1=12$, $\|\delta\|_2^2=\mathrm{room}_{\mathrm{hyp}}/24$ exactly ($\kappa^2=\kappa_{\mathrm{hyp}}$). $\square$

9. **OPEN.** Prove $\|\delta\|_2^2\le\mathrm{room}\cdot((d-1)/d)^2/24$ for all primes $p\ge5$. Then bi-tight via 15.98. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15102.py`, `evidence/e1_gmin_m4_prop15102.json`.

**Proposition 15.103 ($\delta$-bound certified $p=3,5,7$; 16N at $p=5,7$; 2026-07-31).** Continue Prop 15.102.

1. **Budget form (proved).** $\mathrm{room}_{\mathrm{hyp}}/24=4(p^2-9)(p^2-1)^2/(3(p^2-5)(p^2+1))$. $\square$

2. **Census (certified).** Full Max+ at $p=3,5,7$ ($N=12,260,11452$):
   - $p=3$: $\delta^2=0=\mathrm{room}_{\mathrm{hyp}}/24$, $\kappa^2=\kappa_{\mathrm{hyp}}=96n$;
   - $p=5$: $\delta^2=1536/65=\mathrm{room}_{\mathrm{hyp}}/24$, $\kappa^2=\kappa_{\mathrm{hyp}}<96n$;
   - $p=7$: $\delta^2\approx10.424\le\mathrm{room}_{\mathrm{hyp}}/24\approx55.855$ (ratio $\approx0.187$), $\kappa^2\approx3595.6<\kappa_{\mathrm{hyp}}\approx4686<4800=96n$ (strict).
   Equality $\kappa^2=\kappa_{\mathrm{hyp}}$ is **not** universal. $\square$

3. **16N (certified $p=5,7$).** $\lambda_2(P\odot P)\le4/N$ at $p=5$ (exact $11/845<4/260$) and $p=7$ (power method $\lambda_2\approx2.31\cdot10^{-4}<4/11452$). Hence 16N holds; bi-tight empty at these primes via Prop 15.61. Also via $\kappa^2\le96n$+Prop 15.98. $\square$

4. **OPEN.** Prove $\|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24$ (or $\lambda_2\le4/N$) for all primes $p\ge5$. **$\lim\alpha_n$ OPEN.**

Evidence: `src/e1_gmin_m4_prop15103.py`, `evidence/e1_gmin_m4_prop15103.json`; Max+ $p=7$ at `/tmp/e1_p7/maxplus.npy`.




**Proposition 15.22 (liminf controlled by the universal cube/sphere floor).** Write
\[
\rho_{\min}(n)\,:=\,\min_{A\in\mathcal S_n}\rho(A).
\]
Then for every \(n\ge2\),
\begin{equation}
\label{eq:rho-min-lb}
m_n
\;\ge\;
\frac12\,n\sqrt{n-1}\,\rho_{\min}(n),
\qquad
\alpha_n
\;\ge\;
\frac12\sqrt{1-\frac1n}\,\rho_{\min}(n).
\end{equation}
Consequently
\[
\liminf_{n\to\infty}\alpha_n
\;\ge\;
\frac12\,\liminf_{n\to\infty}\rho_{\min}(n).
\]
In particular, any uniform lower bound \(\rho_{\min}(n)\ge\rho_0-o(1)\) upgrades Theorem A's liminf from \(2^{-5/2}\) to \(\rho_0/2\). The Nesterov value \(\rho_0=2/\pi\) would give \(\liminf\alpha_n\ge1/\pi\approx0.3183\).

*Proof.* For each \(A\), \(\Phi(A)=\tfrac12 n\|A\|_{\mathrm{op}}\rho(A)\) and \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\) (Prop 15.5), so \(\Phi(A)\ge\tfrac12 n\sqrt{n-1}\,\rho(A)\). Minimise over \(A\). \(\square\)

**Proposition 15.23 (exhaustive \(\rho_{\min}\) for \(n\le8\)).** Vertex-folded exhaustive search over all Seidel matrices of orders \(n\in\{6,7,8\}\) (respectively \(2^{10}\), \(2^{15}\), \(2^{21}\) matrices; 86-worker process pool, `OMP_NUM_THREADS=1`) yields
\[
\begin{array}{c|ccc}
n&6&7&8\\ \hline
\rho_{\min}(n)&0.745356&0.721996&0.693375\\
\min r(A)&0.745356&1.049781&0.944911\\
\min\Phi=m_n&5&9&10
\end{array}
\]
In all three orders \(\rho_{\min}(n)>2/\pi\approx0.6366\). At \(n=6\), \(\rho_{\min}=\min r=\rho(C)\) (Paley is optimal for both). At \(n=8\), \(\min\Phi=m_8=10\) matches the known exact value. Cross-check: \(\min\Phi=m_n\) on the shipped `exact_m` table for \(n=6,7,8\).

*Remark.* The sequence \(\rho_{\min}(6,7,8)\) is decreasing and still above \(2/\pi\). Whether \(\rho_{\min}(n)\to2/\pi\) (which would give \(\liminf\alpha\ge1/\pi\) by Prop 15.22) is open; a matching universal Nesterov theorem for all Seidel matrices would close it. Cluster-Nesterov numerics support \(\rho(A)\ge2/\pi-o(1)\) but are not a proof.

### §15.3 Why the full transfer still fails

- Prop 15.2 lower-bounds \(\Phi(C)\), not \(m_n\).
- The naive claim \(\mathrm{SDP}_+(A)\ge n\sqrt{n-1}\) for all Seidel \(A\) is **false**: the all-negative matrix \(A=J-I-2(J-I)=I-J\) (off-diagonal \(-1\)) has \(\mathrm{SDP}_+(A)=n\) (dual certificate \(M=J\)) while \(n\sqrt{n-1}>n\). That matrix has huge \(\Phi\) (\(\Phi=\binom{n}{2}\)) and is irrelevant for minimising \(\Phi\). What remains plausible is that \(\max\bigl(\mathrm{SDP}_+(A),\mathrm{SDP}_+(-A)\bigr)\) is minimized by conference matrices; this is not proved.
- Grothendieck / Kashin–Szarek / AMMN block naïve spherical transfer for general (non-Seidel) forms.
- Prop 15.16 kills the global super-linear-\(\min\delta\) repair of Prop 15.14.
- Prop 15.19 kills the *conditional* shell path for large \(n\) (shell \(=\) all of \(\mathcal S_n\) once \(n\gtrsim38\)).
- **Remaining gap for Theorem E(1):** prove \(\rho(A)\|A\|_{\mathrm{op}}\ge\bigl(\rho(C)-o(1)\bigr)\sqrt{n-1}\) for all \(A\in\mathcal S_n\) along Paley orders (product form of Prop 15.9), by a method that does **not** pass through \(\mathbb E[Q^4]\). Natural programme: (i) universal lower bound \(\rho(A)\ge2/\pi-o(1)\) via Nesterov rounding; (ii) rigidity of near-minimal-op Seidel matrices toward a conference class; (iii) prove that the resulting class has \(\rho=1-o(1)\), or compare its cube maximum directly with Paley; and (iv) continuity of cube-max under Frobenius perturbation. Spectral defect alone cannot identify the Paley switching class among inequivalent conference matrices.
- **Remaining gap for Theorem E(2):** prove convergence of the cube-imbalance of \(P_+\) (Prop 15.18). Exact \(\rho(C_n)\) is strictly increasing on \(\{6,14,18,30,38\}\); a monotonicity or Cauchy argument would close E(2).

**Invariant needed for \(L=\rho_*/2\):** asymptotic optimality (Prop 15.9) plus \(\rho(C_k)\to\rho_*\) (Theorem E).

---

## §16. Approach 6: non-existence — also fails

To prove the limit does **not** exist one needs two subsequences with different limit points of \(\alpha_n\). Denseness (Proposition 6.1) forces any such oscillation to be visible along *every* ratio-dense subsequence, including Paley orders and all arithmetic progressions. No construction is known that produces two different asymptotic densities of \(m_n/n^{3/2}\). Log-log periodic abstract sequences satisfy the soft axioms but are not realised by \(\min_A\Phi(A)\). Non-existence is therefore as open as existence.

---

## §17. Main results (rigorous content)

### Theorem A (sandwich) — complete

\[
\frac1\pi
\le\liminf_{n\to\infty}\alpha_n
\le\limsup_{n\to\infty}\alpha_n
\le\frac12.
\]
(The floor \(1/\pi\) is Prop.~5.2; Bohnenblust–Hille \(2^{-5/2}\) remains as a weaker Prop.~5.1.)

### Theorem B (structure) — complete

Equivalence (§1); monotonicity and padding (§3); denseness and the two-ray
conditional convergence theorem (§6); multipartite bounds including reverse
(§7); \(a_n\to\limsup\alpha_n\) (§8).

### Theorem C (obstruction) — complete

No multipartite bound of the form \(\alpha_{kn}\le\alpha_n/\sqrt k+c_k\) with \(c_k\to0\) (uniformly in \(n\)) can exist (§10). Soft comparison inequalities compatible with the sandwich admit abstract non-convergent sequences. Existence of \(\lim\alpha_n\) is **not** a formal consequence of Theorems A–B.

### Theorem D (conference spectral structure) — complete

Propositions 15.1–15.2 and Corollary 15.3: exact cube/sphere formula via \(P_+\); sharp Nesterov expectation \eqref{eq:nest-exact}–\eqref{eq:nest-rho}; sandwich \(1/\pi\le\liminf\Phi(C_n)/n^{3/2}\le\limsup\Phi(C_n)/n^{3/2}\le1/2\) along conference orders.

### Theorem G (switching + min-op + spectral/\(L^4\) extremality + limsup via \(\rho\)) — complete

Propositions 15.4–15.23:
- Seidel switching preserves \(\Phi\); \(\|A\|_{\mathrm{op}}\ge\sqrt{n-1}\) with equality iff conference;
- \(\Phi=\tfrac12 n\rho\|A\|_{\mathrm{op}}\); beaters of conference must have worse \(\rho\);
- \(L^2\)-universality \(\mathbb E[Q^2]=\binom{n}{2}\) for every Seidel \(A\);
- \(\mathrm{tr}(A^4)\ge n(n-1)^2\) with equality iff conference;
- **exact fourth-moment formula** \eqref{eq:Q4}: \(\mathbb E[Q^4]\) is an affine function of \(\mathrm{tr}(A^4)\) and is uniquely minimised at conference matrices;
- **exact optimality criterion** \eqref{eq:delta-bound}: \(\Phi(A)\le\Phi(C)\Rightarrow\delta(A)\le\Delta_*/3\); when the spectral gap of every non-conference matrix exceeds \(\Delta_*/3\), every possible beater is a conference matrix, so \(m_n\) is the minimum cube maximum over all conference classes (equal to \(\Phi(C)\) only when the chosen class attains that minimum);
- **\(m_6=\Phi(C)=5\)** by exhaustive verification of the gap criterion (Cor 15.15);
- asymptotic optimality \(\Leftrightarrow\) conference minimises \(r(A)=\max|x^\top Ax|/(n\sqrt{n-1})\);
- Prop 15.16: global \(\min\delta=\Theta(n)\) (edge flip), so super-linear gap repair is impossible;
- Prop 15.17–15.19: conditional shell criterion is equivalent to optimality at fixed \(n\), but the shell is **vacuous** for large conference orders (\(\Delta_*/3=\Theta(n^5)>\max\delta=O(n^4)\));
- Prop 15.18: \(\rho(C)\) is the cube-imbalance of the spectral projector \(P_+\);
- Prop 15.20: Lipschitz continuity of \(\Phi\) in Frobenius norm;
- Prop 15.21: single-edge local optimality under maximizer balance (verified for Paley \(n\le18\));
- Prop 15.22: \(\alpha_n\ge\tfrac12\sqrt{1-1/n}\,\rho_{\min}(n)\), so \(\liminf\alpha\ge\tfrac12\liminf\rho_{\min}\);
- Prop 15.23: exhaustive \(\rho_{\min}(n)>2/\pi\) for \(n\in\{6,7,8\}\) with \(\min\Phi=m_n\);
- and
\[
\limsup_n\alpha_n\le\tfrac12\limsup_k\rho(C_k)\le\tfrac12
\]
along Paley (strict improvement of Theorem A if \(\limsup\rho<1\)).

Conference matrices are the unique Seidel matrices simultaneously extremal for operator norm, \(\mathrm{tr}(A^4)\), \(\mathbb E[Q^4]\), and the universal cube-\(L^2\) mass, and they are exactly optimal for \(m_n\) at \(n=6\) and *locally* optimal under edge flips for Paley \(n\le18\). The open core of asymptotic optimality (Theorem E(1)) is an \(L^\infty\) comparison that **cannot** pass through fourth moments for large \(n\) (Props 15.16, 15.19) and requires the delocalization/rigidity programme after Prop 15.21. Prop 15.22 remains available as a spectral route to liminf; the dual-Gaussian Prop 5.2 already supplies the constant \(1/\pi\) uniformly, so a separate universal \(\rho_{\min}\ge2/\pi\) theorem is no longer needed for the sandwich floor.

### Existence of \(\lim\alpha_n\)

**Not established.** Approaches 1–5 each fail for a specific structural reason
(§11–§15). Non-existence is equally unproved (§16). Proposition 6.3 gives the
direct gate: Dini-summable amplification at multipliers 2 and 3 would settle
existence.  Alternatively, Proposition 6.2 makes convergence along the
ratio-dense Paley orders necessary and sufficient; proving E(1) there would
settle the stronger value-specific claim \(L=1/2\).

**Conditional Theorem E (existence via conference optimality).** Suppose that along Paley orders \(n_k=q_k+1\):
1. *asymptotic optimality:* \(m_{n_k}=\Phi(C_k)+o(n_k^{3/2})\) (equiv.\ Prop 15.9); and
2. *cube/sphere convergence:* \(\rho(C_k)\to\rho_*\in(0,1]\).

Then \(\lim_{n\to\infty}\alpha_n=\rho_*/2\). In particular \(\rho_*=1\Rightarrow L=\tfrac12\).

*Proof.* By Prop 15.1, \(\Phi(C_k)/n_k^{3/2}=\tfrac12\sqrt{1-1/n_k}\,\rho(C_k)\to\rho_*/2\). Optimality gives \(\alpha_{n_k}\to\rho_*/2\). Apply Prop 6.2. \(\square\)

**Conditional Theorem F (Stolz regularity).** If \(\delta_n/\sqrt n\to\ell\in[0,\infty)\) where \(\delta_n=m_{n+1}-m_n\), then \(\lim\alpha_n=\tfrac23\ell\).

*Proof.* Stolz–Cesàro: \((n+1)^{3/2}-n^{3/2}\sim\tfrac32\sqrt n\), so
\[
\frac{m_{n+1}-m_n}{(n+1)^{3/2}-n^{3/2}}\to\frac{\ell}{3/2}=\frac23\ell,
\]
hence \(m_n/n^{3/2}\to\tfrac23\ell\). \(\square\)

**Status of the two hypotheses of Theorem E.**
- (1) Optimality: proved at \(n=6\); **fails exact optimality at \(n=10\)** (\(m_{10}=13<\Phi=15\)); Prop 15.21 gives *strict local* edge-opt for Paley \(n\in\{6,14,18\}\); multi-core SA on product \(r\) never undercuts \(\rho(C)\) at \(n=6,14\); Prop 15.20 controls Hamming balls of radius \(o(n)\). Gap at \(n=10\) is \(2=o(n^{3/2})\), so asymptotic E(1) remains plausible. Global \(L^4\) path dead (Props 15.16, 15.19). Delocalization/rigidity programme after Prop 15.21 — **not completed**.
- (2) \(\rho(C_k)\to\rho_*\): dual-Gauss / Nesterov \(\liminf\rho\ge2/\pi\); **exact** \(\rho(C_n)\) strictly increasing
  \(\{0.745,0.832,0.889,0.928,0.943,0.959\}\) on \(n\in\{6,14,18,30,38,42\}\).
  Multi-core local search (SCRATCH `attack_E2_rho`): \(\rho_{\mathrm{LB}}(110)\ge0.9665\).
  Constructive **interval signings** on Paley (SCRATCH `interval_rho_large` + `interval_mega`,
  primes \(q\le15000\)):
  \[
  \rho_{\mathrm{int}}(13382)\ge0.9882,\quad
  \rho_{\mathrm{int}}(8762)\ge0.9867,\quad
  \rho_{\mathrm{int}}(3530)\ge0.9826,
  \]
  with \((1-\rho_{\mathrm{int}})\sqrt n=\Theta(1)\) on record-setters — strong constructive evidence for
  \(\limsup\rho=1\), and rigorously
  \[
  \limsup_k\rho(C_{n_k})\;\ge\;0.9882.
  \]
  Full E(2) (\(\rho\to\rho_*\) for *all* large Paley orders, especially \(\rho_*=1\)) is **not proved**:
  an analytic bound \(\rho\ge1-O(n^{-1/2})\) uniform in \(q\) is missing. However, interval quality
  is a poor proxy for true \(\rho\): at \(q=1013\) (\(n=1014\)) one has \(\rho_{\mathrm{int}}\approx0.74\)
  but multi-core local search gives \(\rho_{\mathrm{LB}}\ge0.9259\) (SCRATCH `bad_q1013.json`).
  Thus apparent “dips” in \(\rho_{\mathrm{int}}\) do not certify dips in true \(\rho\).
  Prop 15.8: \(\limsup\alpha\le\tfrac12\limsup\rho\).

**Proposition 15.12 (quantitative near-optimality \(\Rightarrow\) near-conference).** Let \(A\in\mathcal S_n\) satisfy \(r(A)\le R\) and suppose \(\rho(A)\ge\rho_0>0\). Then
\[
\frac{\|A\|_{\mathrm{op}}}{\sqrt{n-1}}\le\frac{R}{\rho_0}.
\]
In particular, if a near-minimiser of \(r\) has cube/sphere ratio bounded below by a constant \(\rho_0\) comparable to \(\rho(C)\), its operator norm is \(O(\sqrt n)\) with leading ratio \(R/\rho_0\). Combined with Prop 15.11,
\[
\mathrm{tr}(A^4)-n(n-1)^2\le n(n-1)^2\Bigl(\Bigl(\frac{R}{\rho_0}\Bigr)^2-1\Bigr),
\]
so \(R\to\rho_*\) and \(\rho_0\to\rho_*\) forces spectral 4th-moment excess \(\to0\), i.e.\ eigenvalue squares concentrate at \(n-1\).

*Proof.* \(r(A)=\rho(A)\,\|A\|_{\mathrm{op}}/\sqrt{n-1}\le R\) rearranges to the op bound. Then \(\max\lambda_i^2\le\|A\|_{\mathrm{op}}^2\le(n-1)(R/\rho_0)^2\), and \(\sum\lambda_i^4\le(\max\lambda_i^2)\sum\lambda_i^2\) yields the trace excess. \(\square\)

This is the natural rigidity route to Theorem E(1): a matching lower bound \(\rho(A)\ge\rho(C)-o(1)\) for every near-minimiser of \(r\), plus \(\rho(C_k)\to\rho_*\), closes asymptotic optimality.

### What would close the problem

1. The two Dini-summable multiplier estimates of Proposition 6.3.
   Proposition 6.5 makes the multiplier-two endpoints equal and reduces the
   ray to the exact mixed-state diamond (6.13); Proposition 6.6 proves it
   outside the explicit residue (6.20). That residue remains, and the
   multiplier-three estimate is independently required. Proposition 6.7
   reduces direct tripling to its tetrahedral diamond, while Proposition 6.8
   reduces the alternative `1:2` composition to (6.42)--(6.43).
2. Asymptotic conference optimality + \(\rho(C_k)\to\rho_*\) (Theorem E) — reduced by Theorem G to product-minimisation of \(\rho\cdot\|A\|_{\mathrm{op}}\).
3. Extension regularity \(\delta_n/\sqrt n\to\ell\) (Theorem F), or the stronger \(\gamma(A^*)=(\tfrac32\alpha_n+o(1))\sqrt n\).
4. Maximizer delocalisation + discrepancy feeding (3).
5. Multipartite rigidity \(\alpha_{kn}\le\alpha_n+o(1)\) blocking \(\lambda\to\Lambda\).
6. Explicit two-density construction for non-existence.

None of (1)–(6) is fully available.  Item 1 is the shortest direct route to
the original existence question; item 2 is the stronger value-specific
Paley route.

---

## §18. Numerical evidence

Shipped library: `src/minmax_quadratic.py`.

| \(n\) | \(m_n\) | \(\alpha_n\) | notes |
|------:|--------:|-------------:|:------|
| 2 | 1 | 0.3536 | |
| 3 | 3 | 0.5774 | |
| 4 | 4 | 0.5000 | |
| 5 | 4 | 0.3578 | flat extension |
| 6 | 5 | 0.3402 | Paley optimal, 12 maximizers |
| 7 | 9 | 0.4860 | |
| 8 | 10 | 0.4419 | |
| 9 | 12 | 0.4444 | exact Gray 86-worker |
| 10 | **13** | **0.4111** | exact Gray 86-worker; Paley-\(q=9\) has \(\Phi=15>13\) so conference **not** exact-optimal |
| 11 | \(\le17\) | \(\le0.4660\) | explicit witness (antipodal cut-code); dual-Gauss LB \(\approx11.07\) |

Optimal extension increments matching \(m_{n+1}-m_n\): \(n=3\to1\), \(4\to0\), \(5\to1\), \(6\to4\), \(7\to1\).

Extension cost \(\gamma(A^*)=\min_s\max_x(|Q|+|s\cdot x|)-\Phi(A^*)\) for one optimal \(A^*\) (exact \(s\)-search; equals \(\delta_n=m_{n+1}-m_n\) on this range):

| \(n\) | \(\alpha_n\) | \(\gamma/\sqrt n\) | \(\tfrac32\alpha_n\) | \(\delta_n\) |
|------:|-------------:|-------------------:|---------------------:|-------------:|
| 3 | 0.577 | 0.577 | 0.866 | 1 |
| 4 | 0.500 | 0.000 | 0.750 | 0 |
| 5 | 0.358 | 0.447 | 0.537 | 1 |
| 6 | 0.340 | 1.633 | 0.510 | 4 |
| 7 | 0.486 | 0.378 | 0.729 | 1 |
| 8 | 0.442 | 0.707 | 0.663 | 2 |

At the local \(\alpha\)-minimum \(n=6\) (Paley), extension is expensive (\(\gamma/\sqrt n\gg\tfrac32\alpha\)); after the climb to \(n=7\), extension is cheap. This is consistent with mean-reversion of \(\alpha\) but is not a proof of uniqueness of the fixed point.

Paley conference spectral diagnostics (shipped `paley_conference_matrix` + `phi`/`phi_local`; 14-way process pool, `OMP_NUM_THREADS=1`; identity Prop 15.1 recon error \(0\); Nesterov closed form Prop 15.2). Here \(\rho=\max|x^\top Cx|/(n\sqrt{n-1})\), \(\alpha_{\mathrm{UB}}=\Phi(C)/n^{3/2}\), and nest LB \(=\frac2\pi\sqrt{n-1}\arcsin(n-1)^{-1/2}\):

| \(n\) | method | \(\rho\) | nest LB | \(\alpha_{\mathrm{UB}}\) | \(\Phi\) |
|------:|:------:|---------:|--------:|-------------------------:|--------:|
| 6 | exact | 0.7454 | 0.6600 | 0.3402 | 5 |
| 14 | exact | 0.8321 | 0.6451 | 0.4009 | 21 |
| 18 | exact | 0.8893 | 0.6430 | 0.4321 | 33 |
| 30 | exact | 0.9285 | 0.6403 | 0.4564 | 75 |
| 38 | exact | 0.9431 | 0.6395 | 0.4653 | 109 |
| 42 | **exact** | **0.9594** | 0.6392 | **0.4739** | **129** |
| 54 | local LB | 0.9615 | 0.6386 | 0.4763 | \(\ge189\) |
| 74 | local LB | 0.9395 | 0.6381 | 0.4666 | \(\ge297\) |
| 98 | local LB | 0.9428 | 0.6377 | 0.4690 | \(\ge457\) |
| 114 | local LB | 0.9490 | 0.6376 | 0.4724 | \(\ge569\) |

**Exact \(\rho\) via 86-worker numba Gray-code half-cube** (`exact_rho_numba.py`, \(2^{n-1}\) patterns): \(n\in\{6,14,18,30,38,42\}\) all match prior local-search lower bounds *exactly* (\(\Phi\in\{5,21,33,75,109,129\}\); \(n=42\) wall \(2087\)s). The exact sequence
\[
\rho(C_n)\in\{0.7454,\,0.8321,\,0.8893,\,0.9285,\,0.9431,\,0.9594\}
\]
is **strictly increasing**. For \(n>42\), multi-core local search + eigenspace sampling (`attack_E2_rho`, 86 workers) recovers
\[
\rho_{\mathrm{LB}}(54)=0.9615,\;
\rho_{\mathrm{LB}}(110)=0.9665,\;
\rho_{\mathrm{LB}}(62)=0.9458,\;
\rho_{\mathrm{LB}}(74)=0.9395,
\]
and **constructive interval signings** on the Paley field (FFT character-sum form; SCRATCH `interval_rho` / large-\(q\) sweep) give
\[
\rho_{\mathrm{int}}(242)\ge0.9545,\;
\rho_{\mathrm{int}}(1010)\ge0.976,\;
\rho_{\mathrm{int}}(1130)\ge0.9785.
\]
Hence rigorously
\[
\limsup_{k\to\infty}\rho(C_{n_k})\;\ge\;0.9785
\]
along Paley (true \(\rho\ge\rho_{\mathrm{int}}\)). This is consistent with \(\rho\to1\) but does **not** prove E(2): interval quality oscillates with \(q\) (some orders only \(\rho_{\mathrm{int}}\approx0.73\)), and local-search dips past \(n=54\) remain uncertified as exact \(\Phi\). Nest / dual-Gauss floor \(\to2/\pi\). Prop 15.19: Q4 shell is all of \(\mathcal S_n\) for \(n\ge38\).

Optima vs random \(\rho\) (86-worker `rho_stats`, brute \(n\le18\), local larger; distinct PIDs \(=85\)):

| class | \(n\) | \(\rho\) mean/value |
|:------|------:|--------------------:|
| exact opt | 3–8 | \(1.00\to0.69\) |
| random | 8–32 | mean \(0.84\to0.78\) |
| Paley | 30–114 | \(0.93\)–\(0.96\) (local LB) |

Product ratio \(r(A)=\max|x^\top Ax|/(n\sqrt{n-1})\) (86 workers, `product_ratio.py`; conference optimality \(\Leftrightarrow\) min \(r\)):

| class | \(n\) | \(r\) min | \(r\) mean | \(\mathrm{tr}(A^4)/n(n-1)^2\) min |
|:------|------:|----------:|-----------:|----------------------------------:|
| Paley | 6 | 0.745 | 0.745 | 1.000 |
| Paley | 14–62 | 0.83–0.96 | =min | 1.000 |
| random | 6 | 1.044 | 1.366 | 1.427 |
| random | 14 | 1.228 | 1.447 | 1.541 |
| random | 24 | 1.321 | 1.464 | 1.749 |

No random sample undercuts Paley in \(r(A)\); Prop 15.11 forces \(\mathrm{tr}(A^4)\) strictly above the conference floor for non-conference matrices.

Paley \(\Phi\) and maximizer counts (half-cube, \(x_1=+1\)):

| \(n\) | \(\Phi\) (Paley) | \(\#\) maximizers | \(\alpha\) upper |
|------:|-----------------:|------------------:|-----------------:|
| 6 | 5 | 12 | 0.340 |
| 14 | 21 | 156 | 0.401 |
| 18 | 33 | 204 | 0.432 |

Exact values where stated, together with older local-search upper bounds
(86-worker SA with **exact** \(\phi\) for \(n\le12\); reconfirmed under
SCRATCH `attack_paley_product`):

| \(n\) | best \(\Phi\) found | notes |
|------:|--------------------:|:------|
| 6 | 5 | hits Paley / exact \(m_6\) |
| 8 | 10 | matches exact \(m_8\) |
| 10 | **13** | matches exact \(m_{10}\) |
| 11 | **17** | exact external cut-code certificate independently replayed in full |
| 12 | **18** | checked witness; exact from \(m_{11}=17\), monotonicity, and even parity |
| 13 | **20** | exact published threshold certificate; local survivor/witness checks passed, billion-record completeness externally trusted |
| 14 | **21** | checked Paley witness; exact from \(m_{13}=20\), monotonicity, and odd parity |
| 15 | **27** | exact published threshold-tower certificate; local top checks passed, tower completeness externally trusted |
| 16 | 34 | prior local \(\phi\) only |

Product ratio \(r\) at Paley orders (86-worker SA): at \(n=6,14\) no sample undercuts \(\rho(C)\); consistent with E(1) at these orders but not a proof.

**Local optimality of Paley under edge flips** (`attack_E1_E2` + `attack_rigidity`, 86 workers): SA edge-flip from Paley at \(q\in\{5,13,17\}\) never improved \(\Phi\) (20 seeds each). Fixed \(k\)-edge flips (\(k=1,\ldots,12\), 15 seeds): **no sample undercuts Paley \(r\)** (`ANY_BEAT=False`). Single edge flip from Paley always yields \(\delta=16(n-2)\) and increases \(\Phi\) by \(+2\) at \(n=6,14,18\). SA minimising \(\|A\|_{\mathrm{op}}\) recovers conference at \(n=6\) (30/30) and never produces \(r<\rho(C)\) when a Paley competitor exists.

**Fourth-moment gap still inconclusive for \(n\ge14\):** random min-\(\delta\) samples grow roughly like \(n^{3.7}\) (log-log fit) but remain \(\ll\Delta_*/3\) (\(\Delta_*/3\approx7341\) at \(n=14\) vs sample min \(\delta=896\); structural single-flip floor \(16(n-2)=192\)). \(E[Q^6]\) polynomial fit unusable for a proof (relative residual \(\gg1\)).

Maximizer cardinalities are consistent with a polynomial bound \(|M|\le n^{O(1)}\) on these instances (not a proof).

```
pytest tests/test_minmax.py -v
# 27+ passed (session-parallel exact_m fixture; never live exact_m(9/10);
# Prop 5.2 dual-Gaussian; Prop 15.1–15.23; m_6 exact opt; m_9=12,m_10=13 recorded)
```

The local Gray harvest gives
\(m_6=5,\;m_7=9,\;m_8=10,\;m_9=12,\;m_{10}=13\).  External certificate
audits extend the exact table to
\(m_{11},\ldots,m_{15}=17,18,20,21,27\), with the large-stream trust
boundaries recorded separately rather than presented as local reruns.

---

## §19. Acceptance checklist

| AC | Status |
|----|--------|
| 1. Exact limit quantity | Yes (Statement) |
| 2. Existence (\(\liminf=\limsup\)) | **Open** — sandwich now \(1/\pi\le\liminf\le\limsup\le1/2\) (Prop 5.2); direct closure can use the two Dini-summable rays of Prop 6.3; the stronger Paley route still needs global E(1)+E(2); certified non-existence is the third alternative |
| 3. Equivalence | §1 complete |
| 4. Numerics | §18 + shipped tests + exact \(m_{\le15}\) with explicit local/external trust boundaries, dual-Gauss checks, and Paley/product SA under SCRATCH |
| 5. Honest obstruction analysis | §9–§17 |

---

## References

1. P. Ivanisvili, MathOverflow 413935 (2022); X status 2081070728422752329 (2026).
2. R.E.A.C. Paley, On orthogonal matrices, *J. Math. Phys.* 12 (1933).
3. A. Defant, M. Mastyło, P. Pérez, *Math. Ann.* (2019).
4. Z. Füredi, I. Z. Ruzsa, Nearly subadditive sequences, arXiv:1810.11723.
5. N. Alon, K. Makarychev, Y. Makarychev, A. Naor, Quadratic forms on graphs, *Invent. Math.* (2006).
6. B. S. Kashin, S. J. Szarek, On the Gram matrices of systems of uniformly bounded functions, *Proc. Steklov Inst.* (2003).
7. Yu. Nesterov, Semidefinite relaxation and nonconvex quadratic optimization, *Optim. Methods Softw.* (1998).
8. H. J. Brascamp and E. H. Lieb, Best constants in Young's inequality, its
   converse, and its generalization to more than three functions,
   *Adv. Math.* **20** (1976), 151--173,
   doi:10.1016/0001-8708(76)90184-5.

## Prop 15.105 (2026-07-31) — Φ variance = κ_orth; 16-criterion

**Proved (Fraction algebra, conference Max+):**
1. ∑_α λ_α(Φ|Z)² = ED4 − 4n².
2. **Variance identity:** ∑_α (λ_α − μ̄)² = ‖κ_orth‖_F² with μ̄ = 8(n−2)/(n−6).
   Hence orth=0 ⇒ Φ is scalar on Z (Fickus residual flat ⇔ Norton operator scalar).
3. **Exact 16-criterion:** if mult(λ_max)≥d and ‖κ‖_F²≤96n, then
   λ_max ≤ μ̄ + 8(p²−9)/(p²−5) = 16, so 16N.

**Certified:** p=3,5,7 (mult(top)=d; 16N; variance match).

**OPEN:** ‖κ‖²≤96n (or direct λ_max≤16) for all primes p≥5. L remains OPEN.

## Prop 15.106 (2026-07-31) — Rest-average-8; kurtosis form of residual

**Proved:**
1. If mult(λ_max(Φ|Z)) ≥ d then λ_max ≤ 16 ⇔ mean(rest spectrum) ≥ 8.
   At λ_max = 16 the rest mean equals **exactly 8** (the Wick baseline):
   (T − 16d)/(m − d) = 8 for conference d = n/2 > 5.
2. ‖κ‖_F² ≤ 96n ⇔ ED4 ≤ wick_hi ⇔ kurtosis κ₄(y·z) ≤ 3 + 12/n.

**Certified** p=3,5,7 (kurtosis bound, 16N, mult=d).

**OPEN:** κ₄ ≤ 3+12/n or λ_max(Φ)≤16 for all primes p≥5. L remains OPEN.

## Prop 15.107 (2026-07-31) — 16N from mult≥d−1 + room_hyp

**Proved (Theorem A):** For every prime \(p\ge 3\),
\[
\mathrm{mult}(\lambda_{\max}(\Phi|Z))\ge d-1
\quad\text{and}\quad
\|\kappa_{\mathrm{orth}}\|^2 \le \mathrm{room}_{\mathrm{hyp}}
\implies
\lambda_{\max}(\Phi|Z)\le 16.
\]
Proof: majorization with mult \(\ge d-1\) gives
\((16-\bar\mu)^2 \ge \mathrm{room}_{\mathrm{hyp}}\cdot(m-m_1)/(m\,m_1)\)
as a Fraction inequality for all primes \(p\ge 3\). Hence 16N and (Prop 15.61) bi-tight empty for \(p\ge 5\).

**Also proved:** Gegenbauer \(Q_4(t)=t^4-\frac{6}{d+4}t^2+\frac{3}{(d+2)(d+4)}\); 4-design ED4 \(\le\) wick_hi.

**Certified:** orth≤room_hyp and 16N at \(p=3,5,7\).

**OPEN residual:** orth ≤ room_hyp for all primes \(p\ge 5\) (equivalently \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\)). L remains OPEN.

## Prop 15.108 (2026-07-31) — Residual-Gram / Schur dual; algebraic Thm A; Parseval \(T_\rho\)

**Proved (Fraction algebra, conference Max+):**

1. **Theorem A\* (polynomial form of Thm A).** For every prime \(p\ge 3\),
   \[
   (16-\bar\mu)^2 - \mathrm{room}_{\mathrm{hyp}}\cdot\frac{m-m_1}{m\,m_1}
   = \frac{128(p-3)(p+3)(p^4-12p^2-5)}{(p^2-5)^2(p^2+1)^2}\ge 0
   \]
   with \(m_1=d-1\), equality at \(p=3\). Hence
   \(\mathrm{mult}(\lambda_{\max})\ge d-1\) and \(\mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}\)
   imply \(\lambda_{\max}(\Phi|Z)\le 16\) for **all** primes \(p\ge 3\) as a polynomial
   identity (not sample checks).

2. **PopP\(\leftrightarrow\Phi\) conversion.** \(\lambda_{\max}(\Phi|Z)=4N\cdot\lambda_2(P\odot P)\).
   Hence 16N \(\Leftrightarrow\) \(\lambda_2(P\odot P)\le 4/N\).

3. **Residual-Gram / Schur dual.** Writing \(R=\mathrm{PopP}_{\mathrm{bulk}}-\lambda_{\mathrm{flat}}\Pi_{\mathrm{bulk}}\),
   \[
   \mathrm{orth}=16N^2\|R\|_F^2=16N^2\Bigl(\sum_{\mathrm{bulk}}\lambda_i^2-\frac{S^2}{m}\Bigr).
   \]
   (Fickus residual-Gram method transfer, arXiv:2605.28738.)

4. **Parseval \(\delta\)-target (single scalar residual).**
   \[
   \mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}
   \;\Leftrightarrow\;
   \|\delta\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24
   \;\Leftrightarrow\;
   \|\rho\|_2^2\le T_\rho(p),
   \]
   where \(T_\rho=\|\rho_{\min}\|_2^2+\mathrm{room}_{\mathrm{hyp}}/24\) is Max+-free:
   \[
   \|\rho_{\min}\|_2^2=\frac{5n(p^2-1)(p^2+3)}{6p^2(p^2-5)},\quad
   \frac{\mathrm{room}_{\mathrm{hyp}}}{24}=\frac{4(p^2-9)(p^2-1)^2}{3(p^2-5)(p^2+1)}.
   \]
   Path C residual is exactly \(\sum_S\rho(S)^2\le T_\rho(p)\).

5. **\(m_4\) expansion.** \(\sum m_4^2=\|\kappa\|_2^2/p^4+(2/p^2)\langle\kappa,\rho\rangle+\|\rho\|_2^2\)
   with \(\|\kappa\|_2^2=n(n-1)(n-2)(n-5)/8\) (Prop 15.71).

**Certified:** \(\sum\rho^2\le T_\rho\) (eq \(p=3,5\); strict ratio \(\approx0.639\) at \(p=7\));
16N via PopP at \(p=3,5,7\).

**OPEN residual:** \(\sum_S\rho(S)^2\le T_\rho(p)\) for all primes \(p\ge 5\)
(equivalently orth\(\le\)room_hyp). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15108.py`, `evidence/e1_gmin_m4_prop15108.json`, `tests/test_prop15108.py`.

## Prop 15.109 (2026-07-31) — Φ–m₄ identity; Aut-invariant δ; PF+rank obstruction

**Proved:**

1. **Φ–m₄ identity.** For every \(B\in Z\),
   \(\mathbb E[(y^\top B y)^2]=6\|B\|_F^2+8\sum_S m_4(S)\,\kappa_B(S)\).
   Hence \(16N\Leftrightarrow \max_{\|B\|=1}\langle m_4,\kappa_B\rangle\le 5/4\).

2. **\(\sum\kappa_B^2\) formula.** For zero-diag symmetric \(B\),
   \(\sum_S\kappa_B(S)^2=\frac14\mathrm{Tr}(B^4)+\frac18(\mathrm{Tr}B^2)^2+\frac12\sum B_{ij}^4-\sum_i(B^2_{ii})^2\).

3. **Aut-invariant reduction.** \(\delta\in E_{4p}^{\mathrm{Aut}}\); residual is finite-dimensional on double-coset orbits.

4. **PF+rank obstruction.** For \(p\ge5\), \(\lambda_2(P\odot P)<d/N\) strictly (nonnegative PopP, bulk sum \(=(d-1)\cdot(d/N)\), full bulk rank \(\binom{d-1}{2}\) forbids mult\(\ge d-1\) at the PF ceiling).

5. **Scalar form.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\) then \(\sum\rho^2=\|\rho_{\min}\|^2+c^2\), so residual \(\Leftrightarrow c^2\le\mathrm{room}_{\mathrm{hyp}}/24\).

**Certified p=5:** \(\dim E_{4p}^{\mathrm{Aut}}=1\), \(\delta=cv_0\), \(c^2=\mathrm{room}_{\mathrm{hyp}}/24\) (equality in \(T_\rho\)).

**OPEN:** \(c^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15109.py`, `evidence/e1_gmin_m4_prop15109.json`, `tests/test_prop15109.py`.

## Prop 15.110 (2026-07-31) — Closed Max+ identities; ρ_min²<budget for p≥7

**Proved:**

1. **∑κ∏ identity.** For boolean \(y\) with \(Cy=py\):
   \(\sum_S\kappa(S)\prod y_i=n(n-1)(n-2)/8\).
   Proof via \((y^\top Cy)^2=p^2n^2=2n(n-1)+8\sum\kappa\prod\) (case \(|3|=0\) from \(C^2=p^2I\)).

2. **e₄ constant.** On Max+, \(|\sum y_i|=p+1\), hence
   \(e_4=-p(p-1)(p+1)(p+4)/12\) by Newton (boolean).

3. **⟨m₄,κ⟩.** \(\sum_S m_4\kappa=n(n-1)(n-2)/8\).

4. **ρ_min² < room_hyp/24 for all primes p≥7.**
   \[
   \|\rho_{\min}\|_2^2-\frac{\mathrm{room}_{\mathrm{hyp}}}{24}
   =-\frac{(p^2-1)(3p^6-105p^4+37p^2-15)}{6p^2(p^2-5)(p^2+1)}<0.
   \]

5. **Sufficient criterion.** For \(p\ge7\): \(\delta^2\le\rho_{\min}^2\Rightarrow\) residual.
   At \(p=5\), equality \(\delta^2=\mathrm{room}_{\mathrm{hyp}}/24\) certified.
   \(c=Q_0(\mathrm{halfspace})\) with \(Q_0\) from Aut-invariant \(4p\)-eigenfunction of \(T\).

**Certified:** \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\); Theorem 4 for primes \(7\le p\le97\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\) (or closed Gauss-sum form of \(c\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15110.py`, `evidence/e1_gmin_m4_prop15110.json`, `tests/test_prop15110.py`.

## Prop 15.111 (2026-07-31) — Pair Schur identity; closed α_κ, α_ρ; Φ residual = 8⟨δ,κ_B⟩

**Proved:**

1. **Zero-diag pairing identity.** For real symmetric zero-diagonal \(C,B\):
   \[
   \sum_S \kappa_C(S)\kappa_B(S)
   =\tfrac14\mathrm{Tr}(CBCB)+\tfrac18(\mathrm{Tr}\,CB)^2
   +\tfrac12\sum_{i,j}C_{ij}^2 B_{ij}^2
   -\tfrac12\sum_i(CB)_{ii}^2-\tfrac12\sum_i(BC)_{ii}^2.
   \]

2. **α_κ on Z.** On \(Z\) (\(CB=pB=BC\), ambient diag\(B=0\)):
   \(\langle\kappa/p^2,\kappa_B\rangle=\alpha_\kappa\|B\|_F^2\) with
   \(\alpha_\kappa=(p^2+2)/(4p^2)\).

3. **Pair target.** \(\mathrm{pair}=(\bar\mu-6)/8=(p^2+11)/(4(p^2-5))\).

4. **Closed α_ρ.** \(\alpha_\rho=\mathrm{pair}-\alpha_\kappa=(7p^2+5)/(2p^2(p^2-5))\).
   Channel form: \(\rho_{\min}=(4p\,b+Tb)/\mathrm{den}\), \(\mathrm{den}=12(p^2-5)\),
   with \(\langle b,\kappa_B\rangle=(6/p)\|B\|^2\) and
   \(\langle Tb,\kappa_B\rangle=6(3p^2+5)/p^2\|B\|^2\) (Schur; certified p=3,5,7)
   recovers \(\alpha_\rho\).

5. **Φ residual is pure δ.** For unit \(B\in Z\):
   \[
   \mathbb E[(y^\top By)^2]=\bar\mu\,\|B\|_F^2+8\langle\delta,\kappa_B\rangle.
   \]
   Hence \(16N\Leftrightarrow\max\langle\delta,\kappa_B\rangle\le(n-10)/(n-6)\).
   The particular solution \(\rho_{\min}\) is absorbed into the flat bulk \(\bar\mu\).

**Certified:** Schur scalarity of \(\kappa/p^2\), \(b\), \(Tb\), \(\rho_{\min}\) on \(Z\) at \(p=3,5,7\) with matching closed forms; zero-diag identity on random pairs \(n=6..10\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) (or \(\max\langle\delta,\kappa_B\rangle\le(n-10)/(n-6)\), or \(c^2\le\mathrm{room}_{\mathrm{hyp}}/24\)) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15111.py`, `evidence/e1_gmin_m4_prop15111.json`, `tests/test_prop15111.py`.

## Prop 15.112 (2026-07-31) — Design moments; conference ‖κ‖²; ED4 residual dictionary

**Proved:**

1. **Conference ‖κ‖².** For conference \(C\) of order \(n=p^2+1\),
   \(\|\kappa\|_2^2=(n p^4/8)(n-6)+n(n-1)/2\) (zero-diag pairing at \(B=C\)).

2. **Antipodality.** \(Cy=py\) boolean \(\Rightarrow C(-y)=p(-y)\).

3. **Design moments.** If \(E[yy^\top]=2P_+\), then for fixed \(y_0\in\mathrm{Max+}\):
   \(E[D]=0\), \(E[D^2]=2n\). Antipodality + 2-design \(\Rightarrow\) spherical 3-design.

4. **ED4 dictionary.**
   \[
   \mathrm{ED4}=\mathrm{ED4}_{\mathrm{flat}}+24\delta^2,
   \qquad
   \delta^2\le\rho_{\min}^2\iff\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}
   \]
   with \(\mathrm{ED4}_{\mathrm{suf}}=\mathrm{ED4}_{\mathrm{flat}}+24\rho_{\min}^2\).
   For \(p\ge7\): \(\mathrm{ED4}_{\mathrm{suf}}<\mathrm{ED4}_{\mathrm{bud}}\).

**Certified:** \(E[yy^\top]=2P_+\), \(E[D^2]=2n\), \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\).

**Attack note:** class_key is not \(m_4\)-equitable at \(p=7\) (do not F19-thrash).

**OPEN:** \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15112.py`, `evidence/e1_gmin_m4_prop15112.json`, `tests/test_prop15112.py`.

## Prop 15.113 (2026-08-01) — ⟨f_y,Tκ⟩; ED4 via W; Q_δ criterion

**Proved:**

1. **⟨κ,Tκ⟩=0** (conference).
2. **⟨f_y,Tκ⟩=2p(p⁴−1)** for every boolean \(Cy=py\) (from ⟨ρ,b⟩=⟨ρ_min,b⟩ and constancy).
3. **ED4 via W:** \(W=\sum_{i<j}y_iy_jz_iz_j\), \(E[W]=n/2\), \(E[D^4]=3n^2+4E[W^2]\).
4. **Criterion:** if \(Q_\delta(y)\le\rho_{\min}^2\) for all \(y\in\mathrm{Max+}\), then \(\delta^2\le\rho_{\min}^2\).

**Certified** at \(p=3,5,7\). **OPEN:** \(Q_\delta\le\rho_{\min}^2\) or ED4≤ED4_suf for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15113.py`, `evidence/e1_gmin_m4_prop15113.json`, `tests/test_prop15113.py`.

## Prop 15.114 (2026-08-01) — γ-calculus; Tf_y multiplicative formula; ∑γ, ∑γ²

**Proved:**

1. **Multiplicative eigenformula.** For boolean \(Cy=py\),
   \((Tf_y)(S)=(4p-2\gamma_y(S))f_y(S)\) where
   \(\gamma_y(S)=\sum_{\{i,j\}\subset S}C_{ij}y_iy_j\).
   Equivalently \(Af_y=2(\gamma_y\odot f_y)\), so \(\langle\delta,\gamma_y\odot f_y\rangle=0\) for all \(\delta\in\ker A\).
2. **∑γ closed form.** \(\sum_S\gamma_y(S)=(6/p)\binom{n}{4}\) (edge double-count + \(y^\top Cy=pn\)).
3. **∑γ² closed form.** \(\sum_S\gamma_y(S)^2=6\binom{n}{4}+n(n-1)(n-2)/4\).
   Adjacent-edge cross terms vanish because \(Cy=py\) and \(n-1=p^2\); matching contribution is \(n(n-1)(n-2)/4\) via \(\sum\kappa\prod=n(n-1)(n-2)/8\).
4. **‖Tf_y‖² closed.** \(\|Tf_y\|_2^2=(16p^2-72)\binom{n}{4}+n(n-1)(n-2)\).
5. **Pair-average residual moment.** The ED4 dictionary uses \(N^{-2}\sum_{y,z}(y\cdot z)^4=E_y[\mathrm{ED4}(y)]\), not a single basepoint.

**Certified:** Thm 1–4 at \(p=3,5\); at \(p=7\), Max+ has **three** ED4(y) types (counts 2352, 8400, 700), all \(\le\mathrm{ED4}_{\mathrm{suf}}\), with global mean \(\le\mathrm{ED4}_{\mathrm{suf}}\) (not 2-point homogeneous).

**OPEN:** pair-average \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) for general \(p\ge5\). Attack: γ-weighted spectral mass of \(f_y\) on \(E_{4p}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15114.py`, `evidence/e1_gmin_m4_prop15114.json`, `tests/test_prop15114.py`.

## Prop 15.115 (2026-08-01) — Max+ residual solves \(A\rho=b\); \(\delta=P_{E_{4p}}m_4\); spectral moments of \(f_y\)

**Proved:**

1. **\(E[\gamma\odot f]=2\kappa/p\).** By the Max+ 2-design \(E[y_iy_j]=C_{ij}/p\) and double-counting matchings on each \(K_4\).
2. **Resolvent.** \(Am_4=4\kappa/p\) and \(A\rho=b\) for the Max+ residual \(\rho=m_4-\kappa/p^2\).
3. **\(\kappa\perp E_{4p}\).** \(T\kappa\in E_\mu\oplus E_{-\mu}\) with \(\mu^2=4(p^2+15)\); \((4p)^2-\mu^2=12(p^2-5)\ne0\) for primes \(p\ge3\). Hence \(\delta=P_{E_{4p}}m_4=E_y P_{E_{4p}}f_y\).
4. **Spectral moments of \(f_y\).** \(m_1=4p-12/p\), \(m_2=16p^2-72+24/(p^2-2)\),
   \(\mathrm{Var}=24(p^2-3)(p^2-4)/(p^2(p^2-2))\).
5. **Jensen.** \(\delta^2\le E_y\|P_{E_{4p}}f_y\|_2^2\) (full \(E_{4p}\) energy is too crude: at \(p=5\) it is \(\gg\rho_{\min}^2\) while residual still holds).

**Certified:** resolvent + moments at \(p=3,5\); \(Q_\delta\) constant (\(=\delta^2\)) at \(p=3,5\); ED4≤ED4_suf at \(p=3,5,7\).

**OPEN:** \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\) via the **coherent** mass \(\|E_y P_{E_{4p}}f_y\|_2\), not full \(E_{4p}\) energy. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15115.py`, `evidence/e1_gmin_m4_prop15115.json`, `tests/test_prop15115.py`.

## Prop 15.116 (2026-08-01) — e₄↔ED4↔δ dictionary; coherent mass; Aut-line criterion

**Proved:**

1. **e₄ poly.** Boolean \(e_4(s)=s^4/24+((4-3n)/12)s^2+n(n-2)/8\).
2. **∑m₄² from ED4.** \(\sum m_4^2=\mathrm{ED4}/24+n(4-3n)/6+n(n-2)/8\).
3. **⟨κ,ρ_min⟩.** \(n(n-1)(n-2)(n-6)/(2p^2(p^2-5))\).
4. **Flat identity.** Pythagoras flat part matches the e₄ constant for all primes \(p\ge3\), \(p^2\neq5\).
5. **Coherent mass.** \(\delta^2=\|E_y P_{E_{4p}}f_y\|_2^2\).
6. **Aut-line.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\) then \(\delta=cv_0\) and residual \(\Leftrightarrow c^2\le\rho_{\min}^2\); when \(Q_0\) constant, \(c=Q_0(x_{\mathrm{hs}})\) (halfspace, Max+-free).
7. **Min-distance envelope.** \(|D|\le p^2-2p-1\) for \(y\neq\pm z\) yields an ED4 UB that **fails** \(\mathrm{ED4}_{\mathrm{suf}}\) for all primes \(p\ge5\) (dead for closing).

**Certified:** dictionary at \(p=3,5,7\); Aut-line at \(p=3,5\); \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\).

**OPEN:** coherent mass \(\delta^2\le\rho_{\min}^2\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15116.py`, `evidence/e1_gmin_m4_prop15116.json`, `tests/test_prop15116.py`.

## Prop 15.117 (2026-08-01) — Path C hyp residual primary; ρ_min pairings

**Proved:**

1. **Path C primary residual.** For all primes \(p\ge5\),
   \[
   \delta^2\le\frac{\mathrm{room}_{\mathrm{hyp}}}{24}
   =\frac{4(p^2-9)(p^2-1)^2}{3(p^2-5)(p^2+1)}
   \]
   \(\Leftrightarrow\) \(\mathrm{orth}\le\mathrm{room}_{\mathrm{hyp}}\) \(\Leftrightarrow\) \(\|\kappa\|_F^2\le\kappa_{\mathrm{hyp}}\) \(\Leftrightarrow\) \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\).
2. **ρ_min² vs hyp.** \(\rho_{\min}^2>\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=5\); \(\rho_{\min}^2<\mathrm{room}_{\mathrm{hyp}}/24\) for all primes \(p\ge7\). Thus \(\delta^2\le\rho_{\min}^2\) is sufficient for Path C when \(p\ge7\); at \(p=5\) the hyp form is the tight target.
3. **Slack.** \(\kappa_{\mathrm{hyp}}-\|\kappa\|_F^2=24(\mathrm{room}_{\mathrm{hyp}}/24-\delta^2)\).
4. **Pairings.** \(\langle b,f_y\rangle=\langle\rho_{\min},b\rangle=2(p^4-1)/p\) for every \(y\in\mathrm{Max+}\); \(\langle\rho_{\min},m_4\rangle=\rho_{\min}^2+\langle\kappa,\rho_{\min}\rangle/p^2\) closed.
5. **γ-channel average.** \(E_y\langle b,\gamma_y\odot f_y\rangle=0\). If pointwise zero, then \(\langle\rho_{\min},f_y\rangle=4(p^4-1)/(3(p^2-5))\) constant on Max+.
6. **Coherent mass.** Path C residual \(\Leftrightarrow\|E_y P_{E_{4p}}f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\).

**Certified:** pointwise \(\langle b,\gamma\odot f\rangle=0\) at \(p=3\); hyp residual \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=3\) (eq), and at \(p=5,7\) when Max+ caches are present (eq at 5, strict at 7).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15117.py`, `evidence/e1_gmin_m4_prop15117.json`, `tests/test_prop15117.py`.

## Prop 15.118 (2026-08-01) — Pointwise ⟨b,γ⊙f⟩=0; T²κ pairing; ρ_min·f closed

**Proved:**

1. **⟨κ,γ⊙f⟩ closed.** For every Max+ \(y\),
   \(\sum_S\kappa(S)\gamma_y(S)f_y(S)=p(p^2+1)(p^2-1)(p^2-4)/4\)
   (edge expansion + \(Cy=py\) + \(C^2=p^2I\)).
2. **⟨T²κ,m₄⟩.** With \(T\kappa\in E_\mu\oplus E_{-\mu}\) and \(\delta\perp T^2\kappa\),
   \(\langle T^2\kappa,m_4\rangle=8p^2(p^4-1)\).
3. **Pointwise criterion.** For each \(y\in\mathrm{Max+}\):
   \(\langle b,\gamma\odot f\rangle=0\Leftrightarrow\langle Tb,f\rangle=4p\langle b,f\rangle\Leftrightarrow\langle T^2\kappa,f\rangle=8p^2(p^4-1)\).
4. **Pointwise vanishing.** \(\langle b,\gamma_y\odot f_y\rangle=0\) and \(\langle T^2\kappa,f_y\rangle=8p^2(p^4-1)\) for all Max+ \(y\) (constancy of the \(T\)-module pairing on Max+ + mean from (2)).
5. **⟨ρ_min,f_y⟩ closed.** \(\langle\rho_{\min},f_y\rangle=4(p^4-1)/(3(p^2-5))\) for all Max+ \(y\).

**Certified:** full Max+ census at \(p=3\) (\(N=12\)) and \(p=5\) (\(N=260\)): bgf=0, \(T^2\) target, ρ_min pairing, and \(\delta^2=\mathrm{room}_{\mathrm{hyp}}/24\) at \(p=5\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (pin \(\langle m_4,f_{\mathrm{hs}}\rangle\) / ED4). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15118.py`, `evidence/e1_gmin_m4_prop15118.json`, `tests/test_prop15118.py`.

## Prop 15.119 (2026-08-01) — Residual budget dictionary; weight enum; halfspace pin

**Proved:**

1. **ED4 budgets (closed rationals).** For primes \(p>\sqrt5\),
   \[
   \mathrm{ED4}_{\mathrm{flat}}=\frac{4(p^2-3)(p^2+1)(3p^2+1)}{p^2-5},\qquad
   \mathrm{ED4}_{\mathrm{bud}}=\frac{4(3p^8+6p^6-104p^4+138p^2-75)}{(p^2-5)(p^2+1)}.
   \]
2. **E[W²] channel.** With \(D=y\cdot z\), \(D^2=n+2W\), \(\mathrm{ED4}=3n^2+4\mathbb E[W^2]\):
   \[
   \mathrm{EW2}_{\mathrm{flat}}=\frac{(p^2+1)(9p^4-20p^2+3)}{4(p^2-5)},\qquad
   \mathrm{EW2}_{\mathrm{bud}}=\frac{9p^8+30p^6-380p^4+594p^2-285}{4(p^2-5)(p^2+1)}.
   \]
   Path C residual \(\Leftrightarrow\mathbb E[W^2]\le\mathrm{EW2}_{\mathrm{bud}}\). Gap:
   \(\mathrm{EW2}_{\mathrm{bud}}-\mathrm{EW2}_{\mathrm{flat}}=\mathrm{room}_{\mathrm{hyp}}/4=6\cdot(\mathrm{room}_{\mathrm{hyp}}/24)\).
3. **⟨m₄,f_y⟩ ↔ ED4.** \(\langle m_4,f_y\rangle=\mathrm{ED4}(y)/24+n(4-3n)/6+n(n-2)/8\). Budget forms:
   \[
   \mathrm{m4f}_{\mathrm{flat}}=\frac{(p-1)(p+1)(p^2+1)(3p^2+17)}{24(p^2-5)},\qquad
   \mathrm{m4f}_{\mathrm{bud}}=\frac{(p-1)(p+1)(3p^2-5)(p^4+20p^2-61)}{24(p^2-5)(p^2+1)}.
   \]
4. **Equivalence.** For \(p\ge5\): \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\Leftrightarrow\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\Leftrightarrow\mathbb E[W^2]\le\mathrm{EW2}_{\mathrm{bud}}\). When \(Q_\delta\) is constant on Max+ (certified \(p=3,5\)): also \(\Leftrightarrow\langle m_4,f_y\rangle\le\mathrm{m4f}_{\mathrm{bud}}\) (halfspace pin).
5. **Weight-enumerator structure.** Max+ dots satisfy \(y\cdot z\equiv2\pmod4\), \(|y\cdot z|\le p^2-2p-1\) off \(\pm\) pairs, antipodal measure. Crude envelope \(\mathrm{ED4}\le2n D_{\max}^2\) is strictly larger than \(\mathrm{ED4}_{\mathrm{bud}}\) for \(p\ge5\) (too weak).

**Certified:** full Max+ weight spectra at \(p=3\) (\(\{\pm10,\pm2\}\)) and \(p=5\) (\(\{\pm26,\pm14,\pm10,\pm6,\pm2\}\) with mults \(\{1,13,20,36,60\}\)); \(\mathrm{ED4}=\mathrm{ED4}_{\mathrm{bud}}\) and \(\langle m_4,f_{\mathrm{hs}}\rangle=\mathrm{m4f}_{\mathrm{bud}}\) at both; \(\mathrm{orth}\cdot N=147456\) at \(p=5\). Prior \(p=7\): \(\delta^2/(\mathrm{room}_{\mathrm{hyp}}/24)=124875/669124\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (independent upper bound on \(\mathbb E[W^2]\) or \(\langle m_4,f_{\mathrm{hs}}\rangle\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15119.py`, `evidence/e1_gmin_m4_prop15119.json`, `tests/test_prop15119.py`.

## Prop 15.120 (2026-08-01) — Pointwise E[W²] factorization; Pythagoras; majorization

**Proved:**

1. **Pointwise residual factorization.** For every \(y\in\mathrm{Max+}\),
   \(\mathbb E_z[W_y(z)^2]=\mathrm{EW2}_{\mathrm{flat}}+6\,Q_\delta(y)\).
   Wick wedges + \(\kappa/p^2+\rho_{\min}\) on disjoint edges give the flat part;
   the \(\delta\)-part of disjoint-edge pairings contributes exactly \(6\langle\delta,f_y\rangle\).
2. **Pythagoras.** \(\sum m_4^2=F(p)+\delta^2\) with
   \(F=\|\kappa/p^2+\rho_{\min}\|_2^2=\mathrm{m4f}_{\mathrm{flat}}\) Max+-free.
3. **Majorization UB.** \(H=G\odot G\succeq0\), \(\lambda_{\max}=2Nn\), \(\mathrm{Tr}=Nn^2\), \(\mathrm{Tr}/\lambda_{\max}=d\)
   \(\Rightarrow\mathrm{ED4}\le 2n^3\). For all primes \(p\ge5\), \(2n^3>\mathrm{ED4}_{\mathrm{bud}}\) (too weak).
4. **Dead independent UBs.** CS via \(\|\tilde\gamma\odot f\|\), discrete moment LP \(2n D_{\max}^2\),
   and majorization \(2n^3\) all exceed residual budgets for \(p\ge5\).

**Certified:** EW2 constant on Max+ and equal to \(\mathrm{EW2}_{\mathrm{flat}}+6\delta^2\) at \(p=3,5\).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (weight enumerator / Gauss sums). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15120.py`, `evidence/e1_gmin_m4_prop15120.json`, `tests/test_prop15120.py`.

## Prop 15.121 (2026-08-01) — Spectral residual dictionary; Frobenius form

**Proved:**

1. **ED4 from FFT.** \(M_{ab}=(y_a\cdot y_b)^2=nJ+2\,\mathrm{FFT}\) (15.93) \(\Rightarrow\)
   \(\mathrm{ED4}=4n^2+4N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\).
2. **E[W²] from FFT.** \(\mathbb E[W^2]=\|\mathrm{FFT}\|_F^2/N^2=d^2+N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\),
   and with 15.120: \(\|\mathrm{FFT}|_{1^\perp}\|_F^2=N^2(\mathrm{EW2}_{\mathrm{flat}}-d^2+6\delta^2)\).
3. **Φ variance.** \(\sum\lambda_\alpha^2=\mathrm{ED4}-4n^2\), \(\sum(\lambda_\alpha-\bar\mu)^2=\mathrm{orth}=24\delta^2\) (15.105);
   at \(\delta=0\), \(\Phi\equiv\bar\mu\) on \(Z\).
4. **EW2_flat closed.** \(\mathrm{EW2}_{\mathrm{flat}}=(n^2+T^2/m)/4\) with \(T=n(n-2)\), \(m=\dim Z\).
5. **Residual ⇔ Frobenius.** Path C residual \(\Leftrightarrow\|\mathrm{FFT}|_{1^\perp}\|_F^2\le N^2(\mathrm{EW2}_{\mathrm{bud}}-d^2)\)
   \(\Leftrightarrow\sum(\lambda_\alpha-\bar\mu)^2\le\mathrm{room}_{\mathrm{hyp}}\).
   Contrast: 16N \(\Leftrightarrow\|\mathrm{FFT}|_{1^\perp}\|_{\mathrm{op}}\le8N\) (operator norm on the same operator).
6. **H/16N·Tr too weak.** \(\|A\|_F^2\le\|A\|_{\mathrm{op}}\mathrm{Tr}(A)\) under H or 16N yields ED4 UBs \(> \mathrm{ED4}_{\mathrm{bud}}\) for all \(p\ge5\).

**Certified:** spectral identities at \(p=3,5\) (full Max+ Gram/FFT).

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15121.py`, `evidence/e1_gmin_m4_prop15121.json`, `tests/test_prop15121.py`.

## Prop 15.122 (2026-08-01) — Max+ disagreement \(u\in V_+\); Aut-line; \(\lambda_{\max}(T)\)

**Proved:**

1. **Disagreement identity.** For \(y,z\in\mathrm{Max}+\), \(D=y\cdot z\), \(k=(n-D)/2\),
   \(u=(y-z)/2\in\{0,\pm1\}^n\) has weight \(k\) and
   \(u^\top Cu=pk\) with \(u\in V_+\) (equivalently \(\|P_+u\|_2^2=k\)).
   Same for \(v=(y+z)/2\) with weight \(n-k\).
2. **Dot support.** Off \(\pm\) pairs, \(|D|\le p^2-2p-1\), \(D\equiv2\pmod4\), and \(k\) is a
   ternary \(V_+\) weight. Cert spectra \(p=3,5\).
3. **\(\lambda_{\max}(T)\) threshold.** \(\lambda_{\max}(T)<4p\Rightarrow\delta=0\) (residual OK).
   Cert: \(p=3\) strict; \(p=5,7\) equality \(\lambda_{\max}=4p\).
4. **Aut-line.** \(\delta\in E_{4p}^{\mathrm{Aut}}\); if \(\dim\le1\) then residual \(\Leftrightarrow c^2\le\mathrm{room}_{\mathrm{hyp}}/24\)
   with \(c=Q_0(x_{\mathrm{hs}})\). Cert equality line at \(p=5\); \(p=7\) has non-constant \(Q_\delta\).
5. **FFT budget.** Residual \(\Leftrightarrow N^{-2}\|\mathrm{FFT}|_{1^\perp}\|_F^2\le B(p)=\mathrm{EW2}_{\mathrm{bud}}-d^2\) (Max+-free).
6. **Dead tight attempts.** Discrete LP with exact \(N\), PGL character sums (F18), and prior majorizations all exceed residual budgets.

**OPEN:** \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for general \(p\ge5\) (ternary \(V_+\) weight enumerator / Gauss \(Q_0\)). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15122.py`, `evidence/e1_gmin_m4_prop15122.json`, `tests/test_prop15122.py`.

## Prop 15.123 (2026-08-01) — Switching; conference srg; dual Krawtchouk residual

**Proved:**

1. **Switching bijection.** After \(C'=D_yCD_y\), \(1\in\mathrm{Max}+'\) and
   \(z\mapsto w=(1-D_yz)/2\) bijects \(\mathrm{Max}+'\leftrightarrow V_+\cap\{0,1\}^n\), with
   \(D(y,z)=n-2\,\mathrm{wt}(w)\).
2. **Conference srg.** \(G\) with edges \(C'=-1\) is
   \(\mathrm{srg}(n,\,p(p-1)/2,\,\mu-1,\,\mu)\), \(\mu=((p-1)/2)^2\);
   \(A\)-eigs \(k_G\), \(\theta=(p-1)/2\) (\(\times d\)), \(\tau=-(p+1)/2\) (\(\times d-1\));
   \(C'\)-eigs \(\pm p\) each of mult \(d\).
3. **Regular sets.** Supports of weight-\(k\) codewords are regular sets with
   \(\alpha=(k-1-p)/2\), \(\beta=k/2\), \(\alpha-\beta=\tau\); allowed
   \(k\in\{0,n\}\cup\mathrm{even}[p+1,p(p-1)]\).
4. **Weight = distance dist.** \(B_i=W_i\), \(|X|=N\), \(\mathrm{ED4}=N^{-1}\sum W_k(n-2k)^4\).
5. **Dual residual.** \(A'_4=\sum m_4^2=\mathrm{m4f}_{\mathrm{flat}}+\delta^2\);
   residual \(\Leftrightarrow A'_4\le\mathrm{m4f}_{\mathrm{bud}}\). (Hamming Delsarte alone too weak.)
6. **Two-valued form.** \(\chi_S-(k/n)\mathbf1\in V_+\cap\mathbf1^\perp\) with only two coordinate values
   and \(\|\cdot\|_2^2=k(n-k)/n\); \(W_k\) counts such vectors.

**Certified:** srg params and full \(W_k\) at \(p=3,5\) (Petersen; srg(26,10,3,4)).

**OPEN:** closed \(W_k\) / \(A'_4\le\mathrm{m4f}_{\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15123.py`, `evidence/e1_gmin_m4_prop15123.json`, `tests/test_prop15123.py`.

## Prop 15.124 (2026-08-01) — Closed weight moments \(j\le3\); \(E[k^4]\) partition; residual as \(R_4\)

Continues 15.123. Isolates the residual as a single 4-wise mass. Does **not** soft-close.

**Proved:**

1. **Moments \(j\le3\) (Max+-free).** After switch, \(E[s^2]=2n\), \(E[s]=0\), \(E[s^3]=0\), and
   \[
   E[k]=\tfrac n2,\quad
   E[k^2]=\tfrac{n(n+2)}{4},\quad
   E[k^3]=\tfrac{n^2(n+6)}{8}.
   \]
2. **Exact \(\le3\) partition of \(E[k^4]\).** Using only pair/triple design averages,
   \[
   \mathrm{exact}_{\le3}
   =\tfrac n2+n^2+\tfrac{3n^2}{4}+\tfrac{3n(n-2)(n+2)}{4},
   \]
   and \(R_4:=E[k^4]-\mathrm{exact}_{\le3}=n(n-1)(n-2)(n-3)\,\overline E[\prod_4 w]\).
3. **Residual dictionary.** \(E[D^4]=\mathrm{ed4\_from\_exact3}(p)+16 R_4\); Path C residual
   \(\Leftrightarrow R_4\le R_{4,\mathrm{bud}}\Leftrightarrow A'_4\le\mathrm{m4f}_{\mathrm{bud}}\).
4. **Hamming Delsarte LP.** Max \(A'_4\) under allowed-weight dual constraints saturates
   \(\mathrm{m4f}_{\mathrm{bud}}\) at \(p=3\) and **strictly exceeds** it at \(p=5,7\) (too weak for \(p\ge5\)).
5. **Hoffman layer.** \(W_{p+1}=d\) (regular cocliques) certified at \(p=3,5\).

**Certified:** moments + partition + \(R_4=\mathrm{bud}\) saturation + \(A'_4=\mathrm{m4f}_{\mathrm{bud}}\) at \(p=3,5\).

**OPEN:** closed \(W_k\) / \(A'_4\le\mathrm{m4f}_{\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15124.py`, `evidence/e1_gmin_m4_prop15124.json`, `tests/test_prop15124.py`.

## Prop 15.125 (2026-08-01) — Perfect 2-colorings; 4-design defect; closed \(R_4\) budget

Continues 15.123–15.124. Does **not** soft-close.

**Proved:**

1. **Perfect 2-colorings.** \(W_k\) equals the number of \(\tau\)-equitable bipartitions
   (perfect 2-colorings) of the conference srg with \(|S|=k\),
   \(\alpha=(k-1-p)/2\), \(\beta=k/2\), \(\alpha-\beta=\tau=-(p+1)/2\).
   Hoffman layer \(k=p+1\) has \(\alpha=0\).
2. **Spherical 2-design / 4-design defect.** After \(V_+\cong\mathbb R^d\), Max+ is a
   spherical 2-design (\(E[uu^\top]=I_d/d\)). It is not a 4-design for \(p\ge5\):
   \(E[s^4]=\mathrm{ED4}\) exceeds \(3n^4/(d(d+2))\). Path C residual is this defect
   controlled by \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{bud}}\).
3. **Closed \(R_4\) budget.**
   \(\mathrm{ed4\_from\_exact3}=-(p^2+1)(p^6+3p^4-25p^2+13)\);
   \(R_{4,\mathrm{bud}}=(\mathrm{ED4}_{\mathrm{bud}}-\mathrm{ed4\_from\_exact3})/16\).
4. **Delsarte + moments \(j\le3\).** Still saturates only at \(p=3\); weak for \(p\ge5\).
5. **Antipodal dual.** \(W_k=W_{n-k}\Rightarrow A'_j=A'_{n-j}\).

**Certified:** algebra + LP weak p=5,7; defect ratios and saturation at \(p=3,5\).

**OPEN:** closed \(W_k\) / \(R_4\le R_{4,\mathrm{bud}}\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15125.py`, `evidence/e1_gmin_m4_prop15125.json`, `tests/test_prop15125.py`.

## Prop 15.126 (2026-08-01) — Geometric Hoffman seed; 1-design; simplex bound

Continues 15.123–15.125. Does **not** soft-close.

**Proved:**

1. **Geometric seed.** The subfield line \(F_p\cup\{\infty\}\subset\mathrm{PG}(1,p^2)\) is a
   Hoffman coclique (\(\tau\)-equitable, \(\alpha=0\), size \(p+1\)) after the standard
   halfspace Seidel switch. Certified \(p=3,5,7\).
2. **1-design algebra.** Hoffman cocliques form an Aut-invariant 1-design:
   \(b(p+1)=nr\), and \(r=(p+1)/2\Longleftrightarrow b=d=n/2\).
3. **Census equality.** At \(p=3,5\): \(W_{p+1}=d\), \(r=(p+1)/2\), and
   \(\{\chi_S\}\) is a basis of \(V_+\).
4. **Simplex bound.** Equal pairwise intersections \(\Rightarrow W_{p+1}\le d\)
   (Gram rank \(b-1\le d-1\)); equality at \(p=3\) (\(\lambda=1\)).
5. **ED4 contribution.** If \(W_{p+1}=W_{n-p-1}=d\), Hoffman layers contribute
   \(2d(p^2-2p-1)^4/N\) to ED4.

**OPEN:** \(W_{p+1}=d\) for all primes \(p\ge5\); full \(W_k\) / 4-design defect bound.
L remains OPEN. **(W=d general disproved in Prop 15.127.)**

Evidence: `src/e1_gmin_m4_prop15126.py`, `evidence/e1_gmin_m4_prop15126.json`, `tests/test_prop15126.py`.

## Prop 15.127 (2026-08-01) — Closed \(W_{p+1}\); inversive plane; \(W=d\) false

Continues 15.126. Does **not** soft-close residual.

**Proved:**

1. **Inversive plane.** \(F_p\)-sublines of \(\mathrm{PG}(1,p^2)\) form the miquelian
   inversive plane of order \(p\): \(S(3,p+1,p^2+1)\), \(b=p(p^2+1)\), \(\lambda_2=p+1\), \(\lambda_3=1\).
2. **Closed Hoffman weight.**
   \[
   W_{p+1}
   =\frac{1+\chi_4(p)}{2}\cdot\frac{p^2+1}{2}
   +\frac{1-\chi_4(p)}{2}\cdot\frac{3p+1}{2},
   \quad
   \chi_4(p)=(-1)^{(p-1)/2}.
   \]
   Equivalently: \(W_{p+1}=d\) if \(p\equiv1\pmod4\), and \(W_{p+1}=(3p+1)/2\) if \(p\equiv3\pmod4\).
3. **Counterexample to \(W_{p+1}=d\).** At \(p=7\), \(W_8=11\neq25=d\).
4. **Census.** Full max-coclique enumeration: \(W=(5,13,11,17)\) at \(p=3,5,7,11\), matching the formula; coincides with the count of regular sublines.
5. **Corrected ED4 part.** Hoffman contribution \(2W_{p+1}D_{\max}^4/N\) (not always \(2d\)).

**OPEN:** full closed \(W_k\) / 4-design defect bound for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15127.py`, `evidence/e1_gmin_m4_prop15127.json`, `tests/test_prop15127.py`.

## Prop 15.128 (2026-08-01) — Full \(W_k\) census \(p=3,5,7\); exact ED4 at \(p=7\)

Continues 15.127. Does **not** soft-close residual.

**Certified:**

1. **Full weight enumerators.**
   - \(p=3\) (\(N=12\)): \(\{0{:}1,4{:}5,6{:}5,10{:}1\}\)
   - \(p=5\) (\(N=260\)): \(\{0{:}1,6{:}13,8{:}20,10{:}36,12{:}60,14{:}60,\ldots\}\)
   - \(p=7\) (\(N=11452\)), free-coord Max+ \(2^{25}\):
     \(\{0{:}1,8{:}11,12{:}112,14{:}159,16{:}280,18{:}728,20{:}1099,22{:}1502,24{:}1834,\ldots\}\)
     with \(W_{10}=W_{40}=0\) (allowed by \(\alpha\ge0\) but empty).
2. **Consistency.** \(\sum W=N\), \(W_k=W_{n-k}\), \(W_{p+1}\) matches 15.127,
   \(E[k^j]\) for \(j\le3\) match 15.124.
3. **ED4 at \(p=7\).** \(\mathrm{ED4}=12835984/409<\mathrm{ED4}_{\mathrm{bud}}=1775728/55\)
   (strict hyp residual); \(\delta^2=82176/4499\).
4. **Saturation.** \(\mathrm{ED4}=\mathrm{ED4}_{\mathrm{bud}}\) at \(p=3,5\).

**OPEN:** closed \(W_k\) for general primes \(p\ge5\) (must allow structural zeros),
or character-sum/PBIBD bound on the 4-design defect. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15128.py`, `evidence/e1_gmin_m4_prop15128.json`, `tests/test_prop15128.py`.

## Prop 15.129 (2026-08-01) — Jensen coherent-mass inequality; Hoffman \(\bar r\)

Continues 15.128. Does **not** soft-close residual.

**Proved:**

1. **Jensen.** \(\delta=E_y[P_{E_{4p}}f_y]\) and
   \(\delta^2\le E\|P f_y\|_2^2\), equality iff \(P f_y\) is a.s. constant on Max+.
   Hence \(E\|P f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\Rightarrow\) Path C residual.
2. **Dictionary.** Linear equivalences ED4 \(\leftrightarrow\sum m_4^2\leftrightarrow\delta^2\)
   certified on full \(W\) at \(p=3,5,7\).
3. **Average replication.** \(\bar r=W_{p+1}(p+1)/n\) with \(W_{p+1}\) from 15.127.
   Integral (1-design possible) for \(p\equiv1\pmod4\) and \(p=3\); at \(p=7\),
   \(\bar r=44/25\notin\mathbb Z\), so the Hoffman layer is **not** a 1-design.
4. **Hoffman geometry at \(p=5\).** 30 disjoint Hoffman pairs, each giving a
   regular 12-set; these are exactly half of \(W_{12}\). Every weight-16 set
   contains exactly two Hoffman cocliques.

**OPEN:** bound \(E\|P_{E_{4p}}f_y\|_2^2\le\mathrm{room}_{\mathrm{hyp}}/24\) for all primes
\(p\ge5\), or closed general \(W_k\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15129.py`, `evidence/e1_gmin_m4_prop15129.json`, `tests/test_prop15129.py`.

## Prop 15.130 (2026-08-01) — \(P m_4=\delta\); \(\rho_{\min}\)-sufficient residual for \(p\ge7\)

Continues 15.129. Does **not** soft-close residual.

**Proved:**

1. **\(P m_4=\delta\).** Since \(m_4=\kappa/p^2+\rho_{\min}+\delta\) with the first two summands
   orthogonal to \(E_{4p}\), \(P m_4=\delta\). With \(m_4=E f_y\), \(\delta=E[P f_y]\).
2. **Jensen.** \(\delta^2\le E\|P f_y\|_2^2\); bound on the right-hand side yields residual.
3. **Gap algebra.** For all primes \(p\ge7\),
   \[
   \frac{\mathrm{room}_{\mathrm{hyp}}}{24}-\rho_{\min}^2
   =\frac{(p-1)(p+1)(3p^6-105p^4+37p^2-15)}{6p^2(p^2-5)(p^2+1)}>0,
   \]
   so \(\delta^2\le\rho_{\min}^2\Rightarrow\) Path C residual. (At \(p=5\), \(\rho_{\min}^2>\mathrm{room}\); use hyp form.)
4. **Census.** \(\delta^2\le\rho_{\min}^2\) at \(p=3,5,7\) (ratio \(\approx0.379\) at \(p=7\)).
5. **Aut-line program.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\), residual \(\Leftrightarrow c^2\le\mathrm{room}/24\)
   with \(c=Q_0(\mathrm{hs})\) when \(Q_0\) is constant.

**OPEN:** \(\delta^2\le\rho_{\min}^2\) (or \(E\|Pf\|^2\le\mathrm{room}\), or Gauss-sum \(Q_0\)) for general
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15130.py`, `evidence/e1_gmin_m4_prop15130.json`, `tests/test_prop15130.py`.

## Prop 15.131 (2026-08-01) — Pair-avg vs basepoint ED4; p=7 three-type \(Q_\delta\)

Continues 15.130. Does **not** soft-close residual.

**Proved / certified:**

1. **Pair vs basepoint.** \(ED4(y)=ED4_{\mathrm{flat}}+24 Q_\delta(y)\) and
   \(ED4_{\mathrm{pair}}=ED4_{\mathrm{flat}}+24\delta^2\). When \(Q_\delta\) is non-constant,
   \(ED4(y_0)\ne ED4_{\mathrm{pair}}\) for a fixed basepoint \(y_0\) (e.g. halfspace);
   the weight-enumerator moment \(E[(n-2k)^4]\) equals \(ED4(y_0)\), not \(\delta^2\).
2. **p=7 spectrum (ProcessPool W=86).** Exactly three types:
   counts \(2352/8400/700\), \(Q_\delta\in\{-124800/4499,\ 82176/4499,\ 200448/4499\}\).
   \(Q_\delta\) can be negative (not 2-point homogeneous for 4th moment).
3. **True \(\delta^2\).** Pair-avg \(\delta^2=19180800/1840091\approx10.424\)
   (prior W-based \(82176/4499\) was \(Q_\delta(\mathrm{hs})\), the middle type).
   Ratios: \(\delta^2/\rho_{\min}^2\approx0.216\), \(\delta^2/\mathrm{room}\approx0.187\).
4. **Pointwise criterion.** \(\max Q_\delta=200448/4499\le\rho_{\min}^2=26000/539\)
   with slack \(812048/220451\); hence residual for \(p=7\) by 15.130.C.
5. **Variance.** \(\mathrm{Var}(Q_\delta)>0\) at \(p=7\), so pointwise is strictly
   stronger than \(\delta^2\le\rho_{\min}^2\).

**OPEN:** \(\max_y Q_\delta(y)\le\rho_{\min}^2\) (or \(\delta^2\le\rho_{\min}^2\)) for all primes
\(p\ge7\); at \(p=5\) use hyp form. Aut-line at \(p=7\) needs care (\(Q\) non-constant).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15131.py`, `evidence/e1_gmin_m4_prop15131.json`, `tests/test_prop15131.py`.

## Prop 15.132 (2026-08-01) — Max+-free residual dictionary; Aut \(\delta\); dead envelopes

Continues 15.131. Does **not** soft-close residual. No full Max+ census for general \(p\).

**Proved / certified:**

1. **Max+-free dictionary.** \(\delta^2=\sum m_4^2-m4f_{\mathrm{flat}}\) with
   \(m4f_{\mathrm{flat}}\) closed Max+-free. Residual
   \(\Leftrightarrow\sum m_4^2\le m4f_{\mathrm{bud}}:=m4f_{\mathrm{flat}}+\mathrm{room}/24\).
   For \(p\ge7\): \(\sum m_4^2\le m4f_{\mathrm{flat}}+\rho_{\min}^2\) suffices.
2. **Aut-invariance.** \(\delta\) is Aut-invariant; \(Q_\delta\) constant on Aut-orbits of Max+.
   At \(p=7\), three \(Q\) types \(\Rightarrow\ge3\) Aut-orbits (not transitive).
3. **\(\gamma\)-parity.** \(\gamma_y(S)\in\{-6,-4,-2,0,2,4,6\}\); formal \(4p\)-fiber is \(\gamma=0\).
4. **\(\gamma=0\) mass.** Constant \(4350\) at \(p=5\); 3-valued in \(p=7\) sample (parallels \(Q\) types).
5. **Dead envelopes.** Moment LP on allowed regular-set \(k\), and pole+\(D_{\max}\) mixes with
   any \(N\ge n\), all strictly exceed \(\mathrm{ED4}_{\mathrm{suf}}\) for primes \(p=5..19\).

**OPEN:** Max+-free bound \(\max Q_\delta\le\rho_{\min}^2\) or \(\delta^2\le\rho_{\min}^2\)
(character sums / Aut-orbit Bose–Mesner / \(\dim E_{4p}^{\mathrm{Aut}}\le1+Q_0\)) for all
primes \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15132.py`, `evidence/e1_gmin_m4_prop15132.json`, `tests/test_prop15132.py`.

## Prop 15.133 (2026-08-01) — class_key Bose–Mesner; F19 quantitative; Aut-line

Continues 15.132. Does **not** soft-close residual. Does **not** thrash class_key (F19).

**Proved / certified:**

1. **Aut-line form.** If \(\dim E_{4p}^{\mathrm{Aut}}\le1\), residual
   \(\Leftrightarrow c^2\le\mathrm{room}/24\) with \(c=Q_0(\mathrm{hs})\) when \(Q_0\) constant.
2. **class_key \(T\)-spectrum (ProcessPool W=86).**
   \(\dim E_{4p}^{\mathrm{ck}}=0,1,0\) at \(p=3,5,7\);
   \(\lambda_{\max}=4p\) only at \(p=5\) among these.
3. **F19 quantitative.** At \(p=7\), \(\dim E_{4p}^{\mathrm{ck}}=0\) but
   \(\delta^2=19180800/1840091>0\), so \(\delta\notin V^{\mathrm{ck}}\)
   (\(m_4\) not class_key-equitable). class_key cannot close residual at \(p=7\).
4. **CR dead.** PGL cross-ratio orbits are not \(\mathrm{Aut}(C)\) orbits
   (\(\kappa\) takes both signs inside CR cells at \(p=5\)); false nullity 0.
5. **Success locus.** Aut-line via class_key works at \(p=5\) (\(c^2=\mathrm{room}=1536/65\));
   ambient \(E_{4p}=0\) at \(p=3\); \(p=7\) needs true \(\mathrm{Aut}(C)\) or character sums.

**OPEN:** true \(\mathrm{Aut}(C)\) Bose–Mesner / Gauss \(Q_0\) / \(\max Q_\delta\le\rho_{\min}^2\)
for general primes. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15133.py`, `evidence/e1_gmin_m4_prop15133.json`,
`evidence/e1_gmin_m4_prop15133_classkey_spectrum.json`, `tests/test_prop15133.py`.

## Prop 15.134 (2026-08-01) — Strict Aut(\(C\)) Bose–Mesner; residual projection

Continues 15.133. Does **not** soft-close residual. No class_key / raw PGL thrash.

**Proved / certified:**

1. **Strict Aut \(G\).** Affine square-semilinear maps
   \(x\mapsto a\cdot\mathrm{Frob}^i(x)+b\) (\(\chi(a)=1\), \(\infty\) fixed) give
   \(|G|=p^2(p^2-1)\) and \(P^\top CP=C\). Inversion is switch-only; adjoining it
   yields all of \(\mathrm{PGL}(2,q)\not\le\mathrm{Aut}(C)\).
2. **Orbits.** Counts \(9/42/128\) at \(p=3,5,7\); \(\kappa\)-constant; \(T\)-equitable;
   strictly finer than class_key at \(p=5,7\).
3. **Spectrum.** \(\dim E_{4p}^{G}=0,2,7\) at \(p=3,5,7\); \(\lambda_{\max}=4p\) for
   \(p=5,7\). Aut-line \(\dim\le1\) **fails** for \(G\) at \(p\ge5\).
4. **Residual projection.** \(\delta=P_{E_{4p}^{G}}m_4\), \(\delta^2=\|Pm_4\|_2^2\)
   recovers \(\delta^2=1536/65\) at \(p=5\) and \(19180800/1840091\) at \(p=7\).
   \(G\) **carries** the residual at \(p=7\) (class_key does not).

**OPEN:** Gauss / character-sum evaluation of \(m_4\) on \(G\)-orbits (Max+-free)
then project to prove \(\delta^2\le\mathrm{room}/24\) for all primes \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15134.py`, `evidence/e1_gmin_m4_prop15134.json`,
`evidence/e1_gmin_m4_prop15134_strict_aut_spectrum.json`,
`evidence/e1_gmin_m4_prop15134_residual_proj.json`, `tests/test_prop15134.py`.

## Prop 15.135 (2026-08-01) — Coherent-mass spectral form; halfspace char sums

Continues 15.134. Does **not** soft-close residual.

**Proved / certified:**

1. **Spectral form.** On an \(L^2\)-ONB \(\{v_j\}\) of \(E_{4p}^{G}\),
   \(\delta=\sum c_j v_j\) with \(c_j=\langle m_4,v_j\rangle=E_y Q_j(y)\), and
   \(\delta^2=\sum c_j^2\). Residual \(\Leftrightarrow\sum c_j^2\le\mathrm{room}/24\).
2. **Halfspace character formula.** \(f_{\mathrm{hs}}(S)\) is Max+-free via
   \(F_p\)-coordinate indicators; \(\sum_S f_{\mathrm{hs}}=e_4\) (same closed form as \(\sum m_4\)).
3. **\(G\cdot\mathrm{hs}\) dead.** \(|O_{\mathrm{hs}}|\in\{60,168\}\) at \(p=5,7\) vs
   \(|\mathrm{Max}^+|\in\{260,11452\}\); \(\delta^2\) from \(m_4^{G\mathrm{hs}}\gg\mathrm{room}\).
4. **Moments don't pin \(\delta\).** \(P_G(\mathbf{1})=P_G(\kappa)=0\) on \(E_{4p}^{G}\);
   \(e_4\) and \(\langle m_4,\kappa\rangle\) do not constrain the free \(c_j\).
5. **Program.** Need character sums over **full** Max+ for the \(c_j\) (Max+ has
   multiple \(G\)-orbits of vectors; \(Q_j\) not constant).

**OPEN:** Gauss/character-sum for \(c_j\) or \(m_4\) on \(G\)-orbits for all primes
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15135.py`, `evidence/e1_gmin_m4_prop15135.json`,
`tests/test_prop15135.py`.

## Prop 15.136 (2026-08-01) — Max+-free flat on \(G\)-orbits; free \(c_j\) remain

Continues 15.135. Does **not** soft-close residual.

**Proved / certified:**

1. **Max+-free flat.** On the \(G\)-orbit space,
   \(\rho_{\min}=A^+(T\kappa/p^2)\) (resolvent on \((\ker A)^\perp\)) and
   \(\mathrm{flat}=\kappa/p^2+\rho_{\min}\) depend only on \(C\). Matches
   \(\|\rho_{\min}\|^2=\rho_{\min}^2\) and \(\|\mathrm{flat}\|^2=m4f_{\mathrm{flat}}\)
   at \(p=3,5,7\); construction available for all primes \(p\ge3\).
2. **Decomposition.** \(m_4=\mathrm{flat}+\delta\) with \(\delta\in E_{4p}^{G}\);
   residual \(\delta^2=\sum c_j^2\). Certified match at \(p=5,7\).
3. **Geometry insufficient.** The invariant
   \((\mathbf{1}_{\infty\in S},\kappa,\dim_{\mathbb F_p}\mathrm{affspan})\)
   does not determine \(m_4\) (split types at \(p=5,7\)).
4. **Character-sum form.** \(c_j=N^{-1}\sum_{y\in\mathrm{Max}^+}Q_j(y)\) with
   Max+-free kernels \(v_j\). Only these \(\nu_G\) coefficients remain.
5. **Partial.** At \(p=5\), some \(\infty\)-orbits have
   \(m_4\in\{-1/5,-21/65\}\).

**OPEN:** Gauss/character-sum evaluation of each \(c_j\) over full Max+ for all
primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15136.py`, `evidence/e1_gmin_m4_prop15136.json`,
`tests/test_prop15136.py`.

## Prop 15.137 (2026-08-01) — \(c_j\) over Max+ \(G\)-orbits; p=5 two-type formula

Continues 15.136. Does **not** soft-close residual.

**Proved / certified:**

1. **G-equivariance.** \(Q_j(g\cdot y)=Q_j(y)\); \(Q_j\) constant on Max+ \(G\)-orbits.
2. **Hemisphere formula.** \(G\) fixes \(\infty\), so
   \(c_j=\sum_t w_t Q_j(y_t)\) over \(G\)-orbits in \(H_+=\{y_\infty=+1\}\).
3. **Census.** \(r=1,2,5\) hemisphere types at \(p=3,5,7\); p=5 weights
   \(3/13\) (hs) and \(10/13\) (other); p=7 sizes \(\{56,84,294,588,1176\}\).
4. **p=5 formula.** \(c_j=(3/13)Q_j(\mathrm{hs})+(10/13)Q_j(y_*)\);
   \(\sum c_j^2=\mathrm{room}=1536/65\).
5. **\(Q_j(\mathrm{hs})\) Max+-free.** Both \(v_j\) and \(f_{\mathrm{hs}}\) are Max+-free.
   Using only \(Q_j(\mathrm{hs})\) for \(c_j\) fails when \(r>1\) (\(\sum Q_j(\mathrm{hs})^2\gg\mathrm{room}\)).

**OPEN:** Max+-free non-hs representatives \(y_t\) and character sums \(Q_j(y_t)\)
for general primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15137.py`, `evidence/e1_gmin_m4_prop15137.json`,
`tests/test_prop15137.py`.

## Prop 15.138 (2026-08-01) — Max+-free non-hs \(y_*\); p=5 residual Max+-free

Continues 15.137. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **hs-switch.** \(C'=D_{\mathrm{hs}}CD_{\mathrm{hs}}\) has \(C'\mathbf1=p\mathbf1\);
   \(y=\mathrm{hs}\odot z'\) with \(C'z'=pz'\) lies in Max+.
2. **Norm circles.** \(S_{t,c}=\{u:N(u-t)=c\}\) (\(N(u)=u^{p+1}\)); lex-first
   \((t^*,c^*)\) making \(S\) a \(C'\)-Hoffman coclique yields Max+-free
   \(y_*\in H_+\) of non-hs type when found (cert \(p=5,7,11\); empty at \(p=13\)).
3. **p=5 residual Max+-free.** With \((t,c)=(0,3)\),
   \(c_j=(3/13)Q_j(\mathrm{hs})+(10/13)Q_j(y_*)\) (both Max+-free) gives
   \(\sum c_j^2=1536/65=\mathrm{room}\) — residual without Max+ census.
4. **Partial p=7.** Norm circles cover hs (size 84) and one non-hs orbit
   (size 588) only — 2 of 5 \(H_+\) types.

**OPEN:** Remaining \(H_+\) G-orbit types for all primes \(p\ge7\); full residual
for all \(p\ge5\). L remains OPEN (do not promote L from p=5 alone).

Evidence: `src/e1_gmin_m4_prop15138.py`, `evidence/e1_gmin_m4_prop15138.json`,
`tests/test_prop15138.py`.

## Prop 15.139 (2026-08-01) — Affine halfspaces + double switch; p=7 size classes

Continues 15.138. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine halfspaces.** For \(L\not\equiv0\) linear and \(|S|=(p+1)/2\),
   \(y_u=+1\Leftrightarrow L(u)\in S\) is Max+ at certified \(p=5,7\) (all such \(S\)).
2. **AP dichotomy at \(p=7\).** Among all \(\binom{7}{4}=35\) sets \(S\):
   4-term AP \(\Rightarrow\) \(H_+\) orbit size 84 (21 sets); non-4-AP \(\Rightarrow\)
   size 56 (14 sets, including QR-half \(\{0,1,2,4\}\)).
3. **Double Seidel–norm-circle.** \(y=y_0\odot z\) with \(C_0=D_{y_0}CD_{y_0}\)
   and \(z\) a norm-circle Hoffman evec of \(C_0\) yields Max+.
4. **All \(p=7\) size classes Max+-free:**
   84 (AP affine), 56 (non-AP affine), 588 (hs\(\odot\)nc),
   1176 (\(y_{56}\odot\)nc), 294 (\(y_{\mathrm{nc}}\odot\)nc).

**OPEN:** Max+-free weights \(w_t\) and character sums \(Q_j(y_t)\) for every
orbit (including four size-1176 orbits); residual for general \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15139.py`, `evidence/e1_gmin_m4_prop15139.json`,
`tests/test_prop15139.py`.

## Prop 15.140 (2026-08-01) — Weights \(|G|/|\mathrm{Stab}|\); character-sum \(Q_j\); p=7 residual form

Continues 15.139. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Orbit–stabiliser weights.** \(|O_t|=|G|/|\mathrm{Stab}^\infty(y_t)|\);
   stabs \(\{42,28,8,4,2\}\) for sizes \(\{56,84,294,588,1176\}\);
   \(|H_+|=5726\); weights sum to 1.
2. **Character-sum residual.** \(c_j=\sum_t w_t Q_j(y_t)\) with \(Q_j\) on
   G-quotient \(v_j\in E_{4p}^{G}\) recovers
   \(\sum c_j^2=19180800/1840091=\delta^2_{\mathrm{pair}}\le\mathrm{room}=3072/55\)
   at \(p=7\).
3. **Coverage.** Seven of eight \(H_+\) orbits have Max+-free geometric reps
   (affine + double Seidel–norm-circle); **one** size-1176 orbit still lacks
   a Max+-free construction (Q from census only for that type).

**OPEN:** Max+-free \(y\) for the remaining size-1176 orbit; full Max+-free
\(c_j\) at \(p=7\); residual for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15140.py`, `evidence/e1_gmin_m4_prop15140.json`,
`tests/test_prop15140.py`.

## Prop 15.141 (2026-08-01) — Size-12 Seidel partner; p=7 residual Max+-free

Continues 15.140. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Size-12 Seidel partner.** For affine halfspace \(y_0\) with
   \(S=\{2,3,4,5\}\) and field set
   \(T=\{10,12,13,16,18,25,29,36,38,42,44,48\}\),
   \(z=-1\) on \(T\), \(C_0=D_{y_0}CD_{y_0}\) has \(C_0z=pz\), and
   \(y_\sharp=y_0\odot z\) is Max+ of \(H_+\) orbit size 1176 with the
   previously missing \(Q_j\) signature — Max+-free.
2. **All eight \(H_+\) types Max+-free** at \(p=7\).
3. **Residual Max+-free at \(p=7\).** Free weights \(w_t=|G|/|\mathrm{Stab}|\)
   and free \(Q_j(y_t)\) give \(\sum c_j^2=\delta^2_{\mathrm{pair}}\le\mathrm{room}\).
4. **Bi-tight at \(p=5,7\).** mult\(\ge d-1\) + residual \(\Rightarrow\) 16N
   \(\Rightarrow\) bi-tight empty (form already proved); residual Max+-free
   at \(p=5,7\).

**OPEN:** General primes \(p\ge5\) (uniform type law beyond p=7-explicit \(T\));
deep ND; Main Theorem. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15141.py`, `evidence/e1_gmin_m4_prop15141.json`,
`tests/test_prop15141.py`.

## Prop 15.142 (2026-08-01) — Uniform affine law; partners; p=11 sample

Continues 15.141. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine all-\(S\).** Every \(S\subset\mathbb F_p\) with \(|S|=(p+1)/2\) gives a
   Max+ affine halfspace at \(p=5,7,11\) (all \(\binom{p}{k}\) sets).
2. **\(k\)-AP split.** \(p=5\): all AP (one affine orbit). \(p=7\): 21 AP \(\to84\),
   14 non-AP \(\to56\). \(p=11\): non-AP further splits (orbits 132, 330, 660).
3. **Fourths-coset partners.** \(z=-1\) on \(t+a\cdot\mathbb F_q^{\times4}\) works
   at \(p=5\) only; **no** evec hits at \(p=7,11\) — not a uniform size-12 law.
4. **\(p=7\) size-12 fibre.** 84 distinct field sets of size 12; explicit \(T\)
   of 15.141 is one; \(k=(q-1)/4\) but \(T\) is not a fourths coset.
5. **\(p=11\) samples.** Max+-free: affine 132/330/660; ystar 3630;
   double-switch 3630/7260. Full residual OPEN.

**OPEN:** Complete type law for \(p\ge11\); free \(Q_j\) residual for all
\(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15142.py`, `evidence/e1_gmin_m4_prop15142.json`,
`tests/test_prop15142.py`.

## Prop 15.143 (2026-08-01) — p=11 affine 6-orbit census; double-switch LB

Continues 15.142. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Affine type census at \(p=11\).** All \(\binom{11}{6}=462\) affine halfspaces
   form **exactly 6** \(H_+\) G-orbits under strict Aut \(G\) (\(|G|=14520\)):
   sizes \(132\times1\), \(330\times2\), \(660\times3\), with constructive samples.
2. **Double Seidel–norm-circle.** ystar orbit 3630; double-switch sizes include
   3630 and 7260; norm-circle count on \(C_0\) depends on affine class
   (size-132 rich, size-660 empty). Lower bound \(|H_+|\ge28182\).
3. **Incomplete non-affine list.** Full \(H_+\) type census at \(p=11\) OPEN.

**OPEN:** Non-affine types at \(p\ge11\); free \(c_j\) residual; general \(p\ge5\).
L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15143.py`, `evidence/e1_gmin_m4_prop15143.json`,
`tests/test_prop15143.py`.

## Prop 15.144 (2026-08-01) — Free orbits; type-enum residual dead for \(p\ge11\)

Continues 15.143. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Free G-orbit.** Affine \(S=\{0,1,2,4,5,7\}\) double-switch path
   \((t,c)=(33,3)\to(69,9)\) yields Max+ with \(|O|=|G|=14520\) (trivial stab).
2. **Size 2420.** ystar chain \((22,4)\to(91,4)\to(25,4)\to(95,5)\) yields
   orbit size 2420. Non-affine sizes include \(\{1210,2420,3630,7260,14520\}\).
3. **Type-enum residual DEAD for \(p\ge11\).** Deep double-switch multiplies
   free orbits — completing a Max+ G-orbit type list for \(c_j=\sum w_t Q_j\)
   is not a viable Max+-free proof path.
4. **Redirect.** Prefer type-free Max+-free residual: \(\delta^2\le\rho_{\min}^2\)
   for \(p\ge7\), or pointwise \(\sum_j Q_j(y)^2\le\mathrm{room}\).

**OPEN:** Type-free residual for general \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15144.py`, `evidence/e1_gmin_m4_prop15144.json`,
`tests/test_prop15144.py`.

## Prop 15.145 (2026-08-01) — Type-free residual package

Continues 15.144. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction algebra, Max+-free):**

1. **Dictionary.** For primes \(p>\sqrt5\),
   \[
   \delta^2\le\rho_{\min}^2
   \;\Longleftrightarrow\;
   \|\rho\|_2^2\le 2\rho_{\min}^2
   \;\Longleftrightarrow\;
   \|m_4\|_2^2\le m4f_{\mathrm{suf}}
   \;\Longleftrightarrow\;
   \mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}},
   \]
   with closed forms
   \[
   m4f_{\mathrm{suf}}=\frac{(p^2-1)(p^2+1)(3p^4+37p^2+60)}{24p^2(p^2-5)},
   \quad
   \mathrm{ED4}_{\mathrm{suf}}=\frac{4(p^2+1)(3p^6-3p^4+7p^2-15)}{p^2(p^2-5)}.
   \]
   For \(p\ge7\), \(\delta^2\le\rho_{\min}^2\) implies Path C residual
   (\(\rho_{\min}^2<\mathrm{room}_{\mathrm{hyp}}/24\)).
2. **Asymptotic.** \(\rho_{\min}^2/(\mathrm{room}_{\mathrm{hyp}}/24)\to 5/8\) as \(p\to\infty\)
   (monotone on primes \(7\ldots97\)).
3. **Type-free targets.** \(E\|P f_y\|^2\le\rho_{\min}^2\), pointwise
   \(Q_\delta\le\rho_{\min}^2\), \(\|m_4\|^2\le m4f_{\mathrm{suf}}\), or
   \(\mathrm{ED4}\le\mathrm{ED4}_{\mathrm{suf}}\) (weight enumerator / Gauss sums).

**Certified:** \(\delta^2\le\rho_{\min}^2\) at \(p=5,7\) only (prior residual closures).

**OPEN:** type-free \(\delta^2\le\rho_{\min}^2\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15145.py`, `evidence/e1_gmin_m4_prop15145.json`,
`tests/test_prop15145.py`.

## Prop 15.146 (2026-08-01) — Type-free R₄ / μ₄ residual channel

Continues 15.145. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **ED4_from_exact3** \(= -n^4+28n^2-40n\).
2. **R₄ dictionary.** \(\delta^2=(2/3)(R_4-R4_{\mathrm{flat}})\);
   \(\delta^2\le\rho_{\min}^2\Leftrightarrow R_4\le R4_{\mathrm{suf}}=R4_{\mathrm{flat}}+(3/2)\rho_{\min}^2\);
   equivalently \(\mu_4\le\mu4_{\mathrm{suf}}\) and
   \(\bar E[\prod_4 w]\le R4_{\mathrm{suf}}/(n)_4\).
3. **Central moments.** \(\mu_2=n/2\), \(\mu_3=0\) Max+-free.
4. **Spectral mass (too weak).** Under \(\lambda_{\max}(T)=4p\),
   \(\|P f_y\|^2\le (p^2+1)(p^2-2)(p^2-3)(p^2-4)/24\)
   with \(w^*=(p^2-3)(p^2-4)/(p^2(p^2-1))\) — far above residual budgets,
   so Jensen+spectral moments cannot close residual.

**Certified:** R₄ channel at \(p=5,7\).

**OPEN:** \(R_4\le R4_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15146.py`, `evidence/e1_gmin_m4_prop15146.json`,
`tests/test_prop15146.py`.

## Prop 15.147 (2026-08-01) — Inclusion-density residual; ULC near-miss

Continues 15.146. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Falling factorial.** \(R_4=E[k^{\underline 4}]\); under design moments,
   \(\mathrm{exact}_{\le3}=6E[k^3]-11E[k^2]+6E[k]\).
2. **Inclusion densities.**
   \(d_1=1/2\), \(d_2=(p^2+1)/(4p^2)\), \(d_3=(p^2+3)/(8p^2)\);
   \(d_4=R_4/(n)_4\); for \(p\ge7\), residual \(\Leftrightarrow d_4\le d4_{\mathrm{suf}}\).
3. **ULC comparison.** \(U=d_3^2/d_2=(p^2+3)^2/(16p^2(p^2+1)) < d4_{\mathrm{suf}}\)
   for all primes \(p\ge5\) via
   \(P(x)=x^5-8x^4+78x^3-4x^2-7x-60=(x-1)(x^4-7x^3+71x^2+67x+60)>0\) at \(x=p^2\ge25\).
   Thus \(d_4\le U\) would close residual.
4. **ULC fails slightly (census).** At \(p=5,7\): \(d_4/U\approx 1.036,1.019\)
   while \(d_4\le d4_{\mathrm{suf}}\) still holds. Near-miss.

**OPEN:** \(d_4\le d4_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15147.py`, `evidence/e1_gmin_m4_prop15147.json`,
`tests/test_prop15147.py`.

## Prop 15.148 (2026-08-01) — Relaxed-ULC residual calculus

Continues 15.147. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Linear defect.** \(C_{\mathrm{act}}=p^2(d_4/U-1)=C_{\mathrm{flat}}+\kappa\delta^2\) with
   \(C_{\mathrm{flat}}=x P_{\mathrm{flat}}(x)/((x-5)(x-1)(x-2)(x+3)^2)\),
   \(P_{\mathrm{flat}}=x^4-8x^3+58x^2-64x+13\),
   \(\kappa=24p^2/((n-2)(n-3)(n+2)^2)\), \(x=p^2\).
2. **C_max.** \(C_{\mathrm{max}}=Q(p^2)/((p^2-5)(p^2-2)(p^2+3)^2)\),
   \(Q=x^4-7x^3+71x^2+67x+60\);
   residual \(\Leftrightarrow C_{\mathrm{act}}\le C_{\mathrm{max}}\);
   \(C_{\mathrm{max}}<1\to1^-\).
3. **Uniform criterion (p≥7).** \(C_{\mathrm{max}}(p)\ge C_7=C_{\mathrm{max}}(7)=79923/87373\);
   thus \(d_4\le U(1+C_7/p^2)\) implies residual for all primes \(p\ge7\).
4. **Census window.** \(C_{\mathrm{act}}(5)\approx0.900\), \(C_{\mathrm{act}}(7)\approx0.907\);
   constant-\(C\) implication window \([C_{\mathrm{act}}(7),C_7]\) nonempty; \(C=1\) fails.

**OPEN:** prove \(d_4\le U(1+C_7/p^2)\) (or \(C_{\mathrm{act}}\le C_{\mathrm{max}}\)) type-free. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15148.py`, `evidence/e1_gmin_m4_prop15148.json`,
`tests/test_prop15148.py`.

## Prop 15.149 (2026-08-01) — Size-bias residual form

Continues 15.148. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction, Max+-free):**

1. **Size-bias.** Let \(\mu\propto W_k\,k^{\underline 3}\). Then
   \(d_4/d_3=E_\mu[k-3]/(n-3)\) and residual for \(p\ge7\) is
   \(E_\mu[k]\le k_{\mathrm{suf}}=3+8\,R4_{\mathrm{suf}}/(n(n-2)(n+2))\).
   Shift: \(E_\mu[k]-k_{\mathrm{flat}}=12\delta^2/(n(n-2)(n+2))\).
2. **\(k_{\mathrm{flat}}\) shift.** \(k_{\mathrm{flat}}-n/2\to 3^-\) with closed
   \(P_k(p^2)/(2(p^2-5)(p^2+1)(p^2-1)(p^2+3))\).
3. **Independence excesses.** \(d_2-1/4=1/(4(n-1))\), \(d_3-1/8=3/(8(n-1))\),
   \(d4_{\mathrm{flat}}-1/16=P_{\mathrm{ind}}(p^2)/(16(p^2-5)n_4)>0\) for \(p\ge5\);
   \(k_{\mathrm{flat}}>(n+3)/2\) (binomial size-bias).
4. **Uniform \(C_7\).** \(E_\mu[k]\le k_{C7}\Rightarrow\) residual for \(p\ge7\).
5. **Gauss program.** \(\mu\) mixes regular-set sizes through srg triple types
   with Aut-constant \(\lambda_\tau\) — finite-type character-sum target.

**OPEN:** \(E_\mu[k]\le k_{\mathrm{suf}}\) for general \(p\ge7\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15149.py`, `evidence/e1_gmin_m4_prop15149.json`,
`tests/test_prop15149.py`.

## Prop 15.150 (2026-08-01) — srg triples; \(\lambda_e\); \(\pi_e\) residual

Continues 15.149. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Triple counts.** Closed \(n_e\) for edge-type \(e=0,1,2,3\) in the
   conference srg after switch:
   \(n_3=p(p^2+1)(p^2-1)(p-3)/48\), etc.
2. **Covering numbers.** Under Aut-constancy on edge-types and affine
   \(\lambda_e=A+Be\), design moments \(j\le3\) force
   \(\lambda_e=N(p+3-2e)/(8p)\). Certified at \(p=5,7\).
3. **Mixture.** \(\pi_e=n_e(p+3-2e)/\mathrm{Tot}\) is Max+-free (\(N\) cancels);
   \(E_\mu[k]=\sum_e\pi_e m_e\) with \(m_e=\) mean regular-set size through a
   type-\(e\) triple.
4. **Residual.** For \(p\ge7\): \(\sum\pi_e m_e\le k_{\mathrm{suf}}\).

**OPEN:** bound \(m_e\) (character sums / regular-set constraints). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15150.py`, `evidence/e1_gmin_m4_prop15150.json`,
`tests/test_prop15150.py`.

## Prop 15.151 (2026-08-01) — \(m_e\) covariance formula

Continues 15.150. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Regular-set identities.** \(\sum t_e=C(k,3)\), \(\sum e\,t_e=e_S(k-2)\),
   \(t_2+3t_3=k\binom{\alpha}{2}\).
2. **\(E[t_e]\) Max+-free.** \(E[t_e]=n_e(p+3-2e)/(8p)\).
3. **Covariance formula.**
   \(m_e=n/2+8p\,\mathrm{Cov}(k,t_e)/(n_e(p+3-2e))\);
   \(E_\mu[k]=n/2+(8p/\mathrm{Tot})\,\mathrm{Cov}(k,C(k,3))\).
4. **Exact \(m_e\) at \(p=5\).** Fractions certified; \(\sum\pi m\le k_{\mathrm{suf}}\);
   \(m_0,m_1<k_{\mathrm{suf}}<m_2,m_3\).
5. **CS/Popoviciu dead** for residual (π-average exceeds \(k_{\mathrm{suf}}\)).

**Note:** weight-constancy of \(t_e\) (Thm B) overstated for \(p=7\); see 15.152.

**OPEN:** closed \(t_e(k)\) or character-sum \(\mathrm{Cov}(k,t_e)\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15151.py`, `evidence/e1_gmin_m4_prop15151.json`,
`tests/test_prop15151.py`.

## Prop 15.152 (2026-08-01) — free-param \(t_3\); multi-orbit; residual \(\equiv R_4\)

Continues 15.151. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Free parameter.** On regular sets,
   \(t_2=R_3-3t_3\), \(t_1=R_2-2R_3+3t_3\), \(t_0=C(k,3)-R_2+R_3-t_3\).
2. **\(p=5\) closed form.** Mono-type per weight;
   \(t_3(\alpha)=0\) (\(\alpha\le3\)), \(t_3(\alpha)=3\alpha^2-21\alpha+40\) (\(\alpha\ge4\)).
3. **\(p=7\) multi-orbit.** Weights \(k\in\{16,18,\ldots,34\}\) have 3–8 distinct
   \(t\)-vectors (full Max+ census \(N=11452\), \(W=86\)). Pure \(t_e(k)\) is **dead**
   for general \(p\). Corrects 15.151.B.
4. **Residual \(\equiv R_4\).** With \(R_4=E[k^{\underline4}]\),
   \(\mathrm{Cov}(k,C(k,3))=(R_4+(3-n/2)E_3)/6\) and
   \(E_\mu[k]=3+8R_4/(n(n-2)(n+2))\). Multi-type structure does not open a new
   residual channel beyond the weight-4 falling moment.
5. **Per-type reduction.** \(\mathrm{Cov}(k,t_e)=\mathrm{Cov}(\mathrm{det}_e(k))+\gamma_e\mathrm{Cov}(k,t_3)\)
   with \(\gamma=(-1,3,-3,1)\); \(E[t_3]=n_3(p-3)/(8p)\).

**OPEN:** \(R_4\le R4_{\mathrm{suf}}\) / character-sum \(m_e\) on fixed triples /
\(\mathrm{Cov}(k,t_3)\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15152.py`, `evidence/e1_gmin_m4_prop15152.json`,
`tests/test_prop15152.py`.

## Prop 15.153 (2026-08-01) — switched \(\mu_4\) residual dictionary

Continues 15.152. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Switched low moments.** \(e_1=e_3=0\), \(e_2=1/(n-1)\) Max+-free.
2. **Four-point expansion.**
   \(d_4=(1+6/(n-1)+\mu_4)/16\) with \(\mu_4=\mathrm{avg}\,m_4^{\mathrm{sw}}\).
3. **Closed budgets.**
   \(\mu4_{\mathrm{flat}}=(3p^2+17)/(p^2(p^2-2)(p^2-5))\),
   \(\mu4_{\mathrm{suf}}=(3p^4+37p^2+60)/(p^4(p^2-2)(p^2-5))\).
   Residual for \(p\ge7\): \(\mu_4\le\mu4_{\mathrm{suf}}\).
4. **Census.** \(\mu_4\le\mu4_{\mathrm{suf}}\) at \(p=5,7\); exact \(m_e\) at \(p=7\).

**OPEN:** Paley/Weil bound \(\mu_4\le\mu4_{\mathrm{suf}}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15153.py`, `evidence/e1_gmin_m4_prop15153.json`,
`tests/test_prop15153.py`.

## Prop 15.154 (2026-08-01) — \(\mathrm{avg}(\chi\kappa)=3/(n-3)\); \(\eta\) residual

Continues 15.153. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Switch.** \(m_4^{\mathrm{sw}}(S)=\chi(S)\,m_4(S)\), \(\chi=\prod z_i\) halfspace.
2. **Combinatorial average (Max+-free).** By Seidel-switched conference \(C_2\)
   row-sum algebra:
   \(\mathrm{avg}(\chi\kappa)=3/(n-3)=3/(p^2-2)\).
3. **Decomposition.**
   \(\mu_4=\kappa_{\mathrm{main}}+\eta\) with
   \(\kappa_{\mathrm{main}}=3/(p^2(p^2-2))\),
   \(\eta=\mathrm{avg}(\chi\,\mathrm{Ext})/(4p)\).
4. **Residual.** \(\mu_4\le\mu4_{\mathrm{suf}}\Leftrightarrow\eta\le\eta_{\mathrm{suf}}\) with
   \(\eta_{\mathrm{suf}}=4(13p^2+15)/(p^4(p^2-2)(p^2-5))\).
   Strictly \(\kappa_{\mathrm{main}}<\mu4_{\mathrm{flat}}<\mu4_{\mathrm{suf}}\).
5. **Census.** \(\eta\le\eta_{\mathrm{suf}}\) at \(p=5,7\).

**OPEN:** Weil/Aut bound on \(\eta\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15154.py`, `evidence/e1_gmin_m4_prop15154.json`,
`tests/test_prop15154.py`.

## Prop 15.155 (2026-08-01) — Aut-line \(e_4/T\chi/Q\); \(\eta=c_1 R_4+c_0\)

Continues 15.154. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **\(e_4(s)\).** For any \(\pm1\) vector with sum \(s\):
   \(e_4=s^4/24+((-3n+4)/12)s^2+n(n-2)/8\).
2. **\(T\chi\).** \(T\chi(S)=\chi(S)(4p-2\sigma_z(S))\) with
   \(\sigma_z=\sum_{uv\subset S}C_{2,uv}\).
3. **\(Q(s)\).** On \(C_2\)-eigenvectors:
   \(Q=(p/4)[s^2(n-4)+n(6-n)]\).
4. **Average.** \(\mathrm{avg}(\chi\,\mathrm{Ext})=E[(4p\,e_4-2Q)/C(n,4)]\).
5. **Affine.** \(\eta=c_1 R_4+c_0\) with
   \(c_1=16/(n)_4\), \(c_0=-(p^4+4p^2-9)/(p^2(p^2-2))\).
   Residual \(\Leftrightarrow R_4\le R4_{\mathrm{suf}}\) (pure \(E[s^4]\)).
6. **Crude dead.** \(E[s^4]\le 2n^3\) exceeds budget.

**OPEN:** Weil/Paley \(E[\langle z,y\rangle^4]\) or spherical 3-design defect. L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15155.py`, `evidence/e1_gmin_m4_prop15155.json`,
`tests/test_prop15155.py`.

## Prop 15.156 (2026-08-01) — \(\kappa_4=E[s^4]-12n^2\) residual dictionary

Continues 15.155. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Cumulant form.** \(\kappa_4:=E[s^4]-12n^2\); residual \(\Leftrightarrow\kappa_4\le\kappa4_{\mathrm{suf}}\).
2. **Closed budgets.**
   \(\kappa4_{\mathrm{flat}}=16(p^2+1)(p^2+3)/(p^2-5)\),
   \(\kappa4_{\mathrm{suf}}=4(p^2+1)(9p^4+22p^2-15)/(p^2(p^2-5))\).
3. **Bridge.** \(\kappa_4=(n)_4\eta-16n\).
4. **Design orientation.** Spherical 4-design value \(3n^4/(d(d+2))\) is a
   **lower** bound among 2-designs; lies below \(\mathrm{ED4}_{\mathrm{flat}}\) for \(p\ge5\).
5. **Dead ends.** Crude \(E[s^4]\le2n^3\); moment LP on allowed weights (factor \(\gtrsim3.5\) at \(p=5\)).
6. **Census.** \(\kappa_4\le\kappa4_{\mathrm{suf}}\) at \(p=5,7\) (ratios \(\approx0.90,0.66\)).

**OPEN:** Weil/Paley bound on \(\kappa_4\), or upper bound on the spherical
3-design defect of Max+ in \(V_+\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15156.py`, `evidence/e1_gmin_m4_prop15156.json`,
`tests/test_prop15156.py`.

## Prop 15.157 (2026-08-01) — Gegenbauer / 3-design defect residual

Continues 15.156. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Expansion.** \(t^4=a_0+a_2 Q_2+a_4 Q_4\) on \(S^{d-1}\) with
   \(a_0=3/(d(d+2))\), \(a_2=6(d-1)/(d(d+4))\),
   \(a_4=(d^2-1)/((d+2)(d+4))\).
2. **2-design reduction.** \(E[s^4]=n^4(a_0+a_4\mu_{G4})\) with
   \(\mu_{G4}=E[Q_4(s/n)]\ge0\).
3. **Residual.** \(\mu_{G4}\le\mu_{G4,\mathrm{suf}}\) where (with \(x=p^2\))
   \(\mu_{G4,\mathrm{suf}}=4(21x^3+19x^2+35x-75)(x+9)/[x(x-5)(x+1)^3(x+3)(x-1)]\).
4. **Census.** Defect positive but inside budget at \(p=5,7\) (ratios \(\approx0.94,0.83\)).
5. **Dead UBs.** \(\mu\le1\); \(\mu\le1/h_4\) (false at \(p=5\)); \(\mu\le d/h_4\) (too weak for \(p\ge7\)).

**OPEN:** Weil/Aut upper bound on \(\mu_{G4}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15157.py`, `evidence/e1_gmin_m4_prop15157.json`,
`tests/test_prop15157.py`.

## Prop 15.158 (2026-08-01) — closed \(Q_4\); Max+ non-scheme; pole bound

Continues 15.157. Does **not** soft-close residual for all \(p\ge5\).

**Proved:**

1. **Closed \(Q_4\).**
   \(Q_4(t)=[(d+2)(d+4)t^4-6(d+2)t^2+3]/(d^2-1)\), \(Q_4(0)=3/(d^2-1)\).
2. **1-homogeneous tight frame.** Absolute distributions identical; Gram
   spectrum \(2N\) (mult \(d\)) + \(0\) (mult \(N-d\)). Certified \(p=5\).
3. **Not IP-scheme.** Intersection numbers fail constancy on several
   inner-product classes at \(p=5\) — Bose–Mesner on pure \(R_s\) blocked.
4. **Pole decomposition.**
   \(\mu_{G4}\le 2/N+P(E)+Q_4(0)P(\mathrm{Eq})\); pure pole bound
   \(\mu_{G4}\le 2/N\) if \(W\) avoids positive-\(Q_4\) weights except poles.
5. **Conditional residual.** If \(N\ge N_*=\lceil 2/\mu_{G4,\mathrm{suf}}\rceil\) and
   nonpositive-\(Q_4\) support, residual holds. \(N\ge N_*\) at \(p=5,7\), but
   support hypothesis fails (equator / Hoffman exterior).
6. **Chebyshev split dead** (\(\mathrm{UB}\sim 0.18\gg\mu_{G4,\mathrm{suf}}\)).

**OPEN:** Weil/Paley \(\sum_y Q_4(\langle z,y\rangle/n)\), Aut-coherent configuration,
or closed \(W_k\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15158.py`, `evidence/e1_gmin_m4_prop15158.json`,
`tests/test_prop15158.py`.

## Prop 15.159 (2026-08-03) — Φ|Z spectrum structure; dual gap \(G=(d/32)(16I-\Phi)\)

Continues 15.158 / strategy reframe. Does **not** soft-close residual for all \(p\ge5\).

**Proved / certified:**

1. **Exact Φ spectrum at \(p=5\) (certified Fraction).** On \(Z\),
   \(\lambda\in\{176/13\,(\times d),\,144/13\,(\times 2d),\,80/13\,(\times 2d)\}\)
   with \(d=13\), sum mult \(=m=65\). In particular
   \(\lambda_{\max}=176/13=16(d-2)/d\), mult\(=d\).
2. **Exact Φ spectrum at \(p=7\) (certified Fraction).** With den \(=N/(4p)=409\),
   \(\lambda\in\{4320,4032,3648,3360,3072\}/409\) with mults
   \((d,2d,2d,4d,2d)\). \(\lambda_{\max}=4320/409<16(d-2)/d\), mult\(=d\).
3. **Design threshold algebra (proved Fraction).** For \(d>2\),
   \(16(d-2)/d<16\). For \(d\ge13\) (\(p\ge5\)),
   \(\bar\mu=8(d-1)/(d-3)\le16(d-2)/d\) iff \((d-3)(d-6)\ge0\).
4. **Dual gap operator.** \(G:=(d/32)(16I-\Phi)\). At \(p=5\): eigs of \(G\) are
   \(\{1,2,4\}\); at \(p=7\): \(G\succeq I\) (strict). Hence \(G\succeq I\Rightarrow\lambda_{\max}\le16(d-2)/d<16\Rightarrow16\mathrm N\).
5. **16N chain predicate (proved, conditional).** mult\(\ge d\) and \(\|\kappa\|_F^2\le96n\) \(\Rightarrow16\mathrm N\) for \(p\ge5\) (15.105 restated).

**OPEN:** dual gap \(G\succeq I\) (or mult\(\ge d\) + \(\|\kappa\|^2\le96n\)) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15159.py`, `evidence/e1_gmin_m4_prop15159.json`,
`tests/test_prop15159.py`.

## Prop 15.160 (2026-08-03) — Dual-gap vs Hypothesis H; \(H\Rightarrow G\succeq I\)

Continues 15.159. Does **not** soft-close residual for all \(p\ge5\).

**Proved (Fraction):**

1. **H vs thr_ray.** \(H(p)=(p+2)^2/d\), thr_ray\(=5-16/d\).
   \(H-\mathrm{thr\_ray}=(p-5)(p+1)/(2d)\) (equivalently checked Fraction form).
   Hence \(H\le\mathrm{thr\_ray}\) for all primes \(p\ge5\), equality only at \(p=5\).
2. **\(H\Rightarrow\) dual gap.** ray_max\(\le H(p)\) and \(p\ge5\) \(\Rightarrow\) ray_max\(\le\mathrm{thr\_ray}\) \(\Rightarrow G\succeq I\Rightarrow16\mathrm N\).
3. **\(H\Rightarrow16\mathrm N\) for \(p\ge3\).** \(H(p)\le5\) (eq only \(p=3\)), so ray\(\le H\Rightarrow\) ray\(\le5\Rightarrow16\mathrm N\) (15.63).
4. **Census.** \(p=5\): ray \(=H=\mathrm{thr\_ray}=49/13\). \(p=7\): ray \(<H<\mathrm{thr\_ray}\).

**OPEN:** Hypothesis H (ray_max\(\le H(p)\)) for all primes \(p\ge5\); equivalent residual forms \(\delta^2\le\mathrm{room}_{\mathrm{hyp}}/24\), orth\(\le\mathrm{room}_{\mathrm{hyp}}\), \(\|\kappa\|^2\le\kappa_{\mathrm{hyp}}\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15160.py`, `evidence/e1_gmin_m4_prop15160.json`,
`tests/test_prop15160.py`.

## Prop 15.161 (2026-08-05) — Φ-frame of Max+; 16N via mult\(\ge d\) + \(\kappa_4\le48n\)

Continues 15.159–15.160. Does **not** soft-close residual / L.

**Proved (Fraction / Max+-free):**

1. **Constant embedding norm.** \(v_y:=P_Z(yy^\top-I)\) satisfies \(\|v_y\|_F^2=n(n-2)\) for all \(y\in\mathrm{Max}_+\) (from \(\mathrm{tr}(\Phi)=m\bar\mu=n(n-2)\) and 1-homogeneity).
2. **Pairwise frame Gram.** With \(v_y=yy^\top-(n/d)P_+\in Z\),
   \(\langle v_y,v_z\rangle_F=(y\cdot z)^2-2n\).
3. **16N budgets under mult\(\ge d\).** \(E[s^4]\le12n(n+4)\) \(\Leftrightarrow\kappa_4\le48n\) \(\Rightarrow\lambda_{\max}\le16\) when mult\(\ge d\); two-level bulk \(b=8\).
4. **Census.** mult\(=d\) and \(\kappa_4\le48n\) at \(p=5,7\).

**OPEN:** mult\((\lambda_{\max})\ge d\) and \(\kappa_4\le48n\) for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15161.py`, `evidence/e1_gmin_m4_prop15161.json`,
`tests/test_prop15161.py`.

## Prop 15.162 (2026-08-05) — Maximizers in \(Z\); mult\(\ge d-1\); \(E[s^4]\) type expansion

Continues 15.161. Does **not** soft-close residual / L.

**Proved:**

1. **Maximizers of \(\Gamma\) lie in \(Z\).** Criticality on Sym\(_0\): \(E[f\,cc^\top]=\lambda A\) (\(\mu=0\) from \(E[f]=0\)); ambient diagonal \(B_{ii}=E[f]/\lambda=0\). Hence mult\((\Gamma\text{ top})=\)mult\((\Phi\text{ top})\).
2. **mult\((\Phi)\ge d-1\) (proved for all primes \(p\ge5\)).** Thm A + Prop 15.97 (mult \(\Gamma=\) mult \(\lambda_2(P\odot P)\)) + Prop 15.98 (PSL min nontrivial irrep dim \(d-1\)).
3. **Type expansion.** \(E[s^4]=C_0+R\) with
   \(C_0=n(3n-2)+2n(n-1)(3n-4)/p^2\) and \(R=24\sum_{4\text{-sets}}m_4^2\ge0\).
4. **16N \(\Leftrightarrow\) mult\(\ge d\) + \(m_4\)-mass.** \(\sum m_4^2\le n(3p^2+61)/24\) \(\Leftrightarrow\kappa_4\le48n\).
5. **Census.** Identities + 16N at \(p=5,7\).

**OPEN:** upgrade mult\(\ge d-1\to d\); prove \(m_4\)-mass / \(\kappa_4\le48n\) for all \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15162.py`, `evidence/e1_gmin_m4_prop15162.json`,
`tests/test_prop15162.py`.

## Prop 15.163 (2026-08-05) — Wick \(m_4\) mass; Aut\(_0\); \(H_C\) split

Continues 15.161–15.162. Does **not** soft-close residual / L.

**Proved / certified:**

1. **Wick 4-set mass (proved Fraction).** \(T=C_{ab}C_{cd}+C_{ac}C_{bd}+C_{ad}C_{bc}\), \(m_4^W=T/p^2\);
   \(\sum T^2=n(n-1)(n-2)(n-5)/8\), \(\sum(m_4^W)^2=(p^4-1)(p^2-4)/(8p^2)\).
2. **η-room after Wick (proved Fraction).** Under \(\sum m_4^W\eta=0\),
   16N \(\Leftrightarrow\sum\eta^2\le n(19p^2-3)/(6p^2)\). Census usage \(\approx95\%,73\%\) at \(p=5,7\).
3. **Aut\(_0\) on \(V_+\) (structure + cert \(p=5\)).** \(V_+\cong\mathbf1\oplus\sigma\) with \(\dim\sigma=d-1\).
4. **\(H_C\) split of Φ-top (cert \(p=5,7\)).** Top mult splits as \(1+(d-1)\).

**OPEN:** mult\(\ge d\) general \(p\); \(\sum\eta^2\le\eta_*\) / 16N for all \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15163.py`, `evidence/e1_gmin_m4_prop15163.json`,
`tests/test_prop15163.py`.

## Prop 15.164 (2026-08-05) — 16N from mult\(\ge d-1\) + \(E[s^4]\le\mathrm{Es4}_*(p)\)

Continues 15.161–15.163. Preferred 16N path (no mult\(\ge d\) needed). Does **not** soft-close L.

**Proved (Fraction):**

1. **Two-level majorization.** Fixed sum \(T\), sum of squares \(Q\), mult\(\ge k\), \(\lambda_i\ge\ell_{\min}\): max top \(L\) is two-level; \(L\) increasing in \(Q\).
2. **Es4\(_*\) budget under mult\(\ge d-1\).** With \(\lambda_{\min}\ge6\), bulk
   \(b_*=8(p^2-1)(p^2-7)/(p^4-8p^2-1)\ge8>6\) for primes \(p\ge5\).
   If \(E[s^4]\le\mathrm{Es4}_*(p)\) then \(\lambda_{\max}\le16\) (16N).
3. **Equivalent forms.** \(E[s^4]\le\mathrm{Es4}_*\Leftrightarrow\kappa_4\le\kappa4_*\Leftrightarrow R\le R_*\Leftrightarrow\sum\eta^2\le\eta_*\) (Wick orth).
4. **Census.** 16N at \(p=5,7\) via spectrum Es4.

**OPEN:** \(E[s^4]\le\mathrm{Es4}_*(p)\) (or \(\sum\eta^2\le\eta_*\)) for all primes \(p\ge5\).
This is the single analytic residual for Path-C 16N given mult\(\ge d-1\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15164.py`, `evidence/e1_gmin_m4_prop15164.json`,
`tests/test_prop15164.py`.

## Prop 15.165 (2026-08-05) — Exact Es4; closed Es4\(_*\)/\(\eta_*\); GoG\(\leftrightarrow\Phi\); \(m_4\) is C-eigen

Continues 15.164. Does **not** soft-close L.

**Proved / certified:**

1. **Moments.** \(E[s]=0\) (central symmetry), \(E[s^2]=2n\) (2-design / \(\sum yy^\top=2N P_+\)).
2. **GoG \(\leftrightarrow\) Φ spectrum.** \(\mathrm{spec}(G\circ G)=\{2nN\}\cup\{N\lambda:\lambda\in\mathrm{spec}(\Phi|_Z)\}\cup\{0\}^{N-1-m}\);
   hence \(E[s^4]=4n^2+\mathrm{tr}(\Phi^2)\).
3. **Closed budgets (Fraction).**
   \[
   \mathrm{Es4}_*(p)=\frac{4(3p^8-6p^6-148p^4-10p^2+129)}{p^4-8p^2-1},\quad
   \eta_*(p)=\frac{(p^2-1)(p^2+1)(19p^4-152p^2-3)}{6p^2(p^4-8p^2-1)}.
   \]
4. **\(m_4\) is C-eigen.** \(p\cdot m_4(a,b,c,d)=\sum_j C_{aj}m_4(j,b,c,d)\) (cert random 4-sets \(p=5,7\)).
5. **Exact Es4 census.** \(p=3,5,7\) via Gram/Φ spectrum (not single-root \(W\); \(p=7\) not 1-homogeneous). H-saturation at \(p=5\).

**OPEN:** Es4\(_*\) / 16N for all primes \(p\ge5\). L remains OPEN.

Evidence: `src/e1_gmin_m4_prop15165.py`, `evidence/e1_gmin_m4_prop15165.json`,
`tests/test_prop15165.py`.

## Prop 15.166 (2026-08-05) — 16N \(\Leftrightarrow\lambda_{\max}(Q_2)\le4N/(d(d-1))\); Wick C-eigen

Continues 15.164–15.165. Does **not** soft-close L.

**Proved:**

1. **Max+ is a spherical 2-design in \(V_+\).** \(\sum\hat u=0\), \((1/N)\sum\hat u\hat u^\top=I_d/d\).
2. **Wick \(m_4\) is C-eigen** (same equation as true \(m_4\)); residual \(\eta=m_4-m_4^W\) is invisible to the C-eigen constraint — \(P_+(\mathrm{Wick})\) does not recover true \(m_4\).
3. **16N \(\Leftrightarrow Q_2\) bound (proved Fraction).**
   \(\lambda_{\max}(\Phi)=4d(d-1)/N\cdot\lambda_{\max}(Q_2)\), hence
   \(\lambda_{\max}(\Phi)\le16\Leftrightarrow\lambda_{\max}(Q_2)\le4N/(d(d-1))\).
4. **Census.** \(Q_2\) thr holds at \(p=5\) (ratio \(\approx0.846\)), \(p=7\) (\(\approx0.660\)).

**OPEN:** \(\lambda_{\max}(Q_2)\le4N/(d(d-1))\) for all primes \(p\ge5\). Dead: Delsarte LP, BM\((C)\), equating Wick to \(m_4\). Preferred: Weil/Jacobi Aut-orbit \(m_4\), Aut\(_0\) isotype, SOS \(Q_4\). L remains OPEN. residual_closed_general=false.

Evidence: `src/e1_gmin_m4_prop15166.py`, `evidence/e1_gmin_m4_prop15166.json`,
`tests/test_prop15166.py`.

## Prop 15.167 (2026-08-05; corrected 2026-08-30) — spectral majorization algebra; bi-tight conclusion retracted

**Valid Fraction algebra (conditional on the unproved floor):**

1. **Majorization UB.** mult\((\lambda_{\max})\ge d-1\) (15.162) + \(\lambda_{\min}(\Phi)\ge6\) + \(\mathrm{tr}(\Phi)=n(n-2)\) \(\Rightarrow\)
   \[
   \lambda_{\max}(\Phi)\le L_*(p)=\frac{p^4+24p^2-1}{2(p^2-1)}.
   \]
2. **\(L_*<2d\).** \(2d-L_*=(p^4-24p^2-1)/(2(p^2-1))\); numerator \(=24\) at \(p=5\) and \(f(x)=x^2-24x-1\) increasing on \(x=p^2\ge25\).
3. **Retracted final arrow.** The implication from a simple top eigenvalue to tight-cover emptiness used the false identity \(\ker(G-(n/2)P_{\mathbf1})=\mathrm{span}\{\mathbf1\}\). The correct kernel also contains \(\ker G\), including all star differences from 15.56. Thus the majorization gap does not prove bi-tight emptiness.
4. **Census.** The spectral inequalities at \(p=5,7\) remain valid, but they have no discrete-cover conclusion by themselves.

The required bi-tight levels are instead closed by Proposition 15.720.

Evidence: `src/e1_gmin_m4_prop15167.py`, `evidence/e1_gmin_m4_prop15167.json`,
`tests/test_prop15167.py`, `src/e1_bitight_chain.py`.

## Prop 15.168 (2026-08-05; corrected 2026-08-30) — E(1) structure with the 15.720 bi-tight gate

Continues 15.167. Does **not** soft-close E(1) or L.

**Proved / checkable predicates (Fraction + prior props):**

1. **Required bi-tight alternatives.** Levels 2 and 3 are empty for every \(p\ge5\) by 15.720. Generic Max\(_+\)-tight covers are not claimed empty.
2. **Deep tight empty** for \(p\ge5\) (15.44.3 + level-2 case of 15.720).
3. **Type I freeness ND** (prior 15.43.1).
4. **Type I freeness-fail \(k=2p-1\)** \(\to\) tight size \(2p\); the master lemma gives no-descent or a level-2 bi-tight alternative, and 15.720 excludes the latter.
5. **Deep auto-freeness** for \(s_+=2\), \(k\le3p-2\): \(N_2/N\) lb \(=2-k/(2p)>(p+1)/(2p)\).
6. **Deep fail-eq \(k=3p-1\)** \(\Rightarrow\) tight \(S\equiv3\) size \(3p\); the level-3 bi-tight alternative is excluded by 15.720.

**OPEN residuals (honest — no soft-close):**

- Type I freeness-fail at \(k=3p-2\) / \(S\in\{1,5\}\) boundary (not reduced to tight \(2p\)).
- Deep non-tight freeness-fail with \(k\ge3p\) (freeze-to-tight sketch not shipped as predicate).

Full \(m_n\ge\Phi-2\) / E(1) / L remain OPEN. residual_closed_general=false. E1_closed_general=false.
Within this Paley route, \(L=1/2\) closes only if bi-tight \(\land\) E(1)
(denseness Prop 6.2) — currently false.  A different convergence proof need
not use this gate.

Evidence: `src/e1_gmin_m4_prop15168.py`, `evidence/e1_gmin_m4_prop15168.json`,
`tests/test_prop15168.py`, `src/e1_main_chain_status.py`.

## Prop 15.169 (2026-08-05) — Type I \(k=3p-2\) ND reduction; deep multi-\(s\) auto-freeness

Continues 15.168. Does **not** soft-close E(1) or L.

**Proved (Fraction / prior props):**

1. **Type I freeness-fail structure at \(k=3p-2\).** At freeness equality \(N_1/N=(p+1)/(2p)\) with \(S\in\{1,5\}\): \(a=(5-k/p)/4\) equals the freeness threshold; affine \(S+2f_e=3\) on Max\(_+\); \(H=G\cup\{e\}\) has size \(3p-1\), scores \(\{2,4\}\), \(s_+^H=2\).

2. **2-Lipschitz of \(\Phi\) under edge flip.** \(Q_y(A\oplus e)=Q_y(A)-2f_e(y)\) with \(f_e=\pm1\), so \(|\Phi(A\oplus e)-\Phi(A)|\le2\). Corollary: if \(\Phi(C\oplus G)\ge\Phi\) then \(\Phi(C\oplus G\oplus e)\ge\Phi-2\) (weak ND for every \(e\)).

3. **Gap-2 undercutter forces \(s_-=-1\).** Type I \(\Rightarrow\Phi\ge\Phi-2\) (15.42.1). Gap-2 means \(\Phi=\Phi-2\), requiring \(s_-\le-1\). With \(k=3p-2\) odd and \(\mathbb E_-[S]=-3+2/p>-3\), odd scores force \(s_-\ge-1\); hence \(s_-=-1\).

4. **ND dichotomy at \(s_-=-1\).** On \(U_-=\{S_G=-1\}\): \(f_e=+1\) gives \(|Q_H|=\Phi\) (strong ND); \(f_e\equiv-1\) on \(U_-\) (bad case) gives only \(\Phi(H)\ge\Phi-4\) from Max\(\pm\) dichotomy.

5. **Deep multi-\(s\) auto-freeness.** For min-level \(s\ge2\) with scores \(\ge s\) step 2: \(N_s/N\) lb \(=(s+2-k/p)/2\); auto-freeness for \(k\le p(s+1)-2\). Recovers \(s=2\Rightarrow k\le3p-2\) (15.168.F). Boundaries checked for \(s=2,3,4,5\).

**OPEN (honest — no soft-close):**

- Residual (i): prove freeness-fail Type I \(k=3p-2\) cannot have \(s_-\le-1\) (or bad case impossible) for all primes \(p\ge5\). **Certified** \(s_-\le-1\) integrally infeasible at \(p=5\) (MILP/HiGHS), not a general proof.
- Residual (ii): deep freeness-fail ND for \(s_+=2\), \(k\ge3p\).

E1_closed_general=false. residual_closed_general=false. L OPEN.

Evidence: `src/e1_gmin_m4_prop15169.py`, `evidence/e1_gmin_m4_prop15169.json`,
`tests/test_prop15169.py`.

### 15.169 addendum — bad-case dual two-level identities

If freeness-fail Type I \(k=3p-2\) occurs with \(s_-=-1\), Max\- two-level \(\{-1,-3\}\), and \(f_e\equiv-1\) on \(U_-\) (bad case), then \(H=G\cup\{e\}\) is dual two-level: mass \(\mathrm{thr}\) at \(S_H=\pm2\), mass \(1-\mathrm{thr}\) at \(\pm4\), with
\[
\mathbb E[S_H^2]=10-\frac6p,\qquad |H|=3p-1=2p+(p-1).
\]
Also \(\sum_{e'<e''\in G}G^+_{e'e''}=\frac{\mathbb E_+[S^2]-k}2=-\frac{3(p-1)(p-4)}{2p}<0\) for \(p\ge5\).

**Status (2026-08-06):** residual (i) **not** closed for general \(p\). Dual-equality Farkas is **conditional** on a disj Gsum LB that is **not** proved in general (see Prop 15.170). Fractional affine+\(s_-\le-1\) remains feasible (integrality gap).

## Prop 15.170 (2026-08-05; status revised 2026-08-06) — Residual (i): Type I \(k=3p-2\), \(s_-\le-1\)

Continues 15.169. **OPEN for general primes \(p\ge5\).** Structure and conditional Farkas algebra shipped; general close blocked by disj Gsum LB hinge.

**Proved / checkable (Fraction / prior props):**

1. **Freeness-fail + gap-2 structure** (15.169): affine \(S=3-2f_e\) on Max\(_+\); gap-2 forces \(s_-=-1\); ND dichotomy at \(s_-=-1\) (good sign \(\Rightarrow\) strong ND; bad case is the only residual risk).

2. **Dual equality correlation (conditional).** If the bad case realises the dual two-level law \(S=-3-2f_e\) on Max\(_-\) with affine Max\(_+\), then
   \[
   (\mathrm{Gsum}\,x)_e = \frac6p-4,
   \]
   where \(\mathrm{Gsum}=\mathbb E_+[ff^\top]+\mathbb E_-[ff^\top]\) and \(x=1_G\).

3. **Gsum identities.** \(\mathrm{Gsum}_{ee}=2\); wedge (share a vertex) \(\mathrm{Gsum}_{ab}=0\) (proved: \(\mathbb E[y_iy_j]_+\!+\!\mathbb E[y_iy_j]_-=0\)); \(\mathrm{Gsum}\,\mathbf1=n\mathbf1\).

4. **Disj Gsum LB — NOT proved for general \(p\).** Candidate bound
   \[
   \mathrm{Gsum}_{ab}\ge -\frac{12}{p(p^2+1)}\quad(a\neq b\text{ disjoint})
   \]
   was previously justified as an association-scheme minimum and is **certified at \(p=5\)** (\(-6/65=-12/(5\cdot26)\)). **Prop 15.158:** Max+ is **not** an IP association scheme, so scheme-min cannot prove the LB for all primes \(p\ge5\). Shipped flag: `gsum_disj_lb_proved_general()=False`.

5. **Conditional box-sum Farkas.** *If* the candidate LB holds, then for \(0\le x\le1\), \(\mathbf1^\top x=k=3p-2\), \(x_e=0\),
   \[
   (\mathrm{Gsum}\,x)_e \ge -\frac{12k}{p(p^2+1)},
   \]
   and dual equality conflicts for all primes \(p\ge5\) because
   \[
   \frac6p-4 < -\frac{12k}{p(p^2+1)} \iff 4p^3-6p^2-32p+18>0.
   \]
   That poly algebra is real (`dual_equality_farkas_algebra_if_lb`); it does **not** alone close residual (i).

6. **Consequence (honest).** residual (i) **OPEN** until disj Gsum LB is proved Max+-free (or an alternate residual-(i) proof). `type_I_k_3p_minus_2_closed_general()=False`. Also \(\mathrm{ES}_2=(13p-12)/p<k\) (integrality seed, independent).

**Addendum (corrected 2026-08-30).** The Gsum path above remains unused
(`gsum_disj_lb_proved_general=False`). The **historical two-level Type-I
residual-(i) slice** is closed by the dual-eq path: 15.272
\(k=1\cup k=3\) span \(\Rightarrow G_+\succ0\Rightarrow\ker=\mathrm{sc}\)
(15.207) \(\Rightarrow\mathrm{cost}_D<2-\alpha\) (15.249)
\(\Rightarrow\) dual-eq empty (15.216). This does not close the current
multi-level Type-I gate. See Prop 15.272 and
`evidence/share/denseness_path_package.md`.

Evidence: `src/e1_gmin_m4_prop15170.py`, `evidence/e1_gmin_m4_prop15170.json`,
`tests/test_prop15170.py`, `evidence/share/denseness_path_package.md`.

## Prop 15.171 (2026-08-05; status revised 2026-08-06) — Residual (ii): deep freeness-fail \(k\ge3p\) ND

Continues 15.168–170. **OPEN for general primes \(p\ge5\).** Structure (parity, auto-freeness, fail-eq empty) shipped; dual two-level Farkas blocked by the same disj Gsum LB hinge as 15.170.

**Proved / checkable (Fraction / prior props):**

1. **Parity.** \(s_+=2\) (even scores) \(\Rightarrow k\) even \(\Rightarrow\) Max\- scores even \(\Rightarrow s_-\ne-1\).

2. **Gap-2 deep undercutter classification.** \(s_+\!=2\) gap-2 undercutter \(\Rightarrow s_-\le-2\) (\(s_-\ge0\) not undercutting; \(s_-=-1\) impossible by parity).

3. **Deep freeness \(\Rightarrow\) weak ND.** If \(f_e\not\equiv+1\) on \(\{S=2\}\), some \(y\) has \(S=2\), \(f_e=-1\), \(S_H=1\), \(Q=\Phi-2S_H=\Phi-2\).

4. **Auto-freeness** \(k\le3p-2\) (15.168) \(\Rightarrow\) freeness \(\Rightarrow\) weak ND.

5. **Fail-eq \(k=3p-1\):** freeness-fail \(+\,S\in\{2,4\}\) \(\Rightarrow H\) tight \(S\equiv3\) size \(3p\) empty under bi-tight/Thm A (15.167–168) when bi-tight holds.

6. **Dual two-level freeness-fail Farkas (conditional).** Two-level \(S\in\{2,4\}/\{-2,-4\}\) with freeness-fail affine \(f_e=3-S\) forces
   \[
   (\mathrm{Gsum}\,x)_e=2\bigl(8-\tfrac{3k}p\bigr).
   \]
   Box-sum needs the same unproved disj LB \(\mathrm{Gsum}_{ab}\ge-12/(pn)\). Under that LB, need \(<\) LB on freeness-fail range (Fraction). Without proved LB: **not closed**.

7. **Consequence (honest).** Residual (ii) **OPEN**. With residual (i) open and Gsum hinge open: **E1_closed_general=false**, **L OPEN** (no soft-close). Denseness cannot force \(L=\tfrac12\) without E(1).

**Addendum (corrected 2026-08-30).** Propositions 15.179, 15.236, and
15.237 close only the bounded even range through \(4p-2\). They do not close
global residual (ii): the non-Walsh multi-level range at even \(k\ge4p\)
remains open, so `residual_ii_full_closed=False`. Proposition 15.272 later
closes the historical two-level Type-I/residual-(i) slice but did not close
multi-level Type I. Proposition 15.750 now closes multi-level Type I for every
prime \(p\ge5\), leaving residual (ii) as the sole false entry in the
then-defined four-unit ledger. Proposition 15.764 later proves that this
ledger is not by itself exhaustive for arbitrary minimal four-gap paths. The
Gsum route remains unused; Path-C / \(16N\) remain optional open.

Evidence: `src/e1_gmin_m4_prop15171.py`, `evidence/e1_gmin_m4_prop15171.json`,
`tests/test_prop15171.py`, `evidence/share/denseness_path_package.md`.

## Prop 15.272 (2026-08-15) — \(k=1\cup k=3\) spans \(\mathcal W_{++}^0\)

**Statement.** For every prime \(p\ge5\), the \(k=1\) cylinders together with the \(k=3\) locked triples span \(\mathcal W_{++}^0\). Hence \(G_+\succ0\), \(\ker(\mathrm{Gsum})=\mathrm{sc}\), dual-eq is empty, and the two-level Type I residual-(i) slice is closed. This does **not** close the multi-level Type-I gate, residual (ii) for even \(k\ge4p\), the spectral floor, E(1), or \(L\).

Aut-Schur is **false** (Jacquet \(\nRightarrow\) PSL-span of \(k=3\) \(F\); \(p=5\) rank \(61/65\)). Gsum disj LB is unused. The pairing \(1^\top K^{-1}v\) is unused.

**Proof (linear).** See `evidence/share/denseness_path_package.md` Lemmas B–J, written for a cold read. In outline:

1. \(V_+\cong\mathbb C^{\{0\}\cup\Omega}\) (FT of \(Cv=pv\); \(\hat z(0)=p v_\infty\); support \(\subseteq\{0\}\cup\Omega\)).
2. Each Aut_∞ isotypic of \(F^\perp\) is the convolution hyperplane on unique pairs of \(S_\mu\) (zero-diag is DFT inversion, not a dimension leap).
3. Bad \(\mu\): \(k=3\) triangles + connected line graph of \(K_m\) span \(1^\perp\). Every triple occurs.
4. Good \(\mu\): Johnson products \(P_S(k)=\hat f(k)\hat f(\alpha-k)\) span the same-line hyperplane of dim \(m-1\) via \(B_{xy}=\omega^{\alpha y}\hat c(x-y)\), 1-swap, and \(\sigma=2m-p=1\). Complementary mixed + through-\(L_0\) Fejer fill the rest. \(p=5\): Veronese rank \(65\).
5. Singer PD of \(k=3\) on \(F\) for \(p\ge7\).
6. 15.207 + 15.249: \(\mathrm{cost}_D<2-\alpha\) is the rational
   \[
   (2-\alpha)-\alpha(S_1+tN)=\frac{2(p^4-3p^3-10p^2+9p+1)}{p^4-4p^2+1}>0
   \]
   for every prime \(p\ge5\) (\(f(5)=46\), \(f'>0\) on \([5,\infty)\)).

Evidence: `src/e1_gmin_m4_prop15272.py`, `src/e1_gmin_m4_prop15270.py`,
`src/e1_gmin_m4_prop15249.py`, `tests/test_prop15272.py`.

## Prop 15.628 (2026-08-24) — Eligible GQR circles close W1, W2, and Walsh

Let \(C\) be the Paley conference matrix on
\(P^1(\mathbb F_{p^2})\), let \(H_0\) be the binary direction code of
Max\(_-\), fix \(e=\{\infty,0\}\), and let
\(U=\{y\in\mathrm{Max}_-:C_{\infty0}y_\infty y_0=-1\}\).

**Theorem A (eligible-circle span).** The nonsquare Miquelian circles that
meet \(e\) in zero or two points span

\[
H_0\cap\ker(e_\infty+e_0)
\]

for every odd prime \(p\).  Tangent-pencil relations identify all circles
meeting \(e\) oddly modulo the eligible span.  The relevant bipartite
tangency graph is the Cayley graph on
\(M=(\mathbb F_{p^2}^{\times})^2\) with connection set

\[
T=\{(1+r\sigma)^{-2}:r\in\mathbb F_p\}.
\]

If \(T\) did not generate \(M\), a nontrivial even character of
\(\mathbb F_{p^2}^{\times}\) would equal one on an affine
\(\mathbb F_p\)-line.  Parameterising that line by the norm-one torus turns
its sum into Katz's \(t=-2\) Soto--Andrade sum.  The two exceptional
character pairs induce the trivial extension-field character; otherwise
the bound \(2\sqrt p<p\) applies.  Thus the graph is connected for
\(p\ge5\); \(p=3\) is direct.

**Theorem B (all affine halfspaces).** Let
\(L:\mathbb F_{p^2}\to\mathbb F_p\) have square kernel \(\mathbb F_p\).
For any \(T\subset\mathbb F_p\) of size \((p+1)/2\), define

\[
h_T(\infty)=1,
\qquad
h_T(u)=s_T(Lu),
\qquad
s_T(b)=2\mathbf1_T(b)-1.
\]

Then \(Ch_T=ph_T\).  Indeed, \(\sum_b s_T(b)=1\).  At a finite point
with \(L(u)=b\), the zero fibre contributes \((p-1)s_T(b)\), while every
nonzero affine square fibre has character sum \(-1\); hence

\[
(Ch_T)_u=1+(p-1)s_T(b)-\sum_{t\ne0}s_T(b+t)=p s_T(b).
\]

For nonsquare \(\sigma\), put
\(z_T(\infty)=-1\) and \(z_T(\sigma u)=h_T(u)\).  Then
\(Cz_T=-pz_T\).  If \(0\in T\), its restriction to
\(S=\{\infty\}\cup\sigma\mathbb F_p\) is the sparse \(-p\) signing of
that nonsquare circle, so flipping \(z_T\) on \(S\) gives another
Max\(_-\) point with binary difference \(\mathbf1_S\).

**Theorem C (every eligible circle is realised in \(U\)).** A pair inside
\(S\) has the \(U\) sign automatically.  For a pair outside \(S\), write
\(x=\sigma u\), \(y=\sigma v\), \(b=L(u)\ne0\), and \(d=L(v)\ne0\).
The \(U\) condition is

\[
s_T(b)s_T(d)=\chi(u-v).
\]

If \(b=d\), the right side is \(+1\).  If \(b\ne d\), prescribe equal
membership of \(b,d\) when the right side is \(+1\), and opposite
membership when it is \(-1\); for \(p\ge5\) this prescription extends to a
set \(T\) of size \((p+1)/2\) containing zero.  PSL transport now realises
every eligible nonsquare circle as a difference of two points in \(U\).

Combining A and C gives

\[
\operatorname{dir}(U)=H_0\cap\ker(e_\infty+e_0)
\]

for every \(p\ge5\).  At \(p=3\), the exact direction rank is \(4\), equal
to the target, and W2 is vacuous.  Therefore W1, W2, and Walsh 15.406 E
hold for every odd prime.  This closes the Walsh slice; it does not by
itself close the unrelated 5+-level or other E1 leftovers.

Evidence: `evidence/NOTE_2026-08-24_w2_gqr_circle_route.md`,
`evidence/w2_gqr_circle_route_2026-08-24.json`,
`scripts/w2_affine_circle_close.py`.

## Prop 15.632 (2026-08-25) — affine slack budget and quadratic parity lifts

This proposition concerns the non-Walsh multi-level residual, independently
of the Paley-lattice Props. 15.629--15.631.  It eliminates the
Eulerian-boundary branch for every odd prime but does not close the residual.

Let `H` be an odd set of signed Paley edges, `h=|H|`, and

\[
 S_H(y)=\sum_{e\in H}C_e y_e.
\]

For each projective \(\mathbb F_p\)-direction \(d\), let \(\epsilon_d\) be
the quadratic type of its kernel.  Its affine halfspaces form a copy of the
middle slice \(J(p,m)\), \(m=(p+1)/2\), inside the
\(\epsilon_dp\)-eigenshell.  If `H` separates the affine shells with margin
three, then

\[
 A_d(y)={\epsilon_dS_H(y)-3\over2}
\]

is a nonnegative integer-valued quadratic on this slice.  Put
\(a_d=2p\mathbb E_dA_d\).

For a fixed direction, the contributions to \(p\epsilon_d\mathbb E S_H\)
of an infinity edge, a parallel finite edge, and a transverse finite edge
are respectively

\[
 1,\qquad p,\qquad-\epsilon_dC_e.                          \tag{15.632.1}
\]

Every contribution is odd modulo two, so every \(a_d\) is even.  There are
\(m\) directions of each quadratic type.  An infinity edge contributes
\(m\) in either type.  A finite edge of sign \(c\) contributes
\(p-(m-1)=m\) in type \(c\), and \(m\) transverse `+1` terms in type
\(-c\).  Therefore

\[
 \sum_{d:\epsilon_d=\tau}a_d
   =m(h-3p),\qquad \tau=\pm1.                             \tag{15.632.2}
\]

Let \(D\) be the odd-degree boundary of `H` and
\(c_H=\prod_{e\in H}C_e\).  The edge-product identity gives

\[
 (-1)^{A_d}
 =\epsilon_d(-1)^{(h-3)/2}c_H\prod_{v\in D}y_v.           \tag{15.632.3}
\]

Let \(B_d\) be the affine fibres in direction \(d\) containing an odd
number of finite points of \(D\).  For middle-slice membership bits \(x_s\),

\[
 \prod_{v\in D}y_v
 =\epsilon_d^{\mathbf1_{\infty\in D}}
  (-1)^{|B_d|+\sum_{s\in B_d}x_s}.
\]

Hence an explicit \(\eta_d\in\{0,1\}\) satisfies

\[
 A_d(x)\equiv\sum_{s\in B_d}x_s+\eta_d\pmod2.            \tag{15.632.4}
\]

Average under `Sym(B_d) x Sym(B_d^c)`.  The result is a quadratic
\(q(t)\), \(t=|X\cap B_d|\), and its nonnegative integer inputs all have
parity \(t+\eta_d\).  Thus

\[
 q(t)\ge(t+\eta_d\bmod2).
\]

Let \(M(p,b,\eta)\) be the minimum expectation of such a quadratic under
the exact hypergeometric law

\[
 \Pr(t)={\binom bt\binom{p-b}{m-t}\over\binom pm}.
\]

This is a three-variable LP.  Its vertices are obtained by interpolating
three active parity values (with the one- and two-point supports immediate),
so `M` is computed exactly over rational numbers.  Since \(a_d\) is even,

\[
 a_d\ge2\lceil pM(p,|B_d|,\eta_d)\rceil.                 \tag{15.632.5}
\]

Combining (15.632.2) and (15.632.5) yields a separate necessary budget in
each quadratic direction type.

At residual size \(h=4p+1\), suppose \(D=\varnothing\).  Then every
\(B_d\) is empty and one quadratic-type half has constant odd slack parity.
Every one of its \(m\) directions costs \(2p\), for total \(p(p+1)\), while
(15.632.2) permits only \((p+1)^2/2\).  The contradiction gap is

\[
 {p^2-1\over2}>0.
\]

Thus the Eulerian-boundary branch is empty for every odd prime.

This is not the full close.  The corrected affine model at \(p=5\) has an
integral twenty-edge `G` for which `H=G union {(0,1)}` has directional
values

\[
 (a_d)=(12,4,0,6,10,4),
\]

all affine slacks are nonnegative pointwise, and the boundary is infinity
plus one affine line.  Nonempty boundary profiles and the full non-affine
shell remain open.

Evidence: `src/e1_gmin_m4_prop15632.py`,
`evidence/e1_gmin_m4_prop15632.json`,
`evidence/NOTE_2026-08-25_affine_slack_parity_budget.md`, and
`tests/test_prop15632.py`.

## Prop 15.633 (2026-08-25) — the complete second Paley-dual shell

Retain

\[
 L=\ker_{\mathbb Z}(C-pI),\qquad
 P={1\over2}(I+C/p),\qquad L^*=P\mathbb Z^n,
\]

with \(n=p^2+1\).  For every odd prime \(p\ge5\), the complete shell of
dual norm \((p-1)/p\) is the disjoint union

\[
 \boxed{
 \{u\in L^*:\|u\|^2=(p-1)/p\}
 =\{\pm P(e_i-C_{ij}e_j):i<j\}
  \mathbin{\dot\cup}
  \{\pm w_S/p:S\in\mathcal B_+\},}               \tag{15.633.1}
\]

where \(\mathcal B_+\) is the PSL orbit of square
\(\mathbb F_p\)-sublines and \(w_S\) is their signed complement satisfying
\(Cw_S=pw_S\).  Hence the signed shell size is

\[
 \boxed{N_2(p)=p(p+1)(p^2+1).}                    \tag{15.633.2}
\]

At \(p=3\), the two descriptions overlap and the signed count is \(30\).

Here is the exhaustion argument.  Use the square-direction circle-word
profiles \(a_{j,s}=\langle x,v_{j,s}\rangle\) and their common sum
\(t=\sum_s a_{j,s}=2px_\infty\).  The tight-frame identity gives

\[
 \sum_{j,s}a_{j,s}^2=p\|x\|^2+{t^2\over2}.       \tag{15.633.3}
\]

For \(\|x\|^2=(p-1)/p\), balancing an integral length-\(p\) profile,
together with the evenness of \(t\), leaves only

\[
 |t|\in\{0,2,p-1\}.                               \tag{15.633.4}
\]

Equality at \(t=p-1\) gives one zero in every profile; degree-one glue
makes those zeros the evaluations of one linear form and yields
\(P(e_\infty-e_u)\).  At \(t=0\), the MDS/Newton equality cases have either
one active direction, whose profile is a translated quadratic character,
or all but one direction active.  The first case gives the complements of
square circles through infinity; degree-one and degree-two glue in the
second case gives \(P(e_u-e_v)\).

At \(t=2\), every profile is \(\delta_\alpha+\delta_\beta\).  Its sum and
product glue to binary forms, so

\[
 (\alpha_j-\beta_j)^2=D(t_j)
\]

for one binary quadratic \(D\), nonzero and square on the selected
quadratic-character half of \(\mathbb P^1(\mathbb F_p)\).  A split
nondegenerate \(D\) has too few positive values.  A rank-one \(D\) gives
\(P(e_u+e_v)\).  If \(D\) is anisotropic and not proportional to the norm
form \(N\), then \(Y^2=DN\) is a smooth genus-one curve with either zero or
\(2(p+1)\) rational points, contradicting Hasse for \(p\ge5\).  Thus
\(D=cN\), which gives exactly the square-circle complements not containing
infinity.  Counting these disjoint equality cases proves
(15.633.1)--(15.633.2).

For an admissible harmonic tensor \(W\), the point-pair part contributes

\[
 {1\over4}\left(1-{(p-1)^2\over d+2}\right)\|W\|_F^2,
                                                               \tag{15.633.5}
\]

while the square-circle part contributes

\[
 {1\over8p^4}\sum_{S\in\mathcal B_+}(w_S^TWw_S)^2
 -{(p-1)^2\over4p(d+2)}\|W\|_F^2.               \tag{15.633.6}
\]

This is an exact decomposition of the complete second shell; its operator
spectrum is determined next.

Evidence: `src/e1_gmin_m4_prop15633.py`,
`evidence/e1_gmin_m4_prop15633.json`,
`evidence/NOTE_2026-08-25_dual_second_shell.md`, and
`tests/test_prop15633.py`.

## Prop 15.634 (2026-08-25) — square-circle operator and second-shell sign

Let \(M\) be point--circle incidence for \(\mathcal B_+\), and let \(A\)
join two square circles when they meet in two points.  The design identity
and the three intersection valencies are

\[
 M^TM={p^2-1\over2}I+{p+1\over2}J,
\]

\[
 k_0={p(p-1)(p-3)\over4},\qquad k_1=p^2-1,\qquad
 k_2={p(p^2-1)\over4}.                              \tag{15.634.1}
\]

A quadratic-character count of common two-secants gives, for distinct
circles with intersection size \(j=0,1,2\), respectively

\[
 {(p-1)^2(p+1)\over8},\qquad {p(p^2-1)\over8},\qquad
 {p^3+p^2-9p-1\over8}.
\]

Equivalently,

\[
 \boxed{A^2+pA={p^2-1\over8}MM^T
              +{(p-1)^2(p+1)\over8}J.}             \tag{15.634.2}
\]

Together with

\[
 AM={p^2-1\over4}J+{(p-1)^2\over4}M,
\]

this yields

\[
 \operatorname{Spec}(A)=
 \left\{k_2^1,
 \left({(p-1)^2\over4}\right)^{n-1},
 (-p)^{n(p-1)/4},0^{n(p-3)/4}\right\}.             \tag{15.634.3}
\]

Let \(b_S\) be the projection of \(w_Sw_S^T\) to the admissible tensor
space \(Z=\{W:PWP=W,\operatorname{diag}W=0\}\).  The exact correlations
\(|w_S^Tw_T|=2p,p,0\) according as \(|S\cap T|=0,1,2\), followed by
orthogonal projection off the diagonal tensors, turns (15.634.2) into

\[
 \boxed{\operatorname{Spec}((\langle b_S,b_T\rangle))
 =\{0^n,[p^3(p-1)]^{n(p-1)/4},
          [p^3(p+1)]^{n(p-3)/4}\}.}                 \tag{15.634.4}
\]

Substituting (15.634.4) into (15.633.5)--(15.633.6) diagonalizes the
complete signed norm-\((p-1)/p\) harmonic shadow shell.  Its eigenvalues
are

\[
\begin{aligned}
 \lambda_0&=-{(p+2)(p^2-4p+1)\over4p(p^2+5)},\\
 \lambda_-&=-{p^3-3p^2-19p+9\over8p(p^2+5)},\\
 \lambda_+&=-{p^3-5p^2-19p-1\over8p(p^2+5)},
\end{aligned}                                      \tag{15.634.5}
\]

with multiplicities
\(n(p-1)(p-3)/8,n(p-1)/4,n(p-3)/4\).  All three are strictly negative
for every odd prime \(p\ge11\).  Thus the second shell is an exact
cancellation channel against Prop. 15.631's positive first shell.

Indeed, retaining only the first two shells is positive semidefinite only
when

\[
 \exp\left(-{\pi(p-2)\over8pt}\right)
 \le {p\over(p+2)(p^2-4p+1)},                       \tag{15.634.6}
\]

which forces \(t=O(1/\log p)\).  This kills a first-shell-only proof but is
not a bound on the omitted shells.  R1, global QVAR, and the limit therefore
remain open.

Evidence: `src/e1_gmin_m4_prop15634.py`,
`evidence/e1_gmin_m4_prop15634.json`,
`evidence/NOTE_2026-08-25_square_circle_operator.md`, and
`tests/test_prop15634.py`.

## Prop 15.635 (2026-08-25) — third dual norm and exact p=11 shell

Scale dual norms by \(s=2p\|x\|^2\).  For a common integral profile sum
\(t\), balancing gives

\[
 g_p(t)=(p+1)f_p(t)-t^2
       =pa^2+2ab+b(p+1-b),\qquad |t|=ap+b.          \tag{15.635.1}
\]

The first odd equality values \(|t|=1,p\) give \(s=p\).  For \(t=1\),
subtract the minimum vector determined by the degree-one profile moment.
If the resulting zero-sum profiles are active in \(h\) directions and have
total positive mass \(M\), write \(s=p+4\Delta\).  Integer energy and the
MDS/Newton floor give

\[
 \Delta\ge h,\qquad \Delta\ge M-h,
 \qquad M\ge h(R-h),\quad R={p+1\over2}.
\]

Thus

\[
 \Delta\ge\max\{h,h(R-h-1)\}\ge R-2.
\]

The \(t=p\) equality case is at least as strong after subtracting
\(Pe_\infty\).  All other odd values in (15.635.1) start at
\(|t|=3,p-2\).  Hence every odd-phase vector outside the minimum shell has

\[
 \boxed{s\ge3p-6.}                                  \tag{15.635.2}
\]

For even \(t\), profile-energy parity makes every increase after the second
norm \(2(p-1)\) a multiple of four, and (15.635.1) has no smaller new
value.  Since \(3p-6>2(p+1)\) for \(p\ge11\), the third dual norm is

\[
 \boxed{s_3=2(p+1),\qquad\|x\|^2=(p+1)/p.}          \tag{15.635.3}
\]

The distinct signed vectors

\[
 \pm P(e_i+C_{ij}e_j),\qquad i<j,                  \tag{15.635.4}
\]

attain (15.635.3) and number \(p^2(p^2+1)\).  Their complete degree-four
harmonic contribution is

\[
 \boxed{-{p^2+4p-3\over4(p^2+5)}\|W\|_F^2.}        \tag{15.635.5}
\]

At \(p=11\), exact `qfminim` enumeration of the saturated dual through
\(s=24\) gives \(31{,}110\) signed vectors.  The first two proved shells
contain \(244+16{,}104\), leaving \(14{,}762=11^2(11^2+1)\).  Therefore
(15.635.4) is the complete third shell at \(p=11\).  For larger primes the
third norm and pair orbit are uniform theorems, but completeness of the
whole shell is not asserted.  This result does not control the remaining
theta tail and does not prove R1.

Evidence: `src/e1_gmin_m4_prop15635.py`,
`evidence/e1_gmin_m4_prop15635.json`,
`evidence/NOTE_2026-08-25_third_dual_norm.md`, and
`tests/test_prop15635.py`.

## Prop 15.636 (2026-08-25) — complete third dual shell

Proposition 15.635 leaves only one possible extra equality profile. Put
(m=(p-1)/2). Its positive root multiset has size (m), with one repeated
root, and its disjoint negative root set has size (m); their power sums
agree through degree (m-1). If (A,B) are the monic root polynomials,
Newton identities give (A-B=\mathrm{constant}). Their roots cover all but
two elements (u,v\in\mathbb F_p), and the repeated root is (alpha), so

\[
 AB={(X^p-X)(X-\alpha)\over(X-u)(X-v)}.             \tag{15.636.1}
\]

Normalize (u=0,v=1,alpha=\lambda). Then (T=(A+B)/2) would satisfy

\[
 T(x)^2=x^{p-1}+a(x+\cdots+x^{p-2})+b,qquad
 a=1-\lambda,\quad a(a-1)\ne0.                    \tag{15.636.2}
\]

Reverse (T), and put (c=a-1). Uniqueness of the formal square root at
zero says that the coefficients in degrees (m+1,ldots,p-2) of

\[
 K(y)=(1+cy)^{m+1}(1-y)^m                         \tag{15.636.3}
\]

vanish. Hence (K=U+q y^{p-1}+r y^p), where (deg U\le m) and
(q\ne0). For (1\le j\le m), the Hasse derivative of order (m+j)
is therefore the monomial (q(-1)^{m+j}y^{m-j}). Evaluating it at the two
roots (1,-1/c) of (15.636.3) and comparing local coefficients gives

\[
 c^{m-j}=-{1\over2j}.                              \tag{15.636.4}
\]

The choices (j=m-1,m-2) force (c=1/3) and (c^2=1/5), hence
(1/9=1/5) and (p\mid4), impossible for (p\ge11). Thus the extra
profile does not exist. The complete third shell for every odd prime
(p\ge11) is exactly (15.635.4), with signed count (p^2(p^2+1)) and
harmonic scalar (15.635.5). The fourth and later shells remain uncontrolled,
so R1 remains open.

Evidence: `src/e1_gmin_m4_prop15636.py`,
`evidence/e1_gmin_m4_prop15636.json`,
`evidence/NOTE_2026-08-25_complete_third_dual_shell.md`,
`evidence/r1_third_shell_profile_obstruction_11_17.json`, and
`tests/test_prop15636.py`.

## Prop 15.637 (2026-08-25) — zero-common-sum gap at the next energy

Put (k=(p-1)/2) and (R=k+1). At the first possible post-third even
energy (E=p+3=2k+4), the profile MDS bound

\[
M\ge h(R-h),\qquad E\ge2M                       \tag{15.637.1}
\]

leaves only (h=1,R-1,R) for (p\ge11). Suppose (h=1). Then
(M\ge k), and integer energy leaves two cases.

If (M=k+1), there is one doubled profile entry and the two degree-(k+1)
root polynomials (A,B) cover all of (\mathbb F_p), with one repeated
root (alpha). Their power sums agree through degree (k-1), so (A-B)
is linear and

\[
AB=(X^p-X)(X-\alpha).                            \tag{15.637.2}
\]

Reversing the square identity for (T=(A+B)/2) forces
(T=(X-\alpha)^{k+1}). Hence
((A-B)^2/4=(X-\alpha)^2), making (A) and (B) share (alpha), a
contradiction.

If (M=k), there are two doubled entries and three omitted roots. Normalize
the latter to \(0,1,\rho\), call the repeated roots \(\alpha,\beta\), and put

\[
N=(1-\alpha y)(1-\beta y),\qquad
D=(1-y)(1-\rho y).
\]

The reversed degree-(k) square root agrees through degree (2k-1) with
(S=\sqrt{N/D}=\sum u_jy^j), so
(u_{k+1}=\cdots=u_{2k-1}=0). The differential equation

\[
2ND S'=(N'D-ND')S                              \tag{15.637.3}
\]

has leading multiplier (alpha\beta\rho\ne0). Its coefficient of degree
(j+3) reduces, whenever the next four coefficients vanish, to
(2\alpha\beta\rho j u_j=0). Since (k\ge5), descending from (j=k)
gives (u_1=\cdots=u_k=0), hence (N=D). This identifies the repeated
roots with omitted roots, the second contradiction.

Therefore no one-active profile occurs at energy \(p+3\). It remains to
exclude the two dense active counts.

Let \(q_d\) denote the degree-\(d\) binary form whose values on the square
directions are the \(d\)-th profile moments. Every signed pair
\(\delta_a-\delta_b\) satisfies

\[
 4q_1q_3-3q_2^2-q_1^4=0.                         \tag{15.637.4}
\]

For \(h=R\), the energy budget gives one energy-four profile and
\(R-1=k\ge5\) ordinary signed pairs. The left side of (15.637.4) is a
binary quartic, so its \(k\) ordinary zeros force it to vanish identically.
On the exceptional profile
\(\delta_a+\delta_b-\delta_c-\delta_d\), however, its value is

\[
 -12(a-c)(a-d)(b-c)(b-d)\ne0,                    \tag{15.637.5}
\]

a contradiction.

Now let \(h=R-1=k\). The unique inactive direction is the root of
\(q_1=L\), and \(q_2=LS\) for a linear form \(S\). At every ordinary
pair,

\[
 4q_3=L(3S^2+L^2).                               \tag{15.637.6}
\]

The excess energy is either carried by two energy-four profiles or by one
energy-six profile. In the first case, the cubic difference in
(15.637.6) vanishes at the \(k-2\) ordinary directions and at the inactive
direction. These \(k-1\ge4\) zeros force an identity, contradicting
(15.637.5) at an exception.

In the second case, the same cubic has \(k\ge5\) zeros and is again an
identity. The energy-six pattern
\(2\delta_a-\delta_b-\delta_c\) is excluded because its defect in
(15.637.4) is

\[
 -12(a-b)^2(a-c)^2\ne0.                          \tag{15.637.7}
\]

The sole remaining pattern has three distinct positive and three distinct
negative roots. Since \(q_4\) also vanishes at the inactive direction,
write \(q_4=LT\), with \(T\) cubic. Ordinary pairs satisfy

\[
 2T=S(S^2+L^2).                                  \tag{15.637.8}
\]

This cubic relation has \(k-1\ge4\) ordinary zeros and is therefore an
identity. At the exceptional direction, (15.637.6) and (15.637.8) say
that its first four power-sum differences equal those of a signed pair
\(\delta_u-\delta_v\). Moving \(u,v\) to opposite sides gives two
four-element multisets with equal first four power sums. Newton identities
make the multisets equal, impossible because the original six roots are
distinct across the two signs.

Thus, for every odd prime \(p\ge11\),

\[
 \boxed{\text{no zero-common-sum profile has energy }p+3.}  \tag{15.637.9}
\]

Balancing leaves only the nonzero common sums
\(|t|\in\{2,p-1,p+1\}\) at scaled norm \(2(p+3)\). These three cases and
the later theta tail remain open, so this is not a fourth-shell
classification or an R1 proof.

Evidence: `src/e1_gmin_m4_prop15637.py`,
`evidence/e1_gmin_m4_prop15637.json`,
`evidence/NOTE_2026-08-25_one_profile_next_energy_gap.md`,
`evidence/r1_two_double_square_shift_11_127.jsonl`, and
`tests/test_prop15637.py`.

## Prop 15.638 (2026-08-25) — empty first post-third even candidate shell

Continue with \(R=(p+1)/2\), the common profile sum \(t\), and scaled norm
\(s=2p\lVert x\rVert^2\). At \(s=2(p+3)\), balancing and parity leave

\[
 |t|\in\{0,2,p-1,p+1\}.                           \tag{15.638.1}
\]

Proposition 15.637 excludes \(t=0\). Up to sign, the remaining values
\(2,p-1,p+1\) have total profile-energy excess \(4,4,2\), respectively.
The profile glue supplies binary moment forms \(q_d\) through degree four
for every \(p\ge11\).

For \(t=2\), ordinary profiles are unsigned pairs and satisfy

\[
\begin{aligned}
2q_3-3q_1q_2+q_1^3&=0,\\
2q_4-q_2^2-2q_1^2q_2+q_1^4&=0.                  \tag{15.638.2}
\end{aligned}
\]

With one energy-six exception, the \(k=(p-1)/2\ge5\) ordinary directions
force both identities. The doubled pattern has nonzero cubic defect
\(6(a-c)^2(b-c)\). For the four-positive/two-negative pattern, (15.638.2)
makes its first four moments those of a two-root multiset; moving the two
negative roots and applying Newton identities to two four-element
multisets contradicts disjointness.

With two energy-four exceptions, the \(k-1\ge4\) ordinary directions force
the cubic identity. Its defect on a three-positive/one-negative profile is

\[
6(a-d)(b-d)(c-d)\ne0,
\]

so both exceptions are doubled points. The quadratic
\(D=2q_2-q_1^2\) is therefore a nonzero square on every ordinary point of
the selected projective half \(T\), and has its two distinct roots in
\(T\). If the half is defined by
\(\eta(N)=\varepsilon\) for an anisotropic norm form \(N\), then

\[
\left|\sum_{\mathbb P^1(\mathbb F_p)}\eta(ND)\right|=p-3. \tag{15.638.3}
\]

Indeed, the projective character sum of split \(D\) is zero, while
\(\sum_T\eta(D)=(p-3)/2\). But \(ND\) is squarefree and
\(Y^2=ND\) is a smooth genus-one curve, so Hasse bounds the left side of
(15.638.3) by \(2\sqrt p\). This contradicts
\(p-3>2\sqrt p\) for \(p\ge11\).

For \(t=p-1\), replace each profile by \(b=1-a\); for \(t=p+1\), use
\(b=a-1\). The transformed profiles have common sum one, and ordinary
profiles are deltas. Their moments obey \(q_d=q_1^d\). The only
energy-three exception has quadratic defect

\[
-2(a-c)(b-c)\ne0,
\]

and the doubled energy-five exception has defect
\(-2(a-b)^2\ne0\). The remaining three-positive/two-negative
energy-five pattern would have the first three moments of one delta.
Moving that delta to the negative side and applying Newton identities to
two three-element multisets again contradicts disjointness. The ordinary
direction counts exceed the degrees of all identities used.

Consequently, for every odd prime \(p\ge11\),

\[
 \boxed{\{x\in L^*:2p\lVert x\rVert^2=2(p+3)\}=\varnothing.}
                                                               \tag{15.638.4}
\]

This does not determine the next nonempty dual norm or bound the full
fourth-and-later harmonic theta tail. R1 and global QVAR remain open.

Evidence: src/e1_gmin_m4_prop15638.py,
evidence/e1_gmin_m4_prop15638.json,
evidence/NOTE_2026-08-25_empty_post_third_even_candidate_shell.md,
scripts/r1_next_shell_half_conic_audit.py,
evidence/r1_next_shell_half_conic_11_43.json, and
tests/test_prop15638.py.

## Prop 15.639 (2026-08-25) — complete first nonminimal odd dual shell

For every odd prime \(p\ge11\), the first possible nonminimal odd scaled
norm from Proposition 15.635 is

\[
 s=2p\lVert x\rVert^2=3p-6.                       \tag{15.639.1}
\]

Write \(r=2px=(pI+C)z\). Odd \(s\) makes every coordinate of \(r\) odd,
while

\[
 \sum_i r_i^2=2ps=6p^2-12p<9(p^2+1).
\]

Thus some \(|r_i|=1\). Signed Paley transport moves it to the profile base,
where the common sum is \(t=1\). Subtract the minimum vector encoded by the
degree-one moment and let \(h\) be the number of active zero-sum profiles.
Equality in the Proposition 15.635 bounds leaves only

\[
 h=1\quad\hbox{or}\quad h=(p+1)/2-2.              \tag{15.639.2}
\]

The one-profile case is the square-circle equality family; the sole other
multiplicity pattern is the exceptional third-shell profile excluded by
Proposition 15.636. In the dense case, active profiles are
\(\delta_\alpha+\delta_\beta-\delta_\gamma\). With
\(A=\alpha-\mu\), \(B=\beta-\mu\), their moment differences are

\[
 Q_2=-2AB,\qquad Q_3=-3AB(2\mu+A+B).              \tag{15.639.3}
\]

The two inactive directions are exactly the roots of \(Q_2\), so
\(Q_2\mid Q_3\). The resulting linear quotient encodes a point \(v\).
Adding \(Pe_v\) produces a common-sum-two vector whose scaled norm is

\[
 2(p-1)+4r,\qquad r\in\{0,1,2\}.                 \tag{15.639.4}
\]

The cases \(r=1,2\) are excluded by the complete third shell and the empty
shell of Proposition 15.638. Reducing \(r=0\) through the complete second
shell proves

\[
\boxed{
\{x\in L^*:2p\lVert x\rVert^2=3p-6\}
=\mathcal T\mathbin{\dot\cup}\mathcal O,}
\]

where \(\mathcal T\) consists of projected unit signed triples whose three
signed conference edges are negative, and

\[
 \mathcal O=\{\mathord\pm(Pe_i+w_S/p):(w_S)_i=-1\}.
\]

The scaled-coordinate signatures have respectively three and one entries
of magnitude \(p-2\), proving disjointness and injectivity. Since
\(\operatorname{tr}(C^3)=0\), exactly half of all coordinate triangles are
negative. Counting circles and their off-circle points gives

\[
\boxed{
N_{\rm odd}(p)=\binom{p^2+1}{3}+p^2(p-1)(p^2+1)
=\frac{p^2(p-1)(p+7)(p^2+1)}6.}                 \tag{15.639.5}
\]

Exact NUKA `qfminim` at \(p=11\) through scaled bound \(28\) returned
\(473,970\) signed vectors and maximum norm \(27\). Removing the
\(31,110\) vectors through the third shell leaves \(442,860\), exactly
(15.639.5). The computation audits but is not used by the uniform proof.

For \(p=11,13\), this is the fourth nonempty shell. For \(p\ge17\), further
even candidates can lie below (15.639.1), so no global ordinal-shell claim
is made. Proposition 15.640 supplies this shell's harmonic operator, but the
intervening and later theta tail remains open; R1 and global QVAR are not
proved.

Evidence: src/e1_gmin_m4_prop15639.py,
evidence/e1_gmin_m4_prop15639.json,
evidence/NOTE_2026-08-25_first_nonminimal_odd_scaled_shell.md,
evidence/r1_dual_shell_count_p11_28.json, and tests/test_prop15639.py.

## Prop 15.640 (2026-08-25) — exact harmonic saddle at scaled norm 3p-6

Let \(\mathcal X\) be the complete shell from Proposition 15.639 and let
\(W\) be admissible, with \(F=\lVert W\rVert_F^2\). For every point \(i\),
the square-circle complement words through it obey

\[
\sum_{S\ni i}w_Sw_S^T
=p^2\left(P-2(Pe_i)(Pe_i)^T\right).              \tag{15.640.1}
\]

After signed Paley transport to infinity, the circles split into
\((p+1)/2\) orthogonal parallel classes. Each class has Gram matrix
\(p(pI-J)\), hence frame eigenvalue \(p^2\), and the class dimensions sum
to \((p^2-1)/2=d-1\). This proves (15.640.1).

For the negative signed triples, put \(B=C\circ W\). The row sums of \(B\)
vanish because \(CW=pW\) and \(\operatorname{diag}W=0\). Expanding over
triangles, the unweighted square sum is \((n-4)F/2\), while weighting by the
conference triangle sign gives \(\operatorname{tr}(CW^2)=pF\). Including
both signs yields

\[
\sum_{x\in\mathcal T}(x^TWx)^2=2(p-3)(p+1)F.     \tag{15.640.2}
\]

For an oriented circle word \(w\), write \(q_S=w^TWw\). Represent the
point--circle vectors by
\(\pm(Pe_i-w_iw/p)\), \(i\notin S\). Using (15.640.1) and
\(\sum_Sw_Sw_S^T=p^2(p-1)P\) gives

\[
\sum_{S,i\notin S}(Ww_S)_i^2=p^2(p-2)F
\]

and therefore

\[
\sum_{x\in\mathcal O}(x^TWx)^2
=8(p-2)F+{2(p-5)\over p^3}\sum_S(w_S^TWw_S)^2.  \tag{15.640.3}
\]

The complete shell has squared norm \(q=3(p-2)/(2p)\), signed size
\(N=p^2(p-1)(p+7)(p^2+1)/6\), and second moment \((Nq/d)P\). The radial
terms of \(H_W\) subtract

\[
{2Nq^2\over d(d+2)}F
={3(p-1)(p+7)(p-2)^2\over p^2+5}F.               \tag{15.640.4}
\]

Combining (15.640.2)--(15.640.4) with Proposition 15.634's square-circle
tensor eigenvalues \(0,p^3(p-1),p^3(p+1)\) gives

\[
\boxed{
\begin{aligned}
\lambda_0&=-{p^4+2p^3-69p^2+136p+26\over p^2+5},\\
\lambda_-&={p^4-14p^3+89p^2-196p+24\over p^2+5},\\
\lambda_+&={p^4-10p^3+69p^2-176p-76\over p^2+5}.
\end{aligned}}                                   \tag{15.640.5}
\]

Writing \(x=p-11\), the numerators of
\(-\lambda_0,\lambda_-,\lambda_+\) become respectively

\[
\begin{aligned}
x^4+46x^3+723x^2+4668x+10476,\\
x^4+30x^3+353x^2+2004x+4644,\\
x^4+34x^3+465x^2+3036x+7668.
\end{aligned}
\]

Thus the signs are \((- ,+,+)\) for every \(p\ge11\): the shell is a
quartic saddle, not a spherical 4-design. At \(p=11\), the spectrum is

\[
(-582/7)^{1220},\qquad(258/7)^{305},\qquad(426/7)^{244}.     \tag{15.640.6}
\]

The norm-parity phase and \(H_W(x/2)=H_W(x)/16\) multiply these values by
\(-1/16\), so the transformed shadow signs are \((+,-,-)\). This exact
channel information does not control the intervening even candidates for
\(p\ge17\) or the later theta tail. R1 and global QVAR remain open.

Evidence: src/e1_gmin_m4_prop15640.py,
evidence/e1_gmin_m4_prop15640.json,
evidence/NOTE_2026-08-25_scaled_norm_3p_minus_6_harmonic_saddle.md, and
tests/test_prop15640.py.

## Prop 15.641 (2026-08-25) — current linear modular data do not determine R1

At \(p=11\), the relevant Kohnen subspace of
\(M_{69/2}(\Gamma_0(44),\chi_{44})\) has dimension 66.  Impose the
currently justified linear information: the infinity coefficients before
the complete second dual shell, the forced gaps at the half, zero, quarter,
and \(1/11\) cusps, and the complete second-shell coefficient.  Exact
rational row reduction gives

\[
 \operatorname{rank}A_{\rm gaps}=29,
 \qquad
 \operatorname{rank}\binom{A_{\rm gaps}}{c_{20}}=30.       \tag{15.641.1}
\]

Thus these data leave a 36-dimensional kernel.  The second-shell row
\(c_{20}\) and the first odd-coset half-cusp target row \(c_*\) have joint
rank two modulo the gap rows.  More explicitly, exact elimination produces
a 66-coordinate rational vector \(w\), with only 21 nonzero coordinates,
such that

\[
 A_{\rm gaps}w=0,
 \qquad c_{20}w=0,
 \qquad c_*w=1.                                  \tag{15.641.2}
\]

Consequently the known shell coefficients and geometric cusp gaps do not
determine, or even sign-determine, the R1 target.  This is a negative result
about the linear coefficient-determination route, not a counterexample to
R1: positivity of complete shell operators and further shell coupling are
not imposed in (15.641.1).

Evidence: src/e1_gmin_m4_prop15641.py,
evidence/e1_gmin_m4_prop15641.json,
evidence/NOTE_2026-08-25_p11_modular_independence.md, and the exact witness
identified there by SHA-256.

## Prop 15.665 (2026-08-27) — conserved positive quartic mass on every R1 dual shell

Let \(X_s\) be a complete dual shell of signed size \(N_s\) and common
squared radius \(r_s\).  Put

\[
 \mathcal Z=\{W=W^T:PWP=W,\ \operatorname{diag}W=0\},
 \qquad z=\dim\mathcal Z={n(n-6)\over8},
\]

and, for \(x\in\operatorname{range}P\), let
\(b_x=\Pi_{\mathcal Z}(xx^T)\).  The unphased quartic shell operator

\[
 R_s=\sum_{x\in X_s}b_x\otimes b_x                    \tag{15.665.1}
\]

is positive semidefinite.  Since the shell is a tight frame in the
irreducible \(d=n/2\) dimensional conference eigenspace,

\[
 \sum_{x\in X_s}xx^T={N_sr_s\over d}P.
\]

Substitution in the degree-four harmonic polynomial of Proposition 15.631
therefore gives the exact operator identity

\[
 \boxed{A_s=R_s-\rho_s I_{\mathcal Z}},
 \qquad
 \boxed{\rho_s={2N_sr_s^2\over d(d+2)}}.             \tag{15.665.2}
\]

It remains to compute the trace of \(R_s\) without constructing a basis of
the \(z\)-dimensional space.  The Gram matrix of the diagonal map on
\(\operatorname{Sym}^2(\operatorname{range}P)\) is

\[
 K=P\circ P={(p^2-1)I+J\over4p^2},
 \qquad
 K^{-1}={4p^2\over p^2-1}I-{2\over p^2-1}J.          \tag{15.665.3}
\]

Writing \(r=\lVert x\rVert^2\), orthogonal projection onto the kernel of
the diagonal map yields

\[
 \boxed{
 \lVert b_x\rVert_F^2
 =r^2-{4p^2\over p^2-1}\sum_i x_i^4
       +{2\over p^2-1}r^2.}                         \tag{15.665.4}
\]

If \(\{W_\alpha\}\) is an orthonormal basis of \(\mathcal Z\),
equivariance and trace give
\(\sum_\alpha W_\alpha^2=(z/d)P\).  Hence the scalar harmonic polynomial

\[
 H_{\rm tr}(x)=\lVert b_x\rVert_F^2
 -{4z\over d(d+4)}r^2
 +{2z\over(d+2)(d+4)}r^2                         \tag{15.665.5}
\]

has shell coefficient \(h_s=\operatorname{tr}A_s\).  The complete raw
quartic mass is consequently

\[
 \boxed{\tau_s=\operatorname{tr}R_s=h_s+z\rho_s}.    \tag{15.665.6}
\]

The real PSL decomposition of \(\mathcal Z\) used above is
multiplicity-free.  If constituent \(c\) has dimension \(m_c\), Schur's
lemma makes \(R_s\) scalar there; call the scalar \(q_{s,c}\).  Equations
(15.665.1) and (15.665.6) prove the conserved inequalities

\[
 \boxed{
 q_{s,c}\ge0,
 \qquad \sum_c m_cq_{s,c}=\tau_s,
 \qquad q_{s,c}\le{\tau_s\over m_c}.}              \tag{15.665.7}
\]

Equivalently, every signed harmonic eigenvalue satisfies
\(a_{s,c}=q_{s,c}-\rho_s\).  The smallest constituent has dimension
\(d=(p^2+1)/2\), so (15.665.7) improves the uncoupled rank-one shell bound
by a factor of order \(p^2\) and couples all channels through one scalar
theta coefficient.

At \(p=11\), the channel dimensions are \(1220,305,244\).  On the four
proved nonempty shells \(s=11,20,24,27\), exact evaluation gives

\[
\begin{array}{c|c|c|c|c}
s&N_s&h_s&\rho_s&\tau_s\\ \hline
11&244&-3538/63&2/63&0\\
20&16104&-5368/21&1600/231&923784/77\\
24&14762&-15921/28&64/7&436943/28\\
27&442860&-527406/7&2430/7&538752
\end{array}                                           \tag{15.665.8}
\]

and the channel-weighted nonnegative \(q_{s,c}\) reproduce every listed
\(\tau_s\).  This supplies the nonlinear shell positivity absent from
Proposition 15.641.  It does not alone prove R1: a modular or multiscale
theta inequality must still transport these conserved masses to the
odd-coset target, uniformly in \(p\).

Evidence: src/e1_gmin_m4_prop15665.py,
evidence/e1_gmin_m4_prop15665.json,
evidence/NOTE_2026-08-27_r1_conserved_quartic_shell_mass.md, and
tests/test_prop15665.py.

## Prop 15.666 (2026-08-27) — every finite p=7 size-eight boundary is excluded

After Propositions 15.662--15.664, the complete finite floor census leaves

\[
154056+1194816+1176+69384=1419432                 \tag{15.666.1}
\]

boundaries per product sign, according as they have \((11,16,24,44)\) exact
directional-mean allocations.  The corresponding number of allocation
leaves is

\[
1694616+19117056+28224+3052896=23892792.          \tag{15.666.2}
\]

The common score system has shape \(282\times1225\), with ranks
\((162,147)\) and left-kernel dimensions \((120,135)\) in characteristics
\((3,7)\).
Conditioning those dependencies to vanish on every raised direction block
and selecting forty independent equations gives complete necessary V100
scans.  Their mod-seven and mod-three survivor counts are respectively
458,822 and 2,671,872; intersecting the same
boundary/stratum/allocation triples leaves

\[
 (77616,0,0,103488),\qquad 181104\text{ leaves}.  \tag{15.666.3}
\]

For variable support \(V\) and a restored subset \(T\), dependencies
vanishing on \(V\setminus T\) give an exact necessary catalog join.  Local,
all-triple, and four-positive joins on twenty-two independent mod-seven rows
reduce (15.666.3) through

\[
181104\longrightarrow124745\longrightarrow78126\longrightarrow62892.
                                                               \tag{15.666.4}
\]

The twenty-two digits are stored losslessly because
\(7^{22}=3909821048582088049<2^{64}\); no probabilistic hash is used.

For each of the last 62,892 leaves and each variable direction \(d\), now
condition on all blocks \(V\setminus\{d\}\).  This filters the exact rows of
catalog \(d\), while retaining each row's signature under one common
twenty-two-row full projection.  An exact shared-memory meet-in-the-middle
join then restores all three or five variable catalogs simultaneously.
There are 1,439,451 distinct isolate/full signature pairs.  Exactly 3,777
leaves have an empty isolated catalog, and the complete join rejects the
other 59,115:

\[
 \boxed{3777+59115=62892,\qquad N_{\rm survivor}=0.}              \tag{15.666.5}
\]

The largest hash-build and probe products are 1,764 and 2,744; the exact
table has 4,096 slots, and the separately instrumented capacity-rejection
counter is zero.  Independent CPU/CUDA prefixes agree at every stage,
including 512 leaves for (15.666.5), and the older full
multi-characteristic engine independently returns zero on representatives
of all three residual catalog classes.

The nonsquare conference anti-isometry audited in Propositions
15.662--15.664 bijects the complete finite floor sets and exchanges the two
product signs.  Hence (15.666.5) proves

\[
 \boxed{\text{every finite (p=7) size-eight boundary is impossible for
 both product signs}.}                                      \tag{15.666.6}
\]

The separate infinity-plus-seven profile, residual (ii), Type I, R1,
global QVAR, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15666.py`,
`evidence/e1_gmin_m4_prop15666.json`,
`evidence/NOTE_2026-08-27_p7_size_eight_complete_exclusion.md`,
`evidence/p7_size8_complete/`, and `tests/test_prop15666.py`.

## Prop 15.656 (2026-08-26) — every p=5 four-point boundary is excluded

For (p=5), let (H) be a putative 21-edge residual graph, containing the
distinguished edge, whose odd-degree boundary (D) has size four. Fix an
eigenshell sign (\varepsilon). Quotienting the complete shell by
(ysim-y) leaves 130 distinct edge-sign rows

\[
 f_y(ab)=\varepsilon y_a y_b C_{ab}\in\{\pm1\}.
\]

Every normalized edge column sums to 26. Writing each required score as

\[
 \sum_{e\in H} f_y(e)=3+2A(y),\qquad A(y)\ge0,
\]

and summing over the shell gives

\[
 \sum_y A(y)={21\cdot26-130\cdot3\over2}=78.       \tag{15.656.1}
\]

The boundary and Paley edge-product sign (c_H) prescribe

\[
 (-1)^{A(y)}=-\varepsilon c_H\prod_{v\in D}y_v.    \tag{15.656.2}
\]

Thus (A=P+2L), where (P\in\{0,1\}^{130}), (0\le L_y\le4), and

\[
 \sum_yL_y={78-|P|\over2}.                         \tag{15.656.3}
\]

The edge-count equation, distinguished-edge equation, and 130 equations
for the number of bad edges form a common (132\times325) zero-one matrix
(M_\varepsilon). Exact elimination gives

\[
 \operatorname{rank}_{\mathbb F_5}M_\varepsilon=67,
 \qquad \dim\ker(M_\varepsilon^T)=65.              \tag{15.656.4}
\]

Substitution of the right sides (9-P_y-2L_y) into the 65 left-null
relations yields a bounded congruence system in the lift variables. A
complete square-semilinear orbit scan applies this necessary condition to
all boundary profiles surviving Proposition 15.632. It excludes 712 direct
orbit cases modulo five. The sole timeout has

\[
 c_H=-1,\quad D=\{2,3,12,13\},\quad |P|=56,
 \quad\sum L=11.
\]

An independent reconstruction over (\mathbb F_7) again has rank 67 and
65 left dependencies, and its bounded lift system is exactly infeasible.
Hence all 713 directly scanned orbit cases are excluded.

It remains to transfer the negative no-infinity cases. If (alpha) is a
nonsquare in (\mathbb F_{25}), multiplication by (alpha), together with
switching only infinity, defines a signed permutation satisfying

\[
 S C[\pi,\pi]S=-C.                                 \tag{15.656.5}
\]

It fixes the distinguished edge, exchanges eigenshells, and preserves
normalized scores. Since (|H|=21) is odd and (D) omits infinity,
(deg_H(\infty)) is even, so the Paley edge product changes sign. This is
an exact bijection between the 489 negative and 489 positive no-infinity
orbits, covering 10,925 boundaries on either side.

The direct scans and transfer therefore exclude all 1,202
floor-surviving orbit/sign cases, representing 26,450 boundary/sign cases.
Proposition 15.632 excludes the other 3,450 cases. Consequently

\[
 \boxed{\text{every size-four residual boundary at }p=5\text{ is impossible}.}
                                                               \tag{15.656.6}
\]

Together with Propositions 15.652--15.655, every size-four boundary is
excluded for every odd prime (p\ge5). Boundary size at least six,
residual (ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: src/e1_gmin_m4_prop15656.py,
evidence/e1_gmin_m4_prop15656.json,
evidence/NOTE_2026-08-26_p5_four_point_full_shell.md,
scripts/p5_size_four_full_shell_mod5_batch.py,
scripts/p5_size_four_full_shell_mod7_exception.py,
scripts/p5_size_four_full_shell_audit.py, and tests/test_prop15656.py.

### Proposition 15.657 — every six-point boundary is impossible for `p>=11`

Let `D` be a six-vertex odd-degree boundary at residual size
`|H|=4p+1`, and let `s` be its number of finite points. In projective
direction `d`, write the finite fibre multiplicities as `n_i` and let
`b_d` count the odd multiplicities. Then

\[
 s-b_d=2\sum_i\lfloor n_i/2\rfloor
       \le 2\sum_i\binom{n_i}{2}.                 \tag{15.657.1}
\]

Every finite pair collides in exactly one projective direction. Therefore

\[
 \sum_d(s-b_d)\le s(s-1).                         \tag{15.657.2}
\]

The deficit budget is 30 for six finite boundary points and 20 for infinity
plus five finite points. Exact positive quadrature extends Proposition
15.652's parity floors through `b=6`. For odd `p>=11`, the phase-zero floor
at `b=5,6` is `2p`; the phase-one floor is `2p-4` at `p=11`, `2p-2` at
`p=13`, and `2p` for primes `p>=17`. The dual quadratures are supported on
`1,3,5`, on `2,3,4` for `p<=15`, and on `0,2,4` for `p>=15`, respectively.

If infinity lies in `D`, use the phase-independent bounds

\[
 f(1)\ge p-1,\qquad f(3)\ge2p-6,\qquad f(5)\ge2p-4.
\]

Starting with cost `(p+1)(2p-4)`, the deficit budget can save at most
`5(p-3)`. The remaining excess over the total budget `(p+1)^2` is

\[
 p^2-9p+10>0\qquad(p\ge11).                       \tag{15.657.3}
\]

If all six points are finite and `p>=13`, use

\[
 f(0)\ge0,\quad f(2)\ge p-1,\quad
 f(4)\ge2p-6,\quad f(6)\ge2p-2.
\]

Starting with `(p+1)(2p-2)`, deficit 30 can save at most `10(p-1)`,
leaving excess

\[
 p^2-12p+7>0\qquad(p\ge13).                       \tag{15.657.4}
\]

At `p=11`, the two quadratic direction types have opposite phases and
budget 72 each. The phase-one type requires pair deficit at least 20; the
phase-zero type requires at least 18. Their total 38 contradicts
(15.657.2). Consequently

\[
 \boxed{\text{every six-point residual boundary is impossible for every
 odd prime }p\ge11.}                              \tag{15.657.5}
\]

The exceptional `p=5,7` size-six cases, boundaries of size at least eight,
residual (ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15657.py`,
`evidence/e1_gmin_m4_prop15657.json`,
`evidence/NOTE_2026-08-26_size_six_boundary_exclusion.md`, and
`tests/test_prop15657.py`.

## Proposition 15.668 — exact p=11 broad-channel theta and finite R1

For a quartic value profile \(a=(a_s)_{s\in\mathbb F_{11}}\), retain in
addition to its scalar moments the marked contraction

\[
 U_4(a)=\sum_{c\in\mathbb F_{11}}
 \left(\sum_s\eta(s-c)a_s\right)^4.              \tag{15.668.1}
\]

Input-affine canonicalization leaves 1,007 profile types, while
output-affine transport reduces their dynamic programs to 20 canonical
ones. The complete ten-dimensional glue code then reduces exactly as

\[
 11^{10}\longrightarrow21,437,340
 \longrightarrow2,584,901                                  \tag{15.668.2}
\]

translation/nonzero-scalar representatives and weighted rich-profile
tuples. Their weights sum to \(11^{10}\). Five independent CRT moduli have
product

\[
31999921744068749461247094447450713426945936557,
\]

which exceeds separate unrestricted bounds for every count, squared excess,
and \(U_4\) coefficient through exponent 120. Hence all recovered integers
are exact, rather than probabilistic modular fingerprints.

For each raw positive shell operator \(R_e\), the scalar trace together with
the two marked contractions \(z^Tz\) and \(z^TA_2z\) linearly recovers its
mass on the square-circle kernel, low, and high eigenspaces. Their dimensions
are

\[
 1220,\qquad305,\qquad244,\qquad1220+305+244=1769. \tag{15.668.3}
\]

After subtracting the universal radial term, all three coefficient sequences
lie in 32-dimensional affine modular spaces with common pivot exponents
\(31,32,35,36,\ldots,91,92\). Thus the prefix through 92 determines every
series uniquely. Exponents 93 through 120 are held out; all 28 coefficients
match exactly in each channel. Evaluation through exponent 800 remains
nonnegative, sums coefficientwise to the aggregate raw trace, and obeys the
dimension-weighted transformed-target identity.

Channelwise shell-mass and transformed-target conservation define eight
rational endpoint LPs. Exact QSopt_ex primal and dual certificates give the
harmonic-target intervals

\[
\begin{array}{c|c}
\text{kernel principal}&[-522.933314,508.493608]\\
\text{low Weil}&[-382.405131,392.954765]\\
\text{low principal}&[-257.123541,268.631365]\\
\text{high principal}&[-219.228926,219.228419].
\end{array}                                                \tag{15.668.4}
\]

Every primal constraint and dual stationarity equation is independently
recomputed over the rationals. Under the exact Poisson conversion, every
interval still contains a target with \(\Phi<6\). Therefore this strict
broad-channel cone does not prove R1. This is a certified limit of the
relaxation, not a counterexample.

Independently, the complete finite \(p=11\) Max+ census gives

\[
 \|\delta\|^2={1382747375360\over583792784981},\qquad
 {n\over12}-\|\delta\|^2
 ={27314875631681\over3502756709886}>0.             \tag{15.668.5}
\]

Consequently strong R1, and hence the exact R1 threshold
\(22143/1682\), holds at \(p=11\). No uniform argument follows: general R1,
global QVAR, residual (ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15668.py`,
`evidence/e1_gmin_m4_prop15668.json`,
`evidence/NOTE_2026-08-28_p11_broad_channel_theta.md`,
`tests/test_prop15668.py`, and the 33-file checked archive at
`/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-28-r1-broad-channel/`.

## Proposition 15.669 — uniform non-Walsh boundary-range exclusion

Let \(p\ge17\) be odd and put \(m=(p+1)/2\). For a set
\(B\subseteq\mathbb F_p\) of size \(b\), let \(X\) be uniform on
\(J(p,m)\) and write \(t=|X\cap B|\). Then

\[
 \mu=\mathbb E t={b(p+1)\over2p},\qquad
 \mathbb E t^2={b(b+1)(p+1)\over4p},
\]

\[
 \sigma^2=\operatorname {Var}(t)
 ={(p+1)b(p-b)\over4p^2}.                         \tag{15.669.1}
\]

For a parity phase \(\eta\in\{0,1\}\), recall that \(M(p,b,\eta)\) is the
minimum hypergeometric expectation of a quadratic \(q\) satisfying

\[
 q(t)\ge (t+\eta\bmod2)
\]

at every point in the support. The constant \(q=1\) is feasible. We prove
that it is optimal whenever \(5\le b\le p-5\).

Complementing \(B\) replaces \(t\) by \(m-t\) and changes the phase by
\(m\bmod2\), so it is enough to take \(5\le b\le(p-1)/2\). Put
\(j=1-\eta\), let \(r=j\), and let \(R\le b\) be the largest integer
congruent to \(j\pmod2\). These are the endpoint contact nodes. At mean
\(\mu\), the largest variance of a measure on the contact nodes is

\[
 v_+=(\mu-r)(R-\mu).
\]

The smallest is \(v_-=(\mu-a)(a+2-\mu)\), where \(a,a+2\) are the adjacent
contact nodes bracketing \(\mu\), with \(v_-=0\) when the mean itself is a
contact. Direct substitution in (15.669.1) gives the following four possible
upper margins after multiplication by \(4p\):

\[
\begin{split}
4p(v_+-\sigma^2)\in\{&b(b-1)(p+1),\\
 &b\{b(p+1)-3p+1\},\\
 &b(b-3)(p+1),\\
 &(b-1)\{b(p+1)-4p\}\}.
\end{split}                                           \tag{15.669.2}
\]

All four are positive for \(b\ge5\). According to \(b\pmod4\) and \(j\),
the lower margin is one of

\[
\begin{split}
\sigma^2-v_-\in\{&{b(p-b-3)\over4p},
 {b(p+1-b)-4p\over4p},\\
 &{b(p-1-b)-3p\over4p},
 {(p-b)(b-3)\over4p}\}.
\end{split}                                           \tag{15.669.3}
\]

The first and fourth quantities are positive immediately on the reduced
range. The second occurs only for even \(b\ge6\), and its numerator is a
concave quadratic whose endpoint lower bounds on
\([6,(p-1)/2]\) are

\[
 2p-30,\qquad {p^2-14p-3\over4}.
\]

The third occurs only for odd \(b\ge5\), and its corresponding bounds are

\[
 2p-30,\qquad {p^2-14p+1\over4}.
\]

All are positive at \(p=17\) and increase thereafter. Thus
\(v_-\le\sigma^2\le v_+\). Mix the adjacent-contact two-point measure with
the endpoint two-point measure in the unique rational proportion that gives
variance \(\sigma^2\). The resulting positive measure is supported entirely
where \(t+\eta\) is odd and matches the hypergeometric moments through degree
two. Every feasible quadratic therefore has expectation at least one. Hence

\[
 \boxed{M(p,b,\eta)=1,
 \qquad 2\lceil pM(p,b,\eta)\rceil=2p
 \quad(5\le b\le p-5).}                            \tag{15.669.4}
\]

Now let \(D\) be a residual boundary and let \(s\) denote its number of
finite points. For each projective direction \(d\), let \(b_d\) be the number
of affine fibres meeting \(D\setminus\{\infty\}\) oddly. If the fibre
multiplicities are \(n_i\), then

\[
 s-b_d=2\sum_i\lfloor n_i/2\rfloor
 \le2\sum_i\binom{n_i}{2}.
\]

Every finite pair lies in exactly one projective direction, so

\[
 \sum_d(s-b_d)\le s(s-1).                          \tag{15.669.5}
\]

First suppose infinity is absent. Then \(s\) is even and the two quadratic
direction types have opposite phases. Each type contains \(m\) directions,
has exact budget \(m(p+1)\), and must save

\[
 R=2pm-m(p+1)={p^2-1\over2}                       \tag{15.669.6}
\]

from the middle-floor baseline. In phase zero, only the low values
\(b=0,2,4\) can improve that baseline in the stated range; their
(saving, deficit) pairs are

\[
 (2p,s),\qquad(p-1,s-2),\qquad(6,s-4).
\]

For \(p\ge17,s\ge6\),

\[
 s(p+1)>4p,\qquad s(p-7)>4(p-4),                  \tag{15.669.7}
\]

so \(b=0\) has the best saving per deficit and \(b=2\) dominates \(b=4\)
afterward. Let \(a=\lfloor R/(2p)\rfloor\) and \(r=R-2pa\). If at least
\(a+1\) zero-fibre directions occur, their deficit is at least
\((a+1)s\). Otherwise the fractional relaxation using \(a\) such directions
and then \(b=2\) gives at least

\[
 as+{r(s-2)\over p-1}.
\]

For \(p\equiv1\pmod4\),
\(a=(p-1)/4,r=(p-1)/2\); for \(p\equiv3\pmod4\),
\(a=(p-3)/4,r=(3p-1)/2\). Both branches imply

\[
 D_0\ge{(p+1)s\over4}-1.                          \tag{15.669.8}
\]

In phase one, only \(b=2\) saves from the middle baseline, by \(p+1\).
Exactly \((p-1)/2\) such savings are needed in the relaxed problem, hence

\[
 D_1\ge{(p-1)(s-2)\over2}.                        \tag{15.669.9}
\]

Subtracting (15.669.5) from (15.669.8)--(15.669.9), four times the
contradiction gap is

\[
 h_p(s)=s(3p+3-4s)-4p.                            \tag{15.669.10}
\]

This quadratic is concave. At the two endpoints of
\([6,3(p-1)/4]\),

\[
 h_p(6)=14p-126>0,\qquad
 h_p(3(p-1)/4)={p-9\over2}>0.
\]

Therefore

\[
 \boxed{\text{every all-finite even boundary with }
 6\le s\le3(p-1)/4\text{ is impossible}.}         \tag{15.669.11}
\]

Suppose instead that infinity lies in \(D\). Then \(s\) is odd and both
quadratic types have the same phase. In phase zero, \(b=1\) saves \(p-1\)
and every other allowed \(b\le p-4\) saves at most six. Since the target
(15.669.6) is exactly \(m(p-1)\), every direction must have \(b=1\). The two
types therefore require total deficit \((p+1)(s-1)\).

In phase one, \(b=1\) saves \(p+1\), while every other allowed value saves at
most six. This statement includes the endpoint exception
\(b=p-4,p\equiv1\pmod4\), whose complemented \(b=4\) floor does save six.
If at most \(m-2\) directions of one type had \(b=1\), even assigning saving
six to each of the other two would give at most

\[
 (m-2)(p+1)+12,
\]

which is short of \(m(p-1)\) by \(p-11>0\). Thus each type has at least
\(m-1=(p-1)/2\) directions with \(b=1\), for total deficit at least
\((p-1)(s-1)\). If \(5\le s\le p-4\), both phase bounds strictly exceed
\(s(s-1)\). Consequently

\[
 \boxed{\text{every infinity-present boundary with odd finite }
 5\le s\le p-4\text{ is impossible}.}             \tag{15.669.12}
\]

For \(p=11,13\), the verifier evaluates the original rational
three-variable parity-majorant LP for every allowed \(b\), then performs an
exact dynamic program over the \((p+1)/2\) directions of each type. Relative
to (15.669.5), it gives phase gaps \(30,18\) for \(p=11\), infinity plus
seven; opposite-phase gap \(4\) for \(p=13\), eight finite; and phase gaps
\((42,30)\), \((40,24)\) for \(p=13\), infinity plus seven and nine. Every
gap is positive, so all these cases are impossible.

This proposition is a range theorem, not a closure of residual (ii). The
first floor-plus-pair survivors are eight finite or infinity plus nine at
\(p=11\), ten finite or infinity plus eleven at \(p=13\), the first even
integer strictly above \(3(p-1)/4\) without infinity for \(p\ge17\), and
\(p-2\) finite points with infinity. Such a survivor is only a directional
count profile, not an affine boundary or residual graph. General residual
(ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: src/e1_gmin_m4_prop15669.py,
evidence/e1_gmin_m4_prop15669.json,
evidence/NOTE_2026-08-28_uniform_boundary_range_exclusion.md, and
tests/test_prop15669.py.

## Proposition 15.670 — every finite p=11 size-eight boundary is impossible

Let \(D\subset\mathbb F_{11^2}\) have eight points. For each projective
\(\mathbb F_{11}\)-direction \(d\), let \(b_d\) be the number of parallel
fibres meeting \(D\) oddly, and let \(\epsilon_d\in\{-1,1\}\) be the
quadratic type of the direction kernel. There are six directions of each
type. At residual size \(4p+1=45\), Proposition 15.632 gives each type the
exact budget

\[
 {(p+1)^2\over2}=72.                              \tag{15.670.1}
\]

For a finite even boundary the phase is
\(\eta_d=\mathbf1_{\epsilon_d=c_H}\). The exact parity-majorant floors are

\[
\begin{array}{c|rrrrr}
b&0&2&4&6&8\\ \hline
f_0(b)&0&12&16&22&16\\
f_1(b)&22&10&22&18&22.
\end{array}                                        \tag{15.670.2}
\]

Choose an ordered pair \(x\ne y\) in \(D\). The unique affine similarity

\[
 z\longmapsto {z-x\over y-x}
\]

sends it to \(0,1\). It therefore suffices for exclusion to test the
\(\binom{119}{6}=3,470,108,187\) eight-sets containing those two points.
This normalization is lossless: multiplication by \(a\ne0\) sends every
direction type to \(\chi(a)\epsilon_d\); transferring \(c_H\) to
\(\chi(a)c_H\) preserves the phase, and the two type budgets are equal. The
pointed-set identity

\[
 \binom{121}{8}\,8\cdot7
 =\binom{119}{6}\,121\cdot120                     \tag{15.670.3}
\]

provides an independent exact coverage ledger.

An exhaustive direct-rank census evaluates both \(c_H\) signs on every
normalized set and accumulates the complete pair of type-cost histograms.
For both signs there are zero pairs with both costs at most 72. More sharply,

\[
 \min_D\max\left\{
  \sum_{\epsilon_d=-1}f_{\eta_d}(b_d),
  \sum_{\epsilon_d=1}f_{\eta_d}(b_d)
 \right\}=76.                                      \tag{15.670.4}
\]

The \(c_H=-1\) first minimizer has costs \((76,66)\), and the
\(c_H=1\) first minimizer has costs \((64,76)\). Thus (15.670.1) is
violated by at least four for every boundary, proving

\[
 \boxed{\text{every finite }p=11\text{ size-eight boundary is
 impossible}.}                                    \tag{15.670.5}
\]

The complete scan was independently replayed on a Tesla V100 under CUDA and
an RX 9070 XT under ROCm/HIP. Both full histograms agree bit for bit, and an
independent CPU combinations traversal agrees on every entry in a 100,000-set
prefix. The verifier additionally checks all 2,892 affine
translation/scalar direction actions needed by the normalization and
recomputes the minimizing costs with the generic rational floor routine.

This closes only the finite-eight \(p=11\) branch. Infinity plus nine and
larger \(p=11\) boundaries remain, as do general residual (ii), R1, global
QVAR, Type I, and the limit.

Evidence: `src/e1_gmin_m4_prop15670.py`,
`scripts/p11_size8_normalized_floor_gpu.py`,
`evidence/e1_gmin_m4_prop15670.json`,
`evidence/p11_size8_normalized_floor_v100.json`,
`evidence/p11_size8_normalized_floor_rx9070xt.json`,
`evidence/NOTE_2026-08-28_p11_size_eight_boundary_exclusion.md`, and
`tests/test_prop15670.py`.

## Proposition 15.671 — rigid-sign collinear near-line exclusion

Suppose the boundary is infinity together with (p-2) finite points on one
affine line. The line direction has (b=1), while each transverse direction
has (b=p-2). In the rigid product-sign branch, the exact endpoint floors
use the full budget of each quadratic direction type. When
(p\equiv1\pmod4), the line type has apparent surplus two, but Proposition
15.642 makes every nonzero lift cost at least four for (p\ge13). Thus every
direction is at its parity baseline.

Put (q=(p-1)/2), let (I) be the infinity-edge count, and let (P_d) be
the finite parallel-edge count. With (\sigma=(-1)^\eta), the special and
transverse targets are

\[
 \epsilon_0S_H=4+\sigma z_j,
 \qquad \epsilon_dS_H=4-z_az_b.
\]

The additive inter-fibre coefficient identity gives

\[
 P_0\equiv4+\sigma-I\pmod q,
 \qquad P_d\equiv4-I\pmod q,                         \tag{15.671.1}
\]

and

\[
 2I+(p+1)P_0\le8p+1+\sigma,
 \qquad2I+(p+1)P_d\le8p+2.                          \tag{15.671.2}
\]

Summing (15.671.1) over all directions and using
(\sum_dP_d=4p+1-I) gives

\[
 I\equiv3+\sigma\pmod q.                            \tag{15.671.3}
\]

If (p\equiv1\pmod4), then (q) is even and (\sigma=-1), so
(15.671.3) makes (I) even, contrary to infinity being in the odd-degree
boundary. If (p\equiv3\pmod4), write

\[
 I=4+qk_0,qquad P_0=1+qa_0,qquad P_d=qa_d.
\]

The edge count gives (k_0+a_0+\sum a_d=8), while (15.671.2) gives

\[
 k_0+(q+1)a_0\le7,qquad k_0+(q+1)a_d\le8.
\]

For (p\ge19), (q\ge9), so every (a) vanishes and the edge count forces
(k_0=8). But (I) odd requires (k_0) odd. Hence the rigid-sign
collinear branch is impossible for (p\equiv1\pmod4,p\ge13) and for
(p\equiv3\pmod4,p\ge19).

Evidence: `src/e1_gmin_m4_prop15671.py`,
`evidence/NOTE_2026-08-28_near_line_rigid_sign.md`, and
`tests/test_prop15671.py`.

## Proposition 15.672 — complete collinear near-line exclusion

Take the opposite product sign on the same collinear boundary. The
transverse (b=p-2) baseline is the xnor of the two omitted fibres and has
scaled mean (p-1). Four nonzero lifts cost at least

\[
 {p^2-1\over p-2}>p+1,
\]

more than either type surplus, so each type retains a transverse baseline.
For a direction of type (\epsilon_d), the exact mean is

\[
 a_d=I+(p+1)P_d-\epsilon_dT-3p.                    \tag{15.672.1}
\]

Same-type means therefore differ by multiples of (p+1). The split budget
leaves exactly one exception per type, each with mean (2p) and parallel
count one above its type baseline. If (x,y) are the two baseline counts
and (m=(p+1)/2), then

\[
 E=m(x+y)+2,qquad I=4p-1-m(x+y),qquad x+y\le7.     \tag{15.672.2}
\]

Every type has a transverse xnor baseline, whose coefficient equation gives

\[
 q\mid I+x-4,qquad q\mid I+y-4.
\]

Substituting (15.672.2), (p=2q+1), and (m=q+1) gives

\[
 q\mid y+1,qquad q\mid x+1.
\]

Thus (x+y\ge2q-2\ge8), contradicting (15.672.2). Combined with
Proposition 15.671, both product signs of every collinear
infinity-plus-((p-2)) boundary are excluded for every prime (p\ge13).

Evidence: `src/e1_gmin_m4_prop15672.py`,
`evidence/NOTE_2026-08-28_near_line_complete.md`, and
`tests/test_prop15672.py`.

## Proposition 15.673 — complete endpoint-only near-line exclusion

Retain infinity and (s=p-2) finite boundary points, but assume only that

\[
 b_d\in\{1,p-2\}\quad\hbox{for every direction }d. \tag{15.673.1}
\]

Equation (15.672.1) implies that the (m=(p+1)/2) means of one quadratic
type have the form

\[
 a_d=r+(p+1)k_d,qquad r=2u,qquad\sum_dk_d=m-u.     \tag{15.673.2}
\]

The two endpoint floors are

\[
\begin{array}{c|c|cc}
p\bmod4&\eta&b=1&b=p-2\\ \hline
3&0&p+1&p+1\\
1&1&p-1&p+1\\
3&1&p-1&p-1\\
1&0&p+1&p-1.
\end{array}                                        \tag{15.673.3}
\]

A nonzero lift costs at least four, so a value only two above its baseline
is impossible. Equations (15.673.2)--(15.673.3) consequently leave only
pure baseline types and, in the nonsaturated rows, one (p+1)-unit mean
jump per relevant type.

The endpoint geometry removes the degenerate combinations. Let (R) count
directions with (b_d=1). The pair-deficit inequality gives (R\le s).
If equality holds, every affine line contains at most two boundary points,
so the finite boundary is a ((p-2))-arc and the other three directions are
undetermined. Adjoin any two of their points at infinity. This gives a
(p)-arc in (\mathrm{PG}(2,p)), hence lies on a conic by Segre's
odd-order (p)-arc theorem. Two choices give conics sharing the
(p-2\ge5) finite points, so they coincide and contain all three collinear
infinity points, impossible. If at most two directions are determined, the
finite set is collinear and Propositions 15.671--15.672 apply.

It remains that \(3\le R\le p-3=2m-4\). The same-type residue equation
has the following exhaustive forms, where \(B\) denotes a type whose
baseline directions have \(b=1\), and \(C\) one whose baseline directions
have \(b=p-2\):

\[
\begin{array}{c|cc|c}
(p\bmod4,\eta)&B\text{-type }R&C\text{-type }R&
 \text{type edge offsets}\\ \hline
(3,0)&m&0&(0,0)\\
(1,1)&m-1\text{ or }m&0&(1,0)\\
(3,1)&m-1\text{ or }m&0\text{ or }1&(1,1)\\
(1,0)&m&0\text{ or }1&(0,1).
\end{array}                                           \tag{15.673.4}
\]

For the equal-floor rows, baseline directions of both endpoint kinds
cannot occur in one type: equal means force equal \(P_d\), while their
coefficient congruences differ by one modulo \(q\). For the unequal-floor
rows, the table follows directly from (15.673.2): a type is either at the
saturated endpoint baseline, or has \(m-1\) low baselines and one value one
period higher. Two \(B\) types violate \(R\le2m-4\), while two \(C\) types
violate \(R\ge3\). Hence every branch has one type of each kind, with total
edge offsets \(0,1,2,1\).

Let (x) and (y) denote the baseline parallel counts in the (b=1) and
(b=p-2) types. Comparing the baseline coefficients and substituting the
global edge count gives

\[
\begin{array}{c|c|c}
(p\bmod4,\eta)&E&\text{congruences}\\ \hline
(3,0)&m(x+y)&q\mid y,\ q\mid x-1\\
(1,1)&m(x+y)+1&q\mid x,\ q\mid y-1\\
(3,1)&m(x+y)+2&q\mid y,\ q\mid x+1\\
(1,0)&m(x+y)+1&q\mid x,\ q\mid y+1,
\end{array}                                        \tag{15.673.5}
\]

where (q=(p-1)/2), (I=4p+1-E\ge1), and hence (x+y\le7). For
(p\equiv3\pmod4), (q\ge9): phase zero forces ((x,y)=(1,0)), but then

\[
 I=7q+4>4q+1=s+2E,
\]

while phase one has no congruence solution. For (p\equiv1\pmod4), phase
one forces ((x,y)=(0,1)), but

\[
 I=7q+3>4q+3=s+2E.
\]

Phase zero has no solution for (q>8). Here (I\le s+2E) is elementary:
the (I) infinity edges toggle (I) finite vertices, the (E) finite
edges toggle at most (2E), and their symmetric difference is the
(s)-point finite boundary.

It remains to treat (p=17,q=8), where (15.673.5) has the sole candidate

\[
 (x,y)=(0,7),\qquad E=64,\qquad I=5.                \tag{15.673.6}
\]

For a complementary baseline direction, let (n_s\ge0) be its infinity
star counts. Then (\sum n_s=5), (P_d=7), and its signed inter-fibre
matrix is

\[
 L_{st}=1-n_s-n_t+\mathbf1_{st=ab}.                 \tag{15.673.7}
\]

If (u) of the 17 counts are positive and (z=17-u), the norm before the
distinguished correction is

\[
 \binom z2+z(5-u)+5(u-1)-\binom u2.
\]

The correction lowers it by at most one. For (u=1,\ldots,5), the lower
bounds are (183,153,125,99,75); the last is attained by five unit counts
including (a,b). Thus the exact (\ell_1) minimum is 75. Only
(E-P_d=57) transverse selected edges are available, contradicting
(15.673.7). Therefore

\[
 \boxed{\text{every endpoint-only infinity-plus-}(p-2)\text{ boundary is
 impossible for both signs and every prime }p\ge17.}
\]

Intermediate odd-fibre counts and the large all-finite range remain open,
so residual (ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15673.py`,
`evidence/e1_gmin_m4_prop15673.json`,
`evidence/NOTE_2026-08-28_endpoint_near_line_complete.md`, and
`tests/test_prop15673.py`.

## Proposition 15.674 — complete infinity-plus-\((p-2)\) shell exclusion

Let the odd-degree boundary consist of infinity and (s=p-2) finite
points, with no restriction on its directional odd-fibre counts. Put

\[
 P=p+1=2m,\qquad q={p-1\over2}.
\]

Every (b_d) is odd. Proposition 15.669's complete floor table gives

\[
\begin{array}{c|c|cc}
p\bmod4&\eta&b=1&b=s\\ \hline
3&0&P&P\\
1&1&P-2&P\\
3&1&P-2&P-2\\
1&0&P&P-2.
\end{array}                                           \tag{15.674.1}
\]

For every intermediate (3\le b\le p-4), the floor is either
(2P-8) or (2P-2). Thus all floors are at least (P-2), and every
intermediate floor is strictly larger than (P).

Within one quadratic type, (15.672.1) and the exact type sum give

\[
 a_d=2u+Pk_d,qquad 0\le u<m,qquad \sum_d k_d=m-u.
                                                            \tag{15.674.2}
\]

If (1\le u\le m-2), every (k_d\ge1), contradicting
(sum k_d=m-u<m). If (u=0), every (k_d\ge1). A floor-((P-2))
direction cannot take value (P), since that would be a forbidden
two-unit lift, and an intermediate direction has floor larger than (P).
Therefore every direction has floor and mean (P). If (u=m-1), then
(sum k_d=1): exactly (m-1) directions have floor and mean (P-2),
and the remaining arbitrary direction has mean (2P-2=2p). Consequently
each type contains at most one intermediate odd-fibre count.

When the two endpoint kinds share one baseline floor, they cannot both occur
among the baseline directions of one type. Their equal means would force
equal (P_d), whereas the coefficient identities give

\[
 q\mid I+P_d-(4+\sigma)\quad(b=1),\qquad
 q\mid I+P_d-4\quad(b=s),                              \tag{15.674.3}
\]

which differ by one modulo (q). Hence every baseline type is homogeneous.

Call the two possible baseline kinds (B) for (b=1) and (C) for
(b=s). Two (B) types contain at least (2(m-1)=p-1>s) baseline
directions. Their contribution to the pair deficit alone exceeds
(s(s-1)), impossible. Two (C) types leave at most their two exceptional
directions determined by pairs of finite boundary points. A noncollinear
set determines at least three directions, so the finite boundary would be
collinear; Propositions 15.671--15.672 exclude that case. Thus one type has
each baseline kind.

The four type pairs and their numbers of exceptional parallel edges are

\[
\begin{array}{c|cc|c}
(p\bmod4,\eta)&B\text{ baselines/exceptions}&
C\text{ baselines/exceptions}&\text{total offset}\\ \hline
(3,0)&m/0&m/0&0\\
(1,1)&(m-1)/1&m/0&1\\
(3,1)&(m-1)/1&(m-1)/1&2\\
(1,0)&m/0&(m-1)/1&1.
\end{array}                                           \tag{15.674.4}
\]

Let (x,y) be the baseline parallel counts of the (B,C) types.
Equations (15.674.3), the edge count, and (15.674.4) reproduce exactly

\[
\begin{array}{c|c|c}
(p\bmod4,\eta)&E&\text{congruences}\\ \hline
(3,0)&m(x+y)&q\mid y,\ q\mid x-1\\
(1,1)&m(x+y)+1&q\mid x,\ q\mid y-1\\
(3,1)&m(x+y)+2&q\mid y,\ q\mid x+1\\
(1,0)&m(x+y)+1&q\mid x,\ q\mid y+1.
\end{array}                                           \tag{15.674.5}
\]

Since (I=4p+1-E\ge1), (x+y\le7). As in Proposition 15.673, the
first and second rows have only ((1,0)) and ((0,1)), respectively,
and both violate (I\le s+2E); the third has no candidate; and the fourth
has no candidate for (q>8).

For (p=17,q=8), the last row leaves

\[
 (x,y)=(0,7),\qquad E=64,\qquad I=5.
\]

The (C) type has at least eight complementary baseline directions. For
any one of them, Proposition 15.673's exact matrix calculation still gives
entrywise norm at least (75), against only (E-P_d=57) transverse
edges. It does not involve the exceptional direction's fibre count. Hence
this final candidate is also impossible, and

\[
 \boxed{\text{every infinity-plus-}(p-2)\text{ boundary is impossible for
 both signs and every odd prime }p\ge17.}
\]

The infinity-plus-(p) shell, the large all-finite range, general residual
(ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15674.py`,
`evidence/e1_gmin_m4_prop15674.json`,
`evidence/NOTE_2026-08-28_full_near_line_shell_complete.md`, and
`tests/test_prop15674.py`.

## Proposition 15.675 — first all-finite survivor in two residue classes

Let (p\ge19), put (P=p+1=2m), and let (s) be the first even integer
strictly above (3(p-1)/4). For one quadratic direction type, write

\[
 a_d=2u+Pk_d,qquad 0\le u<m,qquad \sum_d k_d=m-u.
                                                            \tag{15.675.1}
\]

At this boundary size the phase-one floors are (P-2) at (b=2) and
(2P-2) at (b=0,4), and throughout (6\le b\le s). If
(u\le m-2), every (k_d\ge1), contrary to (15.675.1); at (u=0),
the candidate values over both floor classes are additionally forbidden
two-unit lifts. Hence (u=m-1), with (m-1) directions at (b=2) and
one at (b=s). Its exact minimum deficit is

\[
 D_1=(m-1)(s-2).                                      \tag{15.675.2}
\]

In phase zero, for (u\ge2), the deficit-optimal quotient weights at
(b=0,2,s) are (0,1,2). If (t\in\{0,1\}) is congruent to (m+u)
modulo two, the minimum is

\[
 D_0(u)={s(m+u+t)\over2}-2t.                         \tag{15.675.3}
\]

The successive differences alternate between (s-2) and (2), so
(u=2) is minimal. Therefore

\[
 D_0=\begin{cases}
  s(m+2)/2,&m\text{ even},\\
  s(m+3)/2-2,&m\text{ odd}.
 \end{cases}                                         \tag{15.675.4}
\]

The residue (u=1) is impossible because value two over the zero floor is
forbidden. At (u=0), quotient weights (0,1,2,3) occur at
(b=0,2,4,s). Writing (m=3h+r), (0\le r\le2), its exact minimum is
((m-h)s-2r), strictly larger than (15.675.4) for (m\ge10). Thus
(15.675.2)--(15.675.4) are the exact minima of the relaxed type problem.

The first survivor in the four classes is

\[
\begin{array}{c|cccc}
p\bmod8&1&3&5&7\\ \hline
s&(3p+5)/4&(3p-1)/4&(3p+1)/4&(3p+3)/4.
\end{array}
\]

Subtracting the pair-deficit budget (s(s-1)) gives

\[
 D_0+D_1-s(s-1)=
 \begin{cases}
  -(p-1)/4,&p\equiv1\pmod8,\\
  (p+1)/2,&p\equiv3\pmod8,\\
  (p-1)/2,&p\equiv5\pmod8,\\
  -(p-7)/4,&p\equiv7\pmod8.
 \end{cases}                                         \tag{15.675.5}
\]

The positive middle rows contradict the pair budget. Hence the first even
all-finite boundary size above (3(p-1)/4) is impossible for every prime
(p\ge19) with (p\equiv3) or (5\pmod8). The outer rows remain open.

Evidence: `src/e1_gmin_m4_prop15675.py`,
`evidence/e1_gmin_m4_prop15675.json`,
`evidence/NOTE_2026-08-28_first_all_finite_survivor_half_close.md`, and
`tests/test_prop15675.py`.

## Proposition 15.676 — infinity-plus-p pair equality is impossible

Let the odd boundary contain infinity and (p) finite points, with
(p\ge17). Pair-deficit equality makes every affine line contain at most two
finite boundary points, hence those points form a (p)-arc. Segre's odd-order
(p)-arc theorem puts the arc on a conic.

The line at infinity cannot be secant, because then only (p-1) conic points
are affine. If it is tangent, the affine odd-fibre profile is

\[
 p\text{ copies of }b=1,\qquad 1\text{ copy of }b=p.
\]

If it is external, deleting one of the (p+1) affine conic points gives

\[
 m+1\text{ copies of }b=1,\qquad m-1\text{ copies of }b=3,
 \qquad m=(p+1)/2.
\]

In phase zero, one (b=3) direction already exceeds a type budget; in phase
one each type contains at most one, while the profile has (m-1\ge8)
globally. Thus the external case is impossible. In the tangent case, exact
same-type residues retain a baseline (b=1) direction in each type. The
baseline coefficient congruences force both parallel counts to vanish,
leaving zero or two finite edges and violating boundary support. Therefore
pair-deficit equality is impossible for both signs and every prime
(p\ge17). Strict pair deficit remains open.

Evidence: `src/e1_gmin_m4_prop15676.py`,
`evidence/e1_gmin_m4_prop15676.json`,
`evidence/NOTE_2026-08-28_infinity_plus_p_arc_close.md`, and
`tests/test_prop15676.py`.

## Proposition 15.677 — first all-finite survivor closed from p=19

For the two modulo-eight classes left by Proposition 15.675, exact quotient
arithmetic from (p\ge23) leaves phase-zero residue (u_0=2), and additionally
(u_0=3) when (p\equiv1\pmod8). In either case

\[
 \sum_d k_d=m-u_0<m,
\]

so one phase-zero direction has quotient zero. Its scaled mean is four or
six, below the phase-zero floor at every nonzero even fibre count, and hence
its odd-fibre count is zero. Pointwise its nonnegative slack factors as
(A_d=2B_d), with (B_d) a nonzero nonnegative integer-valued quadratic.
Proposition 15.642 gives

\[
 4p\,\mathbb E[B_d]\ge8,
\]

contradicting the mean four or six. Combined with 15.675, this excludes the
first even all-finite boundary size above (3(p-1)/4) for every prime
(p\ge19). The attempted treatment of the additional (p=17,u_0=0) row in
15.678 is retracted; Proposition 15.721 later excludes that boundary by an
independent signed-transport argument.

Evidence: `src/e1_gmin_m4_prop15677.py`,
`evidence/e1_gmin_m4_prop15677.json`,
`evidence/NOTE_2026-08-28_first_all_finite_survivor_complete_from_p19.md`,
and `tests/test_prop15677.py`.

## Proposition 15.678 — OPEN_RETRACTED_REDUCTION at the exceptional p=17 row

The former statement that exactly two profiles survive at (p=17,s=14) is
false. After retaining the genuine floor-plus-two equality cells, the exact
census has 108 compatible profiles, including 47 arc profiles, with pair-slack
histogram

\[
 \{0:47,\ 4:32,\ 8:18,\ 12:8,\ 16:3\}.
\]

The local exclusions of (u_0=2), (u_0=3), and the later over-budget residues
remain valid. The unique-16-arc/conic argument from the former proof also
remains valid, but applies to only the fourteen arc profiles having the three
required undetermined directions. It therefore excludes fourteen profiles
and leaves 94 compatible profiles uncovered; it does not prove the historical
endpoint theorem.

Accordingly 15.678 has status **OPEN_RETRACTED_REDUCTION**. The failed
two-profile payload is retained as historical evidence rather than used as a
theorem certificate. Proposition 15.721 independently supersedes this
all-finite boundary as a live gate and excludes it by signed transport.

Evidence: `src/e1_gmin_m4_prop15678.py`,
`evidence/e1_gmin_m4_prop15678.json`,
`evidence/e1_gmin_m4_prop15678.historical_retracted.json`,
`evidence/NOTE_2026-08-28_p17_first_all_finite_survivor_exclusion.md`, and
`tests/test_prop15678.py`.

## Proposition 15.679 — next all-finite boundary closed from p=43

Let (p\ge43), put (P=p+1=2m), and let (s) be the second even integer
strictly above (3(p-1)/4). In the four classes modulo eight,

\[
s={3p+13\over4},\ {3p+7\over4},\ {3p+9\over4},\ {3p+11\over4},
\]

respectively, and always (s\le p-5). Phase one has only residue
(u_1=m-1), with exact deficit

\[
D_1=(m-1)(s-2).                                      \tag{15.679.1}
\]

For phase zero and (2\le u\le m-5), put (t=m-u),
(x=\lfloor t/2\rfloor), (y=t\bmod2), and (z=m-x-y). The exact minimum is

\[
D_0(u)=zs+y(s-2),                                    \tag{15.679.2}
\]

using (x) quotient-two directions at (b=s), (y) quotient-one directions
at (b=2), and (z) quotient-zero directions at (b=0). These minima increase
strictly with (u). Residue zero has minimum

\[
D_0(0)=(m-\lfloor m/3\rfloor)s-2(m\bmod3),            \tag{15.679.3}
\]

and is over the pair budget; residue one is infeasible. Substitution at
(u=8) gives positive gaps

\[
{p+11\over2},\ {5p+25\over4},\ {5p+23\over4},\ {p+13\over2}
\]

in the four modulo-eight classes. The final four residues have deficit at
least ((m-4)s), also over budget. Hence only (2\le u\le7) remain.

For each such row, (\sum k_d=m-u<m), so one phase-zero direction has
quotient zero and mean (2u\le14). Its fibre count is (b=0), and its
nonnegative slack is (A_d=2B_d) for a nonzero nonnegative integer-valued
quadratic on the middle slice. Therefore

\[
2u=4p\,\mathbb E B_d\ge {p^2-1\over4(p-2)}.           \tag{15.679.4}
\]

The right side exceeds 14 for (p\ge59). At the only smaller in-scope
primes (43,47,53), exact pair ledgers leave maximum residues (4,6,5), while
the exact lift floors are (12,14,14), respectively. Thus all rows are
impossible, and the second all-finite boundary is excluded for every prime
(p\ge43). This proposition leaves the same boundary at
(p=17,19,23,29,31,37,41) open; Proposition 15.680 next closes (p=37).
The other six endpoints, later sizes, and general residual (ii) remain open.

Evidence: `src/e1_gmin_m4_prop15679.py`,
`evidence/e1_gmin_m4_prop15679.json`,
`evidence/NOTE_2026-08-28_next_all_finite_boundary_p43.md`, and
`tests/test_prop15679.py`.

## Proposition 15.680 — the p=37 next all-finite endpoint is impossible

At (p=37), the second even all-finite size above (3(p-1)/4) is (s=30).
Phase one has its unique profile

\[
 u_1=18,\qquad D_1=504,
\]

and the exact phase-zero quotient/floor replay leaves only

\[
\begin{array}{c|rrrr}
u_0&2&3&4&5\\ \hline
D_0&328&330&358&360\\
\hbox{pair slack}&38&36&8&6.
\end{array}                                             \tag{15.680.1}
\]

Every retained row has quotient sum (19-u_0<19), so some direction has
quotient zero. Its scaled mean is (2u_0\le10), below the least positive
fibre floor (38), and hence its fibre count is zero. The pointwise slack is

\[
 A=2B,\qquad 2u_0=4p\,\mathbb E B,                      \tag{15.680.2}
\]

for a nonzero nonnegative integer-valued quadratic (B) on (J(37,19)).
Proposition 15.642 excludes (u_0=2,3,4). It remains to rule out the sharp
case

\[
 \mathbb E B={5\over74}.                               \tag{15.680.3}
\]

First, a general paired-cube lemma says that every nonzero Boolean
quadratic (f) on (J(p,(p+1)/2)) satisfies

\[
 \mathbb E f\ge {p-3\over4p}.                          \tag{15.680.4}
\]

Indeed, write (p=2m-1). Through a fixed middle set (X), leave one point of
(X) unpaired, biject the other (m-1) points with the complement, and choose
one endpoint of every pair. These choices form a Boolean ((m-1))-cube. If
(T) averages over this construction and (\rho=1/(2m)=1/(p+1)), direct
calculation on (1,x_i,x_ix_j) gives

\[
 T f=\rho f+(1-\rho)\mathbb E f                         \tag{15.680.5}
\]

for every quadratic. If Boolean (f(X)=1), every cube restriction is a
nonzero degree-two cube polynomial and therefore has support density at
least (1/4). Thus (\rho+(1-\rho)\mathbb E f\ge1/4), which is (15.680.4).

Now apply the stabilizer identity of 15.642 at every point where (B=h).
Its endpoint weight at (p=37) is (9/370), so (15.680.3) forces (h\le2).
The degree-two slice-distance floor gives

\[
 \Pr(B\ne0)\ge {171\over2590},
 \qquad
 \Pr(B=2)\le {5\over74}-{171\over2590}={2\over1295}.
                                                               \tag{15.680.6}
\]

If (B=2) anywhere, then (B(B-1)) is a nonzero degree-four polynomial
supported exactly there. The same exact slice-distance lemma gives

\[
 \Pr(B=2)\ge {\binom{29}{15}\over\binom{37}{19}}
 ={1938\over441595}>{2\over1295},                      \tag{15.680.7}
\]

a contradiction. Hence (B) is Boolean. Equation (15.680.4) now gives
(\mathbb E B\ge17/74>5/74), the final contradiction. Therefore every row
in (15.680.1) is impossible and the (p=37,s=30) boundary is closed. The
same boundary at (p=17,19,23,29,31,41), later all-finite sizes, strict
infinity-plus-(p), residual (ii), R1, global QVAR, Type I, and the limit
remain open.

Evidence: `src/e1_gmin_m4_prop15680.py`,
`evidence/e1_gmin_m4_prop15680.json`,
`evidence/NOTE_2026-08-28_p37_next_all_finite_endpoint.md`, and
`tests/test_prop15680.py`.

## Proposition 15.681 — an integral paired-cube lift closes the p=29 endpoint

Let (B) be a nonzero nonnegative integer-valued quadratic on
(J(p,(p+1)/2)), and put (c=4p\mathbb E B). Choose (X) with (B(X)=h\ge1).
The paired-cube operator from 15.680 does not require (B) to be Boolean:
the restriction to every cube through (X) is nonzero of degree at most two,
so its support has density at least (1/4), and integrality makes its average
at least (1/4). Since (\rho=1/(p+1)),

\[
 \rho h+(1-\rho)\mathbb E B\ge {1\over4},
 \qquad c\ge p+1-4h.                                  \tag{15.681.1}
\]

At the same point, the stabilizer identity of 15.642 gives

\[
 c\ge4h\quad(p\equiv3\pmod4),\qquad
 c\ge {4r\over r+1}h\quad(p=4r+1).                   \tag{15.681.2}
\]

A convex combination cancelling (h) between (15.681.1)--(15.681.2) yields

\[
 4p\mathbb E B\ge
 \begin{cases}(p+1)/2,&p\equiv3\pmod4,\\
               (p-1)/2,&p\equiv1\pmod4.
 \end{cases}                                          \tag{15.681.3}
\]

The resulting scaled floors at (p=29,31,37,41) are (14,16,18,20).
Exact second-boundary ledgers leave phase-zero residues

\[
\begin{array}{c|c}
p&u_0\\ \hline
29&0,2,3,4,5\\
31&0,2,3,4,5,6\\
37&2,3,4,5\\
41&0,2,3,4,5,6,7.
\end{array}                                           \tag{15.681.4}
\]

Every positive row forces a quotient-zero (b=0) direction with
(c=2u_0), and (15.681.3) excludes it. Thus (p=37) is reproved and only
residue zero remains at (p=29,31,41).

For (p=29,s=24), the exact residue-zero minima have profiles
(10[b=0]+5[b=24]) and (14[b=2]+[b=24]), with deficits (240) and (308)
against budget (552). Pair slack is a multiple of four. Exhausting the four
available units leaves either a 24-arc with at least four undetermined
directions, or a set with exactly one 3-secant, all other line occupancies
at most two, and six undetermined directions. In the latter case deleting
one triple point gives a 23-arc and preserves all six directions.

Coolsaet--Sticker's exhaustive (\mathrm{PG}(2,29)) classification has ten
classes of 25-arcs and five classes of 26-arcs. Independently,
(\mathrm{PGL}(2,29)) has ten orbits on five-point subsets and five on
four-point subsets of the 30-point projective line. These are exactly the
classes obtained by deleting five or four points from a conic; uniqueness
of the conic through five arc points shows that the orbit counts do not
merge. Hence every 25- and 26-arc is conic-contained.

Adjoin two undetermined infinity points to a 24-arc to obtain a 26-arc, or
to the 23-arc in the one-triple case to obtain a 25-arc. Taking two such
extensions through three distinct undetermined points gives conics sharing
at least 23 affine points, hence the same conic. It would contain three
collinear points at infinity, impossible. Therefore the (p=29,s=24)
boundary is excluded. At this stage (p=31,41) retain only residue zero;
Propositions 15.682--15.683 subsequently exclude both. The endpoints
(p=17,19,23), later all-finite sizes, strict infinity-plus-(p), residual
(ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15681.py`,
`evidence/e1_gmin_m4_prop15681.json`,
`evidence/NOTE_2026-08-28_p29_next_all_finite_endpoint.md`, and
`tests/test_prop15681.py`.

## Proposition 15.682 — the p=31 next all-finite endpoint is impossible

At (p=31,s=26), Proposition 15.681's scaled lift floor is (16). The exact
pair ledger leaves phase-zero residues (u_0=0,2,3,4,5,6). Every positive
row forces a quotient-zero (b=0) lift of scaled mass (2u_0\le12), so only
(u_0=0) remains.

The residue-zero type minima are

\[
10[b=0]+[b=2]+5[b=26],\quad D_0=284,
\qquad
15[b=2]+[b=26],\quad D_1=360,                         \tag{15.682.1}
\]

against pair budget (26\cdot25=650). Exact enumeration through the six
available deficit units leaves fourteen phase-labelled profiles. Since
pair slack is a nonnegative multiple of four, every profile has slack zero
or four. The eleven slack-zero profiles are 26-arcs with at least three
undetermined directions. In each of the three slack-four profiles there is
exactly one 3-secant, every other line has occupancy at most two, and at
least five directions are undetermined. Deleting one point of the unique
triple gives a 25-arc preserving those directions.

Coolsaet's exhaustive classification of complete arcs in
(\mathrm{PG}(2,31)) has no complete arc of size (23) through (31). Every
finite arc extends greedily to a complete arc, so every 27- or 28-arc must
extend to a complete 32-arc. Segre's odd-order theorem identifies that arc
as a nondegenerate conic; hence every 27- and 28-arc is conic-contained.

For a slack-zero profile, adjoining two undetermined infinity points gives
a 28-arc. For a slack-four profile, perform the one deletion and then
adjoin two such points to get a 27-arc. Choose three undetermined infinity
points (D_1,D_2,D_3). The conics supplied by the extensions through
({D_1,D_2}) and ({D_1,D_3}) share at least 25 affine points and therefore
coincide. One nondegenerate conic would contain three collinear infinity
points, a contradiction. Thus (p=31,s=26) is excluded. Proposition 15.683
subsequently excludes the residue-zero (p=41) endpoint. This boundary
remains open at (p=17,19,23); later all-finite sizes, strict
infinity-plus-(p), residual (ii), R1, global QVAR, Type I, and the limit
remain open.

Evidence: `src/e1_gmin_m4_prop15682.py`,
`evidence/e1_gmin_m4_prop15682.json`,
`evidence/NOTE_2026-08-28_p31_next_all_finite_endpoint.md`, and
`tests/test_prop15682.py`.

## Proposition 15.683 — the p=41 next all-finite endpoint is impossible

At (p=41,s=34), Proposition 15.681 excludes all positive phase-zero
residues (u_0=2,ldots,7), leaving only (u_0=0,u_1=20). Its type minima
are

\[
14[b=0]+7[b=34],\quad D_0=476,
\qquad
20[b=2]+[b=34],\quad D_1=640,                         \tag{15.683.1}
\]

against pair budget (34\cdot33=1122). Exact completion-bounded enumeration
leaves nine phase-labelled profiles and four global floor-secant shapes:

\[
\begin{array}{c|l}
\text{pair slack}&\{t:\text{number of directions}\},\quad t=(34-b)/2\\ \hline
4&\{17:14,16:20,1:1,0:7\},\\
0&\{17:14,16:20,3:1,0:7\},\\
0&\{17:14,16:20,2:1,1:1,0:6\},\\
0&\{17:14,16:20,1:3,0:5\}.
\end{array}                                           \tag{15.683.2}
\]

There are seven slack-zero profiles, all 34-arcs, and two slack-four
profiles, each with exactly one 3-secant and no other line of occupancy
above two.

We use Segre's tangent envelope in the polynomial form of
Ball--Lavrauw, *Planar arcs*, Theorem 11. If an arc (A) in odd order has
size (q+2-\tau) and (|A|>2\tau+2), there is a nonzero homogeneous dual
polynomial (\Phi) of degree (2\tau) whose restriction to every point-pencil
line (P^*) is the square of the tangent polynomial at (P). In particular,
each tangent dual point is a double zero on (P^*).

The following component observation will be used twice. Let (D^*) be the
dual line parametrising one direction and suppose that direction contains
(b_D) tangents. If (b_D>2\tau), root counting on (D^*) gives
(D^*\mid\Phi). Write (\Phi=D^*\Psi). At each such tangent point, its
point-pencil line meets (D^*) transversely, while the restriction of
(\Phi) has a double zero. Hence (\Psi) vanishes at all (b_D) points. If
(b_D>2\tau-1), root counting again gives

\[
 (D^*)^2\mid\Phi.                                    \tag{15.683.3}
\]

First take a slack-zero row. The fourteen (t=17) directions are perfect
matchings, the twenty (t=16) directions miss two points each, and the
remaining eight exceptional directions contain three secant edges in
total. Every exceptional direction therefore has at least
(34-2\cdot3=28) tangents. Here (\tau=41+2-34=9), so (\deg\Phi=18) and
(34>20). Since (28>18>17), (15.683.3) makes all eight exceptional
direction lines double components. Their degree is sixteen, leaving a
residual conic (Q).

If (e_P) exceptional secant edges meet an arc point (P), then exactly
(8-e_P) exceptional directions are tangent there. Since (P) has nine
tangents in all, it has (1+e_P) tangents in the remaining directions.
The three exceptional edges are distinct and touch at least three points.
At every touched point (Q|_{P^*}) therefore has at least two distinct
double zeros. Since (4>\deg Q=2), (P^*) must divide (Q). This would give
three distinct line components in a conic, impossible.

Now take a slack-four row. Its exceptional directions have floor-secant
counts (1,0,\ldots,0). Choose a point of the unique triple whose deletion
preserves the exceptional floor-secant as a pair: if that floor-secant is
the triple, any deletion works; otherwise its ordinary pair uses at most
one triple point. The remaining 33 points form an arc. Seven exceptional
directions now contain 33 tangents and the eighth contains 31. Here
(\tau=10), (\deg\Phi=20), and (33>22), so their eight squares again divide
(\Phi), leaving a quartic (Q).

The two endpoints of the surviving exceptional pair have seven exceptional
tangents and hence three other tangents. Each of their point-pencil lines
contains three distinct double zeros of (Q); since (6>4), both lines divide
(Q). Removing them leaves a conic (R). Each of the other 31 arc points has
eight exceptional and two other tangents. Those tangent dual points lie on
neither removed point-pencil line, since a tangent cannot contain a second
arc point. Thus (R) has two distinct double zeros on each remaining
point-pencil, forcing that line to divide (R). Three choices already
contradict (\deg R=2).

All nine residue-zero profiles are impossible, so (p=41,s=34) is excluded.
This second all-finite boundary remains open only at (p=17,19,23); later
all-finite sizes, strict infinity-plus-(p), residual (ii), R1, global QVAR,
Type I, and the limit remain open.

External input: S. Ball and M. Lavrauw, *Planar arcs*, J. Combin. Theory
Ser. A **160** (2018), 261--287, Theorem 11,
doi:10.1016/j.jcta.2018.06.015 (arXiv:1705.10940v4).

Evidence: `src/e1_gmin_m4_prop15683.py`,
`evidence/e1_gmin_m4_prop15683.json`,
`evidence/NOTE_2026-08-28_p41_next_all_finite_endpoint.md`, and
`tests/test_prop15683.py`.

## Proposition 15.684 — OPEN_RETRACTED_REDUCTION at the p=23 next endpoint

At (p=23,s=20), the corrected type arithmetic leaves only (u_1=11) in phase
one and

\[
 u_0\in\{0,2,3,4,5,6,8,9\}.                         \tag{15.684.1}
\]

For (u_0\in\{2,3,4,5,6,8\}), a positive residue forces a quotient-zero
direction and a nonzero nonnegative integer-valued quadratic (B) on
(J(23,12)), with scaled mass (c=4p\mathbb E B=2u_0). Proposition 15.681 gives
(c\ge12), excluding (u_0=2,3,4,5). We now exclude (c=12,16).

Put (H=\max B) and choose a maximum point (X). Stabilizer averaging gives
(c\ge4H), while the paired 11-cube through (X) has mean

\[
 \frac{H+c/4}{24}.                                  \tag{15.684.2}
\]

A nonzero nonnegative integral quadratic on a cube has mean at least
(1/4). Equality cannot occur if a value is two. If its range is contained
in (\{0,1,2,3\}) and it takes value three, then over (\mathbb F_2)

\[
 f=(f\bmod2)+2\left({f\choose2}\bmod2\right).
\]

The two nonzero bits have degrees at most two and four, so the elementary
Reed--Muller distance bound gives (\mathbb E f\ge1/4+2/16=3/8). Comparing
these floors with (15.684.2) excludes every height (H\le3) at (c=12,16).

It remains to exclude (c=16,H=4). The exact stabilizer identity at (X) is

\[
 \mathbb E B=\frac{22}{23}q(6)+\frac{B(X)}{23}.
\]

Both sides apart from (q(6)) equal (4/23), so nonnegativity forces (B=0)
on the shell (|X\cap Y|=6). Now (V_2(J(23,12))) has dimension 253, while
its restriction to (J(12,6)\times J(11,6)) has filtered harmonic dimension

\[
 1+11+10+54+44+11\cdot10=230.
\]

The 23-dimensional kernel is exactly ((|X\cap Y|-6)V_1): the displayed
restriction is onto, and multiplication is injective because an affine
factor vanishing on both neighbouring shells (t=5,7) is zero. Hence

\[
 B(Y)=(|X\cap Y|-6)L(Y)                              \tag{15.684.3}
\]

for affine (L). At (X), (L(X)=2/3). If
(b_{ij}=B(X-i+j)), then (L(X-i+j)=b_{ij}/5), and the affine
parallelogram identity at two replacements yields

\[
B(X-i_1-i_2+j_1+j_2)
 =\frac{4\{3(b_{i_1j_1}+b_{i_2j_2})-10\}}{15}.       \tag{15.684.4}
\]

This cannot be an integer, already modulo three. Thus the treated positive
residues (u_0=2,3,4,5,6,8) are impossible. The corrected floor-plus-two
classifier restores (u_0=9), of scaled mass (c=18), together with an explicit
slack-zero profile. Neither the (c\ge12) bound nor the special (c=12,16)
argument excludes this row. The old assertion that only residue zero remains
is therefore false.

Conditionally on (u_0=0), completion-bounded exact enumeration gives 426 phase-zero
rows, 11 phase-one rows, 1,247 compatible phase-labelled profiles, and 485
global shapes. Their pair-slack profile counts are

\[
\begin{array}{c|rrrrrrrrrrrrrrrr}
\text{slack}&0&4&8&12&16&20&24&28&32&36&40&44&48&52&56&60\\ \hline
\text{profiles}&363&264&189&136&94&68&49&35&21&13&7&4&1&1&1&1.
\end{array}                                           \tag{15.684.5}
\]

All 363 slack-zero rows are 20-arcs. Segre's degree-10 tangent envelope
has each direction with at most four secants as a double component. If any
high direction contains an arc edge, the residual curve is forced to
contain more point-pencil lines than its degree. Otherwise three high
directions are undetermined; adjoining pairs of their infinity points gives
22-arcs, and the complete-arc classification forces one conic through three
collinear infinity points. Thus every arc row is impossible.

For a line occupancy (n), define

\[
 \delta(n)=2\left({n\choose2}-\lfloor n/2\rfloor\right).
\]

Deleting (n-2\le\delta(n)/4) points repairs that line, so slack (4r)
permits deleting at most (r) points to obtain an arc. Coolsaet--Sticker's
complete classification of (PG(2,23)) has no complete arc of size 18
through 23 and has a unique 24-arc, a conic. Hence every repaired arc of
size at least 18 lies on a conic. If (1\le h\le4) original points lie off
that conic, off-conic secant counting gives

\[
 \text{slack}\ge4h(7-h)\ge24.                       \tag{15.684.6}
\]

Consequently all 264 slack-four and 189 slack-eight profiles are excluded.
One undetermined direction closes 135 of 136 slack-twelve profiles, and
two close 93 of 94 slack-sixteen profiles. In total, 1,044 of the 1,247
residue-zero profiles are impossible. Exactly 203 arithmetic profiles remain
inside this conditional residue-zero subledger, with slack histogram

\[
\{12:1,16:1,20:68,24:49,28:35,32:21,36:13,40:7,
44:4,48:1,52:1,56:1,60:1\}.                         \tag{15.684.7}
\]

Because the restored (u_0=9) branch is untouched, the claimed whole-endpoint
reduction is retracted: 15.684 has status **OPEN_RETRACTED_REDUCTION**. Its
residue-zero census and the displayed 1,247-to-203 reduction remain valid only
as a conditional subledger. Proposition 15.721 independently supersedes and
excludes this all-finite boundary as a live gate. Strict infinity-plus-(p),
residual (ii), R1, global QVAR, Type I, and the limit remain open.

External inputs: K. Coolsaet and H. Sticker, *A full classification of the
complete k-arcs of PG(2,23) and PG(2,25)*, J. Combin. Des. **17** (2009),
459--477, doi:10.1002/jcd.20211; and Ball--Lavrauw, *Planar arcs*, Theorem
11, doi:10.1016/j.jcta.2018.06.015.

Evidence: `src/e1_gmin_m4_prop15684.py`,
`evidence/e1_gmin_m4_prop15684.json`,
`evidence/e1_gmin_m4_prop15684.historical_retracted.json`,
`evidence/NOTE_2026-08-28_p23_low_mass_conic_reduction.md`, and
`tests/test_prop15684.py`.

## Proposition 15.685 — the unique p=23 slack-twelve profile is impossible

Within the conditional 203-profile residue-zero subledger retained from
Proposition 15.684, the unique row of pair
slack twelve is

\[
8[b=0]+4[b=18],\quad D_0=168,
\qquad
11[b=2]+[b=18],\quad D_1=200.                       \tag{15.685.1}
\]

Let (S) be a hypothetical realization. The repair lemma in 15.684 deletes
at most three points to obtain an arc (A). If at most two points are
deleted, then (|A|\ge18), so the complete-arc classification puts (A) on
a conic. The off-conic estimate (15.684.6) contradicts the positive slack
(12<24).

Thus write (S=A\mathbin{\dot\cup}D), where (|A|=17) and (|D|=3).
If (A) were incomplete, adjoining one point would give an 18-arc, again
placing (A) on a conic and giving the same contradiction. Hence (A) is a
complete 17-arc.

For (x\notin A), let (\mu_A(x)) be the number of secants of (A) through
(x). If one secant contains (r) points of (D), its occupancy in (S) is
(2+r), and

\[
 \delta(2+r)\ge4r\qquad(r=1,2,3),                 \tag{15.685.2}
\]

with respective values (4,8,16). Summing (15.685.2) over the secants of
(A) gives

\[
 12=\operatorname{slack}(S)
 \ge4\sum_{x\in D}\mu_A(x).                      \tag{15.685.3}
\]

Completeness gives (\mu_A(x)\ge1) for every outside point. Equality is
therefore forced in (15.685.3), and all three points of (D) must satisfy
(\mu_A(x)=1).

Coolsaet--Sticker's exhaustive classification contains exactly five
projective classes of complete 17-arcs in (PG(2,23)). Five explicit
representatives in the evidence record are checked directly: each has 17
points, 136 distinct secants, and all 536 outside points are secant-covered.
Their outside-point secant-multiplicity histograms are

\[
\begin{array}{c|rrrrrrrr}
 &\mu=1&2&3&4&5&6&7&8\\ \hline
1&0&2&6&68&172&190&86&12\\
2&0&1&15&59&159&208&86&8\\
3&1&0&6&69&171&196&78&15\\
4&0&0&14&58&170&206&72&16\\
5&0&1&8&63&185&176&91&12.
\end{array}                                         \tag{15.685.4}
\]

The full histogram is invariant under projective equivalence, and the five
rows are distinct. Thus these representatives are inequivalent and,
because the classification has exactly five classes, exhaustive. No class
has more than one point with (\mu=1), contradicting the three required by
(15.685.3).

The slack-twelve profile is impossible. The conditional residue-zero
remainder decreases from 203 to 202 profiles:

\[
\{16:1,20:68,24:49,28:35,32:21,36:13,40:7,
44:4,48:1,52:1,56:1,60:1\}.                        \tag{15.685.5}
\]

This does not treat the restored (u_0=9) branch. Proposition 15.721, rather
than this historical subledger, independently excludes the all-finite
boundary as a live gate. Strict infinity-plus-(p), residual (ii), R1,
global QVAR, Type I, and the limit remain open.

External input: K. Coolsaet and H. Sticker, *A full classification of the
complete k-arcs of PG(2,23) and PG(2,25)*, J. Combin. Des. **17** (2009),
459--477, doi:10.1002/jcd.20211.

Evidence: `src/e1_gmin_m4_prop15685.py`,
`evidence/e1_gmin_m4_prop15685.json`,
`evidence/NOTE_2026-08-28_p23_slack12_complete17_exclusion.md`, and
`tests/test_prop15685.py`.

## Proposition 15.686 — the unique p=23 slack-sixteen profile is impossible

The unique slack-sixteen row left within the conditional residue-zero
subledger of 15.684--15.685 is

\[
7[b=0]+[b=2]+3[b=18]+[b=20],\quad D_0=164,
\qquad
11[b=2]+[b=18],\quad D_1=200.                       \tag{15.686.1}
\]

Its global floor-secant distribution is
(\{t_0:1,t_1:4,t_9:12,t_{10}:7\}), so it has one undetermined
direction. Let (U) be the corresponding infinity point.

Repair deletes at most four points from a hypothetical realization (S).
If at most three are deleted, the repaired arc together with (U) has size
at least 18. It is conic-contained, and (15.684.6) contradicts slack
(16<24).

Thus (S=A\mathbin{\dot\cup}D), where (|A|=16) and (|D|=4). The set

\[
 K=A\cup\{U\}
\]

is a 17-arc. If it were incomplete, extending it to an 18-arc would give
the same conic-core contradiction, so (K) is complete.

For every (d\in D), the line (Ud) contains no second point of (S), since
(U) is undetermined. Therefore no secant of (K) through (d) uses (U), and

\[
 \mu_A(d)=\mu_K(d)\ge1.                            \tag{15.686.2}
\]

Applying the secant-line slack inequality (15.685.2) to the four deleted
points gives

\[
16=\operatorname{slack}(S)
 \ge4\sum_{d\in D}\mu_A(d)\ge16.                \tag{15.686.3}
\]

Equality forces (\mu_K(d)=1) for all four points. But the exhaustive five
complete-17-arc rows (15.685.4) have respectively (0,0,1,0,0) outside
points of multiplicity one. No class supplies four, a contradiction.

The slack-sixteen profile is impossible. Exactly 201 conditional
residue-zero arithmetic profiles remain, all of slack at least 20:

\[
\{20:68,24:49,28:35,32:21,36:13,40:7,
44:4,48:1,52:1,56:1,60:1\}.                        \tag{15.686.4}
\]

The restored (u_0=9) branch is unaffected. Proposition 15.721 independently
excludes the all-finite boundary as a live gate; strict infinity-plus-(p),
residual (ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15686.py`,
`evidence/e1_gmin_m4_prop15686.json`,
`evidence/NOTE_2026-08-28_p23_slack16_complete17_exclusion.md`, and
`tests/test_prop15686.py`.

## Proposition 15.687 — every p=23 slack-twenty profile is impossible

The 68 exact residue-zero profiles of pair slack twenty have
undetermined-direction
histogram

\[
 \{t_0=2:2,\quad t_0=3:36,\quad t_0=4:30\}.       \tag{15.687.1}
\]

Repair deletes at most five points. We first extend the conic-core estimate
from 15.684 through this fifth point. If (h) of the original twenty points
lie off a conic, each off-conic point has at least eleven full conic
secants. The retained (20-h) conic points omit (4+h), so at least (7-h)
secants retain both endpoints. Hence

\[
 \operatorname{slack}(S)\ge4h(7-h).                \tag{15.687.2}
\]

For (1\le h\le5), the right side has values (24,40,48,48,40). Thus a
positive-slack conic-core set cannot have slack twenty.

For the 66 profiles with at least three undetermined directions, choose
three infinity points (U_1,U_2,U_3) and use the overlapping pairs
(\{U_1,U_2\}) and (\{U_1,U_3\}); the three collinear points are never
adjoined simultaneously. If repair deletes at most four points, each pair
gives an arc of size at least 18. Their conics share the repaired arc of
size at least 16, hence coincide, and would contain all three collinear
infinity points.

If all five deletions are needed, each pair gives a 17-arc. Were either
complete, the argument below would force all five deleted points to have
multiplicity one, contrary to (15.685.4). Thus both pair arcs are
incomplete and extend to 18-arcs. Their conics share the repaired 15-arc,
again coincide, and contain (U_1,U_2,U_3), a contradiction.

It remains to consider the two profiles with exactly two undetermined
directions. If repair deletes at most four points, adjoining both infinity
points again gives an arc of size at least 18. In the five-deletion branch,
write (S=A\mathbin{\dot\cup}D), with (|A|=15), (|D|=5), and adjoin both
infinity points to obtain a 17-arc (K). It must be complete, since otherwise
an 18-arc extension gives the conic contradiction.

No secant of (K) through a deleted point can use either infinity point:
both directions are undetermined for (S). Therefore
(\mu_A(d)=\mu_K(d)\ge1) for all (d\in D). The line-slack inequality
then gives

\[
20\ge4\sum_{d\in D}\mu_A(d)\ge20,               \tag{15.687.3}
\]

so all five deleted points have multiplicity one outside (K). The exhaustive
rows (15.685.4) show that a complete 17-arc has at most one such point.
This contradiction excludes the final two profiles.

All 68 slack-twenty profiles in this subledger are impossible. The
conditional residue-zero remainder is 133 profiles, all of slack at least 24:

\[
\{24:49,28:35,32:21,36:13,40:7,
44:4,48:1,52:1,56:1,60:1\}.                        \tag{15.687.4}
\]

The restored (u_0=9) branch is unaffected. Proposition 15.721 independently
excludes the all-finite boundary as a live gate; strict infinity-plus-(p),
residual (ii), R1, global QVAR, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15687.py`,
`evidence/e1_gmin_m4_prop15687.json`,
`evidence/NOTE_2026-08-28_p23_slack20_complete17_exclusion.md`, and
`tests/test_prop15687.py`.

## Proposition 15.688 — sharp integral lift; p=19 reduced to residue zero

Let \(p\ge5\) be odd and let \(B\not\equiv0\) be a nonnegative,
integer-valued quadratic on \(J(p,(p+1)/2)\). Then, sharply,

\[
4p\,\mathbb E B\ge p-3.                         \tag{15.688.1}
\]

Choose \(X\) with \(B(X)=H=\max B\). The paired-cube operator of
15.680--15.681 satisfies

\[
TB(X)=\frac{H+p\mathbb E B}{p+1}.                \tag{15.688.2}
\]

The cube restriction is a nonzero degree-two polynomial, hence has support
density at least \(1/4\). Because it is integer-valued, its multilinear
coefficients are integers and its cube mean lies in \(\frac14\mathbb Z\).
If \(H=1\), (15.688.2) gives (15.688.1). If \(H\ge2\), mean \(1/4\)
would force value one on a support of density exactly \(1/4\), contradicting
the value at \(X\); the mean is therefore at least \(1/2\), so

\[
4p\mathbb E B\ge2(p+1)-4H.                       \tag{15.688.3}
\]

At the same point, the exact stabilizer identities of 15.642 give

\[
4p\mathbb E B\ge
\begin{cases}
4H,&p\equiv3\pmod4,\\
\dfrac{4r}{r+1}H,&p=4r+1.
\end{cases}                                      \tag{15.688.4}
\]

For \(H\ge2\), the maximum of (15.688.3) and (15.688.4) is at least
\(p+1\) in the first branch and \(p-1\) in the second. Thus the global
minimum occurs at \(H=1\), proving (15.688.1). Equality is attained by
\(B=(1-x_i)(1-x_j)\), whose mean is \((p-3)/(4p)\).

At the live \(p=19\) second all-finite boundary \(s=16\), exact pair
arithmetic leaves the unique phase-one row \(u_1=9\) and phase-zero residues

\[
u_0\in\{0,2,3,4,6\}.                              \tag{15.688.5}
\]

Each positive residue forces a quotient-zero direction. Its scaled mean is
\(2u_0\in\{4,6,8,12\}\), below the least positive fibre floor \(20\),
so it has \(b=0\) and is a nonzero quadratic lift. This contradicts
(15.688.1), whose scaled floor is \(16\). The deficit-minimizing pair is

\[
u_0=0,\quad u_1=9,\qquad
5[b=0]+5[b=16],\quad9[b=2]+[b=16].               \tag{15.688.6}
\]

Its putative slack is \(34\not\equiv0\pmod4\), so this minimum pair is not
realizable. It is not, however, the complete residue-zero row.
Completion-bounded exact enumeration gives 60 phase-zero rows and nine
phase-one rows. Pair compatibility and the exact slack congruence leave
143 phase-labelled profiles, 75 global shapes, with histogram
\(\{0:54,4:37,8:25,12:13,16:7,20:4,24:1,28:1,32:1\}\).

Thus 15.688 removes every positive residue but does not reduce residue zero
to a single profile. The \(p=19\) endpoint, residual (ii), R1, global QVAR,
Type I, and the limit remain open.

Evidence for 15.688: `src/e1_gmin_m4_prop15688.py`,
`evidence/e1_gmin_m4_prop15688.json`,
`evidence/NOTE_2026-08-29_sharp_integral_lift_p19.md`, and
`tests/test_prop15688.py`.

## Proposition 15.689 — p=19 low-slack conic reduction

Of the 143 exact residue-zero profiles in (15.688.7), every profile of pair
slack at most twelve is impossible. Exactly fourteen profiles remain, with
slack histogram `{16:7,20:4,24:1,28:1,32:1}`.

The finite-geometric input is the complete-arc spectrum of
`PG(2,19)`: complete arcs have sizes `10,11,12,13,14,20`, and the unique
complete 20-arc is a conic. Consequently every arc of size at least fifteen
is conic-contained.

At slack zero the boundary `S` is a 16-arc. Three undetermined infinity
points give two overlapping 18-arc extensions. Their conics share `S`, so
coincide and contain three collinear points. With one or two undetermined
directions, the containing conic is tangent or secant to the line at
infinity. Every other direction retains at least six affine conic secants,
so has `b<=4`; all 29 exact one-/two-direction profiles have a
non-undetermined `b>=6`.

For positive slack put

\[
\delta(n)=2\left\{\binom n2-\lfloor n/2\rfloor\right\}.
\]

Deleting at most `slack/4` points repairs the boundary to an arc. Slack
four repairs to a conic-contained 15-arc; the deleted off-conic point has
at least four retained conic secants, forcing slack at least sixteen. At
slack eight or twelve, every profile has two undetermined directions.
Adjoining them after repair gives an arc of size at least sixteen or
fifteen and hence a conic. If `j` deleted boundary points lie off it,
retained-secant counting and `delta(2+r)>=4r` give

\[
\operatorname{slack}(S)\ge4j(5-j).
\]

For `j=1,2,3` this is `16,24,24`, contradicting slack eight or twelve.
Thus `54+37+25+13=129` profiles are excluded. This is a strict reduction,
not closure of the endpoint.

External input: G. Faina, S. Marcugini, A. Milani, and F. Pambianco,
*The spectrum of values k for complete k-arcs in PG(2,q) for q<=23*,
Ars Combinatoria **47** (1997), 3--11; independently tabulated in H.
Sticker's complete-arc classification thesis.

Evidence: `src/e1_gmin_m4_prop15689.py`,
`evidence/e1_gmin_m4_prop15689.json`, and `tests/test_prop15689.py`.

## Proposition 15.690 — exact dilation energy and its method barrier

Let `H=(F_q^*)^2`, `K=H/{±1}`, and let `S_K` be the full-Max+ dilation
energy in the cold strategy note. Exact square-torus character orthogonality
and the affine autocorrelation identity give

\[
\bar L=\frac{q(q-1)(q+5)}{16},\qquad
S_K=\frac{q-1}{2}\frac Vn
=12(q-1)\frac{\|\delta\|^2}{n}.                 \tag{15.690.1}
\]

Consequently

\[
S_K\le q-1\quad\Longleftrightarrow\quad
V\le2n\quad\Longleftrightarrow\quad
\|\delta\|^2\le n/12.                           \tag{15.690.2}
\]

Thus the dilation inequality is exactly strong R1 in new coordinates, not
an auxiliary estimate implied by representation theory. Moreover,
positivity, equivariance, trace, fixed norms, Bochner positivity, and
additive autocorrelation do not imply it. Explicit abstract invariant
spectra and explicit PSD autocorrelations satisfy those relaxations while
violating (15.690.2), the latter by a factor `Theta(q^2)`. These are not
actual full-Max+ counterexamples. Any successful proof must use the Boolean
rank-one identity and exact cancellation among all Max+ orbit types.

Evidence: `evidence/NOTE_2026-08-29_dilation_energy_normalization_and_no_go.md`.

## Proposition 15.691 — the c=2 signed-Eulerian target is false

Let `N=binom(n,2)` and

\[
Y_a(\beta)=\mathbb E_x\cosh(\beta Q_a(x))
=(\cosh\beta)^N P_a(\tanh\beta).
\]

For `0<theta<1`, averaging over independent edge signs gives

\[
\mathbb E_aY_a(\beta)^\theta
\le2^{n(1-\theta)+1}\cosh(\theta\beta)^N.
\]

Hence some deterministic signing satisfies, for `beta=c/sqrt(n)`,

\[
\log P_a(\tanh(c/\sqrt n))
\le-\left(\frac c2-\sqrt{\log2}\right)^2n+o(n). \tag{15.691.1}
\]

At `c=2`, this disproves the proposed uniform lower bound `log P_a>=-o(n)`
by a linear margin. The corrected sufficient target is

\[
\inf_a\log P_a(\tanh(c/\sqrt n))
\ge\left(\frac c2-\frac{c^2}{4}\right)n-o(n).    \tag{15.691.2}
\]

The fractional-moment construction rules out (15.691.2) for
`c<2.0843108...`; the historical candidate left by this proposition was

\[
\inf_a\log P_a(\tanh(3/\sqrt n))\ge-3n/4-o(n).  \tag{15.691.3}
\]

Proposition 6.9 supersedes that status: (15.691.2) is false for every fixed
`c>0`, including (15.691.3), on symmetric conference signings. Only a
growing-temperature version remains logically possible. No top-level gate is
closed.

Evidence: `evidence/NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md`.

## Proposition 15.692 — binary affine-Radon isomorphism at the p=19 endpoint

Let (A) be the binary line-point incidence matrix of
(operatorname{AG}(2,p)), with rows grouped by parallel class. For odd
(p),

\[
A^{\mathsf T}A=I+J\pmod2.                         \tag{15.692.1}
\]

Thus (A^{\mathsf T}A=I) on the even-weight point space. The image has
even weight in each directional block, and the source and target dimensions
both equal

\[
p^2-1=(p+1)(p-1).
\]

Consequently the restricted Radon map is an isomorphism, with inverse

\[
x=A^{\mathsf T}r.                                \tag{15.692.2}
\]

For each of the fourteen profiles left by 15.689, the exact finite problem
is therefore: choose line-parity blocks having the prescribed phase weights
and impose (operatorname{wt}(A^{\mathsf T}r)=16). There are no further
linear compatibility equations. All fourteen pass the resulting mod-four
weight condition. Their fixed first two stripe-count moments also admit
explicit nonnegative distributions supported entirely on multiplicities
({4,6,8}), so pairwise independence and second moments cannot force a
positive odd-parity density.

This is an exact reduction and method barrier, not closure of the endpoint.
Residual (ii), R1, Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15692.py`,
`evidence/e1_gmin_m4_prop15692.json`, and `tests/test_prop15692.py`.

## Proposition 15.693 — the p=19 slack-sixteen block is impossible

Each of the seven slack-sixteen profiles left by Proposition 15.689 has
three or four undetermined directions. The repair lemma deletes at most
four points. With at most three deletions, adjoining two undetermined
infinity points gives an arc of size at least fifteen. The complete-arc
spectrum forces this arc into the conic; the third undetermined infinity
point then has only the line-at-infinity secant although retained-conic
counting forces further secants, a contradiction.

It remains to consider four deletions. Write \(S=A\cup D\), where \(A\) is
a 12-arc and \(|D|=4\), and adjoin two undetermined infinity points:

\[
K=A\cup\{U_1,U_2\},\qquad |K|=14.                \tag{15.693.1}
\]

If \(K\) is incomplete, the gap in the complete-arc spectrum extends it to
the 20-point conic. A third undetermined infinity point has exactly one
\(K\)-secant, the line at infinity. An off-conic point has at least nine
conic secants, while the six omitted conic points destroy at most six, so
at least three \(K\)-secants remain. Hence \(K\) is complete.

For every \(x\in D\), undeterminedness prevents a \(K\)-secant through
\(x\) from using \(U_1\) or \(U_2\). Such a secant is already an
\(A\)-secant and charges four units of slack. Completeness and total slack
sixteen therefore force every point of \(D\) to have secant index one.
Each unused undetermined infinity point also has index one, its unique
\(K\)-secant being the line at infinity. Thus

\[
c_1(K)\ge4+(t-2)\ge5.                            \tag{15.693.2}
\]

Al-Zangana's exhaustive classification of all 83 projective 14-arc classes
in `PG(2,19)` gives \(c_1\le4\) for every class. This contradiction removes
all seven slack-sixteen profiles. The exact remainder is

\[
\{20:4,24:1,28:1,32:1\},
\]

seven profiles in total. The same classification forces every slack-twenty
row to use exactly five repair deletions. The endpoint, residual (ii), R1,
Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15693.py`,
`evidence/e1_gmin_m4_prop15693.json`, and `tests/test_prop15693.py`.

## Proposition 15.694 — exact equality normal form at p=19 slack twenty

Proposition 15.693 forces every slack-twenty witness to use exactly five
repair deletions. Write \(S=A\mathbin{\dot\cup}D\), with \(|A|=11\) and
\(|D|=5\), where \(A\) is an arc. If some \(x\in D\) lay on no secant of
\(A\), then \(A\cup\{x\}\) would be a four-deletion repair, contrary to
15.693. Therefore \(\mu_A(x)\ge1\) for every deleted point, and

\[
20=\operatorname{slack}(S)
 \ge4\sum_{x\in D}\mu_A(x)\ge20.                \tag{15.694.1}
\]

All inequalities are equalities, so every deleted point has
\(\mu_A(x)=1\). On a line \(\ell\), put
\(a=|A\cap\ell|\) and \(d=|D\cap\ell|\). Since \(a\le2\), equality between
the exact line slack and its charged contribution leaves precisely

\[
(a,d)\in\{(0,0),(0,1),(0,2),(1,0),(1,1),
          (2,0),(2,1),(2,2)\}.                 \tag{15.694.2}
\]

Hence \(D\) is also an affine arc, no boundary line has more than four
points, and a line through two deleted points contains zero or two core
points. The five charged incidences have exactly three possible patterns:
five 3-lines; one 4-line and three 3-lines; or two 4-lines and one 3-line.

Adjoin any two of the \(t\) undetermined infinity points to \(A\), obtaining
a 13-arc \(K\). Its five deleted points have secant index one, as do the
other \(t-2\) undetermined infinity points. Thus

\[
c_1(K)\ge5+(t-2)=7\quad(t=4),\qquad
c_1(K)\ge8\quad(t=5).                           \tag{15.694.3}
\]

Al-Zangana's exhaustive 2,733-class 13-arc classification gives
\(c_1\le9\). This is a strict class and search-space filter, but not a
contradiction. All four slack-twenty profiles and the seven-profile p=19
remainder stay open. Bounded exact solver runs returned UNKNOWN and are not
used as evidence.

Evidence: `src/e1_gmin_m4_prop15694.py`,
`evidence/e1_gmin_m4_prop15694.json`, and `tests/test_prop15694.py`.

## Proposition 15.718 — global positive-p7 z7 reduction and Johnson semigroup

After Proposition 15.717, the positive `p=7` infinity-plus-seven remainder
consists of 56 actual line boundaries in two orbits.  Their exact-mean leaves
and the two pointed-star branches give 4,320 pointed cases.  The exact parent
affine-hull sieve over characteristics three and seven gives

\[
4,320\longrightarrow1,296,
\]

with 3,024 rigorous projection rejections.  Exact affine symmetry partitions
the 1,296 survivors into 324 complete four-case classes.  Testing one
representative per class with a joint mod-3/mod-7 catalog join, in which the
same catalog-row index supplies both characteristics, gives

| representative decision | representatives | four-case transfer |
|---|---:|---:|
| rigorous empty global join | 87 | 348 |
| necessary-only join survivor | 159 | 636 |
| explicit side-budget skip | 78 | 312 |
| **total** | **324** | **1,296** |

The empty joins are rigorous: exact byte equality, rather than a digest, is
used for membership, and the audited affine symmetry preserves the joint
same-row decision.  The 159 survivors pass only a relaxation in which high
catalogs are replaced by exact affine hulls; the other 78 representatives are
explicitly unprocessed at the certified state cap.  Most importantly, the
348 transferred rejections count pointed leaf cases.  They do not subtract
348, or any number, from the 56 source line boundaries.

Independently, let

\[
\mathcal S=\{L\in\mathbb Z_{\ge0}^{35}:KL=0\}
\]

be the Johnson incidence semigroup, graded by
\(g=\sum_XL_X/5\).  Its complete Normaliz Hilbert basis has 896 distinct
binary rows and no primitive generator above grade three:

| grade | new Hilbert generators |
|---:|---:|
| 1 | 56 |
| 2 | 168 |
| 3 | 672 |

The 56 grade-one rows are exactly the `S56` catalog.  The complete
`M1764` layer is the disjoint union of 1,596 sums of two grade-one rows and
the 168 primitive grade-two rows.  The exact uncapped semigroup layers are

| grade | rows |
|---:|---:|
| 0 | 1 |
| 1 | 56 |
| 2 | 1,764 |
| 3 | 37,856 |
| 4 | 575,407 |
| 5 | 6,496,938 |
| 6 | 57,232,105 |
| 7 | 410,200,367 |
| 8 | 2,474,264,653 |

Every generator is binary.  Thus every semigroup element through grade six
automatically obeys the catalog coordinate cap \(L_X\le6\).  At grade eight
the displayed count is an uncapped outer count, because repeated generators
can violate that cap.  This proves the high-catalog generation structure, not
an exact edge lift or case feasibility.  All 56 actual `z=7` boundaries, both
orbits and projected profiles, the positive endpoint, and the main theorem
remain open.

Evidence: `src/e1_gmin_m4_prop15718.py`,
`evidence/p7_infinity7_positive_z7_global_semigroup_summary.json`, and
`evidence/NOTE_2026-08-30_p7_infinity7_positive_z7_global_semigroup.md`.

## Proposition 15.719 — finite projected-semigroup stabilization

Fix a finite projected group \(G\).  If \(A_i\subseteq G\) is the image of
the Hilbert generators of grade \(i\), define

\[
T_0=\{0\},\qquad
T_g=\bigcup_{i=1}^3(T_{g-i}+A_i).
\]

If four consecutive supports agree,

\[
T_r=T_{r+1}=T_{r+2}=T_{r+3}=S,
\]

then \(S+A_i\subseteq S\) for every generator family.  Translation is
injective on finite \(G\), so every containment is equality.  It follows that
\(S\) is the subgroup generated by all projected Hilbert generators and
\(T_g=S\) for every \(g\ge r\).

The exact `k=3` and `k=4` certificates have equal raw and anchor-relative
support hashes on grades 3, 4, 5, and 6 in all eight directions.  Their
stabilized subgroup orders are

| projection | directions 0--4 | direction 5 | directions 6--7 |
|---|---:|---:|---:|
| `F_3^6 x F_7^3` | 147 | 3 | 3 |
| `F_3^6 x F_7^4` | 1,029 | 21 | 3 |

Because the Hilbert generators are binary, this identifies the capped
catalog projections exactly through grade six.  Grade eight is still only an
outer support: absence can reject, but presence does not prove the cap
\(L_X\le6\), an integral catalog tuple, or binary edge feasibility.

Two completed campaigns are recorded only as strategy context.  First, the
four `H0_S0_M7` calibration representatives have nonseparating exact torsion
projections for every tested width `k=2,...,6`; every completed target is
present.  Second, bounded case joins over all 51 grade-three-only
representatives found no rejection under eight selected `k=3` projections or
under the completed one-projection `k=4` and `k=5` prefixes.  The even
26-case `k=6` shard likewise found no rejection; the odd shard was cancelled,
so no complete `k=6` statement is made.  A final direct model for one fixed
case imposed all 27 coordinates of the full `F_3^6 x F_7^21` quotient and
returned `UNKNOWN` after 300 seconds.  That status has no mathematical force,
and the semigroup/quotient route is terminated.  These target-presence runs
are necessary-condition diagnostics, not Proposition 15.719, not feasibility
certificates, and not a boundary-count reduction.

Consequently all 56 actual positive `z=7` line boundaries in two orbits, the
positive endpoint, residual (ii), R1, the non-Walsh multi-level remainder,
and the quadratic-minmax-limit theorem remain open.

Evidence: `src/e1_gmin_m4_prop15719.py`,
`evidence/p7_infinity7_positive_z7_projected_stabilization_summary.json`, and
`evidence/NOTE_2026-08-30_p7_infinity7_positive_z7_global_semigroup.md`.

## Proposition 15.720 — degree-congruence obstruction for bi-tight levels 2 and 3

Let \(q=p^2\), \(n=q+1\), and let \(C\) be the symmetric conference
matrix with \(C^2=qI\). Let \(H\) have indicator \(h\), size
\(|H|=sp\), degrees \(d_i\), and density

\[
a=\frac{|H|}{\binom n2}=\frac{2s}{np}.
\]

Suppose \(H\) is bi-tight of level \(s\). Since the full-edge score is
\(pn/2\) on \(\operatorname{Max}_+\) and \(-pn/2\) on
\(\operatorname{Max}_-\), the centered indicator
\(\kappa=h-a\mathbf1\) has zero score on both sets. Thus
\(\kappa\in\ker(G_++G_-)\). Propositions 15.272 and 15.207 give

\[
\ker(G_++G_-)=\mathrm{scheme}\oplus\mathrm{cross}.
\]

Put \(A=C\odot\kappa\) and \(B=C\odot h=A+aC\). Membership in
`scheme+cross` means

\[
A=D_fC+CD_f+X,\qquad \sum_i f_i=0,\qquad CX+XC=0.
\]

The matrix \(XC\) is skew-symmetric, hence has zero diagonal. Comparing the
diagonal of \(AC\) gives

\[
f_i=\frac{d_i-aq}{q-1}.
\]

After absorbing \(aC\), write \(B=D_gC+CD_g+X\), where

\[
g_i=f_i+\frac a2=\frac{d_i-s/p}{q-1}.
\]

The commuting projection

\[
\operatorname{Comm}(M)=\frac12\left(M+\frac1q CMC\right)
\]

kills \(X\) and fixes \(D_gC+CD_g\). Therefore, for \(i\ne j\),

\[
h_{ij}+\frac{C_{ij}(CBC)_{ij}}q=2(g_i+g_j).
\]

The quantity \(C_{ij}(CBC)_{ij}\) is an integer. Clearing denominators and
reducing modulo \(q-1\) therefore gives

\[
q-1\mid2(d_i+d_j)-4ps,
\]

or equivalently

\[
\boxed{d_i+d_j\equiv2ps\pmod{(p^2-1)/2}}. \tag{15.720.1}
\]

Subtracting two equations with one common vertex shows that all degrees have
one common residue modulo \(M=(p^2-1)/2\).

For \(s=2\), \(M>2p\) for every \(p\ge5\), so all degrees are equal. The
handshake identity would then give

\[
d=\frac{4p}{p^2+1}\in(0,1),
\]

which is impossible for an integer degree. For \(s=3\) and \(p\ge7\),
\(M>3p\), and the same argument gives the impossible degree
\(6p/(p^2+1)\in(0,1)\). At \(p=5\), \(M=12\), \(n=26\), and
\(\sum_i d_i=30\). Writing \(d_i=r+12m_i\) forces \(r\in\{0,1\}\),
but neither \(30-26\cdot0\) nor \(30-26\cdot1\) is divisible by \(12\).
Hence

\[
\mathscr C_{2p}\cap(\mathrm{scheme}\oplus\mathrm{cross})=\varnothing,
\qquad
\mathscr C_{3p}\cap(\mathrm{scheme}\oplus\mathrm{cross})=\varnothing
\]

for every prime \(p\ge5\). The required level-2 and level-3 bi-tight
alternatives are therefore empty.

The same argument also excludes bi-tight level 4. For \(p\ge11\),
\(M>4p\), so regularity would force the impossible degree
\(8p/(p^2+1)\). At \(p=5\), total degree \(40\) leaves residues \(4\) and
\(2\) modulo \(12\) for \(r=0,1\); at \(p=7\), total degree \(56\) leaves
residues \(8\) and \(6\) modulo \(24\). This is only a bi-tight corollary:
it does not exclude a generic one-sided Max+- or Max−-tight level-4 cover.

Thus Proposition 15.720 repairs the Type-I \(2p\) and deep \(3p\)
no-descent branches without the spectral floor, QVAR, or R1. At this stage the
one-sided level-4 branch remained inside residual (ii), and the two live E(1)
gates were the multi-level Type-I bad case and non-Walsh residual (ii).
Proposition 15.750 later closes the former, leaving residual (ii) as the sole
false entry in that historical four-unit ledger. Proposition 15.764 later
exposes the separate minimal-four-gap implication bridge. \(\square\)

Evidence: `src/e1_gmin_m4_prop15720.py`,
`tests/test_prop15720.py`, and
`evidence/NOTE_2026-08-30_bitight_degree_congruence.md`.

## Proposition 15.721 — signed boundary transport collapses the all-finite ladder

Let \(R_H\) be the symmetric \(\{\pm1\}\)-mask which is \(-1\) on the
edges of a residual flip set \(H\), so that the switched matrix is
\(C\circ R_H\). Proposition 15.267 gives, for every
\(g\in\operatorname{PSL}(2,p^2)\), a permutation matrix \(P\) and a
diagonal sign matrix \(D_g\) satisfying

\[
 D_gP^TCPD_g=C.                                      \tag{15.721.1}
\]

Signed conjugation therefore acts on the relative flip mask by permutation
alone:

\[
\begin{aligned}
D_gP^T(C\circ R_H)PD_g
  &=(D_gP^TCPD_g)\circ(P^TR_HP)\\
  &=C\circ(P^TR_HP).                                \tag{15.721.2}
\end{aligned}
\]

In particular, it preserves \(|H|\), bijects both Boolean eigenshells,
preserves the two margin-three separation inequalities, and sends the
odd-degree boundary \(D=\partial H\) to its image under the vertex
permutation. The Paley edge-product phase may change, but no exact formula
for that change is needed below because every imported exclusion holds in
both phases.

Choose any \(v\in D\). The fractional linear map

\[
 g_v(z)={1\over z-v},\qquad
 \begin{pmatrix}0&1\\1&-v\end{pmatrix},             \tag{15.721.3}
\]

sends \(v\) to infinity and has determinant \(-1\). Since
\(q=p^2\equiv1\pmod4\),

\[
 \chi_q(-1)=(-1)^{(q-1)/2}=1,                       \tag{15.721.4}
\]

so the projective class of (15.721.3) lies in
\(\operatorname{PSL}(2,q)\). Thus every nonempty residual boundary may be
normalized to contain infinity before any affine direction profile is
analyzed.

Now assume \(p\ge17\), and put \(d=|D|\). The handshake lemma makes \(d\)
even. If

\[
 6\le d\le p-3,
\]

transport one point of \(D\) to infinity. The transformed boundary consists
of infinity and the odd number \(s=d-1\) of finite points, where

\[
 5\le s\le p-4.
\]

This is exactly the complete infinity-present range excluded in both phases
by Proposition 15.669. If \(d=p-1\), the same transport gives infinity plus
\(p-2\) finite points, the complete two-phase shell excluded by Proposition
15.674. The cases \(d=0,2,4\) were already excluded by Propositions 15.632,
15.643/15.647, and 15.652, respectively. Consequently

\[
 \boxed{
 \text{for every prime }p\ge17,
 \quad |\partial H|\le p-1\quad\Longrightarrow\quad H
 \text{ is not residual-compatible}.}              \tag{15.721.5}
\]

Because boundary size is even, the first size not excluded by (15.721.5) is

\[
 \boxed{|D|=p+1.}                                    \tag{15.721.6}
\]

This is a lower floor, not an existence assertion. After normalization it is
the infinity-plus-\(p\) shell. Proposition 15.676 excludes its pair-deficit
equality branch; strict pair deficit remains open.

The boundary-close role of the all-finite ladder in Propositions
15.675--15.712 is therefore redundant. More precisely, the first old shell
is at most \(3(p-1)/4+2\le p-3\) for \(p\ge17\). The second is at most
\(3(p-1)/4+4\le p-3\) for \(p\ge29\); the exceptional pairs are
\((p,d)=(17,14),(17,16),(19,14),(19,16),(23,18),(23,20)\), all covered by
(15.721.5). It does not restore the retracted conclusions of 15.678 or
15.684; only their explicitly retained sublemmas and conditional subledgers
survive. The corrected 15.700--15.712 replay is

\[
2503\to2219\to1744\to1481\to1368\to1228\to1215\to1213
\to1020\to869\to321\to19\to14\to0.
\]

Here 15.705 is partial: it excludes only thirteen historical Orbiter targets
and leaves 74 slack-sixteen rows, all removed later by 15.709. Proposition
15.676 remains load-bearing at the new first shell, and Propositions
15.690--15.691 remain independent optional
no-go results. Small-prime remainders, strict deficit at \(|D|=p+1\),
residual (ii), multi-level Type I, and the limit remain open. \(\square\)

Evidence: `src/e1_gmin_m4_prop15721.py`,
`tests/test_prop15721.py`, and
`evidence/e1_gmin_m4_prop15721.json`.

## Proposition 15.722 — exact phase cocycle and multi-chart p+1 reduction

Let (D=\partial H), (|D|=P=p+1), and let
(g(z)=(az+b)/(cz+d)in\operatorname{PSL}(2,p^2)). The signed Paley
multiplier is

Away from the pole, the multiplier is

\[
 \delta_g(x)=\chi(cx+d).
\]

If \(c\ne0\), its value at the pole and at infinity is \(\chi(c)\).
If \(c=0\), there is no pole and its value on every finite point and at
infinity is \(\chi(d)=\chi(a)\); the last equality follows because
\(ad\) is a square.  This affine case split is necessary: substituting
\(c=0\) into \(\chi(c)\) would give zero rather than a switching sign.
Every transported edge receives its two endpoint multipliers. Even vertex
degrees cancel in their product and the odd degrees leave exactly (D), so

\[
 \boxed{c_{gH}=c_H\prod_{x\in D}\delta_g(x).}       \tag{15.722.1}
\]

If (D) is all finite and (f(X)=\prod_{x\in D}(X-x)), inversion about
(r) gives

\[
 c_r=c_H\chi(f'(r))\quad(r\in D),\qquad
 c_r=c_H\chi(f(r))\quad(r\notin D).                \tag{15.722.2}
\]

The Vandermonde identity and (chi_{p^2}(-1)=1) imply

\[
 \prod_{r\in D}\chi(f'(r))=1.                     \tag{15.722.3}
\]

Thus the boundary-chart phases are coupled and have an even number of
negative transports.

Now send an outside point to infinity and write (b_d) for the odd-fibre
profile of the resulting (P) affine points. If a line has occupancy (n),
the exact normalized pair slack is

\[
 R={\sum_d b_d-P\over4}
 =\sum_{n=2r}r(r-1)+\sum_{n=2r+1}r^2.              \tag{15.722.4}
\]

Hence (R=0) is a (P)-arc. The case (R=1) would have exactly one
trisecant and no other line of occupancy at least three. Delete one point
of that trisecant. The remaining (p)-arc lies on a conic. The deleted
point is off the conic and has at least ((p-1)/2) conic secants; the one
missing conic point destroys at most one. For (p\ge17), more than one
surviving trisecant remains, a contradiction. Thus

\[
 \boxed{R\ne1.}                                    \tag{15.722.5}
\]

Here the finite-geometry input is Segre's precise odd-order `q`-arc theorem:
every `q`-arc in the Desarguesian plane `PG(2,q)`, `q` odd and `q>=5`, is
contained in a nonsingular conic.  This is stronger than, and must not be
replaced by, the statement only about `(q+1)`-arcs.

The same incidence identity also excludes (R=2) and (R=3).  For (R=2),
the positive contributions can only be (2) or (1+1), hence the rich lines
are one 4-secant or two 3-secants.  If the latter share a boundary point,
delete it.  The resulting (p)-arc is a conic minus one point, while the
deleted point is off that conic; it retains at least
((p-1)/2-1=(p-3)/2\ge7) conic secants, contradicting total slack two.  In
the other cases delete two points, leaving a `(p-1)`-arc.  Ball--Lavrauw's
classification of complete `(q-1)`-arcs says that for prime (p\ge17) it
extends to a (p)-arc and hence lies on a conic.  Each deleted point is
off-conic because its rich line retains two conic points, and it retains at
least ((p-1)/2-2=(p-5)/2\ge6) conic secants.  This is again impossible.

For (R=3), the only patterns are one 4-secant plus one 3-secant, or three
3-secants.  Choose two private points from the 4-secant and one private
point from the 3-secant, or one private point from each 3-secant.  The
remaining `(p-2)` points form an arc.  The complete `(q-2)`-arc
classification (whose only exceptions have (q=8,9,11)), followed by the
`(q-1)` result and Segre, puts it on a conic.  Every deleted point is
off-conic and retains at least
((p-1)/2-3=(p-7)/2\ge5) conic secants, more rich lines than (R=3) permits.
There is also a prime-dependent extension.  From every rich line of
occupancy (n) delete (n-2) points, take the union, and make the deletion set
(T) inclusion-minimal.  Then (A=D\setminus T) is an arc and, writing
(t=|T|),

\[
 1\le t\le\sum_\ell(n_\ell-2)\le R,
\]

because for (n=2r,2r+1) the two gaps are respectively
((r-1)(r-2)) and ((r-1)^2).  Minimality says that every (z\in T) lies on a
line retaining exactly two points of (A).  Ball--Lavrauw's prime-field conic
threshold puts (A) on a conic whenever

\[
 R\le r_p:=\max\{r\in\mathbb Z_{\ge0}:(2r+5)^2\le4p\}
       =\left\lfloor\sqrt p-\frac52\right\rfloor .
\]

Indeed, then (|A|=p+1-t\ge p-\sqrt p+7/2).  Every deleted (z) is off the
conic, and after its (t) missing points it retains at least
((p-1)/2-t) conic secants.  But (r_p<(p-1)/4), since
(p-4\sqrt p+9=(\sqrt p-2)^2+5>0).  Hence these secants number more than
(R), although each contributes at least one to (R), a contradiction.
Combining this with the uniform near-complete-arc cases gives

\[
 \boxed{1\le R\le \max\{3,r_p\}\ \hbox{ is impossible};\quad
 R>0\Longrightarrow R\ge\max\{4,r_p+1\}.}             \tag{15.722.5a}
\]

At (R=0), Segre puts the points on a conic with external infinity line,
whose profile is (m[b=0]+m[b=2]). In phase one the two floors are
(2p) and (p-1). If (t) phase-one directions have (b=0), their type
floor is

\[
 t(2p)+(m-t)(p-1)=mP+(t-1)P,
\]

so (t\le1). The conic direction partition and the Paley norm-character
partition therefore disagree at at most two projective directions. Their
character correlation is at least (P-4=p-3). If the two anisotropic
binary quadratic forms `Q,N` were not proportional, a shared root over the
algebraic closure would force the shared Frobenius-conjugate root and hence
proportionality. Thus `QN` is squarefree, (Y^2=QN) is a smooth genus-one
curve, and its projective character sum has absolute value at most
(2\sqrt p), contrary to (p-3>2\sqrt p). Hence the conic is a Miquelian
circle and its direction types align exactly.

A one-point replacement ((C\setminus\{r\})\cup\{z\}), (z\notin C), is
also impossible: send (z) to infinity. Möbius transport leaves (p)
finite points of an affine Miquelian circle, a pair-equality branch excluded
by Proposition 15.676.

Finally normalize a full circle as

\[
 D=\{\infty\}\cup(a+b\mathbb F_p),\qquad
 \epsilon=\chi_{p^2}(b),\qquad m={p+1\over2}.
\]

For an outside point (w), put (u=(w-a)/b). The boundary multiplier is

\[
 \chi\!\left(\prod_{t\in\mathbb F_p}(a+bt-w)\right)
 =\chi(b^p(u-u^p))=(-1)^m\epsilon,                 \tag{15.722.6}
\]

because every nonzero trace-zero element has character ((-1)^m). The
transformed affine circle has tangent directions proportional to
(-b^{-1}/(t-u)^2), so its (b=2) directions have type (\epsilon).
Exact outside-chart alignment forces

\[
 \boxed{c_H=(-1)^m.}                               \tag{15.722.7}
\]

Every circle-point chart therefore has common phase (m\bmod2), profile
one (b=1) plus (p) copies of (b=p), special floor
(P-2(m\bmod2)), and zero transverse floors. Proposition 15.724 closes
this full-circle branch.

Evidence: `src/e1_gmin_m4_prop15722.py`,
`tests/test_prop15722.py`, and
`evidence/e1_gmin_m4_prop15722.json`.
The near-complete-arc inputs are S. Ball and M. Lavrauw, *Planar arcs*,
J. Combin. Theory Ser. A **160** (2018), 261--287, published Corollaries
10--11 (Corollaries 8--9 in the arXiv v4 numbering). The prime-field conic
threshold is published Theorem 5 (Theorem 3 in arXiv v4),
doi:10.1016/j.jcta.2018.06.015.

## Proposition 15.727 — endpoint rigidity and the first four prime closes

Continue with the first integer not excluded by Proposition 15.726:

\[
 R=\left\lfloor{p-1\over3}\right\rfloor,
 \qquad p=3R+c,\quad c\in\{1,2\}.                \tag{15.727.1}
\]

Choose (T) of **minimum cardinality** subject to (A=D\setminus T) being an
arc, and put (t=|T|).  This is stronger than merely choosing (T)
inclusion-minimal and will be used below.  The usual rich-line deletion
still gives (1\le t\le R).

Suppose first that (t<R).  At the worst value (t=R-1),

\[
 |A|-(2(t+1)+2)=p-3t-3=p-3R=c\ge1.              \tag{15.727.2}
\]

Thus the same Ball--Lavrauw tangent envelope used in Proposition 15.726
applies for every (1\le t\le R-1).  Equations (15.726.4) and (15.726.8)
give

\[
 I:=\sum_{z\in T}s_A(z)\le R,
 \qquad I\ge F(t):={t(p-1-3t)\over2}.            \tag{15.727.3}
\]

The quadratic (F) is concave.  Its two endpoint margins on this shorter
interval are

\[
\begin{array}{c|cc}
 &F(1)-R&F(R-1)-R\\ \hline
 p=3R+1&(R-3)/2&(R-3)/2\\
 p=3R+2&(R-2)/2&R-2.
\end{array}                                      \tag{15.727.4}
\]

They are positive for every prime (p\ge17).  Hence (t<R) contradicts
(15.727.3), and

\[
 \boxed{t=R.}                                    \tag{15.727.5}
\]

Minimum cardinality implies inclusion-minimality, so every integer
(s_A(z)\ge1).  Combining (15.727.3) and (15.727.5) forces equality
throughout:

\[
 \boxed{I=R,\qquad s_A(z)=1\quad(z\in T).}       \tag{15.727.6}
\]

This also makes the linewise comparison rigid.  On a line containing
(a\le2) points of (A) and (u) points of (T), the contribution to (I) is
(u) when (a=2), and zero otherwise.  The slack contribution is
(h(a+u)).  Equality in the global sum says that every line has equality
locally.  If (a\le1), this forces (a+u\le2).  If (a=2), the formulas in
Proposition 15.726 show that

\[
 h(2+u)=u\quad\Longleftrightarrow\quad u=0,1,2. \tag{15.727.7}
\]

Consequently every rich line of (D) is either a trisecant with composition
((a,u)=(2,1)) or a 4-secant with composition ((2,2)).

These rich lines are pairwise disjoint as subsets of (D).  Indeed, if two
shared a point (v\in D), choose (v) in the deletion demand (n_\ell-2) for
both lines, and choose arbitrary required points on every other rich line.
Deleting the union makes every original rich line have occupancy at most
two, and hence leaves an arc.  Its size is at most

\[
 \sum_{\ell\ {\rm rich}}(n_\ell-2)-1=R-1,
\]

contrary to (15.727.5).  Thus, if (x,y) count the trisecants and
4-secants,

\[
 \boxed{x+2y=R,\qquad\hbox{all rich lines are }D\hbox{-disjoint}.}
                                                               \tag{15.727.8}
\]

There are (c+1+2y) points outside the rich blocks.  The complete projective
line census is

\[
\begin{aligned}
 N_4&=y,&N_3&=R-2y,&N_2&={p(p+1)\over2}-3R,\\
 N_1&=p+1+3R+2y,&N_0&={p(p-1)\over2}-R-y.        \tag{15.727.9}
\end{aligned}
\]

Removing those (c+1+2y) singleton points and one point from every
4-secant leaves a regular trisecant core of size (3(R-y)): every point lies
on one trisecant and has exactly (c+3+3y) tangents.  This is a genuine
all-prime reduction, not an endpoint exclusion by itself.

For the first four primes, however, (15.727.6) contradicts the already
audited exhaustive arc classifications.  In the standard notation

\[
 c_1(A)=|\{z\notin A:s_A(z)=1\}|,
\]

equation (15.727.6) requires (c_1(A)\ge R).  The four comparisons are

\[
\begin{array}{c|c|c|c}
 p&|A|&\hbox{required }c_1(A)&\hbox{classified maximum}\\ \hline
 17&13&5&4\\
 19&14&6&4\\
 23&17&7&1\\
 29&21&9&0.
\end{array}                                      \tag{15.727.10}
\]

Here is the complete branch check behind the table.

For (p=17), Sticker's eight complete 13-arc classes have index-one counts
(0,0,0,0,0,0,2,3).  If (A) is incomplete, extend it to a 14-arc (K).  If
(K) is complete, it belongs to the unique complete-14 class; auditing all
fourteen subarcs (K\setminus\{x\}) gives (c_1\le4).  If (K) is incomplete,
extend it once more to a 15-arc.  The unique 15-arc class is contained in a
conic.  Relative to a 13-subarc of that conic, a missing conic point has
secant index zero, while an off-conic point retains at least
((17-1)/2-5=3) secants.  Thus this last branch has (c_1=0).

For (p=19), Al-Zangana's exhaustive classification of all 83 projective
14-arc classes, including 13 incomplete and 70 complete classes, gives
(c_1\le4) for every class.

For (p=23), the five complete 17-arc classes have index-one counts
(0,0,1,0,0).  An incomplete 17-arc extends to an 18-arc.  The
Coolsaet--Sticker complete-arc spectrum has no complete size from 18 through
23 and a unique size-24 class, the conic.  A missing conic point again has
index zero, while an off-conic point retains at least
((23-1)/2-7=4) secants.  Hence this branch has (c_1=0).

For (p=29), Coolsaet--Sticker's exhaustive complete-arc classification has
exactly two complete 21-arc classes.  Direct exact audits of their published
representatives give outside secant-index histograms

\[
 \{4:18,5:75,6:190,7:312,8:189,9:63,10:3\}
\]

and

\[
 \{3:3,4:21,5:66,6:187,7:294,8:243,9:27,10:9\},
\]

so both have (c_1=0).  If (A) is incomplete, extend it until it is complete.
Segre's odd-order bound caps an arc at (p+1=30), and the classified spectrum
has no complete sizes 22, 23, or 25 through 29.  Thus the completion has size
24 or 30.  The unique complete 24-arc is the Klein quartic

\[
 x^3y+y^3z+z^3x=0.
\]

Its exact outside secant-index histogram is

\[
 \{6:28,8:126,9:504,10:84,11:84,12:21\}.
\]

Thus for a 21-subarc, a missing Klein point has index zero and a point
outside the Klein arc retains at least (6-3=3) secants.  The unique complete
30-arc is a conic; a 21-subarc has index zero at a missing conic point, while
an off-conic point retains at least ((29-1)/2-9=5) secants.  Every incomplete
branch therefore also has (c_1=0).

It follows that

\[
 \boxed{R=5,6,7,9\text{ is impossible at }p=17,19,23,29,
 \text{ respectively}.}                          \tag{15.727.11}
\]

The current first possible positive slacks at these primes are therefore
(6,7,8,10).  The first prime whose endpoint remains unexcluded here is
(p=31), with (R=10).  No new solver run enters this proposition.  From
(p=31) onward the endpoint is left in the disjoint 3/4-secant normal form;
larger slack, the rest of residual (ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15727.py`,
`tests/test_prop15727.py`,
`evidence/e1_gmin_m4_prop15727.json`, and
`evidence/NOTE_2026-08-30_endpoint_rigidity_small_prime_close.md`.

## Proposition 15.728 — Paley-hard rigidity at the `p=31` endpoint

Specialize the still-open endpoint of Proposition 15.727 to (p=31), so
(|D|=32), (R=10), and (|H|=4p+1=125).  The hypotheses here include the
residual affine separator inequalities of Proposition 15.632; this is not a
claim about an arbitrary 32-point set.  For a direction (d), let (b_d) be
the number of its affine fibres that meet (D) oddly.  Since (|D|) is even,
every (b_d) is even, and the exact pair identity gives

\[
 \boxed{\sum_d b_d=32+4R=72.}                       \tag{15.728.1}
\]

The directional parity sign of Proposition 15.632 is

\[
 \epsilon_d(-1)^{(|H|-3)/2}c_H(-1)^{b_d}
   =-\epsilon_dc_H,                                  \tag{15.728.2}
\]

because ((|H|-3)/2=61) is odd.  Thus the sixteen directions with
(\epsilon_d=c_H) all have phase one, and the other sixteen all have phase
zero.  Both quadratic types have the exact scaled-mean budget

\[
 {32^2\over2}=512.                                   \tag{15.728.3}
\]

For even (b), the exact phase-one floors at (p=31) are

\[
 L_1(b)=
 \begin{cases}
  30,&b=2,30,\\
  62,&b=0,4,6,\ldots,28.
 \end{cases}                                         \tag{15.728.4}
\]

Within one quadratic type the exact means have a common residue modulo 32.
Write

\[
 a_d=2u+32k_d,qquad 0\le u<16.
\]

Summing over the type and using (15.728.3) yields

\[
 \sum_d k_d=16-u.                                    \tag{15.728.5}
\]

For (1\le u\le14), every direction needs (k_d\ge1), already exceeding the
right side of (15.728.5).  At (u=0), quotient one is possible only in the
(b=2,30) cells and gives mean 32, two above their floor.  In either cell,
the parity constraint has a pointwise Boolean baseline (q_0) of scaled mean
30: (q_0=(1-x_i-x_j)^2) for (b=2), and after complementing the 30-set,
(q_0=1-x_i) for (b=30).  A putative slack quadratic (A) of scaled mean 32
would make

\[
 C={A-q_0\over2}\ge0,\qquad 4p\,\mathbb EC=2.
\]

This is a nonzero integral quadratic, whereas Proposition 15.688 gives the
sharp bound (4p\,\mathbb EC\ge p-3=28).  Thus both floor-plus-two lifts are
impossible.  Hence (u=15).  Equation (15.728.5) then has quotient
sum one, and (15.728.4) forces

\[
 \boxed{\{a_d:\epsilon_d=c_H\}=\{30^{15},62^1\}.}   \tag{15.728.6}
\]

The other type is also arithmetically restricted.  If its common residue is
(2u_0), the exact formula
(a_d=I+32P_d-\epsilon_dT-93) gives, after adding the two type residues,

\[
 30+2u_0\equiv2I+6\pmod {32}.
\]

Infinity is not in the boundary, so its degree (I) is even.  Hence
(u_0\equiv I+4\pmod {16}) and

\[
 u_0\in\{0,2,4,6,8,10,12,14\}.                 \tag{15.728.6a}
\]

In particular, fifteen phase-one directions have (b_d\in\{2,30\}) at
their mean-30 baseline; the single mean-62 direction may have any even
(b_d), including an elevated (b_d=2) or (30).  There cannot be two
(b_d=30) directions in this type: together with the other thirteen required
baseline directions at (b_d=2), they would contribute

\[
 2\cdot30+13\cdot2=86>72
\]

to (15.728.1).  Therefore

\[
 \boxed{\text{one Paley type has at least fourteen }b_d=2
        \text{ directions}.}                       \tag{15.728.7}
\]

This also sharpens the disjoint-block geometry.  In a (b_d=2) direction,
let (r_3,r_4) count the trisecants and 4-secants, and let (l_j) count the
(j)-point fibres.  The point and parity identities give

\[
 l_0=14+r_3+r_4,\quad l_1=2-r_3,\quad
 l_2=15-r_3-2r_4,\quad l_3=r_3,\quad l_4=r_4.
                                                               \tag{15.728.8}
\]

If (y) is the total number of 4-secants, Proposition 15.727 has
(x=10-2y), hence only (x+y=10-y) rich lines.  They occupy at most
(10-y) directions.  Equations (15.728.7)--(15.728.8) therefore force at
least

\[
 \boxed{4+y\text{ directions of one Paley type with profile }
        (l_0,l_1,l_2,l_3,l_4)=(14,2,15,0,0).}        \tag{15.728.9}
\]

This is a proved necessary normal form, not an endpoint exclusion.  It uses
no arc classification and no finite configuration search.  At this stage,
the next implication was to show that the disjoint 3/4-secant endpoint cannot
support the same-Paley near-perfect pairing directions in (15.728.9), or that
they force one of the
already-closed circle/conic configurations.  The (p=31) endpoint, the rest
of residual (ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15728.py`,
`tests/test_prop15728.py`,
`evidence/e1_gmin_m4_prop15728.json`, and
`evidence/NOTE_2026-08-31_p31_endpoint_paley_hard_profile.md`.

## Proposition 15.729 — affine unique-trisecant endpoint reduction

Continue the all-prime endpoint normal form of Proposition 15.727.  Thus
(D) is an affine set of (p+1) points, (p=3R+c) with (c\in\{1,2\}), and its
rich lines are pairwise (D)-disjoint: (x) trisecants and (y) 4-secants with

\[
 x+2y=R.                                             \tag{15.729.1}
\]

Choose one rich line.  Retain three points on it, two points on every other
rich line, and every point outside the rich blocks.  Call the resulting set
(U).  If the distinguished line is a trisecant, the number deleted is
((x-1)+2y=R-1).  If it is a 4-secant, the number deleted is
(1+x+2(y-1)=R-1).  Hence

\[
 |U|=p+2-R.                                         \tag{15.729.2}
\]

The set (U) is affine because (U\subset D).  Moreover, every line containing
three points of (U) is already a rich line of (D).  The distinguished rich
line retains three points and every other one retains two.  Therefore

\[
 \boxed{U\text{ is an affine }(p+2-R,3)\text{-arc with exactly one
 trisecant}.}                                      \tag{15.729.3}
\]

Write that trisecant as (\{P,Q,Z\}) and put (B=U\setminus\{P,Q\}).  Then
(B) is an affine arc of size (p-R).  Both (B\cup\{P\}) and (B\cup\{Q\})
are arcs: otherwise the new trisecant would already be a second trisecant
of (U).  The line (PQZ) meets (B) only in (Z), so it is a tangent of (B) at
(Z).  Thus

\[
 \boxed{B\text{ has two distinct affine extension points }P,Q
 \text{ on one tangent}.}                         \tag{15.729.4}
\]

Any of the three points on the unique trisecant may be chosen as the point
left in (B).  The exact sizes in the two residues are

\[
\begin{array}{c|cc}
 &|U|&|B|\\ \hline
 p=3R+1&2R+3&2R+1\\
 p=3R+2&2R+4&2R+2.
\end{array}                                        \tag{15.729.5}
\]

This is a proved structural reduction, not an endpoint exclusion.  It uses
no blocking-set theorem, finite configuration search, or new classification.
At the 15.729 stage, the remaining implication was to exclude or classify
these two families of large affine unique-trisecant 3-arcs, equivalently
their large affine arc subsets with two co-tangent extensions, while
preserving compatibility with the common disjoint-block completion (D).
Propositions 15.730--15.731 supersede that formulation with the full repair
ensemble and its transition coordinates. The endpoint, the rest of residual
(ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15729.py`,
`tests/test_prop15729.py`,
`evidence/e1_gmin_m4_prop15729.json`, and
`evidence/NOTE_2026-08-31_affine_unique_trisecant_reduction.md`.

## Proposition 15.730 — complementary-arc repair ensemble and exact census

Continue the endpoint normal form of Proposition 15.727 and put
(k=p+1-R).  Thus (D=A\sqcup T), where (A) is a (k)-arc, (|T|=R), and
the rich lines of (D) are pairwise (D)-disjoint: (x) trisecants and (y)
4-secants, each containing exactly two points of (A), with

\[
 x+2y=R.                                             \tag{15.730.1}
\]

Let (S) be the points outside the rich blocks.  Then
(|S|=p+1-3x-4y=c+1+2y).  On every rich block (L), choose a two-point
subset (Q_L), and define

\[
 A_Q=S\cup\bigcup_L Q_L,\qquad T_Q=D\setminus A_Q.  \tag{15.730.2}
\]

Any line containing three points of either (A_Q) or (T_Q) would be a rich
line of (D).  But (A_Q) retains exactly two points on every rich block,
while (T_Q) retains one on a trisecant and two on a 4-secant.  Therefore
both sets are arcs, with

\[
 |A_Q|=k=p+1-R,\qquad |T_Q|=R.                      \tag{15.730.3}
\]

Conversely, a (k)-point arc in (D) must delete at least one point on each
trisecant and at least two on each 4-secant.  Those disjoint demands already
total (x+2y=R), so equality forces precisely the choices in
(15.730.2).  Hence

\[
 \boxed{D\text{ has exactly }3^x6^y\text{ maximum arc repairs}.}
                                                               \tag{15.730.4}
\]

For every repair and every (z\in T_Q), the rich block containing (z)
supplies one (A_Q)-secant through it.  A second such secant would be a
second rich block through (z), contrary to disjointness.  Thus all (R)
points of (T_Q) have secant index one, the unique secants form a matching
on (A_Q) with fibre sizes one and two, and each
(A_Q\cup\{z\}) is an affine ((k+1,3))-arc with exactly one trisecant.

There is also an exact two-colour census.  If

\[
 n_{ij}=\#\{\ell:|\ell\cap A_Q|=i,
                    |\ell\cap T_Q|=j\},
\]

then

\[
\begin{array}{c|ccc}
 i\backslash j&0&1&2\\ \hline
0&\frac{p(p-1)}2-R-y&2R+2y&\binom R2-y\\
1&k+2R&R(k-2)&0\\
2&\binom k2-R+y&x=R-2y&y.
\end{array}                                             \tag{15.730.5}
\]

Indeed, both colours are arcs and the rich-line classification forbids type
((1,2)).  Counting pairs of each colour and cross-colour pairs gives

\[
 n_{20}+n_{21}+n_{22}=\binom k2,\quad
 n_{02}+n_{22}=\binom R2,\quad
 n_{11}+2n_{21}+4n_{22}=kR.                         \tag{15.730.6}
\]

Each point of the (k)-arc has (R+1) tangents, and each point of the
(R)-arc has (k+1) tangents.  Hence

\[
 n_{10}+n_{11}=k(R+1),\qquad
 n_{01}+n_{11}+n_{21}=R(k+1),                       \tag{15.730.7}
\]

and (15.730.5) follows from (15.730.1), followed by the total line count.

The co-tangent consequence sharpens as well.  If a rich block contains the
(A_Q)-pair (\{a,b\}), then all points in
(\{a\}\cup(T_Q\cap L)) are individually valid, pairwise incompatible
extensions of (A_Q\setminus\{a\}) on the tangent (L) through (b).
Exchanging (a,b), every fixed repair therefore supplies (2x) bases with two
co-tangent extensions and (2y) bases with **three** co-tangent extensions.

Finally fix an affine direction (d).  Let (\sigma_d,\tau_d) count
(A_Q)- and (T_Q)-secants, let (r_3(d),r_4(d)) count rich blocks, and let
(m_d) count ((1,1))-lines in that direction.  The nine affine line counts
are

\[
\begin{array}{lll}
n_{20}=\sigma_d-r_3-r_4,&n_{21}=r_3,&n_{22}=r_4,\\
n_{10}=k-2\sigma_d-m_d,&n_{11}=m_d,&n_{12}=0,\\
n_{02}=\tau_d-r_4,&n_{01}=R-2\tau_d-r_3-m_d,&
n_{00}=\sigma_d+\tau_d+r_3+r_4+m_d-1.
\end{array}                                             \tag{15.730.8}
\]

In particular the odd-fibre count is

\[
 b_d=p+1-2(\sigma_d+\tau_d+m_d),\qquad
 \sum_d b_d=p+1+4R.                                  \tag{15.730.9}
\]

This is a proved structural reduction, not endpoint closure.  It replaces
the invalid claim that Bartoli--Storme supplies a unique-trisecant size
ceiling: under the theorem's other hypotheses, including
\(d>3+2\sqrt q\) and existence of the configuration, their threshold is the
upper endpoint of a hyperplane-arrangement classification range and is not
such an existence bound.  Proposition 15.731 next
constructs the tangent envelopes directly and reduces compatibility across
one-block swaps to quadratic or cubic transition data.  The endpoint,
larger slack, residual (ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15730.py`,
`tests/test_prop15730.py`,
`evidence/e1_gmin_m4_prop15730.json`, and
`evidence/NOTE_2026-08-31_endpoint_repair_ensemble.md`.

## Proposition 15.731 — endpoint tangent envelopes and low-degree swaps

Fix one of the maximum repairs (A) from Proposition 15.730 and put

\[
 n=|A|=p+1-R,qquad t=p+2-n=R+1,qquad d=2t=2R+2.
                                                               \tag{15.731.1}
\]

For (a\in A), let (f_a(X)) be the product of the (t) tangent-line forms at
(a), scaled by Segre's lemma of tangents so that

\[
 f_a(b)=(-1)^{t+1}f_b(a).                            \tag{15.731.2}
\]

Write (a^*=\{Z:a\cdot Z=0\}) and (L_a(Z)=a\cdot Z).  The map
(X\mapsto X\mathbin\times a) identifies the quotient by (a) with (a^*).
Since every factor of (f_a) vanishes at (a),

\[
 h_a(X\mathbin\times a)=f_a(X)^2                  \tag{15.731.3}
\]

is a well-defined degree-(d) section on (a^*).  At
(a^*\cap b^*) the two prescribed values are (f_a(b)^2) and (f_b(a)^2);
they agree by (15.731.2).  The dual lines have no triple intersection
because (A) is an arc.

We use the following elementary gluing fact.  Compatible degree-(d)
sections on any family of projective lines with no triple intersection
extend to a degree-(d) plane polynomial.  Inductively, the error on a new
line vanishes at all old intersection points.  If there are at most (d) of
them, divide by their distinct linear factors and lift the quotient; if
there are more than (d), the error is identically zero.  The kernel of the
restriction map is

\[
 \left(\prod_i L_i\right)H^0(\mathbb P^2,\mathcal O(d-n)).
                                                               \tag{15.731.4}
\]

Thus (15.731.3) glues to a degree-(d) tangent-envelope polynomial
(\Phi_A(Z)).  The exact endpoint dimensions are

\[
\begin{array}{c|ccc}
 &n&d&\text{envelope space after fixing tangent normalization}\\ \hline
p=3R+1&2R+2&2R+2&
 \Phi_0+\lambda\prod_{a\in A}L_a,\\
p=3R+2&2R+3&2R+2&\text{one polynomial}.
\end{array}                                             \tag{15.731.5}
\]

Without fixing the common tangent normalization, the first row is the
projective pencil spanned by (\Phi_0) and the line product, with the pure
line-product point omitted; the second row is a unique projective envelope.
This applies to all (3^x6^y) repairs.  It is a consequence of the Segre
tangent lemma and line gluing, not an invocation of Ball--Lavrauw's explicit
interpolation theorem below its stated size hypothesis.

There is a low-degree transition law between adjacent repairs.  Suppose

\[
 A=C\cup\{a\},\qquad A'=C\cup\{z\}                 \tag{15.731.6}
\]

differ by one swap on a rich block.  For (u\in C), the products
(\det(u,a,X)f_u^A(X)) and (\det(u,z,X)f_u^{A'}(X)) are tangent products
for (C).  Their Segre normalizations differ by one common scalar; rescale
the primed family to align them.  Restriction to every (u^*) then gives

\[
 L_z^2\Phi_{A'}=L_a^2\Phi_A.
\]

Consequently

\[
 \boxed{L_z^2\Phi_{A'}-L_a^2\Phi_A
   =\left(\prod_{u\in C}L_u\right)Q_{a,z}.}          \tag{15.731.7}
\]

For (p=3R+2), (|C|=d) and (Q_{a,z}) is quadratic.  For (p=3R+1),
(|C|=d-1) and it is cubic; changing the two envelope representatives changes
its class only by an element of (\operatorname{span}\{L_a^3,L_z^3\}).

These relative scalings can be chosen coherently on the whole repair graph.
Every repair contains the common singleton set (S), with
(|S|=c+1+2y\ge2).  Fix (e\in S) and the unique line (\rho_e) through (e)
which avoids (D), and normalize

\[
 f_e^A=\rho_e\prod_{v\in D\setminus A}\det(X,e,v).    \tag{15.731.8}
\]

These are exactly the (R+1) tangent factors at (e).  Under the swap
(a\mapsto z), the right side changes by
(\det(X,e,a)/\det(X,e,z)), so the scalar in (15.731.7) is one on every
edge.  The factors telescope on closed repair-graph walks.  This does not
remove the line-product kernel in the (p=3R+1) residue.

This is a proved algebraic refinement, not endpoint closure.  The repair
graph is the Cartesian product of a (K_3) for every trisecant and a
(J(4,2)) for every 4-secant.  The exact open implication is exclusion of the
common completion under the residual direction, phase, and lift constraints.
At this stage, deriving a nontrivial cycle identity from the quadratic
transitions or cubic classes was the proposed next attack. Proposition
15.732 below proves that the naturally cleared linear circulation is an
exact coboundary and replaces that proposal by the local-jet / simultaneous
baseline fronts; no phase bridge is claimed here.
Endpoint equality, larger slack, residual (ii), Type I, and the limit remain
open.

Evidence: `src/e1_gmin_m4_prop15731.py`,
`tests/test_prop15731.py`,
`evidence/e1_gmin_m4_prop15731.json`, and
`evidence/NOTE_2026-08-31_endpoint_tangent_envelope_gluing.md`.

## Proposition 15.732 — repair-cycle exactness and the surviving local jet

Continue Propositions 15.730--15.731.  For every repair (A), put

\[
 P_A=\prod_{u\in A}L_u,\qquad
 \Theta_A=P_A^2\Phi_A.                              \tag{15.732.1}
\]

Orient an adjacent swap from (A=C\cup\{a\}) to
(A'=C\cup\{z\}).  The coherently normalized transition identity in
15.731 is

\[
 L_z^2\Phi_{A'}-L_a^2\Phi_A=P_CQ_{a,z}.             \tag{15.732.2}
\]

Because (P_A=P_CL_a) and (P_{A'}=P_CL_z), multiplication by (P_C^2)
gives the exact-potential identity

\[
 \boxed{\Theta_{A'}-\Theta_A=P_C^3Q_{a,z}.}          \tag{15.732.3}
\]

The degrees agree: if (p=3R+c), (c\in\{1,2\}), then
(|A|=k=2R+c+1), (\deg\Phi_A=2R+2), and
(\deg Q=4-c), so both sides of (15.732.3) have degree

\[
 (2R+2)+2k=3(k-1)+(4-c).                            \tag{15.732.4}
\]

It follows immediately that on every closed repair-graph walk

\[
 \sum_{A\longrightarrow A'}P_{A\cap A'}^3Q_{A,A'}=0.\tag{15.732.5}
\]

This is an identity of polynomials before any evaluation, coefficient
extraction, derivative, or polarization is applied.  Hence every additive
linear circulation obtained from the 15.731 edge laws is zero.  In the
(c=1) residue the change (\Phi_A\mapsto\Phi_A+\mu_AP_A) changes
(\Theta_A) by the vertex term (\mu_AP_A^3), and changes

\[
 Q_{a,z}\mapsto Q_{a,z}+\mu_{A'}L_z^3-\mu_AL_a^3.   \tag{15.732.6}
\]

Thus (15.732.5) remains an exact coboundary in the envelope-pencil residue;
independently choosing a representative for every edge class would discard
the vertex correlations rather than create holonomy.

For example, let one rich block contain (i,j,k), and write
(A_{ij}=B\cup\{i,j\}) with the outside-block choice (B) fixed.  The
(K_3) cycle, and likewise every triangular face of a (J(4,2)) factor, gives

\[
 L_i^3Q_{j,k}^{,i}+L_k^3Q_{i,j}^{,k}
       +L_j^3Q_{k,i}^{,j}=0.                       \tag{15.732.7}
\]

There are analogous cleared identities on the chordless squares of
(J(4,2)); all of them are instances of (15.732.5), not independent cycle
obstructions.

One nonzero local datum does survive.  Let (q) be the dual point of the rich
line containing the exchanged points (a,z) and the retained point (b).
Write (P_C=P_{C\setminus\{b\}}L_b).  The rich line is a secant of both
repairs, so (\Phi_A(q)) and (\Phi_{A'}(q)) are nonzero.  Taking the lowest
homogeneous term of (15.732.2) at (q) gives

\[
 Q_{a,z}(q)=0,
\]
\[
 P_{C\setminus\{b\}}(q)L_bj_q^1Q_{a,z}
   =\Phi_{A'}(q)L_z^2-\Phi_A(q)L_a^2.               \tag{15.732.8}
\]

The first jet is nonzero, since (L_a,L_z) are nonproportional cotangent
forms and the two displayed scalars are nonzero.  More explicitly, for the
block notation above, let (l_r) be the class of (L_r) in
(\mathfrak m_q/\mathfrak m_q^2), put
(\Delta_{rs}=\det(l_r,l_s)), and set

\[
 K=\Phi_{ij}(q)\Delta_{ij}^2.
\]

The scaled tangent restrictions show that (K) is independent of the
selected pair and is a nonzero square.  Then

\[
 j_q^1Q_i(j,k)=\frac{K}{P_B(q)}
 \frac{l_k^2/\Delta_{ik}^2-l_j^2/\Delta_{ij}^2}{l_i}.
                                                               \tag{15.732.9}
\]

The numerator in (15.732.9) is divisible by (l_i).  In the cubic residue,
the ambiguity (15.732.6) vanishes to order three at (q), so this first jet
is gauge invariant.  Its square character is the repair-coloured character
(\chi_p(P_B(q))).  The residual parity set, however, records the symmetric
difference of the (A)- and complementary-(T) tangent fibres, not either
coloured product separately.  None of 15.688 or 15.722--15.724/15.728
supplies the missing identification.

Two natural attempts to bridge that gap can now be ruled out.  First,
suppose a nonrich direction of (D) has the near-pairing profile

\[
 \frac{p-3}{2}\text{ empty fibres},\quad
 2\text{ singleton fibres},\quad
 \frac{p-1}{2}\text{ double fibres}.               \tag{15.732.10}
\]

Deleting (R) points destroys at most (R) double fibres.  Hence every repair
has at least ((p-1)/2-R) secants in this direction and at most

\[
 |A|-2\left(\frac{p-1}{2}-R\right)=R+2             \tag{15.732.11}
\]

tangents.  These are at most (R+2) distinct known zeros on the
direction-pencil line in the dual plane.  Since the envelope degree is
(2R+2), (15.732.11) is far below the
(2R+3) roots required to force a direction component.  At (p=31,R=10),
every one of 15.728's at least (4+y) nonrich Paley-hard directions gives at
most 12 roots on a degree-22 restriction.

Second, products and quotients of repair contributions reduce modulo
squares to the span of the selected-pair masks on each block.  On a
trisecant those masks are

\[
 110,\quad101,\quad011.
\]

They span only the even-weight subspace of (\mathbb F_2^3), so the full
block mask (111) cannot be recovered.  On a 4-secant the weight-two masks
do contain (1111) in their span.  Therefore a repair-product character
cannot reconstruct a factor which includes an unselected trisecant block.
For the first jet on one rich block this obstruction applies whenever
another trisecant is present; the sole-trisecant jet is the explicit
exception to this mask argument, not a proved phase bridge.

This is a proved method barrier, not endpoint exclusion.  The proposed bare
cycle attack is now closed: a successful continuation must use information
outside its linear circulation.  The exact next possibilities are to relate
the nonzero, gauge-invariant jet (15.732.9) to the signed residual lift, or
to exclude the common completion directly from the many exact Paley-hard
near-pairing directions.  Endpoint equality, larger slack, residual (ii),
Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15732.py`,
`tests/test_prop15732.py`,
`evidence/e1_gmin_m4_prop15732.json`, and
`evidence/NOTE_2026-08-31_repair_cycle_exactness.md`.

## Proposition 15.733 — the (p=31,R=10) endpoint is impossible

Continue Proposition 15.728 with (p=31), (|H|=125), and an all-finite
32-point boundary (D) at outside pair slack (R=10).  Let (I) be the number
of selected infinity edges.  For a direction (d), let (P_d) be the number
of selected finite edges parallel to (d).  The phase-one Paley type
(\epsilon_d=c_H) has fifteen directions of mean 30 and one direction of
mean 62.

We first record the coefficient consequence of equality in a directional
floor.  Put (z_s=2x_s-1), so that (\sum_sz_s=1) on the middle slice.  If a
phase-one (b=2) direction has odd fibres (\{i,j\}), its exact baseline is

\[
 A_d=(1-x_i-x_j)^2,\qquad
 \epsilon_dS_H=4+z_iz_j.                            \tag{15.733.1}
\]

Write (n_s) for the number of selected infinity edges ending in fibre (s)
and (K_{st}) for the signed sum of selected finite edges between fibres
(s,t).  Thus (\sum_sn_s=I), and

\[
 \epsilon_dS_H=P_d+\sum_sn_sz_s
       +\epsilon_d\sum_{s<t}K_{st}z_sz_t.           \tag{15.733.2}
\]

Two multilinear quadratics agreeing on (\sum z_s=1) differ by
((\sum z_s-1)(c+\sum_sa_sz_s)).  Comparing the constant and linear
coefficients in (15.733.1)--(15.733.2) gives

\[
 I+P_d-4=30c.                                      \tag{15.733.3}
\]

Polarization of the integral coefficients makes (2c) integral.  Therefore

\[
 15\mid I+P_d-4.                                   \tag{15.733.4}
\]

The same calculation for a target (4+\sigma z_j) gives

\[
 15\mid I+P_d-(4+\sigma).                          \tag{15.733.5}
\]

All fifteen mean-30 directions have one common (P_d=P), because their
directional mean formula is

\[
 a_d=I+32P_d-\epsilon_dT-93.                       \tag{15.733.6}
\]

Proposition 15.728 already gives at least fourteen (b=2) baselines.  A
phase-one (b=30) baseline is (A_d=1-x_j), hence has target (4-z_j) and,
by (15.733.5), would require (15\mid I+P-3).  This is incompatible with
(15.733.4).  Consequently

\[
 \boxed{\text{all fifteen mean-30 hard directions have }b_d=2.}
                                                               \tag{15.733.7}
\]

Put

\[
 \rho={I+P-4\over15},\qquad s=\rho+P.              \tag{15.733.8}
\]

The hard mean-62 direction has parallel count (P+1) by (15.733.6).  The
number of finite selected edges of the opposite sign is therefore

\[
 E_{\rm opp}=125-I-(16P+1)=15(8-s).                 \tag{15.733.9}
\]

Here (\rho\ge0).  Infinity is outside (D), so (I) is even; from
(I=15\rho+4-P), (s) is even.  Thus

\[
 s\in\{0,2,4,6,8\}.                                \tag{15.733.10}
\]

For an opposite-type direction with (Q_d) parallel selected edges,
(15.733.6) and the mean-30 hard identity give

\[
 a_d=30s-208+32Q_d,qquad
 \sum_{epsilon_d=-c_H}Q_d=15(8-s).                \tag{15.733.11}
\]

Suppose (s<8).  Nonnegativity of (a_d) forces (Q_d\ge7-s) in all sixteen
opposite directions.  Their total excess above that lower bound is only

\[
 15(8-s)-16(7-s)=8+s<16.                            \tag{15.733.12}
\]

Hence some direction has (Q_d=7-s) and

\[
 a_d=16-2s\in\{16,12,8,4\}.                        \tag{15.733.13}
\]

In phase zero every nonzero (b_d) has floor at least 32, so this direction
has (b_d=0).  Its parity is even, hence (A_d=2C) for a nonzero nonnegative
integral quadratic (C), and (a_d=4p\mathbb EC).  Proposition 15.688 gives
(4p\mathbb EC\ge p-3=28), contradicting (15.733.13).  Therefore

\[
 s=8,qquad E_{\rm opp}=0,qquad I=124-16P.          \tag{15.733.14}
\]

Nonnegativity gives (0\le P\le7).  The case (P=0) would have (I=124) and
only one finite selected edge.  If (U) is the set of endpoints of the
infinity edges and (F) the finite edge set, then
(D=U\mathbin\triangle\partial F), so

\[
 I\le |D|+2|F|=34,
\]

a contradiction.  Hence (1\le P\le7).  Equation (15.733.14) says every
finite selected edge has the hard sign.  All sixteen opposite directions
have (Q_d=0), and (15.733.11) gives mean 32 in each.

The exact phase-zero floor table permits only (b_d=0,2,30) at mean 32.  The
(b=2) equality is (A_d=(x_i-x_j)^2), with target (4-z_iz_j), so
(15.733.4) would require

\[
 15\mid I-4.                                       \tag{15.733.15}
\]

The (b=30) equality is (A_d=x_j), with target (4+z_j), so (15.733.5)
would require

\[
 15\mid I-5.                                       \tag{15.733.16}
\]

But (15.733.14) gives (I\equiv4-P\pmod {15}), and neither
(15.733.15) nor (15.733.16) can hold for (1\le P\le7).  Thus all sixteen
opposite directions have (b_d=0).

Finally the exact global identity is (\sum_db_d=72).  The fifteen hard
baselines in (15.733.7) contribute 30 and the opposite type contributes
zero.  The single hard mean-62 direction would therefore need

\[
 b_d=72-30=42,
\]

impossible: among 31 fibres, (b_d) is even and at most 30.  We conclude

\[
 \boxed{\text{the residual }p=31,R=10\text{ endpoint does not exist}.}
                                                               \tag{15.733.17}
\]

This is a symbolic coefficient exclusion, not a finite configuration
search.  In particular, all fifteen baseline hard directions have (b=2),
so before the contradiction at least (5+y) of them are nonrich in a block
row with (y) 4-secants, improving 15.728's (4+y) count.  The stronger
entrywise cell audit in the certificate also reduces the intermediate
possibilities to

\[
 (I,P)=(28,6),(60,4),(92,2),
\]

but that refinement is unnecessary for the final contradiction.  The first
unexcluded endpoint prime is now (p=37,R=12); at (p=31) only larger slack
(R\ge11) remains.  The endpoint for all primes, the full (p+1) shell,
residual (ii), Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15733.py`,
`tests/test_prop15733.py`,
`evidence/e1_gmin_m4_prop15733.json`, and
`evidence/NOTE_2026-08-31_p31_simultaneous_baseline_close.md`.

## Proposition 15.734 — the critical residual-(ii) layer is empty from p=13

Let a critical residual-(ii) witness have original size (k=4p), and let
(H=G\cup\{e\}), so

\[
 |H|=4p+1.                                             \tag{15.734.1}
\]

No hypothesis on the odd-degree boundary (D=\partial H) will be used.
The graph (H) has at most (8p+2) incident vertices among the
(p^2+1) projective vertices. For (p\ge13),

\[
 p^2+1-(8p+2)=p^2-8p-1>0.                            \tag{15.734.2}
\]

Choose an isolated vertex (w). Since (w\notin D), the signed PSL
transport of Proposition 15.721 can send (w) to infinity while preserving
(|H|) and both separator inequalities. In the transported chart,

\[
 \infty\notin D,\qquad I=\deg_H(\infty)=0.             \tag{15.734.3}
\]

The boundary is now all finite. Handshake makes (|D|) even, and hence
every directional odd-fibre count (b_d\equiv |D|\pmod2) is even.

Put (q=(p-1)/2) and (m=q+1). Proposition 15.632 gives the exact budget

\[
 \sum_{\epsilon_d=\tau}a_d=2m^2,qquad
 a_d=2p\,\mathbb E A_d,                               \tag{15.734.4}
\]

for either quadratic direction type. Since
((|H|-3)/2=2p-1) is odd, (15.632.3) and the evenness of (b_d) say that
the type (\epsilon_d=c_H) has phase one. Its means have one common residue
modulo (p+1=2m); write

\[
 a_d=2u+2mk_d,\qquad 0\le u<m,qquad
 \sum_d k_d=m-u.                                      \tag{15.734.5}
\]

The exact phase-one floor is at least (p-1). Thus an interior residue
(1\le u\le m-2) would require all (m) quotients (k_d\ge1), while
their sum is below (m). At (u=0), every mean is (p+1). A (b=2)
cell is a forbidden floor-plus-two lift by Proposition 15.688, so this row
survives only when (p\equiv1\pmod4), with every direction in the exact
(b=p-1) baseline. At (u=m-1), there are (q) low means (p-1) and
one high mean (2p). The low cell is (b=2) when (p\equiv1\pmod4),
and is (b=2) or (b=p-1) when (p\equiv3\pmod4). In the latter case
equal means force one common parallel count, while the two coefficient
offsets below differ by one, so the cells cannot mix.

For (p\ge17) this is Proposition 15.669's exact even-(b) table. At
(p=13), direct evaluation of Proposition 15.632's three-variable LP gives

\[
\begin{array}{c|rrrrrrr}
b&0&2&4&6&8&10&12\\ \hline
\text{phase }0&0&14&20&26&24&26&12\\
\text{phase }1&26&12&26&24&26&20&14,
\end{array}                                           \tag{15.734.6}
\]

so the same endpoint-relevant classification holds.

The numerical LP values are not, by themselves, an equality
classification. Proposition 15.652 supplies the needed rigidity: its
positive quadrature has strictly positive weights at every (b=2)
intersection value, so equality forces the pointwise XNOR polynomial
((1-x_i-x_j)^2). Complementing a (b=p-1) set to its missing fibre and
using the two positive (b=1) weights forces (x_j) when
(p\equiv1\pmod4), and (1-x_j) when (p\equiv3\pmod4). Thus the exact
baseline polynomials below are valid at (p=13) as well as in the range of
Proposition 15.669.

We now reuse the coefficient comparison of Proposition 15.733. With
(z_s=2x_s-1) and (sum_s z_s=1), an exact target
(4+\tau z_iz_j) gives

\[
 q\mid I+P-4,                                         \tag{15.734.7}
\]

while (4+\sigma z_j) gives

\[
 q\mid I+P-(4+\sigma).                               \tag{15.734.8}
\]

The three exhaustive hard branches are therefore

\[
\begin{array}{c|c|c|c}
&A_d&\epsilon_dS_H&\text{offset}\\ \hline
A&(1-x_i-x_j)^2&4+z_iz_j&4\\
B\ (p\equiv1)&x_j&4+z_j&5\\
C\ (p\equiv3)&1-x_j&4-z_j&3.
\end{array}                                           \tag{15.734.9}
\]

Let (P) be the common hard parallel count, write
(\rho=(I+P-C)/q) for the relevant offset (C), and put (s=P+\rho).
Here (\rho\ge0): its numerator is divisible by (q) and is strictly
larger than (-q). Exact finite-edge accounting gives (q(8-s)) opposite
edges in A/B and (q(8-s)+1) in C. Hence (P\le s\le8). Since (I=0)
and (q\ge6), (15.734.7)--(15.734.8) force

\[
 (P,\rho,s)=(4,0,4),\quad(5,0,5),\quad(3,0,3)          \tag{15.734.10}
\]

in A, B, and C respectively.

For A/B, an opposite-type direction with parallel count (Q) has

\[
 a=(p-1)s+(p+1)Q+9-7p,qquad
 \sum Q=q(8-s).                                       \tag{15.734.11}
\]

For C the identities are

\[
 a=(p-1)(s-7)+(p+1)Q,qquad
 \sum Q=q(8-s)+1.                                     \tag{15.734.12}
\]

Nonnegativity forces (Q\ge3,2,4) in A, B, C. The surplus above those
uniform minima is (q-3,q-2,q-3), each below the number (m=q+1) of
opposite directions. Some direction therefore attains the minimum, with

\[
 a=8,\quad6,\quad8                                    \tag{15.734.13}
\]

respectively. That direction has phase zero. Every nonzero even (b_d)
costs at least (p-1\ge12), so (b_d=0). Its parity is even and its mean
is positive, hence (A_d=2B_d) for a nonzero nonnegative integer-valued
quadratic (B_d). Proposition 15.688 gives

\[
 a_d=4p\,\mathbb E B_d\ge p-3\ge10,                   \tag{15.734.14}
\]

contradicting (15.734.13). Thus

\[
 \boxed{\text{residual (ii) at }k=4p\text{ is empty for every prime }
 p\ge13,\text{ for every boundary size}.}             \tag{15.734.15}
\]

At (p=11), the same isolated chart and residue calculation force scaled
mass eight, exactly the sharp value (p-3); no equality classification is
claimed here. Even (k>4p), multi-level Type I, and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15734.py`,
`tests/test_prop15734.py`,
`evidence/e1_gmin_m4_prop15734.json`, and
`evidence/NOTE_2026-08-31_isolated_chart_k_eq_4p_close.md`.

## Proposition 15.735 — the next two residual-(ii) layers are empty from p=13

Let

\[
 k=4p+2t,\qquad |H|=4p+2t+1,\qquad t\in\{1,2\}.
                                                               \tag{15.735.1}
\]

As in Proposition 15.734, at most (2|H|) projective vertices are incident
with an edge.  For (p\ge13),

\[
 p^2+1-2|H|=p^2-8p-4t-1\ge56>0.                \tag{15.735.2}
\]

Transporting an isolated vertex to infinity therefore gives (I=0), an
all-finite even boundary, and even (b_d) in every direction.  Put
(q=(p-1)/2) and (m=q+1).  Each quadratic direction type has exact budget

\[
 2m(m+t).                                           \tag{15.735.3}
\]

The phase-one type is (\epsilon_d=(-1)^tc_H).  Its means have a common
residue,

\[
 a_d=2u+2mk_d,\qquad \sum_d k_d=m+t-u.             \tag{15.735.4}
\]

If (t<u\le m-2), every direction would require (k_d\ge1), although the
right side of (15.735.4) is below (m).  If (0\le u\le t), a direction
with (k_d=1) has mean (p+1+2u).  The only even-(b) cells below that
mean are the explicit (b=2) and (b=p-1) baselines.  Apart from
(u=0, p\equiv1\pmod4, b=p-1), they are nonzero integral lifts of excess
at most (2t+2\le6<p-3), contrary to Proposition 15.688.  Finally,
(u=m-1) forces at least (m-(t+1)) directions of mean (p-1).
Positive quadrature fixes their exact baselines, and offsets four and three
prevent the two (p\equiv3\pmod4) baselines from mixing.

Thus the same three hard branches as in Proposition 15.734 remain.  Their
parallel counts and hard finite-edge totals are

\[
\begin{array}{c|c|c|c}
 &P=s&\text{hard edges}&\text{opposite edges}\\ \hline
A&4&mP+t+1&4q+t\\
B&5&mP+t&3q+t\\
C&3&mP+t+1&5q+t+1.
\end{array}                                         \tag{15.735.5}
\]

Indeed the coefficient congruences are still
(q\mid P-4,P-5,P-3); opposite-edge nonnegativity leaves only the displayed
(\rho=0) rows.  Nonnegativity of an opposite mean gives respectively
(Q\ge3,2,4).  The surplus above these uniform minima is

\[
 q+t-3,\qquad q+t-2,\qquad q+t-3.                  \tag{15.735.6}
\]

For (t\le2), each is below the (m=q+1) opposite directions.  Some
direction attains the minimum and has scaled mean (8,6,8).  A nonzero
even (b) costs at least (p-1); at (b=0), a nonzero integral lift costs
at least (p-3).  Both bounds contradict those means for (p\ge13).  Hence

\[
 \boxed{k\in\{4p,4p+2,4p+4\}\text{ is impossible for every prime }
 p\ge13,\text{ for every boundary size}.}           \tag{15.735.7}
\]

The argument stops honestly at (t=3).  In branch B the surplus in
(15.735.6) equals (m), permitting every opposite direction to have
(Q=3) and mean (p+7).  Nonnegative integral quadratic examples attain
that local mass, and at (p=13) an additional exact (u=3,b=10) row
survives.  These are necessary-condition witnesses, not residual graphs;
they show that another one-direction floor or halving step cannot close
(k=4p+6).

Evidence: `src/e1_gmin_m4_prop15735.py`,
`tests/test_prop15735.py`,
`evidence/e1_gmin_m4_prop15735.json`, and
`evidence/NOTE_2026-08-31_first_three_residual_shells.md`.

## Proposition 15.736 — exact p=11 sharp Boolean-quadratic catalog

Let (\Omega=J(11,6)), so (|\Omega|=462).  The 55 pair monomials span all
degree-at-most-two functions on the slice, because

\[
 \sum_{j\ne i}x_ix_j=5x_i,\qquad
 \sum_{i<j}x_ix_j=15.                              \tag{15.736.1}
\]

Their (462\times55) evaluation matrix has rank 55 modulo 101, hence rank
55 over the reals.  Given a base 3-set and three disjoint swap pairs, take
the alternating sum on the resulting eight 6-sets.  This third difference
annihilates every quadratic.  Deterministic modular elimination examines
8,321 such rows and retains 407 independent rows.  Since

\[
 407=462-55,                                        \tag{15.736.2}
\]

their real nullspace is exactly the quadratic evaluation space.

Proposition 15.688 supplies the equality bridge required by the residual
problem.  At (p=11), sharp scaled mass is eight, whereas every nonnegative
integral lift with maximum at least two has scaled mass at least twelve.
Equality therefore forces maximum one: the lift is Boolean, has mean
(2/11), and has support (84) on (\Omega).

Now impose 462 Boolean variables (f_X), the 407 exact third-difference
identities, and

\[
 \sum_{X\in\Omega}f_X=84.                          \tag{15.736.3}
\]

Exclude the 220 known supports, one at a time: the 55 omitted-pair forms

\[
 (1-x_i)(1-x_j),                                   \tag{15.736.4}
\]

and the 165 all-equal-triple forms

\[
 1-x_i-x_j-x_k+x_ix_j+x_ix_k+x_jx_k.              \tag{15.736.5}
\]

The resulting exact 628-constraint CP-SAT model is infeasible.  Thus
(15.736.4)--(15.736.5) exhaust all sharp Boolean quadratic evaluations on
(J(11,6)).  This is an exhaustive finite certificate, with the linear
space and every no-good independently reconstructed; it does not use the
previously unavailable restriction-extension assertion.

In signed (z_i=2x_i-1) coordinates, the two targets are

\[
 4-z_i-z_j+z_iz_j\quad(\text{offset }2),\qquad
 4+z_iz_j+z_iz_k+z_jz_k\quad(\text{offset }4).      \tag{15.736.6}
\]

At (q=5), the hard-(b=2) residual branch forces a minimum opposite
parallel count (Q=3); neither offset in (15.736.6) is congruent to three,
so that branch is impossible.  The hard-(b=10) branch forces (Q=4): the
omitted-pair form is impossible and only the all-equal triple remains.  The
simultaneous incompatibility of those triple targets is supplied next.

Evidence: `src/e1_gmin_m4_prop15736.py`,
`tests/test_prop15736.py`,
`evidence/e1_gmin_m4_prop15736.json`, and
`evidence/NOTE_2026-08-31_p11_sharp_boolean_quadratic_catalog.md`.

## Proposition 15.737 — a binary quadratic moment closes the first three p=11 layers

Take (p=11), (t\in\{0,1,2\}), and

\[
 k=44+2t,\qquad |H|=45+2t.                         \tag{15.737.1}
\]

The isolated-vertex gaps are (32,28,24), so signed transport again gives
(I=0) and even directional (b).  The phase-one type has budget
(12(6+t)).  Writing (a_d=2u+12k_d) gives

\[
 \sum_d k_d=6+t-u.                                 \tag{15.737.2}
\]

The exact phase-one even-(b) floors are

\[
\begin{array}{c|rrrrrr}
b&0&2&4&6&8&10\\ \hline
f_1(b)&22&10&22&18&22&10.
\end{array}                                         \tag{15.737.3}
\]

For (0\le u\le t), a forced low direction is a positive lift of excess
at most six, below the sharp floor eight.  For (t<u<5), (15.737.2) has too
few quotient units.  Hence (u=5), with at least (5-t\ge3) exact
mean-ten directions.  Positive quadrature makes each either the (b=2)
baseline (4+z_iz_j), or the (b=10) baseline (4-z_j); offsets four and
three forbid mixing.  Proposition 15.736 excludes the first branch.  In the
second, at least (4-t\ge2) opposite minimum directions have an all-equal
triple target.

For an (\mathbb F_{11})-linear fibre functional (L), define the
homogeneous binary quadratic

\[
 M_H(L)=\sum_{\{u,v\}\in H}\chi(u-v)(L(u)-L(v))^2. \tag{15.737.4}
\]

If (K^L_{st}) is the signed selected-edge sum between fibres (s,t), then

\[
 M_H(L)=\sum_{s<t}K^L_{st}(s-t)^2.                 \tag{15.737.5}
\]

In a hard mean-ten (b=10) direction, (P=3) and
(\epsilon_LS_H=4-z_j).  Exact coefficient comparison modulo
(sum_sz_s=1) gives

\[
 \epsilon_LK^L_{st}=
 \begin{cases}-1,&\text{exactly one of }s,t\text{ is }j,\\0,&\text{otherwise}.
 \end{cases}                                        \tag{15.737.6}
\]

Therefore

\[
 M_H(L)=-\epsilon_L\sum_{t\ne j}(j-t)^2=0
 \quad\text{in }\mathbb F_{11}.                    \tag{15.737.7}
\]

There are at least three distinct hard directions of this kind.  A nonzero
homogeneous binary quadratic has at most two projective zeros, so
(M_H\equiv0).

For an opposite all-equal target, (P=4), and coefficient comparison gives
(\epsilon_LK^L_{st}=1) on the triangle ({i,j,k}) and zero elsewhere.
The identity (M_H\equiv0) would then require

\[
 (i-j)^2+(i-k)^2+(j-k)^2=0.                        \tag{15.737.8}
\]

For distinct fibres normalize ((i,j,k)=(0,1,r)).  The left side is
(2(r^2-r+1)), whose discriminant is (-3=8).  Since 8 is a nonsquare
modulo 11, (15.737.8) is impossible.  Thus even one triple target is
forbidden, and

\[
 \boxed{p=11:\quad k=44,46,48\text{ are all impossible}.}   \tag{15.737.9}
\]

Together with Proposition 15.735, the first three even residual layers are
therefore closed for every prime (p\ge11), independently of boundary
size.  At (p=11,t=3), only two hard stars are forced and excess eight
reaches the equality floor, so no (k\ge50) claim is made.  Critical
(p=5,7), all later residual layers, multi-level Type I, and the limit remain
open.

Evidence: `src/e1_gmin_m4_prop15737.py`,
`tests/test_prop15737.py`,
`evidence/e1_gmin_m4_prop15737.json`, and
`evidence/NOTE_2026-08-31_p11_binary_quadratic_moment_close.md`.

## Proposition 15.738 — exact p=13 mass-fourteen residual-cell catalog

Let (B) be a nonzero nonnegative integer-valued quadratic on
(J(13,7)) with

\[
 4p\mathbb E B=14.                                  \tag{15.738.1}
\]

The two inequalities in Proposition 15.688 first make the possible maximum
(H=\max B) completely discrete.  If (H\ge2), then

\[
 14\ge28-4H,\qquad 14\ge3H,
\]

so (H=4).  Otherwise (H=1) and (B) is Boolean.

The residual application supplies more than (15.738.1).  Its relevant
opposite directions have parallel count (Q=0) or (Q=6), and if
(W=\epsilon K) is the normalized signed inter-fibre matrix, then

\[
 \sum W=13(Q-3)-14,qquad
 \sum|W_{st}|\le59-Q,qquad
 \sum_{t\ne s}W_{st}\equiv0\pmod2,                \tag{15.738.2}
\]

and

\[
 B(X)={Q+\sum W-3-2\operatorname{cut}_W(X)\over4}
 \quad( X\in J(13,7)).                             \tag{15.738.3}
\]

For (H=4), permute a maximizing seven-set to a fixed anchor.  The two exact
integer models obtained from (15.738.2)--(15.738.3), one for each value of
(Q), impose (0\le B(X)\le4) on all 1,716 seven-sets and (B=4) at the
anchor.  Both are infeasible.  This is an exhaustive residual-cell
certificate; the permutation anchor is safe because every constraint before
the later field moments is (\operatorname{Sym}(13))-invariant.  Hence every
residual-compatible mass-fourteen cell is Boolean, of support

\[
 {7\over26}{13\choose7}=462.                       \tag{15.738.4}
\]

It remains to classify those Boolean evaluations without importing a slice
restriction theorem.  The 78 pair monomials have rank 78 modulo 101 on
(J(13,7)).  Third differences on a fixed four-set plus three disjoint swap
pairs annihilate every quadratic.  Deterministic elimination selects 1,638
independent eight-term identities, exactly the annihilator dimension.
Consequently their real nullspace is precisely the degree-at-most-two
evaluation space.

There are 1,092 explicit Boolean quadratic supports of size 462:

\[
\begin{array}{c|c|c}
\text{family}&\text{count}&\text{target }3+4B\text{ and offset}\\ \hline
x_ix_j&78&4+z_i+z_j+z_iz_j,\quad6\\
x_i(1-x_j)&156&4+z_i-z_j-z_iz_j,\quad4\\
z_i=z_j=-z_k&858&4+z_iz_j-z_iz_k-z_jz_k,\quad4.
\end{array}                                         \tag{15.738.5}
\]

An exact CP-SAT model has one Boolean variable per seven-set, the 1,638
identities, support 462, and no-goods for (15.738.5).  After the safe anchor
that one fixed support point is present, it is infeasible.  Thus
(15.738.5) is exhaustive.  The coefficient congruence is modulo six, and
both (Q=0) and (Q=6) retain only the first family.  Therefore

\[
 \boxed{B=x_ix_j\text{ in every relevant minimum opposite cell}.}
                                                               \tag{15.738.6}
\]

For either slice gauge, its normalized even moments are

\[
 \sum_{s<t}W_{st}(s-t)^2=(i-j)^2,qquad
 \sum_{s<t}W_{st}(s-t)^4=(i-j)^4.                 \tag{15.738.7}
\]

The two gauges differ by the complete-graph coefficient vector, whose
degree-two and degree-four moments vanish over (\mathbb F_{13}).
Proposition 15.738 is an **exhaustive finite certificate** for a local cell;
the simultaneous residual contradiction is the next proposition.

Evidence: `src/e1_gmin_m4_prop15738.py`,
`tests/test_prop15738.py`,
`evidence/e1_gmin_m4_prop15738.json`, and
`evidence/NOTE_2026-08-31_p13_mass14_boolean_catalog.md`.

## Proposition 15.739 — the exceptional p=13 fourth shell is impossible

Return to (p=13,t=3,u=3), where all seven hard directions have
(A=(2-r)^2) on a three-point complement (C).  The signed target is

\[
 \epsilon S_H
 =5-\sum_{i\in C}z_i+\sum_{\{i,j\}\subset C}z_iz_j. \tag{15.739.1}
\]

The coefficient offset is two, not five.  If (P) is the common hard
parallel count, slice-kernel polarization gives

\[
 6\mid P-2.
\]

Since the seven hard directions contain (7P) of the 59 finite edges,
(P\le8), and hence (P=2) or (P=8).  The exact opposite ledgers are

\[
\begin{array}{c|c|c}
P&\sum Q&a(Q)\\ \hline
2&45&14(Q-5)\\
8&3&14(Q+1).
\end{array}                                         \tag{15.739.2}
\]

For (P=2), a (Q=5,a=0) direction would have target three, whose offset
three is incompatible with (Q=5) modulo six.  Thus all seven (Q\)'s are at
least six and at least one is (Q=6,a=14).  For (P=8), at least one direction
is (Q=0,a=14).

At phase zero and mean fourteen, the (b=2) equality target has offset four,
incompatible with (Q=0,6).  A (b=12) cell would be a two-unit lift of its
pointwise parity baseline, below Proposition 15.688's lift floor ten.  The
remaining cell has (b=0), (A=2B), and is exactly the cell classified in
Proposition 15.738.  Hence one opposite direction has (B=x_ix_j).

For even (d), retain the genuine global homogeneous binary form

\[
 M_d(L)=\sum_{\{u,v\}\in H}\chi(u-v)(L(u)-L(v))^d. \tag{15.739.3}
\]

Let (h) be the fixed sign on the hard type.  The two gauges (P=2,8) in
(15.739.1) have the same even moments as the triangle on (C).  Normalize
(C=\{0,1,r\}) and put (q_0=r^2-r+1).  Then

\[
 S_2=2q_0,\qquad S_4=2q_0^2={1\over2}S_2^2.        \tag{15.739.4}
\]

Since (W=hK) on the hard type, the homogeneous quartic

\[
 G(L)=2hM_4(L)-M_2(L)^2                            \tag{15.739.5}
\]

vanishes on all seven hard projective directions.  A nonzero binary quartic
has at most four projective roots, so (G\equiv0).  On the opposite
(x_ix_j) cell, whose sign is (-h), (15.738.7) instead gives

\[
 G=-3(i-j)^4\ne0\quad\hbox{in }\mathbb F_{13},     \tag{15.739.6}
\]

a contradiction.  Therefore

\[
 \boxed{\text{the exceptional }p=13,t=3,u=3\text{ branch is empty}.}
                                                               \tag{15.739.7}
\]

There remains one generic (t=3) branch.  For (p=4a+1\ge17), at least
(N=m-3=(p-5)/2) hard directions are exact stars.  Applying (15.739.3) to
each even degree (d<N) shows

\[
 \boxed{M_d\equiv0\quad
 d=2,4,\ldots,{p-9\over2}.}                        \tag{15.739.8}
\]

Indeed an exact star has moment
(\sum_{r\in\mathbb F_p^*}r^d=0), and (N>d) roots force the degree-(d)
binary form to vanish.  Thus an opposite normalized matrix must have

\[
 \sum W=-(p+7),\quad \sum|W|\le4p+4,quad
 \operatorname{cut}_W(X)\le-{p+7\over2},quad
 \sum W_{st}(s-t)^d=0\pmod p                      \tag{15.739.9}
\]

for every degree in (15.739.8).  Conditional cut averages further force
(W_{st}\in\{-1,0,1,2,3\}) for every (p\ge17).  At (p=17), row sums are
even in ([-16,16]) and the recorded pair, triple, and four-set inequalities
give a substantially smaller exact model.  These claims are derived in the
certificate from exact conditional cut averages.

One further stabilizer average is useful.  Fix a nine-set (X), and write
(a,b,c) for the averaged coefficients on its 36 internal, 28 external, and
72 crossing edges.  Averaging the cut inequality over nine-sets meeting
(X) in four and five points gives

\[
 36a+28b+72c=-24,\quad
 20a+15b+37c\le-12,\quad
 20a+16b+36c\le-12.
\]

The combination (-9) times the four-intersection inequality, minus
(45/4) times the five-intersection inequality, plus (45/4) times the total
equality yields (72c\ge-27).  Since every row degree, and hence every cut,
is even,

\[
 -26\le\operatorname{cut}_W(X)\le-12.             \tag{15.739.10}
\]

Thus (B(X)=-6-\operatorname{cut}_W(X)/2) is an integer-valued quadratic
with values in (\{0,\ldots,7\}), mean (6/17), and total mass 8,580 on
(J(17,9)).  The full strengthened models are
currently unresolved, so (15.739.8)--(15.739.10) are a **proved open
reduction**, not a close of the generic branch.

The analogous (p=13) generic route genuinely stops at (M_2).  At an elevated
hard direction with (P=6), the signed matrix consisting of a positive
(K_5) on (\{0,1,2,3,5\}) plus the edge (\{0,11\}) has (\sum W=\|W\|_1=11),
odd rows (\{0,11\}), and every seven-cut between zero and seven.  Hence
(A=7-\operatorname{cut}_W\ge0), with scaled mean 28 and (b=2), while its
normalized moments are (S_2=0) and (S_4=5) modulo 13, so the global
moments are (M_2=0) and (M_4=5h) for the hard sign (h).  This is a **counterexample to
the one-direction moment/floor method**, not a common graph.  In particular
the whole (k=4p+6) shell and residual (ii) remain open.

Evidence: `src/e1_gmin_m4_prop15739.py`,
`tests/test_prop15739.py`,
`evidence/e1_gmin_m4_prop15739.json`, and
`evidence/NOTE_2026-08-31_p13_exceptional_quartic_close.md`.

## Proposition 15.740 — translation averages split the generic p=13 branch

In the remaining generic (p=13,t=3) row, write the seven hard means as
(a_L=14k_L).  The exact ledger is

\[
 k_L\ge1,\qquad \sum_L k_L=10.
\]

Thus, up to permutation, the hard quotient partition is one of

\[
 1^6 4,\qquad 1^5 2 3,\qquad 1^4 2^3.             \tag{15.740.1}
\]

Every exact (k_L=1) hard cell is a signed star.  Its even power sums vanish
over (\mathbb F_{13}).  Hence five exact hard directions give more than four
projective roots of both global binary forms and force

\[
 M_2\equiv M_4\equiv0.                            \tag{15.740.2}
\]

Consider any opposite direction.  It has (Q=3), mean 20, and (b=0).  For
its normalized signed matrix (W), coefficient comparison and
nonnegativity give

\[
 \sum W=-20,\quad \sum|W|\le56,\quad
 B(X)=-5-\frac12\operatorname{cut}_W(X)\ge0
 \quad (X\in J(13,7)).                             \tag{15.740.3}
\]

Conditioning a uniform seven-set to contain a fixed pair gives

\[
 \mathbb E[B\mid i,j\in X]={20+12W_{ij}\over44},
\]

so integrality implies (W_{ij}\ge-1).  For (a=1,\ldots,6), let (n_a) be
the sum of (W_{ij}) over the thirteen unordered cyclic pairs of difference
(\pm a).  Equations (15.740.2)--(15.740.3) imply

\[
 -13\le n_a\le18,\quad \sum_a n_a=-20,\quad
 \sum_a|n_a|\le56,
\]

\[
 \sum_a a^2n_a\equiv\sum_a a^4n_a\equiv0
 \pmod {13}.                                      \tag{15.740.4}
\]

For a seven-set (X\subset\mathbb F_{13}), put
(c_a(X)=|X\mathbin\triangle(X+a)|).  A fixed distance-(a) edge is separated
by exactly (c_a(X)) of the thirteen translates of (X).  Summing the thirteen
instances of (15.740.3) therefore gives

\[
 \sum_{a=1}^6c_a(X)n_a\le-130.                    \tag{15.740.5}
\]

The 1,716 seven-sets give exactly 74 distinct vectors (c(X)).  Exact modular
row reduction and bounded enumeration of (15.740.4) leave 32,313 aggregate
vectors; an independent meet-in-the-middle count gives the same number.
Greedily selecting the lexicographically first maximally eliminating vector
from the 74 leaves the exact remainder sequence

\[
 32313\to18091\to8124\to2037\to642\to225\to57\to12\to4\to0.
                                                               \tag{15.740.6}
\]

Thus nine instances of (15.740.5) already make the necessary six-variable
aggregate relaxation infeasible.  An independently encoded 14-variable,
19-constraint CP-SAT model returns exact status `INFEASIBLE`.  Since every
genuine opposite cell maps into this relaxation, the (1^6 4) and
(1^5 2 3) partitions in (15.740.1) are impossible.  Therefore

\[
 \boxed{p=13,t=3:\quad\text{only the hard partition }1^4 2^3
 \text{ remains}.}                                \tag{15.740.7}
\]

This does not close (p=13,k=58).  The remaining theorem must couple four
exact (P=5) stars, three elevated (P=6) hard cells, and seven
(Q=3,b=0) opposite cells through one common 59-edge graph.  The binary
affine-Radon reconstruction used for that coupling is already Proposition
15.692; it is imported, not renumbered.  No further local aggregate census
is an active route.

Evidence: `src/e1_gmin_m4_prop15740.py`,
`tests/test_prop15740.py`,
`evidence/e1_gmin_m4_prop15740.json`, and
`evidence/NOTE_2026-08-31_p13_generic_translation_average.md`.

## Proposition 15.741 — common-graph moments and the difference-Radon gate

Retain the last partition in (15.740.7), and now use that every directional
matrix comes from one 59-edge graph (H) on (V=\mathbb F_{13}^2).  Besides
(M_2) and (M_4), define the orientation-independent homogeneous forms

\[
 T_3(L)=\sum_{\{u,v\}\in H}\chi(u-v)
 (L(u)+L(v))(L(u)-L(v))^2,                         \tag{15.741.1}
\]

\[
 U_4(L)=\sum_{\{u,v\}\in H}\chi(u-v)
 (L(u)+L(v))^2(L(u)-L(v))^2.                      \tag{15.741.2}
\]

They have degrees three and four in (L).  In an exact hard direction the
normalized coefficient matrix is the positive star at a fibre (j).  Writing
(t=j+y), its four local contractions are

\[
 \sum_{y\ne0}y^2,quad
 \sum_{y\ne0}(2j+y)y^2,quad
 \sum_{y\ne0}y^4,quad
 \sum_{y\ne0}(2j+y)^2y^2,
\]

all zero in (\mathbb F_{13}).  The four distinct exact directions therefore
give

\[
 M_2\equiv T_3\equiv0.                            \tag{15.741.3}
\]

Both quartics vanish at those same four points, so they are proportional:

\[
 U_4=\lambda M_4.                                 \tag{15.741.4}
\]

Here (M_4\ne0).  Otherwise (M_2=M_4=0), and every required opposite cell
would enter the infeasible necessary relaxation of Proposition 15.740.
Thus the four exact directions are precisely the projective roots of (M_4),
and all three elevated and seven opposite cells have the same ratio
(\lambda=U_4/M_4).  Equations (15.741.3)--(15.741.4) are invariant under an
affine relabelling (s\mapsto as+c): the extra terms in (T_3) and (U_4) are
multiples of (M_2) and (T_3).

This consequence is sharp locally.  On one elevated fibre set put coefficient
one on (K_5) with vertices ({0,1,4,9,12}) and add one further unit on
({4,9}).  It has sum and (\ell^1)-norm eleven, odd rows ({4,9}), and every
seven-cut is at most seven.  Its contractions are

\[
 (M_2,T_3,M_4,U_4)=(0,0,7,10),\qquad U_4/M_4=7.
\]

For an opposite cell, put (w_{0x}=-1) for (x\ne0), put (w_{1x}=-1) for
(x\notin\{0,1,3,6\}), and put (w_{36}=1).  Its sum is (-20), its
(\ell^1)-norm is 22, every row is even, and every seven-cut lies in
({-14,-12,-10}).  Its contractions are

\[
 (M_2,T_3,M_4,U_4)=(0,0,8,4),\qquad U_4/M_4=7.
\]

These are directional coefficient cells, not a common graph.  They prove
only that the independent cellwise scalar consequences (M_2=T_3=0) and
(U_4=\lambda M_4) do not exclude either cell type.  They do not realize the
values of one common global quartic or the difference transform below.

There is nevertheless a stronger common-graph compression.  Let

\[
 \Omega=(V\setminus\{0\})/\{\pm1\},\qquad |\Omega|=84,
\]

and let (m_\delta) be the number of edges of (H) having unoriented
displacement (\delta).  For a projective functional (L), let
(\epsilon_L) be the character of its kernel and define, for
(a\in\mathbb F_{13}/\{\pm1\}),

\[
 q_L(a)=\epsilon_L\sum_{L\delta=\pm a}\chi(\delta)m_\delta. \tag{15.741.5}
\]

Then (q_L(0)=P_L), and the six nonzero entries are exactly the cyclic
distance aggregates of (W^L).  In particular, the full quartic value code is

\[
 \sum_{a=1}^6a^4q_L(a)=\epsilon_LM_4(L),
 \qquad M_4(L)=c\prod_{i=1}^4\det(L,L_i),\quad c\ne0.
\]

If (B) is the 98-by-84 incidence transform in
(15.741.5), two distinct displacement classes meet in one row when they are
projectively collinear and in two rows otherwise.  Hence

\[
 B^\mathsf TB=13I+2J-G,                            \tag{15.741.6}
\]

where (G) has fourteen diagonal (J_6)-blocks.  Put
(T=\sum_e\chi(e)), let (r) be the projective direction of (\delta), and
write (P_r) for its parallel count.  Multiplying (15.741.5) by
(B^\mathsf T) gives the exact integer inverse

\[
 13m_\delta=P_r-2\epsilon_rT+
 \epsilon_r\sum_L\epsilon_Lq_L(|L\delta|).         \tag{15.741.7}
\]

If (h) is the hard sign and (\sigma_L=\epsilon_L/h), then (T=17h), so the
right side of (15.741.7) lies in (13\mathbb Z_{\ge0}) pointwise.  Equivalently,

\[
 13m_\delta=P_r-34\sigma_r+
 \sigma_r\sum_L\sigma_Lq_L(|L\delta|).
\]

Parseval in (15.741.6), the parallel profile (5^4,6^3,3^7), and the four
exact rows (q=(2,2,2,2,2,2)) give, with
(C=\sum_\delta\binom{m_\delta}{2}),

\[
 \boxed{\sum_{L\ \mathrm{nonexact}}\sum_{a=1}^6q_L(a)^2
       =707+26C.}                                  \tag{15.741.8}
\]

Six multiplicative dilates of the interval seven-set have circulant cut
vector ((2,4,8,10,6,12)); the squared singular values are
(1764,76,84,100,84,76).  Applying their six cut inequalities and the exact
parities bounds an elevated row's energy by 86 and an opposite row's by 106.
Consequently

\[
 0\le C\le11.                                      \tag{15.741.9}
\]

This is not a contradiction.  The fractional assignment which spreads each
parallel count uniformly over its six displacement lengths gives exact rows
(2,\ldots,2), elevated rows ((11/6,\ldots,11/6)), and opposite rows
((-10/3,\ldots,-10/3)); all translated-cut inequalities hold strictly.
This point belongs only to the bare direction-count, exact-row, and
translated-cut relaxation: it is nonintegral, has (M_4=0), and does not test
midpoint identities.  Thus linear averaging of that bare system cannot
finish (15.741.7); integrality and the live quartic constraints are essential.

Finally, conditioning the opposite cut inequality gives, without a census,

\[
 w_{ij}\in\{-1,0,1,2,3\},\qquad
 d_i=\sum_jw_{ij}\in2\mathbb Z\cap[-12,10],qquad
 d_i+d_j\le6w_{ij}.                               \tag{15.741.10}
\]

The lower degree bound uses (w_{ij}\ge-1); the upper bound and pair inequality
come from conditioning on one or two vertices outside the seven-set.

The displacement transform forgets edge midpoints.  Writing an edge as
(\{m\pm\delta/2\}), equations (15.741.1)--(15.741.2) are exactly its first
and second signed midpoint moments.  Therefore the live p13 implication is
layered.  First exclude the 84-class nonnegative integer inverse (15.741.7),
with (15.741.8)--(15.741.10), the exact rows, the full nonzero quartic value
code, and all ten cut/parity lifts.  Only if it survives, lift its
(m_\delta) to the 14,196 binary midpoint variables and impose (T_3=0),
(U_4=\lambda M_4), simplicity, and the exact fibre-pair equations.  Proposition
15.741 is a **proved open reduction and method barrier**.  It does not close
the four-exact partition, (p=13,k=58), residual (ii), Type I, or the limit.

Evidence: `src/e1_gmin_m4_prop15741.py`,
`tests/test_prop15741.py`,
`evidence/e1_gmin_m4_prop15741.json`, and
`evidence/NOTE_2026-08-31_p13_common_graph_moment_transform.md`.

## Proposition 15.742 — six-dilate energy closes the generic p=13 row

Retain the six nonzero distance aggregates (q_L) from (15.741.5).  The
degree-two identity in (15.741.3) gives, in every direction,

\[
 \sum_{a=1}^6 a^2q_L(a)=0\pmod {13}.              \tag{15.742.1}
\]

Use the interval seven-set and its six nonzero multiplicative dilates.  In
the natural distance order their cut vectors are

\[
\begin{pmatrix}
2&4&6&8&10&12\\
12&2&10&4&8&6\\
8&10&2&6&12&4\\
6&12&8&2&4&10\\
10&6&4&12&2&8\\
4&8&12&10&6&2
\end{pmatrix}.                                    \tag{15.742.2}
\]

These are six members of Proposition 15.740's exact 74-vector translated-cut
catalog.  Hence every elevated nonexact row belongs to the integral
relaxation

\[
 \sum_aq_a=11,\quad \lVert q\rVert_1\le53,\quad
 \lVert q\rVert_2^2\le86,\quad (15.742.1),\quad C_6q\le91\mathbf1,
                                                               \tag{15.742.3}
\]

and every opposite row belongs to

\[
 \sum_aq_a=-20,\quad \lVert q\rVert_1\le56,\quad
 \lVert q\rVert_2^2\le106,\quad (15.742.1),\quad C_6q\le-130\mathbf1.
                                                               \tag{15.742.4}
\]

The (\ell^1) bounds count the nonparallel edges, namely (59-6) and
(59-3).  The energy bounds are the six-dilate spectral bounds of 15.741.
They are safe derived bounds: if (r=C_6q), then (r) is even; the elevated
slacks (91-r_i) are positive odd integers summing 84, while the opposite
slacks (-130-r_i) are nonnegative even integers summing 60.  Thus the
rational bounds are (4952/57) and (6050/57), whose integral floors are 86
and 106.

Now enumerate the integer six-tuples in (15.742.3)--(15.742.4).  Five
coordinates determine the sixth.  The energy bounds already give
(|q_a|\le9) and (|q_a|\le10), respectively, so the bounded enumeration is
manifestly exhaustive.  Before applying (C_6), there are 5,844 elevated
rows and 1,704 opposite rows.  Afterwards there are respectively 30 and 24,
and their sharp energy maxima are

\[
 \boxed{\max_E\lVert q\rVert_2^2=31,\qquad
        \max_O\lVert q\rVert_2^2=82.}             \tag{15.742.5}
\]

The six elevated maximizers are the multiplicative-distance images of
((0,3,1,4,1,2)); the six opposite maximizers are the corresponding images of
((-6,-1,-4,-2,-4,-3)).  Equivalently, they are cyclic after ordering the
coordinates as ((1,2,4,5,3,6)).  Direct integer enumeration gives the
complete-row hashes recorded in the evidence.  Independently encoded 19-variable exact
CP-SAT models with the extra requirements (\lVert q\rVert_2^2\ge32) and
(\ge83) both return `INFEASIBLE` with one worker, even after omitting the
prior energy caps 86 and 106 from those audit models.

There are three elevated and seven opposite rows.  Equation (15.742.5)
therefore bounds their total energy by

\[
 3\cdot31+7\cdot82=667.
\]

This contradicts the common-graph identity (15.741.8), whose right side is
(707+26C\ge707).  Consequently the four-exact partition (1^4 2^3) is empty.
Together with Proposition 15.739's exceptional close, this proves

\[
 \boxed{\text{the residual-(ii) row }p=13,k=58\text{ is empty}.}
\]

Proposition 15.742 is an **exhaustive finite certificate** at the six-bin
aggregate level.  It does not classify local coefficient matrices: the
contradiction uses their sharp necessary row-energy maxima only through the
common-graph Parseval identity.  The quartic value code, the four-root orbit
split, and the binary midpoint lift are unnecessary.  Residual (ii), Type I,
and the quadratic-minmax limit remain open.

Evidence: `src/e1_gmin_m4_prop15742.py`,
`tests/test_prop15742.py`,
`evidence/e1_gmin_m4_prop15742.json`, and
`evidence/NOTE_2026-08-31_p13_six_dilate_energy_close.md`.

## Proposition 15.743 — full translated-cut energy closes p=17,k=74

In the generic branch-B row at \(p=17,t=3\), there are nine hard and nine
opposite directions.  The hard quotient variables satisfy

\[
 k_L\ge1,\qquad \sum_{L\ \mathrm{hard}}k_L=12,
\]

so the excess partitions are \((3)\), \((2,1)\), and \((1,1,1)\), where
\(e=k_L-1\).  They leave respectively eight, seven, and six exact hard
stars.  In every case these roots force both homogeneous forms \(M_2\) and
\(M_4\) to vanish identically.

The normalization must be obtained before replacing an exact star by its
eight-bin row.  For an exact hard row \(k_L=1\), the unspecialized local and
common-graph sums give

\[
 17(P_L-3)-18=hT-P_L,
 \qquad\text{hence}\qquad hT=18P_L-69.           \tag{15.743.1}
\]

Thus all exact hard stars have one common \(P_L\).  There are at least six,
so their parallel edges give \(6P_L\le75\), hence \(P_L\le12\).  In the
isolated \(I=0\) chart, literal coefficient comparison gives
\(P_L\equiv5\pmod8\).  Consequently \(P_L=5\), and only now do we obtain
\(hT=21\) and the exact-star row \(q=(2)^8\).  Independently, every opposite
direction has \(Q=3\), so the hard and opposite edge totals are \(48\) and
\(27\), again giving \(hT=48-27=21\).

The parallel count of a general hard row is therefore not a free local
parameter.  Let \(P\) be that count.  The local coefficient identity gives

\[
 \sum_{a=1}^8q_L(a)=17(P-3)-18k_L.                \tag{15.743.2}
\]

Since all rows come from that one signed 75-edge graph and the zero bin of a
hard row is \(P\), the common difference-Radon transform gives independently

\[
 \sum_{a=1}^8q_L(a)=hT-P=21-P.                   \tag{15.743.3}
\]

Equating (15.743.2) and (15.743.3) forces

\[
 \boxed{P=4+k_L=5+e.}                            \tag{15.743.4}
\]

This is the essential cross-direction step.  Local rows with other values
of \(P\) can exist, but cannot be projections of the same signed graph.
Every opposite direction has \(Q=3\), and its off-zero-bin sum is
\(-hT-Q=-24\).

For a nine-set \(X\subset\mathbf F_{17}\), put

\[
 c_X(a)=|X\mathbin\triangle(X+a)|,
 \qquad 1\le a\le8.
\]

As \(X\) ranges over \(J(17,9)\), there are exactly 698 distinct vectors
\(c_X\).  Translation-summing a cut gives

\[
 \sum_t\operatorname{cut}_W(X+t)=c_X\mathbin\cdot q.
\]

For a hard cell \(A=9-\operatorname{cut}_W\), while an opposite cell has
\(B=-6-\operatorname{cut}_W/2\).  Thus every hard row of excess \(e\) lies
in the necessary integral relaxation

\[
 \begin{split}
 &\sum_aq_a=16-e,\qquad \|q\|_1\le70-e,\\
 &\sum_aa^2q_a\equiv\sum_aa^4q_a\equiv0\pmod {17},\\
 &c_X\mathbin\cdot q\le153\qquad(X\in J(17,9)),
 \end{split}                                      \tag{15.743.5}
\]

and every opposite row lies in

\[
 \begin{split}
 &\sum_aq_a=-24,\qquad \|q\|_1\le72,\\
 &\sum_aa^2q_a\equiv\sum_aa^4q_a\equiv0\pmod {17},\\
 &c_X\mathbin\cdot q\le-204\qquad(X\in J(17,9)).
 \end{split}                                      \tag{15.743.6}
\]

Exact one-worker CP-SAT models use the full \(\ell^1\) coordinate domains
and all 698 cut vectors.  They use neither a prior energy cap, the
opposite-entry alphabet, nor the lower cut bound.  The resulting sharp row
certificate is

\[
 \begin{array}{c|c}
 \text{row type}&\text{result}\\ \hline
 e=1&\text{infeasible},\\
 e=2&\max\|q\|_2^2=70,\\
 e=3&\max\|q\|_2^2=119,\\
 \text{opposite}&\max\|q\|_2^2=72.
 \end{array}                                      \tag{15.743.7}
\]

The last three bounds are attained, in natural distance order, by

\[
 (1,-2,5,3,-1,2,5,1),\quad
 (6,4,-1,-3,-2,-1,4,6),\quad
 (-3,-3,-3,-3,-3,-3,-3,-3),
\]

respectively.  Separate broad-domain models prove infeasibility at energies
71, 120, and 73.  No unverified full-row count or spectral truncation enters
the certificate.

The two quotient partitions containing \(e=1\) are already impossible by
(15.743.7).  In the remaining partition \((3)\), the exact
difference-Radon Gram identity gives the nonexact energy

\[
 1211+34C,\qquad
 C=\sum_\delta {m_\delta\choose2}\ge0.            \tag{15.743.8}
\]

Indeed the all-row identity is
\(17\cdot75+2\cdot21^2-2\sum_LP_L^2+34C\), and the eight exact
\((2,\ldots,2)\) rows contribute \(8\cdot32\).  But (15.743.7) bounds the
one \(e=3\) row and nine opposite rows by

\[
 119+9\cdot72=767<1211.
\]

Therefore every quotient partition is impossible and

\[
 \boxed{\text{the residual-(ii) row }p=17,k=74\text{ is empty}.}
\]

Proposition 15.743 is an **exhaustive finite certificate** at the eight-bin
aggregate level.  It does not close any \(p\ge17,t\ge4\) row or the generic
branch-B \(t=3\) range for \(p\equiv1\pmod4\), which now begins at \(p=29\).
Residual (ii), Type I, and the quadratic-minmax limit remain open.

Evidence: `src/e1_gmin_m4_prop15743.py`,
`tests/test_prop15743.py`,
`evidence/e1_gmin_m4_prop15743.json`, and
`evidence/NOTE_2026-08-31_p17_full_translated_cut_energy_close.md`.

## Proposition 15.744 — six roots close the p=13 fifth-shell residue u=3

At \(p=13,t=4\), the odd flip graph has 61 edges.  In the phase-one hard
type write

\[
 a_L=2u+14k_L,\qquad \sum_Lk_L=11-u.              \tag{15.744.1}
\]

The exact even-\(b\) phase-one floors, in the order
\(b=0,2,4,6,8,10,12\), are

\[
 (26,12,26,24,26,20,14).                          \tag{15.744.2}
\]

For \(b=10\), pass to the three-point complement.  Phase one becomes the
phase-zero \(b=3\) parity problem.  Proposition 15.652's exact quadrature
has coefficient polynomial \(r^2-4r+4=(2-r)^2\), contact nodes
\(1,2,3\), and strictly positive weights
\(15/26,3/13,5/26\).  At equality, positivity forces the original
integer quadratic to equal the parity minimum on every point in those
three layers.  To extend across the omitted \(r=0\) layer, restrict the 78
pair monomials to all 1,596 points with \(r=1,2,3\).  Their evaluation
matrix has rank 78 modulo 101.  Since pair monomials span every
degree-at-most-two function on \(J(13,7)\), this proves that the exact cell
is pointwise \(A=(2-r)^2\), not merely that its stabilizer average has that
form.

The \(b=10\) cell at mean 22 is not an ordinary nonnegative lift of this
baseline: \(B=(A-(2-r)^2)/2\) need only satisfy \(B\ge-2\) on \(r=0\),
although \(B\ge0\) on \(r=1,2,3\).  The executable certificate therefore
uses Proposition 15.738's 1,638 independent third-difference identities,
these punctured lower bounds, and \(\sum_XB(X)=66\), equivalently
\(4p\,\mathbb E B=2\).  The resulting deterministic 1,716-variable model
is infeasible.  Thus this floor-plus-two cell is excluded without invoking
the globally nonnegative lift theorem.

For the ordinary \(b=2,12\) baselines, the sharp cost of a nonzero
nonnegative integral lift is ten.  For
\(u\le4\), a quotient-zero row would have mean at most eight, below every
floor, so every quotient is positive and at least \(3+u\) directions have
\(k_L=1\).  Applying (15.744.2) at their low means leaves the exact
\(b=12\) row for \(u=0\), nothing for \(u=1,2\), the exact \(b=10\)
complement triple for \(u=3\), and the sharp mass-ten \(b=2\) lift for
\(u=4\).  At \(u=5\), a forced quotient-zero row has mean ten and is
impossible; at \(u=6\), two such rows can attain the exact \(b=2\) floor
twelve.  Thus the exact residue sieve is

\[
 \boxed{u\in\{0,3,4,6\}}.                          \tag{15.744.3}
\]

Now take \(u=3\).  Equation (15.744.1) forces the hard quotient profile
\(1^6 2\).  The six low rows are exact complement triples of mean 20,
whose signed target has coefficient offset two.  If \(P\) is their common
parallel count, \(R\) is the elevated count, \(h\) is the hard sign, and
\(T\) is the global signed edge total, coefficient normalization gives

\[
 20=14P-hT-39,\qquad 34=14R-hT-39.                 \tag{15.744.4}
\]

Hence \(hT=14P-59\), \(R=P+1\), and the offset congruence is
\(P\equiv2\pmod6\).  Since the hard rows use \(6P+R=7P+1\le61\) edges,

\[
 P\in\{2,8\}.                                     \tag{15.744.5}
\]

For \(P=2\), the seven opposite counts sum to 46 and their means are
\(14(Q-5)\).  A \(Q=5\) zero cell has coefficient offset three and would
require \(Q\equiv3\pmod6\), so every opposite count is at least six; at
least one is \(Q=6\), of mean 14.  For \(P=8\), the opposite counts sum to
four and their means are \(14(Q+1)\), so at least one is \(Q=0\), again of
mean 14.

It remains to classify this forced phase-zero mass-14 cell.  Proposition
15.738's Boolean support classification is independent of \(|H|\), but its
preliminary height-four models used the old bound
\(\sum|W_{st}|\le59-Q\).  They cannot be imported at this layer.  Rebuild
the two necessary coefficient models with

\[
 \sum|W_{st}|\le61-Q,\qquad Q\in\{0,6\}.            \tag{15.744.6}
\]

Exact deterministic one-worker CP-SAT returns `INFEASIBLE` in both cases.
The live height dichotomy is therefore reduced to the Boolean case.  The
edge-independent exhaustive support-462 catalog has offsets six, four, and
four, so at \(Q=0,6\) only

\[
 B=x_i x_j                                             \tag{15.744.7}
\]

survives.

For even \(d\), put

\[
 M_d(L)=\sum_{\{x,y\}\in H}\chi(x-y)(Lx-Ly)^d.
\]

In either hard gauge, an exact complement triple satisfies
\(2S_4=S_2^2\).  Thus the homogeneous binary quartic

\[
 G=2hM_4-M_2^2                                      \tag{15.744.8}
\]

vanishes in six distinct projective directions.  It is identically zero.
The opposite selected-pair cell (15.744.7), whose sign is \(-h\), instead
gives

\[
 G=-3(i-j)^4\ne0\quad\hbox{in }\mathbb F_{13}.      \tag{15.744.9}
\]

This contradiction proves

\[
 \boxed{p=13,t=4,u=3\text{ is empty}.}
\]

This is a **proved branch theorem**, with exhaustive finite certificates for
the contact-layer restriction, the punctured lift, and the changed
height-four premise.  It does not close the other residues or residual (ii).

Evidence: `src/e1_gmin_m4_prop15744.py`,
`tests/test_prop15744.py`,
`evidence/e1_gmin_m4_prop15744.json`, and
`evidence/NOTE_2026-08-31_p13_t4_u3_quartic_close.md`.

## Proposition 15.745 — collision-one energy closes p=13,t=4,u=0

Retain \(p=13,t=4,u=0\).  The seven hard quotients satisfy
\(k_L\ge1\) and \(\sum_Lk_L=11\), so the excesses
\(e_L=k_L-1\) have one of the five partitions

\[
 (4),\ (3,1),\ (2,2),\ (2,1,1),\ (1,1,1,1).       \tag{15.745.1}
\]

Normalize before inserting the exact-star distance row.  If \(P_0\) is an
exact literal star's parallel count, its local and common off-bin sums give

\[
 13(P_0-3)-14=hT-P_0,qquad hT=14P_0-53.           \tag{15.745.2}
\]

For a general hard quotient \(k_L\), the same comparison gives
\(P_L=P_0+k_L-1\).  Summing the hard parallel counts and using the 61-edge
bound gives \(P_0\le8\).  The isolated-chart literal congruence
\(P_0\equiv5\pmod6\) now forces

\[
 P_0=5,qquad hT=17,qquad P_L=4+k_L=5+e_L.        \tag{15.745.3}
\]

The hard rows use 39 parallel edges, so the seven opposite counts sum to
22.  Their scaled mean is \(14Q-22\).  Counts below two give negative mean,
while \(Q=2\) gives mean six, below both the least nonzero phase-zero floor
and the sharp integral-lift floor.  Hence the opposite profile is exactly

\[
 (3,3,3,3,3,3,4).                                 \tag{15.745.4}
\]

An exact literal row is \((2,2,2,2,2,2)\), of energy 24.  Three distinct
exact directions force the binary quadratic \(M_2\) to vanish; five force
the binary quartic \(M_4\) to vanish as well.  The executable certificate
checks all 364 three-direction quadratic evaluation matrices and all 2,002
five-direction quartic matrices over \(\mathbb F_{13}\).

Let \(q=(q_1,\ldots,q_6)\) be a nonzero distance row.  Use all 74 exact
translated seven-set cut vectors.  The necessary integral relaxations have

\[
\begin{array}{lll}
 \text{hard excess }e:&\sum q_a=12-e,&
 \|q\|_1\le56-e,\quad c_Xq\le91,\\
 \text{opposite }Q:&\sum q_a=-17-Q,&
 \|q\|_1\le61-Q,\quad c_Xq\le-130,
\end{array}                                        \tag{15.745.5}
\]

together with the forced moment congruences.  Broad-domain deterministic
one-worker models, independently replayed by a second exact encoding at the
next energy, give

\[
\begin{array}{c|c|c}
 \text{row}&\text{moments}&\text{result}\\ \hline
 Q=3&M_2=M_4=0&\text{infeasible},\\
 e=1&M_2=0&\max\|q\|_2^2=31,\\
 e=2&M_2=0&\max\|q\|_2^2=96,\\
 Q=3&M_2=0&\max\|q\|_2^2=76,\\
 Q=4&M_2=0&\max\|q\|_2^2=111.
\end{array}                                        \tag{15.745.6}
\]

Thus the first three partitions in (15.745.1), which have at least five
exact stars, are empty.  The partition \((1,1,1,1)\) has three exact stars;
its nonexact row energy is at most

\[
 4\cdot31+6\cdot76+111=691,
\]

whereas the exact difference-Radon identity is \(721+26C\ge721\), with
\(C=\sum_\delta\binom{m_\delta}{2}\).

Only \((2,1,1)\) remains.  Its broad row upper bound and exact Radon value
are

\[
 96+2\cdot31+6\cdot76+111=725,
 \qquad 693+26C.                                   \tag{15.745.7}
\]

Consequently \(C\le1\).  Independently, the unique \(e=2\) direction has
\(P=7\) parallel edges in only six unoriented displacement classes, so
\(C\ge1\).  Equality forces its multiplicity profile to be
\((2,1,1,1,1,1)\), places the sole duplicate in that row's zero bin, and
makes every transverse multiplicity zero or one.

For each nonzero bin of this hard direction, the difference-Radon incidence
contains exactly one class from each of the other thirteen projective
directions.  Since \(\chi(-1)=1\), the unoriented signs are well defined;
after hard normalization they split into six positive and seven negative
classes.  The certificate checks all \(14\cdot6=84\) such bins directly in
\(\mathbb F_{13}^2\).  Hence

\[
 -7\le q_a\le6.                                    \tag{15.745.8}
\]

Adding only (15.745.8) to the \(e=2\) relaxation lowers its sharp energy to
66; a separate model makes energy at least 67 infeasible.  Therefore

\[
 66+2\cdot31+6\cdot76+111=695<693+26=719,          \tag{15.745.9}
\]

the final contradiction.  Thus

\[
 \boxed{p=13,t=4,u=0\text{ is empty}.}
\]

Proposition 15.745 is an **exhaustive finite aggregate certificate**.  With
Proposition 15.744, the exact remaining \(p=13,k=60\) residues are
\(u\in\{4,6\}\).  The full row, residual (ii), Type I, and the limit remain
open.

Evidence: `src/e1_gmin_m4_prop15745.py`,
`tests/test_prop15745.py`,
`evidence/e1_gmin_m4_prop15745.json`, and
`evidence/NOTE_2026-08-31_p13_t4_u0_close.md`.

## Proposition 15.746 — sharp mass-ten catalog and the p=13,t=4,u=4 reduction

Let \(\Omega=J(13,7)\).  A hard cell in the surviving residue \(u=4\)
has \(b=2\) and scaled mean 22.  The all-positive quadrature behind the
\(b=2\) floor first forces its equality baseline pointwise, not merely after
stabilizer averaging.  Subtracting that XNOR baseline gives a nonnegative
integral degree-at-most-two lift \(B\) with

\[
 4p\,\mathbb E_\Omega B=10.                       \tag{15.746.1}
\]

Only after this pointwise step does Proposition 15.688 apply.  Its sharp
mass-ten equality gives \(0\le B\le1\), so \(B\) is Boolean and

\[
 |\operatorname{supp}B|={10\over52}{13\choose7}=330. \tag{15.746.2}
\]

Proposition 15.738 already certifies that the quadratic evaluation space on
\(\Omega\) has rank 78 and is exactly the nullspace of 1,638 independent
third-difference identities.  Introduce one Boolean value \(f_X\) for each
of the 1,716 points of \(\Omega\), impose those identities and
\(\sum_Xf_X=330\), and anchor one value \(f_{X_0}=1\).  The full candidate
catalog consists of

\[
 \begin{array}{c|c|c}
 \text{family}&B(x)&\text{count}\\ \hline
 \text{omitted pair}&(1-x_i)(1-x_j)&{13\choose2}=78,\\
 \text{all-equal triple}&
 1-x_i-x_j-x_k+x_ix_j+x_ix_k+x_jx_k&{13\choose3}=286.
 \end{array}                                      \tag{15.746.3}
\]

Every displayed function has support 330 and satisfies all 1,638
identities.  The catalog is invariant under the twelve adjacent coordinate
transpositions, and the induced coordinate action is transitive on
\(\Omega\), so the anchor loses no nonempty support.  Exactly 15 omitted
pairs and 55 all-equal triples contain \(X_0\).  Excluding these 70 anchored
supports gives an exact model with 1,716 Boolean variables and

\[
 1638+1+1+70=1710                                  \tag{15.746.4}
\]

constraints.  Exact CP-SAT returns `INFEASIBLE`; a deterministic one-worker
replay gives the same result.  The full and anchored catalog SHA-256 digests
are, respectively,

```text
4edf1fe1b9c73f05598b667dba121f064807c68421a4df2c8db7090a3e3ff35f
84ce6099dcca66f7cc2792dc60bcbb378672f2e9cac2b19e02812f2f20563c7a
```

Thus (15.746.3) exhausts every Boolean quadratic of support 330 on
\(J(13,7)\).

Return to the seven hard directions.  In \(z_i=2x_i-1\) coordinates the
pointwise target is

\[
 4+z_i z_j+4B.                                    \tag{15.746.5}
\]

For an omitted-pair lift,
\(4B=1-z_a-z_b+z_az_b\), so its total coefficient offset is three.  For an
all-equal triple,
\(4B=1+z_az_b+z_az_c+z_bz_c\), so the offset is five.  The common signed
edge total makes the hard parallel count \(P\) common.  The isolated-chart
coefficient congruence and \(7P\le61\) therefore force

\[
 \begin{array}{c|c|c}
 \text{hard family}&P&hT\\ \hline
 \text{omitted pair}&3&-19,\\
 \text{all-equal triple}&5&9,
 \end{array}
 \qquad hT=14P-61.                                \tag{15.746.6}
\]

In particular the two hard families cannot mix.  If \(Q_L\) is an opposite
parallel count, then

\[
 \sum_LQ_L=61-7P,
 \qquad e_L=P+Q_L-8\ge0,
 \qquad a_L=12+14e_L.                             \tag{15.746.7}
\]

Since \(\sum_Le_L=5\), at least two opposite directions have \(e_L=0\)
and scaled mean 12.  Their minimum parallel count is \(Q=5\) in the
omitted-pair branch and \(Q=3\) in the all-equal-triple branch.

A phase-zero mean-12 cell has exactly the following still-live alternatives:
the exact \(b=12\) literal \(A=1-x_j\), with target
\(3+2A=4-z_j\), or a \(b=0\) cell \(A=2C\), where \(C\) is a nonzero
nonnegative integral quadratic satisfying \(4p\mathbb E C=12\).
Proposition 15.688 makes the latter's height exactly one or four; height one
is Boolean of support 396.  In the omitted-pair branch, a literal at the
minimum \(Q=5\) would require \(6\mid Q-3\), which is false.  Hence every
forced mean-12 cell there is a \(b=0\) mass-12 lift, of height one/support
396 or height four.  At \(P=5,Q=3\), the literal congruence is compatible,
so the all-equal-triple branch retains the literal-or-lift dichotomy.

There is one further common-graph consequence in the omitted-pair branch.
For \(r=1,2,3\), its normalized hard moments have the two-power-sum form
\(hM_{2r}=\alpha^{2r}+\beta^{2r}\).  The degree-six Newton identity gives

\[
 F_6=2hM_6+hM_2^3-3M_2M_4=0.                    \tag{15.746.8}
\]

The complete-graph and star gauge sums vanish in degrees two, four, and six;
the executable certificate also checks all \(78^2=6084\) possible overlaps
between the baseline and omitted-pair coordinates.  Thus the seven hard
directions are seven projective roots of the homogeneous binary sextic
\(F_6\), forcing it to vanish identically.  Every forced
\(P=3,Q=5,b=0\) mass-12 opposite cell must consequently satisfy \(F_6=0\).
For the opposite-cell sign normalization
\(N'_{2r}=(-h)M_{2r}\), this condition must be encoded as

\[
 2N'_6+(N'_2)^3+3N'_2N'_4=0,                    \tag{15.746.9}
\]

not by copying the hard-sign Newton formula verbatim.  For the \(P=5\)
family, the 22,308 baseline-pair/triple patterns give ranks (1,2,3) to the
weighted feature lists
(N_2), ((N_4,N_2^2)), and
((N_6,N_2N_4,N_2^3)).  Thus there is no nonzero universal
weighted-homogeneous polynomial identity in these even moments through
degree six; this is not a claim about every conceivable invariant.

Proposition 15.746 is an **exhaustive finite equality classification and
proved open reduction**.  It does not exclude either \(u=4\) family, so the
exact \(p=13,k=60\) remainder stays \(u\in\{4,6\}\).  The next narrow gate
is one \(Q=5,b=0\) mass-12 cell under (15.746.9), split into height one
(support 396) and height four; a broad support census is not the gate.
Residual (ii), Type I, and the quadratic-minmax limit remain open.

Evidence: `src/e1_gmin_m4_prop15746.py`,
`tests/test_prop15746.py`,
`evidence/e1_gmin_m4_prop15746.json`,
`evidence/p13_support330_boolean_classifier.json`, and
`evidence/NOTE_2026-08-31_p13_support330_u4_reduction.md`.

## Proposition 15.747 — the mass-twelve cut obstruction closes the P=3 branch

Continue in the two `u=4` branches isolated by Proposition 15.746.  For a
phase-zero `b=0` cell write `A=2C`, and let `W` be its normalized integral
coefficient graph.  At parallel count `Q`, averaging the coefficient identity
gives

\[
 S:=\sum_eW_e=13(Q-3)-12=13Q-51,
 \qquad \operatorname{cut}_W(X)=7Q-27-2C(X).       \tag{15.747.1}
\]

Suppose first that `C` has height one.  It is Boolean of support 396.  Take a
uniform six-set, whose cut is the same as that of its complementary point of
`J(13,7)`.  One edge is cut with probability `7/13`; two distinct adjacent
edges are both cut with probability `7/26`, and two disjoint edges with
probability `42/143`.  If

\[
 E_2=\sum_eW_e^2,\qquad
 D_2=\sum_v\left(\sum_{e\ni v}W_e\right)^2,
\]

then exact expansion yields

\[
 \mathbb E[\operatorname{cut}_W^2]
 ={ -7D_2+84E_2+84S^2\over286}.                  \tag{15.747.2}
\]

On the other hand, (15.747.1) has 396 occurrences of `7Q-29` and 1,320
occurrences of `7Q-27`.  Equating the second moments gives

\[
 -7D_2+84E_2+182Q^2-1428Q+2598=0.               \tag{15.747.3}
\]

The last three terms are congruent to one modulo seven for every integral
`Q`, while the first two are divisible by seven.  Thus no Boolean mass-twelve
lift exists at any `Q`.  In particular, the required second moments at
`Q=3,5` would be `552/13` and `748/13`, producing the impossible residuals
`-48` and `8` in (15.747.3).  This argument needs neither the edge `l1`
bound, row parity, third differences, nor the sextic.

It remains to exclude height four at `Q=3,5`.  Here the field moments can be
dropped, restoring the full `S_13` symmetry; a point where `C=4` may therefore
be moved to the first seven-set.  The exact necessary models use 78 integral
weights, 78 absolute values, and 13 row halves.  They impose

\[
 \sum W=-12, \sum|W|\le58, -14\le\operatorname{cut}_W\le-6
 \quad(Q=3),                                      \tag{15.747.4}
\]

or

\[
 \sum W=14, \sum|W|\le56, 0\le\operatorname{cut}_W\le8
 \quad(Q=5),                                      \tag{15.747.5}
\]

together with even coefficient rows and the appropriate anchored lower cut.
Each projected model has 169 variables and 3,526 constraints.  Deterministic
one-worker CP-SAT returns `INFEASIBLE` for both, with model hashes

```text
Q=3  e8404a5684e033b73750b1f36a338aa13038861d6dbfc614cc99b6f0666423d9
Q=5  8f992368fac869f29c23e6ecd20400228c2c10d5bda4d1001b291242dd6e3941
```

Thus no `b=0` mass-twelve lift survives at either minimum count.  The
omitted-pair `P=3` branch forces at least two `Q=5` cells and does not permit
the literal, so it is empty.  In the all-equal-triple `P=5` branch, every
minimum `Q=3` cell is consequently the exact `b=12` literal `A=1-x_j`.

Proposition 15.747 is a **proved branch exclusion with two exhaustive finite
necessary-relaxation certificates**.  It closes the `P=3` half of `u=4`, not
the full residue.

## Proposition 15.748 — literal-root interpolation leaves only excess 1^5

Let `z` be the number of minimum `Q=3` opposite directions in the remaining
`P=5` branch.  Proposition 15.747 makes all of them exact literals, hence
common projective roots of the homogeneous binary forms `M_2,M_4,M_6`.
The seven opposite excesses are nonnegative and sum to five, so `z>=2`.

The complete local moment alphabet of a hard baseline pair plus an all-equal
triple has 69 triples `(N_2,N_4,N_6)`.  Its fourth moment is never zero.  If
`z>=5`, five roots force the binary quartic `M_4` to vanish identically, an
immediate contradiction.  If `z>=3`, the quadratic `M_2` vanishes.  Writing
`w` for the baseline-pair difference and `T_d` for the triangle moments,
`N_2=w^2+T_2=0` and `2T_4=T_2^2` give

\[
 N_4=w^4+T_4={3\over2}w^4=8w^4
 \in\{7,8,11\}.                                  \tag{15.748.1}
\]

For each of the two possible hard signs, exact interpolation now checks all
choices of literal roots.  At `z=4`, `M_4` is a scalar times their four root
factors; at `z=3`, it is their cubic root product times an arbitrary linear
form.  Neither case has a value vector in the alphabet (zero survivors for
both signs).

At `z=2`, write the common quadratic root product as `R_2`.  The exhaustive
parameterization is

\[
 M_2=cR_2,\qquad M_4=R_2Q_2,\qquad M_6=R_2Q_4.    \tag{15.748.2}
\]

After the first two moment filters, 1,554 `(M_2,M_4)` candidates remain per
sign.  Exact degree-four evaluation-code membership checks 2,688 allowed
`N_6` vectors per sign and leaves 336 distinct moment-level survivors per
sign.  These are necessary moment data, not common 61-edge graphs.  The raw
enumeration payload hash is
`894c087d4acae7ff0722ba236b1fac494984b9b331431e6117b2edbde0afbbec`.

Therefore `z=2`.  The other five opposite excesses are positive integers
with sum five, so all are one:

\[
 \boxed{(e_L:e_L>0)=(1,1,1,1,1).}                \tag{15.748.3}
\]

Proposition 15.748 is an **exhaustive finite interpolation certificate and
proved open reduction**.  The `P=5` branch and `u=4` remain open.  The narrow
next gate is to couple the 336 moment survivors per sign to the five
excess-one `Q=4` cells and one common 61-edge graph.

Evidence: `src/e1_gmin_m4_prop15747.py`,
`src/e1_gmin_m4_prop15748.py`, `scripts/p13_p5_literal_interpolation.py`,
`tests/test_prop15747.py`, `tests/test_prop15748.py`,
`evidence/e1_gmin_m4_prop15747.json`, and
`evidence/e1_gmin_m4_prop15748.json`, with the combined derivation in
`evidence/NOTE_2026-08-31_p13_u4_mass12_literal_interpolation.md`.

## Proposition 15.749 — translated-cut moments close p=13,t=4,u=4

Continue with a `z=2` record from Proposition 15.748.  Its two minimum
opposite directions are literal cells, and its other five opposite
directions have `Q=4`.  The common normalization is `hT=9`.  If

\[
 q_a=(-h)\sum_{|L(\delta)|=a}\chi(\delta)m_\delta
 \qquad(1\le a\le6)
\]

is the distance row of one such `Q=4` direction, then

\[
 \sum_aq_a=-13,\qquad \sum_a|q_a|\le57.             \tag{15.749.1}
\]

Moreover its phase-zero cell obeys

\[
 (-h)S_H=Q+\sum_aq_a-2\operatorname{cut}_W
        =-9-2\operatorname{cut}_W=3+2A.
\]

Thus `A=-6-cut_W>=0`.  Summing over all thirteen translates of a
seven-set gives

\[
 c\mathbin\cdot q\le-78                              \tag{15.749.2}
\]

for each of Proposition 15.740's 74 exact translated-cut vectors.

There is a short exact bound before any list recovery.  In the canonical
cut-vector order,

\[
\begin{aligned}
 e_1&={19\over9}{\bf1}-{c_0\over18}-{c_6\over6}-{c_{34}\over18},\\
-e_1&={29\over15}{\bf1}-{c_{63}\over15}-{c_{69}\over30}
                   -{c_{71}\over6}-{c_{73}\over30}.
\end{aligned}                                        \tag{15.749.3}
\]

Equations (15.749.1)--(15.749.3) give
`-52/9<=q_1<=26/15`.  Multiplication by
`F_13^*/{+-1}` preserves the complete cut catalog and acts transitively on
the six distance bins, so integrality gives

\[
                 -5\le q_a\le1\quad(1\le a\le6).     \tag{15.749.4}
\]

Exact enumeration of this six-variable box under (15.749.1)--(15.749.2)
has 522 rows and 492 distinct moment triples

\[
 (N_2,N_4,N_6)=\left(\sum a^2q_a,\sum a^4q_a,
                            \sum a^6q_a\right)\pmod {13}.       \tag{15.749.5}
\]

For each of the 336 `z=2` records per hard sign, reconstruct the full forms
`M2=cR2`, `M4=R2Q2`, and `M6=R2Q4`; the stored seven hard values recover the
quartic quotient exactly.  Opposite rows use the negative of this hard
normalization.  For either hard sign the five nonroot evaluations form the
same 48-element alphabet.  Its intersection with the 492 triples in
(15.749.5) is exactly

\[
\begin{split}
 &(1,0,3),(2,0,1),(3,0,3),(4,0,10),(5,0,1),(6,0,1),\\
 &(7,0,12),(8,0,12),(9,0,3),(10,0,10),(11,0,12),(12,0,10).
\end{split}                                          \tag{15.749.6}
\]

In particular, every admissible `Q=4` direction is a root of the common
quartic `M4`.  The two literal directions are roots as well.  Seven distinct
projective roots force this binary quartic to vanish identically, whereas
the hard moment alphabet of Proposition 15.748 has no zero fourth moment.
This contradiction closes the `P=5` branch.  Proposition 15.747 already
closes `P=3`, hence

\[
             \boxed{p=13,t=4,u=4\text{ is empty}.}    \tag{15.749.7}
\]

This is a **proved branch theorem with an exact aggregate certificate**.
The exact `p=13,k=60` remainder is now `u=6`; residual (ii) globally,
multi-level Type I, and the quadratic-minmax limit remain open.

Evidence: `src/e1_gmin_m4_prop15749.py`, `tests/test_prop15749.py`,
`evidence/e1_gmin_m4_prop15749.json`, and
`evidence/NOTE_2026-09-01_p13_u4_translated_cut_moment_close.md`.

## Proposition 15.750 — isolated-chart parity halving closes multi-level Type I

Let \(p\ge5\) be prime, let \(e\) be a distinguished edge, and suppose that
a Boolean edge set \(G\) is a multi-level Type-I bad case:

\[
 |G|=3p-2,\quad e\notin G,\quad
 S_G=3-2f_e\quad\hbox{on Max+},                    \tag{15.750.1}
\]

while on Max−,

\[
 S_G\le-1,qquad S_G\le-3f_e.                     \tag{15.750.2}
\]

Put \(H=G\cup\{e\}\) and let \(W=G+2e\) be an edge multiset. Then

\[
 |H|=3p-1,\quad |W|=3p,\quad S_W=3\text{ on Max+},
 \quad S_H\le-2\text{ on Max−}.                \tag{15.750.3}
\]

Assume first that \(p\ge11\). At least

\[
 p^2+1-2(3p-1)=p^2-6p+3>0                         \tag{15.750.4}
\]

vertices are isolated in \(H\). Proposition 15.721 transports one to
infinity. Its signed relative-flip identity is linear in nonnegative integer
edge multiplicities, so the same permutation transports \(W=H+e\), including
the doubled edge. Both eigenshell conditions are preserved and every edge is
finite in the resulting affine chart.

Write \(q=(p-1)/2\) and \(m=(p+1)/2\). Fix a square direction, let \(P\)
be the total \(W\)-multiplicity parallel to it, and let \(K_{st}\) be the
signed edge sum between fibres \(s,t\). The affine Max+ cylinders give

\[
 P+\sum_{s<t}K_{st}z_sz_t=3
 \quad\text{whenever }z_s\in\{\pm1\},\ \sum_sz_s=1. \tag{15.750.5}
\]

The Johnson-slice swap lemma makes all \(K_{st}\) equal. Indeed, swap two
coordinates \(a,b\) while fixing a middle subset among the remaining
coordinates. The difference is a linear form with coefficients
\(K_{ak}-K_{bk}\); swapping two remaining coordinates makes these
coefficients equal, and their dot product with a sign vector of sum one
makes the common value zero. Varying \(a,b\) proves \(K_{st}=\kappa\).
Since \(\sum_{s<t}z_sz_t=-q\),

\[
 P=3+q\kappa,\qquad
 3p-P\ge {p\choose2}|\kappa|=pq|\kappa|.          \tag{15.750.6}
\]

For \(p\ge11\), (15.750.6) forces \(\kappa=0\), hence \(P=3\). Summing over
the \(m\) square directions gives positive and negative \(W\)-multiplicities
\(3m\) and \(3q\), and signed total three.

Let \(c=C_e\in\{\pm1\}\). Removing one copy of \(e\) gives
\(\tau:=\sum C_{ab}H_{ab}=3-c\). In a nonsquare direction let \(P_d\)
count parallel \(H\)-edges and define

\[
 T_d(z)={-S_H(z)-2\over2}\ge0.                    \tag{15.750.7}
\]

This is an integer-valued quadratic on \(J(p,m)\). Since
\(\mathbb E[z_sz_t]=-1/p\),

\[
 a_d:=2p\mathbb E[T_d]=(p+1)P_d-2p+\tau.          \tag{15.750.8}
\]

Nonnegativity forces every \(P_d\ge2\). Their sum is \(3q\) if \(c=1\)
and \(3q-1\) if \(c=-1\), both strictly below \(3m\). Thus some direction
has \(P_d=2\), and there

\[
 a_d=4\quad(c=1),\qquad a_d=6\quad(c=-1).         \tag{15.750.9}
\]

The product of the edge features of \(H\) is a constant times the parity of
the fibres meeting the odd-degree boundary of \(H\). Because

\[
 (-1)^{(|H|-S_H)/2}=\prod_{g\in H}f_g,
\]

\(T_d\bmod2\) is an affine fibre parity. For a nonconstant parity supported
on \(b\) fibres, its bias on \(J(p,m)\) is \(K_m(b)/{p\choose m}\). The
central Krawtchouk recurrence

\[
 (p-b)K_m(b+1)=-K_m(b)-bK_m(b-1),\qquad
 K_m(1)=K_m(2)=-{1\over p}{p\choose m},            \tag{15.750.10}
\]

and induction through \(b\le(p-1)/2\), followed by complement symmetry,
give

\[
 |K_m(b)|\le {1\over p}{p\choose m}qquad(1\le b\le p-1). \tag{15.750.11}
\]

Hence either value of a nonconstant affine parity has probability at least
\((p-1)/(2p)\). If \(T_d\) were odd on that class, nonnegativity would give
\(2p\mathbb E[T_d]\ge p-1\); constant odd parity would cost \(2p\).
But (15.750.9) is at most six and \(p-1\ge10\). Thus \(T_d\) is everywhere
even. Then \(B_d=T_d/2\) is a nonzero nonnegative integer-valued quadratic,
so Proposition 15.688 yields

\[
 4p\mathbb E[B_d]\ge p-3.                         \tag{15.750.12}
\]

Yet \(4p\mathbb E[B_d]=2p\mathbb E[T_d]=a_d\in\{4,6\}<p-3\), a
contradiction.

For \(p=5,7\), tracked positive integer Farkas multipliers give exact
identities \(\sum_r\lambda_rA_r=0\) and
\(\sum_r\lambda_rb_r<0\) for the canonical edge. The verifier regenerates
the Paley matrix, checks every stored Boolean eigenvector, and evaluates
both identities with arbitrary-precision integers; it invokes neither an
optimizer nor an eigenshell cache. Finally, \(\operatorname{PSL}(2,p^2)\)
is 2-transitive: determinant-one translate/invert maps send any ordered edge
to \((\infty,0)\), and the signed lift preserves (15.750.1)--(15.750.2).
Therefore

\[
 \boxed{\text{the multi-level Type-I bad case is empty for every prime }
 p\ge5.}                                           \tag{15.750.13}
\]

This closes only Type I. Residual (ii), E1, \(L=1/2\), and the original
convergence question remain open. Evidence:
`src/e1_gmin_m4_prop15750.py`, `src/e1_type_i_small_prime_exact.py`,
`evidence/e1_type_i_badcase_farkas_p5.json`,
`evidence/e1_type_i_badcase_farkas_p7.json`, and
`evidence/e1_gmin_m4_prop15750.json`.

## Proposition 15.751 — influence rigidity closes the fourth residual shell

In generic branch B at (t=3), every opposite direction supplies a nonzero
nonnegative integer-valued quadratic (B) on (J(p,(p+1)/2)) with

\[
                         4p\mathbb EB=p+7.          \tag{15.751.1}
\]

First, any nonnegative integral quadratic (g) on a Boolean cube with
(mathbb Eg=1/2) has (max g\le3). A nonzero quadratic has support density
at least (1/4), restriction means lie in ((1/4)\mathbb Z), and a minimal
counterexample forces every facet away from its maximum to be Boolean. Zero
third differences then force a putative maximum four to have singleton,
pair, and triple values (1,0,1); its value on a four-set is again four,
a contradiction. The bound is sharp for (3-2s+{s\choose2}), whose layer
values are (3,1,0,0,1).

Let (H=\max B). Paired cubes through a maximizer have average mean

\[
 TB(X)={4H+p+7\over4(p+1)}.                         \tag{15.751.2}
\]

If (H\ge2), each cube has mean at least (1/2), so (H\ge(p-5)/4).
The exact stabilizer bound gives

\[
 H\le{(p+7)(p+3)\over4(p-1)},\qquad
 TB(X)\le{p+7\over2(p-1)}<{3\over4}.               \tag{15.751.3}
\]

Some cube therefore has mean (1/2) and (H\le3), contradicting the lower
bound for (p\ge29).

For (H=1), complement the slice and write
(f:J(p,(p-1)/2)\to\{0,1\}), with (mu=(p+7)/(4p)). Put

\[
 I_{ij}={1\over4}\Pr[f(X)\ne f(X^{(ij)})].          \tag{15.751.4}
\]

A relevant derivative is a nonzero ({-1,0,1})-valued affine function on
(J(p-2,(p-3)/2)). Sorting its coefficients shows their total integral
deviation from the median is at most two, so its support density is at least
((p-3)/(2(p-2))), and

\[
 I_{ij}\ge{(p+1)(p-3)\over16p(p-2)}.               \tag{15.751.5}
\]

Zero-influence transpositions form equivalence classes. If the complement
of the largest class has size (L), there are at least (pL/2) relevant
pairs. The correctly normalized Johnson Laplacian identity is

\[
 \sum_{i<j}I_{ij}
 ={1\over2}\sum_{e=1}^2e(p+1-e)\|f_{=e}\|_2^2
 \le(p-1)\mu(1-\mu).                               \tag{15.751.6}
\]

Hence

\[
 L\le {2(p-1)(p-2)(p+7)(3p-7)\over p^2(p+1)(p-3)}<7. \tag{15.751.7}
\]

After (p=x+29), the cleared gap is
(x^4+92x^3+3107x^2+45296x+237300>0); thus (L\le6).
Symmetrizing over the largest class extends (f) to a Boolean quadratic on
the full (L)-cube. Cube derivative support and Parseval reduce it to four
actual coordinates.

The exact four-bit Möbius catalog contains 222 quadratic Boolean tables and
fourteen profiles. Their possible complementary-slice densities are

\[
 0,1,{p-3\over4p},{p+1\over4p},{p-1\over2p},{p+1\over2p},
 {3p-1\over4p},{3(p+1)\over4p},                    \tag{15.751.8}
\]

none equal to ((p+7)/(4p)). CUDA on a V100 and Orin, OpenCL on AMD gfx1201
and Intel Arc A380, and exact CPU replay agree on all 222 tables, fourteen
profiles, and SHA-256
`63c9daf2b117b540a5199b1b007cb4e6997ba01704fbc6017efaaa9735859396`.
This is a fixed four-bit certificate, not a prime census.

Thus generic branch B is empty for every live (p\ge29). With the existing
(p=13,17) certificates and branch-A/C arithmetic,

\[
 \boxed{k=4p+6\text{ is impossible for every prime }p\ge13.} \tag{15.751.9}
\]

Later layers and residual (ii) globally remain open. Evidence:
`src/e1_gmin_m4_prop15751.py`, `tests/test_prop15751.py`,
`scripts/boolean_cube_degree2_gpu_audit.py`, and
`evidence/NOTE_2026-09-01_GENERIC_T3_INFLUENCE_CLOSE.md`.

## Proposition 15.752 — influence rigidity closes a contiguous residual band

Let \(p\ge23\) be prime. There is no nonzero nonnegative integer-valued
quadratic \(B\) on \(J(p,(p+1)/2)\) satisfying

\[
                         4p\,\mathbb EB=p+9.          \tag{15.752.1}
\]

Indeed, put \(H=\max B\). The paired-cube operator through a maximizing
middle set gives

\[
 TB(X)={4H+p+9\over4(p+1)}.                         \tag{15.752.2}
\]

If \(H\ge2\), every paired-cube restriction is a nonzero integral quadratic.
Its mean lies in \((1/4)\mathbb Z\), and mean \(1/4\) would make it Boolean,
contrary to the value \(H\ge2\) at \(X\). Hence every restriction has mean
at least \(1/2\), and

\[
                         H\ge {p-7\over4}.           \tag{15.752.3}
\]

The exact stabilizer inequalities from Proposition 15.688 give

\[
 H\le {p+9\over4}\quad(p\equiv3\pmod4),\qquad
 H\le{(p+9)(p+3)\over4(p-1)}\quad(p\equiv1\pmod4). \tag{15.752.4}
\]

Substitution in (15.752.2) yields respectively

\[
 TB(X)\le {p+9\over2(p+1)}<\frac34,
 \qquad
 TB(X)\le {p+9\over2(p-1)}<\frac34.               \tag{15.752.5}
\]

Thus some paired cube has mean exactly \(1/2\). Proposition 15.751's
dimension-free half-mean theorem bounds its maximum by three, while
(15.752.3) gives \(H\ge4\). This excludes \(H\ge2\).

For \(H=1\), complement the middle slice and write the resulting Boolean
quadratic as \(f:J(p,(p-1)/2)\to\{0,1\}\), with

\[
                         \mu={p+9\over4p}.           \tag{15.752.6}
\]

The relevant-transposition derivative floor in Proposition 15.751 is
unchanged. Combining it with the correctly normalized Johnson Laplacian
identity shows that the complement \(L\) of the largest zero-influence
coordinate class obeys

\[
 L\le {6(p-1)(p-2)(p+9)\over p^2(p+1)}<7.          \tag{15.752.7}
\]

After clearing denominators, the strict inequality is

\[
 p^3-29p^2+150p-108>0.                              \tag{15.752.8}
\]

The left side is \(168\) at \(p=23\), its derivative is already positive
there, and the derivative is increasing thereafter. Hence \(L\le6\).
Symmetrizing over the largest class extends \(f\) to a Boolean cube
quadratic, and the cube influence bound reduces it to four actual
coordinates. The fixed catalog from Proposition 15.751 has only the
density values

\[
 0,1,{p-3\over4p},{p+1\over4p},{p-1\over2p},{p+1\over2p},
 {3p-1\over4p},{3(p+1)\over4p}.                    \tag{15.752.9}
\]

The target in (15.752.6) lies strictly between the fourth and fifth values,
which proves (15.752.1).

Now put

\[
 k=4p+2t,\qquad |H_G|=4p+2t+1,
 \qquad q={p-1\over2},\quad m=q+1,                 \tag{15.752.10}
\]

and assume

\[
\begin{array}{ll}
 p\equiv1\pmod4:&4\le t\le q-4=(p-9)/2,\\
 p\equiv3\pmod4:&4\le t\le q-3=(p-7)/2.
\end{array}                                        \tag{15.752.11}
\]

The isolated-vertex count \(p^2+1-2|H_G|\) is positive throughout this
range. Signed transport therefore gives \(I=0\), an all-finite even
boundary, and even \(b_d\) in every direction. In the phase-one type,

\[
 a_d=2u+(p+1)k_d,qquad \sum_dk_d=m+t-u.            \tag{15.752.12}
\]

The exact phase-one floors and the sharp \(p-3\) integral-lift floor reduce
(15.752.12) to the same three endpoint branches as Propositions
15.734--15.735. Their hard parallel counts, opposite minima, and opposite
parallel-count totals are

\[
\begin{array}{c|c|c|c|c}
 &P&Q_{\min}&a(Q_{\min})&\sum Q\\ \hline
 A&4&3&8&4q+t\\
 B&5&2&6&3q+t\\
 C&3&4&8&5q+t+1.
\end{array}                                        \tag{15.752.13}
\]

Here \(B\) occurs only for \(p\equiv1\pmod4\), and \(C\) only for
\(p\equiv3\pmod4\). The coefficient congruences force the displayed \(P\)
with quotient zero; a positive quotient makes the opposite edge total
negative in (15.752.11).

The minimum means \(8,6,8\) are below both the nonzero phase-zero fibre
floor and the \(p-3\) lift floor, so no direction has \(Q=Q_{\min}\).
After raising every \(Q\) once, (15.752.13) leaves respectively only
\(t-4,t-3,t-4<m\) surplus units. Thus some direction has

\[
 A,C:\quad a=p+9,qquad B:\quad a=p+7.             \tag{15.752.14}
\]

At either mean, every nonzero-\(b_d\) option is an exact parity baseline
plus fewer than \(p-3\) units. Hence \(b_d=0\), the phase-zero slack is
\(2B_d\), and (15.752.14) contradicts (15.752.1) in branches \(A,C\) or
Proposition 15.751 in branch \(B\). Consequently

\[
 \boxed{\begin{array}{ll}
 p\equiv1\pmod4:&4\le t\le(p-9)/2,\\
 p\equiv3\pmod4:&4\le t\le(p-7)/2
 \end{array}
 \quad\Longrightarrow\quad k=4p+2t\text{ is impossible}.}  \tag{15.752.15}
\]

In particular, \(k=4p+8\) is impossible for every prime \(p\ge23\), for
every boundary size. This is a proved infinite-family theorem, not a prime
or configuration census. Its local threshold is real: at \(p=19\), if
\(|R|=4\), \(r=|X\cap R|\), and

\[
 B=3-2r+{r\choose2},                                \tag{15.752.16}
\]

then the layer values are \(3,1,0,0,1\) and
\(4p\mathbb EB=2p-10=p+9\). This local object is not asserted to be a
residual graph. The endpoint \(p=13\), the prior \(p\le11\) gate, and the
layers beyond (15.752.15) remain open. Proposition 15.753 below closes the
two other fifth-shell endpoints \(p=17,19\), and Proposition 15.754 closes
the \(p=13\) endpoint. Residual (ii), E(1), \(L=1/2\), and the original
convergence problem are not closed.

Evidence: `src/e1_gmin_m4_prop15752.py`,
`tests/test_prop15752.py`,
`evidence/e1_gmin_m4_prop15752.json`, and
`evidence/NOTE_2026-09-01_RESIDUAL_BAND_INFLUENCE_CLOSE.md`.

## Proposition 15.753 — the two exceptional fifth-shell endpoints close

There is no residual-(ii) boundary at

\[
 (p,k)=(17,76)\quad\hbox{or}\quad(19,84).          \tag{15.753.1}
\]

This is an exact finite aggregate-row theorem. It enumerates neither graphs
nor coefficient cells. Put

\[
 q={p-1\over2},\qquad m=q+1,\qquad t=4,\qquad
 a_d=2u+(p+1)k_d,\qquad \sum_d k_d=m+4-u.          \tag{15.753.2}
\]

We first derive the complete branch list. For \(0\le u\le4\), a row with
\(k_d=0\) has mean \(2u\), below the least phase-one floor. A row with
\(k_d=1\) has mean \(p+1+2u\). The exact phase-one floor table, followed by
the sharp nonzero integral-lift floor \(p-3\), says that a nonexact such row
would require \(k_d\ge2\). The only exact low-row option is the \(b=16\)
literal at \(p=17,u=0\). If that option is absent, all \(m\) directions have
\(k_d\ge2\), but

\[
                         m+4-u<2m.                 \tag{15.753.3}
\]

Thus every low row is impossible except the branch containing that literal;
at \(p=17,u=0\), (15.753.2) forces at least five exact literals.
For \(4<u<m-1\), \(k_d=0\) is still below the floor, whereas
\(m+4-u<m\), so those rows are impossible. At \(u=m-1\), the exact floor
allows only \(b=2\) for \(p=17\), and \(b\in\{2,18\}\) for \(p=19\). The two
\(p=19\) types cannot mix, since their common parallel count would have to
satisfy both \(P\equiv4\pmod9\) and \(P\equiv3\pmod9\). Hence the exhaustive
branches are

\[
\begin{array}{c|c}
p=17&A:\ b=2\text{ XNOR},\quad B:\ b=16\text{ literal},\\
p=19&A:\ b=2\text{ XNOR},\quad C:\ b=18\text{ complement literal}.
\end{array}                                        \tag{15.753.4}
\]

For a projective direction \(L\), let \(P_L\) be its hard parallel count and
let \(q_L(a)\), \(1\le a\le q\), be the normalized nonzero part of the common
difference-Radon row. Solving the coefficient-offset congruence in the full
integral range first fixes the hard/opposite edge split and
\(hT=|E_h|-|E_{-h}|\). The common graph then gives

\[
 \sum_aq_L(a)=hT-P_L,\qquad
 \sum_aq_{-L}(a)=-hT-Q_L.                          \tag{15.753.5}
\]

On the other hand, the unspecialized local hard-row sum is
\(p(P_L-3)-c-(p+1)k_L\), where \(c=0\) in branch \(17B\) and \(c=p-1\)
otherwise. Equating the two expressions in (15.753.5) gives, with no other
integer solution,

\[
\begin{array}{c|c|c|c|c}
\text{branch}&(|E_h|,|E_{-h}|)&hT&P_L&\text{exact-row energy}\\ \hline
17A&(41,36)&5&4+k_L&1\\
17B&(49,28)&21&4+k_L&32\\
19A&(45,40)&5&4+k_L&1\\
19C&(35,50)&-15&3+k_L&36.
\end{array}                                        \tag{15.753.6}
\]

Thus no locally admissible row may choose an independent normalization.
For all four branches the common difference-Radon Parseval identity is

\[
 \sum_{L,a>0}q_L(a)^2
 =p|E(H)|+2(hT)^2-2\sum_LP_L^2+2pC,\qquad
 C=\sum_\delta {m_\delta\choose2}\ge0.             \tag{15.753.7}
\]

Every row below is constrained by all translation-averaged middle-slice
cuts. There are exactly \(698\) distinct cut rows at \(p=17\) and \(2338\)
at \(p=19\). For each displayed finite maximum, an integral row attains it,
and the same integer model with energy at least one larger is exactly
infeasible under one-worker CP-SAT. Empty systems are replayed directly.

In branch \(17B\), at least five exact literal directions force
\(M_2=M_4=0\). Hard excess one is infeasible, so only the partitions
\((2,2)\) and \((4)\) of the total excess four remain. The exact hard maxima
are \(70,218\), and the opposite \(Q=3,4\) maxima are \(72,101\). After the
known exact-row energies are removed from (15.753.7), the two comparisons
are

\[
\begin{array}{c|c|c|c}
\text{partition}&\text{Parseval lower bound}&\text{row upper bound}&\text{gap}\\ \hline
(2,2)&1245+34C&817&428\\
(4)&1197+34C&895&302.
\end{array}                                        \tag{15.753.8}
\]

Thus branch \(17B\) is empty. In branch \(17A\), a \(Q=3\) opposite row has
mean \(8\). For nonzero \(b\) this is below the phase-zero floor \(16\), and
for \(b=0\) it is a nonzero integral lift below the sharp floor \(14\).
Thus all nine opposite counts are \(Q=4\). The raw hard maxima for excess
\(1,2,3,4\) are
\(28,81,200,289\), and the raw opposite maximum is \(23\). They give strict
gaps

\[
 342,312,282,212,182,138                            \tag{15.753.9}
\]

for the first six partitions of five. For the last partition \((5)\), eight
hard directions are exact XNOR rows. With
\(S_j(L)=\epsilon_LM_j(L)\), where \(\epsilon_L=h\) on hard directions and
\(\epsilon_L=-h\) on opposite directions, the binary quartic

\[
                         G(L)=hM_4(L)-M_2(L)^2      \tag{15.753.10}
\]

has eight distinct projective roots, hence is identically zero. The sign
change on opposite rows is essential:

\[
 \text{hard}:\ S_4=S_2^2,\qquad
 \text{opposite}:\ S_4=-S_2^2\pmod {17}.           \tag{15.753.11}
\]

Under (15.753.11), the hard excess-five maximum is \(384\) and the opposite
\(Q=4\) maximum is \(11\). Consequently

\[
                         483<645+34C,               \tag{15.753.12}
\]

with gap \(162\). This closes branch \(17A\) and proves the \(p=17\) case.

In branch \(19C\), at least five exact complement literals force
\(M_2=M_4=0\). A \(Q=4\) opposite row has mean eight, below the nonzero-\(b\)
phase-zero floor \(20\) and, when \(b=0\), below the sharp integral-lift
floor \(16\). Thus all ten opposite rows have \(Q=5\). Such a row would satisfy

\[
 \sum_aq_a=10,\qquad \sum_a|q_a|\le80,\qquad
 c_X\mathbin\cdot q\le114\quad(2338\text{ cut rows}),\qquad
 M_2=M_4=0\pmod {19}.                              \tag{15.753.13}
\]

The exact integer system (15.753.13) is infeasible, closing branch \(19C\).
Finally, in branch \(19A\), the same two-part floor argument excludes
\(Q=3\), so all ten opposite counts equal four. The hard maxima for excess
\(1,2,3,4,5\) are
\(36,97,194,325,494\), and the opposite maximum is \(23\). The seven
partitions of five give the strict Parseval gaps

\[
                 520,490,460,420,390,312,162.      \tag{15.753.14}
\]

Thus branch \(19A\) is empty as well. All branches in (15.753.4) have been
excluded, proving (15.753.1). Together with Proposition 15.752,

\[
 \boxed{k=4p+8\text{ is impossible for every prime }p\ge17.} \tag{15.753.15}
\]

Among \(p\ge13\), this leaves \(p=13,k=60,u=6\) as the sole fifth-shell
hole; Proposition 15.754 below closes that endpoint. The prior \(p\le11\)
gate still remains. Proposition 15.753 does not close later shells, residual
(ii) globally, E(1), or the original convergence problem. Evidence:
`src/e1_gmin_m4_prop15753.py`,
`tests/test_prop15753.py`, `evidence/e1_gmin_m4_prop15753.json`, and
`evidence/NOTE_2026-09-01_P17_P19_FIFTH_SHELL_CLOSE.md`.

## Proposition 15.754 — the p13 fifth-shell endpoint closes

There is no residual-(ii) boundary in the last fifth-shell branch

\[
                    (p,t,k,u)=(13,4,60,6).          \tag{15.754.1}
\]

This is an exhaustive finite aggregate/common-form certificate and a proved
endpoint theorem. It is not a graph, orbit, cell, support, or common-
realization census.

There are seven hard projective directions. Write their local means as

\[
 a_L=12+14e_L,\qquad e_L\ge0,\qquad \sum_Le_L=5.   \tag{15.754.2}
\]

Thus at least two directions have \(e_L=0\), and the exact floor
classification makes each such row a \(b=2\) XNOR row. Solving the full
coefficient-offset congruence and gluing the unspecialized row sums to the
common signed edge total, before normalizing any individual row, gives

\[
 (|E_h|,|E_{-h}|)=(33,28),\quad hT=5,\quad
 P_L=4+e_L.                                        \tag{15.754.3}
\]

The phase-zero and integral-lift floors exclude opposite parallel count
three, so all seven opposite directions have \(Q_L=4\). If
\(q_L=(q_L(1),\ldots,q_L(6))\) is the nonzero distance row, every hard row
of excess \(e\) and every opposite row respectively satisfy

\[
\begin{array}{c|c|c|c}
 &\sum_aq_L(a)&\sum_a|q_L(a)|&c_X\mathbin\cdot q_L\quad(74\text{ cuts})\\ \hline
 \text{hard }e&1-e&\le57-e&\le13\\
 \text{opposite}&-9&\le57&\le-52.
\end{array}                                        \tag{15.754.4}
\]

An exact XNOR row has energy one. The common difference-Radon Parseval
identity therefore makes the required energy of the nonexact rows

\[
                         B_\lambda+26C,             \tag{15.754.5}
\]

where \(C=\sum_\delta {m_\delta\choose2}\) is the displacement-collision
count. The seven positive partitions of the total excess five, their exact
XNOR root counts, collision minima, and bases are

\[
\begin{array}{c|rrrrrrr}
\lambda&(1^5)&(2,1^3)&(2,2,1)&(3,1,1)&(3,2)&(4,1)&(5)\\ \hline
\#\text{ roots}&2&3&4&4&5&5&6\\
C_{\min}&0&0&0&1&1&2&3\\
B_\lambda&303&298&293&289&284&276&259.
\end{array}                                        \tag{15.754.6}
\]

We exclude these seven partitions in three exact regimes.

For \((1^5)\) and \((2,1^3)\), use the common binary forms

\[
                 U=hM_2,\qquad G=hM_4-M_2^2.       \tag{15.754.7}
\]

The exact XNOR directions are roots of \(G\), so after factoring their root
product only the remaining coefficients are free. For both hard signs,
every exact-root set and every excess assignment, exhaustive finite-field
coefficient enumeration against the exact row-energy tables gives maximum
separable nonexact energies \(293\) and \(290\), respectively. The raw
energy/collision ledger leaves only \(C=0,1,2\); throughout those strata

\[
 293<303+26C,\qquad 290<298+26C.                    \tag{15.754.8}
\]

For the two four-root partitions \((2,2,1)\) and \((3,1,1)\), the common
forms refine to

\[
 U=hM_2,\qquad G=R_4c,\qquad
 J_6=h(M_6-M_2^3)=R_4Q_2.                          \tag{15.754.9}
\]

The opposite normalization is sign-sensitive. With
\(N_{2r}=(-h)M_{2r}\), the hard and opposite table keys are

\[
\begin{array}{c|c}
\text{hard}&(U,G,J_6)=(N_2,N_4-N_2^2,N_6-N_2^3),\\
\text{opposite}&(U,G,J_6)=(-N_2,-N_4-N_2^2,-N_6+N_2^3).
\end{array}                                        \tag{15.754.10}
\]

In particular the exact sixth-moment weight vector is

\[
                         W_6=(1,12,1,1,12,12),      \tag{15.754.11}
\]

computed in ordinary integers before reduction modulo thirteen. The
explicit coefficient join checks both hard signs, all 70 hard-sign/root-set
choices, all excess assignments, and the complete exact row-energy bitsets.
For each partition it starts with \(218{,}320{,}284\) form triples, leaves
336 after all seven opposite rows, and performs 1,008 hard-assignment checks.
No target energy survives: \((2,2,1)\) is excluded for \(C=0,1,2,3\), and
\((3,1,1)\) for \(C=1,2,3,4\). Those are exactly their collision ranges.

Finally, each of \((3,2),(4,1),(5)\) has at least five exact XNOR roots.
The binary quartic \(G\) is therefore identically zero, so every remaining
hard row has \(N_4=N_2^2\), whereas every opposite row has
\(N_4=-N_2^2\). Exact one-worker six-bin models impose this relation, all
74 translated cuts, and the signed six-positive/seven-negative collision
floor. The resulting ledgers are

\[
\begin{array}{c|c}
(3,2)&310>288,\\
(4,1)&\text{gaps }28,50\quad(C=2,3),\\
(5)&\text{gaps }38,42,56,58,82\quad(C=3,4,5,6,7).
\end{array}                                        \tag{15.754.12}
\]

The raw bounds exclude larger collision counts, so (15.754.6) is exhausted
and (15.754.1) follows. Together with Propositions 15.752--15.753,

\[
 \boxed{k=4p+8\text{ is impossible for every prime }p\ge13.} \tag{15.754.13}
\]

This closes the fifth shell, not residual (ii) globally. Propositions
15.768--15.770 later close the next two generic post-band layers and the
\(p=23,t=9,10\) endpoints. Proposition 15.771 adds \(p=23,t=11\), and
15.772 adds \(p\equiv1\pmod4,p\ge29,t=q-1\).
With \(q=(p-1)/2\), the exact live frontier is
\(p=5,7\); \(p=11,t\ge3\) (\(k\ge50\)); \(p=13,17,19,t\ge5\);
\(p=23,t\ge12\); \(p\equiv1\pmod4,p\ge29,t\ge q\);
\(p\equiv3\pmod4,p\ge31,t\ge q\); and the positive \(p=7,z=7\)
boundary systems. Thus residual (ii), E(1), \(L=1/2\), and the original
convergence problem remain open.

Evidence: `src/e1_gmin_m4_prop15754.py`,
`src/p13_u6_cut_equalities.py`, `src/p13_u6_high_root_energy.py`,
`scripts/p13_u6_joint_ug_tables.py`,
`scripts/p13_u6_low_root_ug_bound.py`,
`scripts/p13_u6_four_root_ugj.py`, `tests/test_prop15754.py`,
`evidence/e1_gmin_m4_prop15754.json`,
`evidence/e1_gmin_m4_prop15754_low_root_ug.json`,
`evidence/e1_gmin_m4_prop15754_four_root_221.json`,
`evidence/e1_gmin_m4_prop15754_four_root_311.json`, and
`evidence/NOTE_2026-09-02_P13_U6_COMMON_FORM_CLOSE.md`.

The next five propositions use the separately audited reductions recorded as
Propositions 15.755--15.756.  Those reductions leave a common simple graph
\(H\), and rule out a boundary-only Weil estimate as the missing input.  The
full edge--Radon map of that one graph is therefore the object below.  The
numbering follows the audited source records even though the full proof
sections for 15.755--15.756 are not reproduced here.

## Proposition 15.757 — exact binary edge--Radon image

Let \(V=\mathbf F_p^2\), where \(p\) is odd.  In every one of the \(p+1\)
projective directions \(L\), record modulo two the off-diagonal fibre-pair
counts \(K_L(s,t)\), \(s<t\), and the total parallel count \(P_L\).  Denote
the resulting map on simple edge variables by \(R_2\).  Then

\[
 \operatorname {rank}_{\mathbf F_2}R_2
  =p^2-1+(p+1){p-1\choose2}.                         \tag{15.757.1}
\]

Its image is exactly the data with a common even boundary word and the
ordinary total-parity compatibilities: the row boundaries of the fibre
graphs are the affine Radon transforms of one even word on (V),

\[
 P_L+\sum_{s<t}K_L(s,t)=h\quad\hbox{for every }L,
 \qquad \sum_LP_L=h                                  \tag{15.757.2}
\]

in \(\mathbf F_2\).  There are no further binary cross-direction
conditions.

Indeed, split the off-diagonal \(K_p\)-edge space in direction \(L\) into
its row-boundary quotient and Eulerian cycle space.  Proposition 15.692
identifies the direct sum of all row-boundary quotients with the even point
space, of dimension \(p^2-1\).  Over an algebraic closure of
\(\mathbf F_2\), the cycle space in direction \(L\) has basis

\[
 \zeta^{as+bt}+\zeta^{bs+at},
 \qquad \{a,b\}\subset\mathbf F_p^*,\ a\ne b.       \tag{15.757.3}
\]

It therefore has dimension \({p-1\choose2}\).  Distinct projective
directions have disjoint nonzero Fourier support in \((V^*)^2\), so their
cycle spaces are independent.  The \(p+1\) independent equations in
(15.757.2) give the remaining codimension, proving (15.757.1) and
completeness.

This binary theorem does not eliminate the sharp local survivor.  For
distinct \(a,b,c\),

\[
 B=x_ax_b+x_c-x_ax_c-x_bx_c\in\{0,1\},\qquad
 4B=1+z_az_b-z_az_c-z_bz_c,                          \tag{15.757.4}
\]

and on \(J(p,(p+1)/2)\),

\[
                         4p\,\mathbf EB=p+1.          \tag{15.757.5}
\]

Its coefficient graph is an Eulerian signed triangle modulo two.  The
omitted-pair and all-equal-triple atoms of scaled mass \(p-3\) are Eulerian
as well.  For the compact aggregate rows, the hard boundary words are one
on every fibre except their literal fibre and the opposite boundary words
vanish.  Prescribing the corresponding Fourier coefficients independently
on the \(p+1\) dual lines gives one even point word; (15.757.2) follows from
the common odd edge total and \(\sum_LP_L=|H|\).  Hence every such compact
aggregate target in the recorded `p=1 mod 4` family has a binary edge lift.

This is a proved image theorem and a proved method barrier, not a graph
construction.  It controls neither integral signed multiplicities nor the
nonnegative simple-edge box.  Residual (ii), E(1), and the limit remain
open.

Evidence: `src/e1_gmin_m4_prop15757.py`,
`tests/test_prop15757.py`, and `evidence/e1_gmin_m4_prop15757.json`.

## Proposition 15.758 — sharp coefficient cancellation and two local rays

Put \(m=(p+1)/2\) and \(z_i=2x_i-1\) on \(J(p,m)\), so
\(\sum_i z_i=1\).  Once the offset \(P\) is fixed, the representation

\[
             4B=P+\sum_{i<j}K_{ij}z_iz_j             \tag{15.758.1}
\]

is unique, and

\[
             4p\,\mathbf EB=pP-\sum_{i<j}K_{ij}.     \tag{15.758.2}
\]

If \(a,b\ge0\), \(P=a+b\), and the prescribed scaled mass is
\(a(p-3)+b(p+1)\), then every representation has
\(\sum K_{ij}=3a-b\), and the exact minimum coefficient norm, even among
nonnegative integer-valued quadratics, is

\[
 \min\sum_{i<j}|K_{ij}|=
 \begin{cases}
 |3a-b|,&|3a-b|\ne1,\\
 3,&|3a-b|=1.
 \end{cases}                                        \tag{15.758.3}
\]

The lower bound is \(\sum|K_{ij}|\ge|\sum K_{ij}|\).  If the latter is one,
norm one would make two cube values differ by two modulo four, while norm
two has the wrong parity.  The bounds are attained on three coordinates.
Writing \(T=z_1z_2+z_1z_3+z_2z_3\), the four forms

\[
 1+T,\quad3-T,\quad1+z_1z_2-z_1z_3-z_2z_3,
 \quad2-2z_1z_2                                      \tag{15.758.4}
\]

are four times nonnegative integral quadratics and carry, respectively,
one \(p-3\) unit, three \(p+1\) units, one \(p+1\) unit, and two \(p+1\)
units.  In particular three opposite sharp atoms cancel one same-sign atom
to a constant; no coefficient-capacity lower bound growing with \(a+b\)
can be valid.

There are also two infinite local directional survivor rays.  For a prime
\(p=4r+1\ge17\), branch B with \(u=0\), they exist throughout

\[
       2r^2-5r\le t\le4r^2-6r-3.                    \tag{15.758.5}
\]

There are \(m=2r+1\) hard and \(m\) opposite directions.  The hard rows
have \(e_L\ge0\), \(\sum e_L=t\), cell \(x_j+2B\), and

\[
 4p\,\mathbf EB=e_L(p+1),\qquad P_L=5+e_L.          \tag{15.758.6}
\]

The opposite rows have \(Q_L\ge r\), \(\sum Q_L=6r+t\), and half-lift mass

\[
 (r-1)(p-3)+(Q_L-r)(p+1)=(p+1)Q_L-2p+4.             \tag{15.758.7}
\]

For a prime \(p=4r+3\ge19\), branch C with \(u=m-1\), the interval is

\[
       2r^2-4r-2\le t\le4r^2-2r-5.                  \tag{15.758.8}
\]

Now \(m=2r+2\), \(\sum e_L=t+1\), the hard cell is
\(1-x_j+2B\), and

\[
 4p\,\mathbf E(\text{hard row})=p-1+e_L(p+1),
 \qquad P_L=3+e_L.                                  \tag{15.758.9}
\]

The opposite rows have \(Q_L\ge r+2\),
\(\sum Q_L=10r+6+t\), and half-lift mass

\[
 (r-1)(p-3)+(Q_L-r-2)(p+1)=(p+1)Q_L-4p+4.           \tag{15.758.10}
\]

Balanced allocations satisfy every isolated-chart residue, offset, parity,
directional mean, parallel-edge, and local nonnegative-integral cell
condition, with all \(p+1\) affine directions present and \(I=0\) at the
transported infinity point.  At the lower endpoint of either interval and
for every \(r\ge7\), the exact compact-template row-energy upper value
strictly exceeds the common difference--Radon Parseval lower value.  Thus
scalar Parseval supplies no contradiction on either ray.

These are local directional targets, not the images of a constructed
common graph.  The exact remaining issue is simultaneous integral
midpoint/difference--Radon consistency for one simple \(0/1\) graph.
Residual (ii), E(1), and the limit remain open.

Evidence: `src/e1_gmin_m4_prop15758.py`,
`tests/test_prop15758.py`, and `evidence/e1_gmin_m4_prop15758.json`.

## Proposition 15.759 — complete characteristic-\(p\) moment hierarchy

For the integral edge--Radon rows of one graph on \(V=\mathbf F_p^2\), set

\[
 Q_{d,k}(s,t)=(s-t)^2(st)^k(s+t)^{d-2-2k},\qquad
 C_{d,k}(L)=\sum_{s<t}K_L(s,t)Q_{d,k}(s,t).          \tag{15.759.1}
\]

For \(2\le d\le p-2\), \(0\le k<\lfloor d/2\rfloor\), and
\(0\le j\le p-d-2\), one has

\[
             \sum_{\lambda\in\mathbf F_p}
             \lambda^jC_{d,k}(L_\lambda)=0\pmod p,  \tag{15.759.2}
\]

and the projective endpoint relation is

\[
 \sum_{\lambda\in\mathbf F_p}\lambda^{p-d-1}
 C_{d,k}(L_\lambda)+C_{d,k}(L_\infty)=0\pmod p.     \tag{15.759.3}
\]

Equation (15.759.3) also applies at \(d=p-1\).  If \(p=2m+1\), retain
\(m-1\) of its \(m\) rows: their missing line is the ordinary total row
\((s-t)^{p-1}\), since, with \(D=(s-t)^2\),

\[
 D^m=\sum_{k=0}^{m-1}{m-1\choose k}(-4)^kQ_{p-1,k}, \tag{15.759.4}
\]

and every coefficient is nonzero modulo \(p\).

These moment rows are independent and exhaust the characteristic-\(p\)
left kernel beyond the \(p+1\) ordinary compatibility equations.  Their
number is

\[
 \begin{split}
 S(p)&=\sum_{d=2}^{p-2}\Big\lfloor{d\over2}\Big\rfloor(p-d)+(m-1)\\
     &={(m-1)(4m^2+7m+6)\over6}.                   \tag{15.759.5}
 \end{split}
\]

For the proof, use midpoint and half-difference coordinates
\(a=(s+t)/2\), \(\delta=(s-t)/2\).  Symmetric pair functions with one common
diagonal value have the unique basis

\[
 1;\qquad a^i\delta^{2b},\quad0\le i\le p-1,
 \quad1\le b\le m.                                  \tag{15.759.6}
\]

Pullback along the \(p+1\) projective functionals separates these bidegree
blocks.  For fixed \(i,b\), the directional rank is

\[
                    \min(p+1,i+2b+1).               \tag{15.759.7}
\]

When \(i+2b<p\), this is ordinary Vandermonde rank.  In the remaining range
the only finite-field aliases are exponents differing by \(p-1\), and the
infinity vector separates the final aliased pair.  The constant and the
the \((i,b)=(0,m)\) block are disjoint: a degree-\(p-1\) homogeneous form
constant on every nonzero difference has \(p\) finite-chart roots after
subtracting that constant, and evaluation at infinity forces it to vanish.
Summing (15.759.7) gives the target rank minus \(p+1+S(p)\), proving
exhaustion.  Alternatively, (15.759.2)--(15.759.3) follow edgewise from
the vanishing of every field power sum of degree at most \(p-2\).

In residual normalization, multiply every finite outer coefficient and the
infinity term by its direction sign; the source-edge and direction-block
sign changes are unimodular and do not change the rank.  This is a genuine
common-edge obstruction beyond Proposition 15.757.  It has not yet been
shown that the compact atoms satisfy every row above degree four, so no
compact ray, simple graph, or residual case is excluded here.

Evidence: `src/e1_gmin_m4_prop15759.py`,
`tests/test_prop15759.py`, and `evidence/e1_gmin_m4_prop15759.json`.

## Proposition 15.760 — the integral image is cut out by the moments

Let \(R:\mathbf Z^{\binom V2}\to T\) be the unsigned full edge--Radon map,
and let \({\cal A}\subset T\) be the ordinary compatibility lattice

\[
 {\cal A}=\left\{(P_L,K_L):
 T_L=P_L+\sum_{s<t}K_L(s,t)=h\ (\text{all }L),\quad
 \sum_LP_L=h\right\}.                               \tag{15.760.1}
\]

For every odd prime \(p=2m+1\),

\[
             {\cal A}/R\mathbf Z^{\binom V2}
             \cong(\mathbf Z/p\mathbf Z)^{S(p)},
 \qquad S(p)={(m-1)(4m^2+7m+6)\over6}.              \tag{15.760.2}
\]

Consequently an integral target lies in the image if and only if it
satisfies (15.760.1) and every congruence of Proposition 15.759.  There is
no hidden \(p^2\)-torsion and no obstruction at another prime.

Here is the integral argument.  Taking directional-total differences and
\(\sum_LP_L-T_{L_0}\) defines a split-surjective constraint map, so
\({\cal A}\) is primitive.  Sum the \(K_L(a,\beta)\) coordinates over the
midpoint \(a\), retaining \(P_L\).  This gives an exact sequence

\[
             0\longrightarrow M\longrightarrow{\cal A}
              \longrightarrow{\cal A}_D\longrightarrow0,     \tag{15.760.3}
\]

where \(M\) is generated by
\(K_L(a,\beta)-K_L(a_0,\beta)\).  Forgetting the edge midpoint gives

\[
 0\longrightarrow E_0\longrightarrow E
   \longrightarrow\mathbf Z^\Delta\longrightarrow0,
 \qquad \Delta=(V\setminus\{0\})/\{\pm1\}.          \tag{15.760.4}
\]

The induced pure-difference map \(S\) has one column

\[
 S e_{[\delta]}=P_{L_\delta}+
        \sum_{L\ne L_\delta}e_{(L,L(\delta)^2)}.      \tag{15.760.5}
\]

Exact row intersections give \(p{\cal A}_D\subseteq\operatorname{im}S\),
while reduction modulo \(p\) has rank \((m+1)^2\).  Hence

\[
               \operatorname {coker}S
                \cong(\mathbf Z/p\mathbf Z)^{m^2-1}. \tag{15.760.6}
\]

For the midpoint kernel, if \(L(\delta)^2=\beta\), put

\[
 X_a=\sum_{u:L(u)=a}\{u-\delta,u+\delta\}.
\]

Then all other directions cancel exactly and

\[
 R_0(X_a-X_{a'})=p\bigl(K_L(a,\beta)-K_L(a',\beta)\bigr).      \tag{15.760.7}
\]

Together with the bidegree ranks of Proposition 15.759 this yields

\[
 \operatorname {coker}R_0\cong
 (\mathbf Z/p\mathbf Z)^{m(m-1)(4m+1)/6}.            \tag{15.760.8}
\]

The snake lemma applied to (15.760.3)--(15.760.5) gives

\[
 0\longrightarrow\operatorname {coker}R_0
  \longrightarrow\operatorname {coker}R
  \longrightarrow\operatorname {coker}S\longrightarrow0.     \tag{15.760.9}
\]

Thus the middle group has order \(p^{S(p)}\).  Proposition 15.759 shows
that it needs \(S(p)\) generators modulo \(p\); its invariant factors must
therefore all equal \(p\), proving (15.760.2).  Diagonal source and target
signs transport the result verbatim to the residual normalization.

This exhausts linear divisibility and Smith-form obstructions, but it gives
only an unrestricted signed integral preimage.  If \(\tau_e\) is the source
edge sign and \(z_0\) is one integral lift, a simple graph exists exactly
when

\[
 (z_0+\ker_{\mathbf Z}R)\cap
                  \prod_e\{0,\tau_e\}\ne\varnothing.          \tag{15.760.10}
\]

The ordered remaining gates are: first check every higher moment of the
compact atoms from Proposition 15.758; conditional on their passing, decide
(15.760.10).  Neither step is proved.  In particular an integral lift is
not a nonnegative simple graph, and residual (ii), E(1), and the limit
remain open.

Evidence: `src/e1_gmin_m4_prop15760.py`,
`tests/test_prop15760.py`, and `evidence/e1_gmin_m4_prop15760.json`.

## Proposition 15.761 — full real spectrum and least-norm barrier

Let \(C={p\choose2}=pm\) and \(d=p+1\).  The complete spectrum of
\(RR^t\) on the real ordinary image follows from the exact row
intersections

\[
 |P_L\cap P_M|=0,\qquad |P_L\cap K_M(s,t)|=p,
 \qquad |K_L(a,b)\cap K_M(s,t)|=2                  \tag{15.761.1}
\]

for \(L\ne M\), together with row norms \(p^2m\) and \(p^2\).  If the
common directional total is \(T\), write

\[
 K_L=k_L\mathbf1+w_L,\quad\sum w_L=0,\quad
 k_L={T-P_L\over C},\qquad P_L={T\over d}+q_L,\quad\sum q_L=0. \tag{15.761.2}
\]

The three orthogonal blocks, their eigenvalues, and multiplicities are

\[
\begin{array}{c|c|c}
\text{block}&\text{eigenvalue}&\text{multiplicity}\\ \hline
w_L&p^2&d(C-1)\\
(P_L,K_L)=(q_L,-q_L/C)&p(C+1)&p\\
(P_L,K_L)=(a,a/m)&p^2(m+p)&1.
\end{array}                                                   \tag{15.761.3}
\]

Consequently the exact least squared norm of a real preimage is

\[
 \|R^+y\|^2={1\over p^2}\sum_L\|w_L\|^2
       +{1\over pC}\sum_Lq_L^2
       +{2T^2\over p^2(p^2-1)}.                    \tag{15.761.4}
\]

In residual normalization let \(\eta_L\) be the direction sign,
\(P_L\ge0\) the actual parallel count, and \(W_L\) the normalized
off-diagonal row.  Undoing the orthogonal row signs turns (15.761.4) into

\[
\begin{split}
 {\cal Q}(W,P,T)={}&{1\over p^2}\sum_L
 \left(\|W_L\|^2-{(\eta_LT-P_L)^2\over C}\right)\\
 &+{1\over pC}\sum_L\left(\eta_LP_L-{T\over p+1}\right)^2
 +{2T^2\over p^2(p^2-1)}.                           \tag{15.761.5}
\end{split}
\]

For a simple graph \(H\), the signed source vector is
\(z_e=\tau_e1_{e\in H}\), so \(\|z\|^2=|H|\).  Orthogonal projection onto
\((\ker R)^\perp\) gives the necessary full-midpoint inequality

\[
                         {\cal Q}(W,P,T)\le|H|.      \tag{15.761.6}
\]

This is strictly finer than scalar difference--Radon Parseval, but it still
does not separate either compact ray.  For the balanced allocations of
Proposition 15.758, uniformly in all atom labels,

\[
\begin{array}{c|c|c|c}
 &\max e_L&\max Q_L& (\|W_{\rm hard}\|_1,\|W_{\rm opp}\|_1)\\ \hline
p=4r+1&2r-3&2r-1&(10r-9,14r-10)\\
p=4r+3&2r-2&2r+2&(10r-4,6r-3).
\end{array}                                                   \tag{15.761.7}
\]

Dropping the negative terms in (15.761.5), using
\(\|W\|_2\le\|W\|_1\), and bounding each remaining aggregate term by one
gives

\[
 {\cal Q}<{(2r+1)((10r-9)^2+(14r-10)^2)\over(4r+1)^2}+2,      \tag{15.761.8}
\]

\[
 {\cal Q}<{(2r+2)((10r-4)^2+(6r-3)^2)\over(4r+3)^2}+2.        \tag{15.761.9}
\]

The minimum graph sizes on the two rays are \(4r^2+6r+5\) and
\(4r^2+8r+9\).  After clearing the positive denominators in their gaps and
putting \(u=r-7\), the two numerators are

\[
 64u^4+1328u^3+9796u^2+29864u+30706,
\quad
 64u^4+1744u^3+18108u^2+85374u+154867.              \tag{15.761.10}
\]

Every coefficient is positive.  Hence (15.761.6) holds with strict room
for both compact rays for all \(r\ge7\), uniformly for every choice of atom
labels and hence for any label choice that might satisfy the moment equations.

This is a proved full-spectrum theorem and a strictly stronger proved
method barrier, not a common graph construction.  The least-norm real
preimage need not be integral, sign-compatible, nonnegative, or \(0/1\).
The live gates remain the higher compact-atom moments and then the affine
Boolean-box intersection (15.760.10).  Residual (ii), E(1), and the limit
remain open.

Evidence: `src/e1_gmin_m4_prop15761.py`,
`tests/test_prop15761.py`, and `evidence/e1_gmin_m4_prop15761.json`.

## Proposition 15.762 — integral conference cube gap

Let `C` be a symmetric conference matrix of order `n=p^2+1`, with odd
`p>=5`, and put `Q_C(x)=x^tCx/2` on the Boolean cube. If
`delta=pn/2-Q_C(x)` and `z=(C-pI)x/2`, then `z` is integral,
`Cz=-pz`, and `||z||^2=p delta`. After switching by `x`, the conference
equations force every coordinate of `w=diag(x)z` to have one parity.

If `w=2v`, then `Dv=-pv` and a maximum coordinate gives
`||v||_1 >= (p+1)||v||_infinity`; if `w` is odd, then
`||w||^2>=p^2+1`. Together with the parity of `delta`, these facts exclude
gaps two, four, and six. Hence either `Cx=px` or

\[
 Q_C(x)\le {p(p^2+1)\over2}-8.
\]

Applying this to `-C` shows that a class with no Boolean `+p` or `-p`
eigenvector is already a complete eight-gap certificate. At gap eight the
necessary form is `w=2v`, where `v` has `p-2` entries `+1`, `p+2` entries
`-1`, and all other entries zero. This conference cube gap does not produce
such a class and does not close residual (ii), E1, or the original limit;
all remain OPEN. The full proof and normalization audit are in
`evidence/NOTE_2026-09-03_CONFERENCE_CUBE_GAP.md`.

## Proposition 15.763 — signed affine-alias incidence bound

Continue Proposition 15.755's dangerous shared-maximizer setup and suppose
the active Boolean point is an odd-parameter affine alias. Thus

\[
 \delta_\epsilon(x)=2pr^2,\qquad
 m={p+1\over2}+r\le p,\qquad r\ge1\text{ odd},
\]

and flipping any union \(T_J\) of \(r\) among the \(m\) positive parallel
fibres gives an \(\epsilon p\)-Boolean eigenvector. Put
\(a_h=\epsilon C_hx_ux_v\) for \(h=\{u,v\}\in H\). The active-state identity
is

\[
 \sum_{h\in H}a_h=2-pr^2.                                 \tag{15.763.1}
\]

For every alias \(w_J=x^{T_J}\), the four-unit gap and the odd cardinality of
\(H\) give \(\epsilon S_H(w_J)\ge3\). Hence

\[
 -\sum_{h\text{ crossing }T_J}a_h\ge {pr^2+1\over2}.       \tag{15.763.2}
\]

There are \(N={m\choose r}\) aliases and every edge crosses at most
\(M=2{m-2\choose r-1}\) of them. Retain the signs when summing
(15.763.2). By (15.763.1), the number of negative \(a_h\)'s is
\((|H|+pr^2-2)/2\), whereas positive edges only decrease the signed sum.
Therefore

\[
 {N(pr^2+1)\over2}
 \le M{|H|+pr^2-2\over2},
\]

and

\[
 \boxed{|H|\ge
 { (pr^2+1)m(m-1)\over2r(m-r)}-pr^2+2.}                   \tag{15.763.3}
\]

The exact integral statement rounds the right side up to the next odd
integer. At \(r=1\), (15.763.3) is already the odd integer

\[
 \boxed{|H|\ge {p^2+11\over4}},                            \tag{15.763.4}
\]

strictly improving Proposition 15.755's parity-adjusted
\((p+1)(p+3)/8\) bound.

There is a further exact alternative for the distinguished active edge
\(e\). If an alias has \(\epsilon S_H(w_J)=3\) and its cut does not cross
\(e\), then \(\epsilon S_{H\setminus\{e\}}(w_J)=2\). If no such alias exists,
every noncrossing alias has H-score at least five. The number of noncrossing
aliases and the negative signed contribution of \(e\) cancel in the summed
inequality, and (15.763.3) holds with \(pr^2+1\) replaced by \(pr^2+3\).

This is a proved conditional theorem for the affine-alias subfamily, not a
classification of the full defect-\(2p\) shell. The three-coordinate integral
eigenvector branch can be nonaffine, and deletion-specific affine coordinate
systems need not agree. Thus no all-deletions contradiction follows;
residual (ii), E1, and the original limit remain OPEN. Full proof and exact
arithmetic are in `evidence/NOTE_2026-09-04_SIGNED_AFFINE_ALIAS_BOUND.md` and
`src/e1_gmin_m4_prop15763.py`.

## Proposition 15.764 — exact parity shell bridge for a minimal four-gap set

For \(\epsilon\in\{+1,-1\}\), put

\[
 E_\epsilon=\{y\in\{\pm1\}^n:Cy=\epsilon py\},\qquad
 T_F^\epsilon(y)=\epsilon\sum_{\{u,v\}\in F}C_{uv}y_uy_v,
\]

and let \(m_\epsilon(F)=\min_{E_\epsilon}T_F^\epsilon\). Suppose that

\[
 \Phi(C\mathbin\triangle H)=\Phi-4,
 \qquad
 \Phi(C\mathbin\triangle(H\setminus\{e\}))=\Phi-2
 \quad(e\in H).                                           \tag{15.764.1}
\]

Eigenshell evaluation gives, in both phases,

\[
 T_H^\epsilon\ge2,\qquad T_{H\setminus\{e\}}^\epsilon\ge1,\qquad
 T_F^\epsilon\equiv|F|\pmod2.                             \tag{15.764.2}
\]

Write \(b_e=\epsilon C_{uv}y_uy_v\). Since
\(T_{H\setminus\{e\}}=T_H-b_e\), if \(|H|\) is odd then

\[
 \boxed{
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})=2
 \iff
 \exists\epsilon,y\in E_\epsilon:\ T_H^\epsilon(y)=3.}  \tag{15.764.3}
\]

Indeed a deletion score two makes the odd H-score one or three, and
(15.764.2) rules out one. Conversely, a sum of signs equal to three has a
positive edge, whose deletion has score two. Here \(m_+=s_+\) and
\(m_-=-s_-\), so the left side is exactly the official even-cardinality
residual-(ii) entry level.

This is the full official entry, not merely its numerical level. On every
deletion row with score two, (15.764.2) forces \(b_e=+1\), so the edge freezes
positive on the entire critical level; both deletion phases have shell floor
two. A minus-phase row can be normalized to plus: multiplication by a
nonsquare in \(\mathbb F_{p^2}\) negates the finite Paley block, and switching
at infinity negates the remaining entries, proving \(-C=DP^tCPD\). Finally
the H-score floor three and (15.764.5) give \(|H|\ge3p\); equality would be
bi-tight level three, excluded by Proposition 15.720. Hence
\(|H|\ge3p+2\), so the even deletion has \(k\ge3p+1\), \(s_+=2\),
deep two-sidedness, and freeness failure for the distinguished edge.

If \(|H|\) is even, a deletion score is odd and can never equal two. The
correct parity alternative is

\[
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})\le2
 \iff
 \exists\epsilon,y\in E_\epsilon:\ T_H^\epsilon(y)=2,    \tag{15.764.4}
\]

and the deletion score is one, on the Type-I side of the ledger.
Every critical level-one row again has \(b_e=+1\), and the same phase
normalization gives the plus Type-I convention. The analogous frame argument,
using the excluded bi-tight level two at equality, gives
\(|H|\ge2p+2\) and \(|H|-1\ge2p+1\).

The signed frame mean is

\[
 \mathbb E_{E_\epsilon}T_H^\epsilon={|H|\over p}.         \tag{15.764.5}
\]

Thus failure of (15.764.3) for odd H forces \(|H|\ge5p\). At equality H is
bi-tight of level five. Proposition 15.720's generic degree congruence gives

\[
 d_i+d_j\equiv10p\pmod{(p^2-1)/2}.                        \tag{15.764.6}
\]

For \(p\ge11\) this forces the impossible regular degree
\(10p/(p^2+1)\). For \(p=7\), all degrees have residue 11 or 23 modulo 24,
so their sum exceeds 70. At \(p=5\), the only degree profiles are a full star
and a balanced double star. The star is a vertex switching, hence is not
deeper. For the double star, with centres \(a,b\), twelve leaves on each
side, \(c=C_{ab}\), \(\alpha_i=C_{ai}\), \(\beta_i=C_{bi}\), and
\(r_i=+1\) or \(-1\) according to the side, the `scheme+cross` decomposition
has

\[
 X_{ai}=\alpha_i r_i/2,\qquad X_{bi}=-\beta_i r_i/2,
 \qquad CX+XC=0.
\]

Its leaf-leaf entries make \(\alpha_i\beta_i\) constant on each side;
conference-row orthogonality makes the constants opposite, so
\(r_i=\tau\alpha_i\beta_i\). The \((a,j)\) entry of the anticommutator then
has left side \(-\tau c\alpha_j\) by \((C^2)_{bj}=0\), but right side
\(+\tau c\alpha_j\), a contradiction. Therefore the level-five equality
cannot be deeper, and

\[
 \boxed{|H|\text{ odd},\ |H|\le5p
 \Longrightarrow
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})=2.}    \tag{15.764.7}
\]

For even H, failure of (15.764.4) gives \(|H|\ge4p\); equality is the
bi-tight level four already excluded by Proposition 15.720. Hence the exact
unclosed failure ranges begin at

\[
 \boxed{|H|\ge4p+2\text{ even},\qquad |H|\ge5p+2\text{ odd}.}
                                                                    \tag{15.764.8}
\]

A fully specified non-Paley max-of-affine model at \(p=5,|H|=25\) realizes
all scalar frame, parity, minimal-four-gap, and all-deletions-two-gap
identities while every deletion shell minimum is four. Thus scalar identities
alone cannot bridge (15.764.8); common Paley graph structure is essential.
The construction and the complete level-five calculation are in
`evidence/NOTE_2026-09-04_MINIMAL_GAP4_SHELL_BRIDGE.md`.

This proposition repairs the deletion-to-unit audit only in the range
(15.764.7). It does not prove that an arbitrary minimal four-gap H is odd or
small enough, so residual (ii), E1, and the original limit remain OPEN.
Evidence: `src/e1_gmin_m4_prop15764.py`, `tests/test_prop15764.py`, and
`evidence/e1_gmin_m4_prop15764.json`.

## Proposition 15.765 — a nonaffine point on the first Paley defect shell

Let \(p=11\), identify
\(\mathbb F_{121}=\mathbb F_{11}[a]/(a^2-2)\), and write

\[
 \chi(x+ya)=\left({x^2-2y^2\over11}\right).
\]

Transcribe the right-hand \(11\times11\) Boolean matrix in Section 6,
Figure 4 of Kiss--Somlai, *Special directions on the finite affine plane*,
Designs, Codes and Cryptography 92 (2024), 2587--2597. It gives a
33-point set \(E_0\subset\mathbb F_{11}^2\) whose row \(y=1\) is empty.
Adjoin that full row and apply \(T(x,y)=(x,2y)\):

\[
 E=T\bigl(E_0\mathbin\cup\{(x,1):x\in\mathbb F_{11}\}\bigr),
 \qquad D=\mathbb F_{11}^2\setminus E.
\]

Thus \(|E|=44\) and \(|D|=77\). Direct exact line intersection counts show
that \(E\) is nonconstant in exactly the four directions
\(\infty,0,2,-2\), and is constant of value four on every line in each of
the other eight directions. The four direction representatives
\((0,1),(1,0),(1,2),(1,-2)\) all have Paley character \(+1\). Equivalently,
and checked pointwise over all 121 vertices, the finite Paley character
matrix \(Q_{uv}=\chi(v-u)\) satisfies

\[
 Q\mathbf1_E=11\mathbf1_E-4\mathbf1,
 \qquad
 Q\mathbf1_D=11\mathbf1_D-7\mathbf1.                 \tag{15.765.1}
\]

In particular,

\[
 3+2\sum_{v\in D}\chi(v-u)
   =11\bigl(2\mathbf1_D(u)-1\bigr)
 \quad(u\in\mathbb F_{11}^2).                         \tag{15.765.2}
\]

For the normalized order-122 Paley conference matrix

\[
 C=\begin{pmatrix}0&\mathbf1^t\\[2pt]\mathbf1&Q\end{pmatrix},
\]

define \(y_\infty=3\) and \(y_u=2\mathbf1_D(u)-1\). Since
\(\sum_u y_u=2|D|-121=33\), (15.765.2) proves, coordinate by coordinate,

\[
 Cy=11y,
 \qquad \|y\|_2^2=3^2+121=130,                         \tag{15.765.3}
\]

and infinity is the unique coordinate of magnitude three. Let
\(x=y-2e_\infty\in\{\pm1\}^{122}\). Because \(C_{\infty\infty}=0\),

\[
 x^tCx
 =y^tCy-4e_\infty^tCy
 =11\cdot130-4\cdot33=1298.
\]

Consequently

\[
 q_C(x)=649,
 \qquad
 \Phi-q_C(x)={11\cdot122\over2}-649=22=2p.            \tag{15.765.4}
\]

This shell point is genuinely nonaffine. A nonempty proper union of
parallel affine lines has exactly one nonconstant parallel-class profile:
the selected direction has only values zero and eleven, while every
transverse profile is constant. The exact set \(E\) instead has four
nonconstant directions. Hence \(E\) is not a union of four parallel lines;
by complementation, \(D\) is not a union of seven parallel lines.

Therefore the proposed classification of every one-coordinate-three
integral Paley eigenshell point by parallel affine lines is false, already at
\(p=11\). This does **not** construct a common switching set \(H\), align
different deletion representatives, or satisfy the all-deletions hypotheses
of residual (ii). It blocks only the attempted universalization of
Proposition 15.763; residual (ii), E1, and the original limit remain OPEN.
The cited-coordinate transcription, pointwise convolution, eigenvector,
quadratic value, line profiles, source hashes, and scope guards are in
`src/e1_gmin_m4_prop15765.py`, `tests/test_prop15765.py`, and
`evidence/e1_gmin_m4_prop15765.json`.

## Proposition 15.768 — the first `p=1 mod 4` post-band layer is empty

Put \(q=(p-1)/2\), \(m=q+1\), and \(k=4p+2t\).  If \(p\ge29\) is prime
and \(p\equiv1\pmod4\), then

\[
 t=q-3={p-7\over2},\qquad k=5p-7                         \tag{15.768.1}
\]

is impossible for every boundary size.

Indeed, the isolated-vertex inequality
\(p^2+1-2(5p-6)>0\) permits signed transport to an all-finite chart with
\(I=0\).  Write the hard phase-one means as

\[
 a_L=2u+(p+1)k_L,\qquad \sum_Lk_L=m+t-u.                \tag{15.768.2}
\]

The exact parity floors and Proposition 15.688's sharp \(p-3\) lift floor
leave the two old endpoint branches and one new residue \(u=t\).  The old
branches force the local masses \(p+7\) and \(p+9\), already excluded by
Propositions 15.751--15.752.  In the new residue every \(k_L=1\), every hard
mean is \(2p-6\), and equality in the positive quadrature is pointwise

\[
 A(X)=(2-|X\cap C|)^2,\qquad |C|=3.                    \tag{15.768.3}
\]

To see the last assertion also on the omitted intersection layer, restrict
\(A-(2-|X\cap C|)^2\) to the three-cube formed by three independent swaps
with \(C\).  It vanishes at the seven nonzero vertices, and its third finite
difference is zero, so it vanishes at the origin as well.

In signed coordinates, (15.768.3) has target

\[
 3+2A=5-\sum_{i\in C}z_i+
          \sum_{\{i,j\}\subset C}z_iz_j,
\]

whose coefficient offset is \(5-1-1-1=2\).  Hence the slice-kernel
comparison gives \(q\mid P_L-2\).  The two expressions for the difference
row sum of the one common graph are

\[
 \sum_a q_L(a)=p(P_L-3)-a_L=hT-P_L.                    \tag{15.768.4}
\]

All hard means agree, so (15.768.4) makes \(P_L\) common.  The edge bound
and the congruence force \(P_L=2\), whence

\[
 hT=8-3p,\qquad \sum_LQ_L=4p-7,
 \qquad a(Q)=(p+1)Q-6p+8.                              \tag{15.768.5}
\]

The row \(Q=6\) has mass \(14\), below both relevant nonzero lift floors.
Thus \(Q\ge7\) in every opposite direction, and (15.768.5) forces a
\(Q=7\) row of mass \(p+15\).  Every nonzero-boundary alternative would
leave a positive lift of mass \(14\) or \(16\), again below \(p-3\).
Consequently a survivor would yield a nonzero nonnegative integral
quadratic \(B\) on \(J(p,m)\) with

\[
 4p\,\mathbb E B=p+15.                                 \tag{15.768.6}
\]

We now exclude (15.768.6).  If \(H=\max B\ge2\), paired-cube means through
a maximizer are positive quarter-integers.  For \(p\ge37\), their average
is below \(3/4\), while the paired-cube and stabilizer bounds give
\(H\ge(p-13)/4>3\).  Mean \(1/4\) is support-floor equality and has
maximum at most one, so some cube must have mean \(1/2\), contradicting
the dimension-free half-mean theorem \(\max g\le3\).

The endpoint \(p=29\) uses the following sharp companion:

\[
 g\ge0\text{ integral quadratic on a cube},\quad
 \mathbb E g={3\over4}\quad\Longrightarrow\quad\max g\le6.   \tag{15.768.7}
\]

For completeness, take a minimum-dimensional counterexample with its maximum
at the origin.  Facet means are quarter-integral.  The support floor, the
half-mean theorem, minimality, and the complementary zero-facet case leave
only origin-facet means \(1\) and \(5/4\); every nonorigin vertex then has
value at most three.  On any five-coordinate face, permutation averaging
gives a quadratic \(f(s)\), and interpolation at \(s=1,3,5\) yields

\[
 f(0)={15\over8}f(1)-{5\over4}f(3)+{3\over8}f(5)
       \le {27\over4}<7.
\]

Dimensions at most three have total mass below seven, and in dimension four
the vanishing fourth difference splits the total mass twelve equally between
the two parity classes.  This proves (15.768.7).  At \(p=29\), the half-mean
bound first gives \(H\ge4\), so every paired cube has mean at least \(3/4\);
then \(H\ge12\), while the integral stabilizer bound gives \(H\le12\).
Quarter integrality supplies a mean-\(3/4\) cube containing the height-twelve
maximizer, contradicting (15.768.7).

Finally, if \(H=1\), then \(B\) is Boolean of density \((p+15)/(4p)\).
The corrected Johnson influence bound leaves at most seven slice
coordinates; all patterns extend to the complementary middle slice, and
cube influence leaves at most four active coordinates.  Proposition
15.751's fixed four-bit catalog contains no such density.  This proves
(15.768.1).  No prime, graph, orbit, slice, or residual-candidate census is
used beyond that already fixed four-bit catalog.

## Proposition 15.769 — the first `p=3 mod 4` post-band layer is empty

If \(p\ge31\) is prime and \(p\equiv3\pmod4\), then

\[
 t=q-2={p-5\over2},\qquad k=5p-5                         \tag{15.769.1}
\]

is impossible for every boundary size.  The same exact residue ledger now
leaves the old endpoints and one new all-low residue.  The old endpoints
again force the excluded mass \(p+9\).  In the new residue, subtracting
either exact parity baseline from a hard cell and dividing by two produces a
nonzero nonnegative integral quadratic \(B\) with

\[
 4p\,\mathbb E B=p-3.                                  \tag{15.769.2}
\]

Sharpness in Proposition 15.688 makes \(B\) Boolean.  The corrected Johnson
bound leaves at most five slice coordinates, and the cube influence bound
leaves at most four.  The fixed four-bit catalog has exactly ten tables at
this density: six omitted-pair tables and four all-equal-triple tables.

Combining those two lifts with the XNOR and complementary-literal baselines
gives coefficient offsets \(2,3,4,5\).  The common-row identity first makes
the hard parallel count \(P\) common; the slice congruence and edge bound
then force \(P\) to equal its offset.  In all four families,

\[
 hT=(p+1)P-5p+4,\qquad
 a(Q)=(p+1)Q+hT-3p.                                    \tag{15.769.3}
\]

The row \(Q=8-P\) has mass twelve and lies below both relevant floors.
Moreover

\[
 \sum_L\bigl(Q_L-(9-P)\bigr)=m-9,                      \tag{15.769.4}
\]

so at least nine opposite directions have \(Q=9-P\) and mass \(p+13\).
Any nonzero-boundary alternative again leaves a forbidden mass-twelve lift;
hence a survivor would give a nonzero nonnegative integral \(C\) with

\[
 4p\,\mathbb E C=p+13.                                 \tag{15.769.5}
\]

For \(H=\max C\ge2\), paired-cube averaging gives
\(H\ge(p-11)/4>3\), while the average paired-cube mean is below \(3/4\).
Mean \(1/4\) is excluded by support-floor equality, so quarter integrality
produces a half-mean cube, contradicting
\(\max g\le3\).  For \(H=1\), the Johnson/cube reduction reaches the fixed
four-bit catalog, whose density list omits \((p+13)/(4p)\).  This proves
(15.769.1).

There is one exceptional endpoint in the same layer:

\[
 p=23,\qquad t=9,\qquad k=110.                        \tag{15.769.6}
\]

Here the local mass \(p+13=36\) is attainable, so the preceding local lemma
cannot be extended.  Its equality case is nevertheless rigid.  Every
height-three half-mean paired-cube restriction is, up to dummy coordinates,

\[
 F_r(s)=3-2s+{s\choose2},\qquad r\in\{4,5\}.            \tag{15.769.7}
\]

Comparing overlapping paired cubes gives vanishing additive \(2\times2\)
minors, so (15.769.7) globalizes on \(J(23,12)\).  The hard and opposite
coefficient congruences leave only \(P=4,Q=5,F_5\).  The twelve distinct hard
directions are triangle-minus-full-star roots of the homogeneous forms

\[
 G_4=2hM_4-M_2^2,\qquad
 G_8=24hM_8-32M_2M_6+5M_2^4.                           \tag{15.769.8}
\]

Since their degrees are four and eight, both forms vanish identically.  An
opposite \(F_5\) row is a \(K_5\), so it would have to satisfy both resulting
coefficient equations.  The authoritative exact scan of all
\({23\choose5}=33{,}649\) five-sets finds \(1{,}518\) quartic zeros,
\(2{,}024\) octic zeros, and no simultaneous zero.  An independent
accelerator replay returns the same vector \([33649,1518,2024,0]\).
Therefore (15.769.6) is empty; this is a fixed coefficient certificate, not
a graph census.

## Proposition 15.770 — the next two post-band layers are empty

The newly classified equalities can be advanced by one quotient unit.  An
all-low family becomes \(m-1\) rows of parallel count \(P\) and one high row
of count \(P+1\).  The extra hard edge and the two extra total edges cancel
in the signed difference, so the common value \(hT\) is unchanged.  This
one-row carry is the only new global premise.

For \(p\equiv1\pmod4\), \(p\ge29\), it gives

\[
 t=q-2={p-5\over2},\qquad k=5p-5.                      \tag{15.770.1}
\]

The carried complement-triple family again forces mass \(p+15\).  The new
all-low XNOR residue differs from its baseline by a sharp mass-\(p-3\)
Boolean lift, so the same ten-table classification leaves the omitted-pair
and all-equal-triple offsets \(3,5\); both force mass \(p+13\).  The old
complement-triple boundary at this same mean has floor plus two. Its
exclusion requires the punctured-gap theorem recorded below with 15.772:
the difference from \((r-2)^2\) need not be nonnegative at \(r=0\), so
15.688 alone does not justify that subtraction. The repaired proof leaves
the stated 15.770 range unchanged. The old
endpoints force \(p+7\) or \(p+9\).  The masses \(p+15,p+7,p+9\) were
excluded above or in Propositions 15.751--15.752.  For the remaining
\(p+13\) mass, the height-at-least-two proof uses a half-mean cube for
\(p\ge37\); at \(p=29\), it bootstraps to height twelve and uses
(15.768.7).  The Boolean density \((p+13)/(4p)\) is absent from the fixed
catalog.  Thus (15.770.1)
is empty.

For \(p\equiv3\pmod4\), \(p\ge31\), the carry gives

\[
 t=q-1={p-3\over2},\qquad k=5p-3.                      \tag{15.770.2}
\]

The four carried sharp families again force mass \(p+13\).  The only new
residue is a lift of mass \(p-1\), below the height-at-least-two floor
\(p+1\), hence Boolean.  The Johnson/cube reduction reaches four active
coordinates, and the fixed catalog omits density \((p-1)/(4p)\).  The old
endpoints again force \(p+9\).  Hence (15.770.2) is empty.

Finally, Proposition 15.770 includes the exceptional endpoint

\[
 p=23,\qquad t=10,\qquad k=112.                       \tag{15.770.3}
\]

The new mass-22 residue is Boolean and its density \(11/46\) is absent from
the same four-bit catalog.  In the carried sharp residue there are eleven
low hard rows of count \(P\) and one high row of count \(P+1\).  After the
mass-twelve row is excluded, the opposite surplus is four, so at least eight
opposite rows have mass 36.  The globalized \(F_4/F_5\) classification and
offset compatibility again force \(P=4,Q=5,F_5\).  The eleven low hard
triangle-minus-star directions are eleven distinct roots of both forms in
(15.769.8); since \(11>8\), both forms vanish identically.  The high row is
not needed.  Every forced opposite \(K_5\) would contradict the same exact
33,649-five-set zero-intersection certificate.  This proves (15.770.3).

As of Proposition 15.770, with \(q=(p-1)/2\), the residual-(ii) frontier was

\[
\begin{gathered}
 p=5,7;\qquad p=11, t\ge3\ (k\ge50);\qquad
 p\in\{13,17,19\}, t\ge5;\qquad p=23, t\ge11;\\
 p\equiv1\pmod4, p\ge29, t\ge q-1;\qquad
 p\equiv3\pmod4, p\ge31, t\ge q,
\end{gathered}                                         \tag{15.770.4}
\]

together with the separately tracked positive \(p=7,z=7\) branch.
Residual (ii), E1, \(L=1/2\), and the original limit remain open.  Evidence
and fail-when-wrong replays are in `src/e1_gmin_m4_prop15768.py`,
`src/e1_gmin_m4_prop15769.py`, `src/e1_gmin_m4_prop15770.py`,
`tests/test_prop15768.py`, `tests/test_prop15769.py`,
`tests/test_prop15770.py`, `tests/test_p23_post_band_moment_close.py`,
`tests/test_p23_second_post_band_moment_close.py`,
`evidence/e1_gmin_m4_p23_post_band_moment_close.json`, and
`evidence/e1_gmin_m4_p23_second_post_band_moment_close.json`.

## Proposition 15.771 — the third p23 post-band endpoint is empty

At \(p=23\), residual (ii) is empty for every boundary size at

\[
 t=11,\qquad k=114.                                    \tag{15.771.1}
\]

The set \(H\) has 115 edges and at least 300 isolated vertices. In an
isolated chart, \(I=0\), \(m=12\), \(q=11\), and the exact hard ledger is

\[
 a_L=2u+24k_L,\qquad \sum_L k_L=23-u,\qquad
 hT=24P_L-69-a_L.                                      \tag{15.771.2}
\]

The phase-one floors and sharp mass-20 lift bound leave only
\(u=9,10,11\). At \(u=9\), the two-unit carry has at least ten low hard
roots; ten exceeds the octic degree, and the already proved quartic/octic
five-set certificate excludes the seven forced mass-36 opposite rows.
At \(u=10\), at least eleven low rows force the already excluded mass-22
lift. These are uses of the preceding endpoint theorems, not new censuses.

For \(u=11\), either an exact quotient-zero baseline occurs at \(b=2,22\),
or every hard quotient equals one and every hard row has mean 46.
In the latter case, positive degree-two quadrature forces the local
quadratic \(A\) to equal one on every even boundary-intersection layer.
At \(b=4,20\), the elementary fixed-weight quadratic kernel removes the
pure outside polynomial, and the pair-contact equations make each mixed
coefficient row constant. Thus the function depends only on the four-point
boundary or three-point complement. The complete equality types are
\(4000,2200\), of offset five, and \(000;4,200;2,220;0,400;0\), of
offsets eight, six, four, four, respectively.

At \(b=6,8,\ldots,18\), put \(d=\min(b,23-b)\ge5\). Every twelve-set
lies in a \(d\)-dimensional cross-boundary swap cube: match each small-side
coordinate to a distinct outside coordinate of opposite membership and
fix the other \(12-d\) selected points. The matching capacities have
slacks \(11-d\) and \(12-d\). On the cube, original boundary parity is
the bit-weight parity, including when the small side is the complement.
For degree-two Walsh characters the even-half Gram matrix is
\(2^{d-1}I\), because no two index sets of size at most two are
complements when \(d\ge5\). Hence \(A-1\) vanishes on the whole cube,
contradicting parity at an odd vertex. The remaining \(b=0,2,22\) cases
are constant one or the baseline plus a classified Boolean mass-24 lift.
The exhaustive signed-target offsets are therefore \(4,5,6,7,8\).

Equal means in (15.771.2) first make \(P\) common. Only then do the
coefficient congruence \(P\equiv\mathrm{offset}\pmod{11}\) and
\(12P\le115\) force \(P=\mathrm{offset}\). The opposite ledger is

\[
 \sum Q=115-12P,\qquad a(Q)=24(P+Q)-184.                \tag{15.771.3}
\]

The mass-eight value at \(Q=8-P\) is impossible, so \(Q\ge9-P\).
The surplus above these twelve minima is seven, forcing at least five
mass-32 rows. If a quotient-zero baseline occurs instead, the common
integer \(c=P_L-k_L\) is fixed to three or four, and the ledger is

\[
 \sum Q=103-12c,\qquad a(Q)=24(c+Q)-160.                \tag{15.771.4}
\]

It again has forbidden mass eight, surplus seven, and at least five
mass-32 rows. At phase zero the floors allow mass 32 only at
\(b=0,2,22\). At the latter two boundaries, subtracting the pointwise
Boolean parity minimum (XOR or the omitted bit) and dividing by two gives
a nonnegative integral quadratic of mass \(32-24=8<20\), impossible.
At \(b=0\), write \(A=2L\); the local Proposition 15.752 theorem
excludes \(4pE[L]=32=p+9\). This proves (15.771.1).

The full proof and independently reviewed local lemmas are in
`evidence/NOTE_2026-09-04_P23_THIRD_POST_BAND_CLOSE.md`,
`evidence/NOTE_2026-09-04_P23_SMALL_BOUNDARY_EQUALITY_PROOF.md`, and
`evidence/NOTE_2026-09-04_P23_MIDDLE_BOUNDARY_CUBE_PROOF.md`.
The exact source/test package is `src/e1_gmin_m4_prop15771.py` and
`tests/test_prop15771.py`; complementary four-node checks are recorded in
`evidence/p23_third_post_band_mesh_replay.json`. The latter confirm the
algebra and do not substitute for the proof. The p23 frontier in
(15.770.4) is superseded by \(t\ge12\), equivalently \(k\ge116\).
Every other frontier there was unchanged by 15.771; 15.772 below separately
advances the generic p1 frontier. Residual (ii), E1, \(L=1/2\),
and original convergence remain open.

## Proposition 15.772 — the third generic p1 post-band layer is empty

For every prime \(p\equiv1\pmod4\), \(p\ge29\), residual (ii) is empty at

\[
q={p-1\over2},\qquad t=q-1,\qquad k=5p-3.             \tag{15.772.1}
\]

Here \(|H|=5p-2\) and at least \(p^2-10p+5>0\) vertices are isolated.
Transport one to infinity, put \(m=q+1\), and choose the hard phase-one
sign \(h\). The exact common-row equations are

\[
a_L=2u+(p+1)k_L,\quad 0\le u\le q,\quad
\sum_{L\ {m hard}}k_L=2q-u,\quad
hT=(p+1)P_L-3p-a_L.                                  \tag{15.772.2}
\]

**Punctured complement-triple gap.** Let \(C\) be a three-set,
\(r=|X\cap C|\), and let \(A\) be a nonnegative integral quadratic of
parity \(r\) on \(J(p,m)\). If \(2pE[A]=2p-6+\delta\),
\(0\le\delta\le4\), then precisely the following possibilities occur:

\[
\begin{array}{c|c|c}
\delta&A&\text{signed coefficient offset}\\\hline
0&(r-2)^2&2\\
4&2-x_i-x_j-x_k+2x_ix_j,\quad \{i,j,k\}=C&4.
\end{array}                                         \tag{15.772.3}
\]

There are three labeled forms in the second line. In particular gap two
is impossible, but gap four must not be discarded. To justify the
punctured step, put \(L=(A-(r-2)^2)/2\). This is nonnegative on
\(r=1,2,3\), but only \(L\ge-2\) on \(r=0\).
Positive degree-two quadrature gives

\[
\delta=(p-3)\left(\sum_{|S|=1}\mu_S+\mu_C\right)
       +4\sum_{|S|=2}\mu_S,                         \tag{15.772.4}
\]

where \(\mu_S\) is the section mean with \(X\cap C=S\).
A nonzero integral nonnegative quadratic on an even neighboring slice
\(J(N,N/2\pm1)\) has mean at least
\((N-2)(N-4)/(4N(N-1))\): at most two one-sections (respectively
zero-sections) vanish, by the zero third difference on a three-swap cube,
and every other section has the odd middle-slice 15.688 floor.
At \(N=p-3\), its contribution to (15.772.4) is at least
\((p-5)(p-7)/(4(p-4))>4\). Thus all singleton and triple sections vanish
pointwise. The fixed-weight quadratic kernel and pair-contact sums then
remove every outside coordinate. On the remaining three-cube the pair
values \(a_{ij}\) are nonnegative integers, the zero-pattern value is
\(-\sum a_{ij}\), and \(\delta=4\sum a_{ij}\), proving (15.772.3).
This repairs the floor-plus-two premise used in 15.770 without assuming
global nonnegativity of the punctured difference.

For \(u<q\), the phase-one floors prohibit quotient zero, and at least
\(u+2\) rows have quotient one. At \(u=q\), some quotient is zero
because their sum is \(q<m\). The low means and their possible offsets
are exhaustively

\[
\begin{array}{c|c|c}
u& a_{\rm low}&\text{offsets}\\\hline
0&p+1&5\\
1,\ldots,q-4&p+1+2u&\text{none}\\
q-3&2p-6&2\\
q-2&2p-4&3,5\\
q-1&2p-2&4,6\\
q&p-1&4.
\end{array}                                         \tag{15.772.5}
\]

The pointwise XNOR and complement-literal baselines allow genuine
nonnegative lift subtraction. Sharp mass \(p-3\) is classified by the
fixed four-bit catalog; (15.772.3) handles the other complement-triple
possibilities. Mass \(p-1\) for \(p\equiv1\pmod4\) is excluded as
follows. At height at least two, paired-cube and stabilizer inequalities
force \(H=(p+3)/4\ge8\) and every maximizing cube to have mean one half,
contradicting its maximum bound three. At height one, the corrected
Johnson bound is less than six, reducing through cube influence to the
same fixed four-bit catalog, which omits density \((p-1)/(4p)\).
No classification of mean \(2p\) is used: its residue \(u=q\) already
has a quotient-zero XNOR row.

For a low row of quotient \(\ell\) and count \(P\), common \(hT\) first
gives \(P_L=P+k_L-\ell\) and total hard count
\(m(P-\ell)+2q-u\). Nonnegative opposite count bounds \(0\le P\le9\).
Only now may the congruence \(P\equiv\mathrm{offset}\pmod q\),
\(q\ge14\), give \(P=\mathrm{offset}\). Applying the opposite ledger
\(a_{\rm opp}(Q)=(p+1)Q+hT-3p\) gives

| Low family | Forbidden small mass | Forced next mass | Number of next rows at least |
|---|---:|---|---:|
| complement literal | 6 | `p+7` | 5 |
| complement triple | 14 | `p+15` | 9 |
| XNOR plus sharp lift | 12 | `p+13` | 8 |
| gap-four or literal-plus-sharp equality | 10 | `p+11` | 7 |
| quotient-zero XNOR | 8 | `p+9` | 6 |

The small masses lie below both the nonzero phase-zero boundary floor and
the sharp zero-boundary lift floor. At each next mass, only boundaries
\(b=0,2,p-1\) survive the phase-zero floors. Subtracting the pointwise
XOR or omitted-bit minimum excludes \(b=2,p-1\) by a positive lift mass
below \(p-3\). At \(b=0\), write \(A=2F\). The old local theorems
exclude masses \(p+7,p+9,p+13,p+15\).

For the new mass \(4pE[F]=p+11\), height \(H\ge2\) would give
\(H\ge(p-9)/4\ge5\) while the average maximizing-cube mean is at most
\((p+11)/(2(p-1))<3/4\). Quarter-integrality forces some such cube to
have mean one half and hence maximum at most three, a contradiction.
At height one, corrected Johnson influence is less than eight. The
resulting at-most-seven-coordinate slice junta extends to a cube, where
Boolean quadratic influence leaves at most four active variables. The
fixed catalog omits \((p+11)/(4p)\). This exhausts (15.772.5) and proves
(15.772.1).

The p1 frontier is now \(t\ge q\), matching the p3 frontier; the next
common generic layer \(t=q,k=5p-1\) has not been attacked here. The p23
frontier remains \(t\ge12\), and all small-prime and global open gates
remain unchanged. Residual (ii), E1, \(L=1/2\), and the original limit
remain open.

Full proofs: `evidence/NOTE_2026-09-04_P1_THIRD_POST_BAND_CLOSE.md` and
`evidence/NOTE_2026-09-04_COMPLEMENT_TRIPLE_PUNCTURED_GAP.md`.
Exact packages: `src/e1_gmin_m4_prop15772.py`,
`src/e1_gmin_m4_complement_triple_gap.py`, `tests/test_prop15772.py`, and
`tests/test_complement_triple_gap.py`. The four-node arithmetic record
`evidence/p1_third_post_band_mesh_replay.json` independently checks the
quadrature, contact kernel, equality tables, and row ledger; it is not a
substitute for the analytic proof.

**All-prime branch-C odd--Radon centrality (2026-09-03).**
Let \(p=4r+3\) be prime with \(r\ge7\). Suppose one opposite row on the
branch-C compact ray contains \(b\) arbitrarily labelled compact atoms and
\(r-1\) positive all-equal triangle atoms. If \(3b\le r+2\), its aggregate
signed edge chain \(C\) has coefficient \(\ell _1\) mass and edge occurrence
count at most

\[
 N=3(r+b-1)\le4r-1.
\]

Pair each nonantipodal edge with its negative and use invariant coordinates
\(U=(s+t)^2\) and \(D=(s-t)^2\). If all odd global forms through degree
\(p-2\) vanish, then

\[
 W(U,D)=n_E D(s+t)
\]

is orthogonal on the \((2r+1)\)-by-\((2r+1)\) square grid to every bivariate
polynomial of total degree at most \(2r-1\). The support-isolation theorem says
that a nonzero \(W\) must lie on a maximal line \(U=u_0\), \(D=d_0\), or
\(D=aU\). The edge-orbit differences satisfy

\[
 \lvert n_E\rvert\le (r-1)+2b<2r+1.
\]

The horizontal and diagonal lines require all \(2r+1\) projective nonzero
coefficient classes, so they are impossible. A vertical line instead gives a
constant nonzero integer \(k\) on its \(2r+1\) fixed-sum matching orbits. Since
\(N<2(2r+1)\), the occurrence bound forces \(\lvert k\rvert=1\). After fixing
the sign, each positive all-equal triangle supplies at most one aligned
fixed-sum edge and each compact triangle at most two: two all-equal pair sums
with the same sign, or all three compact signed pair sums, would repeat a
label. The total aligned capacity is therefore

\[
 (r-1)+2b<2r+1,
\]

contradicting the required matching. Hence \(W=0\) and \(C\) is centrally
symmetric.

For the deterministic balanced branch-C allocation, put

\[
 t_0=2r^2-4r-2,\qquad \delta=t-t_0=\sum_L b_L,
 \qquad m=2r+2.
\]

Every \(b_L\) is \(\lfloor\delta/m\rfloor\) or
\(\lceil\delta/m\rceil\). Consequently every opposite row is central whenever

\[
 0\le\delta\le(2r+2)\left\lfloor{r+2\over3}\right\rfloor.
\]

This is an all-prime structural theorem under the zero-odd-form and balanced
allocation hypotheses, not an even-moment exclusion. It does not treat
nonzero odd global forms, unbalanced allocations, coordinate the
degree-six/eight forms, or produce an integral or Boolean common-edge lift.

**Audited \(p=31,t=69\) arbitrary-compact local gate (2026-09-03).**
Consider one arbitrarily labelled compact atom

\[
 K(a,b;c)=\{a,b\}-\{a,c\}-\{b,c\}
\]

and six all-equal triangle atoms, allowing repeated blocks but requiring
three distinct labels within each block. If every odd row
\(d=3,5,\ldots,29\) vanishes, the total signed edge chain is centrally
symmetric: the all-prime theorem applies with \(r=7\), \(b=1\), since
\(3b=3\le9=r+2\). The earlier specialized calculation gives the same
conclusion by putting the resulting dual word on the 15-by-15 square grid;
the general aligned-incidence proof above supersedes its separate vertical
matching bookkeeping.

Centrality does not force \((a,b;c)\) to be centered, so the remaining joint
degree-six/eight problem was classified exactly. The 13,485 labelled compact
atoms form 450 scaling orbits: 449 free noncentered orbits and the unique
centered orbit represented by \((1,30;0)\) at index 435. A sparse exact DFS
tracks all 225 edge-orbit imbalances and all seven degree-six/eight channels.
Its necessary \(\ell_1\), maximum-coordinate, and final one/two-block moment
prunes are exact; the complete central remainder classification through
volume five consists of invariant/negation-pair cores and the Pasch
four-trade. The full run visits 317,916,856 states. Every noncentered orbit is
`UNSAT`, the centered orbit is covered by the preceding independent exact
certificate, and no `SAT` orbit exists. An independent v1/v2 comparison
agrees on all 435 overlapping fibers.

Therefore, at this one row, no arbitrary compact atom plus six all-equal
atoms can make the odd, degree-six, and degree-eight global forms all zero.
The stored parameter \(t=69\) records where that atom profile was first
encountered; the certificate itself depends only on the row.

For \(p=31\), one has \(r=7\), \(t_0=68\), and 16 balanced opposite
directions. The original all-prime theorem above makes every such row central
through \(t=116\); the component-packing theorem below supersedes this by the
guaranteed endpoint \(t=164\). Moreover, exactly when
\(1\le t-68\le31\), at least one balanced
row has compact count one. It therefore has precisely the certified
one-compact/six-all-equal profile. Hence no balanced profile with

\[
 69\le t\le99
\]

can make all odd, degree-six, and degree-eight global forms zero. This is a
balanced zero-global-form band exclusion obtained from a local row
certificate. It does not cover unbalanced allocations, exclude nonzero forms
coordinated across directions, decide the signed Boolean box (15.760.10),
close residual (ii), or prove the original limit.

Evidence and replay are in
`src/e1_gmin_m4_compact_ray_moment_gate.py`,
`tests/test_compact_ray_moment_gate.py`,
`evidence/NOTE_2026-09-02_COMPACT_RAY_HIGHER_MOMENT_GATE.md`,
`evidence/p31_arbitrary_compact_fiber.cpp`,
`evidence/p31_arbitrary_compact_fiber_v2_000_450.log`, and
`evidence/p31_arbitrary_compact_fiber_v2_merge.json`. The pinned SHA-256
values for source, raw log, normalized verdicts, and corrected merge record
are respectively
`1dcfce7b5765630655d049413c4d9138c544a6d05fe19e3308a9a20a2880d1f2`,
`f3f77607181287095aa69644649d14d7b9b5e3a8f24044477b667549ef0512e3`,
`ad3bf3c97b378c9cdebb0b77d486cced544199750ad689060bd2a24f6a2210cb`,
and `c7f5dea5811a8d2aa25d7bd3224b1fceae3fce73bb49fd4c8fe3f335e2e71c2f`.

**Seven-channel algebraic-dominance barrier (2026-09-03).**
For each of the compact and all-equal atom types, sum four atoms with triples
\((a_i,b_i,0)\) and map them to the seven channels
\((6,0),(6,1),(6,2),(8,0),(8,1),(8,2),(8,3)\). At

\[
 ((a_1,b_1),\ldots,(a_4,b_4))=((2,1),(3,2),(4,3),(5,4)),
\]

the Jacobian minor in \((b_1,a_2,b_2,a_3,b_3,a_4,b_4)\) is respectively

\[
\begin{split}
 226534996574208000&=2^{28}3^9 5^3 7^3 &&\text{(compact)},\\
 220242357780480000&=2^{26}3^7 5^4 7^4 &&\text{(all-equal)}.
\end{split}
\]

Both maps are therefore dominant on the distinct-label locus over the
algebraic closure in every characteristic at least 11. Every balanced
branch-C hard row for \(r\ge7\) has at least four compact atoms and every
opposite row has at least six all-equal atoms, so freezing extra atoms
preserves dominance. Dense-open intersection and rescaling by a linear form
with no zero on \(\mathbf P^1(\mathbf F_p)\) then produce common nonzero
degree-six/eight forms over an algebraic closure, and already over some finite
extension.

Consequently there is no universal polynomial identity among the seven
channels and no purely algebraic projective root-count contradiction. This
does **not** provide labels or form coefficients in \(\mathbf F_p\), and hence
does not construct an admissible finite-field row family. It also proves no
compatibility with nonzero odd or higher moments and no signed Boolean lift.
The finite-field rational, unbalanced zero-form, and globally coordinated
gates remain open.

**Full-balanced maximal-line and boundary-cubic exclusions (2026-09-03).**
Let

\[
 p=4r+3,\qquad r\ge7,\qquad h=2r+1,\qquad m=h-2,
\]

and let an opposite balanced row contain $r-1$ positive all-equal atoms
and $0\le b\le r$ compact atoms.  Its signed occurrence and coefficient
budget is

\[
                         N=3(r+b-1)\le3h-6.              \tag{G.1}
\]

Suppose its odd-Radon word is supported on one maximal line of the square
grid $H\times H$.  On a horizontal line or diagonal line, the unique dual
relation runs through every nonzero projective residue class, and therefore

\[
                   \sum_E|n_E|\ge {h(h+1)\over2}>3h-6. \tag{G.2}
\]

On a vertical line all $h$ orbit differences are congruent to one
nonzero residue.  The budget in (G.1) forces the actual integer lifts
to be one constant $k$ with $|k|\in\{1,2\}$.  If $|k|=1$, reducing
modulo two and identifying $x\sim-x$ leaves exactly two odd vertices in
the fixed-sum graph, whereas every atom projects to a triangle boundary.  If
$|k|=2$, the required $2h$ aligned occurrences exceed the capacity

\[
                         (r-1)+2b\le3r-1<2h.             \tag{G.3}
\]

Thus no one-maximal-line support occurs anywhere in the full balanced range.

Couvreur peeling makes the corresponding two-line statement exhaustive.
If a degree-$m$ dual word of support at most $3m$ contains $h=m+2$
points on one line, the remainder is either empty or consists of between
$h-1$ and $h$ points on a second maximal line.  For lines in different
families, the dual space on their union is spanned by the two line relations,
and one component already costs

\[
                         {h(h-1)\over2}>3h-6.             \tag{G.4}
\]

For two vertical lines, either one coefficient fiber is injective and costs
$(r+1)^2>3h-6$, or both coefficients are constant units and quotient
parity again leaves two odd vertices.  For two horizontal lines or two
diagonals, each projective coefficient class occurs at most three times; the
least possible mass

\[
 L_3(h)=3{q(q+1)\over2}+s(q+1)>3h-6,
 \qquad h=3q+s,\quad 0\le s<3,                            \tag{G.5}
\]

is too large.  Hence no realizable word with support at most $3h-6$
contains $h$ collinear points.

At the exact boundary $b=r$, the remaining minimally linked cubic case is
a reduced complete intersection of a cubic $F$ and a degree-$m$ curve,
with $3m=3h-6$ support points.  Equality in the support and
$\ell _1$ budgets makes every orbit difference a unit.  A reducible cubic
is excluded by restricting the Euler--Jacobi residue relation to its line
and conic components; the only apparent vertical-line/tangent-conic case
would force a nonconstant omitted parameter to equal the tangent parameter.
A geometrically integral singular cubic has at most $p+2<3h-6$ rational
points.  For a smooth cubic, adjoining $y^2=U$ gives a connected double
cover of genus at most four.  The support would give at least $3p-15$
rational points on the cover, while Weil gives at most $p+1+8\sqrt p$, and

\[
                         3p-15>p+1+8\sqrt p               \tag{G.6}
\]

for every $p\ge31$.  Therefore every boundary cubic support is excluded
for $p=4r+3\ge31$.  These are support exclusions, not an extension of the
earlier $3b\le r+2$ centrality theorem: the high-intersection irreducible
conic alternative still has to be treated separately.

**Irreducible-conic peeling, classification, and equianharmonic survivor
(2026-09-03).**  Suppose the support contains at least
$2m+2=p-3$ points on an irreducible conic $Q$.  Multiplying the odd-Radon
relation by $Q$ shows that any nonempty remainder would be a degree-$(m-2)$
dual word of support at least $m$, whereas the budget leaves at most
$m-2$ points.  The whole word is therefore supported on $Q$.

The square-coordinate character bounds then force every such
high-intersection conic into the triangle-tangent normal form

\[
 U=uz^2,\qquad D=d(z-1)^2,\qquad u,d\in H,
 \qquad z\in\mathbf F_p\setminus\{0,1\}.                  \tag{C.1}
\]

The restricted dual Reed--Solomon relation is

\[
 W(z)=c z(z-1)(Az+B),\qquad
 n(z)=\alpha+{\beta\over z-1}.                            \tag{C.2}
\]

If $\beta\ne0$, its least-absolute-residue mass exceeds $3m$; hence
$\beta=0$, and the actual orbit differences are one constant unit.  After
scaling $u=1$, write $d=k^2$ and

\[
                         q={1-k\over1+k}.                  \tag{C.3}
\]

For $k=\pm1$, quotient parity leaves two odd vertices, so the star is
impossible.  In the nonstar case, exact alignment scoring shows that an
all-equal atom has score three only on a three-cycle.  If $q^3\ne1$, only
two compact score-three types can occur, and the total aligned/reverse/off
count cannot reach the $p-2$ target coordinates within $b\le r$.
Consequently every nonstar survivor satisfies

\[
 q^3=1,\quad q\ne1
 \quad\Longleftrightarrow\quad k^2=-3,
 \qquad p\equiv7\pmod {12}.                               \tag{C.4}
\]

This survivor is real.  At $p=31,r=b=7,k=11$, take the six positive
all-equal triples

\[
 (1,3,13),(2,7,8),(4,18,26),(5,20,23),(6,14,28),(24,25,30)
\]

and the seven compact atoms, written as $(\{a,b,c\};c_*)$,

\[
\begin{split}
 &(\{0,2,12\};2),(\{0,12,19\};12),(\{4,19,22\};4),\\
 &(\{9,12,16\};16),(\{10,14,17\};14),\\
 &(\{11,19,22\};19),(\{12,16,20\};20).
\end{split}                                               \tag{C.5}
\]

Direct integer replay gives coefficient $+1$ on exactly the 29
constant-conic edges and zero elsewhere.  All 105 odd channels of degrees
$3,5,\ldots,29$ vanish, but

\[
 F_6=(11,19,10),\qquad F_8=(12,11,23,6).                  \tag{C.6}
\]

Thus (C.5) is an exact odd-zero noncentral edge witness, not a
simultaneous zero-six/eight witness and not a common-graph construction.

**Exact $p=31,b=7,k=11$ zero-six/eight finite-fiber exclusion
(2026-09-03).**  In the constant equianharmonic conic fiber above, every
all-equal atom has alignment score at most three and every compact atom has
score at most two.  The target has squared norm 29, while the six/seven-atom
maximum is $6\cdot3+7\cdot2=32$.  Any realization therefore has total
deficit exactly three, in exactly one of the disjoint partitions

\[
                         (3),\qquad(2,1),\qquad(1,1,1).    \tag{M.1}
\]

An independent integer meet-in-the-middle/exact-cover generator exhausts
these cases, including repeated atoms.  The nine maximal all-equal atoms are
the disjoint $q$-cycle triangles and their multiplicities are recovered
uniquely from exact edge equality.  The 111 maximal compact atoms are
completed by their exact 51-off-orbit parameterization.  The exhaustive
ledger is

\[
\begin{array}{c|r|r|r}
\text{deficit case}&\text{maximal completions}&\text{edge hits}&F_6=F_8=0\\ \hline
(3)&13{,}528{,}344&60&0\\
(2,1)&87{,}840{,}508&2{,}160&0\\
(1,1,1),\ \text{at least one AE}&20{,}465{,}801&392&0\\
(1,1,1),\ \text{all compact}&108{,}480{,}057&14{,}464&0\\ \hline
\text{total}&230{,}314{,}710&17{,}076&0.
\end{array}                                               \tag{M.2}
\]

Hence no choice of six positive all-equal atoms and seven compact atoms in
the $p=31,b=7,k=11$ constant-conic fiber realizes the target while all
degree-six and degree-eight channels vanish.  This is an exhaustive finite
fiber theorem and nothing larger; in particular it does not exclude other
equianharmonic primes or parameters and does not close residual (ii).

**Exact signed-Boolean defect, complete-Graver criterion, and ridge quotient
(2026-09-03).**  Retain the signed source convention
$z_e\in\{0,\tau_e\}$.  The edge-sign row is an integral combination of the
parallel Radon rows, so on every integral fiber $Rz=y$,

\[
 H_y=\sum_L\epsilon_LP_L(y)=\tau\mathbin\cdot z           \tag{B.1}
\]

is fixed.  Define

\[
 \beta(z)={1\over2}\sum_ez_e(z_e-\tau_e)
          ={\|z\|_2^2-H_y\over2},\qquad
 \beta_R(y)=\min_{Rz=y}\beta(z).                          \tag{B.2}
\]

Each summand is a nonnegative integer, with equality exactly at
$z_e\in\{0,\tau_e\}$.  Therefore

\[
 (z_0+\ker_{\mathbf Z}R)\cap\prod_e\{0,\tau_e\}\ne\varnothing
 \quad\Longleftrightarrow\quad\beta_R(y)=0.               \tag{B.3}
\]

For a kernel move $g$,

\[
 \beta(z+g)-\beta(z)=z\mathbin\cdot g+{\|g\|_2^2\over2}.
                                                                    \tag{B.4}
\]

Conformal decomposition into the complete Graver basis
$\mathcal G(R)$ makes global optimality equivalent to the exact Voronoi
criterion

\[
 z\text{ minimizes }\beta
 \quad\Longleftrightarrow\quad
 |2z\cdot g|\le\|g\|_2^2\quad(g\in\mathcal G(R)).         \tag{B.5}
\]

This is a finite mathematical alternative: descent terminates either at
$\beta=0$, producing a common graph, or at positive defect certified by
the complete inequalities.  It is not an assertion that the complete
Graver basis has been constructed.

There are explicit linearly supported kernel circuits.  Type P moves compare
two midpoint fibers in one parallel difference class; Type K moves compare
two difference classes with the same nonzero projected square.  If
$K_{\rm ridge}$ is their integer span, then for every odd prime
$p$,

\[
 p\ker_{\mathbf Z}R\subseteq{\cal K}_{\rm ridge}
 \subseteq\ker_{\mathbf Z}R.                              \tag{B.6}
\]

Writing $m=(p-1)/2$ and $d=p+1$, the exact quotient is

\[
 \boxed{
 \ker_{\mathbf Z}R/{\cal K}_{\rm ridge}
 \cong(\mathbf Z/p\mathbf Z)^{\nu_p},\qquad
 \nu_p=dpm^2+{m(m-1)(4m+1)\over6}.}                       \tag{B.7}
\]

A basis of the mod-$p$ dependencies of the ridge synthesis matrix gives
the missing saturating moves and, together with the ridge basis, an exact
parametrization of the full integral fiber.  At a defect minimizer the
elementary ridges give the necessary inequalities

\[
 |S_\delta(\alpha)-S_\delta(\alpha')|\le p,               \tag{B.8}
\]

and, for same-square transverse classes,

\[
 |S_{\delta_1}(\alpha)-S_{\delta_1}(\alpha')
  -S_{\delta_2}(\alpha)+S_{\delta_2}(\alpha')|\le2p.      \tag{B.9}
\]

The quotient in (B.7) is proper and forces at least $2\nu_p$
non-ridge elements in any symmetric complete Graver basis.  Thus the ridge
system is an exact one-step saturation and descent reduction, not the
complete Graver system and not a signed Boolean lift.

**Equianharmonic component-packing upgrade (2026-09-03).**  Continue with
\(p=4r+3\ge31\) and a zero-odd opposite row having \(r-1\) all-equal
atoms and \(b\le r\) compact atoms. In the only unresolved tangent-conic
branch, put

\[
 q^2+q+1=0,\qquad \Phi(x)=qx+1-q.
\]

Pair reverse target occurrences with aligned occurrences, pair non-target
occurrences antipodally, and treat self-antipodal occurrences as caps. For
one connected pairing component let \(K,E,\mu,Z\) denote its compact
vertices, all-equal vertices, cycle rank, and caps. Its exact deficit and
excess are

\[
 \delta=K+2E-2+2\mu+Z,\qquad
 K-2\delta=4-K-4E-4\mu-2Z.                    \tag{E.1}
\]

The complete symbolic classification excludes a distinct-label compact atom
of score three and leaves only a capped singleton (excess one), an HH pair
(excess two), and an O2/R1 three-atom tree (excess one). Their exact unpaired
target supports have weighted pairwise-disjoint packing number three. If
\(C=(p-1)/3\) is the number of nonfixed \(\Phi\)-cycles and
\(L=C-(r-1)=(r+5)/3\) are uncovered, the global deficit identity forces

\[
 \boxed{b\ge2L-1={2r+7\over3}.}               \tag{E.2}
\]

The map \(x\mapsto-x\) negates the chain in a fixed orbit basis, so (E.2)
covers both constant signs. Combining (E.2) with the proved line,
boundary-cubic, and conic dichotomies gives

\[
 \begin{array}{ll}
 p\equiv7\pmod {12}:&3b\le2r+4\Longrightarrow\text{central},\\
 p\equiv11\pmod {12}:&0\le b\le r\Longrightarrow\text{central}.
 \end{array}                                               \tag{E.3}
\]

For \(p=31\), (E.3) guarantees that every balanced opposite row is central
for \(68\le t\le164\). Rows with \(b=7\) can first occur for
\(165\le t\le177\), so the word "guarantees" is essential. At \(p=43\),
an exact nine-all-equal/nine-compact witness attains (E.2), realizes all 41
target coordinates, and kills all 210 odd channels; its nonzero syndromes
\(F_6=(37,19,8)\), \(F_8=(18,17,10,32)\) show only that the threshold is
sharp for the odd fiber.

**Hard compact residual and the closed antisymmetric Boolean half
(2026-09-03).** A balanced hard row is one fixed unit star plus
\(0\le e\le2r-2\) compact atoms. The star is invisible to every moment
through degree \(p-2\). Couvreur peeling, the maximal-line exclusions, and
the strict occurrence bound \(3e<3(h-2)\) eliminate every possible nonzero
odd-Radon support. Hence, under zero odd global forms, the compact residual
is centrally symmetric over the integers. The whole hard row is not
central: the fixed unit-star difference remains.

Let central inversion act by \(J(e)=-e\) on source edges and by \(I\) on
directional cells. With \(p=2h+1\) and \(d=p+1\), the antisymmetric block has

\[
 \operatorname {rank}E^-=d^2h^2,\quad
 \operatorname {rank}{\cal A}^-=dh^2,\quad
 \operatorname {rank}\ker R^-=dph^2,                    \tag{E.4}
\]

and its exact integral cokernel is

\[
 {\cal A}^-/R(E^-)\cong
 (\mathbf Z/p\mathbf Z)^{h(h-1)(h+1)/3},                 \tag{E.5}
\]

precisely the odd-moment rows. There is also a constructive ternary lift.
For a hard functional \(L\), an independent \(M\), and center \(j\ne0\), set

\[
 f(t)={t\over t+1},\quad
 (L,M)(u_t)=j(1,f(t)),\quad (L,M)(v_t)=j(t,t),
 \quad t\in\mathbf F_p\setminus\{-1\}.                  \tag{E.6}
\]

If \(E_{L,M,j}=\{\{u_t,v_t\}:t\ne-1\}\) and
\(z={\bf1}_{-E}-{\bf1}_E\), direct projective involutions give

\[
 R_Lz=S_{-j}-S_j,\qquad R_Nz=0\quad(N\ne L).             \tag{E.7}
\]

Each trade uses exactly \(p-1\) nonfixed inversion orbits. One prior trade
forbids at most \(2p-1\) auxiliary functionals; before the last of at most
\((p+1)/2\) hard trades, the exact union bound still leaves
\((p-1)/2>0\) choices. Thus arbitrary hard centers have pairwise-disjoint
ternary lifts. Restoring the Paley column signs preserves this construction.
The antisymmetric Boolean half is therefore **closed** for the balanced
zero-odd target.

This does not construct a graph. On each used nonfixed orbit the binary pair
total is forced to \(s_e=1\); on an unused nonfixed orbit it is
\(s_e\in\{0,2\}\); and each fixed antipodal edge is an independent binary
variable. Matching those coupled symmetric totals to the full target is the
remaining Boolean gate.

The preceding hard-star support theorem remains useful only as a checked
barrier. It proves the one-active-row ternary floor \(2(p-2)\) and classifies
floor equality by a projective pencil when at least nine hard centers are
active, but its total-edge, fixed-capacity, parallel-count, and Euclidean-norm
ledgers are all feasible. Equation (E.7) supersedes attempts to extract an
antisymmetric contradiction from those scalar diagnostics.

**Equianharmonic threshold even-syndrome barrier (2026-09-03).** At equality
in (E.2), \(b=2L-1\) and the global component excess is exactly one. The only
positive/negative excess totals are \((1,0),(2,1),(3,2)\). Two exact
four-compact blocks \(U(x)\) and \(V(x)\) replace the two all-equal cycles

\[
 A=(\Phi(x),\Phi^2(x),x),\qquad
 B=(\Phi(-x),\Phi^2(-x),-x)                              \tag{E.8}
\]

in every odd edge-orbit channel. For the seven degree-six/eight deviations,
four \(U\)-parameters \(-1,-2,1,0\) and three \(V\)-parameters \(-3,-2,-1\)
give the exact Jacobian

\[
 \det J=4128623683475967290061619200
       =2^{32}3^{26}5^2\,7\,2161.                       \tag{E.9}
\]

Hence the mixed syndrome map is dominant in characteristic zero and this
specialization stays nondegenerate away from the displayed prime factors.
This is a rigorous barrier: component excess alone supplies no universal
affine degree-six/eight obstruction once both trade families occur. It does
not construct a rational finite-field zero, prove simultaneous
\(F_6=F_8=0\), or provide a common or Boolean lift.

**The inversion-symmetric lattice and corrected Mobius capacity
(2026-09-03).** Let

\[
 E^+=\ker_{\mathbf Z}(1-J),\qquad
 {\cal A}^+={\cal A}\cap\ker_{\mathbf Z}(1-I).
\]

The exact symmetric ranks are

\[
 \operatorname {rank}E^+={dh(p^2+1)\over2},\qquad
 \operatorname {rank}{\cal A}^+=dh(h+1),\qquad
 \operatorname {rank}\ker R^+=dph^2.                    \tag{E.10}
\]

Modulo two, the fixed antipodal source edges inject onto the \(dh\)-
dimensional compatible fixed-coordinate space, while the nonfixed source
pairs surject onto the \(dh^2\) paired coordinates. Their target supports
are disjoint, so the full symmetric map is mod-two surjective. The natural
symmetric-cokernel map into the full cokernel is injective: its kernel is
killed by two, while mod-two surjectivity removes two-torsion. Its image is
the even-moment eigenspace. Consequently

\[
 {\cal A}^+/R(E^+)\cong
 (\mathbf Z/p\mathbf Z)^{S_+(p)},\qquad
 S_+(p)={(h-1)(2h^2+5h+6)\over6}.                       \tag{E.11}
\]

Thus zero even moments characterize an unrestricted integral central lift,
and compatible central targets always have a mod-two lift. Neither fact
places the same lift in the restricted Boolean box.

For the Mobius half \(E=\{e_t:t\ne-1\}\) of (E.6), direct projection gives

\[
 P_N(E)=
 \begin{cases}
 1,&N=L,\\
 1,&N=M,\\
 0,&N=L-M,\\
 1+\eta(1+m),&N=L+mM,\quad m\ne0,-1.
 \end{cases}                                             \tag{E.12}
\]

Every transverse row other than \(L\) is central. In a primal basis dual to
\(L,M\), its exact Paley column signs are

\[
 \tau_t=\eta\!\left(Q(e_1-t^2(e_1+e_2))\right),\qquad
 \sum_{t\ne-1}\tau_t=
 \sum_{t\in\mathbf F_p}\eta\!\left(Q(e_1-t^2(e_1+e_2))\right)-\epsilon_L.
                                                               \tag{E.13}
\]

The forced symmetric pair-total chain has parallel vector \(2P_N(E)\).
After an actual antisymmetric selection \(q_U\), its exact remaining central
target is

\[
 T_U=Y-Rq_U={Y+IY-RC_U\over2}.                            \tag{E.14}
\]

Used nonfixed orbits are frozen at zero after subtraction; unused nonfixed
orbits allow \(0\) or \(\tau_O(e+Je)\), and fixed antipodal edges allow
\(0\) or \(\tau_f f\). This restricted affine-box intersection remains
open despite (E.11).

There is an essential overlap correction. For any two distinct hard
directions and arbitrary nonzero centers, one can choose their auxiliaries so
the two localized trades share exactly one origin orbit with opposite signs.
Their sum stays ternary and has

\[
                         2(p-1)-2                         \tag{E.15}
\]

nonzero inversion orbits. Hence the disjoint construction's upper-endpoint
count \(m(p-1)=|H|_{\max}+1\) is not a universal lower bound and supplies no
one-edge capacity contradiction.

**All-active support sharpening (2026-09-03).** Assume additionally that
all \(m=(p+1)/2\) hard-star centers are nonzero, and let \(c\) count source
inversion orbits occupied on exactly one side. The general floor is
\(c\ge p-2\). At equality, the proved pencil theorem gives a common
projective endpoint \([P]\) and \(j_L^2=L(P)^2\) in every hard row. If
\(Q\) is the set of the \(p-2\) other endpoints, hard-row bijectivity gives

\[
 L(Q)=\mathbf F_p\setminus\{L(P),-L(P)\}.                \tag{E.16}
\]

Adjoining \(P,-P\) makes a \(p\)-point set determining at most
\((p+1)/2\) affine directions. The prime-order Redei--Megyesi theorem says a
noncollinear \(p\)-set determines at least \((p+3)/2\), so this set is the
line through \(P,-P\). All \(p-2\) pencil edges are then parallel in its
annihilating opposite row, whose balanced full-ray quota is at most
\(2r+2<p-2\). This contradiction proves

\[
 \boxed{c\ge p-1\quad\text{when every hard center is nonzero}.} \tag{E.17}
\]

The all-active hypothesis is indispensable to the proof. Equation (E.17)
is a support bound, not a restricted central Boolean lift or a residual-(ii)
closure.

**Fixed-edge elimination and the halved symmetric code (2026-09-03).**
Let (U) now denote the actual support of a chosen ternary antisymmetric
lift, after all cancellations.  In bases consisting of fixed antipodal
edges and nonfixed orbit sums on the source, and fixed cells and nonfixed
cell pairs on the target, the restricted symmetric map has block form

\[
                 R^+=\begin{pmatrix}A&2B\\0&C\end{pmatrix}.       \tag{E.18}
\]

Modulo two, (A) is an isomorphism onto the compatible fixed-cell residue
space.  If (g_L(0)) and (g_L(\beta)) are the fixed-cell bits and (L_v)
annihilates (v), its inverse is

\[
 a_{[v]}=g_{L_v}(0)+\sum_L g_L\bigl(L(v)^2\bigr)\pmod2.           \tag{E.19}
\]

Consequently the fixed binary vector (a(T_U)) is uniquely forced, not
relaxed.  Subtracting it and dividing only the fixed block gives the exact
equivalence

\[
 \sum_{O\notin U}b_O\widehat B_O=\widehat T_U,qquad
 b_O\in\{0,1\},qquad
 2\sum_{O\notin U}b_O=|H|-|U|-|a(T_U)|.                          \tag{E.20}
\]

For (p\equiv3\pmod4), the parallel coordinate in each direction (L)
also pins the integer slice

\[
 n_L={P_L-u_L-f_L\over2},qquad 0\le n_L\le dh^2-u_L.            \tag{E.21}
\]

The forced word has a geometric formula.  Write a nonfixed orbit as
({a-\delta,a+\delta}), modulo central inversion.  Then

\[
 \Phi(a,[\delta])=
 \begin{cases}
 0,&a\parallel\delta,\\
 \mathbf1_{\{[\delta+ca]:c\in\mathbf F_p\}},&a\not\parallel\delta,
 \end{cases}qquad
 a(T_U)=a_Y+\sum_{O\in U}\Phi(O).                               \tag{E.22}
\]

The nonzero words in (E.22) are the (dh=|\Delta|) antipodal pairs of
non-origin affine lines.  Their square incidence matrix satisfies
(MM^{\mathsf T}=M^{\mathsf T}M=I) over (mathbf F_2).  Moreover every
block has (h) disjoint (p)-column lifts in (ker C).  It follows that
the full, unpunctured halved map

\[
 D=(C,\Phi):\mathbf F_2^{\Delta^2}longrightarrow
 \mathbf F_2^{dh^2}\oplus\mathbf F_2^\Delta
 \quad\hbox{is onto},\qquad \operatorname {rank}D=dh(h+1).       \tag{E.23}
\]

This does not survive arbitrary puncturing.  For a projective (L) and a
nonzero square (\beta), the row word

\[
 X_{L,\beta}=\{([a],[\delta]):L(a)=0, L(\delta)^2=\beta\}        \tag{E.24}
\]

has weight (ph=|\Delta|-h\); deleting it drops rank.  Thus robustness under
all deletions of size at most (|\Delta|) is false.  The exact criterion for
the actual puncture is that no nonzero word of ({\rm Row}(D)) have support
contained in (U), or, for one target, that every resulting left-kernel
functional annihilate (widehat T_U).  One Möbius half has at most two
midpoints in any spatial direction, so (m=(p+1)/2) halves meet (E.24) in
at most (p+1<ph) columns for (p\ge7).  This excludes only the displayed
rectangle.  There is, however, a sharper projection reduction.  For every
non-origin difference block (C=B_{K,\beta}), put

\[
 r_C([a])=\sum_{[\delta]\in C}w([a],[\delta])\pmod2.
\]

The identity (M^{\mathsf T}M=I) shows that (r_C) is a union of cells in
the partition (A_K=\{K(a)=0\}), (B_{K,\alpha}=\{K(a)^2=\alpha\}), and
that all the (r_C) determine (w).  If every (r_C) uses only (A_K),
then (w) is a disjoint sum of the already excluded rectangles (E.24).
Otherwise, containment ({\rm supp}(w)\subseteq U) forces some full affine
midpoint block (B_{K,\alpha}) over one common (C) to be met by all
(h+1) Möbius halves, with at least (h) halves supplying two midpoint
classes.  This is the nearly saturated cover condition isolated before the
row-code gap theorem below.  That theorem now excludes it for every actual
Hamming-extendable branch-C puncture; it is no longer an open linear gate.

There is also an exact local integer kernel.  On a fixed difference slice
([\delta]), every kernel vector has

\[
 v([a],[\delta])=\gamma_{L_\delta(a)^2},qquad
 \gamma_0=0,qquad \sum_{s\ne0}\gamma_s=0.                       \tag{E.25}
\]

It is the root lattice (A_{h-1}); its primitive moves exchange two whole
(p)-element slabs and preserve every direction weight.  These moves
connect two binary solutions that already agree off that slice, but do not
prove global connectivity, normality, or existence.

Finally, if (N=m(p-1)) is the disjoint Möbius support count and
(|U|=N-2\kappa), then throughout the all-active branch-C ray

\[
 |H|=N-\{2(t_{\max}-t)+1\},qquad
 \kappa\ge t_{\max}-t+1                                      \tag{E.26}
\]

is necessary even before (E.20).  At (p=31) the demand ranges from 110
cancellations at (t=68) to one at (t=177), so the deliberately disjoint
lift is extendable nowhere on that ray.  Two halves share at most two
oppositely oriented orbits, and this bound is sharp, but the sharp
two-cancellation locus is rigid:
(q=r=1/2, A=B=3/2).  It supplies no free parameter for a greedy
multi-pair construction.  The exact first support test is therefore

\[
 \boxed{|U|+\left|a_Y+\sum_{O\in U}\Phi(O)\right|\le |H|},       \tag{E.27}
\]

followed by the punctured divided Boolean fibre (E.20).  Pairwise overlap
counts alone control neither test.

There are four further exact qualifications to this gate.  First, let
\(B_A\) be the span of the \(h\) paired non-origin affine blocks in
direction \(A\).  The orthogonal block incidence matrix decomposes
\(\mathbf F_2^\Delta=\bigoplus_A B_A\), and direct conversion of the raw
halved rows gives

\[
 {\rm Row}(D)=
 (\langle{\bf1}\rangle\otimes\mathbf F_2^\Delta)
 \mathbin{\dot+}\bigoplus_A(B_A\otimes B_A).                    \tag{E.28}
\]

Thus, after the \(M\)-basis change, the exact minimum-support problem is a
matrix \({\bf1}q^{\mathsf T}+T\), with \(T\) block diagonal by direction.
For every nonzero point word \(x\), incidence counting proves the sharp
one-dimensional branch inequality

\[
              {\rm wt}(x)+{\rm wt}(M^{\mathsf T}x)\ge p+1.     \tag{E.29}
\]

Points and affine blocks attain equality.  Equation (E.29) alone does not
sum in the form needed for the matrix in (E.28).  Moreover,
\({\bf1}\otimes e_\delta\) and the scalar graphs
\(\{([a],[\delta]):[a]=[c\delta]\}\) are nonrectangle row words of weight
\(|\Delta|\).  Hence a rectangles-only classification through
\(|\Delta|\) is false.

The sharper group-support inequality is now proved for every odd prime.  Write
a row-code word in the grouped form
\[
 W=\sum_{A,j}c_{A,j}\otimes b_{A,j},\qquad
 c_{A,j}\in\langle{\bf1}\rangle+B_A,\qquad
 S_A=\{a:(c_{A,j}(a))_j\ne0\},
\]
where the common \({\bf1}\)-coefficients come from the unique block
expansion of \(q\), and let \(k\) be the number of nonempty \(S_A\).
Each \(S_A\) is a union of cells in the partition consisting of the radial
line of size \(h\) and the \(h\) affine blocks of size \(p\) in direction
\(A\).  Every nonzero point word \(f\) satisfies

\[
 {\rm wt}(f)+
 \#\{A:(M^{\mathsf T}f)|_{B_A}\ne0\}\ge p+1.                    \tag{E.29a}
\]

For odd support this follows from radial-fibre parity.  For even support
\(s=2n\), choose distinct antipodal representatives \(v_i\) and form
\[
 P_u(X)=\prod_{i=1}^{2n}(X-u(v_i)^2).
\]
The unique monic degree-\(n\) top-half square root leaves homogeneous
remainders of degree at most \(2s\), at least one nonzero because the factors
of \(P_u\) are distinct in \(\mathbf F_p[U,V,X]\).  At every silent
projective direction the roots pair, so every remainder has a double zero.
Equivalently, the canonical remainder vanishes to order at least two there.
Projective root counting gives \(2z\le2s\) for the number \(z\) of silent
directions, proving (E.29a).

It follows that \({\rm wt}(W)\ge ph\).  Indeed, unless one \(S_A=\Delta\), distinct
nonempty \(S_A,S_B\) have
\(|S_A\mathbin\triangle S_B|\ge2h\), whence
\[
 D:=\sum_{A<B}|S_A\mathbin\triangle S_B|\ge hk(d-1).
\]
If \(b_a\) directions are active in row \(a\), (E.29a) gives row weight at
least \(\max(1,d-b_a)\), while
\(b_a(d-b_a)\le k\max(1,d-b_a)\).  Summing gives
\(D\le k\,{\rm wt}(W)\) and hence
\({\rm wt}(W)\ge h(d-1)=ph\).  If some \(S_A=\Delta\), every row is nonzero
and the stronger bound \({\rm wt}(W)\ge|\Delta|>ph\) follows.  The
fixed-transverse rectangles attain \(ph\), hence
\[
                         \boxed{d_{\rm row}(D)=ph}.          \tag{E.29b}
\]

The same counting gives a sharp gap.  If \(k\ge2\) supports are active and
none is full, put \(R=\bigcup_AS_A\),
\(b_a=|\{A:a\in S_A\}|\), and
\(G=\sum_{a\in R}(k-b_a)\).  Pairwise symmetric differences give
\[
 \sum_{A<B}|S_A\mathbin\triangle S_B|
 =\sum_{a\in R}b_a(k-b_a)\ge hk(k-1),\qquad G\ge hk.
\]
For \(k<d\), (E.29a) yields
\[
 {\rm wt}(W)\ge\sum_{a\in R}(d-b_a)
 =(d-k)|R|+G\ge dh=|\Delta|.
\]
For \(k=d\), let \(r=|\{a:b_a=d\}|\).  Then (E.29a) gives
\[
 {\rm wt}(W)\ge\sum_{b_a<d}(d-b_a)+r,
\]
whereas the pair count is at most
\((d-1)\sum_{b_a<d}(d-b_a)\); hence again
\({\rm wt}(W)\ge dh+r\).  A full support makes every row nonzero and gives
the same bound directly.  If \(k=1\), a word below \(|\Delta|\) forces
\(S_A\) to be the radial \(h\)-cell and its nonzero \(B_A\)-coefficient to
be one affine \(p\)-block.  Thus the minimum words are exactly the
fixed-transverse rectangles and there is no weight strictly between \(ph\)
and \(|\Delta|\).

Consequently every puncture below \(ph\) is harmless.  If only
\(q\le h\) hard centres are nonzero, only those \(q\) localized halves are
needed and \(|U|\le q(p-1)<ph\).  In the all-active case, a support capable
of passing (E.27) has \(|U|\le|H|<|\Delta|\), while the Möbius midpoint
bound forbids containment of a minimum rectangle.  Therefore the actual
structured \(D_U\) is onto throughout the balanced zero-odd branch-C
regime, for zero-centre and all-active cases alike.  This
closes the divided **mod-two** gate, not (E.20): its prescribed Hamming
weight, every exact direction slice (E.21), and the integral zero-one
equations remain open.  See
`evidence/NOTE_2026-09-03_GROUPED_UNCERTAINTY_SQUARE.md` and
`evidence/NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE_GAP.md`.

There is also a sharp warning against treating onto plus the scalar quota
bounds as a rounding theorem.  With
\(N=|\Delta|=2h(h+1)\) and
\(\operatorname {rank}D=R=2h(h+1)^2\), any branch-C quota slice of total
weight \(s\le(N-1)/2\) has fewer than \(2^{h(N-1)}<2^{R-d}\) points.
An onto map has exactly \(2^{R-d}\) syndromes over one fixed parallel-parity
vector, so some compatible syndromes with the same feasible quotas have no
binary preimage.  This does not obstruct the actual
\(\widehat T_U\); it proves that the actual transverse target, beyond
surjectivity and (E.21), is indispensable.  See
`evidence/NOTE_2026-09-03_SYMMETRIC_QUOTA_CARDINALITY_BARRIER.md`.

The actual balanced branch-C direction parities do sharpen the support
floor. Put \(m=(p+1)/2\), \(s=(t+1)\bmod(p+1)\), and
\(\kappa=t_{\max}-t+1+j\). A localized Mobius half in hard direction
\(L_i\), with auxiliary \(M_i\), has parallel word
\(e_{L_i}+e_{M_i}\pmod2\); arbitrary ternary cancellations change every
direction count by an even integer. Comparing this word with the balanced
hard and opposite quotas in (E.21) forces

\[
 j\ge
 \begin{cases}
  2,&5\le s\le m,\\
  1,&s\in\{4,m+1\},\\
  0,&\text{otherwise}.
 \end{cases}                                               \tag{E.29c}
\]

Thus \(j=0\) is impossible for \(4\le s\le m+1\), and \(j=1\) is wholly
impossible for \(5\le s\le m\). At \(s\in\{4,m+1\}\), \(j=1\) forces
three fixed edges and no divided column. This parallel-parity theorem
survives triple and higher overlaps, but only excludes these endpoint
slices; it does not solve the remaining target equations. See
`evidence/NOTE_2026-09-03_MOBIUS_PARALLEL_PARITY_ENDPOINT.md`.

Neither equality shortcut suggested by this proof is valid. In the
complementary-profile construction the relative auxiliary scales are fixed,
so a singleton fixed edge on \(\ker F\) would require

\[
                         M_i(x)^2=4j_i^2\qquad\text{for every }i. \tag{E.29d}
\]

But the hard-star contribution is independent of \(j_i\) in every moment
degree retained by Proposition 15.759 (the sole top-degree value is also
independent of nonzero \(j_i\)), and the mass, parallel, centrality, integral,
and binary compatibility data likewise permit the centres to vary
independently after a scaled auxiliary family is fixed. Hence the
centre-coherence condition (E.29d) is not automatic for a preassigned family;
choosing that family adaptively from the centers remains open. Nor
does the local geometry contradict it: the exact four-candidate equations
have a clean one-overlap point \(q=r=2\), \(A=B=3/4\), with the other three
candidates absent. Finally the actual compact source
\(K(v,-v;0)\) has fixed word \(e_{[v]}\), exactly one silent affine-block
group, and automatically satisfies the full common-moment system. Thus
grouped uncertainty is sharp and cannot supply a second silent group. These
are method barriers, not a Boolean completion or exclusion; the remaining
argument must couple the branch atom counts, signs, and quotas to the actual
Möbius support. See
`evidence/NOTE_2026-09-03_MOBIUS_ENDPOINT_BARRIER.md`.

The adaptive part of the centre problem can in fact be solved uniformly.
Normalize a prospective singleton line by \(0\ne x_0\in\ker F\), put

\[
 X_i=L_i/j_i,\qquad \alpha_i=X_i(x_0),\qquad x=cx_0.
\]

For singleton signs \(\epsilon_i,\epsilon_k\), the exact complementary-pair
equations give

\[
 \rho={2\epsilon_i\over c\alpha_k-2\epsilon_k},\qquad
 \nu={4\epsilon_i\epsilon_k\over
 (c\alpha_i-2\epsilon_i)(c\alpha_k-2\epsilon_k)}.       \tag{E.29e}
\]

Colour a valid signed endpoint by
\(\eta(\epsilon(c\alpha-2\epsilon))\). Exactly \(m\) of the \(p-1\)
nonzero evaluations are monochrome. Averaging over \(c\) therefore gives a
scale with at most \(m/2\) monochrome targets; the graph obtained by deleting
the two monochrome cliques has an explicit perfect matching. Choosing
opposite colours on its edges makes every \(\nu\) in (E.29e) nonsquare. Thus
arbitrary nonzero hard centres admit a centre-coherent **target matching**.

This is not the endpoint construction. In the affine chart of directions
other than \(F\), a target pair \((i,k)\) forces the auxiliary coordinates

\[
 U={w_kz_i+(w_i-c)z_k\over w_i+w_k-c},\qquad
 V={(w_k-c)z_i+w_iz_k\over w_i+w_k-c},\qquad
 w_i={2\epsilon_i\over\alpha_i}.                         \tag{E.29f}
\]

All outputs must be distinct and equal the prescribed set containing
\(m-2\) hard and two opposite directions. Equivalently, for a bijection
\(\sigma:H\to A\), an involution \(\tau\), and the cross-assignment
\(\phi(\sigma(\tau i))=i\), put
\(g_V=\alpha_{\phi(V)}(V-z_{\phi(V)})\). Every auxiliary pair must obey

\[
 g_U^2=g_V^2,\qquad {(U-V)^2\over g_U^2}={c^2\over4}.     \tag{E.29g}
\]

No prescribed-set paired SDR satisfying (E.29g) is proved. See
`evidence/NOTE_2026-09-03_ADAPTIVE_MOBIUS_PAIRING.md`.

There is also now an exact target-sensitive fixed-word criterion. Let
\(c_U(D,\beta)\) be the parity of surviving nonzero-\(\Phi\) orbits over the
affine block \((D,\beta)\), let \(\ell\) be the hard literal word, let
\(s_x\) be the singleton word, and let \(z\) record which target triangles
contain an antipodal label pair. Orthogonality of the affine-block incidence
matrix gives

\[
 M^{\mathsf T}a_Y=\ell+z,\qquad
 a_Y+\Phi(U)=e_{[x]}\Longleftrightarrow z=c_U+\ell+s_x.  \tag{E.29h}
\]

For the central fixed-word/odd-moment layer, existentially over atom labels
with the prescribed compact/all-equal counts, the row condition is exactly
\[
 |c_U(D,\cdot)+\ell_D+s_{x,D}|\le n_D,qquad
 |c_U(D,\cdot)+\ell_D+s_{x,D}|\equiv n_D\pmod2.          \tag{E.29i}
\]
It forces \(\Lambda\ge\kappa_0+m+q\); separating zero-\(\Phi\)
cancellations gives the higher-overlap-safe
\(\sigma\ge\kappa_z+m+q\). One half has distinct nonzero block types, and
two distinct halves share at most eight; an explicit disjoint ternary pair
over \(\mathbf F_{31}\) shares three, disproving the proposed bound one.
The valid scalar bounds retain positive room. Equations (E.29h)--(E.29i)
do not solve even moments, nonfixed target cells, the endpoint, or residual
(ii). See
`evidence/NOTE_2026-09-03_MOBIUS_FIXED_WORD_ATOM_COUPLING.md`.

Second, fixed-word parity survives every ternary cancellation at \(p=31\).
Each of the sixteen halves has one zero-\(\Phi\) occurrence and 29
nonzero-\(\Phi\) occurrences.  If \(\kappa_0,\kappa_1\) count cancellation
units of the two types, then

\[
 u_0=16-2\kappa_0,\qquad
 u_{\rm np}=464-2\kappa_1,\qquad
 |U|=480-2(\kappa_0+\kappa_1).                              \tag{E.30}
\]

Here zero-\(\Phi\) parallel means midpoint parallel to difference; it is
not the directionwise parallel-edge count in (E.21).  Summing (E.19) over
\([v]\) cancels every \(P_L\) term and counts each nonzero fixed cell \(p\)
times.  For the full central graph target, the opposite nonfixed cells cancel
in pairs, every row total is \(|H|\), and the parallel totals sum to
\(|H|\).  Consequently

\[
 |a(T_U)|\equiv |H|+u_{\rm np}\equiv1,\qquad
 |H|-|U|-|a(T_U)|\equiv u_0\equiv0\pmod2.                  \tag{E.31}
\]

Thus the Hamming numerator is automatically even for the sixteen-half
construction; parity does not exclude it.  At the minimum
\(\kappa=178-t\), one has \(|H|-|U|=1\), so any completion must have
\(|a(T_U)|=1\) and no selected unused double orbit.  This is conditional
rigidity, not existence.

Third, the source word of the sharp rigid two-cancellation pair is explicit:

\[
 |\Phi_{\rm pair}|=
 4p-26-2\eta(3)-4\eta(6)-4\eta(-2),                         \tag{E.32}
\]

which is 108 at \(p=31\).  This does not determine the full target coset.
The previously used identity \(a_Y=a_{\rm literal}\) is false and is
retracted.  A centrally symmetric compact residual can have odd fixed-cell
coefficients, and the correct decomposition is

\[
 a_Y=a_{\rm literal}+a_{\rm compact},\qquad
 a(T_U)=a_{\rm literal}+a_{\rm compact}+\Phi(U).             \tag{E.33}
\]

Equation (E.31) uses this full target and fixes only its total parity; it
does not control the unknown compact support or turn (E.32) into a Hamming
obstruction.

Fourth, on an equal-square proposed common block \(X=C=B_{K,r^2}\), the
prescribed-center incidence problem has an exact Hall formulation.  For
half \(i\), its doubled pair is precisely

\[
 x_i\in A_i=\{a:L_i(a)=\pm j_i/2\},\qquad
 y_i\in X\setminus T_i,\quad
 T_i=A_i\cup\{a:L_i(a)=0\}.                                \tag{E.34}
\]

A subfamily \(P\) with
\(\left|\bigcup_{i\in P}A_i\right|\le |P|-2\) excludes both the fully
doubled and one-single saturated profiles.  The recorded \(p=31\)
prescribed-center list has such a six-edge/four-anchor deficiency for one
fixed common direction \(K\), and no other \(K\) is excluded by that
witness.

There is also a constructive sufficient condition for branch primes
(pge31), under the explicit hypothesis that every hard center is nonzero.
Choose one prescribed hard pair as the common direction and form, on its affine anchor line
\(\ell\), the simple graph whose \(h\) edges are the other anchor pairs
\(A_i\).  If every component of this link graph is a tree or unicyclic,
then it is a pseudoforest, its edges have an incident-vertex SDR, the
remaining free slots satisfy Hall, and the dependent half supplies the last
point.  Therefore

\[
 G_\ell\ {\rm pseudoforest}
 \quad\Longrightarrow\quad
 \hbox{a saturated equal-square common-block incidence cover}. \tag{E.35}
\]

This is a proved one-way implication.  It is not known that one prescribed
anchor line always satisfies it, and (E.35) proves neither mutual ternarity,
containment of a full row-code support, nor the divided Boolean completion.

Evidence and focused replay for these post-15.761 statements are in
`evidence/NOTE_2026-09-02_COMPACT_RAY_HIGHER_MOMENT_GATE.md`,
`evidence/NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md`,
`evidence/NOTE_2026-09-03_P31_EQUIANHARMONIC_ZERO68_MITM.md`,
`evidence/NOTE_2026-09-03_EDGE_RADON_SIGNED_BOOLEAN_DEFECT.md`,
`evidence/NOTE_2026-09-03_EDGE_RADON_RIDGE_KERNEL.md`,
`evidence/NOTE_2026-09-03_EQUIANHARMONIC_COMPONENT_PACKING.md`,
`evidence/NOTE_2026-09-03_HARD_ROW_COMPACT_ODD_RADON_CENTRALITY.md`,
`evidence/NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md`,
`evidence/NOTE_2026-09-03_HARD_STAR_ANTISYMMETRIC_SUPPORT.md`, and
`evidence/NOTE_2026-09-03_EQUIANHARMONIC_THRESHOLD_EVEN_BARRIER.md`, together
with `evidence/NOTE_2026-09-03_INVERSION_SYMMETRIC_LATTICE.md`,
`evidence/NOTE_2026-09-03_MOBIUS_HALF_SYMMETRIC.md`, and
`evidence/NOTE_2026-09-03_ALL_ACTIVE_PENCIL_SUPPORT.md`, together with
`evidence/NOTE_2026-09-03_SYMMETRIC_FIXED_EDGE_ELIMINATION.md`,
`evidence/NOTE_2026-09-03_SYMMETRIC_HALVED_MOD2.md`,
`evidence/NOTE_2026-09-03_SYMMETRIC_HALVED_MOBIUS_COVER.md`,
`evidence/NOTE_2026-09-03_SYMMETRIC_UNUSED_SLICE_EXCHANGE.md`, and
`evidence/NOTE_2026-09-03_MOBIUS_HALF_INTERSECTIONS.md`, together with
`evidence/NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE.md`,
`evidence/NOTE_2026-09-03_PRESCRIBED_CENTER_COMMON_BLOCK.md`,
`evidence/NOTE_2026-09-03_RIGID_PAIR_FIXED_WORD.md`,
`evidence/NOTE_2026-09-03_GROUPED_UNCERTAINTY_SQUARE.md`,
`evidence/NOTE_2026-09-03_SYMMETRIC_HALVED_ROW_CODE_GAP.md`, and
`evidence/NOTE_2026-09-03_SYMMETRIC_QUOTA_CARDINALITY_BARRIER.md`, together
with `evidence/NOTE_2026-09-03_MOBIUS_PARALLEL_PARITY_ENDPOINT.md`. Replay with

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q \
  tests/test_compact_ray_moment_gate.py tests/test_conic_odd_radon.py \
  tests/test_p31_equi_zero68_mitm.py tests/test_signed_boolean_defect.py \
  tests/test_ridge_kernel.py tests/test_equianharmonic_component_packing.py \
  tests/test_hard_compact_odd_radon.py \
  tests/test_hard_star_antisymmetric_support.py \
  tests/test_inversion_antisymmetric_radon.py \
  tests/test_equianharmonic_threshold_even_barrier.py \
  tests/test_inversion_symmetric_lattice.py \
  tests/test_mobius_half_symmetric.py tests/test_all_active_pencil_support.py \
  tests/test_symmetric_fixed_edge_elimination.py \
  tests/test_symmetric_halved_mod2.py \
  tests/test_symmetric_halved_mobius_cover.py \
  tests/test_symmetric_slice_exchange.py \
  tests/test_mobius_half_intersections.py \
  tests/test_symmetric_halved_row_code.py \
  tests/test_prescribed_center_common_block.py \
  tests/test_rigid_pair_fixed_word.py \
  tests/test_grouped_uncertainty_square.py \
  tests/test_symmetric_halved_row_code_gap.py \
  tests/test_symmetric_quota_cardinality_barrier.py \
  tests/test_mobius_parallel_parity_endpoint.py \
  tests/test_mobius_endpoint_barrier.py \
  tests/test_main_chain_docs.py
```

Residual (ii), E(1), $L=1/2$, and the original MathOverflow problem all
remain **OPEN**.  The finite conic-fiber `UNSAT` result and the exact
Boolean/kernel reductions must not be promoted to any of those closures.
Equivalently: Residual (ii), E1, `L=1/2`, and the original MO limit remain
OPEN.

## Proposition 15.723 — paired-cube obstruction to middle floor-plus-two lifts

Put (p=2m-1\ge17). Let (A\ge0) be an integer-valued quadratic on
(J(p,m)) with

\[
 A(X)\equiv |X\cap B|+\eta\pmod2,qquad b=|B|\ \hbox{odd},
\]

and suppose

\[
 2p\,\mathbb EA=2p+2.                               \tag{15.723.1}
\]

First consider a nonnegative integral quadratic (g) on a Boolean cube
with parity (1+\sum_{i\in R}z_i). Möbius inversion makes its multilinear
coefficients integral and its mean half-integral. If (r=|R|\ge3), Fourier
orthogonality gives

\[
 \mathbb E[g\chi_R]=0,
\]

while (g\ge1) on the (chi_R=1) half, so (\mathbb Eg\ge1). If
(r\ge5) and equality held, (h=g-1) would vanish on that half and satisfy
(h\chi_R=-h). For every (|S|\le2),

\[
 \widehat h(S\mathbin\triangle R)=-\widehat h(S).
\]

The left side vanishes because (|S\mathbin\triangle R|\ge3), so every
Fourier coefficient of (h) vanishes. This would make (g=1) on the
opposite parity half, impossible. Half-integrality therefore sharpens the
strict inequality to

\[
 \boxed{\mathbb Eg\ge{3\over2}\quad(r\ge5).}        \tag{15.723.2}
\]

The bound is sharp for (r=5,6), using
((\sum_{i\in R}z_i-3)^2).

The paired-cube operator through (X\in J(p,m)) leaves one point of (X)
unpaired, bijects the other (m-1) points with (X^c), and averages the
resulting cubes. On degree at most two,

\[
 TA(X)={A(X)+p\mathbb EA\over p+1}.
\]

Under (15.723.1),

\[
 TA(X)=1+{A(X)\over p+1}.                           \tag{15.723.3}
\]

Replace (B) by its complement on the slice when necessary and write

\[
 (k,e)=
 \begin{cases}(b,\eta),&b\le m-1,\\
 (p-b,\eta+m\bmod2),&b\ge m.
 \end{cases}
\]

For (5\le k\le m-1), apply (15.723.2)--(15.723.3) at a contact layer and
then match the first two hypergeometric moments by three positive contact
nodes.  In phase (e=1) the nodes are exactly
`0, 2 floor(k/4), 2 floor(k/4)+2`; the zero node has at least five active
parity coordinates and is upgraded from its parity baseline to at least
`m`.  The excess over (1+1/p) is
((m-1)w_0-1/p), with

\[
 w_0=\begin{cases}
 {p-k-3\over p(k+4)},&k\equiv0\pmod4,\\
 {-k^2+kp-k-3p\over p(k-1)(k+3)},&k\equiv1\pmod4,\\
 {-k^2+kp+k-4p\over p(k-2)(k+2)},&k\equiv2\pmod4,\\
 {p-k\over p(k+1)},&k\equiv3\pmod4.
 \end{cases}                                       \tag{15.723.4}
\]

To make the all-parameter sign step explicit, put `d=p-(2k+1)`, an even
nonnegative integer.  After multiplication by the positive denominator,
the four phase-one numerators for `k mod 4 = 0,1,2,3` factor as

```
(d+k-4)(d+2k+2)
(d+2k+2)(d(k-3)+k^2-8k+3)
(d+2k+2)(d(k-4)+k^2-8k+4)
(d+k-1)(d+2k+2).
```

The only zeros in the allowed range are `(k,d)=(5,6)` and `(6,4)`, both
at `p=17`.  In phase (e=0), (k\ge7), put `s=floor(k/4)`, take the first
node `2s+1` when `k=3 mod 4` and `2s-1` otherwise, the second node two
larger, and the far node `k` for odd `k` or `k-1` for even `k`.  At the
far node every paired cube has at least `k-1` active coordinates for odd
`k` and `k-3` for even `k`, hence at least five.  The four positive gap
numerators become

```
(k-4)d^2+(3k^2-14k-4)d+2k^2(k-7)
(k-3)d^2+(3k^2-8k-3)d+2k^3-6k^2-10k+6
k d^2+(3k^2+6k)d+2k^3+10k^2+8
(k-3)d^2+(3k^2-12k-3)d+2k^3-14k^2-2k+6.
```

Their residue-class minima are positive; the only negative-looking
constant is the last polynomial at `k=7`, where `p>=17` gives `d>=2` and
the value is 128.  For (k=5,6), nodes `1,3,5`, combined endpoint weight
`3(p-5)/(8p)`, and the exact good-cube fraction `(m-5)/m` give numerator
`3(p-5)(p-9)-16` for `p=1 mod 4` and
`3(p-5)(p-11)-16` for `p=3 mod 4`, positive from `p=17` and `p=19`
respectively.  Thus no finite scan supplies the universal quantifier.

Because the original (b) is odd, the two zero-gap cells translate exactly
to

\[
 \boxed{(p,b,\eta)=(17,5,1),\ (17,11,0).}           \tag{15.723.5}
\]

They are genuine. On the smaller parity side (C), set
(t=|X\cap C|) and (A=(t-3)^2). For (k=5,6) at (p=17),

\[
 \mathbb EA={18\over17},\qquad2p\mathbb EA=36=2p+2.
\]

Thus every middle floor-plus-two cell is excluded except (15.723.5), and
those two exceptions must remain in subsequent profile ledgers.

Evidence: `src/e1_gmin_m4_prop15723.py`,
`tests/test_prop15723.py`, and
`evidence/e1_gmin_m4_prop15723.json`.

## Proposition 15.724 — full Miquelian-circle boundary exclusion

Assume (D=\partial H) is a full Miquelian circle, (|D|=P=p+1=2m), and
(|H|=4p+1=8m-3). Normalize a boundary point as in Proposition 15.722.
Every boundary vertex is nonisolated. Since (H) has at most (8p+2)
nonisolated vertices, at most (7p+1) of the (p^2-p) vertices outside
(D) are nonisolated. For (p\ge17), choose an isolated outside point
(w).

Send (w) to infinity. Equations (15.722.6)--(15.722.7) make the new
product sign equal to the circle's (b=2) direction type. Hence the chart
has

\[
 I=0,\qquad m\hbox{ phase-zero }b=0,qquad
 m\hbox{ phase-one }b=2.                            \tag{15.724.1}
\]

For either type, write the exact directional means as

\[
 a_d=I+P P_d-\epsilon_dT-3p=2u+Pk_d,qquad
 \sum_dk_d=m-u.                                    \tag{15.724.2}
\]

The phase-one (b=2) floor is (P-2). For (1\le u\le m-2), every
direction would require (k_d\ge1) although their sum is below (m).
For (u=0), every direction would be a forbidden floor-plus-two lift.
Therefore (u=m-1), with parallel counts

\[
 x,\ldots,x,x+1.                                   \tag{15.724.3}
\]

For the phase-zero type write (P_d=y+k_d). Counting finite edges in
(15.724.2), using (I=0), gives

\[
 8m-3=m(x+y+1)+1-u,qquad
 m(x+y-7)=u-4.                                     \tag{15.724.4}
\]

Since (0\le u<m), this forces (u=4) and (x+y=7). A phase-one xnor
baseline has the exact coefficient congruence

\[
 q={p-1\over2}\mid I+x-4.                          \tag{15.724.5}
\]

This is the sign-independent two-coordinate congruence from Proposition
15.673. That proposition treats both \(4+z_a z_b\) and \(4-z_a z_b\)
(XNOR and XOR in zero-one variables); its symbolic sign \(\tau\) drops out
of the coefficient comparison \((p-1)c=I+P_d-4\). Thus the present XNOR
chart does not rely on identifying the two baselines.

Now (q\ge8), (I=0), and (x,y\ge0), so

\[
 \boxed{(u,x,y)=(4,4,3).}                          \tag{15.724.6}
\]

The phase-zero quotient sum is (m-4). At least four of its (m)
directions therefore have (k_d=0) and (a_d=8). For such a phase-zero
(b=0) direction, parity gives (A_d=2B_d), where (B_d) is a nonzero
nonnegative integer-valued quadratic. Thus

\[
 4p\mathbb EB_d=8.
\]

Proposition 15.688 instead gives

\[
 4p\mathbb EB_d\ge p-3\ge14,                       \tag{15.724.7}
\]

a contradiction. Therefore

\[
 \boxed{\text{no full Miquelian-circle boundary is residual-compatible
 for }p\ge17.}
\]

Together with Proposition 15.722, this excludes the complete outside
pair-slack-zero branch. Proposition 15.722 also excludes every positive
outside slack through `max(3,floor(sqrt(p)-5/2))`. Slack beyond that cutoff,
the whole (p+1) shell, residual (ii), multi-level Type I, and the limit remain
open.

Evidence: `src/e1_gmin_m4_prop15724.py`,
`tests/test_prop15724.py`, and
`evidence/e1_gmin_m4_prop15724.json`.

## Proposition 15.725 — retracted parabola-plus-internal family close

Consider the explicit boundary

\[
 D=\{x+x^2\omega:x\in\mathbf F_p\}\cup\{a\omega\},
 \qquad \omega^2=\nu,\quad \chi(\nu)=\chi(-a)=-1.
\]

Sending (a\omega) to infinity gives exact coordinates

\[
 (A_x,B_x)=\left({x\over Q(x)},-{x^2-a\over Q(x)}\right),
 \qquad Q(x)=x^2-\nu(x^2-a)^2.
\]

For one product-sign orientation, direct modular enumeration at
(p=17,19,23,29,31,37,41,43,47) checks 2,381 parameter cases and 92,664
typed directions; every typed floor sum exceeds its budget. This is exact
finite evidence.

The proposed all-prime continuation is not a proof. It assumes bounds for
three quartic-fibre discriminant/resolvent character sums without deriving
their curves, genera, singular fibres, or points at infinity. In particular,
the admissible locus (4a\nu+1=0) is singular in the asserted generic
discriminant model. The opposite product-sign orientation is also unchecked.
Therefore

\[
 \boxed{\text{the parabola-plus-internal family remains OPEN}.}
\]

Proposition 15.725 has no downstream role.

## Proposition 15.726 — tangent-envelope linear low-slack exclusion

Continue in the normalized first shell of Proposition 15.722: (D) is a
set of (p+1) affine points, where (p\ge17) is prime, and its outside
pair slack is

\[
 R=\sum_\ell h(n_\ell),\qquad
 h(2r)=r(r-1),\quad h(2r+1)=r^2.                 \tag{15.726.1}
\]

Suppose

\[
 1\le R\le\left\lfloor{p-4\over3}\right\rfloor. \tag{15.726.2}
\]

Delete points until the remainder (A=D\setminus T) is an arc, and choose
(T) inclusion-minimal with this property.  Put (t=|T|).  The deletion
construction from Proposition 15.722 and the elementary bounds
(h(n)\ge n-2) give

\[
 1\le t\le R.                                     \tag{15.726.3}
\]

Minimality says that every (z\in T) lies on an (A)-secant.  Let (s_z)
be the number of (A)-secants through (z), and put

\[
 I=\sum_{z\in T}s_z.
\]

For a line (ell), write
(a_\ell=|A\cap\ell|\le2) and (u_\ell=|T\cap\ell|).  Only lines with
(a_\ell=2) contribute to (I), and each such line contributes
(u_\ell).  On the other hand,

\[
 h(2+u)-u=
 \begin{cases}
  r(r-1),&u=2r,\\
  r^2,&u=2r+1,
 \end{cases}
 \quad\ge0.
\]

Consequently

\[
 I=\sum_{a_\ell=2}u_\ell
   \le\sum_\ell h(a_\ell+u_\ell)=R.             \tag{15.726.4}
\]

Now

\[
 |A|=p+1-t=p+2-\tau,\qquad \tau=t+1.
\]

Since (3t\le3R\le p-4), one has

\[
 |A|>2\tau+2.                                     \tag{15.726.5}
\]

The odd-order tangent-envelope theorem of Ball--Lavrauw therefore supplies
a nonzero homogeneous polynomial (Phi) in the dual plane, of degree
(2\tau), such that

\[
 \Phi(X\mathbin\times P)=f_P(X)^2\qquad(P\in A), \tag{15.726.6}
\]

where (f_P) is the product of the (A)-tangent forms at (P).  This is Theorem
11 in arXiv v4 and Theorem 13 in the authors' current manuscript.  Their
common size hypothesis is (|A|\ge2\tau+2), which (15.726.5) satisfies with
room to spare.

Fix (z\in T).  There are exactly

\[
 |A|-2s_z
\]

(A)-tangents through (z).  If this number exceeded (2\tau), their
distinct dual points would give more zeros than the degree of
(Phi|_{z^*}), forcing that restriction to vanish identically.  But
minimality supplies an (A)-secant (zPQ).  At its dual point,

\[
 \Phi(z\mathbin\times P)=f_P(z)^2\ne0,
\]

because (Pz=PQ) is a secant rather than a tangent at (P).  This is a
contradiction.  Hence

\[
 |A|-2s_z\le2\tau,
 \qquad
 s_z\ge {p-1-3t\over2}.                           \tag{15.726.7}
\]

Summing (15.726.7) gives

\[
 I\ge F(t):={t(p-1-3t)\over2}.                    \tag{15.726.8}
\]

The quadratic (F) is concave, so its minimum on (1\le t\le R) occurs
at an endpoint.  Using (3R\le p-4),

\[
 F(1)={p-4\over2}\ge{3R\over2}>R,
 \qquad
 F(R)={R(p-1-3R)\over2}\ge{3R\over2}>R.          \tag{15.726.9}
\]

Thus (I>R), contradicting (15.726.4).  Therefore

\[
 \boxed{1\le R\le\left\lfloor{p-4\over3}\right\rfloor
 \text{ is impossible for every prime }p\ge17.}  \tag{15.726.10}
\]

Equivalently, after the already-closed (R=0) branch, every surviving
positive slack must satisfy

\[
 \boxed{R\ge\left\lfloor{p-1\over3}\right\rfloor.} \tag{15.726.11}
\]

This strictly supersedes Proposition 15.722's square-root cutoff in the
current frontier.  It does not close the rest of the (p+1) shell, larger
boundary shells, residual (ii), multi-level Type I, or the limit.

Evidence: `src/e1_gmin_m4_prop15726.py`,
`tests/test_prop15726.py`,
`evidence/e1_gmin_m4_prop15726.json`, and
`evidence/NOTE_2026-08-30_tangent_envelope_linear_low_slack.md`.
The finite-geometry input is S. Ball and M. Lavrauw, *Planar arcs*,
J. Combin. Theory Ser. A **160** (2018), 261--287, Theorem 11 in arXiv v4,
doi:10.1016/j.jcta.2018.06.015.
