# The complementary Möbius endpoint and its center-coherence barrier

Date: 2026-09-03

Status: proof-grade necessary compatibility and method warning.  A
complementary-profile construction at the top balanced all-active endpoint
requires an additional coherence of the hard-star centers.  That coherence
is invisible to every aggregate, compact-centrality, and common-moment
constraint proved so far, so it is not automatic for a preassigned scaled
auxiliary family.  Choosing the auxiliaries adaptively from arbitrary centers
may still satisfy it; that global construction and the final restricted
Boolean fibre remain open.  No residual case is closed.

## 1. Where center coherence enters

Let

\[
 p=4r+3,\qquad h={p-1\over2},\qquad m={p+1\over2}=h+1,
\]

and consider the top balanced all-active endpoint with cancellation offset
zero.  The exact Hamming ledger leaves one fixed antipodal source edge and no
unused doubled orbit.  In the complementary-profile attempt, the ordinary
parallel vector of a half with target \(L_i\) and auxiliary \(M_i\) is the
vector in (E.12) of `solution.md`.  Pairing its projective square block with
its complement fixes the *relative functional scale* of \(M_i\); independent
rescaling of the auxiliaries no longer preserves the complementary parallel
identity.

Let the remaining fixed edge be represented by

\[
                         \{x_0,-x_0\},\qquad 0\ne x_0\in\ker F,
\]

where \(F\) is its spatial difference direction.  For this complementary
endpoint construction, matching the fixed-cell equations to that singleton
forces, after one common normalization,

\[
                         M_i(x_0)^2=4j_i^2
                         \quad\hbox{for every }i.          \tag{1}
\]

Equivalently, before scaling the representative \(x_0\), it forces

\[
             {j_i^2\over M_i(x_0)^2}
             ={j_k^2\over M_k(x_0)^2}
             \quad\hbox{for all }i,k.                    \tag{2}
\]

Here every denominator is nonzero because the proposed exceptional direction
\(F\) is not one of the auxiliary directions.  Since \(\ker F\) is
one-dimensional, replacing \(x_0\) by \(c x_0\) multiplies every ratio in
(2) by the same square.  Thus (2) is an intrinsic condition, not a choice of
representative.

Equations (1)--(2) are necessary only for this complementary endpoint
construction.  Neither fixed-edge elimination (E.18)--(E.22) nor the general
Möbius lift theorem says that every feasible endpoint must use complementary
profiles or satisfy (2).

## 2. The proved target constraints do not see the centers

A branch-C hard row has the exact form

\[
                        C_L=-S_{j_L}+K_L^{\rm compact}.   \tag{3}
\]

The center \(j_L\) does not occur in the balanced hard mass or parallel
quota: (15.758.9) gives

\[
 4p\,\mathbf E(\text{hard row})=p-1+e_L(p+1),
 \qquad P_L=3+e_L.                                      \tag{4}
\]

It is also invisible to the full characteristic-\(p\) moment hierarchy.
For a fixed star center \(j\), every moment summand contains the diagonal-
vanishing factor \((j-u)^2\); after adjoining the omitted diagonal term, the
sum over \(u\in\mathbf F_p\) is a field power sum of degree at most \(p-2\).
It therefore vanishes in every degree \(2\le d\le p-2\), not merely in the
odd degrees.  This is recorded in
`src/e1_gmin_m4_hard_compact_odd_radon.py`, function
`unit_star_odd_blind_certificate`, and in
`NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md`, Section 4.
At degree \(p-1\), homogeneity and \(j\ne0\) make the star contraction
independent of \(j\): substituting \(u=jx\) factors out
\(j^{p-1}=1\).  Thus the retained top-degree equations also cannot compare
the nonzero centers.

The remaining proved target conditions also introduce no relation among the
centers:

* the ordinary integral compatibility lattice (15.760.1) uses only the
  common row total and the sum of the parallel coordinates;
* the compact-residual centrality theorem is applied after subtracting the
  star and hence is unchanged when \(j_L\) changes; and
* Proposition 15.757 constructs the binary boundary completion for an
  arbitrary list of hard centers by prescribing

  \[
       \widehat f(\lambda n_L)=\omega^{\lambda j_L}
       \quad(L\text{ hard},\ \lambda\ne0),
  \]

  with zero data on the opposite dual lines.

This shows why coherence is not automatic for a *preassigned* complementary
family. Fix its scaled auxiliaries \(M_i\) and a nonzero
\(x_F\in\ker F\). The
center list

\[
                         j_i=\sigma_i cM_i(x_F),
                         \qquad \sigma_i\in\{\pm1\},     \tag{5}
\]

satisfies (2).  Choose any \(\lambda\in\mathbf F_p^*\) with
\(\lambda^2\ne1\), and replace only \(j_k\) by \(\lambda j_k\).  The
literal star in (3), and hence the target row itself, changes.  Nevertheless,
the numerical data in (4), all compact atoms, all common moments, compact
centrality, and the existence of ordinary and binary compatible lifts remain
unchanged.  One ratio in
(2) is multiplied by \(\lambda^2\), while all the others are not.  No other
choice of \(x_0\in\ker F\) can repair this, because it rescales every ratio
simultaneously.

Therefore the current target theory permits both coherent and incoherent
center lists relative to one fixed scaled auxiliary family. It neither
forces nor forbids (2). This argument does not exclude an adaptive rule that
chooses a different complementary family after the centers are known.

## 3. Why this is not a Boolean-fibre obstruction

The conclusion above is deliberately limited to the proved target layer.
Modulo two, the literal star has a \(j_L\)-dependent fixed-cell word, and the
central compact residual can also have odd fixed-cell coefficients.  In the
notation of `NOTE_2026-09-03_RIGID_PAIR_FIXED_WORD.md`, the actual forced word
is

\[
       a(T_U)=a_{\rm literal}+a_{\rm compact}+\Phi(U).    \tag{6}
\]

No theorem presently controls \(a_{\rm compact}\).  Hence a common simple
graph, if one exists, may impose additional center coherence through the
restricted affine box even though the aggregate and moment equations do not.
Conversely, failure of (2) rules out only the complementary-profile attempt,
not every choice of antisymmetric lift or every Boolean completion.

## 4. The local overlap equations are consistent

The complementary construction cannot be rejected by the two-half
intersection equations alone.  In the normal form of
`NOTE_2026-09-03_MOBIUS_HALF_INTERSECTIONS.md`, put \(z=qr\).  The desired
opposite-swapped intersection has the exact locus

\[
 A={(q-1)(z-1)\over z},\qquad
 B={(r-1)(z-1)\over z},                                  \tag{7}
\]

with \(q,r,z\ne0\), \(q,r\ne1\), and \(z\ne1\).  Substitution into the
other three forced candidates shows that both same-orientation candidates
are absent.  The opposite-direct candidate occurs only at
\(q=r=1/2\).  Hence every admissible point of (7) except that one exceptional
point gives exactly one shared orbit, with opposite orientation.

For example,

\[
 q=r=2,qquad A=B={3\over4},qquad t=s=-{1\over2}          \tag{8}
\]

gives a unique clean cancellation in every characteristic at least 31.
With equal singleton signs, choosing \(X(x)=Y(x)=-6\) makes
\(F=X-Y\) annihilate \(x\) and gives \(M_1(x)=M_2(x)=2\), so the local
singleton equations are consistent as well.  This does not embed the pair
in a global complementary family: that step must also preserve its fixed
functional scales, all Paley character placements, (1), and every
cross-pair intersection.

Thus the local four-candidate system supplies neither a contradiction nor a
global construction.  The exact rational replay is in
`src/e1_gmin_m4_mobius_endpoint_barrier.py`.

## 5. Full common moments still do not force two silent groups

There is an exact source-side counterexample to that proposed shortcut.  For
any \(0\ne v\in\mathbf F_p^2\), take the centrally symmetric integral source
chain

\[
 z_v=\{v,-v\}-\{v,0\}-\{-v,0\}.                         \tag{9}
\]

In a row \(L\) with \(L(v)\ne0\), its image is the centered compact atom

\[
 K\bigl(L(v),-L(v);0\bigr)
 =\{L(v),-L(v)\}-\{L(v),0\}-\{-L(v),0\}.               \tag{10}
\]

Because (7) is an actual integral source chain, its full Radon target
automatically satisfies every ordinary compatibility and every common-moment
identity.  It is also central.  Its unique fixed source component modulo two
is the antipodal edge \(\{v,-v\}\); the other two edges form one nonfixed
inversion pair and disappear from the fixed block modulo two.  Fixed-edge
elimination therefore gives the singleton point word

\[
                              a=e_{[v]}.                 \tag{11}
\]

For each projective direction \(L\), the affine-block transform of (9) is
zero exactly when \(L(v)=0\), and otherwise has its single one in the block
\(B_{L,L(v)^2}\).  Thus it has exactly one silent direction and \(p\) active
direction groups:

\[
                 |a|=1,\qquad g=p,\qquad |a|+g=p+1.     \tag{12}
\]

This attains equality in grouped uncertainty.  Consequently centrality plus
the complete common-moment system cannot force two silent groups.  The
counterexample is a method barrier; it is not asserted to have the top
endpoint quotas or to lie in the punctured Möbius support.

## Canonical references

* `solution.md`, (15.758.9), (15.759.1)--(15.759.3),
  (15.760.1), and (E.12), (E.18)--(E.22).
* `evidence/NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md`, Section 4.
* `evidence/NOTE_2026-09-03_HARD_ROW_COMPACT_ODD_RADON_CENTRALITY.md`,
  Sections 1 and 5.
* `evidence/NOTE_2026-09-03_MOBIUS_PARALLEL_PARITY_ENDPOINT.md`, Section 5.
* `evidence/NOTE_2026-09-03_RIGID_PAIR_FIXED_WORD.md`, Section 3.
* `evidence/NOTE_2026-09-03_GROUPED_UNCERTAINTY_SQUARE.md`.
