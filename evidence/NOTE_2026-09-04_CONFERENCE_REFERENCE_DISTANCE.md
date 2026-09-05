# Conference-reference distance and the obstruction to regauging

Date: 2026-09-04.

Classification: **proved theorem and regauging obstruction**. This is not
a minimal-witness localization theorem, residual-(ii) closure, E1 closure,
or a limit theorem. The theorem uses no finite census, solver, or
computational certificate. No attainment of its odd-square distance bound
is claimed.

The scope was checked against [AGENTS.md](../AGENTS.md) and the
[proposition-deduplication audit](PROPOSITION_DEDUP_AUDIT_2026-08-30.md).
The protected
[even collision note](NOTE_2026-09-04_EVEN_MINIMAL_GAP4_TYPE_COLLISION.md)
and [odd collision note](NOTE_2026-09-04_ODD_MINIMAL_GAP4_SHIFTED_COLLISION.md)
remain unchanged. Their bounded-support arguments are not inputs here.
The [older local-stability no-go](NOTE_2026-08-29_global_minimality_and_local_stability_no_go.md)
remains valid: switching minimality and correct norm scale do not imply
subquadratic distance to Paley.

## 1. Disagreement support between two conference signings

Let \(C\ne C'\) be real symmetric conference matrices of order \(n>2\):
their diagonals vanish, their off-diagonal entries are signs, and
\[
C^2=(C')^2=(n-1)I.
\]
Orthogonality of different rows is a sum of \(n-2\) signs, so \(n\)
is even. All distances below count unordered edges, not matrix entries.

Let \(H\) be the disagreement edge set, \(h=|H|\), and
\[
B=\frac{C-C'}2=C\circ{\bf1}_H.
\]
Let \(V\) be its active vertex set, \(v=|V|\), and
\(d_j=\deg_H(j)>0\) for \(j\in V\). Expanding the equality of
conference squares gives
\[
\boxed{CB+BC=2B^2.}                                      \tag{1}
\]
Outside \(V\), the rows of \(B\), \(BC\), and \(B^2\) vanish.
Equation (1) therefore makes those rows of \(CB\) vanish. For each
column \(b_j=Be_j\),
\[
\operatorname{supp}(Cb_j)\subseteq V,\qquad
\|Cb_j\|_2^2=(n-1)\|b_j\|_2^2=(n-1)d_j.                 \tag{2}
\]

The case \(v=2\) is impossible. Its single disagreement edge has
endpoints \(i,j\), with \(b_j=C_{ij}e_i\). Then \(Cb_j\) is nonzero
at every coordinate except \(i\), including a coordinate outside
\(\{i,j\}\) because \(n>2\), contrary to (2). Hence \(v\ge3\).

Every coordinate of \(Cb_j\) has absolute value at most \(d_j\).
Consequently
\[
(n-1)d_j\le v d_j^2,\qquad d_j\ge\frac{n-1}{v}.
\]
Summing over \(V\) gives \(2h\ge n-1\). Since \(2h\) and \(n\)
are even,
\[
\boxed{d_H(C,C')=h\ge n/2.}                              \tag{3}
\]

## 2. Baseline equality forces a perfect matching

Suppose \(h=n/2\). Equality in the crude coordinate estimate above
must not be assumed; a sharper bound is needed.

Fix \(j\in V\), and write \(d=d_j\). At each of its \(d\)
neighbors \(i\in N_H(j)\), the term indexed by \(i\) in
\((Cb_j)_i\) vanishes because \(C_{ii}=0\). Those coordinates
therefore have absolute value at most \(d-1\), while the other
\(v-d\) coordinates in \(V\) have absolute value at most \(d\).
Thus
\[
\begin{split}
(n-1)d&\le(v-d)d^2+d(d-1)^2\\
      &=(v-2)d^2+d,
\end{split}
\qquad
d\ge\frac{n-2}{v-2}.                                    \tag{4}
\]
Division is valid because \(d>0\) and \(v\ge3\). Summing (4),
\[
n=2h\ge\frac{v(n-2)}{v-2}
=n+\frac{2(n-v)}{v-2}.
\]
Since \(v\le n\), this forces \(v=n\). The \(n\) positive integer
degrees sum to \(n\), so all are one:
\[
\boxed{h=n/2\ \Longrightarrow\
H\text{ is a perfect matching and }B^2=I.}               \tag{5}
\]

## 3. Matching equality is impossible at odd square order

Assume now \(n=p^2+1\), with odd \(p\ge3\). Primality is not needed
for this algebraic assertion.

Under the hypothetical equality \(h=n/2\), write \(C=B+D\).
From (1), \(B^2=I\), and \(C^2=(n-1)I\),
\[
BD+DB=0,\qquad D^2=(n-2)I.                              \tag{6}
\]
Indeed \(CB+BC=2I\); substituting \(C=B+D\) gives anticommutation,
and then expanding \(C^2\) gives the second identity.

Conjugate by a diagonal sign matrix so that every matching entry of
\(B\) is \(+1\), and order the vertices in matching pairs. Then \(B\)
has diagonal blocks
\[
J=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]
The diagonal pair blocks of \(D\) are zero. Every off-diagonal pair
block has sign entries, and (6) gives
\[
JM_{ij}+M_{ij}J=0,\qquad
M_{ij}=\begin{pmatrix}a&b\\-b&-a\end{pmatrix},
\qquad a,b\in\{+1,-1\}.                                 \tag{7}
\]

Use the orthonormal vectors
\[
e_i^+=\frac{e_{i,1}+e_{i,2}}{\sqrt2},
\qquad
e_i^-=\frac{e_{i,1}-e_{i,2}}{\sqrt2},
\]
with all plus vectors followed by all minus vectors. In this basis,
\[
B=\begin{pmatrix}I_m&0\\0&-I_m\end{pmatrix},\qquad
D=\begin{pmatrix}0&K\\K^{\mathsf T}&0\end{pmatrix},
\qquad m=n/2.
\]
Anticommutation gives the zero diagonal blocks; symmetry gives the
transpose. Formula (7) explicitly yields
\[
K_{ij}=(e_i^+)^{\mathsf T}D e_j^-=a-b,\qquad
K_{ji}=a+b,\qquad K_{ii}=0.                              \tag{8}
\]
Thus \(K\) has entries \(0,\pm2\). This entry calculation is needed:
the orthonormal change of basis alone would not justify integrality
after division by two.

Equation (6) implies \(KK^{\mathsf T}=(n-2)I_m\).
Set \(W=K/2\), an integral matrix, and \(q=(p-1)/2\). Then
\[
WW^{\mathsf T}
=\frac{p^2-1}{4}I_m
=q(q+1)I_m.                                            \tag{9}
\]
Since \(p\) is odd, \(m=(p^2+1)/2\) is odd. Taking determinants,
\[
(\det W)^2=[q(q+1)]^m.
\]
Every prime valuation on the left is even. Oddness of \(m\) forces
every prime valuation of \(q(q+1)\) to be even, making \(q(q+1)\)
an integer square. But \(q\ge1\) and
\[
q^2<q(q+1)<(q+1)^2.
\]
This contradiction excludes (5). Together with (3),
\[
\boxed{
C\ne C',\quad n=p^2+1,\quad p\ge3\text{ odd}
\ \Longrightarrow\
d_H(C,C')\ge n/2+1.}                                    \tag{10}
\]

No construction attaining \(n/2+1\) is supplied. This is a lower
bound, not an exact determination of the minimum conference distance.

For contrast only, the separate
[one-pair repair record](threshold_valley_probe.json) supplies conference
signings of order six with masks `2393` and `2641`, both squaring exactly
to `5I` and differing on the three-edge matching
`{(0,4),(1,5),(2,3)}`. This attains the general baseline `n/2` in (3)
at order six. That order is not an odd square plus one, so it neither
attains nor contradicts (10), and is not an input to its proof.

## 4. Unique nearest conference reference in the small-support region

Continue with \(n=p^2+1\), odd \(p\ge3\). Let
\(A=C\mathbin\triangle H\) be any complete signing, not necessarily
a conference matrix, and put \(h=|H|\). If
\[
\boxed{h<\frac{n+2}{4},}                                \tag{11}
\]
then for every different symmetric conference matrix \(C'\) of the
same order, the Hamming triangle inequality and (10) give
\[
d_H(A,C')\ge d_H(C,C')-d_H(A,C)
\ge n/2+1-h>h=d_H(A,C).                                 \tag{12}
\]
Hence \(C\) is the unique nearest conference matrix to \(A\).

This includes every conference class, every vertex relabeling, every
diagonal switching, and global negation: all preserve the conference
property. In particular, applying these operations to \(A\) cannot
reduce its disagreement count with fixed \(C\), because applying
the same operation to both matrices preserves Hamming distance.
Operations preserving \(C\) may still permute \(H\); unique nearest
reference does not mean a unique labeled support.

For the three first support sizes left by the relevant 15.774 bounds,
the exact four-times margins in (11) are
\[
\begin{array}{c|c}
h&n+2-4h\\ \hline
5p+6&(p+1)(p-21)\\
6p+6&p^2-24p-21\\
7p+8&(p+1)(p-29).
\end{array}
\]
All are positive for \(p\ge37\). In particular the largest size
\(7p+8\) already satisfies (11), since \((p+1)(p-29)>0\).
Therefore any actual signing with one of these support sizes has
\(C\) as its unique nearest conference reference. This does not assert
that the scalar survivors in
[15.774](NOTE_2026-09-04_SMALL_MASS_TWO_TYPE_BRIDGE.md)
have graph realizations.

## 5. Scope of the obstruction

If a hypothetical first four-gap witness lies in this region,
changing its conference reference cannot make it smaller or cheaply
repair its parity. Every genuinely different reference costs at
least \(n/2+1-h\), quadratic in \(p\) when \(h=O(p)\).

The theorem does not place a first four-gap witness in this region,
prove \(O(p)\) support, or send every minimal witness into an
existing E1 unit. It uses conference orthogonality, not the full
norm or all-deletions constraints of a witness. Those constraints
remain essential to any localization proof.

It also gives no \(o(p^3)\) bound on the improvement over Paley.
Controlling a closest global norm minimizer to Hamming distance
\(o(p^4)\) would suffice by the existing Max-Lipschitz bound,
but a first-four-gap witness is not automatically a global minimizer.
The all-size minimal-four-gap bridge, residual (ii), E1, and the
original limit remain open.
