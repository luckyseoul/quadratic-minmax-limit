# Directed half-cut identity and structure of the two-half orientations

**Status:** proved two all-order identities plus exact structural certificates
for the already stored orders `n=5,6,7,8`.  This note does not extend the
finite census and does not prove the multiplier-two estimate.

Let `A` be a symmetric signing, let `S` be a tournament sign matrix, and put

\[
 R_S=A\circ S,
 \qquad
 B(A,R)=\max_{x,y}\bigl(|Q_A(x)+Q_A(y)|+|x^TRy|\bigr).
\]

The tournament `S=A\circ R` is the gauge-invariant object: simultaneous
Seidel switching `(A,R) -> (DAD,DRD)` leaves `S` unchanged.  Thus an
"orientation by a vertex order" means that `S` itself, not merely a switched
copy of `S`, is a transitive tournament.

## 1. Exact directed half-cut norm identity

For a vertex set `U`, split its cut according to the tournament `S`:

\[
\begin{aligned}
 E_S^+(U)&=\{\{i,j\}:i\in U, j\notin U, S_{ij}=+1\},\\
 E_S^-(U)&=\delta(U)\setminus E_S^+(U).
\end{aligned}                                                    \tag{1}
\]

In the first line the names `i,j` are chosen so that the displayed ordered
pair points from `U` to its complement.  Let `A^E` denote `A` with the signs
on the unordered edge set `E` reversed.  Then, for **every** tournament `S`,

\[
 \boxed{
 {1\over2}B(A,A\circ S)
 =\max_{U\subseteq[n]}\Phi\left(A^{E_S^+(U)}\right).
 }                                                               \tag{2}
\]

To prove this, fix `x,y` and put `U={i:x_i=-y_i}`.  Write `I_y` for the
`A_ij y_i y_j` sum on edges not crossing `U`, and write `F_y,G_y` for the
corresponding sums on `E_S^+(U),E_S^-(U)`.  Directly,

\[
 Q_A(y)=I_y+F_y+G_y,\qquad Q_A(x)=I_y-F_y-G_y,
 \qquad x^T(A\circ S)y=2(G_y-F_y).                 \tag{3}
\]

Therefore

\[
\begin{aligned}
 {1\over2}\bigl(&|Q_A(x)+Q_A(y)|+|x^T(A\circ S)y|\bigr)\\
 &=|I_y|+|G_y-F_y|\\
 &=\max\{|Q_{A^{E_S^+(U)}}(y)|,
          |Q_{A^{E_S^-(U)}}(y)|\}.                \tag{4}
\end{aligned}
\]

Switching all vertices in `U` takes `A^{E_S^+(U)}` to
`A^{E_S^-(U)}`, so the two matrices have the same `Phi` norm.  Maximizing
(4) over `y` and then over `U` proves (2).  Since every skew signing `R`
corresponds uniquely to `S=A\circ R`, (2) also gives

\[
 {1\over2}\min_R B(A,R)
 =\min_{S\ {\rm tournament}}\max_U
   \Phi\left(A^{E_S^+(U)}\right).                  \tag{5}
\]

There is a compact matrix form of the same geometry.  Put
`D_U=diag(r_i)`, where `r_i=-1` on `U` and `r_i=+1` off `U`.  Then the
directed half-cut neighbor in (2) is exactly

\[
 A^{E_S^+(U)}
 ={1\over2}\left(A+D_UAD_U+D_UR-RD_U\right).       \tag{6}
\]

Equivalently, with `P=(I-D_U)/2`, `Q=(I+D_U)/2`,
`Z_U=Q-iP`, and `H=A+iR`,

\[
 A^{E_S^+(U)}
 =\operatorname{Re}(Z_U^*HZ_U).                   \tag{7}
\]

As `U` and a Boolean test vector vary, `Z_Uy` ranges over the entire
fourth-phase cube `{+1,-1,+i,-i}^n`.  Thus (7) is also a direct real-slice
proof of (2): the directed half-cut matrices are precisely the real parts of
all diagonal quarter-phase conjugates of the Hermitian completion `A+iR`.

Thus the unrestricted multiplier-two problem is exactly a simultaneous
neighbor problem: orient the edges so that reversing the outgoing half of
**every** cut keeps the cube norm at most `sqrt(2) Phi(A)` (plus the allowed
Dini error).  Global optimality of `A` gives the reverse inequality
`Phi(A^E)>=Phi(A)` and therefore does not settle (5).

Choosing `S` uniformly at random turns every fixed value in (2) into an
internal energy plus a Rademacher sum.  The resulting partition-function
criterion and its sharp universal entropy threshold are audited separately
in the [random-tournament threshold note](NOTE_2026-09-02_RANDOM_TOURNAMENT_PARTITION_THRESHOLD.md).

For a transitive `S=T_pi`, (2) specializes to the one-sided cut-flip
identity in the separate
[prefix/half-cut audit](NOTE_2026-09-02_ORDERED_SKEW_PREFIX_HALF_CUT.md),
where it also has an exact prefix-total-variation form.

## 2. An all-order edge-stability lemma

If two tournaments `S,T` differ on `f` unordered edges, then

\[
 \boxed{\ |B(A,R_S)-B(A,R_T)|\le 4f.\ }                 \tag{8}
\]

Indeed, changing the orientation of one edge `{i,j}` changes
`x^TRy` by

\[
 2R_{ij}(x_jy_i-x_iy_j),
\]

whose absolute value is at most four.  Summing over the changed edges,
using `||a|-|b||<=|a-b|`, and then taking the maximum over `(x,y)` proves
(8).

Let `tau(S)` be the minimum number of backward edges of `S` over all vertex
orders, and let

\[
 B_{\rm ord}(A)=\min_{\pi}B(A,A\circ T_\pi),
 \qquad B_*(A)=\min_SB(A,A\circ S).
\]

For any unrestricted minimizer `S_*`, (8) gives the exact comparison

\[
 0\le B_{\rm ord}(A)-B_*(A)\le4\tau(S_*).              \tag{9}
\]

Consequently an `O(n)` correction to an ordered orientation costs only
`O(n)`, which is harmless at the required Dini `n^(3/2)` scale.  The missing
theorem is not (8): one would have to control `tau(S_*)` for a minimizer (or
construct an equally good near-transitive `S`).  General tournament theory
does not supply such control; `tau(S)` can have quadratic order.

Together with (2), (9) says precisely that a low-feedback unrestricted
optimizer would transfer the problem to one-sided cut-flip stability, at
additive cost `2 tau(S_*)` on the half-cut side.  Neither identity supplies
the missing low-feedback theorem.

## 3. The `sqrt(2)` multiplier is forced

The target constant in (5) cannot be uniformly improved.  Define

\[
 T_n=\min_{\Phi(A)=m_n}\ \min_{S\ {\rm tournament}}\
       \max_U\Phi\left(A^{E_S^+(U)}\right).                    \tag{10}
\]

Then for every fixed `n_0>=2`,

\[
 \boxed{\quad
 \limsup_{j\to\infty}{T_{2^j n_0}\over m_{2^j n_0}}\ge\sqrt2.
 \quad}                                                        \tag{11}
\]

Indeed, if an eventual dyadic tail satisfied

\[
 T_n\le c m_n+r_n,\qquad c<\sqrt2,\qquad r_n=o(n^{3/2}),       \tag{12}
\]

then the equal-endpoint Hadamard lift and (2) would give

\[
 m_{2n}\le2T_n+n,
 \qquad
 \alpha_{2n}\le{c\over\sqrt2}\alpha_n+o(1).                 \tag{13}
\]

Since `c/sqrt(2)<1`, iteration along that tail forces
`alpha_(2^j n_0)->0`.  This contradicts the uniform Gaussian lower bound

\[
 \alpha_n\ge{\sqrt{1-1/n}\over\pi}.
\]

No Dini summability is needed for this obstruction; a vanishing normalized
error is enough.  Merely having (12) on unrelated infinitely many orders is
not enough, because the contraction needs a complete eventual dyadic tail.
Thus the apparent contact with the `sqrt(2)` line in the finite picture is
consistent with a genuinely critical all-order constant, not slack that can
be removed by aiming for an easier fixed multiplier.

## 4. What the four stored optimizers actually are

For the matrices in `scripts/original_mo_two_half_geometry.py`, the literal
tournaments `S=A circ R` have respectively

\[
 4,\ 6,\ 8,\ 12
\]

directed triangles at `n=5,6,7,8`.  Hence none is a vertex-ordering
orientation.  Explicit vertex orders put them within respectively
`2,3,3,4` edge reversals of a transitive tournament:

```text
n=5: (2,3,0,4,1)
n=6: (0,1,4,2,3,5)
n=7: (1,2,4,5,6,0,3)
n=8: (1,3,0,7,4,6,2,5)
```

Replacing `S` by the displayed transitive tournament gives

```text
                 n=5   n=6   n=7   n=8
stored B          16    18    22    28
ordered B         16    18    26    32
```

Thus ordering happens to retain an optimum at orders five and six, but the
few nontransitive edges lower the sharp finite objective by four at orders
seven and eight.  This is compatible with (8), and it is finite information
only.

There is also no common Paley/conference description of these four stored
orientations.

* At `n=8`, the stored `R` is exactly a skew conference matrix:
  `R^2=-7I`.  Switching it by
  `(1,-1,1,-1,1,-1,-1,1)` and deleting vertex zero leaves the Paley
  tournament on `F_7`; the remaining vertices `1,...,7` receive field labels
  `(0,1,6,2,3,5,4)`, with positive differences `{1,2,4}`.
* At `n=6`, instead, `A^2=5I` and the stored orientation attains equality in
  the Clifford parity floor:
  `||AR+RA||_F^2=48=2n(n-2)`.  The support of `AR+RA` is two disjoint
  triangles, with every nonzero entry of magnitude two.
* Switching and permutation preserve the characteristic polynomial of a
  skew tournament matrix.  At `n=7`, the two relevant polynomials are
  `chi_R(t)=t(t^2+3)(t^2+7)(t^2+11)` and
  `chi_S(t)=t(t^2+7)(t^4+14t^2+1)`; at `n=8` they are
  `chi_R(t)=(t^2+7)^4` and
  `chi_S(t)=(t^2+7)(t^6+21t^4+35t^2+7)`.  Hence `R` and `S` are not even
  switching-permutation equivalent at either order.

The structural conclusion is therefore narrow but rigorous: the stored
solutions support a **near-ordering** hypothesis, not an exact ordering or a
single Paley recipe.  Equations (8)--(9) show precisely what an all-order
version of that hypothesis would buy.  They do not provide the missing
bound on `tau(S_*)`.
