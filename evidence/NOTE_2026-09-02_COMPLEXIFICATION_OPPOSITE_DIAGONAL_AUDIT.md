# Complexification and the opposite-diagonal doubling lift

**Status:** proved the exact all-order diamond, its four-label weighted-slice
form, and an exact counterexample.  The opposite-diagonal block construction
gives a genuinely different open doubling reduction.  Its all-directed
specialization isolates one outgoing half of every signed cut; its coherent
holomorphic specialization is instead a clique-flip problem, and a globally
optimal order-four signing disproves the zero-loss lemma that would close that
specialization.  General `A`-dependent cross blocks remain open.  No finite
census was used.

Let `A` be a symmetric zero-diagonal sign matrix of order `n`, put

\[
 Q_A(x)=\sum_{i<j}A_{ij}x_ix_j,\qquad M=\Phi(A),
\]

and let `C` be an arbitrary `n` by `n` sign matrix.  The order-`2n` signing

\[
 {\cal K}(A,C)=
 \begin{pmatrix}A&C\\ C^T&-A\end{pmatrix}                    \tag{1}
\]

has energy

\[
 Q_A(x)-Q_A(y)+x^TCy.
\]

## 1. Exact opposite-diagonal diamond

Flipping all signs of `x` fixes the two internal quadratic terms and negates
the cross term.  Therefore

\[
 \boxed{
 \Phi({\cal K}(A,C))
 =\max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TCy|\bigr).
 }                                                            \tag{2}
\]

In particular, define

\[
 \Gamma_-(A)=\min_{C\in\{\pm1\}^{n\times n}}
 \max_{x,y}\bigl(|Q_A(x)-Q_A(y)|+|x^TCy|\bigr).               \tag{3}
\]

For every globally optimal `A`, the estimate

\[
 \Gamma_-(A)\le 2\sqrt2\,M+o_{\rm Dini}(n^{3/2})              \tag{4}
\]

would prove the multiplier-two amplification.  This is distinct from the
existing skew completion, whose internal term is `Q_A(x)+Q_A(y)` and whose
cross matrix is skew.  The fractional point `C=0` has objective

\[
 \max_{x,y}|Q_A(x)-Q_A(y)|=P(A)+N(A),                          \tag{5}
\]

where `P(A)=max_x Q_A(x)` and `N(A)=-min_x Q_A(x)`.  Rounding
that point is again a simultaneous, shifted discrepancy problem; (5) alone
contains no integral cross-block information.

## 2. Every cross block is an exact four-label graph

The general cross block has an exact geometric form.  For
`T subseteq [n]`, put `r_i=-1` on `T` and `r_i=1` off `T`, and let
`D_T=diag(r)`.  Define the symmetric zero-diagonal weighted graph

\[
 (G_T)_{ij}
 =A_{ij}{1-r_ir_j\over2}+{C_{ij}r_j+C_{ji}r_i\over2},
 \qquad
 h_T={1\over2}\sum_i C_{ii}r_i.                            \tag{5a}
\]

Every Boolean pair is uniquely `x=s,y=D_Ts`.  Substitution in (1) gives

\[
 Q_{{\cal K}(A,C)}(s,D_Ts)=2\bigl(Q_{G_T}(s)+h_T\bigr),
\]

and therefore

\[
 \boxed{
 \Phi({\cal K}(A,C))=2\max_{T,s}|Q_{G_T}(s)+h_T|.
 }                                                           \tag{5b}
\]

Let `A_delta(T)` retain the `A`-edge on a pair precisely when the pair
crosses `(T,T^c)`, and be zero otherwise.  The first summand of (5a) is
`A_delta(T)`.  Complementing `T` fixes it and negates both the second
summand and `h_T`.  Pairing the two complementary values via
`max(|a+b|,|a-b|)=|a|+|b|` yields the second exact form

\[
 \boxed{
 {1\over2}\Phi({\cal K}(A,C))
 =\max_{T,s}\left(
 |Q_{A_\delta(T)}(s)|+{1\over2}|s^TCD_Ts|\right).
 }                                                           \tag{5c}
\]

This is the precise midpoint--displacement geometry of the
opposite-diagonal lift.  It is not a disk or triangle relaxation.

### The four labels on one edge

For each unordered pair `{i,j}`, the two cross signs have only four
possibilities.

| label | cross signs | weight in `G_T` |
|---|---|---|
| undirected `epsilon` | `C_ij=C_ji=epsilon` | `epsilon` if both vertices are outside, `-epsilon` if both are inside, and `A_ij` if the edge crosses |
| directed `i -> j` | `C_ij=A_ij`, `C_ji=-A_ij` | `2A_ij` if `i in T,j notin T`, and zero otherwise |
| directed `j -> i` | `C_ji=A_ij`, `C_ij=-A_ij` | `2A_ij` if `j in T,i notin T`, and zero otherwise |

The first row has the two choices `epsilon=+1,-1`, so the table contains all
four ordered sign pairs.  In particular, `G_T` is an ordinary signing for
some, equivalently every, `T` if and only if the off-diagonal part of `C` is
symmetric.  Any skew pair has weight zero or `+/-2` for every cut.  Thus the
additional two labels are genuine weighted-graph freedom, and one may not
invoke global minimality of `A` on the resulting slices.

The field term is lower order, so (5b) also gives, with
`L(A,C)=max_T Phi(G_T)`,

\[
 \left|\Phi({\cal K}(A,C))-2L(A,C)\right|\le n.              \tag{5d}
\]

### The displacement cannot be bounded on its own

As `(T,s)` varies, `(s,D_Ts)` ranges over every Boolean pair.  Hence

\[
 \max_{T,s}{1\over2}|s^TCD_Ts|
 ={1\over2}\max_{x,y}|x^TCy|
 ={1\over2}\|C\|_{\infty\to1}.                             \tag{5e}
\]

Averaging over `y` and then choosing the row signs gives

\[
 \|C\|_{\infty\to1}=\max_y\|Cy\|_1
 \ge n\mu_n
 =\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}.                   \tag{5f}
\]

Thus the displacement in (5c) is always at least
`(1/sqrt(2 pi)+o(1)) n^(3/2)`.  The independent allowance in the desired
diamond is only `(sqrt(2)-1)M`, at most
`((sqrt(2)-1)/2+o(1))n^(3/2)` along the known upper construction.  Since
`1/sqrt(2 pi)>(sqrt(2)-1)/2`, an uncoupled estimate is impossible.  The
midpoint and displacement must be anticorrelated state by state.

This formulation also shows that the critical `sqrt(2)` cannot be improved
uniformly.  Put

\[
 \widehat\Lambda_n=\min_{\Phi(A)=m_n}\min_C
        \max_{T,s}|Q_{G_T}(s)+h_T|.
\]

Since every `K(A,C)` is an order-`2n` signing, (5b) gives
`widehat Lambda_n >= m_(2n)/2`.  If an eventual dyadic tail satisfied
`widehat Lambda_n <= c m_n+o(n^(3/2))` for a fixed `c<sqrt(2)`, then
`m_(2n)<=2c m_n+o(n^(3/2))`; normalized iteration would force
`m_n/n^(3/2)` to zero on that tail, contradicting the universal positive
lower bound.  The same conclusion for the unfielded `L(A,C)` follows from
(5d).

### All-directed cross blocks: one outgoing half of every cut

The cleanest nonsymmetric specialization uses only the two directed labels.
Let `S` be a tournament and put `R=A circ S`, so `R` is skew.  Choose any
diagonal signing `D` and use `C=R+D`.  Define

\[
 D_{\to}(A,S)=\max_{T,s}
 \left|\sum_{\substack{u\in T,\ v\notin T\\u\to v}}
       A_{uv}s_us_v\right|.                                 \tag{5g}
\]

The table says exactly that `G_T` has weight `2A_uv` on the arcs leaving
`T` and zero on every other edge.  Taking the skew cross block `C=R` first
gives `h_T=0`; although the matching entries in the resulting order-`2n`
matrix are zero, the energy identity is exact:

\[
 \boxed{
 4D_{\to}(A,S)=\max_{x,y}
 \bigl(|Q_A(x)-Q_A(y)|+|x^TRy|\bigr).
 }                                                           \tag{5h}
\]

Adding the diagonal signing `D` fills the matching and changes every energy
by at most `n`, so

\[
 \left|\Phi({\cal K}(A,R+D))-4D_{\to}(A,S)\right|\le n.     \tag{5h'}
\]

Complementing `T` includes the incoming direction as another outgoing
case.  Therefore the following is now a sufficient all-orders target:

\[
 \boxed{
 \text{orient the edges so that}\quad
 D_{\to}(A,S)\le {M\over\sqrt2}+o_{\rm Dini}(n^{3/2}).
 }                                                           \tag{5i}
\]

It is close in language to the skew half-cut neighbor, but mathematically
distinct: (5i) controls the signed energy of the outward directed half
itself, whereas Proposition 6.5a controls the norm after reversing that
half.  This exact outgoing-half formulation is the main surviving lead from
the four-label audit.

It is also an exact simultaneous paving problem.  For each `T`, let
`B_T(S)` be the rectangular matrix, with rows in `T` and columns in `T^c`,
whose entry is `A_uv` if the arc points from `u` to `v`, and zero if it
points back.  Then

\[
 D_{\to}(A,S)=\max_T\|B_T(S)\|_{\infty\to1}.                 \tag{5i'}
\]

The unpaved cut matrix obeys the exact hereditary bound

\[
 \|A_{T,T^c}\|_{\infty\to1}\le M.                          \tag{5i''}
\]

Indeed arbitrary row and column signs extend to two Boolean states differing
exactly on `T`; their endpoint diamond bounds the cross energy by `M`.
Thus (5i) asks one tournament to give a one-sided `1/sqrt(2)` paving of
every signed cut submatrix simultaneously.  The consistency requirement
between overlapping cuts is what is absent from standard one-matrix paving.

For one fixed face it is elementary.  Write `C_cut=F+G` for its total signed
cut energy and let `k=|T|(n-|T|)`.  If the fixed spin state gives `p`
positive and `q` negative cross-edge contributions, assigning edges to the
outgoing half lets `F` attain every integer in `[-q,p]`, with
`G=C_cut-F`.  Hence

\[
 \min_{\text{orientations of this face}}\max(|F|,|G|)
 =\left\lceil {|C_{\rm cut}|\over2}\right\rceil,            \tag{5j}
\]

whereas the maximum over such assignments is
`(k+|C_cut|)/2`.  The endpoint diamond gives
`|C_cut|<=M-|I|`, so every individual face is easy, at cost at most
`ceil(M/2)`.  The difficulty is solely to make one orientation realize
these balanced choices simultaneously on all faces.  Scalar plotting sees
the strip in the sum direction `F+G`; it supplies no control in the
cancellation direction `F-G`.

There is a concrete necessary best-response gate.  For each `y`, choose
`X_R(y)` coordinatewise so that
`X_R(y)^T R y=||Ry||_1` (with any fixed convention at zero).  The target
(5h)--(5i) would force, pointwise,

\[
 |Q_A(X_R(y))-Q_A(y)|+\|Ry\|_1
 \le2\sqrt2M+o_{\rm Dini}(n^{3/2}).                         \tag{5k}
\]

For uniform `Y`, every row of `R` is a sum of `n-1` independent signs, so

\[
 \mathbb E\|RY\|_1=n\mu_{n-1}.
\]

Consequently any successful orientation necessarily has

\[
 \mathbb E|Q_A(X_R(Y))-Q_A(Y)|
 \le2\sqrt2M-n\mu_{n-1}+o_{\rm Dini}(n^{3/2}).              \tag{5l}
\]

This does not yet construct `R`, but it is a falsifiable joint
energy--skew criterion; an estimate of the two summands separately cannot
reach it.

## 3. Symmetric cross blocks give a two-minimizer interpolation

There is a broader exact Hadamard subclass before imposing holomorphicity.
Let `C_0` be any symmetric zero-diagonal signing and use the cross block
`C=C_0+D`, where `D` is an arbitrary diagonal signing.  For
`T \subseteq [n]`, define the symmetric signing `B_T(A,C_0)` by

\[
 (B_T)_{ij}=\begin{cases}
 (C_0)_{ij},&i,j\notin T,\\
 A_{ij},&|\{i,j\}\cap T|=1,\\
 -(C_0)_{ij},&i,j\in T.
 \end{cases}                                                  \tag{6}
\]

After reordering into two-vertex clouds, the inter-cloud block is

\[
 \begin{pmatrix}A_{ij}&(C_0)_{ij}\\(C_0)_{ij}&-A_{ij}\end{pmatrix}.
\]

For the cloud state `s_i(1,(-1)^{t_i})`, direct multiplication gives

\[
 Q_{{\cal K}(A,C_0+D)}
 =2Q_{B_T(A,C_0)}(s)+\sum_i d_i(-1)^{t_i},
 \qquad T=\{i:t_i=1\}.                                      \tag{7}
\]

Consequently

\[
 \boxed{
 \left|\Phi({\cal K}(A,C_0+D))
 -2\max_T\Phi(B_T(A,C_0))\right|\le n.
 }                                                            \tag{8}
\]

If `A` and `C_0` are both global minimizers, the endpoint slices
`T=\varnothing` and `T=[n]` are `C_0` and `-C_0`, both of norm `M`.
Thus a concrete restricted form of the new live problem is to choose two
minimizers whose every hybrid slice (6) has norm at most
`sqrt(2) M+o_Dini(n^(3/2))`.  Global minimality again gives only the
reverse inequality for those slices.  This two-minimizer interpolation is
not supplied by complexification, but it is an exact `A`-dependent
construction inside (3).

## 4. Holomorphic complexification forces the coherent cross block

Write `z=x+iy`.  The unique holomorphic quadratic extension of `Q_A` is

\[
 \widetilde Q_A(z)=\sum_{i<j}A_{ij}z_iz_j
 =Q_A(x)-Q_A(y)+i\,x^TAy.                                    \tag{9}
\]

This choice is forced.  If a function with real part
`Q_A(x)-Q_A(y)` and imaginary part `x^TCy` is holomorphic, the
Cauchy--Riemann equations give

\[
 Ax=C^Tx\quad\hbox{and}\quad Ay=Cy
\]

for all real `x,y`; hence `C=A`.  Since a cross block also has `n`
diagonal entries whereas `A_ii=0`, the actual signing is uniquely of the
form

\[
 C=A+D,\qquad D=\operatorname{diag}(d_1,\ldots,d_n),
 \quad d_i\in\{\pm1\},                                     \tag{10}
\]

up to an additive matching term of magnitude at most `n`.

This also explains precisely why ordinary complexification does not choose
the imaginary signs in the directed-half-cut formulation.  If `R` is skew,
then

\[
 z^T(A+iR)z=z^TAz                                             \tag{11}
\]

because `z^TRz=0`.  All holomorphic polynomial and symmetric-polarization
theorems are consequently blind to `R`.  The live fourth-phase expression
is instead the nonholomorphic Hermitian form `z^*(A+iR)z`.

## 5. Exact clique-flip form of the coherent lift

For `T \subseteq [n]`, let `A^{K_T}` be obtained by reversing the signs of
all edges with both ends in `T`, and put

\[
 L_{\rm cl}(A)=\max_{T\subseteq[n]}\Phi(A^{K_T}).             \tag{12}
\]

Reorder the vertices of `K(A,A+D)` into the `n` two-vertex clouds.  The
inter-cloud block on `{i,j}` is

\[
 A_{ij}\begin{pmatrix}1&1\\1&-1\end{pmatrix}.               \tag{13}
\]

Every cloud state has the unique form
`s_i(1,(-1)^{t_i})`, with `s_i in {+1,-1}` and
`t_i in {0,1}`.  Direct multiplication in (13) gives

\[
 Q_{{\cal K}(A,A+D)}
 =2Q_{A^{K_T}}(s)+\sum_i d_i(-1)^{t_i},
 \qquad T=\{i:t_i=1\}.                                      \tag{14}
\]

Thus, for every choice of the diagonal matching,

\[
 \boxed{
 2L_{\rm cl}(A)-n
 \le \Phi({\cal K}(A,A+D))
 \le 2L_{\rm cl}(A)+n.
 }                                                            \tag{15}
\]

There is an exact phase interpretation.  For
`z in {+1,-1,+i,-i}^n`, write `z_i=s_i` off `T` and
`z_i=i s_i` on `T`.  If `O,J,X` are respectively the `A`-energies
on `T^c`, on `T`, and across the cut, then

\[
 \widetilde Q_A(z)=O-J+iX.                                   \tag{16}
\]

Changing all `s_i` on `T` negates `X` and preserves `O,J`.
Since `max(|a+b|,|a-b|)=|a|+|b|`, maximizing first over that pair
and then over `s,T` proves

\[
 \boxed{
 L_{\rm cl}(A)=
 \max_{z\in\{\pm1,\pm i\}^n}
 \bigl(|\operatorname{Re}\widetilde Q_A(z)|
      +|\operatorname{Im}\widetilde Q_A(z)|\bigr).
 }                                                            \tag{17}
\]

If

\[
 \Phi_4(A)=\max_{z\in\{\pm1,\pm i\}^n}|\widetilde Q_A(z)|,
\]

then (17) gives only

\[
 \Phi_4(A)\le L_{\rm cl}(A)\le\sqrt2\,\Phi_4(A).           \tag{18}
\]

Consequently the usual norm-only complexification route would need the
lossless special statement `Phi_4(A)<=Phi(A)` for globally optimal
complete signings.  Since the real cube is contained in the fourth-phase
cube, that statement is equality.  General complexification theorems do not
give it: the degree-two case of the classical Visser bound is

\[
 \|\widetilde Q_A\|_{\mathbb T^n}\le2\|Q_A\|_{[-1,1]^n}=2M. \tag{19}
\]

Combining (18)--(19) loses a full factor two relative to (4).  The relevant
primary modern statements are Proposition 18 of Muñoz-Fernández,
Sarantopoulos, and Tonge, [“Complexifications of real Banach spaces,
polynomials and multilinear maps”](https://doi.org/10.4064/sm-134-1-1-33),
and Lemma 1.3(3) of Defant, Mastyło, and Pérez,
[“On the Fourier spectrum of functions on Boolean cubes”](https://arxiv.org/abs/1706.03670).
For `l_infinity^n`, the Taylor complexification norm used there is exactly
the standard complex sup norm.

### Exact two-variable face calculus

The continuous geometry of one phase face can be exhausted exactly.  Fix
`T` and a sign state `s`, and let `O,J,X` be its energies on `T^c`, on `T`,
and across the cut.  Scaling the two signed partial vectors by
`r,t in [-1,1]` gives

\[
 p(r,t)=Or^2+Jt^2+Xrt.
\]

The Boolean norm bounds this polynomial on the entire square, by
multilinearity.  In fact the complete semialgebraic criterion is

\[
\begin{aligned}
 |O+J|+|X|&\le M,\\
 |O-X^2/(4J)|&\le M
   &&\text{if }J\ne0\text{ and }|X|\le2|J|,\\
 |J-X^2/(4O)|&\le M
   &&\text{if }O\ne0\text{ and }|X|\le2|O|.
\end{aligned}                                               \tag{18a0}
\]

These conditions are necessary and sufficient for `|p|<=M` on the square.
Homogeneity moves every nonzero point to one of the four sides.  On
`r=+/-1`, the endpoints give the first line and the only possible interior
quadratic vertex gives the second; on `t=+/-1` one obtains the third.

If the two opposite corners are exact opposite extrema,
`p(1,1)=M,p(1,-1)=-M`, then (18a0) reduces to

\[
 O+J=0,\qquad X=M,qquad |O|=|J|\le M/2.                    \tag{18a1}
\]

On the unit circle,

\[
 p(\cos\theta,\sin\theta)
 ={O+J\over2}+{O-J\over2}\cos2\theta+{X\over2}\sin2\theta.
\]

Exact maximization in `theta` therefore gives

\[
 \boxed{|O+J|+\sqrt{(O-J)^2+X^2}\le2M.}             \tag{18a}
\]

There is a universal but lower-order axial improvement.  If `k=|T|`, then

\[
 |O|\le M-\lfloor k/2\rfloor,
 \qquad
 |J|\le M-\lfloor(n-k)/2\rfloor.                   \tag{18b}
\]

To prove it, expose the vertices of any order-`r` signing sequentially and
choose each new spin so that the new incident sum is nonnegative.  Whenever
the number of earlier vertices is odd that sum has odd parity and magnitude
at least one, so the resulting one-sided maximum is at least
`floor(r/2)`.  Applying the same statement to the negative signing supplies
a one-sided minimum at most `-floor(r/2)`.  Now fix the state giving `O` and
complete the other `k` coordinates.  Their remainder has the form
`L(x)+Q_B(x)`.  Pairing `x` with `-x` removes the linear term in the relevant
maximum/minimum, so the remainder reaches both at least `floor(k/2)` and at
most `-floor(k/2)`.  The full norm bound proves the first inequality in
(18b); the second is symmetric.

This is the complete scalar-calculus gain, and it is insufficient at leading
order.  The polynomial

\[
 p(r,t)=M(r^2-t^2)                                      \tag{18c}
\]

is bounded by `M` on the square and attains equality in (18a0)--(18a).
Equation (18b) removes only `Theta(n)` from its two axial coefficients,
whereas the required scale is `Theta(n^(3/2))`.  Moreover, in the directed
split `X=F+G`, even the full criterion (18a0) contains no information about
`F-G`.  Hence a proof cannot come from plotting or differentiating one face
alone: it must couple different faces or insert new global-minimizer
information.

### An honest sign family saturates the scalar obstruction

The saddle obstruction persists asymptotically inside actual sign matrices,
even after imposing strict one-edge local minimality.  Let `r>=4` admit a
Hadamard matrix `C`, put `k=2r`,

\[
 P=J_k-I_k,\qquad
 K_2=\begin{pmatrix}1&-1\\-1&1\end{pmatrix},\qquad
 H=C\otimes K_2,
\]

and form the order-`4r` signing

\[
 A_r=\begin{pmatrix}P&H\\H^T&-P\end{pmatrix}.               \tag{18d}
\]

Then

\[
 \boxed{\Phi(A_r)=2r^2+2,}                                  \tag{18e}
\]

and flipping any one edge raises the norm to at least `2r^2+4`.

For the calculation, pair the coordinates in each `k`-block and write

\[
 s_h={u_{h,+}+u_{h,-}\over2},\quad
 d_h={u_{h,+}-u_{h,-}\over2},
\]

and similarly `t,e` for the other block.  Exactly one of `s_h,d_h` is a
sign and the other is zero.  With
`a=sum d_h^2`, `b=sum e_h^2`, `S=sum u`, and `T=sum v`, direct expansion is

\[
 Q_{A_r}(u,v)={S^2-T^2\over2}+4d^TCe.                       \tag{18f}
\]

Here `|S|<=2(r-a)`, `|T|<=2(r-b)`, and Hadamard orthogonality gives
`|d^TCe|<=sqrt(rab)`.  Therefore

\[
 Q_{A_r}\le2(r-a)^2+4r\sqrt a\le2r^2+2.                    \tag{18g}
\]

For `a=0,1` this is immediate.  Subtracting the right side for `a>=2`
gives `2a^2-4r(a-sqrt(a))-2`, which is negative directly at `a=2,3` for
`r>=4`, and is at most
`-2a^(3/2)(sqrt(a)-2)-2` for `a>=4`.  Apply the same estimate to `-Q`, with
`a,b` interchanged and `C^T` in place of `C`; the transpose is also
Hadamard, so this gives the lower bound.  Equality follows by taking
`d=e_i`, all other `s_h=1`,
`e=C^T e_i`, and `t=0`; the negative equality is symmetric.

For strict edge stability, the equality families retain two free global
signs.  A positive maximizer centered on pair `i` has

\[
 u_{i,\alpha}=\sigma\alpha,\quad
 u_{h,\alpha}=\eta\ (h\ne i),\quad
 v_{j,\beta}=\sigma C_{ij}\beta.
\]

A symmetric negative family is centered on a pair of the second block.
Choosing the center and the product `sigma eta` makes any prescribed edge
contribute opposite to the sign of the corresponding extremal energy.
Flipping that edge therefore changes `+M` to `M+2`, or `-M` to `-M-2`.

Finally, on the natural two-block face at the all-one state,

\[
 O=2r^2-r,\qquad J=-(2r^2-r),\qquad X=0.                    \tag{18h}
\]

Thus the ratio in (18a) tends to equality.  The essential limitation is that
`Phi(A_r)=Theta((4r)^2)`, whereas global minimizers have
`Theta(n^(3/2))` norm.  This family blocks any proof using only scalar face
calculus plus the inequalities `Phi(A^e)>=Phi(A)`; it does **not** block a
proof using the optimal scale or genuinely nonlocal global minimality.

## 6. A globally optimal signing already has complex and clique inflation

The lossless lemma is false even under global optimality.  Take

\[
 A=\begin{pmatrix}
 0&1&1&1\\
 1&0&-1&1\\
 1&-1&0&1\\
 1&1&1&0
 \end{pmatrix}.                                               \tag{20}
\]

After fixing the first Boolean coordinate to one,

\[
 Q_A(1,a,b,c)=a+b+c-ab+ac+bc.
\]

For `(a,b)=(1,1),(1,-1),(-1,1),(-1,-1)`, this is respectively
`1+3c`, `1+c`, `1+c`, and `-3-c`.  Hence `Phi(A)=4`.
Every order-four quadratic signing is a sum of six signs, so all its values
are even, while Walsh orthogonality gives `E Q^2=6`.  A norm below four
would force `Q^2<=4` everywhere, a contradiction.  Therefore

\[
 \Phi(A)=m_4=4.                                               \tag{21}
\]

But for `z=(1,i,i,1)`,

\[
 \widetilde Q_A(z)=2+4i,
 \qquad |\widetilde Q_A(z)|=2\sqrt5>4.                       \tag{22}
\]

Even more directly, taking `T={2,3}` reverses the unique negative edge in
(20), producing the all-positive signing.  Thus

\[
 L_{\rm cl}(A)=6>\sqrt2\,\Phi(A).                            \tag{23}
\]

This is a counterexample to the zero-error claim that global minimality
controls every coherent clique flip.  It does not refute an asymptotic
`o(n^(3/2))` statement, because it is one fixed order and the matching error
in (15) is lower order only along growing orders.

## 7. Why sharp Grothendieck or vector balancing does not fill the gap

Krivine's sharp rank-two Grothendieck theorem (`K_G(2)=sqrt(2)`) concerns
one fixed coefficient matrix and has quantifiers

\[
 \forall W\quad \exists\hbox{ vertex signs rounding the one objective }W.
\]

The completion problem has the opposite simultaneous pattern

\[
 \exists C\quad\forall(x,y)\quad
 |Q_A(x)-Q_A(y)|+|x^TCy|\le2\sqrt2M+o(n^{3/2}).               \tag{24}
\]

The rounding signs in the former depend on the witness matrix; they do not
produce one sign for each entry of `C` that works for every witness.  The
source for the exact rank-two constant is Krivine,
[“Constantes de Grothendieck et fonctions de type positif sur les sphères”](https://doi.org/10.1016/0001-8708(79)90017-3).

Nor can a general vector-balancing theorem round `C=0` with lower-order
uniform error.  Every `n` by `n` sign matrix satisfies

\[
 \max_{x,y}|x^TCy|=\max_y\|Cy\|_1
 \ge \mathbb E_y\|Cy\|_1=n\mu_n
 =\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}.                     \tag{25}
\]

Thus the exact constraint system itself has leading-order discrepancy.
Spencer/Banaszczyk bounds of order `n^(3/2)` cannot be put into an error
term, and a separate cross-norm budget is too large.  As in the skew route,
any successful general `C` must correlate its large cross values with
states having internal difference below the endpoint maximum.

## 8. Verdict

The opposite-diagonal lift supplies the new exact open reduction (4).
Symmetric cross blocks give the exact two-minimizer interpolation (6)--(8).
The holomorphic specialization is completely rigid and reduces to
(12)--(17); standard complexification loses a factor two, and (20)--(23)
rules out the tempting zero-loss lemma based solely on global optimality.  Sharp
rank-two Grothendieck rounding has the wrong quantifier order, while generic
vector balancing has an unavoidable leading-order discrepancy floor.

What remains live is narrower and genuinely new: either construct a
noncoherent, `A`-dependent sign matrix `C` satisfying (4), or prove the
concrete all-directed estimate (5i).  Both require statewise correlation
between the internal cut midpoint and its cross displacement.  None of the
audited black-box complexification, polarization, Grothendieck, or
vector-balancing theorems provides that correlation.  The multiplier-two
ray and the original MO limit remain open.
