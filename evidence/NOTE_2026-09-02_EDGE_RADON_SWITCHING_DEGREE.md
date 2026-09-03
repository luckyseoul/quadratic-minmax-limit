# Edge-Radon switches are necessarily global

This note concerns exactly the integral edge-Radon map `R` of Proposition
15.760.  It gives a uniform obstruction to a bounded local-switch approach.
It does **not** prove that the signed Boolean fibre is nonempty.

## Row-collision lemma

Write an edge of `V=F_p^2` uniquely as

\[
 e=(a,[\delta])=\{a-\delta,a+\delta\},
 \qquad a\in V,\quad [\delta]\in(V\setminus\{0\})/\{\pm1\}.
\]

For a projective nonzero functional `L`, the row of `R` containing `e` is

\[
 \rho_L(e)=
 \begin{cases}
 P_L,&L(\delta)=0,\\
 K_L\bigl(L(a),L(\delta)^2\bigr),&L(\delta)\ne0.
 \end{cases}                                                     \tag{1}
\]

**Lemma.**  If `e` and `f` are distinct edges, then

\[
 \#\{L:\rho_L(e)=\rho_L(f)\}\le 2.                              \tag{2}
\]

**Proof.**  Put `f=(b,[\epsilon])`.

If `[\delta]=[\epsilon]`, choose `\epsilon=\delta`.  There is one possible
parallel collision, namely the unique projective `L` annihilating `\delta`.
Every nonparallel collision is equivalent, by (1), to `L(a-b)=0`.  Since
`e\ne f`, one has `a-b\ne0`, so this names at most one further projective
direction.

If `[\delta]\ne[\epsilon]` but `\delta,\epsilon` are collinear, write
`\epsilon=c\delta` with `c\ne\pm1`.  Their unique common parallel direction
gives one collision, while outside that direction
`L(\epsilon)^2=c^2L(\delta)^2\ne L(\delta)^2`, so there is no nonparallel
collision.

It remains to consider linearly independent `\delta,\epsilon`.  A parallel
collision is then impossible.  A nonparallel collision requires

\[
 L(a-b)=0,\qquad L(\delta)^2=L(\epsilon)^2.
\]

The second equality says either `L(\delta-\epsilon)=0` or
`L(\delta+\epsilon)=0`.  Both `\delta-\epsilon` and
`\delta+\epsilon` are nonzero.  Each sign therefore supplies at most one
projective direction, and imposing `L(a-b)=0` can only remove it.  This
proves (2).  \(\square\)

## Minimum support and degree

Let `0\ne z\in\ker_{\mathbb Z}R`, and let

\[
 S_+=\{e:z_e>0\},\qquad S_-=\{e:z_e<0\}.
\]

Fix `e\in S_+`.  In each of the `p+1` directions, the sum in the row
containing `e` is zero.  Hence that row also contains some edge of `S_-`.
By (2), any one edge of `S_-` can do this in at most two directions.  Thus

\[
 |S_-|\ge \left\lceil\frac{p+1}{2}\right\rceil=\frac{p+1}{2}.
\]

Interchanging the signs gives the same bound for `S_+`.  Consequently every
nonzero integral edge-Radon switch satisfies

\[
 |\operatorname{supp}z|\ge p+1.                                \tag{3}
\]

For the unsigned matrix, summing all rows in any one direction gives
`\sum_ez_e=0`.  Therefore its Markov degree obeys

\[
 \deg z:=\sum_{z_e>0}z_e=-\sum_{z_e<0}z_e
 \ge \frac{p+1}{2}.                                             \tag{4}
\]

In particular there is no prime-uniform bounded-support or bounded-degree
Markov/switching basis for `R`.  After the source signing of Proposition
15.760, `A=R\operatorname{diag}(\tau)`, support is unchanged.  Thus if two
distinct Boolean graphs `x,x'` have the same normalized edge-Radon data,
then

\[
 d_H(x,x')=|\operatorname{supp}(x-x')|\ge p+1.                  \tag{5}
\]

So a `2x2`, four-edge, or any other constant-size toggle cannot preserve the
full edge-Radon data.  Any exact repair of a Boolean lift must change a
number of graph edges growing at least linearly with `p`.

## Audit of the centered near-pencil symmetry shortcut

There is no nontrivial affine symmetry which fixes the **indexed Radon
rows** pointwise.  Indeed, let `g(a)=Aa+c` and suppose its edge permutation
fixes every set `P_L` and `K_L(alpha,beta)`.  Preservation of every `P_L`
forces `A` to preserve every one-dimensional kernel of a projective
functional, hence `A=lambda I`.  Applying (1) to arbitrary midpoints in any
nonparallel cell then gives

\[
 L(\lambda a+c)=L(a)\quad\hbox{for all }a,L,
\]

so `lambda=1`, `c=0`, and `g` is the identity.  In particular centering the
circle supplies no hidden translation: a translation also moves its unique
centre unless it is zero.

A symmetry may instead permute Radon rows and stabilize one specially
symmetric numerical target.  If `g` does so and `H` is a Boolean lift, then
`gH` is another lift; (5) gives the exact dichotomy

\[
 gH=H\quad\hbox{or}\quad |H\mathbin\triangle gH|\ge p+1.       \tag{6}
\]

The near-pencil surplus does not provide the opposite strict upper bound.
Although all but at most `3s` edges of `H` are isolated matching edges, two
matchings on the same odd-boundary set may differ on an alternating cycle of
arbitrary even length.  Thus `|H symmetric-difference gH|` is not bounded by
`6s` unless
one separately proves that `g` fixes the large isolated matching.  The
boundary and directional quotas currently do not do that.

Likewise, the compact atoms of Proposition 15.758 specify target rows, not
a source lift.  Proposition 15.760 proves existence of some signed integral
preimage only after all moment rows pass; it neither chooses a canonical
preimage nor bounds its Hamming distance from the signed Boolean box.  The
distance theorem does give a useful conditional uniqueness radius: every
Hamming ball of radius at most `(p-1)/2` contains at most one Boolean lift of
a fixed target.  A centered-target application therefore still needs one
of the following genuinely new inputs:

1. an explicit Boolean lift (which would settle existence), or
2. a proof that every Boolean lift lies in such a ball around a specified
   signed integral lift, together with identification of the sole candidate.

Neither follows from the centered circle, the two-outlier boundary, or the
compact target atoms presently recorded.

## What Markov-basis theory does and does not supply

Diaconis--Sturmfels, Theorem 3.1, identifies Markov bases with binomial
generators of the toric ideal and proves connectivity of every *nonnegative
fibre*:

<https://doi.org/10.1214/aos/1030563990>

For upper-bounded fibres, the standard construction is the Lawrence lift;
see Rapallo--Yoshida:

<https://arxiv.org/abs/0905.4841>

These theorems connect two points already known to lie in the nonnegative or
bounded fibre.  They do not turn an unrestricted signed integral preimage
into a Boolean one, nor do they prove that the Boolean fibre is nonempty.
The familiar `2x2` bounded-table theorem applies to one row/column-margin
matrix, not to the simultaneous `p+1` edge projections (1); (3) directly
rules out such local switches here.

Even saturation and absence of Smith torsion are insufficient in general.
For example

\[
 A=\begin{pmatrix}1&1&1\\1&0&0\\0&1&0\end{pmatrix},\qquad
 b=\begin{pmatrix}1\\1\\1\end{pmatrix}
\]

have `det A=1`, while the unique integral preimage is
`(1,1,-1)^t`; no Boolean preimage exists.  Thus Proposition 15.760's
integral-image theorem cannot by itself imply the desired signed Boolean
lift.  What remains is a special box-nonemptiness theorem for this particular
edge-Radon target.  Equations (3)--(5) show that it cannot be proved by a
prime-uniform collection of local toggles.

**Status:** proved kernel-geometry obstruction; residual (ii) remains open.
