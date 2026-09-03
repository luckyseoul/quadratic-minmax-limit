# The exact low-weight gap and the structured Mobius puncture

Date: 2026-09-03

Status: the grouped uncertainty theorem upgrades the halved row code to

\[
 d(\operatorname {Row}(D))=ph,
\]

its minimum words are exactly the fixed-transverse rectangles, and there
are no row-code weights strictly between `ph` and
`N=|Delta|=(p+1)h`.  Consequently every branch-C Mobius puncture which can
pass the elementary physical-size condition leaves the halved **binary** map
onto.  This settles the structured punctured parity map, not the prescribed
Hamming slice, the directionwise integer slices, the divided integral
Boolean fibre, residual (ii), E1, or the original limit.

## 1. Input and a two-partition distance lemma

Put

\[
 p=2h+1,\qquad d=p+1,\qquad
 \Delta=(\mathbf F_p^2\setminus\{0\})/\{\pm1\},\qquad N=dh.
\]

For a projective functional `A`, the partition of `Delta` attached to the
block space `B_A` consists of the radial line `l_A`, of size `h`, and `h`
paired non-origin affine blocks, each of size `p`.  Thus every nonempty union
of cells has size at least `h`, with equality only for `l_A`.

There is a second elementary fact which will be used sharply.  If `A` and
`B` are distinct directions and `S_A,S_B` are nonempty proper unions of
their respective cells, then

\[
                    |S_A\mathbin\triangle S_B|\ge2h.       \tag{1}
\]

To prove it, use the independent coordinates `x=A(v), y=B(v)`.  Write the
two cell unions as even Boolean functions `f(x),g(y)` on `F_p`, and put
`a=|supp(f)|`, `b=|supp(g)|`.  Proper nonemptiness gives
`1<=a,b<=p-1`.  On the full affine plane the number of unequal pairs is

\[
                   F(a,b)=p(a+b)-2ab.                      \tag{2}
\]

The minimum of this bilinear expression on the displayed rectangle occurs
at a corner and is `2(p-1)`.  Removing the origin subtracts
`epsilon=f(0)+g(0) mod 2`.  But

\[
 F(a,b)\equiv a+b\equiv\epsilon\pmod2.
\]

Hence `F-epsilon>=2(p-1)` also when `epsilon=1`.  Dividing by the antipodal
action proves (1).  The same count gives

\[
 |S_A|={pa-f(0)\over2},                                  \tag{3}
\]

which also proves the size assertion and its radial equality case.

## 2. Exact minimum distance and the open interval in the spectrum

Use the already-proved normal form

\[
 \operatorname {Row}(D)=
 (\langle\mathbf1\rangle\otimes\mathbf F_2^\Delta)
 \mathbin{\dot+}\bigoplus_A(B_A\otimes B_A).              \tag{4}
\]

Write a word in grouped form

\[
 W=\sum_{A,j}c_{A,j}\otimes b_{A,j},\qquad
 c_{A,j}\in\langle\mathbf1\rangle+B_A.                  \tag{5}
\]

Let `S_A` be the midpoint support on which the coefficient vector
`(c_{A,j})_j` is nonzero, let `k` be the number of nonempty `S_A`, let
`R` be their union, and for `x in R` put

\[
 b_x=|\{A:x\in S_A\}|.
\]

The affine-block transform of the row `W_x` is active in exactly `b_x`
direction groups.  The all-prime grouped uncertainty theorem therefore
gives

\[
             \operatorname {wt}(W_x)\ge d-b_x,            \tag{6}
\]

and a nonzero row has weight at least one when `b_x=d`.

Assume first that no `S_A` is all of `Delta` and `2<=k<d`.  By (1),

\[
 D_{act}:=\sum_{A<B}|S_A\mathbin\triangle S_B|
          =\sum_x b_x(k-b_x)\ge hk(k-1).                  \tag{7}
\]

Put

\[
 G:=\sum_x(k-b_x)=\sum_A|R\setminus S_A|.
\]

Since `b_x<=k-1` whenever `k-b_x` is nonzero,
`D_act<=(k-1)G`; hence `G>=hk`.  Equations (3) and (6) now give

\[
\begin{aligned}
 \operatorname {wt}(W)
 &\ge\sum_{x\in R}(d-b_x)\\
 &=(d-k)|R|+G\\
 &\ge(d-k)h+hk=dh=N.                                    \tag{8}
\end{aligned}
\]

If `k=d` and no support is full, the same calculation gives

\[
 hd(d-1)\le\sum_xb_x(d-b_x)
            \le(d-1)\sum_{b_x<d}(d-b_x),                 \tag{9}
\]

so (6), together with the positive contribution of any `b_x=d` rows,
again yields

\[
 \operatorname {wt}(W)\ge
 \sum_{0<b_x<d}(d-b_x)+|\{x:b_x=d\}|\ge N+|\{x:b_x=d\}|. \tag{9a}
\]

If some `S_A=Delta`, every midpoint row has a
nonzero transformed component, so directly `wt(W)>=N`.

It remains to take `k=1`.  A nonradial `S_A` has size at least `p`, and
(6) gives `wt(W)>=p^2>N`.  If `S_A=l_A`, every nonzero coefficient word
in (5) supported there is exactly `l_A`.  If `s` right blocks occur, their
supports are disjoint and

\[
                         \operatorname {wt}(W)=sph.        \tag{10}
\]

Weight below `N=(p+1)h` forces `s=1`.  Thus `W` is exactly

\[
                         l_A\otimes b_{A,j}.               \tag{11}
\]

Conversely every word (11) has weight `ph`.  We have proved

\[
 \boxed{d(\operatorname {Row}(D))=ph},                    \tag{12}
\]

the `dh=N` words of form (11) are precisely the minimum words, and

\[
 \boxed{\operatorname {Row}(D)\text{ has no weight }w
                 \text{ with }ph<w<N}.                   \tag{13}
\]

The known vertical fibres and scalar graphs have weight exactly `N`, so
(13) is sharp and does not classify the weight-`N` layer.

## 3. The full structured puncture consequence

Let `q<=h+1` be the number of nonzero hard centers; a zero-center row needs
no localized half.  Let `U` be the actual support of a ternary sum of those
`q` localized Mobius halves, after all cancellations.  The midpoint-conic
theorem gives, for
every minimum rectangle `X_(K,beta)`,

\[
 |U\cap X_{K,\beta}|\le2q\le2(h+1)=p+1<ph              \tag{14}
\]

throughout the branch range `p>=31`.  Thus `U` contains no minimum word.
If additionally `|U|<N`, (13)--(14) show that `U` contains no support of
any nonzero word of `Row(D)`.  By the exact puncture duality criterion,

\[
                  \boxed{D_U\text{ is onto over }\mathbf F_2}. \tag{15}
\]

This covers every number of nonzero hard centers:

* If `q<=h`, even the uncancelled support obeys
  `|U|<=q(p-1)<=h(p-1)=ph-h<ph`, so (15) follows already from
  minimum distance.
* If `q=h+1`, put `|U|=N-2kappa`.  Any cancellation gives `|U|<N` and
  (15).  More importantly, on the balanced branch-C ray
  `|H|<=|H|_max=N-1`; hence the necessary physical-size condition
  `|U|<=|H|` itself implies (15), at every `t` in the ray.

For comparison, minimum distance alone in the all-active case would require

\[
 \kappa\ge r+1,
 \qquad t\le t_{max}-r=4r^2-3r-5                       \tag{16}
\]

when the bare size lower bound on `kappa` is used.  The gap theorem removes
this artificial upper `t` cutoff.  At `p=31`, it upgrades the distance-only
range `t<=170` to the entire branch ray `68<=t<=177`, conditional only on
the same central symmetric reduction and the necessary size condition.  In
particular the punctured parity map is harmless on the proved opposite-row
centrality band `68<=t<=164`.

## 4. Exact Hamming ledger and the boundary of the theorem

For `p=4r+3`, a ternary sum of `q` localized halves has

\[
 |U|=q(p-1)-2\kappa,\qquad |H|=4p+2t+1.                  \tag{17}
\]

The bare Hamming capacity requires

\[
 \kappa\ge\kappa_{size}(q,t):=
 \max\{0,qh-2p-t\}.                                    \tag{18}
\]

Each half has one zero-`Phi` orbit and `p-2` nonzero-`Phi` orbits.
Cancellations remove two occurrences of one intrinsic type.  Therefore

\[
 u_0\equiv u_{np}\equiv q\pmod2,
 \qquad |a(T_U)|\equiv |H|+u_{np}\equiv1+q\pmod2.       \tag{19}
\]

For this particular one-localized-half-per-nonzero-center ansatz, the
Hamming numerator has parity `q`.  Thus odd `q` cannot complete within this
ansatz.  This is not a target obstruction: another antisymmetric preimage or
an antisymmetric kernel move is not ruled out.

For even `q`, the parity is automatic, but (15) still supplies only a
solution of the halved mod-two equations.  A physical common graph must also
satisfy

\[
 2\sum b_O=|H|-|U|-|a(T_U)|,                              \tag{20}
\]

every exact direction weight

\[
 n_L={P_L-u_L-f_L\over2},\qquad
 0\le n_L\le dh^2-u_L,                                  \tag{21}
\]

and the full divided integral target equation.  Surjectivity in (15) does
not prescribe the weight in (20), does not place a parity solution in the
slices (21), and does not prove nonnegative integral realizability.  The
remaining symmetric Boolean fibre, residual (ii), E1, and `L=1/2` therefore
remain open.

## 5. Executable transcription

The formula ledger is in
`src/e1_gmin_m4_symmetric_halved_row_code_gap.py`; its focused tests are in
`tests/test_symmetric_halved_row_code_gap.py`.  They check the exact symbolic
inequalities and branch-C thresholds.  They do not enumerate primes or row
codes.
