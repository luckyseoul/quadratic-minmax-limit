# Proposition 15.764: the exact minimal-gap-four shell bridge

**Status:** proved parity equivalence and a genuine bridge for odd
`|H|<=5p`. The first unclosed failure ranges are even `|H|>=4p+2` and odd
`|H|>=5p+2`. This does not close residual (ii), E1, `L=1/2`, or the original
MathOverflow problem, and no global predicate changes.

## 1. Exact setup

Let `C` be the Paley conference matrix of order `n=p^2+1` and put

\[
 \Phi={p(p^2+1)\over2}.
\]

For \(\epsilon\in\{+1,-1\}\), define

\[
 E_\epsilon=\{y\in\{\pm1\}^n:Cy=\epsilon py\},\qquad
 T_F^\epsilon(y)=
 \epsilon\sum_{\{u,v\}\in F}C_{uv}y_uy_v,
\]

and \(m_\epsilon(F)=\min_{y\in E_\epsilon}T_F^\epsilon(y)\). Thus

\[
 m_+(F)=s_+(F),\qquad m_-(F)=-s_-(F).                    \tag{1}
\]

Suppose `H` is an inclusion-minimal deeper set at the first possible
four-unit gap, in the exact all-deletions form

\[
 \Phi(C\mathbin\triangle H)=\Phi-4,
 \qquad
 \Phi(C\mathbin\triangle(H\setminus\{e\}))=\Phi-2
 \quad(e\in H).                                           \tag{2}
\]

The argument below only needs the corresponding upper bounds. For every
\(y\in E_\epsilon\),

\[
 \epsilon Q_{C\triangle F}(y)=\Phi-2T_F^\epsilon(y).
\]

Consequently (2) gives

\[
 T_H^\epsilon(y)\ge2,
 \qquad
 T_{H\setminus\{e\}}^\epsilon(y)\ge1.                    \tag{3}
\]

Every term in a `T`-score is a sign, so

\[
 T_F^\epsilon(y)\equiv |F|\pmod2.                         \tag{4}
\]

## 2. The parity equivalence

Fix a row and abbreviate

\[
 b_e=\epsilon C_{uv}y_uy_v\in\{+1,-1\},\qquad
 T_{H\setminus\{e\}}=T_H-b_e.                            \tag{5}
\]

Suppose first that `|H|` is odd. By (3)--(4), every H-score is odd and at
least three, while every deletion score is even and at least two. If a
deletion score is two, (5) makes the H-score either one or three; (3) rules
out one. Conversely, if an H-score is three, at least one of its edge signs
is positive, and deleting such an edge gives score two. Therefore

\[
 \boxed{
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})=2
 \iff
 \exists\epsilon,y\in E_\epsilon:\ T_H^\epsilon(y)=3.}    \tag{6}
\]

By (1), the left side is precisely a deletion with `s_+=2` or `s_-=-2`,
the official even-cardinality residual-(ii) entry level. In fact it supplies
the rest of that entry's hypotheses. For every row with deletion score two,
(5) and the H-score floor force `b_e=+1`; hence the distinguished edge freezes
positive on the entire critical level. Both signed phases of the even deletion
have minimum at least two, so it is deep two-sided.

If the critical phase is minus, normalize it to plus by replacing `C` with
`-C`. This stays within the Paley switching class. In the standard
`F_(p^2) union {infinity}` coordinates, multiplication by a nonsquare negates
all finite--finite Paley entries, and switching at infinity negates all
infinity--finite entries. Thus

\[
 -C=DP^tCPD,                                                \tag{6a}
\]

for a permutation `P` and diagonal sign matrix `D`. The edge set is merely
permuted and the signed minus-shell features become plus-shell features.

Finally, the H-score floor three and the frame mean give `|H|>=3p`. Equality
would make H bi-tight of level three, excluded by Proposition 15.720. Since
both quantities are odd,

\[
 |H|\ge3p+2,\qquad |H|-1\ge3p+1.                           \tag{6b}
\]

Thus (6), after phase normalization, gives an even deletion with `s_+=2`,
deep two-sidedness, `f_e=+1` on all of `U_2`, and `k>=3p+1`: the complete
official residual-(ii) entry, not only the numerical level.

If `|H|` is even, every deletion has odd score. It therefore cannot have
signed shell minimum two. The audited parity statement is instead

\[
 \boxed{
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})\le2
 \iff
 \exists\epsilon,y\in E_\epsilon:\ T_H^\epsilon(y)=2,}    \tag{7}
\]

and the deletion score in (7) is exactly one. This is the Type-I parity,
not residual (ii). Every level-one deletion row similarly has `b_e=+1`, and
the same switch/permutation normalization puts it in the plus Type-I
convention. Here the level-two frame mean gives `|H|>=2p`; equality is the
forbidden bi-tight level two, so `|H|>=2p+2` and `|G|>=2p+1`. Thus an argument
that sends every deletion directly into
the even-`k`, `s_+=2` ledger without first proving `|H|` odd has a real gap.

## 3. Frame averaging forces a size dichotomy

Proposition 15.42's signed frame identity says, in either phase,

\[
 \mathbb E_{y\in E_\epsilon}
   [\epsilon C_{uv}y_uy_v]={1\over p}.
\]

Hence

\[
 \mathbb E_{E_\epsilon}T_H^\epsilon={|H|\over p}.         \tag{8}
\]

If odd `H` fails (6), its odd shell scores are at least five in both phases.
Equation (8) forces

\[
 |H|\ge5p.                                                 \tag{9}
\]

At equality, every score is exactly five in both phases: `H` is bi-tight of
level five.

If even `H` fails (7), its even shell scores are at least four in both
phases, so `|H|>=4p`. Equality is bi-tight of level four and is excluded for
every prime `p>=5` by Proposition 15.720. Since `|H|` is even,

\[
 \boxed{\text{even failure requires }|H|\ge4p+2.}          \tag{10}
\]

## 4. The level-five equality cannot be deeper

The generic degree congruence in Proposition 15.720 applies at any bi-tight
level `s`. For `s=5`, with \(M=(p^2-1)/2\), it says

\[
 d_i+d_j\equiv10p\pmod M\quad(i\ne j),                    \tag{11}
\]

and all degrees have one common residue modulo `M`.

For `p>=11`, `M>5p`, so all degrees are equal. The handshake identity would
give

\[
 d={10p\over p^2+1}\in(0,1),
\]

impossible for an integer degree. At `p=7`, (11) gives the common residues
11 or 23 modulo 24. Every degree is therefore at least 11, contradicting
\(\sum_i d_i=70\) on 50 vertices.

At `p=5`, `M=12`, `n=26`, and the total degree is 50. Equation (11) gives
residue 1 or 7; the latter is already too large. Degrees are therefore in
`{1,13,25}`, and distributing the excess 24 above the baseline degree-one
sequence leaves exactly two profiles:

- one degree-25 vertex and twenty-five degree-1 vertices: a full star;
- two degree-13 vertices and twenty-four degree-1 vertices: a balanced
  double star.

A full star is a vertex switch: if `D` changes the sign at its centre, then
`C triangle H=DCD`, so its norm is `Phi`, not `Phi-4`.

It remains to exclude the balanced double star. Let `a,b` be its centres.
Counting incidences forces the edge `ab`, no leaf-leaf edge, and twelve leaves
at each centre. Put

\[
 c=C_{ab},\quad \alpha_i=C_{ai},\quad \beta_i=C_{bi},\quad
 r_i=\begin{cases}+1&i\text{ is attached to }a,\\-1&i\text{ is attached to }b.
 \end{cases}
\]

For the signed support matrix `B=C odot h`, Proposition 15.720's exact
`scheme+cross` decomposition is

\[
 B=D_gC+CD_g+X,\qquad CX+XC=0,
\]

where \(g_i=(d_i-1)/24\). Thus `g_a=g_b=1/2`, every leaf coordinate is
zero, and

\[
 X_{ai}={\alpha_i r_i\over2},\qquad
 X_{bi}={-\beta_i r_i\over2},                              \tag{12}
\]

with every other entry of `X` zero. The leaf-leaf entries of `CX+XC=0`
give

\[
 (\alpha_i\alpha_j-\beta_i\beta_j)(r_i+r_j)=0.            \tag{13}
\]

Therefore \(t_i=\alpha_i\beta_i\) is constant on each twelve-leaf part.
Conference-row orthogonality gives

\[
 \sum_i t_i=(C^2)_{ab}=0,
\]

so the two constants are opposite. There is one sign
\(\tau\in\{\pm1\}\) such that

\[
 r_i=\tau\alpha_i\beta_i.                                 \tag{14}
\]

The `(a,j)` entry of `CX+XC=0` is

\[
 \sum_i\alpha_i r_iC_{ij}=c\beta_jr_j.                   \tag{15}
\]

Using (14) and `(C^2)_{bj}=0`, the left side of (15) is
`-tau*c*alpha_j`; its right side is `+tau*c*alpha_j`. Both are nonzero.
This contradiction excludes the balanced double star.

Thus a bi-tight level-five graph cannot also satisfy the deeper condition in
(2). Combining this with (6) and (9) proves

\[
 \boxed{|H|\text{ odd and }|H|\le5p
 \Longrightarrow
 \exists e,\epsilon:\ m_\epsilon(H\setminus\{e\})=2.}    \tag{16}
\]

The first remaining odd failure size is therefore `5p+2`.

## 5. Exact method barrier for the scalar identities

The size `5p` boundary is not removable from frame means alone. Here is a
fully specified finite max-of-affine model at `p=5`, `Phi=65`, on a formal
25-edge set `H`.

Take as each signed eigenshell all vectors

\[
 b\in\{\pm1\}^{25},\qquad \sum_i b_i=5.
\]

Uniformly, every coordinate has mean `1/5`. Every H-score is five, and a
one-edge deletion has shell score four or six, hence shell minimum four.

For each \(D\subset H\) of size `r=1,2,3`, add an off-shell affine row
`a_D` which is `+1` on `D` and `-1` elsewhere. Give it base defect

\[
 \delta_D=54-4r.
\]

Its base signed score is `65-delta_D`, and at `H` its switched score is

\[
 65-\delta_D-2\sum_{i\in H}(a_D)_i=61=\Phi-4.             \tag{17}
\]

At its own complement `H minus D`, the score is respectively 63, 65, or 67
for `r=1,2,3`. For a complement of size at least four, a fixed-sum shell row
can be chosen whose retained-edge sum is at most one, giving score at least
63. Thus `H` is inclusion-minimal at the four-gap threshold. Each one-edge
deletion has norm exactly 63: its own `r=1` row reaches 63, and the norm can
increase by at most two from the value 61 at `H`. Add the raw negatives of
all rows to make the model two-sided.

This formal model satisfies the flip-affine law, parity, frame first moments,
minimal four-gap condition, and exact all-deletions two-gap condition, yet no
deletion has signed shell minimum two. It is **not** asserted to be induced by
a Boolean Paley quadratic form. Its role is precise: Propositions 15.42 and
15.755's scalar identities cannot prove the global bridge without additional
common-graph structure.

## 6. Exact remaining bridge

To place an arbitrary minimal four-gap path into the existing official units,
one must still prove at least one of the following structural assertions:

1. `H` is odd and `|H|<=5p`; or
2. directly, some signed H-eigenshell row has level three; or
3. separately exclude the two failure regimes
   `even |H|>=4p+2` and `odd |H|>=5p+2`.

Proposition 15.763's affine-alias alternative can produce the row in item 2
for a particular deletion, but it does not classify every off-shell spike or
align the coordinate systems of different deletions. No such global theorem
is currently proved.

Exact replay:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q tests/test_prop15764.py
PYTHONPATH=src /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_prop15764.py
```
