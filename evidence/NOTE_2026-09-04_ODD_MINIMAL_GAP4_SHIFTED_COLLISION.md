# Proposition 15.767: shifted level-five collision for odd minimal gaps

**Status:** proved infinite-band enlargement of the minimal-four-gap bridge.
It does not close the bridge globally, residual (ii), E1, `L=1/2`, or the
original MathOverflow problem. No graph, prime, orbit, solver, or eigenshell
census is used.

## 1. Shift the affine slack from three to five

Let `C` be the Paley conference matrix of order `p^2+1`, and suppose the odd
no-bridge graph left by Proposition 15.764 has

\[
 |H|=h=5p+2t,
 \qquad T_H^\epsilon(y)\ge5
 \quad(y\in E_\epsilon,\ \epsilon=\pm1).                 \tag{1}
\]

In an affine direction `d` of quadratic type `epsilon`, define on
`J(p,(p+1)/2)`

\[
 B_d={T_H^\epsilon-5\over2},\qquad c_d=2p\,\mathbb E B_d. \tag{2}
\]

This is a nonnegative integer-valued quadratic. Let `I` be the infinity
degree, `P_d` the number of finite H-edges parallel to `d`, and

\[
 T=\sum_{uv\in H,\ u,v\ne\infty} C_{uv}.
\]

The exact middle-slice mean used in Proposition 15.632 gives

\[
 c_d=I+(p+1)P_d-\epsilon T-5p.                            \tag{3}
\]

There are `m=(p+1)/2` directions of each type. Summing (3), or equivalently
subtracting `2p` in each direction from Proposition 15.632's level-three
budget, gives

\[
 \sum_{d:\epsilon_d=\epsilon}c_d=(p+1)t.                 \tag{4}
\]

Put `M=p+1=2m`. The `c_d` are nonnegative even integers with one residue in
each type, so uniquely

\[
 c_d=2u_\epsilon+M k_d,qquad
 0\le u_\epsilon<m,quad k_d\ge0,quad
 \sum_d k_d=t-u_\epsilon.                                \tag{5}
\]

In particular `u_epsilon<=t`.

## 2. Every positive shifted slack has mass about p/2

The edge-product identity makes `B_d mod 2` an affine parity on the middle
slice. If it is nonconstant, the central Krawtchouk bound of Proposition
15.750 puts at least `(p-1)/(2p)` of the slice in its odd class, hence
`c_d>=p-1`. If it is constantly odd, `c_d>=2p`. If it is constantly even
and nonzero, write `B_d=2Q_d`; Proposition 15.681 gives

\[
 c_d=4p\,\mathbb E Q_d\ge
 \beta_p:=
 \begin{cases}
 (p+1)/2,&p\equiv3\pmod4,\\
 (p-1)/2,&p\equiv1\pmod4.
 \end{cases}                                               \tag{6}
\]

Both values of `beta_p` are even. Thus universally

\[
 c_d=0\quad\hbox{or}\quad c_d\ge\beta_p.                  \tag{7}
\]

Let `gamma_p=beta_p/2`. While `t<m`, equation (5) has a zero quotient in
each type. Consequently

\[
 u_\epsilon=0\quad\hbox{or}\quad u_\epsilon\ge\gamma_p. \tag{8}
\]

## 3. The low band already contradicts the degree sum

If `t<gamma_p`, (8) forces `u_+=u_-=0` in every affine chart. Adding the two
residues from (3) and dividing by two modulo `m` gives

\[
 I+5\equiv u_++u_-\equiv0\pmod m.                         \tag{9}
\]

Signed PSL transport can send any graph vertex to infinity, so (9) says

\[
 \deg_H(v)\equiv-5\pmod m\qquad\hbox{for every }v.        \tag{10}
\]

For `p>=11`, every degree is therefore at least `m-5=(p-9)/2`. But

\[
 (p^2+1)(m-5)>2(5p+2t)
\]

throughout `1<=t<gamma_p`. This excludes, in particular, `t=1,2` at
`p=11`, `t=1,2` at `p=13`, `t<=3` at `p=17`, and `t<=4` at `p=19`.

## 4. A zero shifted direction is rigid

For the extended band below, `p>=13` and `H` has an isolated vertex. Send
it to infinity, so `I=0`. If one type has `u_epsilon=0`, (5), with `t<m`,
gives a direction with `c_d=0`. Nonnegativity makes `B_d` identically zero,
so the signed H-score is identically five in that direction.

Let `K_st` be the signed sum of nonparallel edges between fibres `s,t`.
Johnson-slice constancy makes every `K_st` one integer `kappa`, and

\[
 P_d-{p-1\over2}\kappa=5.                                \tag{11}
\]

If `kappa>=1`, signed block capacity forces at least

\[
 P_d+{p\choose2}\kappa
 \ge5+{p^2-1\over2}
\]

edges, above the theorem range. If `kappa<=-1`, (11) gives
`P_d<=5-(p-1)/2<0` for `p>=13`. Hence

\[
 \kappa=0,\qquad P_d=5,\qquad \epsilon T=5.              \tag{12}
\]

The two direction types cannot both have zero residue, because (12) would
give simultaneously `T=5` and `T=-5`.

## 5. The isolated-chart type collision

At `I=0`, adding the two residues in (3) gives

\[
 u_++u_-\equiv5\pmod m.                                  \tag{13}
\]

Together with (8) and the exclusion of `(u_+,u_-)=(0,0)`, this has no
solution in the following ranges:

- `p=13,17,19` and `1<=t<=4`;
- every prime `p>=23` and `1<=t<=floor((p+10)/4)`.

For the second line, `gamma_p>5`. If exactly one residue is zero, (13)
would force the other to equal five, contrary to (8). If both are positive,
their sum is at least `2gamma_p>5`; since `t<m`, its first possible value in
(13) is `m+5`. But

\[
 u_++u_-\le2t<m+5
\]

in the displayed range. For `p=13,17,19`, direct interval arithmetic in
(8), not a graph census, gives the same contradiction through `t=4`.

Combining this with Section 3 proves

\[
 \boxed{
 \begin{array}{ll}
 p=11:&1\le t\le2,\\
 p\in\{13,17,19\}:&1\le t\le4,\\
 p\ge23:&1\le t\le\lfloor(p+10)/4\rfloor
 \end{array}
 \quad\Longrightarrow\quad\text{the odd no-bridge branch is empty}.}
                                                                  \tag{14}
\]

The first rows not excluded by this proposition are therefore

\[
 \boxed{
 p=11:t=3;\qquad p=13,17,19:t=5;\qquad
 p\ge23:t=\left\lceil{p+11\over4}\right\rceil.}          \tag{15}
\]

For `p>=13`, the executable certificate records the surviving residue pairs
at the first row. They certify only that this method has reached its exact
barrier; they are not common Paley graphs and not counterexamples.

## 6. Scope and replay

The argument uses only common Paley graph structure, signed PSL chart
transport, exact directional means, boundary parity, Johnson constancy, and
the proved integral quadratic mass floor. It does not use a fixed-prime
search or assume an affine classification of the full defect shell.

Exact replay:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q tests/test_prop15767.py
PYTHONPATH=src /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_prop15767.py
```

Proposition 15.767 enlarges the implication bridge of Proposition 15.764.
It does not address the dangerous off-eigenshell maximizer in an already
given residual-(ii) pair, and it changes no global closure predicate.
