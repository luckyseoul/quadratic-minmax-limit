# Minimum-degree normalization is exact, but aggregate chart equations do not bound the shell

Date: 2026-09-04

Scope: canonical residual (ii), with `G` an even flip graph, `e` the
distinguished edge, and `H=G union {e}`. This note derives the exact general
chart system and tests the proposed minimum-degree normalization. It neither
uses a finite-prime scan nor closes residual (ii).

## 1. Exact residual inequalities

Write `f_a(y)=C_a y_a` for the signed edge feature and
`S_F=sum_(a in F) f_a`. In the deep residual-(ii) bad case,

```
S_G >= 2 on Max+,
S_G=2 implies f_e=+1,
S_H=S_G+f_e <= -3 on Max-.
```

Equivalently, with `T_H^eps=eps*S_H` on the signed eigenshell `E_eps`, the
gap-four box is

```
3 <= T_H^eps <= Phi-2,       eps=+1,-1,
E_eps[T_H^eps]=|H|/p.
```

The upper inequality is essential once `|H|` is no longer automatically
small compared with `Phi`; retaining only the usual nonnegative slack loses
the absolute-value condition.

## 2. Minimum-degree and outside-boundary charts

Signed PSL transport permutes the relative flip mask, its degree sequence,
its odd-degree boundary, and both signed shell inequalities. Therefore a
minimum-degree vertex may be sent to infinity and

```
I <= floor(2|H|/(p^2+1)).
```

In particular `2|H|<p^2+1` forces `I=0`. This includes every existing
`O(p)` shell once `p` is beyond its explicitly handled small thresholds, so
it recovers the isolated chart already used by Propositions 15.734--15.752;
it does not strengthen those floor arguments.

If `D` is the odd-degree boundary, `|D|=s<p^2+1`, then every degree outside
`D` is even and its total is at most `2|H|-s`. Hence some outside vertex has

```
deg(v) <= 2 floor((2|H|-s)/(2(p^2+1-s))).
```

This forces an isolated outside point when `2|H|+s<2(p^2+1)`. If `D` is all
vertices there is no outside point, which is a real obstruction rather than
a chart choice.

The signed PSL action is two-transitive but the minimum vertex is already one
prescribed point. Unless it lies on `e`, one cannot also assume that the two
endpoints of `e` are `(infinity,0)` without an additional stabilizer argument.

## 3. General affine chart equations with `I` retained

For a projective `F_p` direction `L`, let `eps_L` be its quadratic type,
`P_L` the number of finite `H`-edges parallel to it, and
`T=sum_(finite a in H) C_a`. The affine middle-slice slack and its scaled
mean are

```
A_L=(T_H^eps_L-3)/2 >= 0,
a_L=2p E[A_L]=I+(p+1)P_L-eps_L*T-3p.
```

They obey

```
I+sum_L P_L=|H|,
sum_(eps_L=tau) a_L=((p+1)/2)(|H|-3p),
a_L-a_M=(p+1)(P_L-P_M)                 (eps_L=eps_M).
```

Let `D` be the boundary, `b_L` the number of its odd finite fibres, and
`c_H=product_(a in H) C_a`. The exact parity phase is

```
(-1)^eta_L = eps_L^(1+1_(infinity in D))
             (-1)^((|H|-3)/2+b_L) c_H,
```

and the exact quadratic parity-majorant gives

```
a_L >= 2 ceil(p M(p,b_L,eta_L)).
```

For finite displacement multiplicities `m_delta` and normalized difference
rows `q_L(a)`, the common-graph equations additionally give

```
sum_(a>0) q_L(a)=eps_L*T-P_L,
sum_(L,a>0) q_L(a)^2
  = p sum_delta m_delta^2 + 2T^2 - 2sum_L P_L^2,
sum_delta m_delta^2=(|H|-I)+2sum_delta binom(m_delta,2).
```

These are identities and necessary inequalities. The pointwise conditions
`A_L(X)>=0` and the existence of a level-two `G` row are strictly stronger.

## 4. Sharp aggregate counterexample at the isolation threshold

For every odd prime `p>=7`, take the antipodal perfect matching

```
H={{infinity,0}} union {{x,-x}:x in F_(p^2)^*/+/-}.
```

Then

```
|H|=(p^2+1)/2,   deg(v)=1 for every v,   D=all vertices,
I=1,             P_L=(p-1)/2 for every L,   T=0.
```

Every finite antipodal displacement occurs once. In every direction the
normalized difference row consists of `(p-1)/2` copies of `-1`. Therefore
the common energy is exactly

```
p(|H|-I)-2 sum_L P_L^2=(p^2-1)/2=|H|-I.
```

The boundary has `b_L=p`, phase zero in every direction, and

```
a_L=(p^2-6p+1)/2.
```

For `p=3 mod 4` the phase floor is zero. For `p=1 mod 4` it is `2p`, which
is at most `a_L` for every relevant prime `p>=13`. Thus this superlinear
family passes degree averaging, both type budgets, every phase floor, the
same-type congruences, and the common energy identity.

It is not a residual graph. If `z_t` denotes the affine fibre sign and `r`
counts split nonzero antipodal fibre pairs, direct evaluation gives

```
T_H^eps(X)=z_0+2r.
```

Its minimum is `-1` for `p=3 mod 4` and `1` for `p=1 mod 4`, both below the
required value three. Deleting `{infinity,0}` gives minimum `G`-score zero,
not two.

This is a precise method barrier: minimum degree plus the aggregated
phase/floor/common-energy system cannot prove `|H|=O(p)` or reduce every high
shell to the known band. Any successful continuation must use the full
pointwise shell inequalities together with the distinguished level-two row,
or a new consequence of inclusion-minimality that is not visible in those
aggregates.

## Replay

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_residual_min_degree_barrier.py
```
