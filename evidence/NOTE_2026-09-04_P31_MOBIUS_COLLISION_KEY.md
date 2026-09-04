# A centre-free collision key for p31 localized-Mobius halves

Date: 2026-09-04

**Status:** proved a projective necessary-and-sufficient criterion for two
chosen nonorigin Mobius skeleton edges to be centred onto one inversion
orbit.  An exact class join excludes the one frozen sixteen-half hard-fixed
witness in
`/tmp/resii_p31_hard_fixed_fixed_pairs_signature_seed920260927.json` without
enumerating centres.  This is a physical prefilter for that frozen witness,
not a closure of the hard-fixed branch, branch C, or residual (ii).

## 1. Collision-key theorem

Let `L,M` be independent linear functionals over `F_p`, let the Mobius
centre be `c != 0`, and let `t != -1` be its parameter.  In `(L,M)`
coordinates the edge at unit centre is

\[
 E_{L,M,t}(1)=\left\{
 T_{L,M}\left(1,{t\over t+1}\right),
 T_{L,M}(t,t)
 \right\},                                             \tag{1}
\]

and the edge at centre `c` is `c E_{L,M,t}(1)`.

For `t != 0,-1`, define three projective functional directions

\[
 U_t=[tL-(t+1)M],\qquad V=[L-M],\qquad
 D_t=[t^2L+(1-t^2)M].                                  \tag{2}
\]

They annihilate respectively the first endpoint, the second endpoint, and
the edge displacement in (1).  If `Delta=det(L,M)`, direct expansion gives

\[
 \det(U_t,V)=\Delta,\qquad
 \det(V,D_t)=\Delta,\qquad
 \det(U_t,D_t)=t(t+1)\Delta.                            \tag{3}
\]

Thus `U_t,V,D_t` are pairwise distinct.  The **collision key** is

\[
                 K(L,M,t)=(D_t,\{U_t,V\}).              \tag{4}
\]

It is independent of the centre.

**Theorem.** Two chosen nonorigin skeleton edges `E_1(1),E_2(1)` can be
centred onto the same inversion orbit if and only if their keys (4) agree.
More precisely, when the keys agree there is a unique `lambda != 0` with

\[
                         E_2(1)=\lambda E_1(1).          \tag{5}
\]

For nonzero centres `c_1,c_2`, their physical edges occupy the same
inversion orbit exactly when

\[
                  c_1=\lambda c_2\quad\hbox{or}\quad
                  c_1=-\lambda c_2.                     \tag{6}
\]

The second sign gives negative physical edges and therefore opposite
normalized antisymmetric coefficients; it is exactly the cancellation
choice.

Necessity follows because one inversion orbit fixes the unordered endpoint
rays and the displacement ray.  For sufficiency, match the endpoint rays and
write `x_2=a x_1`, `y_2=b y_1`.  Equality of displacement directions says
`a x_1-b y_1` is proportional to `x_1-y_1`.  The two endpoint rays are
independent by (3), so coefficient comparison forces `a=b=lambda`, proving
(5).  The centre relation (6) follows immediately.

The exceptional parameter `t=0` is the origin edge.  Its only nonzero
endpoint and its spatial direction both have class `[M]`.  Hence two origin
edges can meet only when their auxiliary directions agree; an auxiliary SDR
excludes such a cancellation.  The parameter `t=-1` is absent.

## 2. Boundary and parallel corrections

Removing one nonorigin inversion pair with key `(D,{U,V})` changes the two
coupled ledgers by

\[
 \hbox{boundary signature}=e_U\mathbin\oplus e_V,
 \qquad \hbox{parallel loss}=2e_D.                       \tag{7}
\]

A fixed antipodal edge in direction `F` contributes boundary `e_F` and
parallel gain `e_F`.  Consequently the hard-fixed `j=0` equations have the
form

\[
 R=P+2e_D-e_F,\qquad G\mathbin\oplus e_F=e_U\mathbin\oplus e_V. \tag{8}
\]

In particular, the two-bit word `h=G xor e_F` fixes the unordered endpoint
pair `{U,V}`.  It does **not** in general force `D=F`; the physical key must
supply `D`.

## 3. One-candidate-per-half hash join

The second endpoint class of every edge in a labelled half is fixed:

\[
                              V=[L-M].                   \tag{9}
\]

If the required correction support is `{A,B}`, reject the half unless (9)
is one of `A,B`.  Let `H` be the other class.  The equation
`[t(L-M)-M]=H` then has the unique solution

\[
                 t={\det(M,H)\over\det(L-M,H)}.          \tag{10}
\]

Reject `t=0,-1`; otherwise (2) gives one and only one direction bucket `D`
for this half.  A pair or `2:1` triple collision with this boundary support
is possible only if at least two selected halves enter the same bucket.
Conversely, any two entries in one bucket can be aligned with the centre
relation (6).  This is a hash join over projective keys; no centre scan is
needed.  It is only a local collision criterion: a full sixteen-centre
assignment must still avoid or account for all other simultaneous overlaps.

There is also an exact inverse that avoids scanning the 29 parameters or the
available auxiliary scales.  Fix a target `L`, choose which prescribed
endpoint is `K=V`, call the other one `H=U`, and let the requested spatial
direction be `D`.  Then `t` is projectively invariant and

\[
 t={\det(D,H)\det(K,L)\over
        \det(K,D)\det(L,H)},\qquad
 L-M=sK,\quad
 s={\det(L,H)\over\det(K,H)(t+1)}.                       \tag{11}
\]

The scalar `s` depends covariantly on the chosen representative of `K`, while
the functional `sK=L-M` does not.  There are only two choices for `K`, so a
fixed `(D,{U,V})` shard contains at
most two eligible labelled options per target, at most 32 over all sixteen
targets.  The implementation validates the reconstructed option by feeding
it back through (2).

An exhaustive p31 forward/inverse audit covers all
`14,430*29=418,470` hard-fixed labelled key memberships.  The inverse catalog
equals the direct formula catalog exactly: both set differences are empty.
The ordered membership catalog has SHA-256
`2775cd4ad86a1844834dba1981a8396a92fb2f6da9b840a739e9099abef32897`.
This audit uses no centre variables.

For an unpinned correction word there are

\[
 (p+1)\binom p2=14,880\quad (p=31)                      \tag{12}
\]

possible `(D,{U,V})` keys, because `U,V` are distinct and both differ from
`D`.  Sharding the simultaneous solver by `{U,V}` reduces each labelled
half from 29 collision memberships to at most one candidate from (10).

## 4. Frozen hard-fixed witness replay

The frozen choices `(target, auxiliary direction, relative scale)` are

```text
(0,22,12) (1,27,27) (2,15,11) (3,28,7)
(7,16,9)  (9,4,5)   (10,7,20) (15,26,16)
(16,24,28) (21,29,15) (22,9,19) (24,10,15)
(28,21,21) (29,31,27) (30,2,1) (31,3,25)
```

The exact replay recovers the raw profile

```text
15 15 14 14 15 16 16 14 16 14 14 16 16 16 16 14
14 16 16 16 16 14 14 16 14 16 15 15 14 14 15 14
```

and aggregate selector signature `0x00800005`.  With the fixed direction
`F=0`, the correction is `0x00800004`, supported on `{2,23}`, while the
parallel ledger requires cancellation direction `D=0`.

Applying (9)--(10), only one selected half contains the prescribed endpoint
pair:

```text
(target, auxiliary, scale) = (29,31,27),  t=12,
K = (D,{U,V}) = (8,{2,23}).
```

There is no second owner, and its direction is `8`, not the required `0`.
The frozen family is therefore physically impossible.

As an independent cross-check, joining all `16*29=464` labelled keys gives
exactly the five shared classes

```text
(4,{3,9}) (9,{16,31}) (20,{28,29})
(21,{27,31}) (29,{7,16}).
```

Each shared class gives `2(p-1)=60` nonzero centre pairs on the same
inversion orbit, exactly `p-1=30` with cancelling orientation.  The key join
therefore derives 300 shared-orbit incidences and 150 cancelling incidences,
agreeing with the prior 108,000-pair physical replay; none has direction
zero.  The ordered labelled-key catalog has SHA-256
`f3512567d9655427212ec0383b6d1e0f2f15b364031fa348980a44cb5724d80b`.

## 5. Reproduction

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_p31_mobius_collision_key.py \
  tests/test_residual_branch_c_hard_fixed_center_obstruction.py
```

The implementation is
`src/e1_gmin_m4_p31_mobius_collision_key.py`.  It verifies (2) directly
against every physical nonorigin edge in the frozen family, tests the exact
necessary-and-sufficient centre relation, replays the five-class join, and
checks the one-candidate `{2,23}` obstruction.

The optional full 418,470-membership inverse audit is

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -c \
  'from e1_gmin_m4_p31_mobius_collision_key import full_p31_inverse_catalog_replay as f; print(f())'
```
