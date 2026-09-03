# Adaptive center coherence for complementary Mobius pairs

Date: 2026-09-03

Status: proved target-pairing theorem for every branch prime
\(p=4r+3\ge31\) and every list of nonzero hard-star centers.  After an
opposite fixed-edge direction is selected, one can choose the magnitude of
that edge, singleton signs, and a perfect matching of the hard targets so
that every matched pair has the required nonsquare complementary parameter.
This does **not** prove that the resulting auxiliary directions are globally
distinct or have the quota-required Paley types.  It does not construct the
full parallel target, the singleton fixed word, a common graph, or residual
(ii).

## 1. Exact scaled pair equations

For one half write \(X=L/j\), \(N=M/j\).  A complementary mate with
parameter \(\nu\) has

\[
 L'=\nu L+(1-\nu)M,\qquad M'=\nu(L-M).
\]

If its prescribed center is \(j'\), put \(\rho=j'/j\),
\(Y=L'/j'\), and \(N'=M'/j'\).  Then

\[
 \rho Y=\nu X+(1-\nu)N,
 \qquad N'={\nu\over\rho}(X-N).                         \tag{1}
\]

Choose a prospective fixed antipodal edge \(\{x,-x\}\), put
\(a=X(x)\), \(b=Y(x)\), and prescribe singleton signs
\(\epsilon,\epsilon'\in\{\pm1\}\):

\[
 N(x)=2\epsilon,\qquad N'(x)=2\epsilon'.                \tag{2}
\]

Solving (1)--(2) gives

\[
 \boxed{\rho={2\epsilon\over b-2\epsilon'},\qquad
 \nu={4\epsilon\epsilon'
 \over(a-2\epsilon)(b-2\epsilon')}}.                   \tag{3}
\]

The two projective auxiliary directions are therefore

\[
 [N]=[\rho Y-\nu X],\qquad [N']=[X-\rho Y].             \tag{4}
\]

All exceptions are visible in (3): the first sign is invalid at
\(a=2\epsilon\), the second at \(b=2\epsilon'\), and
\(\rho\ne0\), \(\nu\ne0\).  Complementarity requires
\(\eta(\nu)=-1\), so \(\nu=1\) is then automatically absent.  In the
application the fixed direction is opposite and every target is hard, hence
\(a,b\ne0\).

Define the color of a valid signed endpoint by

\[
 c_\epsilon(a)=\eta\bigl(\epsilon(a-2\epsilon)\bigr).
                                                               \tag{5}
\]

Equation (3) says exactly

\[
                         \eta(\nu)
 =c_\epsilon(a)c_{\epsilon'}(b).                         \tag{6}
\]

Thus a pair is complementary precisely when its selected endpoint colors
are opposite.

## 2. Flexible and monochrome targets

For \(a\ne\pm2\), both signs are valid and

\[
 c_+(a)c_-(a)=-\eta(a^2-4).
\]

Hence the endpoint is flexible precisely when
\(\eta(a^2-4)=+1\), and monochrome precisely when
\(\eta(a^2-4)=-1\).  At \(a=2\) only \(\epsilon=-1\) is valid; at
\(a=-2\) only \(\epsilon=+1\) is valid.  Both exceptional endpoints are
monochrome of color \(-1\).

Put \(h=(p-1)/2\), \(m=h+1\).  The standard quadratic-character sum
\(\sum_a\eta(a^2-4)=-1\), with the value at \(a=0\) and the two zeros
handled separately, gives among \(a\in\mathbf F_p^*\)

\[
 \#\{\text{flexible}\}=h-1,qquad
 \#\{\text{monochrome}\}=m.                             \tag{7}
\]

More precisely, if \(\eta(2)=+1\), the two monochrome colors each occur
\(r+1\) times.  If \(\eta(2)=-1\), color \(+1\) occurs \(r\) times and
color \(-1\) occurs \(r+2\) times.

## 3. Arbitrary centers admit a target perfect matching

Fix \(0\ne x_0\in\ker F\) and write

\[
 \alpha_i=(L_i/j_i)(x_0)\ne0
\]

for the \(m\) hard targets.  The magnitude \(x=cx_0\),
\(c\in\mathbf F_p^*\), gives \(a_i=c\alpha_i\).  For every fixed \(i\),
this runs once through all nonzero field values.  By (7), summing the number
of monochrome vertices over all magnitudes gives exactly

\[
                              m^2.                       \tag{8}
\]

Since \(p-1=2h\) and

\[
 {m^2\over p-1}={m^2\over2h}< {m\over2}+1,              \tag{9}
\]

some magnitude has at most \(m/2\) monochrome vertices in total.

At that magnitude the feasible graph on hard targets is the complete graph
with only the edges internal to either monochrome color class deleted.
Each class has size at most \(m/2\).  Match the two monochrome classes to
each other as far as possible, match the excess of the larger class to
flexible vertices of the opposite selected color, and pair the remaining
flexible vertices with opposite selected colors.  The remainder is even
because \(m\) is even.  This constructs a perfect matching, with signs, for
which every \(\nu\) in (3) is a nonsquare.

This is a genuine arbitrary-center coherence theorem.  It is only a
perfect matching of the **hard targets**.

## 4. The forced auxiliary map and the exact remaining gate

Normalize every projective direction other than \(F\) to take value one at
\(x_0\), and let \(z_i\) be the affine coordinate of target \(i\).  For a
chosen sign set

\[
 w_i={2\epsilon_i\over\alpha_i},\qquad
 \mu_i={w_i\over w_i-c},\qquad
 H_i(z)=z_i+\mu_i(z-z_i).                                \tag{10}
\]

For a matched pair \((i,k)\), its auxiliary coordinates \(U,V\) are the
unique two-cycle

\[
                         V=H_i(U),\qquad U=H_k(V),        \tag{11}
\]

namely

\[
 \boxed{
 U=z_k+{w_k(z_i-z_k)\over w_i+w_k-c},\qquad
 V=z_i+{w_i(z_k-z_i)\over w_i+w_k-c}.}                   \tag{12}
\]

The excluded denominators are

\[
 c=w_i\ (a_i=2\epsilon_i),\quad
 c=w_k\ (a_k=2\epsilon_k),\quad
 c=w_i+w_k\ (\nu=1).                                    \tag{13}
\]

For a valid complementary pair, \(U,V,z_i,z_k\) are four distinct finite
directions; in particular each pair is locally clean.  Nothing in the
matching proof prevents an auxiliary from one pair equaling an auxiliary
from another pair.

At the top parallel endpoint, choosing the singleton direction \(F\)
opposite would require the union of all auxiliary directions to be exactly
\(m\) distinct directions, of which \(m-2\) are hard and two are opposite.
Equivalently, one must select a perfect matching in the labeled feasible
edge set while also avoiding every cross-edge intersection of the forced
two-point blocks \(\{U,V\}\) and imposing that exact Paley-type weight.
This is a colored induced-matching/paired-SDR problem, not the ordinary
target perfect matching proved in Section 3.  It remains open.

There is a useful exact invariant for that last problem.  Let \(\tau\) be
the fixed-point-free involution pairing the targets and let
\(\sigma:H\to A\) assign to each target its own auxiliary direction.  The
full affine endpoint equations are

\[
 c\alpha_i\bigl(z_{\sigma(\tau i)}-z_i\bigr)
 =2\epsilon_i\bigl(z_{\sigma(\tau i)}-z_{\sigma i}\bigr)
 \quad(i\in H).                                          \tag{14}
\]

Cross-assign the target labels to the auxiliary endpoints by
\(\phi(\sigma(\tau i))=i\), and put

\[
             g_V=\alpha_{\phi(V)}(z_V-z_{\phi(V)}).
\]

For every auxiliary pair \(\{U,V\}\), equation (14) gives

\[
 \boxed{g_U^2=g_V^2,\qquad
 { (U-V)^2\over g_U^2}={c^2\over4}.}                     \tag{15}
\]

Thus a prescribed quota-compatible auxiliary set needs a bijection for
which the cross-assigned \(g^2\) multiset pairs into equal values, with the
same chord ratio in every pair.  This condition is necessary; no theorem
here proves that such a bijection exists.  The focused test also records a
single \(p=31\) formula replay where the target matching from Section 3 has
only nine distinct auxiliary directions among sixteen occurrences.  That
shows why target matching does not itself prove the SDR.  It is not a
search, and it does not say another matching for those data cannot work.

Even making every square fibre even is not sufficient. For a fixed
cross-assignment \(\phi\), the condition \(g_U^2=g_V^2\) alone is equivalent
to \(\prod_{V\in A}(T-g_V^2)\) being a square polynomial. For one prospective
square \(R=c^2/4\), each fibre \(A_q=\{V:g_V^2=q\}\) must additionally have
a perfect matching using only chords with \((U-V)^2=Rq\), followed by the
nonsquare-complement filter. The unfiltered graph has degree at most two, so
this is an exact component-parity test, not a product or discriminant test.

The distinction is visible on the abstract set \(A=\{0,1,2,4\}\) with all
\(g_V^2=1\). Its polynomial is \((T-1)^4\), while the squared chord pairs in
its three perfect matchings are \((1,4)\), \((4,9)\), and \((16,1)\). They
are unequal over every branch prime. Thus no common \(R\) exists. This is a
counterexample to sufficiency of the paired-square invariant, not a branch-C
target counterexample.

## 5. Exact-magnitude warning

If at one exact representative \(x\) all normalized evaluations happen to
equal \(a_i=2\), then only \(\epsilon_i=-1\) is valid and every
\(\mu_i=1/2\).  Thus every putative pair has \(\nu=1/4\), a square.  This
does not obstruct the direction \(F\): replacing \(x\) by \(cx\) changes
all evaluations to \(2c\), and Section 3 guarantees that some magnitude
works.  The example is only a scope guard against silently freezing the
fixed-edge magnitude.

## Reproduction

```text
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
  tests/test_adaptive_mobius_pairing.py

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
  /home/nick/.venvs/mo-exact/bin/python \
  src/e1_gmin_m4_adaptive_mobius_pairing.py
```

These are formula replays at one field, not a prime or target census.
