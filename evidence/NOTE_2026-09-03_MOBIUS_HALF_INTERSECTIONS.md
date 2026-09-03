# Exact intersections of two Mobius halves

Date: 2026-09-03

Status: proved a uniform, sharp two-half intersection theorem. Two localized
halves aimed at distinct hard directions share at most four inversion orbits
in total and at most two with opposite orientation. Hence a pair of trades
can cancel at most two orbits, and this is attained over every odd prime at
least five. The induced multitrade pair bound is too weak to settle branch C.
The forced fixed-edge objective and the symmetric Boolean fibre remain open.

## 1. Normal form

Take two distinct target functionals L1,L2 with nonzero star centers j1,j2,
and normalize coordinates by

\[
                     X=L_1/j_1,\qquad Y=L_2/j_2.
\]

Every auxiliary independent of its target has unique parameters
\(q,r,A,B\in\mathbf F_p\), with \(A,B\ne0\), such that

\[
 {M_1\over j_1}=(1-q/A)X+(1/A)Y,\qquad
 {M_2\over j_2}=(1/B)X+(1-r/B)Y.                 \tag{1}
\]

The closed Mobius parameterization then gives the two selected halves

\[
\begin{aligned}
 E_1(t)&=\{(1,q-A/(t+1)),(t,qt)\},\\
 E_2(s)&=\{(r-B/(s+1),1),(rs,s)\},
\end{aligned}
\qquad t,s\ne-1.                                  \tag{2}
\]

Each \(E_i\) contains one representative of every inversion orbit used by
that half. The full signed trade is
\(\mathbf1_{-E_i}-\mathbf1_{E_i}\).

## 2. Four forced candidates

A common inversion orbit satisfies

\[
                         E_1(t)=\sigma E_2(s),
                         \qquad \sigma\in\{+1,-1\}.       \tag{3}
\]

For each sign there are only two endpoint matchings.

For the direct matching, the first coordinate of the first endpoint and the
second coordinate of the second endpoint force

\[
             t={A\over q-\sigma}-1,\qquad
             s={B\over r-\sigma}-1.                       \tag{4}
\]

If either denominator vanishes, no direct solution exists. The other two
coordinates accept or reject this single candidate.

For the swapped matching, the fixed coordinates instead force

\[
                         t={\sigma\over q},\qquad
                         s={\sigma\over r}.                \tag{5}
\]

If \(q=0\) or \(r=0\), no swapped solution exists. Again the remaining
coordinates only accept or reject the forced candidate. Excluding
\(t=-1\) and \(s=-1\) completes the classification. Therefore

\[
 \boxed{\#(E_1/\{\pm1\}\cap E_2/\{\pm1\})\le4,\qquad
        \#\{\hbox{oppositely oriented common orbits}\}\le2.}       \tag{6}
\]

The distinction in (6) is essential. A same-orientation overlap adds the
two trade coefficients to plus or minus two, so the sum of two trades is not
ternary. An opposite-orientation overlap cancels the entire inversion orbit.
Thus two trades whose sum is ternary cancel at most two orbits.

## 3. Sharpness and the old one-origin choice

The earlier common-origin construction is

\[
                         q=r=0,\qquad A=B=1.
\]

Equations (4)--(5) give exactly one accepted opposite-orientation candidate,
at \(t=s=0\).

The bound two is uniformly sharp. For every odd prime \(p\ge5\), take

\[
                  q=r={1\over2},\qquad A=B={3\over2}.
\]

Both normalized auxiliaries equal \((2/3)(X+Y)\). There are exactly two
opposite-orientation intersections and no same-orientation intersections:

\[
\begin{array}{c|c|c}
(t,s)&E_1(t)&-E_2(s)\\ \hline
(0,0)&\{(1,-1),(0,0)\}&\{(1,-1),(0,0)\}\\
& &\\[-2ex]
(-2,-2)&\{(1,2),(-2,-1)\}&\{(1,2),(-2,-1)\}.
\end{array}                                               \tag{7}
\]

The two-trade sum is consequently ternary and loses four selected-orbit
occurrences, or two cancellation units.

### The two-cancellation locus is rigid

There is no free auxiliary parameter behind (7). To impose both opposite
matchings, let \(t_d,s_d\) denote the direct candidate. Its two remaining
endpoint equations are

\[
                         t_d+r s_d=0,\qquad q t_d+s_d=0.   \tag{8}
\]

The opposite swapped candidate exists only when \(q,r\ne0,1\), and its
remaining equations give

\[
 A={(q-1)(qr-1)\over qr},\qquad
 B={(r-1)(qr-1)\over qr}.                                 \tag{9}
\]

In particular \(qr\ne1\). Equation (8) therefore forces
\(t_d=s_d=0\), so the direct candidate also gives
\(A=q+1\), \(B=r+1\). Equating these with (9) yields

\[
                    q=r=1-2qr,\qquad
                    (2q-1)(q+1)=0.                        \tag{10}
\]

The root \(q=r=-1\) would force \(A=B=0\), contrary to independence of the
auxiliaries. The only admissible point is therefore

\[
                    \boxed{q=r=1/2,\qquad A=B=3/2.}       \tag{11}
\]

In the original coordinates this forces

\[
 M_1={2\over3}\left(L_1+{j_1\over j_2}L_2\right),\qquad
 M_2={2\over3}\left({j_2\over j_1}L_1+L_2\right).
\]

At this point the two same-orientation candidates are rejected, as (7)
already verifies. Thus the sharp locus is a single normalized point, not a
family. There is no parameter available for a greedy choice avoiding the
supports of previously paired trades, so the proposed free-locus pairing
argument stops here.

For scale only, the unproved disjoint-pair scenario would have
\(\kappa=m\) and \(|U|=m(p-3)\). At \(p=31\), this would mean
\(\kappa=16\), \(|U|=448\), and would pass the support count only for
\(t\ge162\). Even then the exact fixed-word gate would be

\[
 \left|a_Y+\sum_{O\in U}\Phi(O)\right|
       \le |H|-|U|=2t-323,
\]

an odd allowance ranging from one at \(t=162\) to 31 at \(t=177\).
Rigidity proves neither mutual disjointness nor this fixed-word inequality.

## 4. The induced multitrade bound

Suppose \(m\) localized trades aimed at distinct hard directions sum to a
ternary antisymmetric lift. At an inversion orbit \(O\), let \(k_+(O)\) and
\(k_-(O)\) be the numbers of selected halves using the two orientations.
Ternarity says
\(\lvert k_+(O)-k_-(O)\rvert\le1\). If \(\kappa_O\) denotes the number of
cancellation units at that orbit, then

\[
 \kappa_O={k_+(O)+k_-(O)
              -|k_+(O)-k_-(O)|\over2}
          =\min(k_+(O),k_-(O))
          \le k_+(O)k_-(O).                               \tag{12}
\]

Charge the last product to oppositely oriented pairs of trades. By (6), each
pair contributes at most two charges. Hence

\[
                         \boxed{\kappa\le
                         2{m\choose2}=m(m-1).}             \tag{13}
\]

For branch C, \(p=4r+3\), \(m=2r+2\), and

\[
 t_{\min}=2r^2-4r-2,\qquad
 t_{\max}=4r^2-2r-5.
\]

The raw number of selected-orbit occurrences and the target edge count are

\[
 N=m(p-1)={p^2-1\over2},\qquad
 |H|=4p+2t+1=N-\{2(t_{\max}-t)+1\}.                       \tag{14}
\]

Thus mere size feasibility requires

\[
                         \kappa\ge t_{\max}-t+1.           \tag{15}
\]

At \(p=31\), (15) ranges from 110 at \(t=68\) to 1 at \(t=177\), while
(13) gives 240. Therefore the exact pair theorem does not rule out the needed
cancellation anywhere on this ray. It identifies the first place where a
genuinely global compatibility theorem, rather than another two-trade
estimate, is needed.

## 5. Why support still is not the feasibility objective

After cancellations, \(|U|=N-2\kappa\). Fixed-edge elimination strengthens
the scalar test to

\[
 \boxed{|U|+
   \left|a_Y+\sum_{O\in U}\Phi(O)\right|\le |H|,}          \tag{16}
\]

with the same parity as \(|H|\). Here \(\Phi(O)\) is zero for a parallel
midpoint/difference pair and otherwise is its paired affine-line word in the
fixed-edge coordinates. The intersection classification determines which
orbits leave \(U\), and therefore determines \(|U|\), but it does not
determine the binary affine-line sum in (16) for an arbitrary collection of
trades. In particular, (13) cannot be promoted to a symmetric-completion
theorem.

The next exact object is therefore a multitrade construction or obstruction
that controls the affine-line word simultaneously with its cancellations.
Even success at (16) would only pass the first feasibility test; the
punctured halved Boolean fibre would still remain to be solved.

## Reproduction

The implementation evaluates only the four symbolically forced candidates.
The p=31 checks are fail-when-wrong witnesses, not a prime census.

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python -m pytest -q -n 0 \
      tests/test_mobius_half_intersections.py

    PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. \
      /home/nick/.venvs/mo-exact/bin/python \
      src/e1_gmin_m4_mobius_half_intersections.py
