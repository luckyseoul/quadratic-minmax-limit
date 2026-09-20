# Superlevel discrete derivatives (direction 2)

**Status:** identities D1–D4 proved; Hamming-1 superlevel corollary proved.
The plus-I / non-plus-I fork is **not** proved. Original limit OPEN.
This is not a matrix census and not Paley E(1).

OpenAI `suggest_direction` (2026-09-19): a stability theorem for the
one-sided supremum of \(Q_A\) at scale \(n^{3/2}\), from discrete
derivatives and local-optimality of superlevel points, not Kindler–Safra.

Convention: \(Q_A(x)=\sum_{i<j}A_{ij}x_ix_j=\tfrac12 x^\top Ax\),
\(\Phi(A)=\max|Q_A|\). Write \(w_i(x)=x_i(Ax)_i\).

## D1 (coordinate flip)

Let \(x^{(k)}\) be \(x\) with coordinate \(k\) flipped. Then
\[
Q_A(x^{(k)})=Q_A(x)-2w_k(x).
\]
*Proof.* The pairs incident to \(k\) contribute \(x_k(Ax)_k=w_k(x)\) to
\(Q_A(x)\). Flipping \(k\) negates that block, so the change is
\(-2w_k(x)\). \(\square\)

## D2 (1-opt / maximizer fields)

If \(x\) maximises \(Q_A\) then \(w(x)\ge0\) coordinatewise, and
\(\Phi(A)=Q_A(x)=\tfrac12\sum_i w_i(x)=\tfrac12\|Ax\|_1\).

## D3 (two-coordinate flip)

\[
Q_A(x^{(ij)})=Q_A(x)-2w_i(x)-2w_j(x)+4A_{ij}x_ix_j.
\]
If \(x\) maximises \(Q_A\) then \(A_{ij}x_ix_j\le(w_i+w_j)/2\).

## D4 (edge flip of the signing)

Let \(A^{(uv)}\) flip the single pair \(\{u,v\}\). Then
\[
Q_{A^{(uv)}}(x)=Q_A(x)-2A_{uv}x_ux_v.
\]
If \(A\) is a local minimiser of \(\Phi\) in the edge-flip graph, then
\(\Phi(A^{(uv)})\ge\Phi(A)\). Combined with D1 this is the maximizer-balance
gate of Prop 15.21, not a new covering argument.

## Hamming-1 superlevel

Let \(x\) maximise \(Q_A\), \(\Phi=Q_A(x)\), and
\(S(t)=\{y:Q_A(y)\ge\Phi-t\}\). Then the Hamming-1 neighbours of \(x\)
that lie in \(S(t)\) are exactly the flips of coordinates with
\(w_k(x)\le t/2\).

Call those coordinates **light** at width \(t\). The light set at
width \(O(1)\) is the discrete-derivative kernel of the maximizer.

This is the §12 superlevel invariant specialised to Hamming-1, not a
junta theorem.

## D5 (plus-I unfolding)

Let \(B\in\mathcal S_m\) with \(m\) odd, \(K=K(B)\), and \((x,y)\in\{\pm1\}^m\times\{\pm1\}^m\).
Set \(s=x\circ y\) and \(T=\{i:s_i=-1\}\). Then
\[
Q_{K}(x,y)=2Q_B(x)-4Q_{B[T]}(x|_T)+m-2|T|,
\]
where \(Q_{B[T]}=0\) if \(|T|<2\). In particular
\[
\Phi\bigl(K(B)\bigr)=\max_{x\in\{\pm1\}^m,\,T\subseteq[m]}
\bigl|2Q_B(x)-4Q_{B[T]}(x|_T)+m-2|T|\bigr|.
\]

*Proof.* Expand \(Q_K=Q_B(x)-Q_B(y)+x^TBy+x\cdot y\). Substitute
\(y_i=s_ix_i\). The coefficient of \(B_{ij}x_ix_j\) becomes
\(1-s_is_j+s_i+s_j\), which equals \(2\) unless both endpoints lie in
\(T\), in which case it equals \(-2\). The sum of those contributions is
\(2(Q_B(x)-Q_{B[T]})-2Q_{B[T]}\). The remaining inner-product term is
\(\sum s_i=m-2|T|\). \(\square\)

**Corollary (T-extremes).** \(T=\emptyset\) gives \(Q_K(x,x)=2Q_B(x)+m\).
\(T=[m]\) gives \(Q_K(x,-x)=-2Q_B(x)-m\). Hence
\[
\Phi\bigl(K(B)\bigr)\ge 2\Phi(B)+m.
\]
Equality holds at \(n=10\) (\(B=C_5\), \(\Phi(B)=4\), \(\Phi(K)=13\)).
At \(n=26\) one has \(2\Phi(B)+m=53<61=\Phi(K)\): a proper \(T\) supplies
the extra \(8\).

This is the plus-I transfer identity. It does **not** by itself bound
\(\Phi(K)\) by \(\Phi(C_{2m})-o(n^{3/2})\); that still needs a bound on
the max over \(T\).

The width-2 light set of a \(K\)-maximizer is **not** identical to \(T\)
in general (false at the n=26 maximizer). D5 does not claim that.

## D9 (T-gain form of D5)

The second-block fields on \(T\) are algebraic, for every \((x,y)\):
if \(j\in T\) then \(w^y_j=-2d_j^{T,\mathrm{int}}-1\), where
\(d_j^{T,\mathrm{int}}=\sum_{k\in T\setminus\{j\}}B_{jk}x_jx_k\).
Substituting into D5 gives the identity (all \(x,y\), not only maximizers)
\[
Q_K(x,y)=2Q_B(x)+m+\sum_{i\in T}(w^y_i-1).
\]
Write \(\mathrm{Extra}(x,y)=\sum_{i\in T}(w^y_i-1)\). At a maximizer of
\(Q_K\) with \(w\ge1\) one has \(\mathrm{Extra}\ge0\), recovering
\(\Phi(K)\ge2\Phi(B)+m\).

If a maximizer of \(K\) has \(|Q_B(x)|=\Phi(B)\), then
\(\mathrm{Extra}=\Phi(K)-2\Phi(B)-m\). Such an \(x\) exists among the
n=26 maximizers (\(|Q_B(x)|\in\{4,8,12,20\}\), and \(20=\Phi(B)\)).

On Paley orders \(n=2m=p^2+1\), \(\Phi(C)=mp\), and
\[
\Phi(K)=\Phi(C)-(p-1)
\quad\Longleftrightarrow\quad
\mathrm{Extra}=(m-1)(p-1)-2\Phi(B).
\]
This Extra value is attained at the three exact lifts (0, 8, 24 for
\((m,p,\Phi(B))\in\{(5,3,4),(13,5,20),(25,7,60)\}\)). That is a
certificate of the criterion, not a proof for general \(p\).

The remaining plus-I inequality is to bound Extra so that
\(\Phi(K)\ge\Phi(C)-o(n^{3/2})\) (cannot undercut Paley by a definite
\(c\)), which the operator-norm bound already gives for near-min-op
\(B\); or Extra exactly equal to that display, which would make the
gap \(p-1\) a theorem along this family.

## D10 (Extra as induced plus-count)

Let \(x\) maximise \(Q_B\) and \(\beta_{ij}=B_{ij}x_ix_j\). For
\(T\subseteq[m]\) with \(t=|T|\ge 2\), let \(e_+(T)\) be the number of
pairs in \(T\) with \(\beta_{ij}=+1\), and set \(s_i=-1\) on \(T\),
\(+1\) off \(T\), \(y=s\circ x\). Then
\[
\mathrm{Extra}(x,y)=2t(t-2)-8e_+(T),
\]
and therefore
\[
Q_K(x,y)=2\Phi(B)+m+2t(t-2)-8e_+(T).
\]
*Proof.* \(Q_{B[T]}(x|_T)=\sum_{i<j,\,i,j\in T}\beta_{ij}=e_+-e_-\) and
\(e_++e_-=\binom t2\), so \(Q_T=2e_+-\binom t2\). D9/D5 give
\(\mathrm{Extra}=-4Q_T-2t=4\binom t2-8e_+-2t=2t(t-1)-8e_+-2t=2t(t-2)-8e_+\). \(\square\)

Hence, for every \(B\),
\[
\Phi\bigl(K(B)\bigr)\ge 2\Phi(B)+m+\Gamma(B),\qquad
\Gamma(B):=\max_{T\subseteq[m]}\bigl(2|T|(|T|-2)-8e_+(T)\bigr),
\]
the max taken in the switching that makes a maximizer of \(B\) into \(1\).
(The \(t\le 1\) terms contribute \(0\).)

## D11 (Case A: large \(\Phi(B)\) needs no Extra)

Let \(n=2m=p^2+1\) with \(p\) an odd prime, and let \(C\) be Paley
conference of order \(n\), so \(\Phi(C)=mp\). If
\[
\Phi(B)\ge\frac{(m-1)(p-1)}{2},
\]
then already \(2\Phi(B)+m\ge mp-(p-1)=\Phi(C)-(p-1)\), hence
\(\Phi(K(B))\ge\Phi(C)-(p-1)\). *Proof.* Rearrange
\(2\Phi(B)+m-(mp-(p-1))=2\Phi(B)-(m-1)(p-1)\). \(\square\)

## D12 (\(p=3\) is entirely Case A)

For \(p=3\), \(m=5\), \(\frac{(m-1)(p-1)}{2}=4\). Exact \(m_5=4\), so
every \(B\in\mathcal S_5\) has \(\Phi(B)\ge4\). D11 gives
\(\Phi(K(B))\ge 13=\Phi(C_{10})-2\). Equality holds for \(B\) the
5-cycle. Thus every plus-I lift of order 10 satisfies
\(\Phi(K)\ge\Phi(C_{10})-(p-1)\), and the 5-cycle lift matches \(m_{10}\).

## D13 (bipartite minus graph is Case A)

Let \(B\) be switched so that \(1\) maximises \(Q_B\). If the minus graph
of \(B\) is bipartite with parts \(X,Y\), then every minus edge lies in
the cut \((X,Y)\). The 1-maximizer cut condition gives
\(2e_+^{\mathrm{cut}}\ge|X||Y|\), hence \(e_-\le|X||Y|/2\) and
\[
\Phi(B)=\binom{m}{2}-2e_-\ge\binom{m}{2}-|X||Y|\ge\frac{(m-1)^2}{4}
\]
(the last step uses \(|X||Y|\le\lfloor m^2/4\rfloor\)). For
\(m=(p^2+1)/2\) one has \((m-1)^2/4\ge(m-1)(p-1)/2\) iff \(p\ge3\).
Thus a bipartite minus graph forces Case A, so it is **not** a Case B
counterexample. (This is the missing 1-maximizer input that Mantel lacked.)

## D14 (\(p=5\), \(M_+=22\): Extra bound holds)

D10 is evaluated at a maximizer of \(Q_B\), not of \(|Q_B|\). Write
\(M_+=\max Q_B\). Switching that vector to \(1\) gives \(M_+=Q(1)=78-2e_-\)
and Hamming-1 minus \(\Delta\le 6\). Same-sign / \(|Q|\) cuts are not used.
Case A is \(M_+\ge 24\); Case B is \(M_+\le 22\), and D9/D10 give
\(\Phi(K)\ge 2M_++13+\Gamma\) so the Paley-gap target \(61\) needs
\(\Gamma\ge 48-2M_+\). Extra never equals \(2\) or \(4\) (\(t\) even
\(\Rightarrow\) Extra \(\equiv 0\pmod 8\); \(t\) odd \(\Rightarrow\)
Extra \(\equiv 6\pmod 8\)). A minus triangle has Extra \(=6\). Mantel:
triangle-free minus \(\Rightarrow\Gamma=0\).

nauty `geng -t -D6 13 28:39` enumerates every unlabeled triangle-free
13-vertex graph with \(\Delta\le 6\) and \(e\ge 28\) (163477 graphs:
100457 with 28 edges, then 39920, 15057, 5355, 1840, 587, 186, 54, 18, 3).
Each was tested against the \(2^{12}\) Q-max cut inequalities
\(e_-^{\mathrm{cut}}\le |U||U^c|/2\). Zero passed. Self-check: the empty
minus graph is accepted as a Q-maximizer, the complete minus graph is
rejected, and the switched \(m_{13}=20\) block (\(M_+=20\), \(e_-=29\),
a minus triangle, \(\Gamma=8\)) is accepted. OpenAI `math_review` PASS
on this reduction (after a BLOCK on the false \(|Q|\)-maximizer form).

Hence there is no triangle-free Case B Q-maximizer at \(m=13\). Every
Q-maximizer with \(M_+=22\) has a minus triangle, so \(\Gamma\ge 6\ge 4\),
and \(\Phi(K)\ge 2\cdot 22+13+6=63\ge 61=\Phi(C_{26})-(p-1)\).

The case \(M_+\le 20\) (need \(\Gamma\ge 8\)) is open. The known order-13
minimizer attains equality (\(\Gamma=8\)). A counterexample would have
to be diamond-free (\(K_4-e\)-free) with a triangle. Hamming-2 forces
the deg-6 vertices to be a minus clique; a deg-6 triangle is
diamond-free-impossible (12 disjoint outer neighbours, 10 vertices
left), so a counterexample has at most two deg-6 vertices. A 75-worker
search of that class checked 17.4 million 29-edge leaves and found no
Q-maximizer, but 36/75 shards hit a node cap, so that search is not a
proof. Receipts: `evidence/case_b_gamma_p5/tf_qmax_census.json`,
`scripts/case_b_tf_geng_check.py`. Limit OPEN.

## Case B (open for \(p\ge 5\) except \(p=5\), \(M_+=22\))

If \(M_+<(m-1)(p-1)/2\), D11 does not apply and \(\Gamma\) must supply
the deficit \((m-1)(p-1)-2M_+\). Taking \(|T|=p-1\), D10
says this holds as soon as some \((p-1)\)-set has
\[
e_+(T)\le\frac{2(p-1)(p-3)-(m-1)(p-1)+2M_+}{8}.
\]
At the n=26 minimizer this threshold is \(1\), and a 4-set with
\(e_+=1\) exists (Extra \(=8\)). That is a certificate for one matrix,
not a proof for every \(B\in\mathcal S_m\). At \(p=5\), D14 kills the
triangle-free / \(M_+=22\) slice. Remaining: \(M_+\le 20\) at \(m=13\),
and all Case B at \(p\ge 7\).

## Fork (unproved)

Let \(n_k=p_k^2+1\) and \(\Delta(A)=\Phi(C_{n_k})-\Phi(A)\) with \(C\)
Paley conference. A fixed-deficit sequence \(\Delta(A_k)=O(n_k^{1/2})\)
(the plus-I scale) must, after switching and permutation, either

1. admit a plus-I decomposition \(K(B)\) together with a sharp identity
   relating \(\Delta(K(B))\) to data of \(B\), or
2. contain a stable non-plus-I light-set template that lifts along a
   ratio-dense subsequence.

Neither arm is proved. Arm 1 is a transfer identity, not another
\(\Phi(K_{50})\) pin. Arm 2 is forbidden to become a Hadamard/family
scan or an E(1) Max+ cover.

## Certificates (existing matrices, not a search)

On already-certified maximizers, D1 holds and the light sets at width 2
(i.e. \(w_k\le1\)) are:

| matrix | \(\Phi\) | \(w\)-multiset | light \(\{w\le1\}\) |
|---|---|---|---|
| Paley \(n=10\) | 15 | \(3^{10}\) | empty |
| plus-I \(C_5\) | 13 | \(1^6,5^4\) | 6 |
| plus-I \(n=26\) | 61 | \(1^6,5^{16},9^4\) | 6 |

Paley maximizers have no width-2 light coordinates (\(w\equiv p\)).
The known beaters do. That is the discrete-derivative distinction the
fork has to explain, not a new finite record.

Tests: `tests/test_superlevel_discrete_derivative.py`.
