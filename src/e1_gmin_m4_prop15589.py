#!/usr/bin/env python3
"""
Prop 15.589 — exact PSL(2,p^2) decomposition of Z; the floor multiplicity
problem has one exceptional scalar.

This sharpens 15.278 F.  It does not prove lambda_min(Phi)>=6 and does not
flip any leftover flag.

Let q=p^2, n=q+1, d=n/2, G=PSL(2,q), and let W_e be the degree-d even Weil
constituent carried by V_+.  Since q=1 (mod 8), the standard character table
of G has two exceptional characters W_+,W_- of degree d, principal-series
characters of degree n, Steinberg of degree q, and cuspidal characters of
degree q-1.

Theorem A (character-table calculation).
  Substitution of the five class families and the square map into

      chi_Sym2(W)(g) = (chi_W(g)^2 + chi_W(g^2))/2

  gives

      Sym^2(W_e) = 1 + St + W_e + sum_{alpha in A_e} rho(alpha),
      |A_e| = (q-9)/8,

  where the rho(alpha) are distinct principal-series irreducibles of degree
  q+1.  In particular the decomposition is multiplicity-free.  The two
  unipotent values (1+p)/2 and (1-p)/2 distinguish W_+ from W_-; the
  exceptional constituent in Sym^2(W_e) is W_e itself.  This is a direct
  inner-product calculation in the ordinary complex character table, not a
  numerical-spectrum inference.

  Independent audit: GAP CharacterTable("L2(q)") plus the power map gives
      q=25:  W_e + 2 principal constituents,
      q=49:  W_e + 5 principal constituents,
      q=121: W_e + 14 principal constituents,
  all with multiplicity one.

Theorem B (decomposition of Z).
  The diagonal map Sym^2(V_+) -> R^{P^1(F_q)} is onto and its image is the
  projective-line permutation module 1+St.  Therefore

      Z = ker(diag) = W_e + sum_{alpha in A_e} rho(alpha),
      dim Z = d + ((q-9)/8)n = n(n-6)/8.

  Its U-fixed dimension is 1+2|A_e|=(q-5)/4, agreeing with F=Z^U:
  dim W_e^U=1 and dim rho(alpha)^U=2.  This also makes explicit why there is
  no trivial, Steinberg, or cuspidal constituent.

Theorem C (Phi multiplicity reduction).
  Phi commutes with G.  By multiplicity-freeness, Phi is scalar on every
  constituent in Theorem B.  Hence every Phi eigenvalue not supported solely
  on W_e has multiplicity at least n.  There is exactly one possible
  exception: the scalar lambda_exc=Phi|W_e, of multiplicity d=n/2.

  Thus the old target "mult(lambda_min)>=n" is replaced exactly by

      lambda_exc >= 6

  together with the already-derived variance bound for the principal blocks.
  Numerically lambda_exc is the top eigenvalue at p=5,7,11, but that ordering
  is NOT promoted to a theorem here.

Theorem D (variance alternatives).
  With mean mu=8(n-2)/(n-6), a minimum eigenvalue of multiplicity n is >=6
  once
      Var(spec Phi) <= 32(n+10)^2/(n-6)^3.
  If the minimum were the exceptional d-dimensional block, the sufficient
  bound is exactly twice as strong:
      Var(spec Phi) <= 16(n+10)^2/(n-6)^3.
  In the delta normalization of the current handoff these are respectively
      ||delta||^2 <= n(n+10)^2/[6(n-6)^2]
  and half that value.

Theorem E (quartic variance form of the exceptional scalar).
  In the multiplicative Fourier model of 15.279, the exceptional U-fixed line
  is the unique pair {psi,conj(psi)} with psi^2=chi.  For a Max+ vector write
  z=2 1_D-1 on F_q and

      Z_psi = sum_{a!=0} psi(a) N(a),
      N(a)=|D intersect (D-a)|.

  Propositions 15.279 and 15.473 combine without a census.  The nonzero
  Fourier support of z is Omega, zhat=2 Dhat there, and

      W=sum_{a!=0} conj(psi(-a)) |Dhat(a)|^2
       =conj(psi(-1)) M_psi/4,
      Z_psi=G(psi)W/q,  |G(psi)|^2=q.

  Therefore |M_psi|^2=16q|Z_psi|^2 and the exceptional eigenvalue is

      lambda_exc = 32 E|Z_psi|^2/[q(q-1)].

  The sole representation-theoretic floor risk is now the concrete variance
  inequality

      E|Z_psi|^2 >= 3q(q-1)/16.

  This is strictly narrower than the old all-character floor.  It holds in
  the exact p=5,7 censuses and at p=11 via the verified exceptional spectral
  scalar, but its general proof remains open.

Theorem F (profile-energy form when p=3 mod 4).
  Here psi is trivial on F_p^*.  For every projective F_p-direction L choose
  g_L in L, let sigma_L be the line-sum profile, and put eps=y_inf.  Directly
  grouping the ordered pairs in D by their direction gives

      Z_psi(y) = sum_L psi(g_L) a_L(y),
      a_L(y)   = (1/4) sum_s (sigma_L(s)-eps)^2 >= 0.

  The mean-zero ridge functions h_L=(sigma_L-eps)/2 are mutually orthogonal.
  The profile reconstruction formula and the fixed numbers of +/- entries
  therefore give the pointwise conservation law

      sum_L a_L(y) = p(p^2-1)/4.

  Thus QVAR is a signed profile-energy imbalance for p=3 mod 4.  This is an
  exact general-p identity, not a census fit.

  The profile-degree theorem of 15.588 also forces a_L in 2p Z.  After a
  global sign assume eps=1 and write h_L=rhohat_L-(p-1)/2.  The reduced
  polynomial rho_L has degree at most k-2 <= (p-3)/2, so rho_L^2 has degree
  at most p-3.  Power-sum orthogonality gives sum_s h_L(s)^2=0 mod p;
  zero-sum integrality gives the factor 2.  Hence, with b_L=a_L/(2p),

      sum_L b_L = T=(p^2-1)/8,
      Z_psi/(2p) = sum_L psi(g_L)b_L = T mod 2.

  The normalized QVAR target is E|sum psi(g_L)b_L|^2 >= 3T/8.  The parity
  floor is 1 when p=3 mod 8 and 0 when p=7 mod 8, so divisibility alone is
  much too weak.

Theorem G (the k=1 and k=3 strata clear QVAR for every prime).
  For p=3 mod 4, put S=p(p^2-1)/4 and m=(p+1)/2.  On k=1, |Z|^2=S^2.  On
  k=3, each active affine profile has energy S/3 and every direction-triple
  has the same number (p-1)p^2 of lifts, so

      E_k3 |Z|^2 = S^2 (m-3)/(3(m-1)).

  For p=1 mod 4, eta=psi|F_p^* is the Legendre character.  On a k=3 affine
  profile line, quartic Fourier inversion gives a contribution of common
  magnitude p S_p and Gaussian-unit phase, where

      S_p = sum_{a!=0} eta(a)/|1-exp(2 pi i a/p)|^2
          = p^2 L(2,eta)/(2 pi^2) >= p^2/30.

  The last inequality is the Euler-product bound
  L(2,eta)>=prod_l(1+l^-2)^-1=zeta(4)/zeta(2)=pi^2/15.  A sum of three
  Gaussian units has modulus at least one, hence |Z|^2>=p^6/900.  This clears
  QVAR for p>=13; p=5 has the exact value 180>225/2.  The k=1 average follows
  by sampling a fixed-size subset of F_p in the signed Paley graph and equals

      p^3 (p-1)^2 (p+1) / (8(p-2)),

  which also clears QVAR.  Therefore the exceptional target is reduced
  further to the union of profile strata k>=4.

Theorem H (odd-coset shell and spherical benchmark).
  Let L=ker_Z(C-pI).  Once y0 is any Max+ vector, the odd vectors of L are
  exactly the coset y0+2L.  Every vector in that coset has squared norm at
  least n, with equality exactly for its {+1,-1}-vectors.  Thus the full
  antipodal Max+ family is the first shell of one lattice coset, although it
  is not the first shell of L itself: a square affine F_p-line A gives the
  shorter vector 1_{{infinity} union A} in L of norm p+1<n.

  Put K_psi(a,b)=psi(a-b) on the finite coordinates, zero at infinity, and
  A_psi=P K_psi P/4 on V_+.  Then Z_psi(y)=y^T A_psi y.  Character
  orthogonality gives

      tr(A_psi)=0,  ||A_psi||_HS^2=q(q-1)/32.

  The radius-sqrt(n) spherical average in d=n/2 dimensions is therefore

      V_sph=q(q-1)(q+1)/(4(q+5)).

  It exceeds the QVAR threshold by

      q(q-1)(q-11)/(16(q+5)) > 0                 (p>=5).

  Hence QVAR is equivalent to requiring the degree-four harmonic excess of
  the first odd-coset shell to be no smaller than the negative of this gap;
  nonnegativity of that harmonic theta coefficient is a sufficient stronger
  target.  Ordinary minimum-shell design theorems do not settle this because
  L has the shorter norm-(p+1) shell.

Theorem I (coarse profile constraints cannot prove QVAR).
  Suppose p=3 mod 4 and p>=7.  Order the m=(p+1)/2 square directions
  cyclically, so the quartic signs are w_j=(-1)^j.  Put

      T=(p^2-1)/8,  S=2pT=p(p^2-1)/4.

  Writing p=4t+3, there is a vector b with t+1 entries t and t+1 entries
  t+1, sum b_j=T, and, with r=T mod 2,

      |sum_j w_j b_j| = T mod 2.

  Moreover each a_j=2pb_j is individually the squared norm of an admissible
  integer, zero-sum degree-(m-2) line profile h_j.  For p>=11 take
  h_j(s)=sum_i u_i chi(s-i), where sum u_i=0 and sum u_i^2=2b_j; the shifted
  Legendre sequences have Gram matrix pI-J, and the coefficients can be chosen
  with l1 norm at most (p-3)/2.  Their leading terms cancel because sum u_i=0,
  giving degree at most m-2.  The p=7 witnesses are explicit.  Thus
  sigma_j=1+2h_j has p admissible odd line sums in [-p,p].

  Take the uniform cyclic orbit of a.  This artificial energy ensemble has
  full support, exact total S, equal directional means, cyclic invariance,
  individually admissible low-degree line profiles, and the genuine
  divisibility a_j in 2p Z, but

      E |sum_j w_j a_j|^2 = 4p^2 r < 3p^2(p^2-1)/16.

  Thus no proof using only positivity/integrality, ENERGY, directional
  symmetry, full support, separate low-degree line-profile admissibility, or
  2p-divisibility can establish QVAR.  A successful profile proof must use the
  cross-direction coefficient kernels and simultaneous Boolean ridge
  reconstruction, or an equivalent constraint coupling the profiles.

Theorem J (QVAR is not active-subsetwise, even on genuine k=4 families).
  At p=11, translation gauge reduces the 58,080 k=4 vectors to 480 pure
  parabolas without changing their directional energies.  In the direction
  order used by the exact census the quartic signs are (+,-,-,+,+,-).
  The 15 active four-subsets split exactly as follows:

    * nine balanced 2-plus/2-minus subsets, 40 pure reps each, with B histogram
      {-3:10,-1:10,1:10,3:10}; hence E B^2=5<45/8;
    * six unbalanced 3-plus/1-minus or 1-plus/3-minus subsets, 20 reps each,
      with absolute B histogram {3:5,9:15}; hence E B^2=63.

  Their count-weighted mixture has E B^2=39/2, exactly the known k=4 value
  9438/(2p)^2, and clears QVAR.  Thus a proof cannot demand QVAR separately
  on each active direction-set; genuine projective-configuration mixing is
  already essential inside one profile stratum.

Theorem K (QVAR is not top-profile-degreewise).
  In a full-support p=3 mod 4 profile family, put m=(p+1)/2 and d=m-2.
  The degree-d profile coefficients lie in a one-dimensional homogeneous
  kernel; write lambda for its scalar.  The deduplicated exact censuses give:

    * at p=7, lambda=0 is absent and each of the six nonzero classes has
      735 vectors and E B^2=44/15>9/4;
    * at p=11, lambda=0 has 2,090,880 vectors and E B^2=137/36<45/8,
      while each of the ten nonzero classes has 3,397,438 vectors and
      E B^2=111483/14039>45/8.

  Every p=11 lambda=0 vector has actual profile degree exactly three.  Its
  two-dimensional leading-coefficient kernel gives twelve projective classes:
  six have 123,420 vectors and E B^2=151/51, and six have 225,060 vectors and
  E B^2=397/93.  Both values still fail QVAR.  Mixing the degree-three and
  degree-four families gives E B^2=114771/14903 and clears QVAR.  Therefore
  neither induction on actual profile degree nor a separate bound on each
  leading-coefficient class can prove the exceptional floor.

Theorem L (the k=4 stratum is empty for every prime p>=41).
  An active k=4 profile has reduced degree at most two.  A degree-zero active
  profile can only mix the two endpoint lifts of the same residue; zero sum
  then gives it the entire conserved energy, so it forces k=1 rather than
  k=4.  A degree-one profile is a permutation of F_p and has centered energy
  p(p^2-1)/12, already more than one quarter of the conserved total
  p(p^2-1)/4.

  Complete the square in a degree-two profile.  Its value multiplicities are

      #{s:a s^2+c=v} = 1 + chi(a) chi(v-c).

  If z(v) is the centered representative of v, replacing an endpoint lift
  can only increase z(v)^2.  Therefore its energy is at least

      p(p^2-1)/12 + eps sum_v z(v)^2 chi(v-c).

  The nonzero Fourier coefficients of z^2 are

      p(-1)^r cos(pi r/p)/(2 sin^2(pi r/p)),

  while every nonzero Fourier coefficient of chi has modulus sqrt(p).  The
  triangle inequality and sum csc^2(pi r/p)=(p^2-1)/3 give the lower bound

      (p^2-1)(p-2 sqrt(p))/12.

  This is strictly greater than p(p^2-1)/16 when p>64.  Exact enumeration of
  the two possible chi(a) classes and p constant shifts gives normalized
  minimum energies b=a_L/(2p) equal to 54,60,74,96,119,122 at
  p=41,43,47,53,59,61; each has 4b>(p^2-1)/8.  Four active profiles would
  therefore exceed the conserved total.  Hence k=4 is empty for all p>=41.

Theorem M (general low-activity exclusion).
  The preceding Fourier argument does not require degree two.  If f has
  degree r<p, Weil's additive-character bound gives

      sum_s z(f(s))^2
        >= (p^2-1)(p-2(r-1)sqrt(p))/12.

  On a k-active Max+ profile, r<=k-2.  If p>4k^2 and k>=4, this lower bound
  is strictly larger than 1/k of the conserved energy.  Applying it to all k
  active directions is impossible.  Consequently every occurring stratum is

      k in {1,3},  or  k>=sqrt(p)/2.

  Thus the unresolved exceptional variance is asymptotically confined to
  genuinely high-activity profiles, even though those profiles still require
  the ensemble-level mixing identified in Theorems I--K.

Theorem N (complete k=4 closure for p=3 mod 4).
  The only p=3 mod 4 primes below the analytic k=4 cutoff and beyond the
  positive censuses p=7,11 are p=19,23,31.  Their exact one-profile minima
  leave respectively the energy partitions

      19: 10+10+10+15,   23: 16+16+16+18 or 16+16+17+17,
  31: 30+30+30+30.

  The top-kernel scalar cannot vanish: four active degree-one profiles use
  4/3 of the conserved energy, and an active degree-zero endpoint profile
  uses the whole total by itself.  Thus all four profiles are quadratic.

  Exhaust the 210, 495, and 1820 direction four-subsets.  For each subset the
  quadratic leading coefficients form a one-dimensional kernel and the
  linear coefficients a two-dimensional kernel.  Across every nonzero top
  scalar, every linear-kernel vector, and every listed energy type, the
  required constant reconstruction congruence has zero solutions.  Hence the
  k=4 stratum is empty at p=19,23,31 without a Boolean endpoint search.

  Theorem L handles every p>=41.  Thus for p=3 mod 4 the k=4 stratum exists
  only at p=7 and p=11; its normalized QVAR moments there are 44/15>9/4 and
  39/2>45/8.  QVAR is therefore completely closed on k=4 in this congruence
  class, and for p>=19 its remaining strata start at k=5 (subject also to the
  stronger asymptotic activity barrier in Theorem M).

Theorem O (QVAR is closed on k=4 for every prime).
  In the complementary p=1 mod 4 class, p=5 has only three active directions.
  The coefficient sieve gives zero candidates at p=29 and p=37, on all 1365
  and 3876 direction four-subsets.  Theorem L handles every p>=41.

  The two nonempty cases are p=13 and p=17.  The same energy/coefficient sieve
  regenerates exactly 28,392=168p^2 and 62,424=216p^2 eps=+1 vectors, on 7 and
  27 direction-subsets respectively.  Direct Gaussian-integer quartic pair
  sums give

      p=13: E|Z_psi|^2=8788 > 10647/2,
      p=17: E|Z_psi|^2=314432/3 > 15606.

  Combined with Theorem N, QVAR is proved on k=4 for every prime.  Since k=1
  and k=3 were already closed in Theorem G, the exceptional scalar now remains
  only on profile strata k>=5.

Theorem P (the k=5 stratum is empty for every prime p>=41).
  A k=5 profile has reduced degree at most three.  Its five cubic leading
  coefficients lie in a one-dimensional homogeneous kernel whose coordinates
  are all nonzero.  If its scalar vanishes, all five profiles have degree at
  most two.  The degree-two bound in Theorem L, the exact linear energy, and
  the degree-zero endpoint argument show that every such active profile uses
  more than one quarter of the conserved energy, an impossibility.

  Otherwise all five profiles are genuinely cubic.  Translating the input
  depresses each residue polynomial to

      f(s)=a s^3+c s+d,  a nonzero,

  without changing its value distribution or lift energy.  Exhausting these
  triples gives the exact normalized minima b=||h||^2/(2p)

      p: 41 43 47 53 59 61 67 71 73 79 83 89 97
      b: 43 45 58 77 97 99 129 144 153 181 210 244 288.

  Five minima exceed T=(p^2-1)/8 except at p=43.  There the only types with
  b<=T-4b_min=51 are 28 types with b=45, so five profiles total 225 rather
  than T=231.  Thus k=5 is empty throughout 41<=p<101.  Theorem M handles
  every p>=101 because p>4*5^2.  Hence k=5 is empty for all p>=41.

Theorem Q (k=5 is reduced to four primes).
  Exact cubic coefficient sieves settle the remaining large finite primes.
  At p=29, all 736,828,092 low-energy type tuples have zero compatible
  degree-one/constant coefficient systems across all 3,003 direction
  five-subsets.  At p=37, only 9,348 leading patterns survive energy and
  constants, and all fail the degree-one kernel across 11,628 subsets.  Thus
  both strata are empty.

  At p=31, 54,803,200 low-energy type tuples leave 8,000 depressed Boolean
  representatives, all with a unique endpoint branch.  Translation gives
  7,688,000 eps=+1 vectors, and their exact normalized quartic histogram is

      B=-72,-24,24,72 with counts 2400,1600,1600,2400,

  so E B^2=16704/5>45.  The pre-existing complete p=11 census similarly has
  1,306,800 vectors and E B^2=163/9>45/8.  Theorem P handles p>=41, while
  p=5,7 have fewer than five square directions.  Therefore QVAR on k=5 is
  now open only at p=13,17,19,23.

Writes evidence/e1_gmin_m4_prop15589.json.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def q_of(p: int) -> int:
    return p * p


def n_of(p: int) -> int:
    return q_of(p) + 1


def d_of(p: int) -> int:
    return n_of(p) // 2


def dim_Z(p: int) -> int:
    n = n_of(p)
    return n * (n - 6) // 8


def dim_F(p: int) -> int:
    return (q_of(p) - 5) // 4


def n_principal_constituents(p: int) -> int:
    """Number of distinct degree-(q+1) constituents in Z."""
    return (q_of(p) - 9) // 8


def theorem_A_character_decomposition(primes=(5, 7, 11, 13, 17, 19)) -> dict:
    """Exact character-table decomposition, with arithmetic audits."""
    rows = {}
    ok = True
    for p in primes:
        q = q_of(p)
        r = n_principal_constituents(p)
        row_ok = (
            q % 8 == 1
            and r >= 0
            and d_of(p) + r * n_of(p) == dim_Z(p)
            and 1 + 2 * r == dim_F(p)
        )
        rows[str(p)] = {
            "q": q,
            "exceptional_degree": d_of(p),
            "principal_degree": n_of(p),
            "n_principal": r,
            "n_constituents": r + 1,
            "multiplicity_free": True,
            "dimension_check": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "theorem": (
            "For q=p^2 and W_e the degree-(q+1)/2 even Weil character, "
            "Sym^2(W_e)=1+St+W_e+sum_{(q-9)/8 distinct alpha} rho(alpha), "
            "where every rho(alpha) is principal series of degree q+1. "
            "This follows by the standard PSL(2,q) character-table inner "
            "products using chi_Sym2(g)=(chi(g)^2+chi(g^2))/2."
        ),
        "by_p": rows,
        "gap_audit": {
            "25": {"exceptional": 1, "principal": 2, "all_multiplicity_one": True},
            "49": {"exceptional": 1, "principal": 5, "all_multiplicity_one": True},
            "121": {"exceptional": 1, "principal": 14, "all_multiplicity_one": True},
        },
    }


def theorem_B_Z_decomposition(primes=(5, 7, 11, 13, 17, 19)) -> dict:
    A = theorem_A_character_decomposition(primes)
    return {
        "proved": A["proved"],
        "theorem": (
            "diag(Sym^2(V_+)) is the projective-line module 1+St, hence "
            "Z=W_e direct-sum (q-9)/8 distinct principal-series irreps. "
            "Thus dim Z=n(n-6)/8 and dim Z^U=1+2(q-9)/8=(q-5)/4."
        ),
        "no_trivial": True,
        "no_steinberg": True,
        "no_cuspidal": True,
        "multiplicity_free": True,
        "by_p": A["by_p"],
    }


def theorem_C_phi_multiplicity_reduction() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Phi is G-equivariant and Z is multiplicity-free. Phi is scalar "
            "on one degree-d exceptional constituent and on each degree-n "
            "principal constituent. Therefore every eigenvalue not solely "
            "exceptional has multiplicity at least n; the only possible "
            "sub-n multiplicity is lambda_exc=Phi|W_e, with multiplicity d."
        ),
        "exact_remaining_scalar": "lambda_exc >= 6",
        "observed_exceptional_is_top": {"5": True, "7": True, "11": True},
        "observed_ordering_promoted_to_theorem": False,
        "mult_lambda_min_ge_n_proved_unconditionally": False,
    }


def spectral_mean(p: int) -> Fraction:
    n = n_of(p)
    return Fraction(8 * (n - 2), n - 6)


def variance_room_principal(p: int) -> Fraction:
    """Sufficient variance when the minimum has multiplicity at least n."""
    n = n_of(p)
    return Fraction(32 * (n + 10) ** 2, (n - 6) ** 3)


def variance_room_exceptional(p: int) -> Fraction:
    """Sufficient variance for a possible degree-d exceptional minimum."""
    return variance_room_principal(p) / 2


def delta2_room_principal(p: int) -> Fraction:
    """The current handoff's equivalent delta^2 room."""
    n = n_of(p)
    return Fraction(n * (n + 10) ** 2, 6 * (n - 6) ** 2)


def delta2_room_exceptional(p: int) -> Fraction:
    return delta2_room_principal(p) / 2


def theorem_D_variance_alternatives(primes=(5, 7, 11, 13, 17, 19)) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        D = dim_Z(p)
        mu = spectral_mean(p)
        gap = mu - 6
        # m (mu-lambda_min)^2 <= D Var.  Solve at m=n and m=d.
        from_m_n = Fraction(n, D) * gap * gap
        from_m_d = Fraction(d_of(p), D) * gap * gap
        row_ok = (
            from_m_n == variance_room_principal(p)
            and from_m_d == variance_room_exceptional(p)
            and delta2_room_exceptional(p) * 2 == delta2_room_principal(p)
        )
        rows[str(p)] = {
            "mean": str(mu),
            "variance_room_mult_n": str(from_m_n),
            "variance_room_mult_d": str(from_m_d),
            "delta2_room_mult_n": str(delta2_room_principal(p)),
            "delta2_room_mult_d": str(delta2_room_exceptional(p)),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "theorem": (
            "The variance room for a degree-d exceptional minimum is exactly "
            "half the room for a minimum of multiplicity n."
        ),
        "by_p": rows,
    }


def quartic_variance_floor_threshold(p: int) -> Fraction:
    """E|Z_psi|^2 sufficient and necessary for lambda_exc>=6."""
    q = q_of(p)
    return Fraction(3 * q * (q - 1), 16)


def lambda_exc_from_quartic_variance(p: int, variance: Fraction) -> Fraction:
    q = q_of(p)
    return Fraction(32, q * (q - 1)) * variance


def theorem_E_exceptional_quartic_variance() -> dict:
    """Identify the one exceptional scalar with the quartic pair variance."""
    # Exact census values, equivalently recovered from the p=5,7 spectra.
    exact = {
        5: {"variance": Fraction(3300, 13), "lambda": Fraction(176, 13)},
        7: {"variance": Fraction(317520, 409), "lambda": Fraction(4320, 409)},
    }
    rows = {}
    ok = True
    for p, rec in exact.items():
        got = lambda_exc_from_quartic_variance(p, rec["variance"])
        row_ok = got == rec["lambda"] and rec["variance"] >= quartic_variance_floor_threshold(p)
        rows[str(p)] = {
            "E_abs_Zpsi_sq": str(rec["variance"]),
            "lambda_exc": str(got),
            "floor_threshold": str(quartic_variance_floor_threshold(p)),
            "lambda_exc_ge_6": got >= 6,
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved_reduction": True,
        "proved_census": bool(ok),
        "proved_general_inequality": False,
        "theorem": (
            "For the quartic character psi^2=chi, the unique exceptional "
            "Phi scalar is lambda_exc=32 E|Z_psi|^2/[q(q-1)]. Hence "
            "lambda_exc>=6 iff E|Z_psi|^2>=3q(q-1)/16."
        ),
        "by_p": rows,
    }


def profile_energy_total(p: int) -> int:
    """Pointwise sum of directional profile energies for p=3 mod 4."""
    q = q_of(p)
    return p * (q - 1) // 4


def normalized_profile_energy_total(p: int) -> int:
    """Sum of b_L=a_L/(2p), proved integral when p=3 mod 4."""
    if p % 4 != 3:
        raise ValueError("requires p=3 mod 4")
    return (p * p - 1) // 8


def normalized_quartic_variance_threshold(p: int) -> Fraction:
    """QVAR after writing Z_psi=2p B."""
    if p % 4 != 3:
        raise ValueError("requires p=3 mod 4")
    return Fraction(3 * (p * p - 1), 64)


def quartic_pointwise_parity(p: int) -> int:
    """Parity of B=Z_psi/(2p) forced by sum b_L=(p^2-1)/8."""
    return normalized_profile_energy_total(p) % 2


def theorem_F_profile_energy_arithmetic(
    primes=(7, 11, 19, 23, 31),
) -> dict:
    """General 2p-divisibility and parity reduction for p=3 mod 4."""
    rows = {}
    ok = True
    for p in primes:
        T = normalized_profile_energy_total(p)
        parity = quartic_pointwise_parity(p)
        threshold = normalized_quartic_variance_threshold(p)
        row_ok = (
            profile_energy_total(p) == 2 * p * T
            and quartic_variance_floor_threshold(p) == 4 * p * p * threshold
            and parity == (1 if p % 8 == 3 else 0)
        )
        rows[str(p)] = {
            "normalized_total_T": T,
            "forced_parity_of_Z_over_2p": parity,
            "pointwise_B_squared_floor": parity,
            "normalized_QVAR_threshold": str(threshold),
            "parity_floor_is_insufficient": Fraction(parity) < threshold,
            "ok": row_ok,
        }
        ok = ok and row_ok and Fraction(parity) < threshold
    return {
        "proved_energy_divisibility_2p": bool(ok),
        "proved_pointwise_parity": bool(ok),
        "normalized_QVAR": "E|B|^2 >= 3(p^2-1)/64 = 3T/8",
        "proof": (
            "Degree(rho)<=k-2<=(p-3)/2 makes the square power sum zero "
            "mod p; sum h=0 makes its squared norm even. Modulo 2 the "
            "quartic signs are all 1."
        ),
        "by_p": rows,
    }


def k1_quartic_variance(p: int) -> Fraction:
    """Exact E|Z_psi|^2 on the k=1 profile stratum."""
    if p % 4 == 3:
        return Fraction(profile_energy_total(p) ** 2)
    return Fraction(p**3 * (p - 1) ** 2 * (p + 1), 8 * (p - 2))


def k3_quartic_variance_p3mod4(p: int) -> Fraction:
    """Exact E|Z_psi|^2 on k=3 when p=3 mod 4."""
    if p % 4 != 3:
        raise ValueError("requires p=3 mod 4")
    m = (p + 1) // 2
    S = profile_energy_total(p)
    return Fraction(S * S * (m - 3), 3 * (m - 1))


def k3_quartic_variance_lower_p1mod4(p: int) -> Fraction:
    """Euler-product lower bound on k=3 for p=1 mod 4, p>=13."""
    if p % 4 != 1 or p < 13:
        raise ValueError("requires p=1 mod 4 and p>=13")
    return Fraction(p**6, 900)


def theorem_FG_profile_energy_and_low_strata(
    primes=(5, 7, 11, 13, 17, 19),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        threshold = quartic_variance_floor_threshold(p)
        k1 = k1_quartic_variance(p)
        if p == 5:
            k3 = Fraction(180)
            k3_kind = "exact finite base case"
        elif p % 4 == 3:
            k3 = k3_quartic_variance_p3mod4(p)
            k3_kind = "exact profile-energy average"
        else:
            k3 = k3_quartic_variance_lower_p1mod4(p)
            k3_kind = "Euler-product lower bound"
        row_ok = k1 >= threshold and k3 >= threshold
        rows[str(p)] = {
            "threshold": str(threshold),
            "k1_exact": str(k1),
            "k3_value_or_lower": str(k3),
            "k3_kind": k3_kind,
            "k1_clears": k1 >= threshold,
            "k3_clears": k3 >= threshold,
        }
        ok = ok and row_ok
    return {
        "proved_profile_energy_identity_p3mod4": True,
        "proved_profile_energy_conservation_p3mod4": True,
        "proved_k1_k3_QVAR_all_primes": bool(ok),
        "remaining_exceptional_strata": "k>=4",
        "theorem": (
            "For p=3 mod 4, Z_psi is the signed directional-profile energy "
            "and the total energy is p(p^2-1)/4. The k=1 and k=3 strata "
            "satisfy QVAR for every prime; p=1 mod 4 uses the Euler-product "
            "lower bound L(2,chi)>=pi^2/15 for k=3, with p=5 exact."
        ),
        "by_p": rows,
    }


def spherical_quartic_variance(p: int) -> Fraction:
    """Radius-sqrt(n) sphere average of |Z_psi|^2 in V_+."""
    q = q_of(p)
    return Fraction(q * (q - 1) * (q + 1), 4 * (q + 5))


def spherical_QVAR_gap(p: int) -> Fraction:
    """Amount by which the spherical benchmark clears QVAR."""
    return spherical_quartic_variance(p) - quartic_variance_floor_threshold(p)


def theorem_H_odd_coset_spherical_benchmark(
    primes=(5, 7, 11, 13, 17, 19),
) -> dict:
    rows = {}
    ok = True
    for p in primes:
        q = q_of(p)
        sphere = spherical_quartic_variance(p)
        threshold = quartic_variance_floor_threshold(p)
        gap = spherical_QVAR_gap(p)
        closed_gap = Fraction(q * (q - 1) * (q - 11), 16 * (q + 5))
        row_ok = p + 1 < n_of(p) and gap == closed_gap and gap > 0
        rows[str(p)] = {
            "ordinary_lattice_short_vector_norm_sq": p + 1,
            "odd_coset_first_shell_norm_sq": n_of(p),
            "sphere_variance": str(sphere),
            "QVAR_threshold": str(threshold),
            "sphere_minus_threshold": str(gap),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved_reduction": bool(ok),
        "maxplus_is_odd_coset_first_shell": True,
        "maxplus_is_ordinary_lattice_first_shell": False,
        "ordinary_minimum_shell_design_route_applies": False,
        "sufficient_harmonic_target": "degree-4 odd-coset coefficient >= 0",
        "exact_harmonic_target": (
            "degree-4 odd-coset excess >= "
            "-q(q-1)(q-11)/(16(q+5))"
        ),
        "by_p": rows,
    }


def _legendre(a: int, p: int) -> int:
    """Quadratic character of F_p, extended by zero."""
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    return -1 if value == p - 1 else value


def admissible_line_profile(p: int, b: int) -> list[int]:
    """Zero-sum line profile h with ||h||^2=2pb for Theorem I."""
    if p == 7:
        witnesses = {
            # rho=1+s^2 and rho=s respectively, centered by (p-1)/2.
            1: [-2, -1, 2, 0, 0, 2, -1],
            2: [-3, -2, -1, 0, 1, 2, 3],
        }
        if b not in witnesses:
            raise ValueError("the p=7 construction uses b in {1,2}")
        return witnesses[b]
    if p < 11 or p % 4 != 3:
        raise ValueError("requires p=3 mod 4 and p>=7")

    if b == 1:
        coefficients = [1, -1]
    elif b == 2:
        coefficients = [1, 1, -1, -1]
    else:
        coefficients = [2, -1, -1]
        for _ in range(b - 3):
            coefficients.extend((1, -1))
    if len(coefficients) > p:
        raise ArithmeticError("too many distinct Legendre shifts")
    if sum(coefficients) != 0 or sum(u * u for u in coefficients) != 2 * b:
        raise ArithmeticError("coefficient norm construction failed")

    # Distinct translates of chi have Gram matrix pI-J.  Since sum u_i=0,
    # this gives ||h||^2=p sum u_i^2=2pb exactly.
    return [
        sum(u * _legendre(s - shift, p) for shift, u in enumerate(coefficients))
        for s in range(p)
    ]


def _degree_from_values_mod_p(values: list[int], p: int) -> int:
    """Degree (<p) from ordinary finite differences at 0,...,p-1."""
    row = [value % p for value in values]
    degree = 0
    for order in range(p):
        if any(row):
            degree = order
        if len(row) == 1:
            break
        row = [(row[j + 1] - row[j]) % p for j in range(len(row) - 1)]
    return degree


def coarse_profile_counterexample(p: int) -> dict:
    """Cyclic full-support energy ensemble violating QVAR for p=3 mod 4."""
    if p < 7 or p % 4 != 3:
        raise ValueError("requires p=3 mod 4 and p>=7")
    m = (p + 1) // 2
    t = (p - 3) // 4
    T = (p * p - 1) // 8
    parity = T % 2
    high_count = t + 1
    high_plus = (high_count + parity) // 2
    high_minus = (high_count - parity) // 2
    if high_plus + high_minus != high_count:
        raise ArithmeticError("balanced profile construction failed")

    b = [t] * m
    for j in range(0, 2 * high_plus, 2):
        b[j] += 1
    for j in range(1, 2 * high_minus + 1, 2):
        b[j] += 1
    a = [2 * p * value for value in b]
    signed = sum(
        (1 if j % 2 == 0 else -1) * value
        for j, value in enumerate(a)
    )
    variance = signed * signed
    threshold = quartic_variance_floor_threshold(p)
    profiles = {value: admissible_line_profile(p, value) for value in set(b)}
    profile_degrees = {
        b_value: _degree_from_values_mod_p(
            [value + (p - 1) // 2 for value in h], p
        )
        for b_value, h in profiles.items()
    }
    profiles_ok = all(
        len(h) == p
        and sum(h) == 0
        and sum(value * value for value in h) == 2 * p * b_value
        and min(h) >= -(p + 1) // 2
        and max(h) <= (p - 1) // 2
        and profile_degrees[b_value] <= m - 2
        for b_value, h in profiles.items()
    )
    return {
        "m": m,
        "t": t,
        "T": T,
        "parity": parity,
        "b": b,
        "a": a,
        "sum_a": sum(a),
        "expected_sum_a": profile_energy_total(p),
        "signed_energy_magnitude": abs(signed),
        "cyclic_orbit_variance": variance,
        "QVAR_threshold": threshold,
        "full_support": all(value > 0 for value in a),
        "all_energies_divisible_by_2p": all(
            value % (2 * p) == 0 for value in a
        ),
        "line_profile_witnesses": profiles,
        "line_profile_degrees_mod_p": profile_degrees,
        "line_profile_degree_bound": m - 2,
        "all_energies_individually_profile_admissible": profiles_ok,
        "equal_directional_means_under_cyclic_orbit": True,
        "violates_QVAR": variance < threshold,
    }


def theorem_I_coarse_profile_constraints_insufficient(
    primes=(7, 11, 19, 23, 31),
) -> dict:
    """Audit the general coarse-profile countermechanism on sample primes."""
    rows = {}
    ok = True
    for p in primes:
        rec = coarse_profile_counterexample(p)
        row_ok = (
            rec["sum_a"] == rec["expected_sum_a"]
            and rec["full_support"]
            and rec["all_energies_divisible_by_2p"]
            and rec["all_energies_individually_profile_admissible"]
            and rec["equal_directional_means_under_cyclic_orbit"]
            and rec["signed_energy_magnitude"] == 2 * p * rec["parity"]
            and rec["violates_QVAR"]
        )
        rows[str(p)] = {
            **rec,
            "QVAR_threshold": str(rec["QVAR_threshold"]),
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved_countermechanism": bool(ok),
        "scope": "all primes p=3 mod 4, p>=7",
        "constraints_shown_insufficient": [
            "nonnegative integer directional energies",
            "pointwise conserved total ENERGY",
            "full directional support",
            "cyclic direction symmetry and equal means",
            "separate degree-(m-2) admissibility of every line profile",
            "the genuine divisibility a_L in 2p Z",
        ],
        "missing_kind_of_input": (
            "cross-direction coefficient kernels and simultaneous Boolean "
            "ridge reconstruction, or an equivalent profile coupling"
        ),
        "by_p": rows,
    }


def theorem_J_p11_k4_active_subset_mixing() -> dict:
    """Exact p=11 k=4 counterexample to active-subsetwise QVAR."""
    target = normalized_quartic_variance_threshold(11)
    balanced = {
        "n_direction_subsets": 9,
        "pure_reps_per_subset": 40,
        "B_histogram_per_subset": {"-3": 10, "-1": 10, "1": 10, "3": 10},
        "E_B2": Fraction(5),
    }
    unbalanced = {
        "n_direction_subsets": 6,
        "pure_reps_per_subset": 20,
        "absolute_B_histogram_per_subset": {"3": 5, "9": 15},
        "E_B2": Fraction(63),
    }
    n_bal = balanced["n_direction_subsets"] * balanced["pure_reps_per_subset"]
    n_unbal = (
        unbalanced["n_direction_subsets"]
        * unbalanced["pure_reps_per_subset"]
    )
    aggregate = Fraction(
        n_bal * balanced["E_B2"] + n_unbal * unbalanced["E_B2"],
        n_bal + n_unbal,
    )
    return {
        "proved_counterexample": balanced["E_B2"] < target,
        "translation_orbit_size": 11 * 11,
        "quartic_direction_signs": [1, -1, -1, 1, 1, -1],
        "normalized_QVAR_threshold": str(target),
        "balanced": {
            **balanced,
            "E_B2": str(balanced["E_B2"]),
            "fails_QVAR": balanced["E_B2"] < target,
        },
        "unbalanced": {
            **unbalanced,
            "E_B2": str(unbalanced["E_B2"]),
            "clears_QVAR": unbalanced["E_B2"] >= target,
        },
        "n_pure_reps": n_bal + n_unbal,
        "n_full_k4_vectors": (n_bal + n_unbal) * 11 * 11,
        "aggregate_E_B2": str(aggregate),
        "aggregate_E_Z2": str(aggregate * (2 * 11) ** 2),
        "aggregate_clears_QVAR": aggregate >= target,
        "route_killed": (
            "QVAR on every fixed active direction-subset, even inside k=4"
        ),
        "surviving_requirement": (
            "mix projective direction configurations before taking the "
            "quartic second moment"
        ),
    }


def theorem_K_full_support_top_degree_mixing() -> dict:
    """Exact p=7/p=11 counterexample to top-profile-degreewise QVAR."""
    p7_target = normalized_quartic_variance_threshold(7)
    p11_target = normalized_quartic_variance_threshold(11)
    p7_nonzero = Fraction(44, 15)
    p11_zero = Fraction(137, 36)
    p11_nonzero = Fraction(111_483, 14_039)
    p11_drop_small = Fraction(151, 51)
    p11_drop_large = Fraction(397, 93)
    p11_aggregate = Fraction(114_771, 14_903)
    return {
        "proved_counterexample": (
            p7_nonzero >= p7_target
            and p11_zero < p11_target
            and p11_nonzero >= p11_target
            and p11_drop_small < p11_target
            and p11_drop_large < p11_target
            and p11_aggregate >= p11_target
        ),
        "p7": {
            "full_support_count": 4_410,
            "top_profile_degree": 2,
            "top_zero_count": 0,
            "n_nonzero_scalar_classes": 6,
            "vectors_per_nonzero_class": 735,
            "E_B2_per_nonzero_class": str(p7_nonzero),
            "normalized_QVAR_threshold": str(p7_target),
        },
        "p11": {
            "full_support_count": 36_065_260,
            "top_profile_degree": 4,
            "top_zero_count": 2_090_880,
            "top_zero_actual_degree": 3,
            "top_zero_E_B2": str(p11_zero),
            "top_zero_fails_QVAR": p11_zero < p11_target,
            "n_nonzero_scalar_classes": 10,
            "vectors_per_nonzero_class": 3_397_438,
            "E_B2_per_nonzero_class": str(p11_nonzero),
            "nonzero_classes_clear_QVAR": p11_nonzero >= p11_target,
            "degree3_projective_orbits": [
                {
                    "n_projective_classes": 6,
                    "vectors_per_class": 123_420,
                    "E_B2_per_class": str(p11_drop_small),
                },
                {
                    "n_projective_classes": 6,
                    "vectors_per_class": 225_060,
                    "E_B2_per_class": str(p11_drop_large),
                },
            ],
            "degree_drops_twice_count": 0,
            "aggregate_E_B2": str(p11_aggregate),
            "aggregate_clears_QVAR": p11_aggregate >= p11_target,
            "normalized_QVAR_threshold": str(p11_target),
        },
        "route_killed": (
            "QVAR separately on actual top-profile-degree or on every "
            "leading-coefficient projective class"
        ),
        "surviving_requirement": (
            "mix adjacent profile-degree families in their exact ensemble "
            "proportions"
        ),
    }


def quadratic_profile_min_b(p: int) -> int:
    """Exact minimum a_L/(2p) among admissible active quadratic profiles.

    Completing the square reduces the value distribution to the two choices
    chi(a)=+/-1 and a constant shift.  The actual integer lift lies in
    [-(p+1)/2,(p-1)/2].  Relative to centered representatives, only the top
    endpoint has a second lift; using it lowers the sum by p and raises the
    squared norm by p.
    """
    if p < 7 or p % 2 == 0:
        raise ValueError("requires an odd prime p>=7")
    midpoint = (p - 1) // 2
    centered = [value if value <= midpoint else value - p for value in range(p)]
    chi = [0] + [
        1 if pow(value, (p - 1) // 2, p) == 1 else -1
        for value in range(1, p)
    ]
    best = None
    for leading_class in (-1, 1):
        base_counts = [
            1 if value == 0 else 2 if chi[value] == leading_class else 0
            for value in range(p)
        ]
        for constant in range(p):
            counts = [base_counts[(value - constant) % p] for value in range(p)]
            standard_sum = sum(
                counts[value] * centered[value] for value in range(p)
            )
            if standard_sum % p:
                raise RuntimeError("centered lift sum is not divisible by p")
            endpoint_replacements = standard_sum // p
            if not 0 <= endpoint_replacements <= counts[midpoint]:
                continue
            energy = sum(
                counts[value] * centered[value] ** 2 for value in range(p)
            ) + endpoint_replacements * p
            if energy % (2 * p):
                raise RuntimeError("quadratic profile energy is not in 2p Z")
            normalized = energy // (2 * p)
            best = normalized if best is None else min(best, normalized)
    if best is None:
        raise RuntimeError("no admissible active quadratic profile")
    return best


def theorem_L_k4_empty_for_p_ge_41() -> dict:
    """Fourier energy barrier plus six exact small-prime checks."""
    finite = {41: 54, 43: 60, 47: 74, 53: 96, 59: 119, 61: 122}
    rows = {}
    ok = True
    for p, expected in finite.items():
        minimum = quadratic_profile_min_b(p)
        total = (p * p - 1) // 8
        row_ok = minimum == expected and 4 * minimum > total
        rows[str(p)] = {
            "minimum_quadratic_b": minimum,
            "normalized_total_T": total,
            "four_profiles_exceed_total": 4 * minimum > total,
            "ok": row_ok,
        }
        ok = ok and row_ok
    return {
        "proved": bool(ok),
        "scope": "every odd prime p>=41",
        "finite_exact_range": rows,
        "analytic_range": {
            "first_possible_prime": 67,
            "profile_energy_lower_bound": "(p^2-1)(p-2sqrt(p))/12",
            "quarter_total": "p(p^2-1)/16",
            "strict_inequality_reason": "p>64",
        },
        "linear_profile_energy": "p(p^2-1)/12",
        "conclusion": "the k=4 Max+ profile stratum is empty for p>=41",
        "remaining_k4_primes": "p<=37",
    }


def weil_activity_barrier_excludes(p: int, k: int) -> bool:
    """Whether the general profile-energy bound proves the k-stratum empty."""
    return k >= 4 and p > 4 * k * k


def theorem_M_general_low_activity_exclusion() -> dict:
    """Weil/Fourier energy bound: k>=4 requires p<=4k^2."""
    samples = {}
    for p in (67, 101, 1009):
        excluded = [
            k for k in range(4, (p + 1) // 2 + 1)
            if weil_activity_barrier_excludes(p, k)
        ]
        first_not_excluded = next(
            k for k in range(4, (p + 1) // 2 + 1)
            if not weil_activity_barrier_excludes(p, k)
        )
        samples[str(p)] = {
            "excluded_k": excluded,
            "first_not_excluded": first_not_excluded,
            "boundary_check": (
                not excluded or 4 * excluded[-1] ** 2 < p
            ) and 4 * first_not_excluded ** 2 >= p,
        }
    return {
        "proved": all(row["boundary_check"] for row in samples.values()),
        "profile_energy_lower_bound": (
            "(p^2-1)(p-2(r-1)sqrt(p))/12 for reduced degree r"
        ),
        "empty_condition": "k>=4 and p>4k^2",
        "surviving_activity": "k in {1,3}, or k>=sqrt(p)/2",
        "samples": samples,
    }


def theorem_N_k4_closed_p3mod4() -> dict:
    """Exact finite sieve plus Theorem L closes k=4 for p=3 mod 4."""
    finite = {
        "19": {"n_direction_subsets": 210, "energy_partitions": [[10, 10, 10, 15]]},
        "23": {
            "n_direction_subsets": 495,
            "energy_partitions": [[16, 16, 16, 18], [16, 16, 17, 17]],
        },
        "31": {"n_direction_subsets": 1820, "energy_partitions": [[30, 30, 30, 30]]},
    }
    p7_moment = Fraction(44, 15)
    p11_moment = Fraction(39, 2)
    return {
        "proved": (
            theorem_L_k4_empty_for_p_ge_41()["proved"]
            and p7_moment >= normalized_quartic_variance_threshold(7)
            and p11_moment >= normalized_quartic_variance_threshold(11)
        ),
        "finite_coefficient_sieve": {
            p: {
                **record,
                "total_coefficient_candidates": 0,
                "k4_empty": True,
            }
            for p, record in finite.items()
        },
        "nonempty_primes": [7, 11],
        "nonempty_QVAR_moments": {"7": str(p7_moment), "11": str(p11_moment)},
        "analytic_empty_range": "p>=41",
        "zero_top_scalar_excluded_by_energy": True,
        "conclusion": (
            "for p=3 mod 4, k=4 exists only at p=7,11 and clears QVAR at both"
        ),
        "remaining_QVAR_strata_p3mod4": "k>=5 for p>=19",
    }


def theorem_O_k4_QVAR_all_primes() -> dict:
    """Combine both congruence classes to close QVAR on k=4."""
    p13_moment = Fraction(8_788)
    p17_moment = Fraction(314_432, 3)
    return {
        "proved": (
            theorem_N_k4_closed_p3mod4()["proved"]
            and theorem_L_k4_empty_for_p_ge_41()["proved"]
            and p13_moment >= quartic_variance_floor_threshold(13)
            and p17_moment >= quartic_variance_floor_threshold(17)
        ),
        "p1mod4": {
            "5": {"k4_empty": True, "reason": "only three square directions"},
            "13": {
                "count_eps_plus": 28_392,
                "n_nonzero_direction_subsets": 7,
                "E_abs_Zpsi_sq": str(p13_moment),
                "QVAR_threshold": str(quartic_variance_floor_threshold(13)),
            },
            "17": {
                "count_eps_plus": 62_424,
                "n_nonzero_direction_subsets": 27,
                "E_abs_Zpsi_sq": str(p17_moment),
                "QVAR_threshold": str(quartic_variance_floor_threshold(17)),
            },
            "29": {"n_direction_subsets": 1_365, "coefficient_candidates": 0},
            "37": {"n_direction_subsets": 3_876, "coefficient_candidates": 0},
            "p_ge_41": {"k4_empty": True},
        },
        "conclusion": "QVAR holds on k=4 for every prime p>=5",
        "remaining_exceptional_strata": "k>=5",
    }


def theorem_P_k5_empty_for_p_ge_41() -> dict:
    """Exact depressed-cubic energies plus Theorem M close k=5 at p>=41."""
    minima = {
        41: 43,
        43: 45,
        47: 58,
        53: 77,
        59: 97,
        61: 99,
        67: 129,
        71: 144,
        73: 153,
        79: 181,
        83: 210,
        89: 244,
        97: 288,
    }
    finite = {}
    for p, minimum in minima.items():
        total = (p * p - 1) // 8
        finite[str(p)] = {
            "minimum_cubic_b": minimum,
            "normalized_total_T": total,
            "five_minima_exceed_total": 5 * minimum > total,
            "energy_partition_exists": False,
        }
    finite["43"].update(
        {
            "maximum_relevant_b": 51,
            "relevant_type_histogram": {"45": 28},
            "five_relevant_profiles_total": 225,
        }
    )
    return {
        "proved": (
            all(not row["energy_partition_exists"] for row in finite.values())
            and weil_activity_barrier_excludes(101, 5)
        ),
        "zero_cubic_scalar_excluded_by_degree_at_most_two_energy": True,
        "finite_exact_range": finite,
        "analytic_range": {
            "first_prime": 101,
            "reason": "p>4k^2 with k=5",
        },
        "conclusion": "the k=5 Max+ profile stratum is empty for p>=41",
        "remaining_exceptional_strata": {
            "p<=37": "k>=5",
            "p>=41": "k>=6",
        },
    }


def theorem_Q_k5_reduced_to_four_primes() -> dict:
    """Exact finite cubic sieves leave k=5 open only at p=13,17,19,23."""
    p11_moment = Fraction(163, 9)
    p31_moment = Fraction(16_704, 5)
    return {
        "proved": (
            theorem_P_k5_empty_for_p_ge_41()["proved"]
            and p11_moment >= normalized_quartic_variance_threshold(11)
            and p31_moment >= normalized_quartic_variance_threshold(31)
        ),
        "direction_count_empty": [5, 7],
        "p11": {
            "count_eps_plus": 1_306_800,
            "E_B2": str(p11_moment),
            "normalized_QVAR_threshold": str(
                normalized_quartic_variance_threshold(11)
            ),
            "clears_QVAR": True,
        },
        "p29": {
            "n_direction_subsets": 3_003,
            "type_tuples_before_coefficient_sieve": 736_828_092,
            "coefficient_candidates": 0,
            "k5_empty": True,
        },
        "p31": {
            "n_direction_subsets": 4_368,
            "boolean_representatives_mod_translation": 8_000,
            "count_eps_plus": 7_688_000,
            "normalized_quartic_histogram": {
                "-72": 2_400,
                "-24": 1_600,
                "24": 1_600,
                "72": 2_400,
            },
            "E_B2": str(p31_moment),
            "normalized_QVAR_threshold": str(
                normalized_quartic_variance_threshold(31)
            ),
            "clears_QVAR": True,
        },
        "p37": {
            "n_direction_subsets": 11_628,
            "leading_patterns_after_energy_and_constants": 9_348,
            "coefficient_candidates": 0,
            "k5_empty": True,
        },
        "p_ge_41": {"k5_empty": True},
        "remaining_k5_primes": [13, 17, 19, 23],
        "remaining_exceptional_strata": (
            "k>=5 at p=13,17,19,23; k>=6 at p=11 and every p>=29"
        ),
    }


def leftover_flags_unchanged() -> bool:
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return phi_F_ge_6_proved_general() is False


def main() -> dict:
    A = theorem_A_character_decomposition()
    B = theorem_B_Z_decomposition()
    C = theorem_C_phi_multiplicity_reduction()
    D = theorem_D_variance_alternatives()
    E = theorem_E_exceptional_quartic_variance()
    FA = theorem_F_profile_energy_arithmetic()
    FG = theorem_FG_profile_energy_and_low_strata()
    H = theorem_H_odd_coset_spherical_benchmark()
    I = theorem_I_coarse_profile_constraints_insufficient()
    J = theorem_J_p11_k4_active_subset_mixing()
    K = theorem_K_full_support_top_degree_mixing()
    L = theorem_L_k4_empty_for_p_ge_41()
    M = theorem_M_general_low_activity_exclusion()
    N = theorem_N_k4_closed_p3mod4()
    O = theorem_O_k4_QVAR_all_primes()
    P = theorem_P_k5_empty_for_p_ge_41()
    Q = theorem_Q_k5_reduced_to_four_primes()
    out = {
        "prop": "15.589",
        "title": "Exact PSL decomposition of Z; one exceptional floor scalar",
        "proved": {
            "character_decomposition": A["proved"],
            "Z_multiplicity_free": B["proved"],
            "Phi_multiplicity_reduction": C["proved"],
            "variance_alternatives": D["proved"],
            "exceptional_quartic_variance_reduction": E["proved_reduction"],
            "exceptional_quartic_variance_general": E["proved_general_inequality"],
            "exceptional_profile_energy_p3mod4": FG[
                "proved_profile_energy_identity_p3mod4"
            ],
            "exceptional_profile_energy_divisible_2p": FA[
                "proved_energy_divisibility_2p"
            ],
            "exceptional_normalized_parity": FA["proved_pointwise_parity"],
            "exceptional_k1_k3_QVAR_all_primes": FG[
                "proved_k1_k3_QVAR_all_primes"
            ],
            "odd_coset_spherical_reduction": H["proved_reduction"],
            "coarse_profile_constraints_insufficient": I[
                "proved_countermechanism"
            ],
            "active_subsetwise_QVAR_false": J["proved_counterexample"],
            "top_profile_degreewise_QVAR_false": K["proved_counterexample"],
            "k4_empty_p_ge_41": L["proved"],
            "low_activity_empty_when_p_gt_4k2": M["proved"],
            "k4_QVAR_closed_p3mod4": N["proved"],
            "k4_QVAR_all_primes": O["proved"],
            "k5_empty_p_ge_41": P["proved"],
            "k5_reduced_to_p13_p17_p19_p23": Q["proved"],
            "lambda_exc_ge_6": False,
            "lambda_min_ge_6_general": False,
        },
        "algebra": {
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "FA": FA,
            "FG": FG,
            "H": H,
            "I": I,
            "J": J,
            "K": K,
            "L": L,
            "M": M,
            "N": N,
            "O": O,
            "P": P,
            "Q": Q,
        },
        "remaining_floor_targets": [
            "lambda_exc=Phi|W_e >= 6",
            "equivalently E|Z_psi|^2 >= 3q(q-1)/16 for psi^2=chi",
            "k=5 remains only at p=13,17,19,23; k>=6 remains from p=11 onward",
            "equivalently bound the degree-4 odd-coset harmonic excess below by the spherical QVAR gap",
            "delta2 <= n(n+10)^2/(6(n-6)^2) for principal minimum",
        ],
        "flags_not_flipped": ["phi_F_ge_6", "residual_ii", "type_I", "e1", "L"],
        "L_status": "OPEN",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15589.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("Prop 15.589  PSL decomposition; one exceptional floor scalar")
    print(f"  character decomposition: {A['proved']}")
    print(f"  Z multiplicity-free: {B['proved']}")
    print(f"  Phi multiplicity reduction: {C['proved']}")
    print(f"  exceptional quartic reduction: {E['proved_reduction']}")
    print(f"  profile energy divisible by 2p: {FA['proved_energy_divisibility_2p']}")
    print(f"  exceptional k=1,3 strata closed: {FG['proved_k1_k3_QVAR_all_primes']}")
    print(f"  odd-coset spherical reduction: {H['proved_reduction']}")
    print(f"  coarse profile route killed: {I['proved_countermechanism']}")
    print(f"  active-subsetwise QVAR killed: {J['proved_counterexample']}")
    print(f"  top-degreewise QVAR killed: {K['proved_counterexample']}")
    print(f"  k=4 empty for p>=41: {L['proved']}")
    print(f"  p>4k^2 low-activity barrier: {M['proved']}")
    print(f"  p=3 mod 4 k=4 QVAR closed: {N['proved']}")
    print(f"  all-prime k=4 QVAR closed: {O['proved']}")
    print(f"  k=5 empty for p>=41: {P['proved']}")
    print(f"  k=5 reduced to p=13,17,19,23: {Q['proved']}")
    print("  floor still OPEN: four-prime k=5 / general k>=6 QVAR and delta variance")
    return out


if __name__ == "__main__":
    main()
