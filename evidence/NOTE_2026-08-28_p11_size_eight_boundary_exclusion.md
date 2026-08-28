# Finite p=11 size-eight boundary exclusion

Date: 2026-08-28. This is Proposition 15.670. It closes the first
all-finite survivor left by Proposition 15.669: every eight-point boundary
in the affine plane over \(\mathbb F_{11}\) violates the exact split
parity-floor budget.

This is not a closure of residual (ii). The separate infinity-plus-nine
profile, all larger \(p=11\) boundaries, the remaining primes, and the full
residual graph constraints remain open.

## 1. Exact boundary condition

Let \(D\subset\mathbb F_{11^2}\) have eight points. For each of the twelve
projective \(\mathbb F_{11}\)-directions \(d\), let \(b_d\) be the number of
parallel affine fibres containing an odd number of points of \(D\). Since
\(|D|\) is even,

\[
 b_d\in\{0,2,4,6,8\}.
\]

There are six directions of each quadratic kernel type
\(\epsilon_d\in\{-1,1\}\). At residual size \(|H|=4p+1=45\), each type has
the exact scaled slack budget

\[
 { (p+1)^2\over2}=72.                              \tag{15.670.1}
\]

For a finite even boundary, the parity phase is
\(\eta_d=\mathbf 1_{\epsilon_d=c_H}\). Proposition 15.632's exact rational
parity-majorant LP gives

\[
\begin{array}{c|rrrrr}
b&0&2&4&6&8\\ \hline
f_0(b)&0&12&16&22&16\\
f_1(b)&22&10&22&18&22.
\end{array}                                         \tag{15.670.2}
\]

Thus a necessary condition for a residual graph with boundary \(D\) is,
for its sign \(c_H\in\{-1,1\}\),

\[
 \sum_{d:\epsilon_d=-1} f_{\mathbf1_{\epsilon_d=c_H}}(b_d)\le72,
 \qquad
 \sum_{d:\epsilon_d=1} f_{\mathbf1_{\epsilon_d=c_H}}(b_d)\le72. \tag{15.670.3}
\]

The census below tests exactly (15.670.3). It is a boundary-level necessary
condition, not a search for edge lifts.

## 2. Lossless affine normalization

Choose an ordered pair \(x\ne y\) in \(D\). The unique affine similarity

\[
 z\longmapsto {z-x\over y-x}                       \tag{15.670.4}
\]

sends that pair to the field points \(0,1\). Consequently it suffices for
exclusion to inspect all

\[
 \binom{119}{6}=3,470,108,187                       \tag{15.670.5}
\]

eight-sets containing \(0,1\), instead of all

\[
 \binom{121}{8}=899,749,479,915                     \tag{15.670.6}
\]

finite eight-sets. The exact pointed-set identity is

\[
 \binom{121}{8}\,8\cdot7
 =\binom{119}{6}\,121\cdot120.                     \tag{15.670.7}
\]

Multiplication by a scalar \(a\ne0\) sends every direction type
\(\epsilon_d\) to \(\chi(a)\epsilon_d\). Transfer
\(c_H\) simultaneously to \(\chi(a)c_H\); then the phase test
\(\epsilon_d=c_H\) is unchanged. The two type budgets are equal, so a
nonsquare type swap changes no feasibility condition. Testing both values
of \(c_H\) on every normalized set is therefore lossless.

The proposition verifier checks this finite-field action directly for all
120 nonzero scalars, all 121 translations, and all twelve directions:
1,440 scalar-direction pairs, 1,452 translation-direction pairs, and 2,880
phase-transfer cases.

## 3. Exact exhaustive census

`scripts/p11_size8_normalized_floor_gpu.py` ranks the six free points in
lexicographic combination order and decodes every rank directly on the GPU.
For each set it forms the twelve odd-fibre masks, computes both exact type
cost pairs from (15.670.2), atomically accumulates their complete
\(133\times133\) histograms, and atomically retains the first minimizer.
Output is written through a temporary file and an atomic rename.

The complete result for each sign is

| \(c_H\) | floor survivors | \(\min_D\max(C_{-},C_{+})\) | first minimizing rank | first type-cost pair |
|---:|---:|---:|---:|---:|
| \(-1\) | 0 | 76 | 108,601,023 | \((76,66)\) |
| \(+1\) | 0 | 76 | 135,919,443 | \((64,76)\) |

The exact budget is 72, so every normalized boundary has at least one type
cost at least 76. The contradiction gap is four:

\[
 \boxed{\min_{D,c_H}\max(C_{-},C_{+})-72=4.}        \tag{15.670.8}
\]

Together with the lossless normalization, this proves

\[
 \boxed{\text{every finite }p=11\text{ size-eight residual boundary is
 impossible}.}                                     \tag{15.670.9}
\]

## 4. Independent replays and checks

The full rank interval was replayed on two GPU stacks and architectures:

| device | stack | checked normalized sets | elapsed | sets/second |
|---|---|---:|---:|---:|
| Tesla V100-SXM2-16GB | CUDA | 3,470,108,187 | 4.852482 s | 715,120,320 |
| AMD Radeon RX 9070 XT | ROCm/HIP | 3,470,108,187 | 1.239626 s | 2,799,317,607 |

Both records reproduce every nonzero entry of both complete cost-pair
histograms. Their semantic little-endian uint64 hashes are

```text
c_H=-1  106de27687d372ce083a3fcfb1adb3fcf50830f14cd314842eebacd9488fbf01
c_H=+1  8df9472aa3d7bb1c24db7240a967d4a861b6e63343e07ea48adf07a473958e39
```

Before each full run, an independent CPU `itertools.combinations` traversal
checked the first 100,000 sets. Every histogram entry, survivor count, and
minimum key matched the GPU result. The theorem verifier also reconstructs
the two minimizing boundaries with the generic rational floor routine and
pins both raw evidence-file hashes.

## 5. Literature and OEIS context check

The closest targeted finite-geometry hit was Ball--Csajbók,
[On sets of points with few odd secants](https://arxiv.org/abs/1711.10876),
which studies lower bounds for odd secants of \(q+2\)-point sets in a
projective plane. Kiermaier--Kurz,
[Maximal integral point sets in affine planes over finite fields](https://arxiv.org/abs/1401.2825),
studies prescribed directions and a Paley-clique connection. Neither source
states the eight-point split-type floor budget (15.670.3), the exact minimum
76, or the exclusion (15.670.9). This is a duplicate/context check, not an
unqualified priority claim.

Individual OEIS searches found 3,470,108,187 in
[A004379](https://oeis.org/A004379) and 899,749,479,915 in the triangle
[A126450](https://oeis.org/A126450). Those occurrences are unrelated
binomial-table values; here the numbers are simply \(\binom{119}{6}\) and
\(\binom{121}{8}\). Searches for the two first-minimizer ranks returned no
relevant sequence context. No sequence submission is proposed.

## 6. Exact scope after Proposition 15.670

At \(p=11\), the first remaining floor-plus-pair profiles are now infinity
plus nine finite points and all-finite boundaries of size at least ten.
They are only relaxed count profiles, not known affine boundaries or known
residual graphs. Proposition 15.670 does not prove general residual (ii),
R1, global QVAR, Type I, or existence of the limit.

## 7. Reproduction

~~~bash
# Generic theorem audit after the two records are present.
python src/e1_gmin_m4_prop15670.py
python -m pytest -q tests/test_prop15670.py

# Soulkiller / Tesla V100 (CuPy CUDA 12.8).
CUDA_PATH=/usr/local/cuda-12.8 \
LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64 \
python scripts/p11_size8_normalized_floor_gpu.py \
  --output evidence/p11_size8_normalized_floor_v100.json

# NUKA / RX 9070 XT (CuPy ROCm 7.2).
/home/nick/.venvs/rocm72/bin/python \
  scripts/p11_size8_normalized_floor_gpu.py \
  --output evidence/p11_size8_normalized_floor_rx9070xt.json
~~~

The compact theorem record is
`evidence/e1_gmin_m4_prop15670.json`.
