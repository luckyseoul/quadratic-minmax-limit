# Exact relative-gauge occupancy: balanced Paley `C18`

This package is a finite diagnostic of the active labeled relative-gauge
composition state. It does **not** prove a cross-order estimate,
good-fiber abundance, or convergence.

## Question tested

Split the Paley symmetric conference matrix of order 18 (prime 17) into two
consecutive 9-vertex blocks. For the resulting internal blocks and rectangular
cross block, form the exact relative-gauge occupancy on

\[
G=P_9\times P_9\times\{\pm1\},\qquad |G|=2^{17}=131072.
\]

The diagnostic uses the first admissible odd energy strictly above the formal
balanced ceiling from the known finite value `m_9=12`:
\[
2^{3/2}m_9=24\sqrt2<35.
\]
The selected principal blocks are not asserted to be order-9 optimizers.
This is therefore a structured next-size calibration, not a multiplier-two
proof.

## Exact result

The local maxima are `(14,14,33)`, giving independent ceiling 61 and deficit
cutoff `61-35=26`.

| quantity | value |
| --- | ---: |
| subthreshold product triples | 85,986,322 |
| included shell triples | 168 |
| empty fibers | 1 |
| least positive occupancy | 140 |
| largest occupancy | 1,203 |

Thus the original `C18` gauge supplies an empty fiber at this finite target,
but it is only one fiber among 131,072. This repeats the sparse-vacancy
behavior in the older `C14` calibration rather than supplying the desired
good-fiber abundance theorem.

The transform uses signed 64-bit integer arithmetic. Exact integrality after
the Walsh inverse was checked, and independent CPU runs on Soulkiller (88
threads) and NUKA (16 threads) returned the same fields in
[results.json](results.json).

## Reproduce

```bash
g++ -O3 -std=c++20 -Wall -Wextra -Wpedantic -Werror -fopenmp \
  scripts/relative_gauge_paley18.cpp -o /tmp/relative_gauge_paley18
OMP_NUM_THREADS=88 /tmp/relative_gauge_paley18
```

This is an explicit relative-gauge fiber calculation, not a Boolean-norm
search. It does not claim that `m_18=33`.
