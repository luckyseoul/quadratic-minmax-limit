/*
 * Exhaustive enumeration of switching-normalised symmetric zero-diagonal +-1 matrices.
 *
 * Phi(A) = max_x |sum_{i<j} A_ij x_i x_j| is invariant under A -> DAD for a diagonal
 * +-1 matrix D, so every matrix can be normalised to have first row and column +1.  The
 * program runs over all 2^((n-1)(n-2)/2) normalised matrices of order n and prints, for
 * each value v <= T, the number of normalised matrices with Phi = v.  Running with
 * T = m_n shows that no matrix has Phi < m_n (no smaller value is printed) and gives the
 * number of normalised minimisers.
 *
 * Build:  gcc -O2 -fopenmp -o normalized_exhaustive_phi_counts normalized_exhaustive_phi_counts.c
 * Run:    ./normalized_exhaustive_phi_counts n T
 *         n=5 T=4, n=6 T=5, n=7 T=9, n=8 T=10 take under a second; n=9 T=12 takes about
 *         15 seconds on 8 cores (2^28 matrices, early abort once |Q| exceeds T).
 * Expected output (number of normalised matrices with Phi = m_n):
 *         n=5: 12   n=6: 12   n=7: 3240   n=8: 4200   n=9: 1607760
 * Order n <= 11 is supported by the array sizes; only the C standard library is used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: %s n T\n", argv[0]); return 1; }
  int n = atoi(argv[1]), T = atoi(argv[2]);
  if (n < 3 || n > 11) { fprintf(stderr, "3 <= n <= 11\n"); return 1; }
  int pi[64], pj[64], m = 0;
  for (int i = 1; i < n; i++)
    for (int j = i + 1; j < n; j++) { pi[m] = i; pj[m] = j; m++; }
  uint64_t tot = 1ULL << m;
  long long cnt[256] = {0};
#pragma omp parallel
  {
    long long loc[256] = {0};
#pragma omp for schedule(dynamic, 4096)
    for (uint64_t b = 0; b < tot; b++) {
      int A[12][12];
      for (int j = 1; j < n; j++) A[0][j] = A[j][0] = 1;
      for (int i = 0; i < n; i++) A[i][i] = 0;
      for (int e = 0; e < m; e++) {
        int v = ((b >> e) & 1) ? -1 : 1;
        A[pi[e]][pj[e]] = A[pj[e]][pi[e]] = v;
      }
      int x[12], h[12], Q = 0;
      for (int i = 0; i < n; i++) {
        x[i] = 1;
        h[i] = 0;
        for (int j = 0; j < n; j++) h[i] += A[i][j];
      }
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) Q += A[i][j];
      int mx = abs(Q), ok = (mx <= T);
      /* classes {x,-x} with x_0 = +1, Gray order over x_1..x_{n-1} */
      uint32_t tt = 1u << (n - 1);
      for (uint32_t g = 1; g < tt && ok; g++) {
        int k = __builtin_ctz(g) + 1;
        Q -= 2 * x[k] * h[k];
        for (int i = 0; i < n; i++)
          if (i != k) h[i] -= 2 * A[i][k] * x[k];
        x[k] = -x[k];
        int a = abs(Q);
        if (a > mx) mx = a;
        if (mx > T) ok = 0;
      }
      if (ok) loc[mx]++;
    }
#pragma omp critical
    for (int i = 0; i < 256; i++) cnt[i] += loc[i];
  }
  for (int i = 0; i < 256; i++)
    if (cnt[i]) printf("n=%d Phi=%d normalised matrices: %lld\n", n, i, cnt[i]);
  return 0;
}
