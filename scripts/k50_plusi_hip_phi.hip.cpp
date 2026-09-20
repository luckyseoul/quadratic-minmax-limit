// Exact Phi of plus-I lift K=[[B,B+I],[B+I,-B]] on nuka RX 9070 XT (gfx1201).
// Projective cube: x[0]=+1, |x|=m, |y|=m, n=2m. Q = Q_B(x)-Q_B(y)+x^T(B+I)y.
// hipBLAS GEMM + two-stage absmax. ROCm 10: hipcc --offload-arch=gfx1201.
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

#define HIP_CHECK(call)                                                        \
  do {                                                                         \
    hipError_t e = (call);                                                     \
    if (e != hipSuccess) {                                                     \
      fprintf(stderr, "HIP %s:%d %s\n", __FILE__, __LINE__,                    \
              hipGetErrorString(e));                                           \
      return 1;                                                                \
    }                                                                          \
  } while (0)

#define BLAS_CHECK(call)                                                       \
  do {                                                                         \
    hipblasStatus_t s = (call);                                                \
    if (s != HIPBLAS_STATUS_SUCCESS) {                                         \
      fprintf(stderr, "hipBLAS %s:%d status %d\n", __FILE__, __LINE__, (int)s);\
      return 1;                                                                \
    }                                                                          \
  } while (0)

// X: NX x M, x[0]=+1, NX=1<<(M-1). Column-major M x NX with lda=M.
__global__ void fill_X(float *X, int M, int NX) {
  int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (idx >= NX)
    return;
  X[0 + (size_t)idx * M] = 1.f;
  for (int c = 1; c < M; ++c)
    X[c + (size_t)idx * M] = ((idx >> (c - 1)) & 1) ? -1.f : 1.f;
}

// Y panel: Py columns starting at global y0. Column-major M x Py.
__global__ void fill_Y_panel(float *Y, int M, int Py, int y0) {
  int i = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (i >= Py)
    return;
  int idx = y0 + i;
  for (int c = 0; c < M; ++c)
    Y[c + (size_t)i * M] = ((idx >> c) & 1) ? -1.f : 1.f;
}

// QY[idx] = 0.5 y^T B y from bits; B is column-major M x M.
__global__ void fill_QY(const float *B, float *QY, int M, int NY) {
  int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (idx >= NY)
    return;
  float s = 0.f;
  for (int i = 0; i < M; ++i) {
    float yi = ((idx >> i) & 1) ? -1.f : 1.f;
    for (int j = i + 1; j < M; ++j) {
      float yj = ((idx >> j) & 1) ? -1.f : 1.f;
      s += B[i + j * M] * yi * yj;
    }
  }
  QY[idx] = s;
}

// Qvec[i] = 0.5 * sum_c X[c,i] * (B X)[c,i]; BX is M x N, same layout.
__global__ void quad_from_BX(const float *X, const float *BX, float *Q, int M,
                             int N) {
  int i = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  if (i >= N)
    return;
  float s = 0.f;
  for (int c = 0; c < M; ++c)
    s += X[c + (size_t)i * M] * BX[c + (size_t)i * M];
  Q[i] = 0.5f * s;
}

// Q is Py x Px column-major (ldc=Py): Q[j + c*Py]. Threads along j coalesce.
__global__ void absmax_coalesced(const float *Q, const float *QX, const float *QY,
                                 float *best, int Px, int Py) {
  int j = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  float v = 0.f;
  if (j < Py) {
    float qy = QY[j];
    for (int c = 0; c < Px; ++c) {
      float q = Q[j + (size_t)c * Py] + QX[c] - qy;
      v = fmaxf(v, fabsf(q));
    }
  }
  __shared__ float sh[256];
  sh[threadIdx.x] = v;
  __syncthreads();
  for (int s = 128; s > 0; s >>= 1) {
    if ((int)threadIdx.x < s)
      sh[threadIdx.x] = fmaxf(sh[threadIdx.x], sh[threadIdx.x + s]);
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    unsigned int *t = (unsigned int *)best;
    unsigned int assumed, old = *t;
    do {
      assumed = old;
      float nv = fmaxf(__uint_as_float(assumed), sh[0]);
      old = atomicCAS(t, assumed, __float_as_uint(nv));
    } while (old != assumed);
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: k50_plusi_hip_phi B.txt [chunk]\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  std::vector<float> Bh;
  float z;
  while (in >> z)
    Bh.push_back(z);
  int M = (int)std::lround(std::sqrt((double)Bh.size()));
  if (M * M != (int)Bh.size() || M < 3 || (M % 2) == 0) {
    fprintf(stderr, "B must be odd-order square, got %zu entries\n", Bh.size());
    return 2;
  }
  int PX = (argc > 2) ? atoi(argv[2]) : 256;
  if (PX < 1)
    PX = 256;
  int PY = (argc > 3) ? atoi(argv[3]) : (1 << 20);
  if (PY < 1)
    PY = 1 << 20;

  const int NX = 1 << (M - 1);
  const int NY = 1 << M;
  if (PX > NX)
    PX = NX;
  if (PY > NY)
    PY = NY;
  printf("{\"m\":%d,\"n\":%d,\"NX\":%d,\"NY\":%d,\"Px\":%d,\"Py\":%d}\n", M,
         2 * M, NX, NY, PX, PY);
  fflush(stdout);

  std::vector<float> Bcm(M * M), Dcm(M * M);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < M; ++j) {
      float v = Bh[i * M + j];
      Bcm[i + j * M] = v;
      Dcm[i + j * M] = v + (i == j ? 1.f : 0.f);
    }

  float *dB = nullptr, *dD = nullptr, *dX = nullptr, *dBX = nullptr;
  float *dQX = nullptr, *dQY = nullptr, *dY = nullptr, *dM2 = nullptr;
  float *dQ = nullptr, *dBest = nullptr;
  HIP_CHECK(hipMalloc(&dB, sizeof(float) * M * M));
  HIP_CHECK(hipMalloc(&dD, sizeof(float) * M * M));
  HIP_CHECK(hipMemcpy(dB, Bcm.data(), sizeof(float) * M * M, hipMemcpyHostToDevice));
  HIP_CHECK(hipMemcpy(dD, Dcm.data(), sizeof(float) * M * M, hipMemcpyHostToDevice));
  HIP_CHECK(hipMalloc(&dX, sizeof(float) * (size_t)M * NX));
  HIP_CHECK(hipMalloc(&dBX, sizeof(float) * (size_t)M * NX));
  HIP_CHECK(hipMalloc(&dQX, sizeof(float) * NX));
  HIP_CHECK(hipMalloc(&dQY, sizeof(float) * NY));
  HIP_CHECK(hipMalloc(&dY, sizeof(float) * (size_t)M * PY));
  HIP_CHECK(hipMalloc(&dM2, sizeof(float) * (size_t)M * PY));
  HIP_CHECK(hipMalloc(&dQ, sizeof(float) * (size_t)PX * PY));
  HIP_CHECK(hipMalloc(&dBest, sizeof(float)));
  HIP_CHECK(hipMemset(dBest, 0, sizeof(float)));
  size_t free_b = 0, total_b = 0;
  HIP_CHECK(hipMemGetInfo(&free_b, &total_b));
  printf("gpu_mem_ready free_MiB=%.0f total_MiB=%.0f\n", free_b / 1048576.0,
         total_b / 1048576.0);
  fflush(stdout);

  fill_X<<<(NX + 255) / 256, 256>>>(dX, M, NX);
  fill_QY<<<(NY + 255) / 256, 256>>>(dB, dQY, M, NY);
  HIP_CHECK(hipDeviceSynchronize());

  hipblasHandle_t h;
  BLAS_CHECK(hipblasCreate(&h));
  const float alpha = 1.f, beta = 0.f;
  BLAS_CHECK(hipblasSgemm(h, HIPBLAS_OP_N, HIPBLAS_OP_N, M, NX, M, &alpha, dB,
                          M, dX, M, &beta, dBX, M));
  quad_from_BX<<<(NX + 255) / 256, 256>>>(dX, dBX, dQX, M, NX);
  HIP_CHECK(hipDeviceSynchronize());
  (void)hipFree(dBX);
  dBX = nullptr;

  hipEvent_t t0, t1;
  HIP_CHECK(hipEventCreate(&t0));
  HIP_CHECK(hipEventCreate(&t1));
  HIP_CHECK(hipEventRecord(t0));

  const int nxc = (NX + PX - 1) / PX;
  const int nyc = (NY + PY - 1) / PY;
  const int ntile = nxc * nyc;
  float best = 0.f;
  int tile = 0;
  for (int y0 = 0; y0 < NY; y0 += PY) {
    int Py = std::min(PY, NY - y0);
    fill_Y_panel<<<(Py + 255) / 256, 256>>>(dY, M, Py, y0);
    BLAS_CHECK(hipblasSgemm(h, HIPBLAS_OP_N, HIPBLAS_OP_N, M, Py, M, &alpha, dD,
                            M, dY, M, &beta, dM2, M));
    for (int x0 = 0; x0 < NX; x0 += PX) {
      int Px = std::min(PX, NX - x0);
      // Q = M2^T * Xpanel : Py x Px, col-major ldc=Py (coalesced in j)
      BLAS_CHECK(hipblasSgemm(h, HIPBLAS_OP_T, HIPBLAS_OP_N, Py, Px, M, &alpha,
                              dM2, M, dX + (size_t)x0 * M, M, &beta, dQ, Py));
      absmax_coalesced<<<(Py + 255) / 256, 256>>>(dQ, dQX + x0, dQY + y0, dBest,
                                                  Px, Py);
      ++tile;
      if ((tile % 2048) == 0 || tile == ntile) {
        HIP_CHECK(hipMemcpy(&best, dBest, sizeof(float), hipMemcpyDeviceToHost));
        float ms = 0.f;
        HIP_CHECK(hipEventRecord(t1));
        HIP_CHECK(hipEventSynchronize(t1));
        HIP_CHECK(hipEventElapsedTime(&ms, t0, t1));
        printf("progress %d/%d best=%.1f elapsed_s=%.1f\n", tile, ntile, best,
               ms / 1000.f);
        fflush(stdout);
      }
    }
  }

  HIP_CHECK(hipMemcpy(&best, dBest, sizeof(float), hipMemcpyDeviceToHost));
  float ms = 0.f;
  HIP_CHECK(hipEventRecord(t1));
  HIP_CHECK(hipEventSynchronize(t1));
  HIP_CHECK(hipEventElapsedTime(&ms, t0, t1));
  printf("{\"phi\":%.1f,\"paley50\":175,\"m\":%d,\"n\":%d,\"wall_s\":%.3f,"
         "\"backend\":\"hipblas-gfx1201\"}\n",
         best, M, 2 * M, ms / 1000.f);
  fflush(stdout);

  hipblasDestroy(h);
  (void)hipFree(dB);
  (void)hipFree(dD);
  (void)hipFree(dX);
  (void)hipFree(dQX);
  (void)hipFree(dQY);
  (void)hipFree(dY);
  (void)hipFree(dM2);
  (void)hipFree(dQ);
  (void)hipFree(dBest);
  return 0;
}
