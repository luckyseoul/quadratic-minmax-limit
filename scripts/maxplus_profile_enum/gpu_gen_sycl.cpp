// SYCL emit_cands — same algorithm as gpu_gen_device.py CUDA/HIP kernel.
// Build on jellyfin: source /opt/intel/oneapi/setvars.sh
//   icpx -fsycl -O3 -shared -fPIC -o gpu_gen_sycl.so gpu_gen_sycl.cpp
#include <sycl/sycl.hpp>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

static sycl::device pick_gpu() {
    for (auto const& d : sycl::device::get_devices(sycl::info::device_type::gpu)) {
        auto name = d.get_info<sycl::info::device::name>();
        if (name.find("Arc") != std::string::npos)
            return d;
    }
    return sycl::device{sycl::gpu_selector_v};
}

static sycl::queue& qref() {
    static sycl::queue q{pick_gpu()};
    return q;
}

extern "C" const char* sycl_device_name() {
    static std::string name =
        qref().get_device().get_info<sycl::info::device::name>();
    return name.c_str();
}

extern "C" int64_t emit_chunk_sycl(
    int p, int k, int u_lo, int nci, int c0, int thi, int tlo, int npr,
    int two_p, uint64_t ncombo,
    const int32_t* av, const int32_t* af, const int32_t* an,
    const uint8_t* aull, const int32_t* UU, const int32_t* vmap,
    const int16_t* probes, const int16_t* cprobes,
    int64_t* codes, int64_t* fsums, int64_t cap
) {
    if (nci <= 0 || k < 2 || k > 8 || p <= 0)
        return 0;
    try {
        auto& q = qref();
        const int pp = p * p;
        const size_t nUU = (size_t)(u_lo + nci);  // UU indexed by absolute ci
        const size_t avn = (size_t)k * (size_t)p * (size_t)p;
        const size_t ann = (size_t)k * (size_t)p;
        const size_t uun = nUU * (size_t)k;
        const size_t prn = (size_t)k * (size_t)pp * (size_t)npr;

        int32_t *d_av = sycl::malloc_device<int32_t>(avn, q);
        int32_t *d_af = sycl::malloc_device<int32_t>(avn, q);
        int32_t *d_an = sycl::malloc_device<int32_t>(ann, q);
        uint8_t *d_aull = sycl::malloc_device<uint8_t>(ann, q);
        int32_t *d_UU = sycl::malloc_device<int32_t>(uun, q);
        int32_t *d_vmap = sycl::malloc_device<int32_t>(avn, q);
        int16_t *d_pr = sycl::malloc_device<int16_t>(prn, q);
        int16_t *d_cpr = sycl::malloc_device<int16_t>(prn, q);
        int64_t *d_codes = sycl::malloc_device<int64_t>((size_t)cap, q);
        int64_t *d_fsums = sycl::malloc_device<int64_t>((size_t)cap, q);
        uint64_t *d_ctr = sycl::malloc_device<uint64_t>(1, q);

        q.memcpy(d_av, av, avn * sizeof(int32_t));
        q.memcpy(d_af, af, avn * sizeof(int32_t));
        q.memcpy(d_an, an, ann * sizeof(int32_t));
        q.memcpy(d_aull, aull, ann * sizeof(uint8_t));
        q.memcpy(d_UU, UU, uun * sizeof(int32_t));
        q.memcpy(d_vmap, vmap, avn * sizeof(int32_t));
        q.memcpy(d_pr, probes, prn * sizeof(int16_t));
        q.memcpy(d_cpr, cprobes, prn * sizeof(int16_t));
        q.memset(d_ctr, 0, sizeof(uint64_t));
        q.wait();

        const uint64_t ntot = (uint64_t)nci * ncombo;
        const size_t tpb = 256;
        size_t nblocks = (size_t)((ntot + tpb - 1) / tpb);
        if (nblocks > 262144)
            nblocks = 262144;
        if (nblocks < 1)
            nblocks = 1;

        q.parallel_for(
             sycl::nd_range<1>(sycl::range<1>(nblocks * tpb), sycl::range<1>(tpb)),
             [=](sycl::nd_item<1> it) {
                 const uint64_t gid = (uint64_t)it.get_global_linear_id();
                 const uint64_t stride = (uint64_t)it.get_global_range(0);
                 int idx[8];
                 for (uint64_t t = gid; t < ntot; t += stride) {
                     const int ci_off = (int)(t / ncombo);
                     uint64_t combo = t - (uint64_t)ci_off * ncombo;
                     const int ci = u_lo + ci_off;
                     int ok = 1;
                     for (int j = 0; j < k; ++j) {
                         const int u = d_UU[ci * k + j];
                         if (!d_aull[j * p + u]) {
                             ok = 0;
                             break;
                         }
                     }
                     if (!ok)
                         continue;
                     int lens_ok = 1;
                     for (int j = 0; j < k - 1; ++j) {
                         idx[j] = (int)(combo % (uint64_t)p);
                         combo /= (uint64_t)p;
                         const int u = d_UU[ci * k + j];
                         const int ln = d_an[j * p + u];
                         if (idx[j] >= ln) {
                             lens_ok = 0;
                             break;
                         }
                     }
                     if (!lens_ok)
                         continue;
                     int vs = 0, fs = 0;
                     for (int j = 0; j < k - 1; ++j) {
                         const int u = d_UU[ci * k + j];
                         const int base = (j * p + u) * p + idx[j];
                         vs += d_av[base];
                         fs += d_af[base];
                     }
                     int vlast = c0 - vs;
                     vlast %= p;
                     if (vlast < 0)
                         vlast += p;
                     const int ul = d_UU[ci * k + (k - 1)];
                     const int pos = d_vmap[((k - 1) * p + ul) * p + vlast];
                     if (pos < 0)
                         continue;
                     const int ftot = fs + d_af[((k - 1) * p + ul) * p + pos];
                     int good = 1;
                     for (int a = 0; a < npr; ++a) {
                         int ps = 0, cv = 0;
                         for (int j = 0; j < k; ++j) {
                             const int u = d_UU[ci * k + j];
                             const int v = (j < k - 1)
                                               ? d_av[(j * p + u) * p + idx[j]]
                                               : vlast;
                             const int uv = u * p + v;
                             const int po = (j * pp + uv) * npr + a;
                             ps += (int)d_pr[po];
                             cv += (int)d_cpr[po];
                         }
                         if (ftot == 0) {
                             if (ps != thi && ps != tlo) {
                                 good = 0;
                                 break;
                             }
                         } else {
                             const int d = ps - thi;
                             int g = d / two_p;
                             if (d < 0 && (d % two_p) != 0)
                                 g -= 1;
                             if (g < -1 || g > ftot || g > cv) {
                                 good = 0;
                                 break;
                             }
                         }
                     }
                     if (!good)
                         continue;
                     int64_t code = (int64_t)ci * 16777216LL;
                     int64_t place = 1;
                     for (int j = 0; j < k - 1; ++j) {
                         code += (int64_t)idx[j] * place;
                         place *= 16;
                     }
                     code += (int64_t)pos * place;
                     sycl::atomic_ref<uint64_t, sycl::memory_order::relaxed,
                                      sycl::memory_scope::device,
                                      sycl::access::address_space::global_space>
                         ctr_ref(*d_ctr);
                     const uint64_t slot = ctr_ref.fetch_add(1);
                     if ((int64_t)slot < cap) {
                         sycl::atomic_ref<int64_t, sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>
                             c_ref(d_codes[slot]);
                         sycl::atomic_ref<int64_t, sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>
                             f_ref(d_fsums[slot]);
                         c_ref.store(code);
                         f_ref.store((int64_t)ftot);
                     }
                 }
             })
            .wait();

        uint64_t nout = 0;
        q.memcpy(&nout, d_ctr, sizeof(uint64_t)).wait();
        if ((int64_t)nout > cap) {
            sycl::free(d_av, q);
            sycl::free(d_af, q);
            sycl::free(d_an, q);
            sycl::free(d_aull, q);
            sycl::free(d_UU, q);
            sycl::free(d_vmap, q);
            sycl::free(d_pr, q);
            sycl::free(d_cpr, q);
            sycl::free(d_codes, q);
            sycl::free(d_fsums, q);
            sycl::free(d_ctr, q);
            return (int64_t)nout;
        }
        if (nout) {
            q.memcpy(codes, d_codes, nout * sizeof(int64_t));
            q.memcpy(fsums, d_fsums, nout * sizeof(int64_t));
            q.wait();
        }
        sycl::free(d_av, q);
        sycl::free(d_af, q);
        sycl::free(d_an, q);
        sycl::free(d_aull, q);
        sycl::free(d_UU, q);
        sycl::free(d_vmap, q);
        sycl::free(d_pr, q);
        sycl::free(d_cpr, q);
        sycl::free(d_codes, q);
        sycl::free(d_fsums, q);
        sycl::free(d_ctr, q);
        return (int64_t)nout;
    } catch (sycl::exception const& e) {
        std::fprintf(stderr, "sycl: %s\n", e.what());
        return -1;
    }
}

// Host NumpyTester is 90%+ of a380 wall. Same gather+filter as GpuTester,
// on the Arc so emit is not GPU-theater.

static inline int floor_div_pos(int n, int d) {
    int q = n / d;
    int r = n % d;
    if (r != 0 && n < 0)
        q -= 1;
    return q;
}

struct TestCtx {
    int p = 0, k = 0, qn = 0, thi = 0, tlo = 0, two_p = 0, nUU = 0;
    int32_t *d_UU = nullptr;
    int32_t *d_av = nullptr;
    int32_t *d_af = nullptr;
    int32_t *d_Tm = nullptr;
    int16_t *d_cont = nullptr;
    int8_t *d_cov = nullptr;
};

static TestCtx g_test;

static void test_release() {
    auto& q = qref();
    auto rel = [&](auto*& p) {
        if (p) {
            sycl::free(p, q);
            p = nullptr;
        }
    };
    rel(g_test.d_UU);
    rel(g_test.d_av);
    rel(g_test.d_af);
    rel(g_test.d_Tm);
    rel(g_test.d_cont);
    rel(g_test.d_cov);
}

extern "C" int sycl_test_load(
    int p, int k, int qn, int thi, int tlo, int nUU,
    const int32_t* UU, const int32_t* av, const int32_t* af, const int32_t* Tm,
    const int16_t* cont, const int8_t* cov
) {
    if (p <= 0 || k < 2 || k > 8 || qn <= 0 || qn > 256 || nUU <= 0 || k * p > 128)
        return -2;
    try {
        auto& q = qref();
        test_release();
        const int pp = p * p;
        g_test.p = p;
        g_test.k = k;
        g_test.qn = qn;
        g_test.thi = thi;
        g_test.tlo = tlo;
        g_test.two_p = 2 * p;
        g_test.nUU = nUU;
        const size_t uun = (size_t)nUU * (size_t)k;
        const size_t avn = (size_t)k * (size_t)p * (size_t)p;
        const size_t tmn = (size_t)k * (size_t)qn;
        const size_t cn = (size_t)k * (size_t)pp * (size_t)qn;
        g_test.d_UU = sycl::malloc_device<int32_t>(uun, q);
        g_test.d_av = sycl::malloc_device<int32_t>(avn, q);
        g_test.d_af = sycl::malloc_device<int32_t>(avn, q);
        g_test.d_Tm = sycl::malloc_device<int32_t>(tmn, q);
        g_test.d_cont = sycl::malloc_device<int16_t>(cn, q);
        g_test.d_cov = sycl::malloc_device<int8_t>(cn, q);
        q.memcpy(g_test.d_UU, UU, uun * sizeof(int32_t));
        q.memcpy(g_test.d_av, av, avn * sizeof(int32_t));
        q.memcpy(g_test.d_af, af, avn * sizeof(int32_t));
        q.memcpy(g_test.d_Tm, Tm, tmn * sizeof(int32_t));
        q.memcpy(g_test.d_cont, cont, cn * sizeof(int16_t));
        q.memcpy(g_test.d_cov, cov, cn * sizeof(int8_t));
        q.wait();
        return 0;
    } catch (sycl::exception const& e) {
        std::fprintf(stderr, "sycl_test_load: %s\n", e.what());
        return -1;
    }
}

extern "C" int sycl_test_batch(
    int64_t B,
    const int64_t* codes, const int64_t* fsums,
    int32_t* f0_idx, int32_t* fl_idx,
    int32_t* n_f0, int32_t* n_fl
) {
    if (B <= 0) {
        *n_f0 = 0;
        *n_fl = 0;
        return 0;
    }
    if (!g_test.d_UU || B > 4000000)
        return -2;
    try {
        auto& q = qref();
        const int p = g_test.p, k = g_test.k, qn = g_test.qn;
        const int thi = g_test.thi, tlo = g_test.tlo, two_p = g_test.two_p;
        const int pp = p * p;
        const int nUU = g_test.nUU;
        int32_t *d_UU = g_test.d_UU, *d_av = g_test.d_av, *d_af = g_test.d_af;
        int32_t *d_Tm = g_test.d_Tm;
        int16_t *d_cont = g_test.d_cont;
        int8_t *d_cov = g_test.d_cov;

        int64_t *d_codes = sycl::malloc_device<int64_t>((size_t)B, q);
        int64_t *d_fsums = sycl::malloc_device<int64_t>((size_t)B, q);
        int32_t *d_f0 = sycl::malloc_device<int32_t>((size_t)B, q);
        int32_t *d_fl = sycl::malloc_device<int32_t>((size_t)B, q);
        uint32_t *d_nf0 = sycl::malloc_device<uint32_t>(1, q);
        uint32_t *d_nfl = sycl::malloc_device<uint32_t>(1, q);
        q.memcpy(d_codes, codes, (size_t)B * sizeof(int64_t));
        q.memcpy(d_fsums, fsums, (size_t)B * sizeof(int64_t));
        q.memset(d_nf0, 0, sizeof(uint32_t));
        q.memset(d_nfl, 0, sizeof(uint32_t));
        q.wait();

        const size_t tpb = 256;
        size_t nblocks = ((size_t)B + tpb - 1) / tpb;
        q.parallel_for(
             sycl::nd_range<1>(sycl::range<1>(nblocks * tpb), sycl::range<1>(tpb)),
             [=](sycl::nd_item<1> it) {
                 const int64_t b = (int64_t)it.get_global_linear_id();
                 if (b >= B)
                     return;
                 const int64_t code = d_codes[b];
                 const int ftot = (int)d_fsums[b];
                 const int ci = (int)(code / 16777216LL);
                 if (ci < 0 || ci >= nUU)
                     return;
                 int vidx[8];
                 int uvj[8];
                 int64_t c = code;
                 for (int j = 0; j < k; ++j) {
                     vidx[j] = (int)(c % 16);
                     c /= 16;
                     const int u = d_UU[ci * k + j];
                     const int vv = d_av[(j * p + u) * p + vidx[j]];
                     uvj[j] = u * p + vv;
                 }
                 if (ftot == 0) {
                     int ok = 1;
                     for (int x = 0; x < qn; ++x) {
                         int s = 0;
                         for (int j = 0; j < k; ++j)
                             s += (int)d_cont[(j * pp + uvj[j]) * qn + x];
                         if (s != thi && s != tlo) {
                             ok = 0;
                             break;
                         }
                     }
                     if (ok) {
                         sycl::atomic_ref<uint32_t, sycl::memory_order::relaxed,
                                          sycl::memory_scope::device,
                                          sycl::access::address_space::global_space>
                             nref(*d_nf0);
                         const uint32_t slot = nref.fetch_add(1);
                         d_f0[slot] = (int32_t)b;
                     }
                     return;
                 }
                 int Gl[128];
                 for (int i = 0; i < k * p; ++i)
                     Gl[i] = 0;
                 int gok = 1;
                 for (int x = 0; x < qn; ++x) {
                     int s = 0, v = 0;
                     for (int j = 0; j < k; ++j) {
                         const int off = (j * pp + uvj[j]) * qn + x;
                         s += (int)d_cont[off];
                         v += (int)d_cov[off];
                     }
                     const int gx = floor_div_pos(s - thi, two_p);
                     if (gx < -1 || gx > ftot || gx > v) {
                         gok = 0;
                         break;
                     }
                     for (int j = 0; j < k; ++j)
                         Gl[j * p + d_Tm[j * qn + x]] += gx;
                 }
                 if (!gok)
                     return;
                 int line_ok = 1;
                 for (int j = 0; j < k && line_ok; ++j) {
                     const int u = d_UU[ci * k + j];
                     const int fj = d_af[(j * p + u) * p + vidx[j]];
                     const int Fj = ftot - fj;
                     for (int s = 0; s < p; ++s) {
                         const int lo = Gl[j * p + s] - Fj;
                         const int hi = Gl[j * p + s] + p - Fj;
                         const int can0 = (lo <= 0 && hi >= 0);
                         const int can1 = (lo <= p && hi >= p);
                         if (!can0 && !can1) {
                             line_ok = 0;
                             break;
                         }
                     }
                 }
                 if (!line_ok)
                     return;
                 sycl::atomic_ref<uint32_t, sycl::memory_order::relaxed,
                                  sycl::memory_scope::device,
                                  sycl::access::address_space::global_space>
                     nref(*d_nfl);
                 const uint32_t slot = nref.fetch_add(1);
                 d_fl[slot] = (int32_t)b;
             })
            .wait();

        uint32_t nf0 = 0, nfl = 0;
        q.memcpy(&nf0, d_nf0, sizeof(uint32_t));
        q.memcpy(&nfl, d_nfl, sizeof(uint32_t));
        q.wait();
        if (nf0)
            q.memcpy(f0_idx, d_f0, nf0 * sizeof(int32_t));
        if (nfl)
            q.memcpy(fl_idx, d_fl, nfl * sizeof(int32_t));
        q.wait();
        *n_f0 = (int32_t)nf0;
        *n_fl = (int32_t)nfl;
        sycl::free(d_codes, q);
        sycl::free(d_fsums, q);
        sycl::free(d_f0, q);
        sycl::free(d_fl, q);
        sycl::free(d_nf0, q);
        sycl::free(d_nfl, q);
        return 0;
    } catch (sycl::exception const& e) {
        std::fprintf(stderr, "sycl_test_batch: %s\n", e.what());
        return -1;
    }
}
