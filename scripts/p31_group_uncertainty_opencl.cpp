// Targeted p=31 OpenCL hunt for the grouped block-incidence inequality.
//
// This is deliberately not a prime sweep.  It searches for a nonzero set S
// of antipodal point classes for which more projective direction groups of
// M^T 1_S vanish than |S|.  Such a witness would disprove the proposed
// all-prime lemma
//
//     #{A : (M^T 1_S)|_{B_A}=0} <= |S|.
//
// Odd support is already handled by a symbolic injection argument.  Every
// chain here therefore starts even and uses two-bit moves, so GPU time is
// spent only on the genuinely open even-support case.  Any reported support
// is reconstructed and checked independently on the host before it is printed.

#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int P = 31;
constexpr int N = (P * P - 1) / 2;
constexpr int D = P + 1;
constexpr int H = (P - 1) / 2;
constexpr int WORDS = (N + 63) / 64;

struct Pair {
  int x;
  int y;
};

int modp(int x) {
  x %= P;
  return x < 0 ? x + P : x;
}

std::vector<Pair> antipodal_representatives() {
  std::vector<Pair> out;
  for (int x = 0; x < P; ++x) {
    for (int y = 0; y < P; ++y) {
      if (x == 0 && y == 0) continue;
      const Pair neg{modp(-x), modp(-y)};
      if (x < neg.x || (x == neg.x && y < neg.y)) out.push_back({x, y});
    }
  }
  if (static_cast<int>(out.size()) != N) {
    throw std::runtime_error("antipodal representative count changed");
  }
  return out;
}

std::vector<Pair> projective_functionals() {
  std::vector<Pair> out;
  for (int slope = 0; slope < P; ++slope) out.push_back({1, slope});
  out.push_back({0, 1});
  return out;
}

std::vector<int8_t> point_direction_block_map(
    const std::vector<Pair>& points, const std::vector<Pair>& directions) {
  std::array<int, P> square_index{};
  square_index.fill(-1);
  std::vector<int> squares;
  for (int value = 1; value < P; ++value) {
    const int square = modp(value * value);
    if (square_index[square] == -1) {
      square_index[square] = -2;
      squares.push_back(square);
    }
  }
  std::sort(squares.begin(), squares.end());
  if (static_cast<int>(squares.size()) != H) {
    throw std::runtime_error("nonzero-square count changed");
  }
  for (int i = 0; i < H; ++i) square_index[squares[i]] = i;

  std::vector<int8_t> map(N * D, -1);
  for (int point = 0; point < N; ++point) {
    for (int direction = 0; direction < D; ++direction) {
      const int value = modp(directions[direction].x * points[point].x +
                             directions[direction].y * points[point].y);
      if (value != 0) {
        map[point * D + direction] =
            static_cast<int8_t>(square_index[modp(value * value)]);
      }
    }
  }
  return map;
}

void check(cl_int code, const char* where) {
  if (code != CL_SUCCESS) {
    std::ostringstream message;
    message << where << " failed with OpenCL error " << code;
    throw std::runtime_error(message.str());
  }
}

std::string json_escape(const std::string& value) {
  std::ostringstream out;
  for (char c : value) {
    if (c == '\\' || c == '"') out << '\\';
    out << c;
  }
  return out.str();
}

const char* KERNEL = R"CLC(
#define N 480
#define D 32
#define H 15
#define WORDS 8

uint next_u32(uint *state) {
  uint x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x ? x : 0x9e3779b9U;
  return *state;
}

kernel void hunt(
    global const char *point_map,
    const uint seed,
    const uint iterations,
    global int *best_margin_out,
    global ushort *best_weight_out,
    global ushort *best_silent_out,
    global ulong *best_state_out) {
  const uint gid = get_global_id(0);
  uint rng = seed ^ (0x9e3779b9U * (gid + 1U));
  ushort masks[D];
  ulong state[WORDS];
  ulong best_state[WORDS];
  for (int a = 0; a < D; ++a) masks[a] = 0;
  for (int w = 0; w < WORDS; ++w) state[w] = 0;

  ushort weight = 0;
  const uint mode = gid & 3U;
  if (mode == 0U || mode == 1U) {
    int first = (int)(next_u32(&rng) % N);
    int second = (int)(next_u32(&rng) % (N - 1));
    if (second >= first) ++second;
    state[first >> 6] ^= ((ulong)1) << (first & 63);
    state[second >> 6] ^= ((ulong)1) << (second & 63);
    weight = 2;
    if (mode == 1U) {
      // Add another pair; repeated points merely cancel and preserve parity.
      for (int j = 0; j < 2; ++j) {
        const int point = (int)(next_u32(&rng) % N);
        const ulong bit = ((ulong)1) << (point & 63);
        const bool was_set = (state[point >> 6] & bit) != 0;
        state[point >> 6] ^= bit;
        weight = (ushort)(weight + (was_set ? -1 : 1));
      }
      if (weight == 0) {
        state[0] = 1;
        state[1] = 1;
        weight = 2;
      }
    }
  } else if (mode == 2U) {
    const int direction = (int)(next_u32(&rng) % D);
    const int block = (int)(next_u32(&rng) % H);
    for (int point = 0; point < N; ++point) {
      if ((int)point_map[point * D + direction] == block) {
        state[point >> 6] ^= ((ulong)1) << (point & 63);
        ++weight;
      }
    }
    // Toggle one point so the odd p-point block becomes an even state.
    const int point = (int)(next_u32(&rng) % N);
    const ulong bit = ((ulong)1) << (point & 63);
    const bool was_set = (state[point >> 6] & bit) != 0;
    state[point >> 6] ^= bit;
    weight = (ushort)(weight + (was_set ? -1 : 1));
  } else {
    const int count = 2 + 2 * (int)(next_u32(&rng) % 15U);
    for (int j = 0; j < count; ++j) {
      const int point = (int)(next_u32(&rng) % N);
      const ulong bit = ((ulong)1) << (point & 63);
      const bool was_set = (state[point >> 6] & bit) != 0;
      state[point >> 6] ^= bit;
      weight = (ushort)(weight + (was_set ? -1 : 1));
    }
    if (weight == 0) {
      state[0] = 1;
      state[1] = 1;
      weight = 2;
    }
  }

  for (int point = 0; point < N; ++point) {
    if ((state[point >> 6] & (((ulong)1) << (point & 63))) == 0) continue;
    for (int a = 0; a < D; ++a) {
      const int block = (int)point_map[point * D + a];
      if (block >= 0) masks[a] ^= (ushort)(1U << block);
    }
  }
  ushort silent = 0;
  for (int a = 0; a < D; ++a) silent += masks[a] == 0;
  int best_margin = (int)silent - (int)weight;
  ushort best_weight = weight;
  ushort best_silent = silent;
  for (int w = 0; w < WORDS; ++w) best_state[w] = state[w];

  for (uint iteration = 0; iteration < iterations; ++iteration) {
    const int first = (int)(next_u32(&rng) % N);
    int second = (int)(next_u32(&rng) % (N - 1));
    if (second >= first) ++second;
    const ulong first_bit = ((ulong)1) << (first & 63);
    const ulong second_bit = ((ulong)1) << (second & 63);
    const bool removing_first =
        (state[first >> 6] & first_bit) != 0;
    const bool removing_second =
        (state[second >> 6] & second_bit) != 0;
    const ushort new_weight = (ushort)(
        weight + (removing_first ? -1 : 1) +
        (removing_second ? -1 : 1));
    if (new_weight == 0) continue;

    ushort new_silent = silent;
    for (int a = 0; a < D; ++a) {
      const ushort old_mask = masks[a];
      ushort new_mask = old_mask;
      const int first_block = (int)point_map[first * D + a];
      const int second_block = (int)point_map[second * D + a];
      if (first_block >= 0) new_mask ^= (ushort)(1U << first_block);
      if (second_block >= 0) new_mask ^= (ushort)(1U << second_block);
      if (old_mask == 0 && new_mask != 0) --new_silent;
      if (old_mask != 0 && new_mask == 0) ++new_silent;
    }
    const int old_score = (int)weight - (int)silent;
    const int new_score = (int)new_weight - (int)new_silent;
    const int uphill = new_score - old_score;
    bool accept = uphill <= 0;
    if (!accept) {
      const uint phase = iteration & 2047U;
      const uint shift = (uint)min(15, 3 + uphill + (phase > 1535U ? 0 : 2));
      accept = (next_u32(&rng) & ((1U << shift) - 1U)) == 0U;
    }
    if (!accept) continue;

    state[first >> 6] ^= first_bit;
    state[second >> 6] ^= second_bit;
    for (int a = 0; a < D; ++a) {
      const int first_block = (int)point_map[first * D + a];
      const int second_block = (int)point_map[second * D + a];
      if (first_block >= 0) masks[a] ^= (ushort)(1U << first_block);
      if (second_block >= 0) masks[a] ^= (ushort)(1U << second_block);
    }
    weight = new_weight;
    silent = new_silent;
    const int margin = (int)silent - (int)weight;
    if (margin > best_margin ||
        (margin == best_margin && weight > best_weight)) {
      best_margin = margin;
      best_weight = weight;
      best_silent = silent;
      for (int w = 0; w < WORDS; ++w) best_state[w] = state[w];
    }
  }

  best_margin_out[gid] = best_margin;
  best_weight_out[gid] = best_weight;
  best_silent_out[gid] = best_silent;
  for (int w = 0; w < WORDS; ++w) {
    best_state_out[gid * WORDS + w] = best_state[w];
  }
}
)CLC";

}  // namespace

int main(int argc, char** argv) try {
  const size_t chains = argc > 1 ? std::stoull(argv[1]) : 8192;
  const uint32_t iterations = argc > 2 ? std::stoul(argv[2]) : 10000;
  const uint32_t seed = argc > 3 ? std::stoul(argv[3]) : 157521U;
  if (chains == 0 || iterations == 0) {
    throw std::runtime_error("chains and iterations must be positive");
  }

  const auto points = antipodal_representatives();
  const auto directions = projective_functionals();
  const auto map = point_direction_block_map(points, directions);

  cl_uint platform_count = 0;
  check(clGetPlatformIDs(0, nullptr, &platform_count), "clGetPlatformIDs(count)");
  if (platform_count == 0) throw std::runtime_error("no OpenCL platform");
  std::vector<cl_platform_id> platforms(platform_count);
  check(clGetPlatformIDs(platform_count, platforms.data(), nullptr),
        "clGetPlatformIDs(list)");

  cl_device_id device = nullptr;
  for (cl_platform_id platform : platforms) {
    if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr) ==
        CL_SUCCESS) break;
    device = nullptr;
  }
  if (!device) throw std::runtime_error("no OpenCL GPU device");

  size_t device_name_size = 0;
  check(clGetDeviceInfo(device, CL_DEVICE_NAME, 0, nullptr, &device_name_size),
        "clGetDeviceInfo(size)");
  std::string device_name(device_name_size, '\0');
  check(clGetDeviceInfo(device, CL_DEVICE_NAME, device_name.size(),
                        device_name.data(), nullptr),
        "clGetDeviceInfo(name)");
  while (!device_name.empty() && device_name.back() == '\0') device_name.pop_back();

  cl_int error = CL_SUCCESS;
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &error);
  check(error, "clCreateContext");
  cl_command_queue queue = clCreateCommandQueue(context, device, 0, &error);
  check(error, "clCreateCommandQueue");
  const size_t source_size = std::strlen(KERNEL);
  cl_program program = clCreateProgramWithSource(context, 1, &KERNEL, &source_size, &error);
  check(error, "clCreateProgramWithSource");
  error = clBuildProgram(program, 1, &device, "-cl-std=CL1.2", nullptr, nullptr);
  if (error != CL_SUCCESS) {
    size_t log_size = 0;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
    std::string log(log_size, '\0');
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log.size(), log.data(), nullptr);
    throw std::runtime_error("OpenCL build failed: " + log);
  }
  cl_kernel kernel = clCreateKernel(program, "hunt", &error);
  check(error, "clCreateKernel");

  std::vector<int> margins(chains);
  std::vector<uint16_t> weights(chains), silents(chains);
  std::vector<uint64_t> states(chains * WORDS);
  cl_mem map_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     map.size(), const_cast<int8_t*>(map.data()), &error);
  check(error, "clCreateBuffer(map)");
  cl_mem margin_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                        margins.size() * sizeof(int), nullptr, &error);
  check(error, "clCreateBuffer(margins)");
  cl_mem weight_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                        weights.size() * sizeof(uint16_t), nullptr, &error);
  check(error, "clCreateBuffer(weights)");
  cl_mem silent_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                        silents.size() * sizeof(uint16_t), nullptr, &error);
  check(error, "clCreateBuffer(silents)");
  cl_mem state_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                       states.size() * sizeof(uint64_t), nullptr, &error);
  check(error, "clCreateBuffer(states)");

  check(clSetKernelArg(kernel, 0, sizeof(map_buffer), &map_buffer), "clSetKernelArg(0)");
  check(clSetKernelArg(kernel, 1, sizeof(seed), &seed), "clSetKernelArg(1)");
  check(clSetKernelArg(kernel, 2, sizeof(iterations), &iterations), "clSetKernelArg(2)");
  check(clSetKernelArg(kernel, 3, sizeof(margin_buffer), &margin_buffer), "clSetKernelArg(3)");
  check(clSetKernelArg(kernel, 4, sizeof(weight_buffer), &weight_buffer), "clSetKernelArg(4)");
  check(clSetKernelArg(kernel, 5, sizeof(silent_buffer), &silent_buffer), "clSetKernelArg(5)");
  check(clSetKernelArg(kernel, 6, sizeof(state_buffer), &state_buffer), "clSetKernelArg(6)");

  const size_t global = chains;
  check(clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel");
  check(clFinish(queue), "clFinish");
  check(clEnqueueReadBuffer(queue, margin_buffer, CL_TRUE, 0,
                            margins.size() * sizeof(int), margins.data(), 0, nullptr, nullptr),
        "clEnqueueReadBuffer(margins)");
  check(clEnqueueReadBuffer(queue, weight_buffer, CL_TRUE, 0,
                            weights.size() * sizeof(uint16_t), weights.data(), 0, nullptr, nullptr),
        "clEnqueueReadBuffer(weights)");
  check(clEnqueueReadBuffer(queue, silent_buffer, CL_TRUE, 0,
                            silents.size() * sizeof(uint16_t), silents.data(), 0, nullptr, nullptr),
        "clEnqueueReadBuffer(silents)");
  check(clEnqueueReadBuffer(queue, state_buffer, CL_TRUE, 0,
                            states.size() * sizeof(uint64_t), states.data(), 0, nullptr, nullptr),
        "clEnqueueReadBuffer(states)");

  size_t best = 0;
  for (size_t i = 1; i < chains; ++i) {
    if (margins[i] > margins[best] ||
        (margins[i] == margins[best] && weights[i] > weights[best])) best = i;
  }

  std::array<uint16_t, D> host_masks{};
  std::vector<int> support;
  for (int point = 0; point < N; ++point) {
    if ((states[best * WORDS + (point >> 6)] &
         (uint64_t{1} << (point & 63))) == 0) continue;
    support.push_back(point);
    for (int a = 0; a < D; ++a) {
      const int block = map[point * D + a];
      if (block >= 0) host_masks[a] ^= uint16_t{1} << block;
    }
  }
  std::vector<int> silent_directions;
  for (int a = 0; a < D; ++a) {
    if (host_masks[a] == 0) silent_directions.push_back(a);
  }
  const int verified_margin =
      static_cast<int>(silent_directions.size()) - static_cast<int>(support.size());
  const bool verified = support.size() == weights[best] &&
                        silent_directions.size() == silents[best] &&
                        verified_margin == margins[best] && !support.empty();
  const bool even_support = support.size() % 2 == 0;

  std::cout << "{\"p\":" << P
            << ",\"device\":\"" << json_escape(device_name) << "\""
            << ",\"chains\":" << chains
            << ",\"iterations\":" << iterations
            << ",\"seed\":" << seed
            << ",\"scope\":\"even support only\""
            << ",\"best_margin\":" << verified_margin
            << ",\"weight\":" << support.size()
            << ",\"silent_group_count\":" << silent_directions.size()
            << ",\"device_reported_margin\":" << margins[best]
            << ",\"device_reported_weight\":" << weights[best]
            << ",\"device_reported_silent\":" << silents[best]
            << ",\"counterexample\":" << (verified_margin > 0 ? "true" : "false")
            << ",\"verified\":" << (verified && even_support ? "true" : "false")
            << ",\"support_indices\":[";
  for (size_t i = 0; i < support.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << support[i];
  }
  std::cout << "],\"silent_direction_indices\":[";
  for (size_t i = 0; i < silent_directions.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << silent_directions[i];
  }
  std::cout << "]}\n";

  clReleaseMemObject(state_buffer);
  clReleaseMemObject(silent_buffer);
  clReleaseMemObject(weight_buffer);
  clReleaseMemObject(margin_buffer);
  clReleaseMemObject(map_buffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
  return verified && even_support ? 0 : 2;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
