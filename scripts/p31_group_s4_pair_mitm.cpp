// Exact pair-signature meet-in-the-middle decision procedure for the p=31,
// weight-four grouped uncertainty case.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

constexpr int P = 31;
constexpr int POINTS = 480;
constexpr int DIRECTIONS = 32;
constexpr int FIBRES = 16;
constexpr std::uint64_t EXPECTED_PAIRS = 114960ULL;

struct Point {
    int x;
    int y;
    friend bool operator<(const Point& a, const Point& b) {
        return std::tie(a.x, a.y) < std::tie(b.x, b.y);
    }
};

struct PairRecord {
    std::uint64_t signature;
    std::uint32_t code;
};

int mod(int value) {
    value %= P;
    return value < 0 ? value + P : value;
}

std::vector<Point> antipodal_classes() {
    std::set<Point> classes;
    for (int x = 0; x < P; ++x) {
        for (int y = 0; y < P; ++y) {
            if (x == 0 && y == 0) continue;
            const Point point{x, y};
            const Point negative{mod(-x), mod(-y)};
            classes.insert(std::min(point, negative));
        }
    }
    return {classes.begin(), classes.end()};
}

std::array<Point, DIRECTIONS> directions() {
    std::array<Point, DIRECTIONS> out{};
    for (int slope = 0; slope < P; ++slope) out[slope] = {1, slope};
    out[P] = {0, 1};
    return out;
}

std::uint32_t encode_pair(int a, int b) {
    if (!(0 <= a && a < b && b < POINTS)) throw std::runtime_error("invalid pair");
    return static_cast<std::uint32_t>(a) | (static_cast<std::uint32_t>(b) << 9);
}

std::array<int, 2> decode_pair(std::uint32_t code) {
    return {static_cast<int>(code & 511U), static_cast<int>((code >> 9) & 511U)};
}

bool disjoint(const std::array<int, 2>& a, const std::array<int, 2>& b) {
    return a[0] != b[0] && a[0] != b[1] && a[1] != b[0] && a[1] != b[1];
}

double seconds_since(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start
    ).count();
}

void print_int_array(const auto& values, std::size_t count) {
    std::cout << "[";
    for (std::size_t i = 0; i < count; ++i) {
        if (i) std::cout << ",";
        std::cout << values[i];
    }
    std::cout << "]";
}

}  // namespace

int main() {
    const auto started = std::chrono::steady_clock::now();
    const auto points = antipodal_classes();
    const auto dirs = directions();
    if (points.size() != POINTS) throw std::runtime_error("point count changed");

    std::array<int, P> square_index{};
    square_index.fill(-1);
    std::vector<int> squares;
    for (int value = 1; value < P; ++value) squares.push_back(value * value % P);
    std::sort(squares.begin(), squares.end());
    squares.erase(std::unique(squares.begin(), squares.end()), squares.end());
    if (squares.size() != FIBRES - 1) throw std::runtime_error("square count changed");
    for (int i = 0; i < FIBRES - 1; ++i) square_index[squares[i]] = i;

    std::array<std::array<std::uint8_t, POINTS>, DIRECTIONS> category{};
    std::array<std::array<int, FIBRES>, DIRECTIONS> fibre_sizes{};
    for (int d = 0; d < DIRECTIONS; ++d) {
        for (int i = 0; i < POINTS; ++i) {
            const int value = mod(dirs[d].x * points[i].x + dirs[d].y * points[i].y);
            const int fibre = value == 0 ? FIBRES - 1 : square_index[value * value % P];
            if (fibre < 0 || fibre >= FIBRES) throw std::runtime_error("bad fibre");
            category[d][i] = static_cast<std::uint8_t>(fibre);
            ++fibre_sizes[d][fibre];
        }
    }
    for (int d = 0; d < DIRECTIONS; ++d) {
        for (int f = 0; f < FIBRES - 1; ++f) {
            if (fibre_sizes[d][f] != P) throw std::runtime_error("affine fibre size changed");
        }
        if (fibre_sizes[d][FIBRES - 1] != 15) throw std::runtime_error("kernel size changed");
    }

    std::array<std::uint64_t, POINTS> point_signature{};
    for (int i = 0; i < POINTS; ++i) {
        point_signature[i] =
            (1ULL << category[0][i])
            ^ (1ULL << (16 + category[1][i]))
            ^ (1ULL << (32 + category[31][i]));
    }

    std::vector<PairRecord> pairs;
    pairs.reserve(EXPECTED_PAIRS);
    for (int i = 0; i < POINTS; ++i) {
        for (int j = i + 1; j < POINTS; ++j) {
            pairs.push_back({point_signature[i] ^ point_signature[j], encode_pair(i, j)});
        }
    }
    if (pairs.size() != EXPECTED_PAIRS) throw std::runtime_error("pair count changed");
    std::sort(pairs.begin(), pairs.end(), [](const PairRecord& a, const PairRecord& b) {
        if (a.signature != b.signature) return a.signature < b.signature;
        return a.code < b.code;
    });

    std::uint64_t collision_groups = 0;
    std::uint64_t equal_signature_pairs = 0;
    std::uint64_t overlapping_pair_pairs = 0;
    std::uint64_t disjoint_pair_pairs = 0;
    std::uint64_t max_group_size = 0;
    int best_silent = -1;
    std::array<int, 4> best_support{};
    std::array<int, DIRECTIONS> best_silent_dirs{};
    int best_silent_count = 0;

    for (std::size_t first = 0; first < pairs.size();) {
        std::size_t last = first + 1;
        while (last < pairs.size() && pairs[last].signature == pairs[first].signature) ++last;
        const std::uint64_t group_size = last - first;
        max_group_size = std::max(max_group_size, group_size);
        if (group_size >= 2) {
            ++collision_groups;
            equal_signature_pairs += group_size * (group_size - 1) / 2;
            for (std::size_t left = first; left + 1 < last; ++left) {
                const auto a = decode_pair(pairs[left].code);
                for (std::size_t right = left + 1; right < last; ++right) {
                    const auto b = decode_pair(pairs[right].code);
                    if (!disjoint(a, b)) {
                        ++overlapping_pair_pairs;
                        continue;
                    }
                    ++disjoint_pair_pairs;
                    std::array<int, 4> support{a[0], a[1], b[0], b[1]};
                    std::sort(support.begin(), support.end());
                    int silent_count = 0;
                    std::array<int, DIRECTIONS> silent_dirs{};
                    for (int d = 0; d < DIRECTIONS; ++d) {
                        std::uint16_t parity = 0;
                        for (const int point : support) {
                            parity ^= static_cast<std::uint16_t>(1U << category[d][point]);
                        }
                        if (parity == 0) silent_dirs[silent_count++] = d;
                    }
                    for (const int fixed : {0, 1, 31}) {
                        if (!std::binary_search(
                                silent_dirs.begin(), silent_dirs.begin() + silent_count, fixed
                            )) {
                            throw std::runtime_error("signature collision replay failed");
                        }
                    }
                    if (silent_count > best_silent) {
                        best_silent = silent_count;
                        best_support = support;
                        best_silent_dirs = silent_dirs;
                        best_silent_count = silent_count;
                    }
                    if (silent_count >= 5) {
                        std::cout << "{\n  \"status\":\"COUNTEREXAMPLE_FOUND\",\n"
                                  << "  \"support_indices\":";
                        print_int_array(support, support.size());
                        std::cout << ",\n  \"silent_direction_indices\":";
                        print_int_array(silent_dirs, silent_count);
                        std::cout << ",\n  \"wall_seconds\":" << seconds_since(started) << "\n}\n";
                        return 2;
                    }
                }
            }
        }
        first = last;
    }

    std::cout << "{\n  \"p\":31,\n  \"point_weight\":4,\n"
              << "  \"status\":\"EXHAUSTIVE_NO_COUNTEREXAMPLE\",\n"
              << "  \"fixed_silent_direction_indices\":[0,1,31],\n"
              << "  \"pairs_generated\":" << pairs.size() << ",\n"
              << "  \"collision_groups\":" << collision_groups << ",\n"
              << "  \"max_collision_group_size\":" << max_group_size << ",\n"
              << "  \"equal_signature_pair_pairs\":" << equal_signature_pairs << ",\n"
              << "  \"overlapping_pair_pairs_skipped\":" << overlapping_pair_pairs << ",\n"
              << "  \"disjoint_pair_pairs_checked\":" << disjoint_pair_pairs << ",\n"
              << "  \"four_sets_counted_with_partition_multiplicity_removed\":"
              << disjoint_pair_pairs / 3 << ",\n"
              << "  \"best_silent_direction_count\":" << best_silent << ",\n"
              << "  \"best_support_indices\":";
    print_int_array(best_support, best_support.size());
    std::cout << ",\n  \"best_silent_direction_indices\":";
    print_int_array(best_silent_dirs, best_silent_count);
    std::cout << ",\n  \"coverage_argument\":\"every four-set silent in directions 0,1,31 has zero total 48-bit fibre signature, so each of its three unordered 2+2 partitions is a disjoint equal-signature pair pair\",\n"
              << "  \"wall_seconds\":" << seconds_since(started) << "\n}\n";
    return 0;
}
