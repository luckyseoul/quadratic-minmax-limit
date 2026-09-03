// Exact meet-in-the-middle decision procedure for the p=31, weight-six
// grouped uncertainty case after normalizing three silent directions.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

constexpr int P = 31;
constexpr int POINTS = 480;
constexpr int DIRECTIONS = 32;
constexpr int BLOCKS_PER_DIRECTION = 15;
constexpr std::uint64_t EXPECTED_TRIPLES = 18316960ULL;

struct Point {
    int x;
    int y;
    friend bool operator<(const Point& a, const Point& b) {
        return std::tie(a.x, a.y) < std::tie(b.x, b.y);
    }
    friend bool operator==(const Point& a, const Point& b) {
        return a.x == b.x && a.y == b.y;
    }
};

struct TripleRecord {
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
            Point point{x, y};
            Point negative{mod(-x), mod(-y)};
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

std::uint32_t encode_triple(int i, int j, int k) {
    return static_cast<std::uint32_t>(i)
        | (static_cast<std::uint32_t>(j) << 9)
        | (static_cast<std::uint32_t>(k) << 18);
}

std::array<int, 3> decode_triple(std::uint32_t code) {
    return {
        static_cast<int>(code & 511U),
        static_cast<int>((code >> 9) & 511U),
        static_cast<int>((code >> 18) & 511U),
    };
}

bool disjoint(const std::array<int, 3>& a, const std::array<int, 3>& b) {
    int i = 0;
    int j = 0;
    while (i < 3 && j < 3) {
        if (a[i] == b[j]) return false;
        if (a[i] < b[j]) ++i;
        else ++j;
    }
    return true;
}

std::array<int, 6> merged_support(
    const std::array<int, 3>& a,
    const std::array<int, 3>& b
) {
    std::array<int, 6> out{};
    std::merge(a.begin(), a.end(), b.begin(), b.end(), out.begin());
    return out;
}

double seconds_since(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start
    ).count();
}

void print_int_array(const auto& values) {
    std::cout << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
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
    if (squares.size() != BLOCKS_PER_DIRECTION) {
        throw std::runtime_error("square count changed");
    }
    for (int i = 0; i < BLOCKS_PER_DIRECTION; ++i) square_index[squares[i]] = i;

    std::array<std::array<std::int8_t, POINTS>, DIRECTIONS> block{};
    std::array<std::array<int, BLOCKS_PER_DIRECTION>, DIRECTIONS> block_sizes{};
    for (int d = 0; d < DIRECTIONS; ++d) {
        for (int i = 0; i < POINTS; ++i) {
            const int value = mod(dirs[d].x * points[i].x + dirs[d].y * points[i].y);
            if (value == 0) {
                block[d][i] = -1;
            } else {
                const int index = square_index[value * value % P];
                if (index < 0) throw std::runtime_error("nonsquare squared value");
                block[d][i] = static_cast<std::int8_t>(index);
                ++block_sizes[d][index];
            }
        }
    }
    for (int d = 0; d < DIRECTIONS; ++d) {
        for (int b = 0; b < BLOCKS_PER_DIRECTION; ++b) {
            if (block_sizes[d][b] != P) throw std::runtime_error("block size changed");
        }
    }
    for (int i = 0; i < POINTS; ++i) {
        int nonkernel = 0;
        for (int d = 0; d < DIRECTIONS; ++d) nonkernel += block[d][i] >= 0;
        if (nonkernel != P) throw std::runtime_error("point incidence changed");
    }

    constexpr std::array<int, 3> FIXED{0, 1, 31};
    std::array<std::uint64_t, POINTS> point_signature{};
    for (int i = 0; i < POINTS; ++i) {
        std::uint64_t signature = 0;
        for (int slot = 0; slot < 3; ++slot) {
            const int b = block[FIXED[slot]][i];
            if (b >= 0) signature ^= (1ULL << (slot * BLOCKS_PER_DIRECTION + b));
        }
        point_signature[i] = signature;
    }

    std::vector<TripleRecord> records;
    records.reserve(EXPECTED_TRIPLES);
    for (int i = 0; i < POINTS - 2; ++i) {
        for (int j = i + 1; j < POINTS - 1; ++j) {
            const std::uint64_t pair_signature = point_signature[i] ^ point_signature[j];
            for (int k = j + 1; k < POINTS; ++k) {
                records.push_back({
                    pair_signature ^ point_signature[k],
                    encode_triple(i, j, k),
                });
            }
        }
    }
    if (records.size() != EXPECTED_TRIPLES) {
        throw std::runtime_error("triple count changed");
    }
    std::cerr << "generated " << records.size() << " triples in "
              << seconds_since(started) << " seconds; record bytes="
              << sizeof(TripleRecord) * records.size() << "\n";

    std::sort(records.begin(), records.end(), [](const TripleRecord& a, const TripleRecord& b) {
        if (a.signature != b.signature) return a.signature < b.signature;
        return a.code < b.code;
    });
    const double sorted_seconds = seconds_since(started);
    std::cerr << "sorted in " << sorted_seconds << " seconds\n";

    std::uint64_t collision_groups = 0;
    std::uint64_t equal_signature_pairs = 0;
    std::uint64_t disjoint_triple_pairs = 0;
    std::uint64_t overlapping_triple_pairs = 0;
    std::uint64_t max_group_size = 0;
    int best_silent = -1;
    std::array<int, 6> best_support{};
    std::array<int, DIRECTIONS> best_silent_directions{};
    int best_silent_direction_count = 0;

    for (std::size_t first = 0; first < records.size();) {
        std::size_t last = first + 1;
        while (last < records.size() && records[last].signature == records[first].signature) ++last;
        const std::uint64_t group_size = last - first;
        max_group_size = std::max(max_group_size, group_size);
        if (group_size >= 2) {
            ++collision_groups;
            equal_signature_pairs += group_size * (group_size - 1) / 2;
            for (std::size_t left = first; left + 1 < last; ++left) {
                const auto a = decode_triple(records[left].code);
                if (!(a[0] < a[1] && a[1] < a[2] && a[2] < POINTS)) {
                    throw std::runtime_error("bad decoded triple");
                }
                for (std::size_t right = left + 1; right < last; ++right) {
                    const auto b = decode_triple(records[right].code);
                    if (!disjoint(a, b)) {
                        ++overlapping_triple_pairs;
                        continue;
                    }
                    ++disjoint_triple_pairs;
                    const auto support = merged_support(a, b);
                    int silent_count = 0;
                    std::array<int, DIRECTIONS> silent_directions{};
                    for (int d = 0; d < DIRECTIONS; ++d) {
                        std::uint16_t parity = 0;
                        for (const int point : support) {
                            const int index = block[d][point];
                            if (index >= 0) parity ^= static_cast<std::uint16_t>(1U << index);
                        }
                        if (parity == 0) silent_directions[silent_count++] = d;
                    }
                    if (!std::binary_search(
                            silent_directions.begin(),
                            silent_directions.begin() + silent_count,
                            0
                        ) || !std::binary_search(
                            silent_directions.begin(),
                            silent_directions.begin() + silent_count,
                            1
                        ) || !std::binary_search(
                            silent_directions.begin(),
                            silent_directions.begin() + silent_count,
                            31
                        )) {
                        throw std::runtime_error("signature collision failed fixed-direction replay");
                    }
                    if (silent_count > best_silent) {
                        best_silent = silent_count;
                        best_support = support;
                        best_silent_directions = silent_directions;
                        best_silent_direction_count = silent_count;
                    }
                    if (silent_count >= 7) {
                        std::cout << "{\n"
                                  << "  \"p\":31,\n"
                                  << "  \"point_weight\":6,\n"
                                  << "  \"status\":\"COUNTEREXAMPLE_FOUND\",\n"
                                  << "  \"support_indices\":";
                        print_int_array(support);
                        std::cout << ",\n  \"support_points\":[";
                        for (std::size_t q = 0; q < support.size(); ++q) {
                            if (q) std::cout << ",";
                            std::cout << "[" << points[support[q]].x << "," << points[support[q]].y << "]";
                        }
                        std::cout << "],\n  \"silent_direction_indices\":[";
                        for (int q = 0; q < silent_count; ++q) {
                            if (q) std::cout << ",";
                            std::cout << silent_directions[q];
                        }
                        std::cout << "],\n"
                                  << "  \"silent_direction_count\":" << silent_count << ",\n"
                                  << "  \"active_direction_count\":" << DIRECTIONS - silent_count << ",\n"
                                  << "  \"group_branch_weight\":" << 6 + DIRECTIONS - silent_count << ",\n"
                                  << "  \"triples_generated\":" << records.size() << ",\n"
                                  << "  \"disjoint_triple_pairs_checked_before_witness\":" << disjoint_triple_pairs << ",\n"
                                  << "  \"wall_seconds\":" << seconds_since(started) << "\n"
                                  << "}\n";
                        return 2;
                    }
                }
            }
        }
        first = last;
    }

    std::cout << "{\n"
              << "  \"p\":31,\n"
              << "  \"point_weight\":6,\n"
              << "  \"status\":\"EXHAUSTIVE_NO_COUNTEREXAMPLE\",\n"
              << "  \"fixed_silent_direction_indices\":[0,1,31],\n"
              << "  \"triples_generated\":" << records.size() << ",\n"
              << "  \"collision_groups\":" << collision_groups << ",\n"
              << "  \"max_collision_group_size\":" << max_group_size << ",\n"
              << "  \"equal_signature_pairs\":" << equal_signature_pairs << ",\n"
              << "  \"overlapping_triple_pairs_skipped\":" << overlapping_triple_pairs << ",\n"
              << "  \"disjoint_triple_pairs_checked\":" << disjoint_triple_pairs << ",\n"
              << "  \"best_silent_direction_count\":" << best_silent << ",\n"
              << "  \"best_support_indices\":";
    print_int_array(best_support);
    std::cout << ",\n  \"best_support_points\":[";
    for (std::size_t q = 0; q < best_support.size(); ++q) {
        if (q) std::cout << ",";
        std::cout << "[" << points[best_support[q]].x << "," << points[best_support[q]].y << "]";
    }
    std::cout << "],\n  \"best_silent_direction_indices\":[";
    for (int q = 0; q < best_silent_direction_count; ++q) {
        if (q) std::cout << ",";
        std::cout << best_silent_directions[q];
    }
    std::cout << "],\n"
              << "  \"coverage_argument\":\"every six-set silent in directions 0,1,31 has zero total 45-bit signature; each of its ten 3+3 partitions is therefore a disjoint equal-signature triple pair\",\n"
              << "  \"sorted_seconds\":" << sorted_seconds << ",\n"
              << "  \"wall_seconds\":" << seconds_since(started) << "\n"
              << "}\n";
    return 0;
}
