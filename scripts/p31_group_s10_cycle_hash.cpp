// Exact scalar-orbit/cycle-space hash decision procedure for the p=31,
// weight-ten grouped uncertainty case.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

constexpr int P = 31;
constexpr int POINTS = 480;
constexpr int DIRECTIONS = 32;
constexpr int FIBRES = 16;
constexpr std::uint64_t EXPECTED_HALVES = 23909340ULL;
constexpr std::uint32_t NO_GROUP = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t EMPTY_KEY = std::numeric_limits<std::uint64_t>::max();

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

struct HalfRecord {
    std::uint64_t signature;
    std::uint64_t code;
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

std::uint64_t encode_half(std::array<int, 4> values) {
    std::sort(values.begin(), values.end());
    if (!(values[0] < values[1] && values[1] < values[2]
          && values[2] < values[3] && values[3] < POINTS)) {
        throw std::runtime_error("invalid half support");
    }
    return static_cast<std::uint64_t>(values[0])
        | (static_cast<std::uint64_t>(values[1]) << 9)
        | (static_cast<std::uint64_t>(values[2]) << 18)
        | (static_cast<std::uint64_t>(values[3]) << 27);
}

std::array<int, 4> decode_half(std::uint64_t code) {
    return {
        static_cast<int>(code & 511ULL),
        static_cast<int>((code >> 9) & 511ULL),
        static_cast<int>((code >> 18) & 511ULL),
        static_cast<int>((code >> 27) & 511ULL),
    };
}

std::uint32_t encode_pair(int a, int b) {
    if (a > b) std::swap(a, b);
    if (!(0 <= a && a < b && b < POINTS)) throw std::runtime_error("invalid pair");
    return static_cast<std::uint32_t>(a) | (static_cast<std::uint32_t>(b) << 9);
}

std::array<int, 2> decode_pair(std::uint32_t code) {
    return {static_cast<int>(code & 511U), static_cast<int>((code >> 9) & 511U)};
}

bool halves_disjoint(const std::array<int, 4>& a, const std::array<int, 4>& b) {
    int i = 0;
    int j = 0;
    while (i < 4 && j < 4) {
        if (a[i] == b[j]) return false;
        if (a[i] < b[j]) ++i;
        else ++j;
    }
    return true;
}

bool avoids_pair(const std::array<int, 4>& half, const std::array<int, 2>& pair) {
    return !std::binary_search(half.begin(), half.end(), pair[0])
        && !std::binary_search(half.begin(), half.end(), pair[1]);
}

std::array<int, 10> merged_support(
    const std::array<int, 4>& a,
    const std::array<int, 4>& b,
    const std::array<int, 2>& pair
) {
    std::array<int, 10> out{a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3],
                            pair[0], pair[1]};
    std::sort(out.begin(), out.end());
    return out;
}

std::uint64_t hash64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
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

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: p31_group_s10_cycle_hash FOURTH_SILENT [PAIR_ORBIT_LIMIT]\n";
        return 64;
    }
    const int lambda = std::atoi(argv[1]);
    const int orbit_limit = argc == 3 ? std::atoi(argv[2]) : 0;
    if (lambda < 2 || lambda > 30 || orbit_limit < 0) return 64;
    const auto started = std::chrono::steady_clock::now();
    const auto points = antipodal_classes();
    const auto dirs = directions();
    if (points.size() != POINTS) throw std::runtime_error("point count changed");

    std::map<Point, int> point_index;
    for (int i = 0; i < POINTS; ++i) point_index[points[i]] = i;

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
    std::array<std::uint64_t, POINTS> four_direction_signature{};
    for (int i = 0; i < POINTS; ++i) {
        point_signature[i] =
            (1ULL << category[1][i])
            ^ (1ULL << (16 + category[31][i]))
            ^ (1ULL << (32 + category[lambda][i]));
        four_direction_signature[i] =
            (1ULL << category[0][i])
            ^ (1ULL << (16 + category[1][i]))
            ^ (1ULL << (32 + category[31][i]))
            ^ (1ULL << (48 + category[lambda][i]));
    }

    std::vector<std::uint64_t> all_pair_signatures;
    all_pair_signatures.reserve(POINTS * (POINTS - 1) / 2);
    for (int i = 0; i < POINTS; ++i) {
        for (int j = i + 1; j < POINTS; ++j) {
            all_pair_signatures.push_back(
                four_direction_signature[i] ^ four_direction_signature[j]
            );
        }
    }
    std::sort(all_pair_signatures.begin(), all_pair_signatures.end());
    if (std::adjacent_find(all_pair_signatures.begin(), all_pair_signatures.end())
        != all_pair_signatures.end()) {
        throw std::runtime_error("four-direction pair signatures are not injective");
    }
    all_pair_signatures.clear();
    all_pair_signatures.shrink_to_fit();

    std::array<std::vector<int>, FIBRES> rows;
    for (int i = 0; i < POINTS; ++i) rows[category[0][i]].push_back(i);

    std::vector<HalfRecord> halves;
    halves.reserve(EXPECTED_HALVES);
    for (const auto& row : rows) {
        for (std::size_t a = 0; a + 3 < row.size(); ++a) {
            for (std::size_t b = a + 1; b + 2 < row.size(); ++b) {
                for (std::size_t c = b + 1; c + 1 < row.size(); ++c) {
                    const std::uint64_t abc = point_signature[row[a]]
                        ^ point_signature[row[b]] ^ point_signature[row[c]];
                    for (std::size_t d = c + 1; d < row.size(); ++d) {
                        halves.push_back({
                            abc ^ point_signature[row[d]],
                            encode_half({row[a], row[b], row[c], row[d]}),
                        });
                    }
                }
            }
        }
    }
    for (int r = 0; r < FIBRES; ++r) {
        for (int s = r + 1; s < FIBRES; ++s) {
            const auto& left = rows[r];
            const auto& right = rows[s];
            for (std::size_t a = 0; a + 1 < left.size(); ++a) {
                for (std::size_t b = a + 1; b < left.size(); ++b) {
                    const std::uint64_t left_signature =
                        point_signature[left[a]] ^ point_signature[left[b]];
                    for (std::size_t c = 0; c + 1 < right.size(); ++c) {
                        const std::uint64_t prefix =
                            left_signature ^ point_signature[right[c]];
                        for (std::size_t d = c + 1; d < right.size(); ++d) {
                            halves.push_back({
                                prefix ^ point_signature[right[d]],
                                encode_half({left[a], left[b], right[c], right[d]}),
                            });
                        }
                    }
                }
            }
        }
    }
    if (halves.size() != EXPECTED_HALVES) throw std::runtime_error("half count changed");
    std::sort(halves.begin(), halves.end(), [](const HalfRecord& a, const HalfRecord& b) {
        if (a.signature != b.signature) return a.signature < b.signature;
        return a.code < b.code;
    });

    std::vector<std::uint64_t> unique_signatures;
    std::vector<std::uint32_t> group_starts;
    std::vector<std::uint32_t> group_counts;
    unique_signatures.reserve(halves.size());
    group_starts.reserve(halves.size());
    group_counts.reserve(halves.size());
    for (std::size_t first = 0; first < halves.size();) {
        std::size_t last = first + 1;
        while (last < halves.size() && halves[last].signature == halves[first].signature) ++last;
        unique_signatures.push_back(halves[first].signature);
        group_starts.push_back(static_cast<std::uint32_t>(first));
        group_counts.push_back(static_cast<std::uint32_t>(last - first));
        first = last;
    }

    std::size_t table_capacity = 1;
    while (table_capacity < 2 * unique_signatures.size()) table_capacity <<= 1;
    std::vector<std::uint64_t> table_keys(table_capacity, EMPTY_KEY);
    std::vector<std::uint32_t> table_groups(table_capacity, NO_GROUP);
    const std::size_t table_mask = table_capacity - 1;
    for (std::uint32_t group = 0; group < unique_signatures.size(); ++group) {
        const std::uint64_t key = unique_signatures[group];
        std::size_t slot = hash64(key) & table_mask;
        while (table_keys[slot] != EMPTY_KEY) slot = (slot + 1) & table_mask;
        table_keys[slot] = key;
        table_groups[slot] = group;
    }
    auto lookup = [&](std::uint64_t key) -> std::uint32_t {
        std::size_t slot = hash64(key) & table_mask;
        while (table_keys[slot] != EMPTY_KEY) {
            if (table_keys[slot] == key) return table_groups[slot];
            slot = (slot + 1) & table_mask;
        }
        return NO_GROUP;
    };
    for (std::uint32_t group = 0; group < unique_signatures.size(); ++group) {
        if (lookup(unique_signatures[group]) != group) throw std::runtime_error("hash replay failed");
    }

    std::set<std::uint32_t> unseen_pairs;
    for (const auto& row : rows) {
        for (std::size_t a = 0; a + 1 < row.size(); ++a) {
            for (std::size_t b = a + 1; b < row.size(); ++b) {
                unseen_pairs.insert(encode_pair(row[a], row[b]));
            }
        }
    }
    if (unseen_pairs.size() != 7080) throw std::runtime_error("row-pair count changed");
    std::vector<std::uint32_t> pair_representatives;
    while (!unseen_pairs.empty()) {
        const std::uint32_t seed = *unseen_pairs.begin();
        const auto pair = decode_pair(seed);
        std::set<std::uint32_t> orbit;
        for (int scalar = 1; scalar <= 15; ++scalar) {
            std::array<int, 2> image{};
            for (int q = 0; q < 2; ++q) {
                Point scaled{mod(scalar * points[pair[q]].x), mod(scalar * points[pair[q]].y)};
                const Point negative{mod(-scaled.x), mod(-scaled.y)};
                scaled = std::min(scaled, negative);
                image[q] = point_index.at(scaled);
            }
            orbit.insert(encode_pair(image[0], image[1]));
        }
        if (orbit.size() != 15) throw std::runtime_error("scalar pair action is not free");
        const std::uint32_t representative = *orbit.begin();
        pair_representatives.push_back(representative);
        for (const auto code : orbit) {
            if (!unseen_pairs.erase(code)) throw std::runtime_error("pair orbit overlap");
        }
    }
    if (pair_representatives.size() != 472) throw std::runtime_error("pair orbit count changed");
    std::sort(pair_representatives.begin(), pair_representatives.end());
    const std::size_t reps_to_check = orbit_limit == 0
        ? pair_representatives.size()
        : std::min<std::size_t>(orbit_limit, pair_representatives.size());

    std::cerr << "lambda=" << lambda << " halves=" << halves.size()
              << " unique_signatures=" << unique_signatures.size()
              << " hash_capacity=" << table_capacity
              << " scalar_pair_orbits=" << pair_representatives.size()
              << " checking=" << reps_to_check
              << " setup_seconds=" << seconds_since(started) << "\n";

    std::atomic<bool> found{false};
    std::mutex result_mutex;
    std::uint64_t signature_relations = 0;
    std::uint64_t disjoint_ten_candidates_checked = 0;
    int best_silent = -1;
    std::array<int, 10> best_support{};
    std::array<int, DIRECTIONS> best_silent_directions{};
    int best_silent_count = 0;

#pragma omp parallel
    {
        std::uint64_t local_relations = 0;
        std::uint64_t local_candidates = 0;
        int local_best = -1;
        std::array<int, 10> local_best_support{};
        std::array<int, DIRECTIONS> local_best_dirs{};
        int local_best_count = 0;

#pragma omp for schedule(dynamic, 1)
        for (std::size_t rep_index = 0; rep_index < reps_to_check; ++rep_index) {
            if (found.load(std::memory_order_relaxed)) continue;
            const auto pair = decode_pair(pair_representatives[rep_index]);
            const std::uint64_t target = point_signature[pair[0]] ^ point_signature[pair[1]];
            if (target == 0) throw std::runtime_error("zero row-pair signature");
            for (std::uint32_t left_group = 0;
                 left_group < unique_signatures.size(); ++left_group) {
                if (found.load(std::memory_order_relaxed)) break;
                const std::uint64_t left_signature = unique_signatures[left_group];
                const std::uint64_t right_signature = left_signature ^ target;
                if (left_signature >= right_signature) continue;
                const std::uint32_t right_group = lookup(right_signature);
                if (right_group == NO_GROUP) continue;
                ++local_relations;
                const std::uint32_t left_start = group_starts[left_group];
                const std::uint32_t right_start = group_starts[right_group];
                for (std::uint32_t ai = 0; ai < group_counts[left_group]; ++ai) {
                    const auto a = decode_half(halves[left_start + ai].code);
                    if (!avoids_pair(a, pair)) continue;
                    for (std::uint32_t bi = 0; bi < group_counts[right_group]; ++bi) {
                        const auto b = decode_half(halves[right_start + bi].code);
                        if (!halves_disjoint(a, b) || !avoids_pair(b, pair)) continue;
                        ++local_candidates;
                        const auto support = merged_support(a, b, pair);
                        int silent_count = 0;
                        std::array<int, DIRECTIONS> silent_dirs{};
                        for (int d = 0; d < DIRECTIONS; ++d) {
                            std::uint16_t parity = 0;
                            for (const int point : support) {
                                parity ^= static_cast<std::uint16_t>(1U << category[d][point]);
                            }
                            if (parity == 0) silent_dirs[silent_count++] = d;
                        }
                        for (const int fixed : {0, 1, 31, lambda}) {
                            if (!std::binary_search(
                                    silent_dirs.begin(), silent_dirs.begin() + silent_count, fixed
                                )) {
                                throw std::runtime_error("candidate failed fixed-direction replay");
                            }
                        }
                        if (silent_count > local_best) {
                            local_best = silent_count;
                            local_best_support = support;
                            local_best_dirs = silent_dirs;
                            local_best_count = silent_count;
                        }
                        if (silent_count >= 11) {
                            std::lock_guard<std::mutex> lock(result_mutex);
                            if (!found.exchange(true)) {
                                best_silent = silent_count;
                                best_support = support;
                                best_silent_directions = silent_dirs;
                                best_silent_count = silent_count;
                            }
                            break;
                        }
                    }
                }
            }
        }

#pragma omp critical
        {
            signature_relations += local_relations;
            disjoint_ten_candidates_checked += local_candidates;
            if (!found.load() && local_best > best_silent) {
                best_silent = local_best;
                best_support = local_best_support;
                best_silent_directions = local_best_dirs;
                best_silent_count = local_best_count;
            }
        }
    }

    const bool complete = reps_to_check == pair_representatives.size();
    const char* status = found.load() ? "COUNTEREXAMPLE_FOUND"
        : (complete ? "EXHAUSTIVE_NO_COUNTEREXAMPLE" : "PARTIAL_BENCHMARK");
    std::cout << "{\n  \"p\":31,\n  \"point_weight\":10,\n"
              << "  \"fourth_silent_direction_index\":" << lambda << ",\n"
              << "  \"status\":\"" << status << "\",\n"
              << "  \"row_even_halves\":" << halves.size() << ",\n"
              << "  \"unique_half_signatures\":" << unique_signatures.size() << ",\n"
              << "  \"scalar_pair_orbits_total\":" << pair_representatives.size() << ",\n"
              << "  \"scalar_pair_orbits_checked\":" << reps_to_check << ",\n"
              << "  \"half_signature_relations\":" << signature_relations << ",\n"
              << "  \"disjoint_ten_point_candidates_checked\":"
              << disjoint_ten_candidates_checked << ",\n"
              << "  \"best_silent_direction_count\":" << best_silent << ",\n"
              << "  \"best_support_indices\":";
    print_int_array(best_support, best_support.size());
    std::cout << ",\n  \"best_silent_direction_indices\":";
    print_int_array(best_silent_directions, best_silent_count);
    std::cout << ",\n  \"coverage_argument\":\"remove a same-direction-0-fibre pair from any row-even ten-set; scalar normalization leaves 472 pair orbits, and the row-even eight-point complement splits into two enumerated four-halves whose signatures xor to the pair signature\",\n"
              << "  \"wall_seconds\":" << seconds_since(started) << "\n}\n";
    return found.load() ? 2 : 0;
}
