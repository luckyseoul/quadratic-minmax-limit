// Exact cycle-space meet-in-the-middle decision procedure for the p=31,
// weight-eight grouped uncertainty case.
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
constexpr std::uint64_t EXPECTED_HALVES = 23909340ULL;

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

bool disjoint(const std::array<int, 4>& a, const std::array<int, 4>& b) {
    int i = 0;
    int j = 0;
    while (i < 4 && j < 4) {
        if (a[i] == b[j]) return false;
        if (a[i] < b[j]) ++i;
        else ++j;
    }
    return true;
}

std::array<int, 8> merged_support(
    const std::array<int, 4>& a,
    const std::array<int, 4>& b
) {
    std::array<int, 8> out{};
    std::merge(a.begin(), a.end(), b.begin(), b.end(), out.begin());
    return out;
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
    if (argc != 2) {
        std::cerr << "usage: p31_group_s8_cycle_mitm FOURTH_SILENT_DIRECTION\n";
        return 64;
    }
    const int lambda = std::atoi(argv[1]);
    if (lambda < 2 || lambda > 30) {
        std::cerr << "fourth direction must lie in [2,30]\n";
        return 64;
    }
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

    // Category 15 is the kernel fibre.  For an even support, parity zero in
    // all sixteen categories is equivalent to silence in the fifteen
    // nonzero-square blocks used by the grouped map.
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
        if (fibre_sizes[d][FIBRES - 1] != (P - 1) / 2) {
            throw std::runtime_error("kernel fibre size changed");
        }
    }

    std::array<std::uint64_t, POINTS> half_point_signature{};
    std::array<std::uint64_t, POINTS> four_direction_point_signature{};
    for (int i = 0; i < POINTS; ++i) {
        half_point_signature[i] =
            (1ULL << category[1][i])
            ^ (1ULL << (16 + category[31][i]))
            ^ (1ULL << (32 + category[lambda][i]));
        four_direction_point_signature[i] =
            (1ULL << category[0][i])
            ^ (1ULL << (16 + category[1][i]))
            ^ (1ULL << (32 + category[31][i]))
            ^ (1ULL << (48 + category[lambda][i]));
    }

    // Independent generic-case sanity check: all unordered pair syndromes
    // for the four fixed directions must be distinct.
    std::vector<std::uint64_t> pair_signatures;
    pair_signatures.reserve(POINTS * (POINTS - 1) / 2);
    for (int i = 0; i < POINTS; ++i) {
        for (int j = i + 1; j < POINTS; ++j) {
            pair_signatures.push_back(
                four_direction_point_signature[i] ^ four_direction_point_signature[j]
            );
        }
    }
    std::sort(pair_signatures.begin(), pair_signatures.end());
    if (std::adjacent_find(pair_signatures.begin(), pair_signatures.end())
        != pair_signatures.end()) {
        throw std::runtime_error("four-direction pair signatures are not injective");
    }
    pair_signatures.clear();
    pair_signatures.shrink_to_fit();

    std::array<std::vector<int>, FIBRES> rows;
    for (int i = 0; i < POINTS; ++i) rows[category[0][i]].push_back(i);
    std::array<int, FIBRES> sorted_row_sizes{};
    for (int r = 0; r < FIBRES; ++r) sorted_row_sizes[r] = rows[r].size();
    std::sort(sorted_row_sizes.begin(), sorted_row_sizes.end());
    if (sorted_row_sizes[0] != 15 || sorted_row_sizes[1] != 31
        || sorted_row_sizes[15] != 31) {
        throw std::runtime_error("row-size profile changed");
    }

    std::vector<HalfRecord> halves;
    halves.reserve(EXPECTED_HALVES);
    // Four selected points in one direction-zero fibre.
    for (const auto& row : rows) {
        for (std::size_t a = 0; a + 3 < row.size(); ++a) {
            for (std::size_t b = a + 1; b + 2 < row.size(); ++b) {
                for (std::size_t c = b + 1; c + 1 < row.size(); ++c) {
                    const std::uint64_t abc = half_point_signature[row[a]]
                        ^ half_point_signature[row[b]] ^ half_point_signature[row[c]];
                    for (std::size_t d = c + 1; d < row.size(); ++d) {
                        halves.push_back({
                            abc ^ half_point_signature[row[d]],
                            encode_half({row[a], row[b], row[c], row[d]}),
                        });
                    }
                }
            }
        }
    }
    // Two selected points in each of two distinct direction-zero fibres.
    for (int r = 0; r < FIBRES; ++r) {
        for (int s = r + 1; s < FIBRES; ++s) {
            const auto& left = rows[r];
            const auto& right = rows[s];
            for (std::size_t a = 0; a + 1 < left.size(); ++a) {
                for (std::size_t b = a + 1; b < left.size(); ++b) {
                    const std::uint64_t left_signature =
                        half_point_signature[left[a]] ^ half_point_signature[left[b]];
                    for (std::size_t c = 0; c + 1 < right.size(); ++c) {
                        const std::uint64_t prefix =
                            left_signature ^ half_point_signature[right[c]];
                        for (std::size_t d = c + 1; d < right.size(); ++d) {
                            halves.push_back({
                                prefix ^ half_point_signature[right[d]],
                                encode_half({left[a], left[b], right[c], right[d]}),
                            });
                        }
                    }
                }
            }
        }
    }
    if (halves.size() != EXPECTED_HALVES) throw std::runtime_error("half count changed");
    std::cerr << "lambda=" << lambda << " generated " << halves.size()
              << " row-even halves in " << seconds_since(started)
              << " seconds; record bytes=" << sizeof(HalfRecord) * halves.size() << "\n";

    std::sort(halves.begin(), halves.end(), [](const HalfRecord& a, const HalfRecord& b) {
        if (a.signature != b.signature) return a.signature < b.signature;
        return a.code < b.code;
    });
    const double sorted_seconds = seconds_since(started);
    std::cerr << "lambda=" << lambda << " sorted in " << sorted_seconds << " seconds\n";

    std::uint64_t collision_groups = 0;
    std::uint64_t equal_signature_pairs = 0;
    std::uint64_t disjoint_half_pairs = 0;
    std::uint64_t overlapping_half_pairs = 0;
    std::uint64_t max_group_size = 0;
    int best_silent = -1;
    std::array<int, 8> best_support{};
    std::array<int, DIRECTIONS> best_silent_directions{};
    int best_silent_direction_count = 0;

    for (std::size_t first = 0; first < halves.size();) {
        std::size_t last = first + 1;
        while (last < halves.size() && halves[last].signature == halves[first].signature) ++last;
        const std::uint64_t group_size = last - first;
        max_group_size = std::max(max_group_size, group_size);
        if (group_size >= 2) {
            ++collision_groups;
            equal_signature_pairs += group_size * (group_size - 1) / 2;
            for (std::size_t left = first; left + 1 < last; ++left) {
                const auto a = decode_half(halves[left].code);
                if (!(a[0] < a[1] && a[1] < a[2] && a[2] < a[3]
                      && a[3] < POINTS)) {
                    throw std::runtime_error("bad decoded half");
                }
                for (std::size_t right = left + 1; right < last; ++right) {
                    const auto b = decode_half(halves[right].code);
                    if (!disjoint(a, b)) {
                        ++overlapping_half_pairs;
                        continue;
                    }
                    ++disjoint_half_pairs;
                    const auto support = merged_support(a, b);
                    int silent_count = 0;
                    std::array<int, DIRECTIONS> silent_directions{};
                    for (int d = 0; d < DIRECTIONS; ++d) {
                        std::uint16_t parity = 0;
                        for (const int point : support) {
                            parity ^= static_cast<std::uint16_t>(1U << category[d][point]);
                        }
                        if (parity == 0) silent_directions[silent_count++] = d;
                    }
                    for (const int fixed : {0, 1, 31, lambda}) {
                        if (!std::binary_search(
                                silent_directions.begin(),
                                silent_directions.begin() + silent_count,
                                fixed
                            )) {
                            throw std::runtime_error("half collision failed fixed-direction replay");
                        }
                    }
                    if (silent_count > best_silent) {
                        best_silent = silent_count;
                        best_support = support;
                        best_silent_directions = silent_directions;
                        best_silent_direction_count = silent_count;
                    }
                    if (silent_count >= 9) {
                        std::cout << "{\n  \"p\":31,\n  \"point_weight\":8,\n"
                                  << "  \"fourth_silent_direction_index\":" << lambda << ",\n"
                                  << "  \"status\":\"COUNTEREXAMPLE_FOUND\",\n"
                                  << "  \"support_indices\":";
                        print_int_array(support, support.size());
                        std::cout << ",\n  \"support_points\":[";
                        for (std::size_t q = 0; q < support.size(); ++q) {
                            if (q) std::cout << ",";
                            std::cout << "[" << points[support[q]].x << ","
                                      << points[support[q]].y << "]";
                        }
                        std::cout << "],\n  \"silent_direction_indices\":";
                        print_int_array(silent_directions, silent_count);
                        std::cout << ",\n  \"silent_direction_count\":" << silent_count
                                  << ",\n  \"active_direction_count\":" << DIRECTIONS - silent_count
                                  << ",\n  \"group_branch_weight\":" << 8 + DIRECTIONS - silent_count
                                  << ",\n  \"disjoint_half_pairs_checked_before_witness\":"
                                  << disjoint_half_pairs
                                  << ",\n  \"wall_seconds\":" << seconds_since(started) << "\n}\n";
                        return 2;
                    }
                }
            }
        }
        first = last;
    }

    std::cout << "{\n  \"p\":31,\n  \"point_weight\":8,\n"
              << "  \"fourth_silent_direction_index\":" << lambda << ",\n"
              << "  \"status\":\"EXHAUSTIVE_NO_COUNTEREXAMPLE\",\n"
              << "  \"fixed_silent_direction_indices\":[0,1,31," << lambda << "],\n"
              << "  \"pair_signatures_injective\":true,\n"
              << "  \"row_even_halves_generated\":" << halves.size() << ",\n"
              << "  \"collision_groups\":" << collision_groups << ",\n"
              << "  \"max_collision_group_size\":" << max_group_size << ",\n"
              << "  \"equal_signature_pairs\":" << equal_signature_pairs << ",\n"
              << "  \"overlapping_half_pairs_skipped\":" << overlapping_half_pairs << ",\n"
              << "  \"disjoint_half_pairs_checked\":" << disjoint_half_pairs << ",\n"
              << "  \"best_silent_direction_count\":" << best_silent << ",\n"
              << "  \"best_support_indices\":";
    print_int_array(best_support, best_support.size());
    std::cout << ",\n  \"best_support_points\":[";
    for (std::size_t q = 0; q < best_support.size(); ++q) {
        if (q) std::cout << ",";
        std::cout << "[" << points[best_support[q]].x << ","
                  << points[best_support[q]].y << "]";
    }
    std::cout << "],\n  \"best_silent_direction_indices\":";
    print_int_array(best_silent_directions, best_silent_direction_count);
    std::cout << ",\n  \"coverage_argument\":\"an eight-set silent in direction 0 has one of the even row-degree types 8, 6+2, 4+4, 4+2+2, or 2+2+2+2; each type partitions into two enumerated row-even four-sets, and silence in directions 1,31,lambda makes their 48-bit signatures equal\",\n"
              << "  \"sorted_seconds\":" << sorted_seconds << ",\n"
              << "  \"wall_seconds\":" << seconds_since(started) << "\n}\n";
    return 0;
}
